"""Derive and validate the minimal config contract for an EM-seq Snakemake workflow.

Static analysis of a wrapper .smk: parse the plain-Python preamble (before the
first rule) for config[...] / config.get(...) accesses and alias assignments,
follow include: directives, and scan module bodies for alias subscripts to
recover nested per-entry schema. See
docs/superpowers/specs/2026-06-09-snakemake-config-contract-design.md.

Usage:
    python scripts/emseq_config_contract.py list workflows/test-analysis.smk [--format yaml]
    python scripts/emseq_config_contract.py validate workflows/test-analysis.smk config/test.yaml
"""
from __future__ import annotations

import ast
import os
import pathlib
import re
from dataclasses import dataclass, field
from typing import Any, Optional

_RULE_RE = re.compile(r"^(rule|checkpoint)\b")


@dataclass(frozen=True)
class ConfigPath:
    path: tuple[str, ...]
    requiredness: str  # "mandatory" | "optional" | "conditional"
    default: Any = None
    condition: Optional[tuple[str, ...]] = None
    source: str = ""

    # Equality intentionally on (path, requiredness) so _dedupe collapses same-path entries.
    def __eq__(self, other: object) -> bool:
        return isinstance(other, ConfigPath) and self.path == other.path and \
            self.requiredness == other.requiredness

    def __hash__(self) -> int:
        return hash((self.path, self.requiredness))


@dataclass
class Contract:
    paths: list[ConfigPath] = field(default_factory=list)
    incomplete: bool = False
    warnings: list[str] = field(default_factory=list)
    baked_in: dict[str, Any] = field(default_factory=dict)


class _DynamicKey(Exception):
    pass


def slice_preamble(text: str) -> str:
    """Return the source lines before the first rule/checkpoint, dropping the
    Snakemake `configfile:` directive (not valid plain Python)."""
    out = []
    for line in text.splitlines():
        if _RULE_RE.match(line):
            break
        if line.lstrip().startswith("configfile:"):
            continue
        out.append(line)
    return "\n".join(out)


def _literal_subscript_chain(node: "ast.Subscript") -> Optional[list[str]]:
    """Return the literal key path for a config[...] subscript chain (possibly
    nested), or None if the chain does not bottom out at the `config` name.
    Raise _DynamicKey only when the base IS `config` but a key in the chain is
    non-literal (a genuinely dynamic config access)."""
    segments: list[str] = []
    dynamic = False
    cur: Any = node
    while isinstance(cur, ast.Subscript):
        key = cur.slice
        if isinstance(key, ast.Constant) and isinstance(key.value, str):
            segments.append(key.value)
        else:
            dynamic = True
        cur = cur.value
    if isinstance(cur, ast.Name) and cur.id == "config":
        if dynamic:
            raise _DynamicKey()
        return list(reversed(segments))
    return None


def _config_get_call(node: "ast.Call"):
    """If node is config.get('k', default) or config.get('a',{}).get('b',d),
    return (path_tuple, requiredness, default, condition) or None."""
    if not (isinstance(node.func, ast.Attribute) and node.func.attr == "get"):
        return None
    if not node.args:
        return None
    key = node.args[0]
    if not (isinstance(key, ast.Constant) and isinstance(key.value, str)):
        return None
    default = None
    if len(node.args) > 1 and isinstance(node.args[1], ast.Constant):
        default = node.args[1].value
    recv = node.func.value
    if isinstance(recv, ast.Name) and recv.id == "config":
        return ((key.value,), "optional", default, None)
    # chained: config.get('a', {}).get('b', d)
    if isinstance(recv, ast.Call) and isinstance(recv.func, ast.Attribute) \
            and recv.func.attr == "get" and isinstance(recv.func.value, ast.Name) \
            and recv.func.value.id == "config" and recv.args \
            and isinstance(recv.args[0], ast.Constant):
        parent = recv.args[0].value
        return ((parent, key.value), "conditional", default, (parent,))
    return None


def _comprehension_wildcards(node, source_name: str) -> list[ConfigPath]:
    """For {... for k,v in config['P'].items()} or .values(), find v[<lit>]
    subscripts and emit ConfigPath ('P','*',<lit>) mandatory."""
    out: list[ConfigPath] = []
    for gen in node.generators:
        it = gen.iter
        if not (isinstance(it, ast.Call) and isinstance(it.func, ast.Attribute)
                and it.func.attr in ("items", "values")):
            continue
        recv = it.func.value
        if not isinstance(recv, ast.Subscript):
            continue
        try:
            parent = _literal_subscript_chain(recv)
        except _DynamicKey:
            parent = None
        if not parent:
            continue
        # bind the value loop variable
        if it.func.attr == "items" and isinstance(gen.target, ast.Tuple) \
                and len(gen.target.elts) == 2 and isinstance(gen.target.elts[1], ast.Name):
            valvar = gen.target.elts[1].id
        elif it.func.attr == "values" and isinstance(gen.target, ast.Name):
            valvar = gen.target.id
        else:
            continue
        for sub in ast.walk(node):
            if isinstance(sub, ast.Subscript) and isinstance(sub.value, ast.Name) \
                    and sub.value.id == valvar and isinstance(sub.slice, ast.Constant) \
                    and isinstance(sub.slice.value, str):
                out.append(ConfigPath(
                    path=tuple(parent) + ("*", sub.slice.value),
                    requiredness="mandatory",
                    source=f"{source_name}:{getattr(node, 'lineno', '?')}",
                ))
    return out


def extract_from_preamble(preamble_src: str, source_name: str,
                          return_aliases: bool = False,
                          ) -> Contract | tuple[Contract, dict[str, tuple[str, ...]]]:
    contract = Contract()
    aliases: dict[str, tuple[str, ...]] = {}
    tree = ast.parse(preamble_src)
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            got = _config_get_call(node)
            if got:
                path, req, default, cond = got
                contract.paths.append(ConfigPath(
                    path=path, requiredness=req, default=default, condition=cond,
                    source=f"{source_name}:{getattr(node, 'lineno', '?')}",
                ))
            continue
        if isinstance(node, (ast.DictComp, ast.ListComp, ast.SetComp, ast.GeneratorExp)):
            contract.paths.extend(_comprehension_wildcards(node, source_name))
            # do not 'continue' — inner config[...] subscripts still register the parent
        if isinstance(node, ast.Assign) and len(node.targets) == 1 \
                and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            val = node.value
            if isinstance(val, ast.Subscript):
                try:
                    p = _literal_subscript_chain(val)
                except _DynamicKey:
                    p = None
                if p:
                    aliases[name] = tuple(p)
            elif isinstance(val, (ast.List, ast.Constant)):
                try:
                    contract.baked_in[name] = ast.literal_eval(val)
                except (ValueError, SyntaxError):
                    pass
        if isinstance(node, ast.Subscript):
            try:
                path = _literal_subscript_chain(node)
            except _DynamicKey:
                contract.incomplete = True
                warn_msg = (
                    f"{source_name}: dynamic config key near line "
                    f"{getattr(node, 'lineno', '?')} — contract may be incomplete"
                )
                if warn_msg not in contract.warnings:
                    contract.warnings.append(warn_msg)
                continue
            if path:
                contract.paths.append(ConfigPath(
                    path=tuple(path), requiredness="mandatory",
                    source=f"{source_name}:{getattr(node, 'lineno', '?')}",
                ))
    # Fix 3: drop spurious optional('fastp',) that the chained .get walk emits as a
    # standalone optional when the real entry is the conditional ('fastp','extra').
    _conds = {p.condition for p in contract.paths
              if p.requiredness == "conditional" and p.condition}
    contract.paths = [p for p in contract.paths
                      if not (p.requiredness == "optional" and p.default is None
                              and p.path in _conds)]
    _dedupe(contract)
    return (contract, aliases) if return_aliases else contract


def _dedupe(contract: Contract) -> None:
    seen: dict[tuple, ConfigPath] = {}
    for p in contract.paths:
        # Prefer mandatory over optional/conditional if both seen for same path.
        key = p.path
        if key not in seen or (p.requiredness == "mandatory"
                               and seen[key].requiredness != "mandatory"):
            seen[key] = p
    contract.paths = list(seen.values())


def build_contract(wrapper_path) -> Contract:
    """Orchestrate preamble extraction + include scanning into a single Contract."""
    wrapper_path = pathlib.Path(wrapper_path)
    text = wrapper_path.read_text()
    preamble = slice_preamble(text)
    contract, aliases = extract_from_preamble(preamble, wrapper_path.name,
                                              return_aliases=True)
    includes, inc_warnings = find_includes(text, wrapper_path.parent)
    contract.warnings.extend(inc_warnings)
    for mod in includes:
        if not mod.exists():
            contract.warnings.append(f"included module not found: {mod}")
            contract.incomplete = True
            continue
        contract.paths.extend(
            extract_alias_subkeys(mod.read_text(), aliases, mod.name))
    _dedupe(contract)
    return contract


_INCLUDE_RE = re.compile(r"""^\s*include:\s*['"]([^'"]+)['"]""")


def find_includes(text: str, wrapper_dir) -> tuple[list[pathlib.Path], list[str]]:
    """Parse include: directives; resolve paths relative to wrapper_dir.
    Indented (conditional) includes are skipped and generate a warning."""
    includes, warnings = [], []
    for line in text.splitlines():
        m = _INCLUDE_RE.match(line)
        if not m:
            continue
        if line[: len(line) - len(line.lstrip())]:  # indented => inside a block
            warnings.append(f"conditional include skipped: {m.group(1)}")
            continue
        includes.append((pathlib.Path(wrapper_dir) / m.group(1)).resolve())
    return includes, warnings


def extract_alias_subkeys(module_text: str, aliases: dict[str, tuple[str, ...]],
                          source_name: str) -> list[ConfigPath]:
    """Find <alias>[...]['literal'] in module text; emit (alias_path + ('*', literal))."""
    out: list[ConfigPath] = []
    for alias, base in aliases.items():
        pat = re.compile(
            r'(?<![A-Za-z0-9_])' + re.escape(alias)
            + r'(?![A-Za-z0-9_])(?:\s*\[[^\]]+\])+\s*\[\s*[\'"]([A-Za-z0-9_\-]+)[\'"]\s*\]'
        )
        for m in pat.finditer(module_text):
            out.append(ConfigPath(
                path=base + ("*", m.group(1)),
                requiredness="mandatory",
                source=source_name,
            ))
    return out


# ---------------------------------------------------------------------------
# Section: List rendering
# ---------------------------------------------------------------------------

def _dotted(path: tuple[str, ...]) -> str:
    return ".".join(path)


def render_list(contract: Contract, fmt: str = "human") -> str:
    """Return a formatted string describing the contract.

    fmt="human": grouped text with MANDATORY / OPTIONAL / CONDITIONAL sections.
    fmt="yaml":  a YAML skeleton with <REQUIRED> / default placeholders.
    """
    if fmt == "yaml":
        return _render_yaml_skeleton(contract)
    lines = []
    groups: dict[str, list[ConfigPath]] = {"mandatory": [], "optional": [], "conditional": []}
    for p in sorted(contract.paths, key=lambda x: x.path):
        groups[p.requiredness].append(p)
    lines.append("MANDATORY:")
    for p in groups["mandatory"]:
        lines.append(f"  {_dotted(p.path)}")
    lines.append("OPTIONAL (default shown):")
    for p in groups["optional"]:
        lines.append(f"  {_dotted(p.path)} = {p.default!r}")
    if groups["conditional"]:
        lines.append("CONDITIONAL (required if parent present):")
        for p in groups["conditional"]:
            lines.append(f"  {_dotted(p.path)} = {p.default!r} "
                         f"(if {_dotted(p.condition)} set)")
    if contract.baked_in:
        lines.append("BAKED-IN (hardcoded in wrapper, not user-settable):")
        for k, v in sorted(contract.baked_in.items()):
            lines.append(f"  {k} = {v!r}")
    if contract.warnings:
        lines.append("WARNINGS:")
        for w in contract.warnings:
            lines.append(f"  ! {w}")
        if contract.incomplete:
            lines.append("  (contract may be INCOMPLETE)")
    return "\n".join(lines) + "\n"


def _render_yaml_skeleton(contract: Contract) -> str:
    import yaml
    tree: dict = {}
    # Sort so shorter (parent) paths are processed before longer (child) paths.
    sorted_paths = sorted(contract.paths, key=lambda x: (len(x.path), x.path))
    for p in sorted_paths:
        cur = tree
        segs = p.path
        for i, seg in enumerate(segs):
            last = i == len(segs) - 1
            key = "<ENTRY>" if seg == "*" else seg
            if last:
                if p.requiredness == "mandatory":
                    # Only set if not already a dict (child paths may have
                    # promoted this key to a mapping already).
                    if not isinstance(cur.get(key), dict):
                        cur[key] = "<REQUIRED>"
                else:
                    if not isinstance(cur.get(key), dict):
                        cur.setdefault(key, p.default)
            else:
                # Promote scalar placeholder to dict if a deeper path needs it.
                if not isinstance(cur.get(key), dict):
                    cur[key] = {}
                cur = cur[key]
    return yaml.safe_dump(tree, default_flow_style=False, sort_keys=True)


# ---------------------------------------------------------------------------
# Section: Path resolution + validation
# ---------------------------------------------------------------------------


@dataclass
class Violation:
    level: str   # "error" | "warning" | "info"
    path: tuple
    message: str


def resolve_config_paths(config_dict: dict) -> None:
    """Expand ~ and $VAR in string values, recursively, in place.

    Mirrors the resolve_config_paths helper in wrapper .smk files.
    """
    for k, v in config_dict.items():
        if isinstance(v, str):
            config_dict[k] = os.path.expandvars(os.path.expanduser(v))
        elif isinstance(v, dict):
            resolve_config_paths(v)
        elif isinstance(v, list):
            config_dict[k] = [os.path.expandvars(os.path.expanduser(i))
                               if isinstance(i, str) else i for i in v]


def _present(cfg, path: tuple[str, ...]) -> bool:
    cur = cfg
    for seg in path:
        if not isinstance(cur, dict) or seg not in cur:
            return False
        cur = cur[seg]
    return True


def validate_config(contract: Contract, cfg: dict) -> list[Violation]:
    """Check cfg against contract; return ALL violations without short-circuiting.

    Wildcard handling: if a parent mapping is absent, report once and skip
    per-entry checks to avoid double-counting.
    """
    violations: list[Violation] = []
    for p in contract.paths:
        if "*" not in p.path:
            if p.requiredness == "mandatory":
                if not _present(cfg, p.path):
                    violations.append(Violation("error", p.path,
                        f"missing required key: {_dotted(p.path)}"))
            elif p.requiredness == "conditional":
                # Fix 1: when parent present and child absent, error if no default.
                if p.condition is not None and _present(cfg, p.condition):
                    if not _present(cfg, p.path):
                        if p.default is None:
                            violations.append(Violation("error", p.path,
                                f"missing required key: {_dotted(p.path)} "
                                f"(required when {_dotted(p.condition)} is set)"))
                        else:
                            violations.append(Violation("info", p.path,
                                f"optional {_dotted(p.path)} absent; default {p.default!r}"))
            else:
                # optional
                if not _present(cfg, p.path):
                    violations.append(Violation("info", p.path,
                        f"optional {_dotted(p.path)} absent; default {p.default!r}"))
            continue
        # wildcard path: split at first '*'
        idx = p.path.index("*")
        parent, leaf = p.path[:idx], p.path[idx + 1:]
        if not _present(cfg, parent):
            violations.append(Violation("error", parent,
                f"missing required mapping: {_dotted(parent)}"))
            continue
        container = cfg
        for seg in parent:
            container = container[seg]
        if not isinstance(container, dict):
            violations.append(Violation("error", parent,
                f"{_dotted(parent)} must be a mapping"))
            continue
        for entry_name, entry_val in container.items():
            full = parent + (entry_name,) + leaf
            if not (isinstance(entry_val, dict) and _present(entry_val, leaf)):
                violations.append(Violation("error", full,
                    f"missing required key: {_dotted(full)}"))
    # unknown top-level keys (downgraded to info when contract is incomplete)
    known_top = {p.path[0] for p in contract.paths} | set(contract.baked_in)
    for k in cfg:
        if k not in known_top:
            violations.append(Violation(
                "info" if contract.incomplete else "warning", (k,),
                f"config key not read by workflow (typo/dead?): {k}"))
    # Fix 2: dedupe by (level, path), preserving first-occurrence order.
    seen_viol: set[tuple] = set()
    deduped: list[Violation] = []
    for v in violations:
        key = (v.level, tuple(v.path))
        if key not in seen_viol:
            seen_viol.add(key)
            deduped.append(v)
    return deduped


# ---------------------------------------------------------------------------
# Section: CLI
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    """Entry point for list/validate subcommands."""
    import argparse
    import yaml

    parser = argparse.ArgumentParser(
        prog="emseq_config_contract",
        description="Derive/validate the minimal config contract of an EM-seq workflow.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    pl = sub.add_parser("list", help="print the minimal required config")
    pl.add_argument("wrapper")
    pl.add_argument("--format", choices=["human", "yaml"], default="human")

    pv = sub.add_parser("validate", help="validate a config YAML against the contract")
    pv.add_argument("wrapper")
    pv.add_argument("config")
    pv.add_argument("--dry-run", action="store_true",
                    help="also run `snakemake -n` as a cross-check (stub)")

    args = parser.parse_args(argv)
    # Fix 4: friendly error on missing wrapper file.
    try:
        contract = build_contract(args.wrapper)
    except FileNotFoundError:
        print(f"error: wrapper not found: {args.wrapper}")
        return 2

    if args.cmd == "list":
        print(render_list(contract, fmt=args.format), end="")
        return 0

    # Fix 4: friendly error on missing config file.
    if not pathlib.Path(args.config).exists():
        print(f"error: config not found: {args.config}")
        return 2
    cfg = yaml.safe_load(open(args.config).read()) or {}
    resolve_config_paths(cfg)
    violations = validate_config(contract, cfg)
    errors = [v for v in violations if v.level == "error"]
    for v in violations:
        prefix = {"error": "ERROR", "warning": "WARN", "info": "info"}[v.level]
        print(f"[{prefix}] {v.message}")
    if contract.incomplete:
        print("NOTE: contract may be incomplete (dynamic keys detected); "
              "consider --dry-run cross-check")
    if args.dry_run:
        print("(--dry-run cross-check not yet wired; run `snakemake -n` manually)")
    print(f"{len(errors)} error(s), "
          f"{sum(1 for v in violations if v.level == 'warning')} warning(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
