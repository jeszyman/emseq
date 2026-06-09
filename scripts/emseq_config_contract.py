"""Derive and validate the minimal config contract for an EM-seq Snakemake workflow.

Static analysis of a wrapper .smk: parse the plain-Python preamble (before the
first rule) for config[...] / config.get(...) accesses and alias assignments,
follow include: directives, and scan module bodies for alias subscripts to
recover nested per-entry schema. See
docs/superpowers/specs/2026-06-09-snakemake-config-contract-design.md.
"""
from __future__ import annotations

import ast
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


def _literal_subscript_chain(node: "ast.Subscript") -> Optional[list[str]]:
    """Walk a chained subscript and return the key path as a list of strings,
    or None if the chain does not bottom out at the `config` name. Raises
    _DynamicKey on any non-literal (non-string-constant) key."""
    segments: list[str] = []
    cur: Any = node
    while isinstance(cur, ast.Subscript):
        key = cur.slice
        if isinstance(key, ast.Constant) and isinstance(key.value, str):
            segments.append(key.value)
        else:
            raise _DynamicKey()
        cur = cur.value
    if isinstance(cur, ast.Name) and cur.id == "config":
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
