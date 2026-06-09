"""Derive and validate the minimal config contract for an EM-seq Snakemake workflow.

Static analysis of a wrapper .smk: parse the plain-Python preamble (before the
first rule) for config[...] / config.get(...) accesses and alias assignments,
follow include: directives, and scan module bodies for alias subscripts to
recover nested per-entry schema. See
docs/superpowers/specs/2026-06-09-snakemake-config-contract-design.md.
"""
from __future__ import annotations

import ast
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


def _literal_subscript_chain(node: "ast.Subscript") -> Optional[list[str]]:
    """If `node` is config[...] (optionally chained, all-literal), return the
    key path as a list of strings. Return None if it does not bottom out at the
    `config` name, and [] sentinel handling is via the caller. A non-literal key
    raises _DynamicKey."""
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


class _DynamicKey(Exception):
    pass


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


def extract_from_preamble(preamble_src: str, source_name: str) -> Contract:
    contract = Contract()
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
        if isinstance(node, ast.Subscript):
            try:
                path = _literal_subscript_chain(node)
            except _DynamicKey:
                contract.incomplete = True
                contract.warnings.append(
                    f"{source_name}: dynamic config key near line "
                    f"{getattr(node, 'lineno', '?')} — contract may be incomplete"
                )
                continue
            if path:
                contract.paths.append(ConfigPath(
                    path=tuple(path), requiredness="mandatory",
                    source=f"{source_name}:{getattr(node, 'lineno', '?')}",
                ))
    _dedupe(contract)
    return contract


def _dedupe(contract: Contract) -> None:
    seen: dict[tuple, ConfigPath] = {}
    for p in contract.paths:
        # Prefer mandatory over optional/conditional if both seen for same path.
        key = p.path
        if key not in seen or (p.requiredness == "mandatory"
                               and seen[key].requiredness != "mandatory"):
            seen[key] = p
    contract.paths = list(seen.values())


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
