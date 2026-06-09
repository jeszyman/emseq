"""Derive and validate the minimal config contract for an EM-seq Snakemake workflow.

Static analysis of a wrapper .smk: parse the plain-Python preamble (before the
first rule) for config[...] / config.get(...) accesses and alias assignments,
follow include: directives, and scan module bodies for alias subscripts to
recover nested per-entry schema. See
docs/superpowers/specs/2026-06-09-snakemake-config-contract-design.md.
"""
from __future__ import annotations

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
