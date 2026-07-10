from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal


ClaimClass = Literal[
    "exact_math",
    "diagnostic",
    "conditional_physics",
    "coherence_requirement",
]


@dataclass(slots=True)
class ProbeResult:
    claim_id: str
    claim_class: ClaimClass
    verified_fact: str
    premises: list[str] = field(default_factory=list)
    disputed_premises: list[str] = field(default_factory=list)
    conditional_inference: str | None = None
    non_claims: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)
    caveats: list[str] = field(default_factory=list)
    pinned_expressions: dict[str, Any] = field(default_factory=dict)
    symbolic_verdicts: dict[str, Any] = field(default_factory=dict)
    sources: list[str] = field(default_factory=list)
    raw_data: Any = None
