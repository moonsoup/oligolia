"""The one conversion between a Biopython location and the model (#94).

A location is not two integers. Besides its coordinates it carries the `<` / `>`
position class of each boundary, the `join` vs `order` operator over its parts,
and — on a part like `J00194.1:100..202` — the accession of the entry the bases
actually live in. The reader used to flatten all of that to `(int, int)` pairs
and the writer used to rebuild a plain `FeatureLocation`, so none of it could be
recovered at export (#94, the unfixed remainder of #57).

Both directions live here, as one pair, so there is a single place where a
location becomes a model and a single place where it becomes a location again.
Every INSDC reader (`read_genbank`, `read_embl`) and the writer go through them.
"""

from __future__ import annotations

from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    FeatureLocation,
)

from ..models.sequence import Annotation, LocationPart

#: The two position classes that carry INSDC partiality. Both subclass `int`,
#: which is why `int(p.start)` silently erased them.
_CLASS_NAMES = {BeforePosition: "before", AfterPosition: "after"}
_CLASS_TYPES = {"before": BeforePosition, "after": AfterPosition}


def _class_of(position) -> str:
    return _CLASS_NAMES.get(type(position), "exact")


def _position(value: int, class_name: str):
    """A boundary of the given class — an `int` (exact) unless `<` or `>`."""
    factory = _CLASS_TYPES.get(class_name)
    return factory(value) if factory else value


def read_location(location) -> tuple[list[tuple[int, int]], list[LocationPart], str]:
    """`(parts, part_details, operator)` for a parsed Biopython location.

    Parts are in the location's own 5'-to-3' order, as before. `part_details`
    is index-aligned with them; `operator` is `"join"` for a simple location,
    since that is what a bare list of parts has always meant here.
    """
    parts: list[tuple[int, int]] = []
    details: list[LocationPart] = []
    for part in location.parts:
        parts.append((int(part.start), int(part.end)))
        details.append(LocationPart(
            start_class=_class_of(part.start),
            end_class=_class_of(part.end),
            ref=getattr(part, "ref", None) or None,
        ))
    operator = getattr(location, "operator", "") or "join"
    return parts, details, operator


def write_location(ann: Annotation, strand_val: int):
    """The Biopython location an annotation denotes, partiality and all.

    More than one part means a compound location under the annotation's own
    operator — `order()` stays `order()`, and `join()` (the default) stays
    `join()`.
    """
    parts = list(ann.parts) or [(ann.start, ann.end)]
    locs = [
        FeatureLocation(
            _position(start, detail.start_class),
            _position(end, detail.end_class),
            strand=strand_val,
            ref=detail.ref or None,
        )
        for (start, end), detail in zip(parts, ann.details_per_part())
    ]
    if len(locs) == 1:
        return locs[0]
    return CompoundLocation(locs, operator=ann.location_operator or "join")
