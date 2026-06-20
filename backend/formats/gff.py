"""GFF3 / GTF genome annotation format parser."""

from __future__ import annotations
from typing import TextIO
from ..models.sequence import Annotation, Strand


def _parse_attributes_gff3(attr_str: str) -> dict:
    result: dict = {}
    for item in attr_str.strip().split(";"):
        item = item.strip()
        if "=" in item:
            k, v = item.split("=", 1)
            result[k.strip()] = v.strip()
    return result


def _parse_attributes_gtf(attr_str: str) -> dict:
    result: dict = {}
    for item in attr_str.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(None, 1)
        if len(parts) == 2:
            k, v = parts
            result[k] = v.strip('"')
    return result


def _strand(s: str) -> Strand:
    if s == "+":
        return Strand.PLUS
    if s == "-":
        return Strand.MINUS
    return Strand.BOTH


def parse_gff3(source: str | TextIO) -> list[Annotation]:
    lines = source.splitlines() if isinstance(source, str) else source
    annotations = []
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        seqid, source_, feature, start, end, score, strand_str, phase = parts[:8]
        attrs_str = parts[8] if len(parts) > 8 else ""
        attrs = _parse_attributes_gff3(attrs_str)
        attrs.update({"seqid": seqid, "source": source_, "score": score, "phase": phase})
        annotations.append(Annotation(
            feature_type=feature,
            start=int(start) - 1,  # convert 1-based to 0-based
            end=int(end),
            strand=_strand(strand_str),
            qualifiers=attrs,
        ))
    return annotations


def parse_gtf(source: str | TextIO) -> list[Annotation]:
    lines = source.splitlines() if isinstance(source, str) else source
    annotations = []
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        seqid, source_, feature, start, end, score, strand_str, frame = parts[:8]
        attrs_str = parts[8] if len(parts) > 8 else ""
        attrs = _parse_attributes_gtf(attrs_str)
        attrs.update({"seqid": seqid, "source": source_, "score": score, "frame": frame})
        annotations.append(Annotation(
            feature_type=feature,
            start=int(start) - 1,
            end=int(end),
            strand=_strand(strand_str),
            qualifiers=attrs,
        ))
    return annotations


def write_gff3(annotations: list[Annotation], seqid: str = "unknown") -> str:
    lines = ["##gff-version 3"]
    for ann in annotations:
        attrs = ";".join(f"{k}={v}" for k, v in ann.qualifiers.items()
                         if k not in ("seqid", "source", "score", "phase"))
        seq = ann.qualifiers.get("seqid", seqid)
        src = ann.qualifiers.get("source", "oligolia")
        score = ann.qualifiers.get("score", ".")
        phase = ann.qualifiers.get("phase", ".")
        strand = ann.strand.value
        lines.append("\t".join([
            seq, src, ann.feature_type,
            str(ann.start + 1), str(ann.end),
            score, strand, phase, attrs or ".",
        ]))
    return "\n".join(lines) + "\n"
