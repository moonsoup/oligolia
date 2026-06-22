"""File upload/download/format conversion endpoints."""

from fastapi import APIRouter, HTTPException, UploadFile, File
from fastapi.responses import PlainTextResponse, Response
from ..models.sequence import Sequence
from ..formats import (
    read_fasta, read_fastq, read_genbank, read_embl,
    write_fasta, write_fastq, write_genbank,
    parse_vcf, write_vcf, parse_gff3, parse_gtf,
)
from ..models.variant import Variant

router = APIRouter(prefix="/files", tags=["files"])

FORMAT_READERS = {
    "fasta": read_fasta,
    "fa": read_fasta,
    "fna": read_fasta,
    "ffn": read_fasta,
    "faa": read_fasta,
    "fastq": read_fastq,
    "fq": read_fastq,
    "gb": read_genbank,
    "gbk": read_genbank,
    "genbank": read_genbank,
    "embl": read_embl,
}


def _decode(content: bytes) -> str:
    """Decode bytes trying UTF-8, Latin-1, and Windows-1252 in order."""
    for enc in ("utf-8", "latin-1", "cp1252"):
        try:
            return content.decode(enc)
        except (UnicodeDecodeError, ValueError):
            continue
    raise HTTPException(422, "File encoding not supported — please re-save as UTF-8")


@router.post("/parse/sequence", response_model=list[Sequence])
async def parse_sequence_file(file: UploadFile = File(...)) -> list[Sequence]:
    """Upload any supported sequence file (FASTA, FASTQ, GenBank, EMBL) and get parsed sequences."""
    ext = (file.filename or "").rsplit(".", 1)[-1].lower()
    if ext not in FORMAT_READERS:
        raise HTTPException(415, f"Unsupported format: .{ext}. Supported: {', '.join(FORMAT_READERS)}")
    content = await file.read()
    try:
        return FORMAT_READERS[ext](_decode(content))
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(422, f"Parse error: {e}")


@router.post("/parse/vcf", response_model=list[Variant])
async def parse_vcf_file(file: UploadFile = File(...)) -> list[Variant]:
    content = await file.read()
    try:
        return parse_vcf(_decode(content))
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(422, f"VCF parse error: {e}")


@router.post("/parse/gff", response_model=list[dict])
async def parse_gff_file(file: UploadFile = File(...)) -> list[dict]:
    content = await file.read()
    ext = (file.filename or "gff3").rsplit(".", 1)[-1].lower()
    try:
        if "gtf" in ext:
            annotations = parse_gtf(_decode(content))
        else:
            annotations = parse_gff3(_decode(content))
        return [a.model_dump() for a in annotations]
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(422, f"GFF parse error: {e}")


@router.post("/convert")
async def convert_file(file: UploadFile = File(...), target_format: str = "fasta") -> Response:
    """Convert between sequence formats (fasta, fastq, genbank)."""
    ext = (file.filename or "").rsplit(".", 1)[-1].lower()
    if ext not in FORMAT_READERS:
        raise HTTPException(415, f"Unsupported input format: .{ext}")
    content = await file.read()
    try:
        sequences = FORMAT_READERS[ext](_decode(content))
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(422, f"Parse error: {e}")

    target = target_format.lower()
    writers = {"fasta": write_fasta, "fastq": write_fastq, "genbank": write_genbank, "gb": write_genbank}
    if target not in writers:
        raise HTTPException(400, f"Unsupported target format: {target}. Supported: {', '.join(writers)}")

    output = writers[target](sequences)
    ext_map = {"fasta": "fa", "fastq": "fq", "genbank": "gb", "gb": "gb"}
    return PlainTextResponse(
        output,
        headers={"Content-Disposition": f'attachment; filename="converted.{ext_map.get(target, target)}"'},
    )


@router.post("/download/fasta")
def download_fasta(sequences: list[Sequence]) -> PlainTextResponse:
    content = write_fasta(sequences)
    return PlainTextResponse(
        content, headers={"Content-Disposition": 'attachment; filename="sequences.fa"'}
    )


@router.post("/download/genbank")
def download_genbank(sequences: list[Sequence]) -> PlainTextResponse:
    content = write_genbank(sequences)
    return PlainTextResponse(
        content, headers={"Content-Disposition": 'attachment; filename="sequences.gb"'}
    )


@router.post("/download/vcf")
def download_vcf(variants: list[Variant]) -> PlainTextResponse:
    content = write_vcf(variants)
    return PlainTextResponse(
        content, headers={"Content-Disposition": 'attachment; filename="variants.vcf"'}
    )
