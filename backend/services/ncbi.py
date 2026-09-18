"""NCBI Entrez E-utilities + Datasets v2 REST API client."""

import time
import io
import json
import re
import zipfile
import httpx
from typing import Any

ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
# Rate limit: 3 req/s without API key, 10 req/s with key
_last_request: float = 0.0


def _throttle() -> None:
    global _last_request
    elapsed = time.monotonic() - _last_request
    if elapsed < 0.34:
        time.sleep(0.34 - elapsed)
    _last_request = time.monotonic()


class BlastFailed(RuntimeError):
    """A BLAST search that cannot produce a result, with the reason."""


def parse_blast_rid(text: str) -> str | None:
    """The RID from a PUT response's QBlastInfo block.

    The old parser required the literal prefix `"    RID = "` — exactly four
    spaces — so any change in NCBI's spacing silently produced "no RID".
    """
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("RID"):
            _, _, value = stripped.partition("=")
            rid = value.strip()
            if rid:
                return rid
    return None


def parse_blast_status(text: str) -> str:
    """WAITING | READY | FAILED | UNKNOWN, from the QBlastInfo block only.

    The old code tested `"Status=READY" in status_r.text` against the whole body.
    That is wrong twice over: the body for a JSON2 request is zip bytes, in which
    the substring never appears, and for a text body any hit description
    containing the phrase would satisfy it.

    Anything not positively understood is UNKNOWN, never READY — defaulting the
    other way is how a zipped payload looked like a finished search (#68).
    """
    inside = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped == "QBlastInfoBegin":
            inside = True
            continue
        if stripped == "QBlastInfoEnd":
            break
        if inside and stripped.startswith("Status"):
            _, _, value = stripped.partition("=")
            status = value.strip().upper()
            if status in ("WAITING", "READY", "FAILED", "UNKNOWN"):
                return status
    return "UNKNOWN"


def unpack_blast_json2(content: bytes) -> dict:
    """The BLAST report out of a JSON2 payload.

    NCBI returns JSON2 as a **zip archive** — Biopython asserts the `PK\x03\x04`
    header for exactly this case (Bio/Blast/__init__.py:1265). The archive holds a
    top-level `<RID>.json` index whose `BlastJSON` entries point at numbered
    `<RID>_1.json` reports.

    Plain JSON is still accepted, in case NCBI ever serves it unzipped. Anything
    that is neither raises rather than being guessed at — an HTML error page used
    to be scanned for a `{` and parsed from there.
    """
    if not content:
        raise BlastFailed("BLAST returned an empty response")

    if content.startswith(b"PK\x03\x04"):
        try:
            with zipfile.ZipFile(io.BytesIO(content)) as zf:
                names = [n for n in zf.namelist() if n.lower().endswith(".json")]
                if not names:
                    raise BlastFailed(
                        f"BLAST returned a zip with no JSON member: {zf.namelist()}"
                    )

                # Prefer a numbered report; the bare <RID>.json is just an index.
                reports = sorted(n for n in names if re.search(r"_\d+\.json$", n))
                for name in reports or sorted(names):
                    payload = json.loads(zf.read(name))
                    if isinstance(payload, dict) and "BlastJSON" in payload:
                        continue  # the index, not a report
                    return payload
                raise BlastFailed("BLAST zip contained only an index, no report")
        except zipfile.BadZipFile as e:
            raise BlastFailed(f"BLAST returned a corrupt zip: {e}") from None

    try:
        payload = json.loads(content)
    except ValueError:
        head = content[:120].decode("utf-8", "replace")
        raise BlastFailed(
            f"BLAST returned neither a zip nor JSON; response began: {head!r}"
        ) from None
    if not isinstance(payload, dict):
        raise BlastFailed(f"BLAST returned JSON that is not an object: {type(payload).__name__}")
    return payload


class NCBIClient:
    def __init__(self, api_key: str = "", email: str = "") -> None:
        self.api_key = api_key
        self.email = email
        self._params: dict[str, str] = {}
        if api_key:
            self._params["api_key"] = api_key
        if email:
            self._params["email"] = email

    def _get(self, url: str, params: dict) -> httpx.Response:
        _throttle()
        merged = {**self._params, **params}
        with httpx.Client(timeout=30) as client:
            r = client.get(url, params=merged)
            r.raise_for_status()
            return r

    # ── ESearch ──────────────────────────────────────────────────────────────

    def esearch(self, db: str, term: str, retmax: int = 20) -> dict[str, Any]:
        """Return list of UIDs matching term in db."""
        r = self._get(
            f"{ENTREZ_BASE}/esearch.fcgi",
            {"db": db, "term": term, "retmax": retmax, "retmode": "json"},
        )
        return r.json()

    # ── EFetch ───────────────────────────────────────────────────────────────

    def efetch(self, db: str, ids: list[str], rettype: str = "gb", retmode: str = "text") -> str:
        """Fetch full records for given UIDs."""
        r = self._get(
            f"{ENTREZ_BASE}/efetch.fcgi",
            {"db": db, "id": ",".join(ids), "rettype": rettype, "retmode": retmode},
        )
        return r.text

    def efetch_json(self, db: str, ids: list[str]) -> dict[str, Any]:
        """Fetch records in JSON summary format."""
        r = self._get(
            f"{ENTREZ_BASE}/esummary.fcgi",
            {"db": db, "id": ",".join(ids), "retmode": "json"},
        )
        return r.json()

    # ── ELink ────────────────────────────────────────────────────────────────

    def elink(self, dbfrom: str, db: str, ids: list[str]) -> dict[str, Any]:
        r = self._get(
            f"{ENTREZ_BASE}/elink.fcgi",
            {"dbfrom": dbfrom, "db": db, "id": ",".join(ids), "retmode": "json"},
        )
        return r.json()

    # ── Convenience ──────────────────────────────────────────────────────────

    def gene_uid_to_refseq_rna(self, gene_uid: str) -> list[str]:
        """Convert an NCBI Gene UID to RefSeq mRNA nuccore UIDs via elink."""
        result = self.elink("gene", "nuccore", [gene_uid])
        for linkset in result.get("linksets", []):
            for lslink in linkset.get("linksetdbs", []):
                if lslink.get("linkname") == "gene_nuccore_refseqrna":
                    return [str(i) for i in lslink.get("links", [])]
        # Fallback: any nuccore link
        for linkset in result.get("linksets", []):
            for lslink in linkset.get("linksetdbs", []):
                if lslink.get("linkname") == "gene_nuccore":
                    links = [str(i) for i in lslink.get("links", [])]
                    return links[:1]  # just the first to avoid huge genomic records
        return []

    def fetch_fasta_for_gene(self, gene_uid: str) -> str:
        """Fetch RefSeq mRNA FASTA for a Gene UID (handles the Gene UID → nuccore conversion)."""
        nuccore_ids = self.gene_uid_to_refseq_rna(gene_uid)
        if not nuccore_ids:
            raise ValueError(
                f"No RefSeq mRNA records found for Gene UID {gene_uid}. "
                "The gene may not have an mRNA sequence in NCBI."
            )
        return self.efetch("nucleotide", nuccore_ids[:1], rettype="fasta", retmode="text")

    def search_genes(self, query: str, organism: str = "Homo sapiens", max_results: int = 20) -> list[dict]:
        term = f"{query}[Gene Name] AND {organism}[Organism]"
        result = self.esearch("gene", term, retmax=max_results)
        ids = result.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        summary = self.efetch_json("gene", ids)
        records = summary.get("result", {})
        uids = records.get("uids", [])
        return [records[uid] for uid in uids if uid in records]

    def fetch_nucleotide(self, accession: str) -> str:
        """Fetch GenBank flat file for a nucleotide accession."""
        return self.efetch("nucleotide", [accession], rettype="gb", retmode="text")

    def fetch_fasta(self, accession: str, db: str = "nucleotide") -> str:
        return self.efetch(db, [accession], rettype="fasta", retmode="text")

    #: Seconds between status polls, and how many to make. NCBI asks callers not
    #: to poll more often than every 10 s for a search of any size.
    BLAST_POLL_SECONDS = 10
    BLAST_MAX_POLLS = 30

    def blast_search(self, sequence: str, program: str = "blastn", database: str = "nt", max_hits: int = 10) -> dict:
        """Submit a BLAST search and return the JSON2 report.

        Follows NCBI's documented two-step pattern, which the previous version
        conflated (#68):

        1. poll `CMD=Get&FORMAT_OBJECT=SearchInfo`, a small TEXT response whose
           QBlastInfo block carries Status=WAITING|READY|FAILED;
        2. only once READY, fetch `FORMAT_TYPE=JSON2` and unzip it.

        Before, the status check was `"Status=READY" in body` against a body that
        was a zip archive, so the substring never matched and every search burned
        the full poll budget before raising TimeoutError.
        """
        url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
        with httpx.Client(timeout=60) as client:
            r = client.post(url, data={
                "CMD": "Put",
                "PROGRAM": program,
                "DATABASE": database,
                "QUERY": sequence,
                "FORMAT_TYPE": "JSON2",
                "HITLIST_SIZE": max_hits,
            })
            r.raise_for_status()

            rid = parse_blast_rid(r.text)
            if not rid:
                raise BlastFailed(
                    "BLAST did not return an RID; response began: "
                    f"{r.text[:200]!r}"
                )

            for _ in range(self.BLAST_MAX_POLLS):
                time.sleep(self.BLAST_POLL_SECONDS)
                info = client.get(url, params={
                    "CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo",
                })
                status = parse_blast_status(info.text)
                if status == "WAITING":
                    continue
                if status == "FAILED":
                    raise BlastFailed(f"BLAST search {rid} failed server-side")
                if status == "UNKNOWN":
                    raise BlastFailed(
                        f"BLAST search {rid} expired or is not recognised by the server"
                    )
                # READY
                result = client.get(url, params={
                    "CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2",
                })
                result.raise_for_status()
                return unpack_blast_json2(result.content)

            raise TimeoutError(
                f"BLAST search {rid} still running after "
                f"{self.BLAST_POLL_SECONDS * self.BLAST_MAX_POLLS}s"
            )

    def search_clinvar(self, gene: str, max_results: int = 50) -> list[dict]:
        term = f"{gene}[Gene Name]"
        result = self.esearch("clinvar", term, retmax=max_results)
        ids = result.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []
        summary = self.efetch_json("clinvar", ids)
        records = summary.get("result", {})
        uids = records.get("uids", [])
        return [records[uid] for uid in uids if uid in records]
