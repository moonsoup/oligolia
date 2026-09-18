"""#68: BLAST asked for JSON2, which NCBI returns as a zip.

Biopython asserts exactly this at Bio/Blast/__init__.py:1265 —
`assert data.startswith(b"PK\\x03\\x04")  # zipped file` for XML2/JSON2.

So `.json()` failed on the response, `"Status=READY" in status_r.text` never
matched binary zip bytes, and every search burned the full 150 s poll.

The parsing and status logic are pure functions, so this runs with no network:
the fixtures are zips built in the test, in NCBI's actual JSON2 layout (a
top-level `<RID>.json` index pointing at `<RID>_1.json` reports).
"""

from __future__ import annotations

import io
import json
import zipfile

import pytest

from backend.services.ncbi import (
    BlastFailed,
    parse_blast_rid,
    parse_blast_status,
    unpack_blast_json2,
)

REPORT = {
    "BlastOutput2": {
        "report": {
            "program": "blastn",
            "results": {"search": {"query_id": "Query_1", "hits": [
                {"num": 1, "description": [{"title": "Homo sapiens BRCA1"}]},
            ]}},
        }
    }
}


def _zip(members: dict[str, str]) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for name, text in members.items():
            zf.writestr(name, text)
    return buf.getvalue()


# ── the zip really is a zip ──────────────────────────────────────────────────

def test_a_json2_payload_starts_with_the_zip_magic() -> None:
    """The premise, asserted so the fixtures cannot drift from reality."""
    data = _zip({"ABCD1234_1.json": json.dumps(REPORT)})
    assert data.startswith(b"PK\x03\x04")


def test_the_report_is_extracted_from_a_single_member_zip() -> None:
    data = _zip({"ABCD1234_1.json": json.dumps(REPORT)})
    assert unpack_blast_json2(data) == REPORT


def test_the_report_is_extracted_past_the_index_member() -> None:
    """NCBI's real layout: a top-level index pointing at numbered reports."""
    index = {"BlastJSON": [{"File": "ABCD1234_1.json"}]}
    data = _zip({
        "ABCD1234.json": json.dumps(index),
        "ABCD1234_1.json": json.dumps(REPORT),
    })
    assert unpack_blast_json2(data) == REPORT


def test_the_first_report_is_returned_for_a_multi_query_zip() -> None:
    second = {"BlastOutput2": {"report": {"program": "blastn", "results": {}}}}
    index = {"BlastJSON": [{"File": "R_1.json"}, {"File": "R_2.json"}]}
    data = _zip({
        "R.json": json.dumps(index),
        "R_1.json": json.dumps(REPORT),
        "R_2.json": json.dumps(second),
    })
    assert unpack_blast_json2(data) == REPORT


def test_plain_json_is_still_accepted() -> None:
    """If NCBI ever returns unzipped JSON, do not insist on a zip."""
    assert unpack_blast_json2(json.dumps(REPORT).encode()) == REPORT


def test_something_that_is_neither_is_refused_not_guessed() -> None:
    with pytest.raises(BlastFailed):
        unpack_blast_json2(b"<!DOCTYPE html><html>error page</html>")
    with pytest.raises(BlastFailed):
        unpack_blast_json2(b"")


def test_a_zip_with_no_json_member_is_refused() -> None:
    with pytest.raises(BlastFailed):
        unpack_blast_json2(_zip({"readme.txt": "nothing useful"}))


# ── status comes from SearchInfo, not from a substring of the payload ────────

WAITING = """<p><!--
QBlastInfoBegin
\tStatus=WAITING
QBlastInfoEnd
--><p>
"""

READY = """<p><!--
QBlastInfoBegin
\tStatus=READY
\tThereAreHits=yes
QBlastInfoEnd
--><p>
"""

FAILED = """<p><!--
QBlastInfoBegin
\tStatus=FAILED
QBlastInfoEnd
--><p>
"""

UNKNOWN = """<p><!--
QBlastInfoBegin
\tStatus=UNKNOWN
QBlastInfoEnd
--><p>
"""


@pytest.mark.parametrize("text,expected", [
    (WAITING, "WAITING"),
    (READY, "READY"),
    (FAILED, "FAILED"),
    (UNKNOWN, "UNKNOWN"),
])
def test_status_is_parsed_from_the_qblastinfo_block(text: str, expected: str) -> None:
    assert parse_blast_status(text) == expected


def test_an_unparseable_status_is_unknown_not_ready() -> None:
    """Defaulting to READY is how a zip body looked like a finished search."""
    assert parse_blast_status("") == "UNKNOWN"
    assert parse_blast_status("<html>server error</html>") == "UNKNOWN"


def test_status_is_not_matched_inside_unrelated_text() -> None:
    """The old check was `"Status=READY" in text`, which a hit title could satisfy."""
    sneaky = READY.replace("Status=READY", "Status=WAITING") + "\nhit title: Status=READY gene"
    assert parse_blast_status(sneaky) == "WAITING"


# ── the RID ──────────────────────────────────────────────────────────────────

PUT_RESPONSE = """<p><!--
QBlastInfoBegin
\tRID = ABCD1234567
\tRTOE = 27
QBlastInfoEnd
--><p>
"""


def test_the_rid_is_parsed() -> None:
    assert parse_blast_rid(PUT_RESPONSE) == "ABCD1234567"


def test_a_missing_rid_is_none() -> None:
    assert parse_blast_rid("<html>no rid here</html>") is None


def test_an_rid_with_unexpected_spacing_is_still_found() -> None:
    """The old parser required exactly four leading spaces before 'RID ='."""
    assert parse_blast_rid("QBlastInfoBegin\n  RID=XYZ999\nQBlastInfoEnd") == "XYZ999"
    assert parse_blast_rid("\tRID = XYZ999 \n") == "XYZ999"
