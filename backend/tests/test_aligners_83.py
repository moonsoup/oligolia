"""#83 — `GET /alignment/aligners`, the detection the MSA panel asks up front.

Before this, a user learned that multiple alignment cannot run only *after*
pasting sequences and pressing Run, when `/alignment/multiple` returned its 503.
The endpoint here exposes the same decision the run path makes, from the same
function, so the panel and the 503 cannot disagree.

Availability is always controlled here — `shutil.which` patched to find nothing,
or a fake executable written into a tmp dir that is put on PATH. Nothing in this
file depends on whether MUSCLE is really installed on the machine running it,
which is the only way the "available" and "unavailable" halves can both run
everywhere.
"""

import os
import shutil
import stat

import pytest
from fastapi.testclient import TestClient

from backend.routers.alignment import aligner_statuses

HEMOGLOBINS = [
    {"id": "HBB_human", "seq": "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCC"},
    {"id": "HBA1_human", "seq": "ATGGTGCTGTCTCCTGCCGACAAGACCAACGTC"},
]


@pytest.fixture()
def no_aligners(monkeypatch: pytest.MonkeyPatch) -> None:
    """`shutil.which` finds nothing at all."""
    monkeypatch.setattr(shutil, "which", lambda *a, **k: None)


def _fake_executable(directory, name: str) -> str:
    """A real, executable file `name` in `directory`; returns its path."""
    path = os.path.join(str(directory), name)
    with open(path, "w") as fh:
        fh.write("#!/bin/sh\nexit 0\n")
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


@pytest.fixture()
def fake_muscle(tmp_path, monkeypatch: pytest.MonkeyPatch) -> str:
    """A fake `muscle` first on PATH, and nothing else found anywhere."""
    path = _fake_executable(tmp_path, "muscle")
    monkeypatch.setenv("PATH", str(tmp_path))
    return path


def test_aligners_endpoint_reports_both_unavailable_with_install_hints(
    client: TestClient, no_aligners: None
) -> None:
    """Acceptance 1 (endpoint half): nothing installed, both reported, both hinted."""
    r = client.get("/alignment/aligners")
    assert r.status_code == 200, r.text
    data = r.json()
    assert data["any_available"] is False

    by_name = {a["name"]: a for a in data["aligners"]}
    assert set(by_name) == {"muscle", "clustalw"}
    for info in by_name.values():
        assert info["available"] is False
        assert info["path"] is None
        assert info["hint"], "an unavailable aligner must say what to install"
        assert info["url"].startswith("http")

    assert "brew install muscle" in by_name["muscle"]["hint"]
    assert "clustal" in by_name["clustalw"]["hint"].lower()


def test_aligners_endpoint_reports_a_fake_muscle_on_path_as_available(
    client: TestClient, fake_muscle: str
) -> None:
    """Acceptance 2 (endpoint half): the resolved path is the one on PATH."""
    r = client.get("/alignment/aligners")
    assert r.status_code == 200, r.text
    data = r.json()
    assert data["any_available"] is True

    by_name = {a["name"]: a for a in data["aligners"]}
    assert by_name["muscle"]["available"] is True
    assert by_name["muscle"]["path"] == fake_muscle
    # Only muscle was planted; clustalw is still missing and still says so.
    assert by_name["clustalw"]["available"] is False
    assert by_name["clustalw"]["path"] is None


def test_detection_function_is_what_the_endpoint_returns(
    client: TestClient, no_aligners: None
) -> None:
    """One detection, shared: the endpoint is a view of `aligner_statuses()`.

    The GUI calls the function; HTTP callers and the QA pipeline call the
    endpoint. If these two ever diverge, the panel and the 503 can disagree,
    which is the whole point of #83.
    """
    direct = [a.model_dump() for a in aligner_statuses()]
    over_http = client.get("/alignment/aligners").json()["aligners"]
    assert over_http == direct


def test_msa_503_hint_is_the_same_string_the_endpoint_reports(
    client: TestClient, no_aligners: None
) -> None:
    """The 503 and the up-front notice quote the same install hint.

    Not merely "both mention muscle": the hint text itself is shared, so editing
    the hint table cannot leave the panel advising one thing and the error
    another.
    """
    hint = {a["name"]: a["hint"] for a in client.get("/alignment/aligners").json()["aligners"]}

    r = client.post("/alignment/multiple", json={
        "sequences": HEMOGLOBINS, "algorithm": "muscle",
    })
    assert r.status_code == 503, r.text
    assert hint["muscle"] in r.json()["detail"]


def test_msa_503_behaviour_is_unchanged_when_which_finds_nothing(
    client: TestClient, no_aligners: None
) -> None:
    """The existing contract of `/alignment/multiple`, pinned from this side too.

    `backend/tests/test_alignment.py` asserts this only when no aligner happens
    to be installed on the host. With `shutil.which` patched it holds
    everywhere — including on a machine that does have MUSCLE.
    """
    r = client.post("/alignment/multiple", json={
        "sequences": HEMOGLOBINS, "algorithm": "muscle",
    })
    assert r.status_code == 503, r.text
    detail = r.json()["detail"]
    assert "muscle" in detail.lower()
    assert "not installed" in detail.lower()
    assert "install" in detail.lower()
    assert "#58" in detail, "the refusal still explains why nothing is approximated"


def test_msa_still_rejects_an_unknown_algorithm_before_detection(
    client: TestClient, no_aligners: None
) -> None:
    """A bad algorithm name is still a 400, not the new 503 path."""
    r = client.post("/alignment/multiple", json={
        "sequences": HEMOGLOBINS, "algorithm": "not-an-aligner",
    })
    assert r.status_code == 400, r.text
    assert "not-an-aligner" in r.json()["detail"]


def test_msa_too_few_sequences_still_beats_detection(
    client: TestClient, no_aligners: None
) -> None:
    """400 for one sequence, even with no aligner installed — order unchanged."""
    r = client.post("/alignment/multiple", json={
        "sequences": HEMOGLOBINS[:1], "algorithm": "muscle",
    })
    assert r.status_code == 400, r.text
