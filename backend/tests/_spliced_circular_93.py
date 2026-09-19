"""A tiny circular record whose only interesting features are compound (#93).

Reverse Complement mapped every interval correctly but then reversed the part
list, so a spliced CDS read its exons back to front and encoded a different
protein. Only a feature with a `join()` can show that, and the two shapes that
matter are both here:

* `CDS join(11..25,51..86)` — a spliced ORF. The 51 bases it splices are a real
  reading frame, `ATG…TAA`, translating to `MARKEPGIHRLQALGT*`, so "the exons
  came back in the wrong order" is visible as "the protein changed" rather than
  as a coordinate diff.
* `rep_origin join(115..120,1..6)` — origin-spanning on a circular molecule, the
  one case where a descending join *is* legitimate, so the fix must not simply
  sort parts by coordinate.

Shared by `backend/tests/test_reverse_complement_part_order_93.py` and
`gui/panels/test_sequence_panel_reverse_complement_93.py` so the backend
arithmetic and the widget the user clicks are checked against the same bytes.
"""

from __future__ import annotations

#: The protein `CDS join(11..25,51..86)` encodes on the unflipped record, and
#: therefore the protein it must still encode on the flipped one.
SPLICED_PROTEIN = "MARKEPGIHRLQALGT*"

#: The bases that CDS splices together, in reading order.
SPLICED_CDS = "ATGGCACGTAAAGAACCTGGAATTCATCGACTGCAAGCATTGGGAACGTAA"

#: The bases `rep_origin join(115..120,1..6)` reads, in reading order.
ORIGIN_SPANNING = "TGCATGACGTTG"

SPLICED_CIRCULAR = """LOCUS       SPLICE93                 120 bp    DNA     circular SYN 01-JAN-2024
DEFINITION  a spliced CDS and an origin-spanning feature, for #93.
ACCESSION   SPLICE93
FEATURES             Location/Qualifiers
     source          1..120
     rep_origin      join(115..120,1..6)
                     /label="wraps"
     CDS             join(11..25,51..86)
                     /gene="demo"
                     /label="spliced"
ORIGIN
        1 acgttgcatg atggcacgta aagaagcatg acgttgcatg acgttgcatg cctggaattc
       61 atcgactgca agcattggga acgtaacatg acgttgcatg acgttgcatg acgttgcatg
//
"""
