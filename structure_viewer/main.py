"""Oligolia Structure Viewer — standalone companion app.

Kept separate from the core Oligolia installer so the main app never bundles a
Chromium/WebEngine runtime (see gui/panels/structure_panel.py's docstring for
why). Launched by the core app's "Open in 3D Viewer" button
(gui/plugins/structure_viewer_launcher.py) as a one-shot subprocess handoff —
not a live two-way session:

    oligolia-structure-viewer --pdb <path> [--interactions <path.json>]

Renders the structure via a vendored, MIT-licensed copy of 3Dmol.js
(assets/3Dmol-min.js) inlined directly into the page HTML, coloring by
per-residue B-factor/pLDDT confidence and highlighting any residues flagged
as putative interaction points.
"""

import argparse
import json
import sys
from pathlib import Path

from PyQt6.QtWidgets import QApplication, QMainWindow
from PyQt6.QtWebEngineWidgets import QWebEngineView

ASSETS_DIR = Path(__file__).resolve().parent / "assets"

try:
    from version import VERSION
except ImportError:
    VERSION = "0.0.0"

_HTML_TEMPLATE = """<!doctype html>
<html>
<head>
<meta charset="utf-8">
<style>html, body, #viewer { width: 100%; height: 100%; margin: 0; padding: 0; background: #0b1220; }</style>
<script>__THREEDMOL_JS__</script>
</head>
<body>
<div id="viewer"></div>
<script>
  var viewer = $3Dmol.createViewer("viewer", {backgroundColor: "0x0b1220"});
  var pdbText = __PDB_JSON__;
  var points = __POINTS_JSON__;
  viewer.addModel(pdbText, "pdb");
  viewer.setStyle({}, {cartoon: {colorscheme: {prop: "b", gradient: "roygb", min: 50, max: 100}}});
  var interactionResidues = points.filter(function (p) { return p.is_putative_interaction_point; })
                                   .map(function (p) { return p.residue_index; });
  if (interactionResidues.length) {
    viewer.addStyle({resi: interactionResidues}, {sphere: {color: "magenta", radius: 0.6}});
  }
  viewer.zoomTo();
  viewer.render();
</script>
</body>
</html>"""


def build_html(pdb_text: str, points: list[dict]) -> str:
    threedmol_js = (ASSETS_DIR / "3Dmol-min.js").read_text()
    return (
        _HTML_TEMPLATE
        .replace("__THREEDMOL_JS__", threedmol_js)
        .replace("__PDB_JSON__", json.dumps(pdb_text))
        .replace("__POINTS_JSON__", json.dumps(points))
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Oligolia Structure Viewer")
    parser.add_argument("--pdb", required=True, help="Path to a .pdb structure file")
    parser.add_argument("--interactions", help="Path to an interaction-points JSON file (optional)")
    args = parser.parse_args()

    pdb_text = Path(args.pdb).read_text()
    points = json.loads(Path(args.interactions).read_text()) if args.interactions else []

    app = QApplication(sys.argv)
    window = QMainWindow()
    window.setWindowTitle(f"Oligolia Structure Viewer {VERSION}")
    window.resize(1000, 800)

    view = QWebEngineView()
    view.setHtml(build_html(pdb_text, points))
    window.setCentralWidget(view)
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
