"""Every panel method the main window calls must exist.

#60 was one instance of a class: `gui/main_window.py` called
`self._seq_panel._refresh_active()`, which exists nowhere. Python only finds that at
the moment the menu item is clicked — and in that case the click had already mutated
the user's sequence, so the failure surfaced as "Optimization failed" *after* the
damage, with no undo entry.

Nothing in the suite could catch it, because the only cost of a wrong attribute name
is a crash in a slot nobody drives. This test drives none of them; it reads the call
sites and checks the names resolve, which is enough for this defect class and costs
no Qt event loop.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

GUI = Path(__file__).resolve().parent

# attribute -> the class that attribute holds, for the panels the main window drives.
PANEL_ATTRS = {
    "_seq_panel": ("gui.panels.sequence_panel", "SequencePanel"),
}


def _called_attrs(source: str, obj_attr: str) -> set[str]:
    """Names invoked as `self.<obj_attr>.<name>(...)` anywhere in the source."""
    tree = ast.parse(source)
    found: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        # self._seq_panel.NAME(...)
        if (
            isinstance(fn, ast.Attribute)
            and isinstance(fn.value, ast.Attribute)
            and fn.value.attr == obj_attr
            and isinstance(fn.value.value, ast.Name)
            and fn.value.value.id == "self"
        ):
            found.add(fn.attr)
    return found


def test_every_panel_method_the_main_window_calls_exists() -> None:
    source = (GUI / "main_window.py").read_text()
    missing: list[str] = []

    for obj_attr, (module_name, class_name) in PANEL_ATTRS.items():
        called = _called_attrs(source, obj_attr)
        assert called, f"no calls found on self.{obj_attr} — has main_window.py changed?"

        module = __import__(module_name, fromlist=[class_name])
        cls = getattr(module, class_name)

        for name in sorted(called):
            if not hasattr(cls, name):
                missing.append(f"main_window.py calls self.{obj_attr}.{name}(), absent from {class_name}")

    assert not missing, "\n".join(missing)


def test_the_refresh_active_typo_specifically_stays_gone() -> None:
    """The exact name from #60, pinned so a revert is loud.

    Checked against called attribute names rather than raw text, so the comment in
    main_window.py that explains the history does not trip it.
    """
    source = (GUI / "main_window.py").read_text()
    tree = ast.parse(source)
    called = {
        node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    }
    assert "_refresh_active" not in called, (
        "_refresh_active() does not exist on SequencePanel; use _commit_edit(new_seq, msg), "
        "which pushes an undo checkpoint and re-renders (#60)"
    )


def test_codon_optimization_commits_through_the_undo_path() -> None:
    """#60's second half: the mutation must go through a path that records undo.

    Assigning `_active.seq` directly skips the undo stack, so a user who optimised
    codons could not get back. `_commit_edit` is the only method that pushes onto
    the per-sequence UndoStack before mutating.
    """
    source = (GUI / "main_window.py").read_text()
    optimize = re.search(r"def _optimize_codons\b.*?(?=\n    def |\Z)", source, re.S)
    assert optimize, "could not find _optimize_codons in main_window.py"
    body = optimize.group(0)

    assert "_commit_edit" in body, "codon optimization must commit through _commit_edit (#60)"
    assert not re.search(r"_active\.seq\s*=", body), (
        "codon optimization assigns _active.seq directly, which bypasses the undo stack (#60)"
    )
