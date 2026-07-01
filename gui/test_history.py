"""Tests for the sequence-editor undo/redo stack (gui/history.py)."""

import pytest

from gui.history import UndoStack


def test_empty_stack_cannot_undo_or_redo() -> None:
    s = UndoStack()
    assert not s.can_undo()
    assert not s.can_redo()
    assert s.undo("x") is None
    assert s.redo("x") is None


def test_single_undo_redo_roundtrip() -> None:
    s = UndoStack()
    s.push("ATCG")          # checkpoint the original before editing to ATCGGG
    assert s.can_undo()
    restored = s.undo("ATCGGG")
    assert restored == "ATCG"
    assert not s.can_undo()
    assert s.can_redo()
    reapplied = s.redo("ATCG")
    assert reapplied == "ATCGGG"
    assert s.can_undo()
    assert not s.can_redo()


def test_multi_step_linear_history() -> None:
    s = UndoStack()
    states = ["A", "AB", "ABC"]
    for st in states:
        s.push(st)
    cur = "ABCD"
    assert s.undo(cur) == "ABC"
    assert s.undo("ABC") == "AB"
    assert s.undo("AB") == "A"
    assert not s.can_undo()
    # Redo back up
    assert s.redo("A") == "AB"
    assert s.redo("AB") == "ABC"
    assert s.redo("ABC") == "ABCD"
    assert not s.can_redo()


def test_new_edit_clears_redo() -> None:
    s = UndoStack()
    s.push("A")
    s.undo("AB")            # now redo has "AB"
    assert s.can_redo()
    s.push("A")             # a fresh edit forks the timeline
    assert not s.can_redo()


def test_limit_drops_oldest() -> None:
    s = UndoStack(limit=3)
    for st in ["s1", "s2", "s3", "s4", "s5"]:
        s.push(st)
    assert s.depth == 3
    # Only the three most recent checkpoints survive.
    assert s.undo("cur") == "s5"
    assert s.undo("s5") == "s4"
    assert s.undo("s4") == "s3"
    assert not s.can_undo()


def test_invalid_limit_rejected() -> None:
    with pytest.raises(ValueError):
        UndoStack(limit=0)


def test_clear_resets_both_stacks() -> None:
    s = UndoStack()
    s.push("A")
    s.undo("B")
    assert s.can_redo()
    s.clear()
    assert not s.can_undo()
    assert not s.can_redo()
