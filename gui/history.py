"""Linear undo/redo history for the sequence editor.

Pure Python (no Qt) so it can be unit-tested independently of the GUI. Holds
snapshots of an editable value (the sequence string). A single ``UndoStack``
is kept per sequence so switching sequences never crosses histories.
"""

from __future__ import annotations

DEFAULT_LIMIT = 50


class UndoStack:
    """A bounded linear undo/redo stack of value snapshots.

    Usage: call :meth:`push` with the *current* value immediately before
    replacing it with an edited value. :meth:`undo` returns the value to
    restore (and records the value you pass so :meth:`redo` can reinstate it).
    """

    def __init__(self, limit: int = DEFAULT_LIMIT) -> None:
        if limit < 1:
            raise ValueError("limit must be >= 1")
        self._limit = limit
        self._undo: list[str] = []
        self._redo: list[str] = []

    def push(self, state: str) -> None:
        """Record ``state`` as an undoable checkpoint before a new edit.

        Pushing invalidates any redo history (a new edit forks the timeline)
        and drops the oldest checkpoint once the limit is exceeded.
        """
        self._undo.append(state)
        if len(self._undo) > self._limit:
            self._undo.pop(0)
        self._redo.clear()

    def can_undo(self) -> bool:
        return bool(self._undo)

    def can_redo(self) -> bool:
        return bool(self._redo)

    def undo(self, current: str) -> str | None:
        """Return the previous checkpoint, or ``None`` if nothing to undo.

        ``current`` (the value being undone away from) is pushed onto the redo
        stack so :meth:`redo` can restore it.
        """
        if not self._undo:
            return None
        self._redo.append(current)
        return self._undo.pop()

    def redo(self, current: str) -> str | None:
        """Return the next checkpoint, or ``None`` if nothing to redo."""
        if not self._redo:
            return None
        self._undo.append(current)
        return self._redo.pop()

    def clear(self) -> None:
        self._undo.clear()
        self._redo.clear()

    @property
    def depth(self) -> int:
        """Number of available undo steps (for tests / diagnostics)."""
        return len(self._undo)
