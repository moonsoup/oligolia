"""Drive the real WorkflowPanel offscreen.

#55's second confirmed crash path: entering `[]` as a step's parameters. `_run`
caught only `json.JSONDecodeError`, so `[]` parsed fine as JSON and then pydantic's
`ValidationError` escaped the slot — which PyQt6 turns into exit -6.

The crash guard now catches anything that gets this far, but a user who typed the
wrong thing should get a message about their parameters, not a generic "unexpected
error" dialog. So the specific handler is tested here, and the guard is tested in
gui/test_crash_guard.py.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402
from PyQt6.QtWidgets import QApplication, QTableWidgetItem  # noqa: E402

from gui.panels import workflow_panel as wp  # noqa: E402


@pytest.fixture(scope="module")
def app() -> QApplication:
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def recorded(monkeypatch) -> list[tuple[str, str]]:
    """Capture QMessageBox calls instead of blocking on a modal dialog."""
    seen: list[tuple[str, str]] = []

    class FakeBox:
        @staticmethod
        def critical(_parent, title, text, *a, **k):
            seen.append((title, text))

        @staticmethod
        def warning(_parent, title, text, *a, **k):
            seen.append((title, text))

        @staticmethod
        def information(_parent, title, text, *a, **k):
            seen.append((title, text))

    monkeypatch.setattr(wp, "QMessageBox", FakeBox)
    return seen


def _panel_with_params(text: str) -> wp.WorkflowPanel:
    panel = wp.WorkflowPanel()
    panel._add_step()
    panel._steps_table.setItem(0, 1, QTableWidgetItem(text))
    return panel


def test_a_json_list_as_parameters_is_reported_not_fatal(app, recorded) -> None:
    """#55: `[]` is valid JSON and an invalid params value — the exact repro."""
    panel = _panel_with_params("[]")

    panel._run()  # must not raise

    assert recorded, "the user was told nothing"
    title, text = recorded[0]
    assert "param" in (title + text).lower(), recorded


def test_malformed_json_is_still_reported(app, recorded) -> None:
    """The case that already worked must keep working."""
    panel = _panel_with_params("{not json")

    panel._run()

    assert recorded, "malformed JSON should still be reported"
    assert "param" in (recorded[0][0] + recorded[0][1]).lower(), recorded


def test_a_json_string_as_parameters_is_reported(app, recorded) -> None:
    panel = _panel_with_params('"just a string"')
    panel._run()
    assert recorded, recorded


def test_a_json_number_as_parameters_is_reported(app, recorded) -> None:
    panel = _panel_with_params("42")
    panel._run()
    assert recorded, recorded


def test_no_steps_is_reported_rather_than_run(app, recorded) -> None:
    panel = wp.WorkflowPanel()
    panel._run()
    assert recorded
    assert "step" in (recorded[0][0] + recorded[0][1]).lower()
