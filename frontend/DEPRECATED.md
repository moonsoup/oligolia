# This directory is dead code

**Superseded by `gui/` (PyQt6). Last touched 2026-06-20, before v0.3.0.**

Oligolia was briefly a Tauri + React + Vite app talking to a FastAPI server over
HTTP. It is now a PyQt6 desktop app that imports the FastAPI router functions and
calls them **in-process** — no server, no socket, no port. `src-tauri/` is gone;
this directory is what is left of the front end that went with it.

Nothing in the shipped app builds, imports, bundles or serves anything here. The
PyInstaller specs do not reference it and `oligolia.py` does not know it exists.

## Do not port from it by inspection

If you want behaviour that existed here, check it against the current code and a
real run before reproducing it. The React panels were written against a different
API shape, a different set of endpoints and a different threading model, and
several of the endpoints they call have since changed meaning — the CRISPR guide
convention, the primer Tm calculation and the MSA fallback all changed on
2026-09-17 alone. Copying a component from here and assuming it still describes
the product is how a corrected defect comes back.

## Why it is still in git

Deleting 26 tracked files is the owner's call, not an agent's, and the history is
worth keeping. It is marked rather than removed (#71). If it is retired properly
later, the project convention is to move it out rather than delete it.
