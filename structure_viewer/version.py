"""Independent version for the Oligolia Structure Viewer companion app.

Deliberately not coupled to the core app's version.py / per-push patch-bump
cadence (CLAUDE.md rule 3 is about the core app) — this only bumps when the
viewer itself actually changes, since it triggers its own 3-platform,
Chromium-bundled build.
"""

VERSION = "0.1.0"
APP_NAME = "Oligolia Structure Viewer"
GITHUB_OWNER = "moonsoup"
GITHUB_REPO = "oligolia"
