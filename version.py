"""Single source of truth for the application version."""

VERSION = "0.3.1"
APP_NAME = "Oligolia"
GITHUB_OWNER = "moonsoup"
GITHUB_REPO = "oligolia"

# URL checked on startup for update availability
# Uses GitHub Releases API — works once the repo is public or releases are published
RELEASES_API_URL = f"https://api.github.com/repos/{GITHUB_OWNER}/{GITHUB_REPO}/releases/latest"
RELEASES_PAGE_URL = f"https://github.com/{GITHUB_OWNER}/{GITHUB_REPO}/releases"
