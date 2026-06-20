"""
Post the weekly digest as a comment on a pinned 'Weekly Digest' tracking issue.
Creates the tracking issue if it doesn't exist.
"""

from __future__ import annotations
import os
import datetime
import httpx

OWNER = "moonsoup"
REPO  = "oligolia"
GH_API = "https://api.github.com"
TRACKING_LABEL = "digest"
TODAY = datetime.date.today().isoformat()

GH_HEADERS = {
    "Authorization": f"Bearer {os.environ['GITHUB_TOKEN']}",
    "Accept": "application/vnd.github+json",
    "X-GitHub-Api-Version": "2022-11-28",
}


def ensure_label() -> None:
    with httpx.Client() as c:
        r = c.get(f"{GH_API}/repos/{OWNER}/{REPO}/labels/{TRACKING_LABEL}", headers=GH_HEADERS)
        if r.status_code == 404:
            c.post(f"{GH_API}/repos/{OWNER}/{REPO}/labels", headers=GH_HEADERS,
                   json={"name": TRACKING_LABEL, "color": "4ade80",
                         "description": "Weekly issue digest"})


def find_tracking_issue() -> int | None:
    with httpx.Client() as c:
        r = c.get(f"{GH_API}/repos/{OWNER}/{REPO}/issues",
                  headers=GH_HEADERS,
                  params={"labels": TRACKING_LABEL, "state": "open", "per_page": 1})
        issues = r.json()
        if issues:
            return issues[0]["number"]
    return None


def create_tracking_issue() -> int:
    with httpx.Client() as c:
        r = c.post(f"{GH_API}/repos/{OWNER}/{REPO}/issues", headers=GH_HEADERS,
                   json={"title": "📋 Weekly Issue Digest (auto-updated)",
                         "body": "This issue is auto-updated every Monday with a digest of all open issues "
                                 "and AI-suggested actions. Reply here with your decisions or email back.",
                         "labels": [TRACKING_LABEL]})
        return r.json()["number"]


def post_comment(issue_number: int) -> None:
    try:
        with open("digest_output.html"):
            pass  # just verify it exists
        body = (
            f"## Digest — {TODAY}\n\n"
            "<details><summary>View full digest</summary>\n\n"
            "```\nSee email for formatted version\n```\n\n</details>\n\n"
            "📧 Formatted digest sent to moonsoup@gmail.com"
        )
    except FileNotFoundError:
        body = f"## Digest — {TODAY}\n\n(Digest file not found — check workflow logs)"

    with httpx.Client() as c:
        c.post(f"{GH_API}/repos/{OWNER}/{REPO}/issues/{issue_number}/comments",
               headers=GH_HEADERS, json={"body": body})
    print(f"Posted comment on issue #{issue_number}")


def main() -> None:
    ensure_label()
    issue_number = find_tracking_issue() or create_tracking_issue()
    post_comment(issue_number)
    print(f"Digest posted to #{issue_number}")


if __name__ == "__main__":
    main()
