"""
Weekly issue digest generator.

1. Fetches all open GitHub issues via API
2. Scans for injection/spam/phishing — silently drops and closes bad ones
3. Uses Claude API to analyse each legitimate issue and suggest action
4. Writes digest_output.html (emailed to owner)
5. Sets GitHub Actions output: date, issue_count

Expects env vars:
  GITHUB_TOKEN       — from GH Actions secrets
  ANTHROPIC_API_KEY  — from GH Actions secrets
"""

from __future__ import annotations
import os
import json
import re
import datetime
import httpx
import anthropic

OWNER = "moonsoup"
REPO = "oligolia"
OWNER_EMAIL = "moonsoup@gmail.com"
GH_API = "https://api.github.com"
TODAY = datetime.date.today().isoformat()

GITHUB_TOKEN = os.environ["GITHUB_TOKEN"]
ANTHROPIC_KEY = os.environ.get("ANTHROPIC_API_KEY", "")

GH_HEADERS = {
    "Authorization": f"Bearer {GITHUB_TOKEN}",
    "Accept": "application/vnd.github+json",
    "X-GitHub-Api-Version": "2022-11-28",
}


# ── Security filter ────────────────────────────────────────────────────────────

SPAM_PATTERNS = [
    r"(earn|make).{0,20}\$[\d,]+.{0,20}(day|week|month|hour)",
    r"(click here|visit now|limited offer|act now).{0,30}http",
    r"(phish|cred(ential)?s?|passw(or)?d|log\s?in).{0,30}(steal|harvest|grab)",
    r"(bitcoin|crypto|nft).{0,40}(invest|profit|earn)",
    r"<script|javascript:|on(load|click|error)\s*=",
    r"(union\s+select|drop\s+table|insert\s+into|exec\s*\()",
    r"(xss|csrf|sql\s+inject)",
    r"(sex|porn|xxx|nude|naked).{0,20}(video|photo|link|click)",
    r"\b(hack|exploit|bypass|crack).{0,20}(password|account|system)",
]
_SPAM_RE = re.compile("|".join(SPAM_PATTERNS), re.IGNORECASE)


def is_spam(issue: dict) -> bool:
    text = f"{issue.get('title','')} {issue.get('body','')}"
    return bool(_SPAM_RE.search(text))


def close_spam_issue(issue_number: int, reason: str) -> None:
    with httpx.Client() as c:
        c.patch(
            f"{GH_API}/repos/{OWNER}/{REPO}/issues/{issue_number}",
            headers=GH_HEADERS,
            json={"state": "closed", "state_reason": "not_planned"},
        )
        c.post(
            f"{GH_API}/repos/{OWNER}/{REPO}/issues/{issue_number}/comments",
            headers=GH_HEADERS,
            json={"body": "🚫 Closed automatically: flagged as spam/injection/off-topic."},
        )
    print(f"Closed spam issue #{issue_number}")


# ── Fetch issues ───────────────────────────────────────────────────────────────

def fetch_issues() -> list[dict]:
    issues = []
    page = 1
    with httpx.Client() as c:
        while True:
            r = c.get(
                f"{GH_API}/repos/{OWNER}/{REPO}/issues",
                headers=GH_HEADERS,
                params={"state": "open", "per_page": 50, "page": page,
                        "sort": "created", "direction": "desc"},
            )
            r.raise_for_status()
            batch = r.json()
            if not batch:
                break
            # Exclude pull requests
            issues.extend(i for i in batch if "pull_request" not in i)
            page += 1
    return issues


# ── Claude analysis ────────────────────────────────────────────────────────────

def analyse_issues(issues: list[dict]) -> list[dict]:
    if not ANTHROPIC_KEY or not issues:
        return [{"number": i["number"], "title": i["title"],
                 "labels": [lb["name"] for lb in i.get("labels", [])],
                 "suggestion": "No Claude API key — manual review needed.",
                 "action": "review", "priority": "medium"} for i in issues]

    client = anthropic.Anthropic(api_key=ANTHROPIC_KEY)

    issue_list = "\n\n".join(
        f"Issue #{i['number']}: {i['title']}\n"
        f"Labels: {', '.join(lb['name'] for lb in i.get('labels', [])) or 'none'}\n"
        f"Body: {(i.get('body') or '')[:400]}"
        for i in issues
    )

    prompt = f"""You are the maintainer of Oligolia, an open-source gene editing desktop app.
Review these open GitHub issues and for each provide:
1. A one-line assessment
2. Recommended action: implement | investigate | close-duplicate | close-wontfix | needs-more-info
3. Priority: high | medium | low
4. Estimated effort: small (hours) | medium (days) | large (weeks)

Issues:
{issue_list}

Respond with a JSON array:
[{{"number": 1, "assessment": "...", "action": "implement", "priority": "high", "effort": "small"}}]
Only the JSON array, no other text."""

    msg = client.messages.create(
        model="claude-sonnet-4-6",
        max_tokens=2000,
        messages=[{"role": "user", "content": prompt}],
    )
    try:
        analyses = json.loads(msg.content[0].text)
    except Exception:
        analyses = []

    # Merge with original issue data
    analysis_map = {a["number"]: a for a in analyses}
    result = []
    for issue in issues:
        n = issue["number"]
        a = analysis_map.get(n, {})
        result.append({
            "number": n,
            "title": issue["title"],
            "url": issue["html_url"],
            "labels": [lb["name"] for lb in issue.get("labels", [])],
            "created": issue["created_at"][:10],
            "author": issue["user"]["login"],
            "assessment": a.get("assessment", "Needs review"),
            "action": a.get("action", "review"),
            "priority": a.get("priority", "medium"),
            "effort": a.get("effort", "unknown"),
        })
    return result


# ── HTML report ────────────────────────────────────────────────────────────────

ACTION_COLOR = {
    "implement": "#16a34a", "investigate": "#2563eb",
    "close-duplicate": "#6b7280", "close-wontfix": "#6b7280",
    "needs-more-info": "#d97706", "review": "#7c3aed",
}
PRIORITY_EMOJI = {"high": "🔴", "medium": "🟡", "low": "🟢"}


def build_html(issues: list[dict], spam_count: int) -> str:
    rows = ""
    for i in issues:
        action_color = ACTION_COLOR.get(i["action"], "#6b7280")
        rows += f"""
        <tr>
          <td><a href="{i['url']}" style="color:#4ade80">#{i['number']}</a></td>
          <td><a href="{i['url']}" style="color:#e2e8f0">{i['title']}</a></td>
          <td>{PRIORITY_EMOJI.get(i['priority'], '⚪')} {i['priority']}</td>
          <td><span style="background:{action_color};color:white;padding:2px 8px;border-radius:4px;font-size:12px">{i['action']}</span></td>
          <td style="color:#94a3b8">{i['effort']}</td>
          <td style="color:#94a3b8;font-size:13px">{i['assessment']}</td>
        </tr>"""

    return f"""<!DOCTYPE html>
<html>
<head><meta charset="UTF-8">
<style>
  body{{font-family:'Segoe UI',system-ui,sans-serif;background:#030712;color:#e2e8f0;padding:32px;max-width:900px;margin:0 auto}}
  h1{{color:#4ade80;margin-bottom:4px}}
  .meta{{color:#94a3b8;font-size:14px;margin-bottom:32px}}
  table{{width:100%;border-collapse:collapse;margin-bottom:32px}}
  th{{background:#0f172a;color:#94a3b8;font-size:12px;text-transform:uppercase;padding:10px 12px;text-align:left;border-bottom:1px solid #1e293b}}
  td{{padding:10px 12px;border-bottom:1px solid #1e293b;vertical-align:top}}
  tr:hover td{{background:#0a1628}}
  .footer{{color:#475569;font-size:12px;border-top:1px solid #1e293b;padding-top:16px}}
  .respond-box{{background:#0f172a;border:1px solid #1e293b;border-radius:8px;padding:20px;margin-bottom:24px}}
</style>
</head>
<body>
<h1>🧬 Oligolia Issue Digest</h1>
<div class="meta">{TODAY} · {len(issues)} open issues · {spam_count} spam removed this week</div>

<div class="respond-box">
  <strong style="color:#fbbf24">📧 To take action:</strong> Reply to this email with your decisions, e.g.:<br>
  <code style="color:#4ade80">#42: implement · #18: close-wontfix · #99: needs-more-info</code><br><br>
  Or visit <a href="https://github.com/{OWNER}/{REPO}/issues" style="color:#60a5fa">GitHub Issues</a> directly.
</div>

<table>
  <thead>
    <tr>
      <th>#</th><th>Title</th><th>Priority</th><th>Suggested Action</th><th>Effort</th><th>Assessment</th>
    </tr>
  </thead>
  <tbody>{rows}</tbody>
</table>

<div class="footer">
  Generated automatically by Oligolia maintainer bot ·
  <a href="https://github.com/{OWNER}/{REPO}" style="color:#4ade80">GitHub</a> ·
  <a href="https://github.com/{OWNER}/{REPO}/actions" style="color:#4ade80">Actions</a>
</div>
</body>
</html>"""


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    print(f"Fetching issues from {OWNER}/{REPO}…")
    all_issues = fetch_issues()
    print(f"  {len(all_issues)} open issues")

    # Security filter — close spam silently
    spam_count = 0
    clean_issues = []
    for issue in all_issues:
        if is_spam(issue):
            close_spam_issue(issue["number"], "spam")
            spam_count += 1
        else:
            clean_issues.append(issue)
    print(f"  {spam_count} spam/injection issues removed")

    # Claude analysis
    print("Analysing with Claude…")
    analysed = analyse_issues(clean_issues)

    # Sort: high priority first
    order = {"high": 0, "medium": 1, "low": 2}
    analysed.sort(key=lambda x: order.get(x["priority"], 3))

    # Write HTML report
    html = build_html(analysed, spam_count)
    with open("digest_output.html", "w") as f:
        f.write(html)
    print("Wrote digest_output.html")

    # GitHub Actions outputs
    with open(os.environ.get("GITHUB_OUTPUT", "/dev/null"), "a") as f:
        f.write(f"date={TODAY}\n")
        f.write(f"issue_count={len(analysed)}\n")


if __name__ == "__main__":
    main()
