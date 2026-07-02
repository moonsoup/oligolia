/**
 * Oligolia workflow-integration-request proxy.
 *
 * Receives the JSON payload from docs/integrate.html, files it as a GitHub
 * issue via the REST API using a repo-scoped token held only as a Worker
 * secret, and returns the created issue's URL/number. The browser never
 * sees the token and never navigates to github.com — see docs/js/integrate.js.
 *
 * Required Worker secret (set via `wrangler secret put GITHUB_TOKEN`,
 * never committed): a fine-grained GitHub PAT scoped to ONLY
 * moonsoup/oligolia, with "Issues: Read and write" and nothing else.
 * See README.md in this directory for full setup steps.
 */

const REPO_OWNER = "moonsoup";
const REPO_NAME = "oligolia";
const ALLOWED_ORIGIN = "https://moonsoup.github.io";
const LABELS = ["workflow-integration", "needs-review"];

// Defensive limits — this is a public, unauthenticated endpoint.
const MAX_FIELD_LEN = 4000;
const MAX_TITLE_LEN = 200;
const MAX_STEPS = 10;

export default {
  async fetch(request, env) {
    if (request.method === "OPTIONS") return corsPreflight();
    if (request.method !== "POST") {
      return json({ ok: false, error: "method not allowed" }, 405);
    }

    let body;
    try {
      body = await request.json();
    } catch {
      return json({ ok: false, error: "invalid JSON" }, 400, true);
    }

    // Server-side honeypot check — client JS already does this, but a bot
    // posting directly to this endpoint skips the browser entirely.
    if (typeof body.website === "string" && body.website.trim() !== "") {
      // Pretend success; create nothing.
      return json({ ok: true, honeypot: true }, 200, true);
    }

    const err = validate(body);
    if (err) return json({ ok: false, error: err }, 422, true);

    const title = `Workflow Integration: ${clip(body.title, MAX_TITLE_LEN)}`;
    const bodyMd = renderBody(body);

    const ghRes = await fetch(
      `https://api.github.com/repos/${REPO_OWNER}/${REPO_NAME}/issues`,
      {
        method: "POST",
        headers: {
          Authorization: `Bearer ${env.GITHUB_TOKEN}`,
          Accept: "application/vnd.github+json",
          "X-GitHub-Api-Version": "2022-11-28",
          "User-Agent": "oligolia-issue-proxy",
        },
        body: JSON.stringify({ title, body: bodyMd, labels: LABELS }),
      }
    );

    if (!ghRes.ok) {
      const detail = await ghRes.text().catch(() => "");
      console.error("GitHub API error", ghRes.status, detail);
      return json({ ok: false, error: "GitHub API error" }, 502, true);
    }

    const issue = await ghRes.json();
    return json(
      { ok: true, issue_url: issue.html_url, issue_number: issue.number },
      200,
      true
    );
  },
};

function validate(b) {
  if (!isNonEmptyStr(b.title)) return "title is required";
  if (!isNonEmptyStr(b.org)) return "org is required";
  if (!isNonEmptyStr(b.system_type)) return "system_type is required";
  if (!isNonEmptyStr(b.input_shape)) return "input_shape is required";
  if (!isNonEmptyStr(b.output_shape)) return "output_shape is required";
  if (b.steps_needed !== undefined) {
    if (!Array.isArray(b.steps_needed) || b.steps_needed.length > MAX_STEPS) {
      return "steps_needed must be an array (max 10 items)";
    }
  }
  for (const field of ["title", "org", "system_type", "input_shape", "output_shape",
                        "example_workflow", "constraints", "priority"]) {
    if (typeof b[field] === "string" && b[field].length > MAX_FIELD_LEN) {
      return `${field} is too long (max ${MAX_FIELD_LEN} chars)`;
    }
  }
  return null;
}

function isNonEmptyStr(v) {
  return typeof v === "string" && v.trim().length > 0;
}

function clip(s, n) {
  return (s || "").slice(0, n);
}

function renderBody(b) {
  const steps = Array.isArray(b.steps_needed) && b.steps_needed.length
    ? b.steps_needed.map((s) => `- [x] ${s}`).join("\n")
    : "_none selected_";

  return `## Organization / lab name
${b.org}

## What are you integrating with?
${b.system_type}

## Which existing Workflow steps does this involve?
${steps}

## What input format do you need Oligolia to accept?
${b.input_shape}

## What output format do you need Oligolia to produce?
${b.output_shape}

## Example end-to-end workflow
${isNonEmptyStr(b.example_workflow) ? b.example_workflow : "_not provided_"}

## Constraints
${isNonEmptyStr(b.constraints) ? b.constraints : "_not provided_"}

## Priority
${isNonEmptyStr(b.priority) ? b.priority : "_not provided_"}

---
_Submitted via the Integrate Your Workflow form on the Oligolia site, not the raw GitHub issue form._`;
}

function corsPreflight() {
  return new Response(null, {
    status: 204,
    headers: {
      "Access-Control-Allow-Origin": ALLOWED_ORIGIN,
      "Access-Control-Allow-Methods": "POST, OPTIONS",
      "Access-Control-Allow-Headers": "Content-Type",
    },
  });
}

function json(data, status, cors = false) {
  const headers = { "Content-Type": "application/json" };
  if (cors) headers["Access-Control-Allow-Origin"] = ALLOWED_ORIGIN;
  return new Response(JSON.stringify(data), { status, headers });
}
