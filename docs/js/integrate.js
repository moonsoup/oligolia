// Workflow Integration Request form — submits to a small Cloudflare Worker
// (cloudflare-worker/issue-proxy.js in this repo) that holds a repo-scoped
// GitHub token server-side and files the issue via the GitHub API. The
// browser never sees a token and never navigates to github.com — the whole
// flow stays on this page, including the success/error message.
//
// WORKER_URL must point at your deployed Worker (see
// cloudflare-worker/README.md for setup). Until it's deployed, submissions
// will fail gracefully to the error box, which links to the plain GitHub
// template as a fallback.
const WORKER_URL = "https://oligolia-issue-proxy.YOUR-SUBDOMAIN.workers.dev/submit";

const form = document.getElementById("integrate-form");
const submitBtn = document.getElementById("submit-btn");
const submitStatus = document.getElementById("submit-status");
const successBox = document.getElementById("success-box");
const errorBox = document.getElementById("error-box");
const successLink = document.getElementById("success-link");
const successNum = document.getElementById("success-num");

let submitted = false; // guards against double-submit (double-click, Enter+click)

form.addEventListener("submit", async (e) => {
  e.preventDefault();
  if (submitted) return;

  // Honeypot: real visitors never fill this hidden field in. If it's
  // filled, silently pretend to succeed rather than tipping off a bot.
  const honeypot = document.getElementById("website").value.trim();
  if (honeypot) {
    submitted = true;
    showSuccess(null, null);
    return;
  }

  if (!form.checkValidity()) {
    form.reportValidity();
    return;
  }

  const steps = Array.from(
    document.querySelectorAll('#steps-needed input[type="checkbox"]:checked')
  ).map((cb) => cb.value);

  const payload = {
    title: document.getElementById("title").value.trim(),
    org: document.getElementById("org").value.trim(),
    system_type: document.getElementById("system-type").value,
    steps_needed: steps,
    input_shape: document.getElementById("input-shape").value.trim(),
    output_shape: document.getElementById("output-shape").value.trim(),
    example_workflow: document.getElementById("example-workflow").value.trim(),
    constraints: document.getElementById("constraints").value.trim(),
    priority: document.getElementById("priority").value,
  };

  submitBtn.disabled = true;
  submitStatus.textContent = "Submitting…";
  errorBox.hidden = true;

  try {
    const res = await fetch(WORKER_URL, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const data = await res.json();
    if (!data.ok) throw new Error(data.error || "unknown error");

    submitted = true;
    showSuccess(data.issue_url, data.issue_number);
  } catch (err) {
    submitBtn.disabled = false;
    submitStatus.textContent = "";
    errorBox.hidden = false;
    console.error("Workflow integration submission failed:", err);
  }
});

function showSuccess(issueUrl, issueNumber) {
  form.querySelectorAll("input, select, textarea, button").forEach((el) => {
    el.disabled = true;
  });
  submitStatus.textContent = "";
  if (issueUrl) {
    successLink.href = issueUrl;
    successNum.textContent = issueNumber;
  } else {
    // Honeypot path — no real issue was created; keep the message generic.
    successLink.href = "https://github.com/moonsoup/oligolia/issues";
    successNum.textContent = "";
  }
  successBox.hidden = false;
  successBox.scrollIntoView({ behavior: "smooth", block: "center" });
}
