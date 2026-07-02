# Oligolia issue-proxy Worker

Backs the "Integrate Your Workflow" form (`docs/integrate.html`). Holds a
GitHub token server-side and files a real issue on `moonsoup/oligolia` when
someone submits — the browser never sees the token and never leaves the
Oligolia site.

This is the only piece of live infrastructure this project runs; everything
else is a static Pages site + a desktop app. Free tier covers this easily
(Cloudflare Workers free tier: 100,000 requests/day).

## One-time setup (you do this — Claude can't create accounts or hold your token)

### 1. Generate a GitHub fine-grained token, scoped to only this repo

1. GitHub → Settings → Developer settings → **Personal access tokens →
   Fine-grained tokens** → Generate new token.
2. **Repository access:** "Only select repositories" → `moonsoup/oligolia`.
3. **Permissions:** Repository permissions → **Issues: Read and write**.
   Leave everything else at "No access."
4. Generate, copy the token. You'll paste it directly into Cloudflare in
   step 3 below — not into a chat, not into any file in this repo.

### 2. Create a Cloudflare account + install Wrangler (if you haven't already)

```bash
npm install -g wrangler
wrangler login   # opens a browser to authenticate
```

### 3. Deploy the Worker and set the secret

```bash
cd cloudflare-worker
wrangler deploy
wrangler secret put GITHUB_TOKEN
# paste the token from step 1 when prompted — it's stored encrypted by
# Cloudflare, never written to disk here, never visible in `wrangler.toml`
```

`wrangler deploy` prints the Worker's URL, something like:

```
https://oligolia-issue-proxy.<your-subdomain>.workers.dev
```

### 4. Wire the URL into the form

Edit `docs/js/integrate.js`, replace the `WORKER_URL` constant's placeholder
with `https://oligolia-issue-proxy.<your-subdomain>.workers.dev/submit`,
commit, push. That's it — no further changes needed here.

## Testing without spamming the real repo

```bash
curl -X POST https://oligolia-issue-proxy.<your-subdomain>.workers.dev/submit \
  -H "Content-Type: application/json" \
  -d '{"title":"test","org":"test org","system_type":"LIMS (lab information management system)","input_shape":"x","output_shape":"y"}'
```

This **will** file a real issue (labeled `workflow-integration` +
`needs-review`) — close/delete it afterward, or point `REPO_OWNER`/
`REPO_NAME` in `issue-proxy.js` at a scratch repo temporarily while testing.

## Abuse posture

This is a public, unauthenticated POST endpoint that creates content on a
real public repo. Current mitigations: a honeypot field (checked both
client- and server-side), field-length limits, and Cloudflare's own
edge-level bot/DDoS protection. If spam becomes a real problem, the natural
next step is adding Cloudflare Turnstile (free CAPTCHA alternative) to the
form and verifying the token server-side before creating the issue — not
implemented yet since there's no evidence it's needed.

## Rotating or revoking the token

Fine-grained tokens can be revoked instantly from GitHub Settings →
Developer settings → Fine-grained tokens, with zero impact on anything else
in the project (it has no access beyond this one repo's issues).
