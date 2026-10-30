#!/usr/bin/env python3
"""
upload_to_zenodo.py: REST API upload for the CSPI Zenodo deposition.

Why a script: the 30 m raster is ~30 GB and lives on the remote server (it is
not practical to pull it locally first). Run this script on firebreather (the
default host; Cardinal only if Aaron names it) to upload everything directly
via the Zenodo REST API. Draft uploads only with Aaron's explicit yes; publishing
is Aaron's action.

Authentication: token is sent in the Authorization header as a Bearer token,
never in the URL query string. This means the token cannot leak into request
logs, server access logs, error tracebacks, or shell history.

Setup (one time):
  1. Create or log in to your Zenodo account at https://zenodo.org
  2. Go to Applications -> Personal access tokens, create a token with the
     scopes `deposit:write` and `deposit:actions`. Copy the token.
  3. Save the token to a file readable only by you (see SKILL.md for the search
     order; on firebreather, fire_upload it into the job folder, chmod 600, and
     shred -u it after the verified upload):
        printf '%s' 'YOUR_TOKEN_HERE' > ~/.zenodo_token && chmod 600 ~/.zenodo_token
  4. Edit `zenodo_metadata.json` and replace any remaining PLACEHOLDER_* fields.

Run (on firebreather, in the job folder; requests is available to python3 there.
On Cardinal, first: module load python/3.11 && pip install --user requests):
  python3 upload_to_zenodo.py \
      --token-file ./.zenodo_token \
      --metadata zenodo_metadata.json \
      --files-list files_to_upload.txt \
      --publish        # autopublish: Aaron's call only; omit for a draft

Sandbox vs production: --sandbox uses https://sandbox.zenodo.org (test instance).
Default is production https://zenodo.org.

On error, the script prints the server response body so you can see the actual
validation message. The Authorization header is never echoed.
"""

import argparse, json, sys, time
from pathlib import Path
import requests

def _auth_headers(token):
    return {"Authorization": f"Bearer {token}"}

def _bail(resp, what):
    # Print body but never the request headers (which contain the token)
    print(f"\nERROR during {what}: HTTP {resp.status_code}", file=sys.stderr)
    try:
        body = resp.json()
        print(json.dumps(body, indent=2), file=sys.stderr)
    except Exception:
        print(resp.text[:2000], file=sys.stderr)
    sys.exit(2)

def post_deposition(api, token, metadata):
    r = requests.post(
        f"{api}/deposit/depositions",
        headers={**_auth_headers(token), "Content-Type": "application/json"},
        json={"metadata": metadata["metadata"]},
    )
    if not r.ok:
        _bail(r, "deposit creation")
    return r.json()

def upload_file(api, token, bucket_url, path):
    name = path.name
    size_mb = path.stat().st_size / 1e6
    print(f"  uploading {name} ({size_mb:.1f} MB) ...", flush=True)
    with path.open("rb") as fp:
        r = requests.put(
            f"{bucket_url}/{name}",
            headers=_auth_headers(token),
            data=fp,
        )
    if not r.ok:
        _bail(r, f"upload of {name}")
    return r.json()

def publish_deposition(api, token, dep_id):
    r = requests.post(
        f"{api}/deposit/depositions/{dep_id}/actions/publish",
        headers=_auth_headers(token),
    )
    if not r.ok:
        _bail(r, "publish")
    return r.json()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--token-file", required=True, help="path to file containing your Zenodo personal access token")
    ap.add_argument("--metadata", required=True, help="path to zenodo_metadata.json")
    ap.add_argument("--files-list", required=True, help="text file with one absolute path per line")
    ap.add_argument("--sandbox", action="store_true", help="use sandbox.zenodo.org instead of zenodo.org")
    ap.add_argument("--publish", action="store_true", help="publish the deposition after upload (otherwise leave as draft)")
    args = ap.parse_args()

    token = Path(args.token_file).read_text().strip()
    if not token:
        print("ERROR: token file is empty.", file=sys.stderr); sys.exit(2)

    api = "https://sandbox.zenodo.org/api" if args.sandbox else "https://zenodo.org/api"
    metadata = json.loads(Path(args.metadata).read_text())

    print(f"creating deposition at {api} ...", flush=True)
    dep = post_deposition(api, token, metadata)
    dep_id = dep["id"]
    bucket_url = dep["links"]["bucket"]
    reserved = dep.get("metadata", {}).get("prereserve_doi", {}).get("doi", "(none)")
    print(f"  deposition id: {dep_id}")
    print(f"  edit URL:      {dep['links']['html']}")
    print(f"  reserved DOI:  {reserved}   (citable once published; safe to put in the manuscript, withheld in the blinded file)")
    print()

    base = Path(args.files_list).resolve().parent
    seen = set()
    for raw in Path(args.files_list).read_text().splitlines():
        line = raw.split("#", 1)[0].strip()   # strip inline comments and full-comment lines
        if not line:
            continue
        p = Path(line).expanduser()
        path = (p if p.is_absolute() else (base / p)).resolve()
        if path in seen:
            continue
        seen.add(path)
        if not path.exists():
            print(f"  MISSING: {path} (skipping)", flush=True)
            continue
        upload_file(api, token, bucket_url, path)
        time.sleep(0.5)

    if args.publish:
        print("publishing deposition ...", flush=True)
        rec = publish_deposition(api, token, dep_id)
        print(f"  published. DOI: {rec.get('doi', '(pending)')}")
        print(f"  record URL:    {rec.get('links', {}).get('record_html', dep['links']['html'])}")
    else:
        print(f"\nFinished. Deposition is a DRAFT.")
        print(f"Review and publish at: {dep['links']['html']}")

if __name__ == "__main__":
    main()
