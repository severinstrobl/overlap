#!/usr/bin/env python3

"""Script to convert pre-commit config repo tags to commit hashes.

This script reads the .pre-commit-config.yaml file, replaces repo tag references with their
corresponding commit hashes using the GitHub API, and updates the file if changes are made.
"""

import os
import re
import sys
from pathlib import Path

import requests

PATTERN = r"repo:\s+(https://github\.com/[^\s]+)\s+rev:\s+([^\s]+)(?:\s*#\s*frozen:\s*([^\s]+))?"


def get_commit_hash_for_tag(repo_url: str, tag: str) -> str | None:
    """Get commit hash for a specific tag from GitHub API."""
    match = re.search(r"github\.com/([^/]+)/([^/]+?)(?:\.git)?$", repo_url)
    if not match:
        return None

    owner, repo = match.groups()

    headers = {}
    github_token = os.getenv("GITHUB_TOKEN")
    if github_token:
        headers["Authorization"] = f"token {github_token}"

    api_url = f"https://api.github.com/repos/{owner}/{repo}/git/refs/tags/{tag}"
    try:
        response = requests.get(api_url, headers=headers)
        response.raise_for_status()
        tag_data = response.json()
        if tag_data["object"]["type"] == "commit":
            return tag_data["object"]["sha"]
        elif tag_data["object"]["type"] == "tag":
            tag_response = requests.get(tag_data["object"]["url"], headers=headers)
            tag_response.raise_for_status()
            return tag_response.json()["object"]["sha"]
    except Exception:
        print(f"❌ Error fetching hash for {repo_url} tag {tag}", file=sys.stderr)
        return None


def replace_with_hash(match: re.Match) -> str:
    """Replace the rev value in a pre-commit config repo entry with its commit hash."""
    repo_url = match.group(1)
    current_rev = match.group(2)
    frozen_version = match.group(3) if match.group(3) else current_rev

    if re.match(r"^[a-f0-9]{40}$", current_rev):
        new_hash = get_commit_hash_for_tag(repo_url, frozen_version)
        if new_hash and new_hash != current_rev:
            print(f"🔄 Updating hash for {repo_url} {frozen_version}")
            return f"repo: {repo_url}\n    rev: {new_hash}  # frozen: {frozen_version}"
        return match.group(0)
    else:
        commit_hash = get_commit_hash_for_tag(repo_url, current_rev)
        if commit_hash:
            print(f"✅ Converting {repo_url} tag {current_rev} to hash")
            return f"repo: {repo_url}\n    rev: {commit_hash}  # frozen: {current_rev}"
        else:
            print(
                f"❌ Failed to get hash for {repo_url} tag {current_rev}",
                file=sys.stderr,
            )
            return match.group(0)


config_path = Path(".pre-commit-config.yaml")
if not config_path.exists():
    print("🚫 No .pre-commit-config.yaml found", file=sys.stderr)
    sys.exit(1)

content = config_path.read_text()
updated_content = re.sub(PATTERN, replace_with_hash, content, flags=re.MULTILINE)

if updated_content != content:
    config_path.write_text(updated_content)
    print("✨ Pre-commit config updated with commit hashes")
else:
    print("🟢 No changes needed")
