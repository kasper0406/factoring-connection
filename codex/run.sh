#!/bin/bash
# Usage: codex/run.sh <name> [effort]   (prompt in codex/<name>.md; writes <name>.out/.log)
# Repo is read-only for codex; scratch writes go to codex/scratch/<name>/ (git-ignored) and /tmp.
set -u
here="$(cd "$(dirname "$0")" && pwd)"
name="$1"; effort="${2:-high}"
mkdir -p "$here/scratch/$name"
cd "$here/scratch/$name"
perl -e 'alarm shift; exec @ARGV' 4800 \
  codex exec -s workspace-write -m gpt-6.1-sol -c model_reasoning_effort="\"$effort\"" \
  -c tools.web_search=true --skip-git-repo-check \
  -o "$here/$name.out" - < "$here/$name.md" > "$here/$name.log" 2>&1
