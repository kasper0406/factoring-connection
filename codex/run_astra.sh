#!/bin/bash
# Usage: codex/run_astra.sh <name> [effort]  (prompt codex/<name>.md -> codex/<name>_codex.out/.log; scratch codex/scratch/<name>_c)
set -u
here="$(cd "$(dirname "$0")" && pwd)"; name="$1"; effort="${2:-xhigh}"
mkdir -p "$here/scratch/${name}_c"; cd "$here/scratch/${name}_c"
perl -e 'alarm shift; exec @ARGV' 9000 \
  codex exec -s workspace-write -m gpt-6.2-astra -c model_reasoning_effort="\"$effort\"" \
  -c tools.web_search=true -c 'sandbox_workspace_write.network_access=true' --skip-git-repo-check \
  -o "$here/${name}_codex.out" - < "$here/$name.md" > "$here/${name}_codex.log" 2>&1
