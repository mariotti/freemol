#!/bin/bash
# After a build, verify `make cleanall`/`make cleanconfig` (needs csh) leave
# the working tree clean: git status --porcelain --ignored must show nothing
# beyond files these tests themselves may create.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FREEMOL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_DIR="$(cd "$FREEMOL_DIR/.." && pwd)"

cd "$FREEMOL_DIR"
make cleanall
make cleanconfig

cd "$REPO_DIR"
leftovers="$(git status --porcelain --ignored)"
if [[ -n "$leftovers" ]]; then
    echo "check_clean: working tree not clean after cleanall/cleanconfig:" >&2
    echo "$leftovers" >&2
    exit 1
fi

echo "check_clean: PASS (working tree clean)"
