#!/bin/bash
# `config/autoconfigure -r` is meant to auto-detect the machine/compiler and
# then actually run `config/configure`. It has a bash-ism (`[ $1 == "-r" ]`)
# that silently no-ops under a POSIX /bin/sh (e.g. Ubuntu's dash) -- the
# script calls plain `sh`, not bash, so that's what matters here.
#
# On macOS, /bin/sh is itself bash-derived and tolerates `==`, so this won't
# reproduce the bug locally unless `dash` is installed and used explicitly:
#   dash Freemol/config/autoconfigure -r </dev/null
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FREEMOL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$FREEMOL_DIR"
out="$(sh ./config/autoconfigure -r </dev/null 2>&1)"
rc=$?
echo "$out"

if [[ $rc -eq 0 ]] && [[ "$out" == *"FREEMOL configure NORMAL TERMINATION"* ]]; then
    echo "check_autoconfigure: PASS"
    exit 0
else
    echo "check_autoconfigure: FAIL (exit $rc, or missing NORMAL TERMINATION)"
    exit 1
fi
