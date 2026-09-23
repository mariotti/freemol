#!/bin/bash
# Configure and build freemol, stopping at the last working target.
#
# Usage: build_all.sh [machine]   (default: m_generic_linux)
#
# Does NOT run `make freemol`: its helpdocs step is broken (Freemol/help has
# no Makefile) -- see CLAUDE.md. This builds everything that step depends on.
set -euo pipefail

MACHINE="${1:-m_generic_linux}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FREEMOL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$FREEMOL_DIR"
./config/configure "$MACHINE" gfortran "$PWD"
make others includes utilities moduledata modules programs
