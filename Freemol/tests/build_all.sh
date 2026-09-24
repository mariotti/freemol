#!/bin/bash
# Configure and build freemol via the documented top-level target.
#
# Usage: build_all.sh [machine]   (default: m_generic_linux)
set -euo pipefail

MACHINE="${1:-m_generic_linux}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FREEMOL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$FREEMOL_DIR"
./config/configure "$MACHINE" gfortran "$PWD"
make freemol
