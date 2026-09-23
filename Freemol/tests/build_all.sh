#!/bin/sh
# Configure and build the Freemol libraries and programs with gfortran.
# Usage: Freemol/tests/build_all.sh [machine]   (default: m_generic_linux)
# adfrom is not built: it needs the ADF libraries.
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
MACHINE=${1:-m_generic_linux}
cd "$ROOT"
./config/configure "$MACHINE" gfortran "$ROOT"
for target in others includes utilities moduledata modules programs; do
    echo "=== make $target"
    make "$target"
done
ls -l bin/*.exe
