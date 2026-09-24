#!/bin/bash
# Proves the README "add a program" workflow still works: create a program
# dir with a Makelocal, regenerate programs/Makefile, build it, run it.
#
# Always cleans up afterwards (even on failure): removes the program dir,
# its built binary and docs entry, and regenerates programs/Makefile so the
# working tree is left exactly as it was found.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FREEMOL_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PROG="zz_devtest"
PROG_DIR="$FREEMOL_DIR/programs/$PROG"

cleanup() {
    rm -rf "$PROG_DIR"
    rm -f "$FREEMOL_DIR/bin/$PROG.exe"
    rm -f "$FREEMOL_DIR/docs/$PROG.exe.txt"
    (cd "$FREEMOL_DIR" && ./config/makeprograms programs > programs/Makefile)
}
trap cleanup EXIT

mkdir -p "$PROG_DIR"

cat > "$PROG_DIR/Makelocal" <<'EOF'
#
# Makelocal for zz_devtest (Freemol/tests/run_devflow.sh)
#
HERELIBS =
HEREPROG = zz_devtest.exe
HEREOBJS = zz_devtest.o
USEDLIBS = -L$(FMLIBDIR) -lmodules -lmoduledata -lutils -lincludes -lfmlapack -lfmblas
EOF

cat > "$PROG_DIR/zz_devtest.F90" <<'EOF'
program zz_devtest
  use vartypes
  use messages
  use pcmdline
  use baseio
  use osec
  use molecule
  implicit none

  integer(FINT) :: irc
  integer(FINT) :: nargs
  character(FLCHARS), dimension(:), allocatable :: cargs
  character(FLCHARS) :: params
  character(FLCHARS) :: finput, foutput
  integer(FINT) :: iunin, iunout

  call messages_init()
  irc = baseio_init()

  call pcmd_iargc(nargs)
  allocate(cargs(nargs), STAT=irc)
  call pcmd_getarg(nargs, cargs)
  call pcmd_getio(nargs, cargs, finput, foutput)

  irc = baseio_open(iunin, finput, stat='OLD')
  if (irc.lt.0) then
     call message(MESERRO, '[zz_devtest] Cannot open input file.')
     stop 1
  end if
  iunout = 6

  rewind(iunin)
  irc = osec_set(iunin, 'molecule', params)
  if (irc.lt.0) then
     call message(MESERRO, '[zz_devtest] No molecule section.')
     stop 1
  end if

  irc = molecule_init()
  if (irc.lt.0) then
     call message(MESERRO, '[zz_devtest] Cannot initialize the molecule.')
     stop 1
  end if

  call molecule_read(iunin, params)
  call molecule_print(iunout)
  call message_value(MESOUT, '[zz_devtest] Atom count:', molecule_getnat())

end program zz_devtest
EOF

cd "$FREEMOL_DIR"
./config/makeprograms programs > programs/Makefile
# `make programs` alone doesn't reliably rebuild the MODFILES chain (the
# legacy submake for one lib can leave another lib's .mod stale/missing),
# so rebuild the full library chain explicitly (skips helpdocs, unneeded here).
make others includes utilities moduledata modules programs

exe="$FREEMOL_DIR/bin/$PROG.exe"
if [[ ! -x "$exe" ]]; then
    echo "run_devflow: FAIL ($exe was not built)"
    exit 1
fi

out="$("$exe" -i "$FREEMOL_DIR/data/ftemplate/C3v.mld" -o /dev/stdout 2>&1)"
echo "$out"

count="$(echo "$out" | sed -n 's/.*Atom count:[[:space:]]*//p' | tail -1)"
if [[ "$count" != "4" ]]; then
    echo "run_devflow: FAIL (expected 4 atoms, got '$count')"
    exit 1
fi

echo "run_devflow: PASS (4 atoms)"
