# freemol — notes for Claude Code

Historical Fortran 90 chemistry framework (2003, CVS 2009, GitHub 2015).
The tag Legacy-2016 marks the untouched historical state: never move or delete it,
never force-push master.

## Rules
- Do not delete files unless a task explicitly says so. stubs/ and
  Freemol/others/pFUnit are kept on purpose (historical), even though unused.
- Do not change Fortran code except in a task that explicitly fixes a bug.
  One bug = one commit, with a test that fails before and passes after.
- Keep the existing build system (configure + makemake + make). No CMake/fpm.
- Every change goes through a branch + PR; CI must be green before merge.

## Build (Linux; macOS uses m_generic_osx)
    cd Freemol
    ./config/configure m_generic_linux gfortran $PWD
    make others includes utilities moduledata modules programs
Requirements: gfortran, perl (bin/makemake*), csh (bin/cleanup, used by
make cleanall), python3 (tests). adfrom is skipped by design without the
commercial ADF libraries.
Known: `make freemol` fails at the last step (helpdocs: Freemol/help has no
Makefile). New program dirs are only picked up after re-running configure or
`./config/makeprograms programs > programs/Makefile`.

## Programs and test data
(paths below are relative to Freemol/)
- CSMG: data/CSMG/tests/*.mld, reference outputs data/CSMG/outs/osx/*.out
  (macOS 2016, produced reading stdin).
- fit1Dpol: data/fit1Dpol/examples/*.inp
- ch4sym2cart, XY4Coord, XY4PolySphere: no test data yet.
- tools/CSMD/*.m are Octave (not MATLAB) scripts.
