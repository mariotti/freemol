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
    make freemol
Requirements: gfortran, perl (bin/makemake*), csh (bin/cleanup, used by
make cleanall), python3 (tests), octave (tools/CSMD only). adfrom is
skipped by design without the commercial ADF libraries.
New program dirs are only picked up after re-running configure or
`./config/makeprograms programs > programs/Makefile`.

## Tests (Freemol/tests/, wired into .github/workflows/ci.yml)
build_all.sh, check_autoconfigure.sh, run_csmg_regression.sh (+
numdiff.py), run_smoke.sh, run_devflow.sh, check_clean.sh,
check_csmd_tools.sh. Run the full suite locally before pushing; see
README.md's "Build and test" section for what each one actually checks.
`config/configure` wipes and regenerates Freemol/MODFILES on every run,
and stale-but-unchanged .o files won't get recompiled to refill it -- if
a rebuild inexplicably fails to find a .mod file, `make cleanall` and
rebuild from scratch rather than debugging further.

## Programs and test data
(paths below are relative to Freemol/)
- CSMG: data/CSMG/tests/*.mld (59 cases, reads stdin). Reference outputs:
  data/CSMG/outs/osx/*.out (8, macOS 2016 -- independent cross-decade
  validation) and data/CSMG/outs/generated/*.out (the other 51 -- a
  self-generated baseline, only proves reproducibility going forward,
  see run_csmg_regression.sh's header comment).
- fit1Dpol: data/fit1Dpol/examples/*.inp
- ch4sym2cart: data/ch4sym2cart/tests/ch4_s1.inp + {s2x,s2y,s2z,s2a,s2b,
  s4x,s4y,s4z,sr}_only.inp (one displacement dimension each -- note
  ch4sym2cart's r/angle grouping of these 10 fields is NOT the same as
  XY4Coord's, despite identical field names/order).
- XY4Coord: data/XY4Coord/tests/ch4_s1.inp, equilibrium.inp +
  {s2x,s2y,s2z}_only.inp (pass directly) and
  {s2a,s2b,s4x,s4y,s4z}_only.inp (pass via the Sr redundancy solve, see
  run_smoke.sh's xy4coord_redundancy_selfcheck).
- XY4PolySphere: data/XY4PolySphere/tests/ch4_poly.inp, ch4_poly_random.inp
- tools/CSMD/*.m are Octave (not MATLAB) scripts.

## Known issues (see README.md for the full evidence, PRECISION_NOTES.md
for the numerical mechanism behind the do_checks() one and other
compiler-noise findings from CI)
- XY4Coord: angular-type displacements (S2a, S2b, S4x, S4y, S4z) all pass
  now -- were broken (eval_sr's redundancy solve used a cosine-based Sr
  from an unavailable paper, while get_ra/get_symc use a
  radian-displacement-based Sr), fixed by root-finding Sr directly
  against do_checks()'s own Gamma-Sum closure condition instead
  (XY4coord.F90, see PRECISION_NOTES.md section 3). Radial-type
  displacements (S1, S2x, S2y, S2z) pass too, as before. `Sr` alone is
  still expected to fail (geometrically impossible, not a bug -- raising
  all six angles together violates the bond-vectors' Gram matrix
  positive-semidefiniteness). A related copy-paste bug (H4 test checking
  H2's formula) was fixed without a test -- proven to have no effect on
  do_checks()'s verdict for any reachable input, see README.
- ch4sym2cart's "Unchenged Coordina..." diagnostic lines are
  geometrically wrong (missing a degrees->radians conversion);
  diagnostic-only output, not fixed.

## mcp/ (Python MCP server, not Fortran)
Wraps XY4PolySphere/XY4Coord/ch4sym2cart's built binaries as MCP tools
(subprocess + parse output; no numerics reimplemented), plus
list_freemol_sections/read_freemol_section: a pure-Python mirror of
osec_set's [section-name] scanning convention (no binary involved --
text scanning isn't numerics -- so these two work without Freemol/bin
built). The "don't change Fortran except to fix a bug" rule above
doesn't apply here -- it's a separate Python package. Has its own tests
(`mcp/tests/`, pytest) run in CI as the `mcp` job; the three
binary-wrapping tools need `Freemol/bin/*.exe` built first, the two
section tools don't. See `mcp/README.md`.

## Versioning and releases
Semver via git tags (`v1.1.1` is current). `Freemol/config/printversion`'s
`distnum` (and `vernum`/`vernumdist`) should be bumped to match before
tagging a new release. Pushing a `v*.*.*` tag (or a manual
`workflow_dispatch`) runs `.github/workflows/release.yml`, publishing
Linux x86_64 + macOS arm64 binaries as a GitHub Release.
