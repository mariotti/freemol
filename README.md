# freemol

[![CI](https://github.com/mariotti/freemol/actions/workflows/ci.yml/badge.svg)](https://github.com/mariotti/freemol/actions/workflows/ci.yml)

**A Fortran 90 chemistry toolkit from 2003, brought back to life in 2026 --
and a working case study in keeping scientific code reproducible.**

freemol began in 2003 as a small framework plus a handful of tools:
continuous symmetry measures for molecules (CSMG), coordinate
transformations for XY4 molecules such as methane, and a Minuit-based
fitter. The code as it stood before the refresh is preserved in the
[`Legacy-2016`](https://github.com/mariotti/freemol/tree/Legacy-2016) tag.
This repository keeps it building, tested, and honest about what it can
and cannot do.

## Why it is still worth a look

**It reproduces, and CI proves it.** CSMG results match reference
outputs generated on macOS in 2016 to the last printed digit -- today, on
Linux with gfortran 12, 13 and 14 and on current macOS. The comparison
runs on every commit. Eight of these references are independent,
decade-old outputs; the remaining 51 test cases have a baseline captured
from the current code, guarding against future regressions.

**It shows what rescuing legacy scientific code really involves.** The
refresh went from zero tests to CI coverage of every program, and in the
process found bugs that had been sitting in the code for years: a tool
that crashed on every input, a self-check that rejected valid
geometries, a copy-paste error in a geometric validation, and a build
target broken since 2015. Not everything is fixed -- open issues are
documented with their evidence, including one root cause we first got
wrong.

**It treats numerics as evidence, not intent.** Double precision
throughout via kind parameters, with platform-dependent noise,
tolerances and numerical findings written up in
[PRECISION_NOTES.md](PRECISION_NOTES.md).

## What it is not

freemol is not a current production tool. For new work on continuous
symmetry measures, actively developed libraries such as
[cosymlib](https://cosymlib.readthedocs.io) are the better starting point.
The XY4 coordinate tools serve a narrow vibrational-spectroscopy niche;
see [Known issues](#known-issues) for what's still open.

## Start here

- [Build and test](#build-and-test) -- build everything and run the full test suite
- [Known issues](#known-issues) -- what doesn't work yet, and why
- [PRECISION_NOTES.md](PRECISION_NOTES.md) -- what testing taught us about floating point
- [`Legacy-2016`](https://github.com/mariotti/freemol/tree/Legacy-2016) -- the original, untouched

A Fortran 90 toolkit for molecular geometry and symmetry: computing the
Continuous Symmetry Measure (CSM) of a set of weighted points, converting
between Cartesian and internal/polyspherical coordinates for XY4-type
(methane-like) molecules, and fitting 1D data with Minuit -- all driven
by a shared, sectioned INI-style input format that predates molden's own
(see [Sectioned input files](#sectioned-input-files-an-old-ini-style-format-from-before-molden)
below). Started in 2003
(as "Freemol2000"), moved through CVS in 2009 and onto GitHub in 2015 --
the [`Legacy-2016`](https://github.com/mariotti/freemol/tree/Legacy-2016) tag marks that untouched
historical state. Since 2026 it also has a CI-verified build and tagged
releases; see [Build and test](#build-and-test) and
[Releases](#releases) below.

## Sectioned input files: an old INI-style format, from before molden

Every program above reads the same underlying file convention: a plain
ASCII file made of `[section-name] inline params` header lines, each
followed by whatever data that section expects, read in free format --
not fixed columns -- until the next `[...]` header or end of file. It
reads like a Windows-style `.ini` file with data attached to each
section rather than just `key=value` pairs. `[molecule] nrec=5
format=nscxyz`, `[x-xy4-symmcoord]`, `[x-csmg-cgauss]` in the examples
throughout this README are all this same syntax -- just different
section names, understood by different programs.

The convention is older than it looks: it predates molden's now-familiar
`[Atoms]`/`[GTO]`-style section format, and molden and the chemistry
file formats that followed its lead effectively inherited the idea from
tools like this one, not the other way around. The code still carries
that family resemblance -- `molecule.F90` has a whole "MOLDEN format
compatibility" branch for writing the same geometry block as an
`[Atoms]` section instead of freemol's own `[molecule]`, right next to
the routine that writes freemol's native form.

Mechanically (`osec_set`, `modules/osec.F90`): given a section name, scan
forward from the file's current read position for a line starting with
`[`; if the name up to the closing `]` matches (case-insensitively),
stop there and hand back whatever trails the `]` on that same line as
the section's own inline parameters (e.g. `nrec=5 format=nscxyz`). The
caller then reads that section's data lines itself with `linetools`, in
free format -- mixed integers/reals/strings on one line, adapted from
CCL's FREEREAD (see Includes below). Nothing about the format enforces
a fixed row shape per section, which is why one input file can carry a
`[molecule]` block of per-atom xyz rows next to an `[x-csmg-symop]`
block of symmetry-operation rows without either reader knowing anything
about the other's row shape ahead of time.

See [Adding a program](#adding-a-program) below for reading a section
from new Fortran code, and [Libraries](#libraries) below (`osec`,
`linetools`) for the modules involved.

## Careful, portable, reproducible double precision

This is old code, but precision was taken seriously from the start, and
that's worth calling out rather than losing in a pile of `real*8`:

- **Kind-parametrised types.** `Freemol/includes/vartypes.F90` defines
  `FREAL = kind(1.0d0)`, used consistently as `..._FREAL` throughout the
  framework and the main programs, so changing precision there is a
  one-line decision, not a search-and-replace. (`SREAL` and `BREAL` are
  defined alongside `FREAL` for future differentiation -- today all
  three happen to resolve to the same double precision.) Double
  precision holds throughout the codebase, but not every file goes
  through `FREAL` to get there: ~230 bare `d0`-style literals remain,
  mostly in `programs/adfrom` (which does `use vartypes` but still
  writes raw `d0`) and in the vendored NSWC quartic solver
  (`programs/XY4Coord/qtcrtmods.F90`, which defines its own `dp` kind via
  `SELECTED_REAL_KIND(15, 60)` rather than using `FREAL` at all).
- **Machine-aware tolerances.** `F_EPS = 10 * EPSILON(1.0_FREAL)` (and
  `S_EPS`, `B_EPS` likewise) instead of a hardcoded epsilon -- the
  tolerance always tracks whatever `FREAL` actually is on the build
  machine, rather than assuming it.
- **pi to 44 significant digits**
  (`3.1415926535897932384626433832795028841971694_FREAL`, exact well
  beyond double precision), used consistently as `LPI` across CSMG,
  ch4sym2cart, XY4Coord and XY4PolySphere -- one constant, never a source
  of rounding error on its own.
- **Evidence, not just intent.** `Freemol/tests/run_csmg_regression.sh`
  compares CSMG's CSM value (printed to 11 decimal digits) against stored
  reference output, on every CI run -- currently Linux with gfortran
  12/13/14 and macOS with current Homebrew gfortran. 8 of the 59 covered
  cases are reference output generated on macOS in 2016: those match to
  the last printed digit across a decade, two operating systems and three
  compiler versions. The other 51 (previously untested entirely) are a
  self-generated baseline captured from this codebase's own current
  output, so they only prove reproducibility going forward from here, not
  across a decade -- still worth having, since 51 of 59 CSMG inputs had
  zero regression coverage before. Either way, it's not a claim, it's a
  passing test you can rerun yourself. (See [Known issues](#known-issues)
  below for the places numerics *aren't* solid yet -- being honest about
  those is part of the same discipline. [PRECISION_NOTES.md](PRECISION_NOTES.md)
  goes deeper: real compiler-noise findings from CI and the numerical
  mechanism behind the do_checks() issue below.)

## Copyright notice/Licences: Please read this
Please note that some code lines might be a copy of other sources,
usually free source code, as this project started many years ago.

Minuit in particular: I started it from an original F77 (free from the 90s) code
and translated into F90. I contacted CERN (current minuit copyright holder) for a confirmation but I could not get a clear answer. I think it is obvious as it is a really old code which was free on the net many years ago.

There is also some free code from CCL (Computational Chemistry List), it is mentioned as comment in the code.

If code, packages or else might not fit with this "main" repository licence they are present in the folder "others/". Please read each package or code licence. If a package/code is there, then the redistribution is granted but under the given package licence which you should read.

### Why so complicated?

There is a different solution: get the required packages from other/external sources.

But I want to create a "self-consistent" package which will run without external dependencies. For these reasons:

- It might need to compile for systems without a network connection.
- If a hardware system (for example) is not supported by the main distribution we can create patches which apply only to the given distribution.
- Security: it is better to download, check, compile, run, check, put in "others/", check, compile, run, check than a simple download from the net.
- Then the obvious: I promote Fortran as the most portable source code. I need all code in a unique package.

Warning: for example BLAS and LAPACK are definitely more optimised if you use your system's own.

Also linetools (described under Libraries) has this top comment:

    !H----------------------------------------------------------------------
    !H      Some Routine are adapted from FREEREAD package taken from CCL.
    !H      CCL: Computational Chemistry List http://www.ccl.org/
    !H      A statment from Jan Labanowski:
    !H      This software was taken from anon. ftp on rani.chem.yale.edu
    !H      It is a set of routines which allow format free input in FORTRAN
    !H
    !H      Jan Labanowski
    !H      jkl@osc.edu
    !H----------------------------------------------------------------------

All the rest should be considered within the latest GPL.

## Build and test

### Requirements
- gfortran
- perl (used by `bin/makemake*`)
- csh (used by `bin/cleanup`, i.e. `make cleanall`)
- python3 (used by the test scripts)
- octave, for `tools/CSMD` (not needed to build or run anything else)

### Build
    cd Freemol
    ./config/configure m_generic_linux gfortran $PWD   # macOS: m_generic_osx
    make freemol

`adfrom` is skipped by design: it needs the commercial ADF libraries, which
aren't part of this repository.

Run any built program with `-i`/`-o` for its input/output files, e.g.
`bin/CSMG.exe -i input.txt -o output.txt` (see [Programs](#programs) below
for each program's input format and example inputs).

### Test scripts
All in `Freemol/tests/`, runnable after a build (see `.github/workflows/ci.yml`
for how CI wires them together):
- `build_all.sh [machine]` -- configure + build (default machine:
  `m_generic_linux`).
- `check_autoconfigure.sh` -- `config/autoconfigure -r` must actually run
  configure (a regression test for a bash-ism that silently no-ops under
  `/bin/sh` on Ubuntu).
- `run_csmg_regression.sh` -- replays every `data/CSMG/tests/*.mld` that has a
  reference output under `data/CSMG/outs/{osx,generated}/` (59 of 59
  inputs, as of this writing) through `bin/CSMG.exe`, and compares CSM
  values (not optimizer intermediates, which are allowed to differ across
  compilers) via `numdiff.py`.
- `run_smoke.sh` -- checks `fit1Dpol`, `ch4sym2cart`, `XY4PolySphere`,
  `XY4Coord` and `Frimol` against fixtures that have no reference output to
  regress against.
- `run_devflow.sh` -- proves the "add a program" workflow below still works
  end to end, then cleans up after itself.
- `check_clean.sh` -- after a build, `make cleanall`/`make cleanconfig` (and
  `.gitignore`) must leave `git status` clean.
- `check_csmd_tools.sh` -- checks `tools/CSMD` (Octave) against the values
  used in the CSMG example below.

## Known issues

- **XY4Coord's `Sr` redundancy solve was broken; fixed.**
  `[x-xy4-symmcoord]`'s 10 values split into two groups by what they
  actually move: `S1, S2x, S2y, S2z` change bond lengths unevenly while
  leaving all 6 angles at the tetrahedral 109.4712 (all four pass);
  `S2a, S2b, S4x, S4y, S4z` leave all 4 bonds at the reference 1.0900
  while perturbing angles. Every linear angle displacement leaves a
  second-order residual (expected, not a defect -- see
  [PRECISION_NOTES.md](PRECISION_NOTES.md#3-xy4coord-a-redundancy-solver-that-solved-for-the-wrong-sr-fixed)
  for the full numerical write-up), which XY4Coord's redundant
  coordinate `Sr` exists to absorb via "Generate Redundancies" ->
  `eval_sr()`. That solver used to fail: it solved a quartic (from a
  paper not present in this repo, `wang-carrington-2003.pdf`) for a
  *cosine-based* `Sr`, while the rest of the program (`get_ra`/`get_symc`)
  uses a *radian-displacement-based* `Sr` -- two different quantities
  sharing one variable, so the physical root (`Sr=0` at equilibrium) was
  never found, and every one of the five angular displacements above
  failed its own self-check.

  Fixed without the paper: `eval_sr` (`XY4coord.F90`) now root-finds `Sr`
  directly against `do_checks()`'s own already-correct closure condition
  (each of its four Gamma-Sum checks, `XY4coord.F90:577-589,610-639`,
  equal to `2*pi` exactly at a valid geometry) via bisection from `Sr=0`,
  instead of the old quartic/`qtcrt` call. Verified end to end against
  the real binary: `equilibrium`, `s2a_only`, `s2b_only`, `s4x_only` and
  two new fixtures (`s4y_only`, `s4z_only`, added since none existed
  before) all now print `Check OK at Sr solution` unmodified, and the
  resulting Cartesian coordinates reproduce every requested bond/angle
  displacement, cross-checked independently in Python. Wired into
  `Freemol/tests/run_smoke.sh` as a real regression test
  (`xy4coord_redundancy_selfcheck`), replacing the previous
  `xy4coord_known_issue` entries for these fixtures.

  Investigating the original bug also found a real, separate one: the
  "H4 test" block in `do_checks()` used the identical condition already
  used for the "H2 test" two blocks above (`acagam(2)+acagam(3)+acagam(12)`
  instead of `acagam(7)+acagam(8)+acagam(10)`, a copy-paste typo -- the
  accompanying debug message right next to it already used the correct
  formula). Fixed in an earlier PR. No test was added for it: an
  exhaustive search (800,000+ sampled `Sda` combinations) found no input
  where it changes `do_checks()`'s overall pass/fail verdict, for this
  symmetric reference geometry.

  **`Sr` alone is a different, expected failure -- not a bug, unaffected
  by the fix above.** Raising all six X-C-X angles together (`Sr`'s own
  totally-symmetric direction) is geometrically impossible for four
  fixed-length bonds from one center: the four bond-direction unit
  vectors' Gram matrix stays positive-semidefinite only while their
  common pairwise cosine stays >= -1/3 (the tetrahedral value); pushing
  it more negative -- exactly what widening every angle does -- has no
  real solution. Confirmed directly: `Sr = 0.03` alone overshoots
  `do_checks()`'s `2*pi` Gamma Sum closure by about 5.2 degrees, and the
  excess scales *linearly* in `Sr`, not quadratically the way `S4x`'s
  does -- a first-order-forbidden direction, categorically different from
  the second-order residual `eval_sr`'s redundancy solve now correctly
  absorbs for the other five displacements.

- **ch4sym2cart: the "Unchenged Coordina..." diagnostic lines for H3/H4
  are geometrically wrong** (H2-C-H3 comes out around 33 degrees instead
  of the tetrahedral 109.47 degrees). `ch4sym2cart.F90` lines ~303-304 and
  ~330-331 pass a literal `120.0_FREAL` straight into `cos()`/`sin()` as
  if it were already radians; everywhere else in the file (e.g. `fa109` at
  line 239) converts degrees via `/180.0*LPI` first. This block is a
  diagnostic-only reconstruction of the reference geometry (it doesn't
  affect `ch4sym2cart`'s actual `NEW Coordinates` output), so it hasn't
  been touched.

## Releases

Versions follow [semver](https://semver.org/) (`vMAJOR.MINOR.PATCH`) as git
tags, starting at `v1.0.0`. Pushing a tag matching `v*.*.*` (or a manual
`workflow_dispatch`) runs `.github/workflows/release.yml`, which builds and
tests the Linux (x86_64) and macOS (arm64) binaries, then -- for an actual
tag push -- publishes a GitHub Release with a `freemol-<tag>-<platform>.tar.gz`
per platform (the built `.exe` programs plus `README.md` and `COPYING`).
`adfrom` is never included: it isn't built without the commercial ADF
libraries. There is no Windows build: no `m_generic_*` config exists for it.

When cutting a release, update `distnum` (and `vernum`/`vernumdist`) in
`Freemol/config/printversion` to match the new tag first.

## MCP tools

`mcp/` is a Python [MCP](https://modelcontextprotocol.io/) server exposing
the coordinate-transformation programs as tools for LLM clients: a thin
subprocess wrapper around the same built binaries CI already regression-
tests, not a reimplementation. Every response includes a citation (exact
routine, file:lines, freemol version, and a GitHub permalink pinned to the
commit that ran) so results trace back to the original Fortran. Currently
wraps `XY4PolySphere` (polyspherical -> Cartesian), `XY4Coord` and
`ch4sym2cart` (symmetric displacement -> Cartesian); `fit1Dpol` and `CSMG`
aren't wrapped yet. Two further tools, `list_freemol_sections` and
`read_freemol_section`, give generic access to the
[sectioned input format](#sectioned-input-files-an-old-ini-style-format-from-before-molden)
itself rather than any one program -- these are a Python mirror of
`osec_set`'s scanning convention (there's no numerics in listing/reading
section text, so no binary is involved, and they work without freemol
even being built). See [`mcp/README.md`](mcp/README.md) to build and run
it.

### Example prompts

Once the server is connected to an MCP client (Claude Desktop, Claude
Code, ...), these are the kind of plain-language requests that route to
the tools above -- each one verified against the real tools, not just
written down:

- *"I just want to build an ideal methane molecule -- one carbon atom in
  the middle, four hydrogen atoms spaced out evenly around it like the
  corners of a tiny pyramid, all the same distance away. Can you give me
  the actual 3D coordinates?"*
  No chemistry or coding background needed for this one. Behind the
  scenes it's the same computation as the next example below --
  `xy4polysphere_to_cartesian` with the bond length and the angles that
  make a perfect tetrahedron (109.4712 and 60 degrees) -- but the person
  asking never has to know those numbers or their names; the AI works
  them out from "evenly spaced, same distance away."

- *"Using freemol, convert a methane-like XY4 molecule with all four
  bonds at 1.10, polar angles 109.4712 degrees, and both azimuthal angles
  at 60 degrees into Cartesian coordinates. Where in the code did that
  computation actually happen?"*
  Calls `xy4polysphere_to_cartesian`; returns 4x1.1000 bonds, 6x109.4712
  degree angles, Metpot4=107.6015 (vs. 9.5797 at the equilibrium
  reference), full XYZ coordinates, and a citation linking straight to
  `poly2cart` in `XY4PolySphere.F90`.

- *"Apply a symmetric stretch of S1=0.08 to methane's reference geometry
  with ch4sym2cart and show me the new hydrogen positions."*
  Calls `ch4sym2cart_apply_displacement`; returns H1=(1.1000, 0, 0) and
  the other three hydrogens at the same 1.1000 distance from the origin,
  tetrahedral angles preserved.

- *"Do that same symmetric stretch with XY4Coord instead, and check the
  two tools agree."*
  Calls `xy4coord_apply_displacement`; the bonds/angles it reports
  (1.1000 / 109.4712, same as above) match `ch4sym2cart`'s result -- a
  useful sanity check across two independently-implemented programs.

- *"Now try a pure angle displacement, S2a=0.03, on XY4Coord instead of a
  bond stretch -- does that work?"*
  Calls `xy4coord_apply_displacement` again; the direct displacement
  alone doesn't close (a second-order residual -- see
  [PRECISION_NOTES.md](PRECISION_NOTES.md#3-xy4coord-a-redundancy-solver-that-solved-for-the-wrong-sr-fixed)),
  so the tool falls back to XY4Coord's own `Sr` redundancy solve and
  returns the geometry it finds, with `sr_used` reporting which `Sr`
  that actually was.

- *"What sections does `Freemol/data/CSMG/tests/G_C1.mld` have, and what's
  in the symmetry-operations one?"*
  Calls `list_freemol_sections` (finds `molecule`, `x-csmg-cgauss`,
  `x-csmg-symop`) then `read_freemol_section` on `x-csmg-symop`; no need
  to know CSMG's own format in advance, or for freemol to even be built.

## Adding a program

The framework lives under `./Freemol`. To add a program, make a directory:

    mkdir programs/mycode

and add `.F90` files in there. A `makemake` utility generates the Makefile
for you; you'll typically also need a `Makelocal` file like this:

    USEDLIBS = -L$(FMLIBDIR) -lmodules -lmoduledata -lutils -lincludes \
               -lfmlapack -lfmblas

Read the command line like this:

    use pcmdline
    !
    ! We read the command line for I/O file names
    !--------------------------------------------
    call pcmd_iargc(nargs)
    allocate(cargs(nargs),STAT=irc)
    !
    ! Get the current arguments
    call pcmd_getarg(nargs,cargs)
    !
    ! Check if there is -i or -o for input or output
    call pcmd_getio(nargs,cargs,finput,foutput)
    !
    ! Open the input
    irc = baseio_open(iunin,finput,stat='OLD')
    !
    ! Check a complex flag: -p noshift,onestep
    irc = pcmd_checkarg(nargs,cargs,'-p',params)
    irc = index(params,'noshift')
    if(irc.gt.0) then
       demo_noshift = .true.
    endif
    !
    irc = pcmd_checkarg(nargs,cargs,'-p',params)
    irc = index(params,'onestep')
    if(irc.gt.0) then
       demo_onestep = .true.
    end if

Or read a molecule section (INI-file-like, molden-inspired) like this:

    ! We need the molecule section
    !-----------------------------
    rewind(iunin)
    irc = osec_set(iunin,'molecule',params)
    if(irc.lt.0) then
       call message(MESERRO,'[demo_rdinput] No molecule section.')
       return
    end if
    !
    call message_comment(iunbin,'Got molecule section. With Parameters:')
    call message_comment(iunbin,adjustl(params))
    !
    call message(MESLOG,"[demo_rdinput] molecule_init called.")
    irc=molecule_init()
    if(irc.lt.0) then
       call message(MESERRO,'[demo_rdinput] Cannot Initialize the molecule.')
       return
    end if
    !
    call molecule_read(iunin,params)
    call molecule_print(iunbin)

`Freemol/tests/run_devflow.sh` builds and runs a throwaway program using
exactly this pattern on every CI run, so this walkthrough stays true.

## Programs

### CSMG
Computes the CSM ("Continuous Symmetry Measure") of a set of weighted
points using Gaussian functions -- e.g. how close a molecule's geometry is
to belonging to a given point group. See the paper
([PDF](https://fabiomariotti.files.wordpress.com/2013/07/csm_paper.pdf))
and a [short blog post](https://blog.techottis.ch/2016/04/01/continuous-symmetry-measure-an-old-work/)
for background. Example input:

    # from symmetry_mod file
    [molecule] nrec=3 format=nscxyz
    o  1 8 1.000000 2.000000 3.000000
    h  2 1 2.000000 2.000000 3.000000
    h  3 1 1.000000 3.000000 3.000000
    [x-csmg-cgauss]
    1 2.883843 1.591247
    2 3.838395 7.701635
    3 3.838395 7.701635
    [x-csmg-symop]
     E    1.0000  0.0000  0.0000    0.0  0.0000  0.0000  0.0000  F
     C2_2   0.707107 0.707107 0.000000 180.000000  0.0000  0.0000  0.0000  F
     S1   -0.707107 0.707107 0.000000 180.0000  0.0000  0.0000  0.0000  I
     S2   0.000000 0.000000 1.000000 180.0000  0.0000  0.0000  0.0000  I

The `x-csmg-cgauss` values (Gaussian exponent and normalization per atom)
come from `tools/CSMD`'s `genGTOcoeff` -- see Tools below;
`check_csmd_tools.sh` regenerates exactly the values used here.

### XY4Coord
Converts symmetry-adapted internal displacement coordinates (`S1`,
`S2x/y/z`, `S2a/b`, `S4x/y/z`, `Sr`) to Cartesian for XY4-type
(methane-like, tetrahedral AB4) molecules -- the kind of displacement a
vibrational/normal-mode or potential-energy-surface calculation moves
along -- checking its own result against the requested displacement
before generating the redundant `Sr` solutions. See
[Known issues](#known-issues) for the one displacement direction that's
expected to fail by design (raising `Sr` alone). Example input (moving
`S1` by 0.08, all other displacements zero):

    [molecule] nrec=5 format=nscxyz
    c  1 6 0.0 0.0 0.0
    h  2 1  0.62931179  0.62931179  0.62931179
    h  3 1  0.62931179 -0.62931179 -0.62931179
    h  4 1 -0.62931179  0.62931179 -0.62931179
    h  5 1 -0.62931179 -0.62931179  0.62931179
    [x-xy4-symmcoord]
    0.08 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0

### XY4PolySphere
The polyspherical-coordinate counterpart to XY4Coord: converts a
polyspherical description (4 bond lengths, 3 polar and 2 azimuthal
angles) to Cartesian for the same XY4 molecule family. Polyspherical
coordinates are the representation quantum-dynamics calculations tend to
prefer (their kinetic energy operator separates more cleanly than in
Cartesian or simple internal coordinates), so this and XY4Coord cover the
same molecules through two different coordinate conventions rather than
one replacing the other. Can also sample a Metpot4 potential over ranges
of those coordinates (`x-xy4-genrandom`), with `-R <seed>` for a
reproducible sampling run. Example input (the reference/equilibrium
geometry -- 4 equal bonds, all six angles at the tetrahedral 109.4712):

    [molecule] nrec=5 format=nscxyz
    c  1 6 0.0 0.0 0.0
    h  2 1  0.62931179  0.62931179  0.62931179
    h  3 1  0.62931179 -0.62931179 -0.62931179
    h  4 1 -0.62931179  0.62931179 -0.62931179
    h  5 1 -0.62931179 -0.62931179  0.62931179
    [x-xy4-polyspherical]
    1.10 1.10 1.10 1.10 109.4712 109.4712 109.4712 60.0 60.0

### ch4sym2cart
Converts symmetric-coordinate displacements to Cartesian, specifically for
CH4 (methane)-like molecules -- an earlier, CH4-only counterpart to
XY4Coord's more general XY4 treatment. Despite sharing the same field
names and order, its radial/angular grouping of those fields is *not* the
same as XY4Coord's (see [Known issues](#known-issues)), so inputs aren't
interchangeable between the two. Example input (same displacement as the
XY4Coord example above -- note it moves a *different* geometric feature
here, per that known convention mismatch):

    [molecule] nrec=5 format=nscxyz
    c  1 6 0.0 0.0 0.0
    h  2 1  0.62931179  0.62931179  0.62931179
    h  3 1  0.62931179 -0.62931179 -0.62931179
    h  4 1 -0.62931179  0.62931179 -0.62931179
    h  5 1 -0.62931179 -0.62931179  0.62931179
    [x-ch4-symmcoord]
    0.08 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0

### fit1Dpol
Fits data in 1 dimension using Minuit, with a polynomial setup built in
(add your own function for anything else) -- e.g. fitting a scanned
potential-energy curve to a polynomial. Works with multi-column ASCII
files, and can script which data to fit with ranged commands like copy,
add, sub, shift, scale, etc. Bundles its own copy of Minuit
(`fit1Dpol_minuit.F90`), not the shared one under Libraries below -- see
Libraries/Modules for why. Example input (fitting `y = x^2` with a
degree-2 polynomial):

    [x-fit1dpol-data] numrec=5 numcol=2
    1.0  1.0
    2.0  4.0
    3.0  9.0
    4.0 16.0
    5.0 25.0
    [x-fit1dpol-pars] numpars=3 funct=polynom eps=1.0D-14
    #name idx start_value init_step  Fixed
    f0    1   0.1         0.01         N
    f1    2   0.0         0.0001       Y
    f2    3   2.2         0.001        N

### Freemol
Just a placeholder at present (`Frimol.exe` prints "Not Yet Read" and
exits) -- no framework demo program has been built out here yet. Still
useful as-is: `run_smoke.sh` builds and runs it on every CI run as a
minimal end-to-end check that the build plumbing itself works,
independent of any real program logic.

### adfrom
Reads ADF (Amsterdam Density Functional)'s binary `Tape21` output and
produces an ASCII-readable format. Needs the commercial ADF libraries,
which aren't part of this repository, so it isn't built by `make freemol`
(see Build above).

## Framework

- **configure / autoconfigure**: `config/configure <machine> <compiler>
  <FMHOME>` sets up `Makeflags`/`Make.machine` for a given machine
  (`m_generic_linux`, `m_generic_osx`, ...) and compiler.
  `config/autoconfigure -r` auto-detects both and runs `configure` for you.
- **makemake**: `bin/makemake*` (perl) generates each directory's Makefile
  from its `Makelocal`, so adding a program is "add files + Makelocal",
  not "hand-write a Makefile".
- **help extraction**: comment lines starting `!H` are extracted from each
  program's sources into `docs/<program>.exe.txt` during its build (see
  `Freemol/help/Makefile` for where that fits into `make freemol`).
- **manuals**: longer-form LaTeX documentation lives under `manuals/` (see
  e.g. `manuals/CSMD`).

## Libraries
These libraries and tools are used by the programs above. Each item below
has a short line showing its shape; for a fuller working example, check
the programs (particularly "Adding a program" above), and check the
generated Makefiles for how to link against them.

### Includes
- **linetools**: format-free input (adapted from CCL's FREEREAD, see
  Licences above). E.g. `call line_read(line, iwrd, rwrd, swrd, niwrd,
  nrwrd, nswrd, form)` splits one free-form input line into typed
  integer/real/string words.
- **messages**: a logging facility. E.g. `call message_value(MESWARN,
  'current too low:', current)`.
- **pcmdline**: command-line handling. E.g. `call pcmd_getio(nargs,
  cargs, finput, foutput)` pulls `-i`/`-o` filenames out of argv (see
  "Adding a program" above).
- **strtools**: string manipulation routines. E.g. `call
  str_upcase(string)`.
- **vartypes**: kind-parametrised variable types (see "Careful, portable,
  reproducible double precision" above). E.g. `real(FREAL) :: x`.
- **chemconst**: constants defined at different precisions. E.g.
  `real(FREAL), parameter :: bohr_to_angstrom = cc_au2ang`.

### Utils
- **mathtools**: a few simple, common matrix tools. E.g. `call
  mathtools_masscenter(xyz, we, rvec)` returns the (optionally
  weighted) center of a set of Cartesian points.
- **qnumbers**: a small tool for quantum numbers. E.g. `call qn_ms(m,
  msign, idmn)` fills the `m`/sign(`m`) pairs for angular momentum
  `idmn`.

### Modules
- **baseio**: a file manager. E.g. `irc = baseio_open(iunin, finput,
  stat='OLD')` (see "Adding a program" above).
- **extio**: extensions to baseio for easy file manipulation. E.g. `irc =
  extio_open(lefh, file, status, frmt, acc)` opens a file and returns a
  typed `extfile` handle -- see the module's own commented-out
  `extio_openefh` helper for the intended wrapping pattern.
- **osec / sections**: a tool to work with sectioned files, INI-style /
  molden-lineage (see [Sectioned input files](#sectioned-input-files-an-old-ini-style-format-from-before-molden)
  above). `osec` is the layer every program actually uses (`irc =
  osec_set(iunin, 'molecule', params)`, see "Adding a program" above),
  built directly on `linetools`/`strtools`, not on `sections`. `sections`
  (`sections_init`, `section_openfile`, `section_next`) is a separate,
  higher-level file-open API for the same idea that was never finished --
  `section_openfile` currently always returns "not implemented".
- **minuit**: a multidimensional minimization tool, a CERN-derived
  F77-to-F90 port. Nothing currently links against this shared copy --
  CSMG and fit1Dpol each bundle their own local copy instead
  (`csm_minuit.F90`, `fit1Dpol_minuit.F90`), same historical-code
  lineage, just not wired up as a shared dependency.

### Moduledata
- **molecule**: a molden-based structure and tools to store molecule data
  and read it in "free" format. E.g. `irc = molecule_init(); call
  molecule_read(iunin, params); call molecule_print(iunbin)` (see
  "Adding a program" above).

### Tools
- **CSMD** (`tools/CSMD`, Octave): supporting tools for CSMG, notably
  `genGTOcoeff` for the Gaussian exponent/normalization values used in
  `x-csmg-cgauss` sections

### others
- local BLAS and LAPACK (see "Why so complicated?" above for why they're
  vendored)
