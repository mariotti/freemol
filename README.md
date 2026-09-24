# freemol

[![CI](https://github.com/mariotti/freemol/actions/workflows/ci.yml/badge.svg)](https://github.com/mariotti/freemol/actions/workflows/ci.yml)

A Fortran 90 toolkit for molecular geometry and symmetry: computing the
Continuous Symmetry Measure (CSM) of a set of weighted points, converting
between Cartesian and internal/polyspherical coordinates for XY4-type
(methane-like) molecules, and fitting 1D data with Minuit. Started in 2003
(as "Freemol2000"), moved through CVS in 2009 and onto GitHub in 2015 --
the [`Legacy-2016`](https://github.com/mariotti/freemol/tree/Legacy-2016) tag marks that untouched
historical state. Since 2026 it also has a CI-verified build and tagged
releases; see [Build and test](#build-and-test) and
[Releases](#releases) below.

## Careful, portable, reproducible double precision

This is old code, but precision was taken seriously from the start, and
that's worth calling out rather than losing in a pile of `real*8`:

- **Kind-parametrised types.** `Freemol/includes/vartypes.F90` defines
  `FREAL = kind(1.0d0)`, and every real literal in the codebase is written
  `..._FREAL` rather than a bare double. Changing precision is a one-line
  decision, not a search-and-replace across the source. (`SREAL` and
  `BREAL` are defined alongside `FREAL` for future differentiation --
  today all three happen to resolve to the same double precision.)
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
  compares CSMG's CSM value (printed to 11 decimal digits) against
  reference output generated on macOS in 2016, on every CI run --
  currently Linux with gfortran 12/13/14 and macOS with current Homebrew
  gfortran. It matches to the last printed digit, across a decade,
  two operating systems and three compiler versions. That's not a claim,
  it's a passing test you can rerun yourself. (See
  [Known issues](#known-issues) below for the places numerics *aren't*
  solid yet -- being honest about those is part of the same discipline.)

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

### Test scripts
All in `Freemol/tests/`, runnable after a build (see `.github/workflows/ci.yml`
for how CI wires them together):
- `build_all.sh [machine]` -- configure + build (default machine:
  `m_generic_linux`).
- `check_autoconfigure.sh` -- `config/autoconfigure -r` must actually run
  configure (a regression test for a bash-ism that silently no-ops under
  `/bin/sh` on Ubuntu).
- `run_csmg_regression.sh` -- replays every `data/CSMG/tests/*.mld` that has a
  reference output under `data/CSMG/outs/osx/` through `bin/CSMG.exe`, and
  compares CSM values (not optimizer intermediates, which are allowed to
  differ across compilers) via `numdiff.py`.
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

- **XY4Coord: pure angle displacement fails its own Cartesian self-check.**
  With only `S2a` set (e.g. `0.0 0.0 0.0 0.0 0.03 0.0 0.0 0.0 0.0 0.0` in
  `[x-xy4-symmcoord]`, no bond-length change), the self-check at
  `XY4coord.F90:756-763` reports "Error in Cartesian routine". The
  evidence points at one specific angle: of the 6 computed angles, 5
  match the expected value exactly, and only `ryxy(6)` (the angle between
  H3 and H4, i.e. atoms `mxyz(:,4)` and `mxyz(:,5)`) is off (~112.5 deg
  computed vs ~110.5 deg expected). H4's angles to H1 and H2 are both
  correct, which points specifically at the ad-hoc sign disambiguation
  for H4's position (`XY4coord.F90:697-736`, the block that ends with the
  code's own `"Assumed signed 'sin'. It can be inconsistent with input
  sym data."` warning at line 747) rather than a wholesale failure of the
  angle path. Not fixed: the root cause is plausible but not proven, and
  the ~25-year-old trigonometric derivation in that block would need a
  real re-derivation to fix with confidence rather than a guess.

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
aren't wrapped yet. See [`mcp/README.md`](mcp/README.md) to build and run
it.

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
Converts symmetry-adapted internal displacement coordinates (`S1`, `S2x/y/z`,
`S2a/b`, `S4x/y/z`, `Sr`) to Cartesian for XY4-type (methane-like)
molecules, checking its own result against the requested displacement
before generating the redundant `Sr` solutions. See
[Known issues](#known-issues) for the one displacement path that isn't
fully correct yet.

### XY4PolySphere
The polyspherical-coordinate counterpart to XY4Coord: converts a
polyspherical description (4 bond lengths, 3 polar and 2 azimuthal
angles) to Cartesian, and can also sample a Metpot4 potential over
ranges of those coordinates (`x-xy4-genrandom`), with `-R <seed>` for a
reproducible sampling run.

### ch4sym2cart
Converts symmetric-coordinate displacements to Cartesian, specifically for
CH4 (methane)-like molecules.

### fit1Dpol
Fits data in 1 dimension using Minuit, with a polynomial setup built in
(add your own function for anything else). Works with multi-column ASCII
files, and can script which data to fit with ranged commands like copy,
add, sub, shift, scale, etc.

### Freemol
Just a placeholder at present (`Frimol.exe` prints "Not Yet Read" and
exits).

### adfrom
Works with ADF (Amsterdam Density Functional). Needs the commercial ADF
libraries, which aren't part of this repository, so it isn't built by
`make freemol` (see Build above).

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
These libraries and tools are used by the programs above. For a working
example of their usage, check the programs; to link against them, check
the generated Makefiles.

### Includes
- **linetools**: format-free input (adapted from CCL's FREEREAD, see
  Licences above)
- **messages**: a logging facility
- **pcmdline**: command-line handling
- **strtools**: string manipulation routines
- **vartypes**: kind-parametrised variable types (see "Careful, portable,
  reproducible double precision" above)
- **chemconst**: constants defined at different precisions

### Utils
- **mathtools**: a few simple, common matrix tools
- **qnumbers**: a small tool for quantum numbers

### Modules
- **baseio**: a file manager
- **extio**: extensions to baseio for easy file manipulation
- **osec / sections**: a tool to work with sectioned files, molden-format-like
- **minuit**: a multidimensional minimization tool

### Moduledata
- **molecule**: a molden-based structure and tools to store molecule data
  and read it in "free" format

### Tools
- **CSMD** (`tools/CSMD`, Octave): supporting tools for CSMG, notably
  `genGTOcoeff` for the Gaussian exponent/normalization values used in
  `x-csmg-cgauss` sections

### others
- local BLAS and LAPACK (see "Why so complicated?" above for why they're
  vendored)
