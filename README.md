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

- **XY4Coord: every angular displacement fails its own self-checks.**
  `[x-xy4-symmcoord]`'s 10 values split into two groups by what they
  actually move: `S1, S2x, S2y, S2z` change bond lengths unevenly while
  leaving all 6 angles at the tetrahedral 109.4712 (all four pass);
  `S2a, S2b, S4x, S4y, S4z, Sr` leave all 4 bonds at the reference 1.0900
  while perturbing angles -- and every one of those six currently fails,
  confirmed individually for all six, not just `S2a` as originally found.
  Two distinct failure paths:
  - `S2a`, `S2b` fail the Cartesian self-check at `XY4coord.F90:756-763`
    ("Error in Cartesian routine"). Both show the *same specific*
    signature: of the 6 computed angles, 5 match exactly, and only
    `ryxy(6)` (the angle between H3 and H4, atoms `mxyz(:,4)` and
    `mxyz(:,5)`) is off. That repeatability across two different inputs
    reinforces the same diagnosis: the ad-hoc sign disambiguation for
    H4's position (`XY4coord.F90:697-736`, the block that ends with the
    code's own `"Assumed signed 'sin'. It can be inconsistent with input
    sym data."` warning at line 747), not a wholesale failure of the
    angle path.
  - `S4x`, `S4y`, `S4z`, `Sr` fail an *earlier* check, `do_checks()`'s
    "Gamma Sum bigger than 180 for Y atom N" validation
    (`XY4coord.F90:611-639`), before `get_cart` is even reached. Root
    cause now understood, evidence-backed: each Gamma Sum must equal
    exactly 2*pi for a mathematically valid vertex closure, and the
    check compares against `2.0*LPI` with **zero tolerance**. Two
    compounding effects push real inputs past that exact boundary: (1)
    a few ULPs of floating-point noise in the cos -> divide -> acos ->
    sum chain, present even at equilibrium; (2) a genuine,
    displacement-*squared* defect inherent to representing angle
    changes with a **linear** symmetry-coordinate formula (`get_ra`'s
    `da(1:6)`, `XY4coord.F90:390-395`) -- confirmed empirically:
    excess-over-2*pi for a range of `S4x` values fits `~38*S4x^2`
    degrees almost exactly (37.4-39.7 across a 30x range in `S4x`),
    vanishing into the floating-point noise floor once `S4x` drops
    below ~0.001. `S1/S2x/S2y/S2z` never trip this because their `da()`
    contribution is always zero -- angles never move at all, so there's
    no defect to accumulate. A real fix would need a deliberately-chosen
    tolerance (a design decision -- how large a displacement should the
    linear approximation still be trusted for -- not a one-line
    correction), so it hasn't been attempted here. See
    [PRECISION_NOTES.md](PRECISION_NOTES.md#3-xy4coord-do_checks-a-zero-tolerance-check-on-a-value-that-must-be-exact-undermined-by-a-real-quadratic-defect)
    for the full numerical write-up.

    Investigating this also found a real, separate bug: the "H4 test"
    block used the identical condition already used for the "H2 test"
    two blocks above (`acagam(2)+acagam(3)+acagam(12)` instead of
    `acagam(7)+acagam(8)+acagam(10)`, a copy-paste typo -- the
    accompanying debug message right next to it already used the
    correct formula). Fixed. No test was added: an exhaustive search
    (800,000+ sampled `Sda` combinations, wide and boundary-focused,
    every single displacement dimension alone) found no input where
    this specific bug changes `do_checks()`'s overall pass/fail verdict
    -- for this symmetric reference geometry, Y2's and Y4's Gamma Sums
    take different numeric values but always land on the same side of
    the 2*pi boundary, so Y1/Y2/Y3 already independently catch anything
    Y4 alone would have. Fixed anyway on correctness grounds (it matches
    the debug message beside it and the pattern of the other three
    blocks); flagging the missing test explicitly rather than skipping
    it silently.

  Not fixed: the `S2a`/`S2b` `get_cart` root cause is plausible but not
  proven, and the ~25-year-old trigonometric derivation in that block
  would need a real re-derivation to fix with confidence rather than a
  guess; the `do_checks()` zero-tolerance issue needs a deliberate
  tolerance choice, not attempted here.

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
  Calls `xy4coord_apply_displacement` again; this is the documented Known
  Issue above, and the tool returns a clear error explaining the
  self-check failed, rather than silently returning wrong coordinates.

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
[Known issues](#known-issues) for the displacement paths that aren't
fully correct yet. Example input (moving `S1` by 0.08, all other
displacements zero):

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
- **osec / sections**: a tool to work with sectioned files,
  molden-format-like. `osec` is the friendly layer every program actually
  uses (`irc = osec_set(iunin, 'molecule', params)`, see "Adding a
  program" above); `sections` (`sections_init`, `section_openfile`,
  `section_next`) is the lower-level reader it's built on.
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
