# Precision notes

This file collects, in one place, what testing this codebase has actually
taught us about floating-point behaviour: real compiler-dependent noise
found in CI, and one genuine numerical defect root-caused down to a
formula. It's a companion to two things that already exist but don't fully
cover this:

- README's [Careful, portable, reproducible double precision](README.md#careful-portable-reproducible-double-precision)
  section, which documents the *design* (kind-parametrised reals,
  machine-aware epsilons, a 44-digit pi) -- intent, not evidence.
- README's [Known issues](README.md#known-issues) section, which documents
  *correctness* bugs (wrong results). Some of those bugs turned out to have
  a precision angle too; where that's true, this file has the numerical
  detail and README has the short version plus a pointer here.

Nothing below is news to whoever fixed it at the time -- this is a
deliberate second write-up, collected after the fact, so a future
investigation into "why did CI fail on this one platform" or "why does
this check reject a seemingly-fine input" has one place to start instead
of re-deriving it.

## 1. CSMG regression: what's allowed to differ, and what isn't

`Freemol/tests/run_csmg_regression.sh` compares `CSMG.exe`'s output
against stored references via `Freemol/tests/numdiff.py`. A full-text diff
isn't usable here: CSMG's Minuit-derived simplex optimizer takes a
legitimately different path to the same answer depending on
compiler/platform (e.g. flat directions in the fit that the simplex drifts
along with no effect on the result), so `numdiff.py` only ever compares
lines matching `CSM|Minimization at step` -- the CSM value itself, and the
optimizer's step count and `Value:` field, not its raw intermediate
parameters. That filtering is what makes the reproducibility claim in
README possible at all: 8 of 59 CSMG cases match reference output
generated on a different machine, OS and decade (macOS, 2016) to the last
printed digit, specifically because everything that's allowed to be
noise has already been filtered out before comparison.

## 2. numdiff.py: the glued `Tol(Act/Req)` token

Found while extending CSMG regression coverage from 8/59 to 59/59 inputs
(the other 51 references were freshly generated from current code on
macOS). Several of those 51 cases initially failed on Linux CI only --
not locally, not on the macOS CI job -- despite the CSM values matching.

The cause: `Minimization at step` lines print an achieved-tolerance
diagnostic glued directly onto the previous number with **no whitespace**,
e.g. `0.0000Tol(Act/Req):0.6835E-05/...`. `numdiff.py` tokenizes on
whitespace, so this glued blob was one token, failed to parse as a number,
and fell through to exact-string comparison instead of the intended
numeric tolerance comparison -- silently defeating the tolerance logic for
that one token on every matched line. The tolerance-achieved value is
optimizer-internal state and genuinely differs in its last couple of
digits across compilers, same as everything else Section 1 describes; it
had just never been exercised by a Linux build before, because the
generated baseline was captured on the same macOS machine that was
re-comparing against itself locally.

Fix (`Freemol/tests/numdiff.py`): truncate each matched line at `Tol(`
before comparing, so the step count and `Value:` field are still checked
but the noisy suffix isn't. Verified the fix still catches real
divergences with a synthetic mismatched-step-number test before trusting
it.

**Lesson**: a Fortran `write` statement with no explicit field separator
between two values can produce a token that looks numeric to a human but
isn't whitespace-delimited for a naive parser -- and a tolerant-comparison
tool that falls back to exact-string matching on unparseable tokens will
silently stop being tolerant for exactly that token, with no error, just
occasional platform-dependent failures.

## 3. XY4Coord `do_checks()`: a zero-tolerance check on a value that must be exact, undermined by a real quadratic defect

The deepest root-cause investigation in this codebase so far. Full
narrative and the file:line references are in README's Known issues
entry; this is the numerical mechanism.

`do_checks()` (`Freemol/programs/XY4Coord/XY4coord.F90:611-639`) checks
that four per-vertex "Gamma Sum" values each equal `2*LPI` (2*pi) --
a geometric closure identity that is true by construction for a valid
tetrahedral vertex. The check is `.gt.(2.0_FREAL*LPI)`: **zero tolerance**,
not even one ULP. Two effects independently push a real computed value
past that exact boundary:

1. **Floating-point noise**, present even at equilibrium (zero
   displacement): the cos -> divide -> acos -> sum chain that builds each
   Gamma Sum accumulates a few ULPs of rounding error before comparing
   against the boundary.
2. **A genuine displacement-squared defect**, independent of noise: angle
   displacements are parametrised through `get_ra`'s `da(1:6)`
   (`XY4coord.F90:390-395`), a **linear** symmetry-coordinate formula. A
   linear approximation to an inherently nonlinear (angular) quantity has
   error that grows with the *square* of the displacement, not the
   displacement itself. Measured directly: the excess of a Gamma Sum over
   `2*pi`, as a function of the `S4x` displacement, fits `~38 * S4x^2`
   degrees closely (37.4-39.7 across a 30x range of `S4x` values from
   0.03 down to 0.001), and only drops into the floating-point noise floor
   of effect (1) below `S4x ~= 0.001`.

Radial displacements (`S1, S2x, S2y, S2z`) never trigger this at all --
their contribution to `da()` is structurally zero, so angles never move
and there's nothing for the quadratic defect to act on. This is why the
known-issue split in README is exactly radial-vs-angular, not
case-by-case.

A real fix needs a deliberately chosen tolerance (how large a displacement
should the linear approximation still be trusted for is a design
decision, not a one-line correction), so none has been attempted.

### A related but separate bug, found along the way

While tracing this, the "H4 test" block turned out to reuse the "H2
test"'s condition verbatim (`acagam(2)+acagam(3)+acagam(12)` where it
should have been `acagam(7)+acagam(8)+acagam(10)`) -- an ordinary
copy-paste typo, not a precision issue, fixed in commit `6f75ef8`. Worth
recording here anyway because of how it was verified safe to fix without
a regression test: an exhaustive search over 800,000+ sampled `Sda`
combinations (wide-range and boundary-focused, every displacement
dimension in isolation and combined) found no input where the bug changes
`do_checks()`'s overall pass/fail verdict, because for this reference
geometry Y2's and Y4's Gamma Sums always land on the same side of the
`2*pi` boundary even though their numeric values differ. Fixed anyway on
correctness grounds; the absence of a test is deliberate and explained in
the commit message, not an oversight.

## Lessons for future numerical debugging here

- **Test on more than one compiler/platform before trusting a passing
  local run.** Every noise-related finding above (Sections 1 and 2) was
  invisible on the machine that generated the reference output, and only
  showed up once CI ran the same comparison on a different gfortran
  version. This is the actual reason the CI matrix builds Linux
  (gfortran 12/13/14) and macOS separately rather than just one of them.
- **Watch for glued, non-whitespace-separated tokens in Fortran text
  output** when writing a comparison tool against it -- a naive
  whitespace tokenizer will silently stop being tolerant for a token like
  that (Section 2).
- **A check against a value that's mathematically exact by construction
  (sums to exactly `2*pi`, etc.) needs a deliberately-chosen tolerance,
  never zero.** Zero tolerance doesn't mean "exact" in floating point; it
  means "will eventually fail on some legitimate input" (Section 3).
- **A linear approximation to a nonlinear quantity has error that grows
  with displacement squared**, not linearly -- "the displacement is small"
  is not by itself a reason to expect a downstream zero-tolerance check to
  pass (Section 3).
