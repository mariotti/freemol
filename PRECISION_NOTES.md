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

## 3. XY4Coord: a redundancy solver that solves for the wrong Sr

The deepest root-cause investigation in this codebase so far, and one
that needed correcting once already: an earlier version of this section
diagnosed the angular-displacement known issue (README's Known issues)
as a "zero-tolerance boundary check defeated by a linear
symmetry-coordinate defect," and recommended a tolerance change as the
fix. That was wrong on both counts -- found by trying to actually apply
the fix it suggested. This is the corrected mechanism, reverified against
the built binary and the actual source.

### The quadratic residual is expected, not a defect

`XY4Coord`'s angle displacements (`S2a, S2b, S4x, S4y, S4z`) are
parametrised through `get_ra`'s `da(1:6)` (`XY4coord.F90:390-395`), a
**linear** symmetry-coordinate formula. Any linear parametrisation of an
inherently nonlinear (angular) quantity leaves a second-order residual --
that's standard, not a bug, and it's exactly what a redundant coordinate
is *for*. XY4Coord has one: after computing a displacement, it calls
`get_symc` -> `eval_sr` (`XY4coord.F90:786-828`, `1361-1417`), which
solves a quartic (`qtcrt`, coefficients "from JCP 118 (2003) 6260 -
Wang, Carrington") for the redundant coordinate `Sr` (`Sda(6)`), then
retries the displacement with each root substituted in. That machinery
exists precisely to absorb the second-order inconsistency a linear angle
parametrisation produces -- the `~38 * S4x^2` degrees of excess-over-
`2*pi` in `do_checks()` this section used to call a "defect" (still
correctly measured: 37.4-39.7 across a 30x range of `S4x`, vanishing into
floating-point noise below `S4x ~= 0.001`) is exactly the residual that
solver is supposed to correct.

### The solver doesn't work, because it solves for a different Sr

Verified directly (`bin/XY4coord.exe -i data/XY4Coord/tests/
equilibrium.inp`): at equilibrium every requested displacement is zero,
`Sr = 0` is exactly correct (the main run prints `Check OK at Sr input
value: 0`) -- and yet `eval_sr`'s own quartic is not zero at `Sr = 0`.
With all `Sda(1:5) = 0` it collapses to
`A(x) = -x^4/12 + (2*sqrt(6)/9)*x^3 - x^2 + 1`, and `A(0) = 1`: `Sr = 0`
is never a root. The roots the program actually prints are `-0.8165` and
`2.4495` (a numerical triple root). Both have a clean closed form,
checked in Python against the real coefficients in `eval_sr`
(`XY4coord.F90:1389-1410`):
- `-0.8165 = sum(cos(109.4712 deg))/sqrt(6)` for the six *equilibrium*
  bond angles -- confirmed to solve `A(x) = 0` to `1e-15`.
- `2.4495 = sqrt(6)`, the degenerate case where every pairwise angle has
  collapsed to 0 (every cosine = 1).

That is: `eval_sr`'s `Sr` is a symmetric combination of angle *cosines*
-- the Wang-Carrington paper's own coordinate, matching its citation --
while `get_ra`/`get_symc`'s `Sr` (`Sda(6)`, `XY4coord.F90:812`) is a
symmetric combination of *radian angle displacements from equilibrium*.
Two different physical quantities share one variable name and array slot,
and a root computed under the cosine definition gets assigned straight
into the radian-displacement one (`XY4coord.F90:295`,
`Sda(6)=ZSda(isZ)`) -- so the root the program actually needs is never
among the candidates it tries. `get_symc`'s own last line agrees this was
never finished: `call message(MESERRO,"TODO: perform an internal check!
redundant should be consistent!")`.

Consequence, confirmed for every angular fixture that has an XY4Coord
input (`s2a_only`, `s2b_only`, `s4x_only`): none of the roots `eval_sr`
returns pass the self-check either. For `s2a_only` (`S2a = 0.03`), the
Cartesian-reconstructed H3-C-H4 angle at `Sr = 0` is `1.9629` rad against
the `1.9280` rad requested -- `0.035` rad (2 degrees) off. That's a real
geometric inconsistency, not floating-point noise (which would show up at
the `1e-8`-`1e-10` scale), so **loosening `do_checks()`'s tolerance is
not the fix** -- it would just accept structures that don't actually
match the requested displacement. The fix is reconciling `eval_sr`'s
coordinate definition with `get_ra`/`get_symc`'s, against the cited
paper -- a real re-derivation, not attempted here.

(The `S2a`/`S2b` `get_cart` self-check failure -- "Error in Cartesian
routine," only `ryxy(6)` off, previously diagnosed as an H4
sign-disambiguation bug on its own -- is consistent with this same root
cause: at `Sr = 0`, the requested angle set for those inputs is only
realizable with a nonzero `Sr` correction the solver can't supply. The
sign-disambiguation block (`XY4coord.F90:697-736`) stays a plausible
*secondary* contributor, not the diagnosis.)

### `Sr` alone is a different, expected failure -- not a bug

Raising all six X-C-X angles together (`Sr`'s own totally-symmetric
direction, tested directly: all `Sda = 0` except `Sr = 0.03`) is
geometrically impossible for four fixed-length bonds from one center,
independent of any solver bug: the four bond-direction unit vectors'
Gram matrix (`(1-c)*I + c*J` for common pairwise cosine `c`) stays
positive-semidefinite only for `c >= -1/3`, the tetrahedral value itself
-- pushing `c` more negative, exactly what widening every angle does, has
no real solution for *any* nonzero `Sr` in that direction. This shows up
sharply in `do_checks()`'s own numbers: `Sr = 0.03` overshoots the
`2*pi` Gamma Sum closure by about 5.2 degrees, and the excess scales
*linearly* in `Sr` (checked at `0.0001`, `0.001`, `0.03`) rather than
quadratically the way `S4x`'s does -- confirming a first-order-forbidden
direction, categorically different from the second-order residual
`S2a/S2b/S4x/S4y/S4z` leave for the (broken) redundancy solver to absorb.
`Sr` failing on its own is expected by design, not a known issue.

### A related but separate bug, found along the way

While tracing this, the "H4 test" block turned out to reuse the "H2
test"'s condition verbatim (`acagam(2)+acagam(3)+acagam(12)` where it
should have been `acagam(7)+acagam(8)+acagam(10)`) -- an ordinary
copy-paste typo, not a precision issue, fixed in commit `0a0f687`. Worth
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
- **An expected second-order effect may point to a missing/broken
  compensation step rather than a precision defect.** A linear
  coordinate parametrisation leaving a quadratic-in-displacement residual
  is normal, not a flaw; the question worth asking first is whether the
  codebase already has a mechanism meant to absorb that residual (here,
  XY4Coord's `eval_sr` redundancy solver) before reaching for a tolerance
  change (Section 3).
