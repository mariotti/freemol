MODULE constants_NSWC
! Contains the NSWC functions IPMPAR, SPMPAR, DPMPAR, EPSLN, DEPSLN,
! EXPARG & DXPARG
!-----------------------------------------------------------------------
!     WRITTEN using F90 intrinsics by
!        Alan Miller
!        amiller @ bigpond.net.au
!     Latest revision - 1 February 1997
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)

CONTAINS

FUNCTION ipmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
!     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
!     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...

!  INTEGERS.

!     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM

!               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )

!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.

!     IPMPAR(1) = A, THE BASE (radix).

!     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS (digits).

!     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE (huge).

!  FLOATING-POINT NUMBERS.

!     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
!     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
!     NONZERO NUMBERS ARE REPRESENTED IN THE FORM

!               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)

!               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
!               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.

!     IPMPAR(4) = B, THE BASE.

!  SINGLE-PRECISION

!     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.

!     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.

!     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.

!  DOUBLE-PRECISION

!     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.

!     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.

!     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.

!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
INTEGER             :: fn_val

SELECT CASE(i)
  CASE( 1)
    fn_val = RADIX(i)
  CASE( 2)
    fn_val = DIGITS(i)
  CASE( 3)
    fn_val = HUGE(i)
  CASE( 4)
    fn_val = RADIX(1.0)
  CASE( 5)
    fn_val = DIGITS(1.0)
  CASE( 6)
    fn_val = MINEXPONENT(1.0)
  CASE( 7)
    fn_val = MAXEXPONENT(1.0)
  CASE( 8)
    fn_val = DIGITS(1.0D0)
  CASE( 9)
    fn_val = MINEXPONENT(1.0D0)
  CASE(10)
    fn_val = MAXEXPONENT(1.0D0)
  CASE DEFAULT
    RETURN
END SELECT

RETURN
END FUNCTION ipmpar



FUNCTION spmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
REAL                :: fn_val

! Local variable
REAL                :: one = 1.0

SELECT CASE (i)
  CASE (1)
    fn_val = EPSILON(one)
  CASE (2)
    fn_val = TINY(one)
  CASE (3)
    fn_val = HUGE(one)
END SELECT

RETURN
END FUNCTION spmpar



FUNCTION dpmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     DPMPAR PROVIDES THE DOUBLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     DOUBLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        DPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        DPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
REAL (dp)           :: fn_val

! Local variable
REAL (dp)    :: one = 1._dp

SELECT CASE (i)
  CASE (1)
    fn_val = EPSILON(one)
  CASE (2)
    fn_val = TINY(one)
  CASE (3)
    fn_val = HUGE(one)
END SELECT

RETURN
END FUNCTION dpmpar


FUNCTION epsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.0 + EPS .GT. 1.0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
IMPLICIT NONE
REAL                :: fn_val

! Local variable
REAL                :: one = 1.0

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION epsln


FUNCTION exparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     EXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: l
REAL                :: fn_val

! Local variable
REAL                :: one = 1.0

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION exparg


FUNCTION depsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.D0 + EPS .GT. 1.D0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
IMPLICIT NONE
REAL (dp)           :: fn_val

! Local variable
REAL (dp)    :: one = 1._dp

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION depsln


FUNCTION dxparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  DXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     DEXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  DXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF DEXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR DXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: l
REAL (dp)           :: fn_val

! Local variable
REAL (dp)    :: one = 1._dp

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION dxparg

END MODULE constants_NSWC



SUBROUTINE aord (a, n)
!-----------------------------------------------------------------------
! THE AORD SORTING PROCEDURE IS USED TO REORDER THE ELEMENTS OF A SO THAT
! ABS(A(I)) .LE. ABS(A(I+1)) FOR I = 1,...,N-1.  IT IS ASSUMED THAT N >= 1.
!-----------------------------------------------------------------------
USE constants_NSWC
IMPLICIT NONE

REAL (dp), INTENT(IN OUT) :: a(:)
INTEGER, INTENT(IN)       :: n

! Local variables
INTEGER   :: k(10) = (/ 1, 4, 13, 40, 121, 364, 1093, 3280, 9841, 28524 /), &
             imax, i, ii, ki, jmax, j, l, ll
REAL (dp) :: s
!------------------------

!             SELECTION OF THE INCREMENTS K(I) = (3**I-1)/2

IF (n < 2) RETURN
imax = 1
DO i = 3,10
  IF (n <= k(i)) EXIT
  imax = imax + 1
END DO

!            STEPPING THROUGH THE INCREMENTS K(IMAX),...,K(1)

i = imax
DO ii = 1,imax
  ki = k(i)
  
!             SORTING ELEMENTS THAT ARE KI POSITIONS APART
!                 SO THAT ABS(A(J)) .LE. ABS(A(J+KI))
  
  jmax = n - ki
  DO j = 1,jmax
    l = j
    ll = j + ki
    s = a(ll)
    DO
      IF (ABS(s) >= ABS(a(l))) EXIT
      a(ll) = a(l)
      ll = l
      l = l - ki
      IF (l <= 0) EXIT
    END DO
    a(ll) = s
  END DO
  
  i = i - 1
END DO
RETURN
END SUBROUTINE aord


FUNCTION cbrt (x) RESULT(fn_val)
!-----------------------------------------------------------------------
!                   CUBE ROOT OF A REAL NUMBER
!-----------------------------------------------------------------------
USE constants_NSWC
IMPLICIT NONE

REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

!     Local variable
REAL (dp) :: r, zero = 0.0_dp, three = 3.0_dp

IF (x < zero) THEN
  r = LOG(-x) / three
  fn_val = -EXP(r)
ELSE IF (x > zero) THEN
  r = LOG(x) / three
  fn_val = EXP(r)
ELSE
  fn_val = zero
END IF

RETURN
END FUNCTION cbrt


SUBROUTINE qdcrt (a, z)
!-----------------------------------------------------------------------

!        QDCRT COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z + A(3)*Z**2
!     AND STORES THE RESULTS IN Z.  IT IS ASSUMED THAT A(3) IS NONZERO.
!-----------------------------------------------------------------------
!     Converted to be compatible with ELF90 by
!        Alan Miller
!        amiller @ bigpond.net.au
!     WWW-page: http://users.bigpond.net.au/amiller
!     Latest revision - 27 February 1997
!-----------------------------------------------------------------------

USE constants_NSWC
IMPLICIT NONE

REAL (dp), INTENT(IN)     :: a(:)
COMPLEX (dp), INTENT(OUT) :: z(:)
!-----------------------------------------------------------------------

! Local variables
REAL (dp) :: d, eps, r, w, x, y, zero = 0.0_dp

!     ***** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!           SMALLEST NUMBER SUCH THAT 1.0 + EPS .GT. 1.0.

eps = dpmpar(1)

!-------------------
IF (a(1) == zero) GO TO 40
d = a(2)*a(2) - 4.0D0*a(1)*a(3)
IF (ABS(d) <= 2.0D0*eps*a(2)*a(2)) GO TO 20
r = SQRT(ABS(d))
IF (d < zero) GO TO 30

!                 DISTINCT REAL ROOTS

IF (a(2) /= zero) GO TO 10
x = ABS(0.5D0*r/a(3))
z(1) = CMPLX(x, zero, dp)
z(2) = CMPLX(-x, zero, dp)
RETURN

10 w = -(a(2) + SIGN(r,a(2)))
z(1) = CMPLX(2.0D0*a(1)/w, zero, dp)
z(2) = CMPLX(0.5D0*w/a(3), zero, dp)
RETURN

!                  EQUAL REAL ROOTS

20 z(1) = CMPLX(-0.5D0*a(2)/a(3), zero, dp)
z(2) = z(1)
RETURN

!                   COMPLEX ROOTS

30 x = -0.5D0*a(2)/a(3)
y = ABS(0.5D0*r/a(3))
z(1) = CMPLX(x, y, dp)
z(2) = CMPLX(x,-y, dp)
RETURN

!                 CASE WHEN A(1) = 0

40 z(1) = CMPLX(zero, zero, dp)
z(2) = CMPLX(-a(2)/a(3), zero, dp)
RETURN
END SUBROUTINE qdcrt


SUBROUTINE cbcrt (a, z)
!-----------------------------------------------------------------------

!        CBCRT COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z + A(3)*Z**2 + A(4)*Z**3
!     AND STORES THE RESULTS IN Z. IT IS ASSUMED THAT A(4) IS NONZERO.

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS
!        NAVAL SURFACE WEAPONS CENTER
!        DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
!     Converted to be compatible with ELF90 by
!        Alan Miller
!        amiller @ bigpond.net.au
!     WWW-page: http://users.bigpond.net.au/amiller
!     Latest revision - 27 February 1997
!-----------------------------------------------------------------------
USE constants_NSWC
IMPLICIT NONE

REAL (dp), INTENT(IN)     :: a(:)
COMPLEX (dp), INTENT(OUT) :: z(:)

INTERFACE
  SUBROUTINE qdcrt (a, z)
    USE constants_NSWC
    IMPLICIT NONE
    REAL (dp), INTENT(IN)     :: a(:)
    COMPLEX (dp), INTENT(OUT) :: z(:)
  END SUBROUTINE qdcrt
  FUNCTION cbrt (x) RESULT(fn_val)
    USE constants_NSWC
    IMPLICIT NONE
    REAL (dp), INTENT(IN) :: x
    REAL (dp)             :: fn_val
  END FUNCTION cbrt
END INTERFACE

!-------------------
REAL (dp) :: aq(3), arg, c, cf, d, eps, p, p1, q, q1, r, ra, rb, rq, rt, &
             rt3 = 1.7320508075689D0, r1, s, sf, sq, sum, t, tol, t1, w, &
             w1, w2, x, x1, x2, x3, y, y1, y2, y3, zero = 0.0_dp
!-----------------------------------------------------------------------

!     ***** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!           SMALLEST NUMBER SUCH THAT 1.0 + EPS .GT. 1.0.

eps = dpmpar(1)

!-------------------
IF (a(1) == zero) GO TO 100
p = a(3)/(3.0D0*a(4))
q = a(2)/a(4)
r = a(1)/a(4)
tol = 4.0D0*eps

c = zero
t = a(2) - p*a(3)
IF (ABS(t) > tol*ABS(a(2))) c = t/a(4)

t = 2.0D0*p*p - q
IF (ABS(t) <= tol*ABS(q)) t = zero
d = r + p*t
IF (ABS(d) <= tol*ABS(r)) GO TO 110

!           SET  SQ = (A(4)/S)**2 * (C**3/27 + D**2/4)

s = MAX(ABS(a(1)), ABS(a(2)), ABS(a(3)))
p1 = a(3)/(3.0D0*s)
q1 = a(2)/s
r1 = a(1)/s

t1 = q - 2.25D0*p*p
IF (ABS(t1) <= tol*ABS(q)) t1 = zero
w = 0.25D0*r1*r1
w1 = 0.5D0*p1*r1*t
w2 = q1*q1*t1/27.0D0
IF (w1 < zero) GO TO 10
w = w + w1
sq = w + w2
GO TO 12
10 IF (w2 < zero) GO TO 11
w = w + w2
sq = w + w1
GO TO 12
11 sq = w + (w1 + w2)
12 IF (ABS(sq) <= tol*w) sq = zero
rq = ABS(s/a(4))*SQRT(ABS(sq))
IF (sq >= zero) GO TO 40

!                   ALL ROOTS ARE REAL

arg = ATAN2(rq, -0.5D0*d)
cf = COS(arg/3.0D0)
sf = SIN(arg/3.0D0)
rt = SQRT(-c/3.0D0)
y1 = 2.0D0*rt*cf
y2 = -rt*(cf + rt3*sf)
y3 = -(d/y1)/y2

x1 = y1 - p
x2 = y2 - p
x3 = y3 - p
IF (ABS(x1) <= ABS(x2)) GO TO 20
t = x1
x1 = x2
x2 = t
20 IF (ABS(x2) <= ABS(x3)) GO TO 30
t = x2
x2 = x3
x3 = t
IF (ABS(x1) <= ABS(x2)) GO TO 30
t = x1
x1 = x2
x2 = t

30 w = x3
IF (ABS(x2) < 0.1D0*ABS(x3)) GO TO 70
IF (ABS(x1) < 0.1D0*ABS(x2)) x1 = - (r/x3)/x2
z(1) = CMPLX(x1, zero, dp)
z(2) = CMPLX(x2, zero, dp)
z(3) = CMPLX(x3, zero, dp)
RETURN

!                  REAL AND COMPLEX ROOTS

40 ra = cbrt(-0.5D0*d - SIGN(rq,d))
rb = -c/(3.0D0*ra)
t = ra + rb
w = -p
x = -p
IF (ABS(t) <= tol*ABS(ra)) GO TO 41
w = t - p
x = -0.5D0*t - p
IF (ABS(x) <= tol*ABS(p)) x = zero
41 t = ABS(ra - rb)
y = 0.5D0*rt3*t

IF (t <= tol*ABS(ra)) GO TO 60
IF (ABS(x) < ABS(y)) GO TO 50
s = ABS(x)
t = y/x
GO TO 51
50 s = ABS(y)
t = x/y
51 IF (s < 0.1D0*ABS(w)) GO TO 70
w1 = w/s
sum = 1.0D0 + t*t
IF (w1*w1 < 0.01D0*sum) w = - ((r/sum)/s)/s
z(1) = CMPLX(w,zero, dp)
z(2) = CMPLX(x, y, dp)
z(3) = CMPLX(x,-y, dp)
RETURN

!               AT LEAST TWO ROOTS ARE EQUAL

60 IF (ABS(x) < ABS(w)) GO TO 61
IF (ABS(w) < 0.1D0*ABS(x)) w = - (r/x)/x
z(1) = CMPLX(w, zero, dp)
z(2) = CMPLX(x, zero, dp)
z(3) = z(2)
RETURN
61 IF (ABS(x) < 0.1D0*ABS(w)) GO TO 70
z(1) = CMPLX(x, zero, dp)
z(2) = z(1)
z(3) = CMPLX(w, zero, dp)
RETURN

!     HERE W IS MUCH LARGER IN MAGNITUDE THAN THE OTHER ROOTS.
!     AS A RESULT, THE OTHER ROOTS MAY BE EXCEEDINGLY INACCURATE
!     BECAUSE OF ROUNDOFF ERROR.  TO DEAL WITH THIS, A QUADRATIC
!     IS FORMED WHOSE ROOTS ARE THE SAME AS THE SMALLER ROOTS OF
!     THE CUBIC.  THIS QUADRATIC IS THEN SOLVED.

!     THIS CODE WAS WRITTEN BY WILLIAM L. DAVIS (NSWC).

70 aq(1) = a(1)
aq(2) = a(2) + a(1)/w
aq(3) = -a(4)*w
CALL qdcrt(aq, z)
z(3) = CMPLX(w, zero, dp)

IF (AIMAG(z(1)) == zero) RETURN
z(3) = z(2)
z(2) = z(1)
z(1) = CMPLX(w, zero, dp)
RETURN
!-----------------------------------------------------------------------

!                  CASE WHEN A(1) = 0

100 z(1) = CMPLX(zero, zero, dp)
CALL qdcrt(a(2:), z(2:))
RETURN

!                   CASE WHEN D = 0

110 z(1) = CMPLX(-p, zero, dp)
w = SQRT(ABS(c))
IF (c < zero) GO TO 120
z(2) = CMPLX(-p, w, dp)
z(3) = CMPLX(-p,-w, dp)
RETURN

120 IF (p /= zero) GO TO 130
z(2) = CMPLX(w, zero, dp)
z(3) = CMPLX(-w, zero, dp)
RETURN

130 x = -(p + SIGN(w,p))
z(3) = CMPLX(x, zero, dp)
t = 3.0D0*a(1)/(a(3)*x)
IF (ABS(p) > ABS(t)) GO TO 131
z(2) = CMPLX(t, zero, dp)
RETURN
131 z(2) = z(1)
z(1) = CMPLX(t, zero, dp)
RETURN
END SUBROUTINE cbcrt


SUBROUTINE qtcrt (a, z)
!-----------------------------------------------------------------------

!         QTCRT COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!               A(1) + A(2)*Z + ... + A(5)*Z**4
!         AND STORES THE RESULTS IN Z.  IT IS ASSUMED THAT A(5)
!         IS NONZERO.

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS
!        NAVAL SURFACE WEAPONS CENTER
!        DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
!     Converted to be compatible with ELF90 by
!        Alan Miller
!        amiller @ bigpond.net.au
!     WWW-page: http://users.bigpond.net.au/amiller
!     Latest revision - 27 February 1997
!-----------------------------------------------------------------------
USE constants_NSWC
IMPLICIT NONE

REAL (dp), INTENT(IN)     :: a(:)
COMPLEX (dp), INTENT(OUT) :: z(:)

INTERFACE
  SUBROUTINE cbcrt (a, z)
    USE constants_NSWC
    IMPLICIT NONE
    REAL (dp), INTENT(IN)     :: a(:)
    COMPLEX (dp), INTENT(OUT) :: z(:)
  END SUBROUTINE cbcrt
  SUBROUTINE aord (a, n)
    USE constants_NSWC
    IMPLICIT NONE
    REAL (dp), INTENT(IN OUT) :: a(:)
    INTEGER, INTENT(IN)       :: n
  END SUBROUTINE aord
END INTERFACE

!     Local variables
COMPLEX (dp) :: w
REAL (dp)    :: b, b2, c, d, e, h, p, q, r, t, temp(4), u, v, v1, v2,  &
                x, x1, x2, x3, y, zero = 0.0_dp

IF (a(1) == zero) GO TO 100
b = a(4)/(4.0D0*a(5))
c = a(3)/a(5)
d = a(2)/a(5)
e = a(1)/a(5)
b2 = b*b

p = 0.5D0*(c - 6.0D0*b2)
q = d - 2.0D0*b*(c - 4.0D0*b2)
r = b2*(c - 3.0D0*b2) - b*d + e

!          SOLVE THE RESOLVENT CUBIC EQUATION. THE CUBIC HAS AT LEAST ONE
!          NONNEGATIVE REAL ROOT.  IF W1, W2, W3 ARE THE ROOTS OF THE CUBIC
!          THEN THE ROOTS OF THE ORIGINIAL EQUATION ARE

!             Z = -B + CSQRT(W1) + CSQRT(W2) + CSQRT(W3)

!          WHERE THE SIGNS OF THE SQUARE ROOTS ARE CHOSEN SO
!          THAT CSQRT(W1) * CSQRT(W2) * CSQRT(W3) = -Q/8.

temp(1) = -q*q/64.0D0
temp(2) = 0.25D0*(p*p - r)
temp(3) =  p
temp(4) = 1.0D0
CALL cbcrt(temp,z)
IF (AIMAG(z(2)) /= zero) GO TO 60

!               THE RESOLVENT CUBIC HAS ONLY REAL ROOTS
!                REORDER THE ROOTS IN INCREASING ORDER

x1 = DBLE(z(1))
x2 = DBLE(z(2))
x3 = DBLE(z(3))
IF (x1 <= x2) GO TO 10
t = x1
x1 = x2
x2 = t
10 IF (x2 <= x3) GO TO 20
t = x2
x2 = x3
x3 = t
IF (x1 <= x2) GO TO 20
t = x1
x1 = x2
x2 = t

20 u = zero
IF (x3 > zero) u = SQRT(x3)
IF (x2 <= zero) GO TO 41
IF (x1 >= zero) GO TO 30
IF (ABS(x1) > x2) GO TO 40
x1 = zero

30 x1 = SQRT(x1)
x2 = SQRT(x2)
IF (q > zero) x1 = -x1
temp(1) = (( x1 + x2) + u) - b
temp(2) = ((-x1 - x2) + u) - b
temp(3) = (( x1 - x2) - u) - b
temp(4) = ((-x1 + x2) - u) - b
CALL aord (temp,4)
IF (ABS(temp(1)) >= 0.1D0*ABS(temp(4))) GO TO 31
t = temp(2)*temp(3)*temp(4)
IF (t /= zero) temp(1) = e/t
31 z(1) = CMPLX(temp(1), zero, dp)
z(2) = CMPLX(temp(2), zero, dp)
z(3) = CMPLX(temp(3), zero, dp)
z(4) = CMPLX(temp(4), zero, dp)
RETURN

40 v1 = SQRT(ABS(x1))
v2 = zero
GO TO 50
41 v1 = SQRT(ABS(x1))
v2 = SQRT(ABS(x2))
IF (q < zero) u = -u

50 x = -u - b
y = v1 - v2
z(1) = CMPLX(x, y, dp)
z(2) = CMPLX(x,-y, dp)
x =  u - b
y = v1 + v2
z(3) = CMPLX(x, y, dp)
z(4) = CMPLX(x,-y, dp)
RETURN

!                THE RESOLVENT CUBIC HAS COMPLEX ROOTS

60 t = DBLE(z(1))
x = zero
IF (t < zero) THEN
  GO TO 61
ELSE IF (t == zero) THEN
  GO TO 70
ELSE
  GO TO 62
END IF
61 h = ABS(DBLE(z(2))) + ABS(AIMAG(z(2)))
IF (ABS(t) <= h) GO TO 70
GO TO 80
62 x = SQRT(t)
IF (q > zero) x = -x

70 w = SQRT(z(2))
u = 2.0D0*DBLE(w)
v = 2.0D0*ABS(AIMAG(w))
t =  x - b
x1 = t + u
x2 = t - u
IF (ABS(x1) <= ABS(x2)) GO TO 71
t = x1
x1 = x2
x2 = t
71 u = -x - b
h = u*u + v*v
IF (x1*x1 < 0.01D0*MIN(x2*x2,h)) x1 = e/(x2*h)
z(1) = CMPLX(x1, zero, dp)
z(2) = CMPLX(x2, zero, dp)
z(3) = CMPLX(u, v, dp)
z(4) = CMPLX(u,-v, dp)
RETURN

80 v = SQRT(ABS(t))
z(1) = CMPLX(-b, v, dp)
z(2) = CMPLX(-b,-v, dp)
z(3) = z(1)
z(4) = z(2)
RETURN

!                         CASE WHEN A(1) = 0

100 z(1) = CMPLX(zero, zero, dp)
CALL cbcrt(a(2:), z(2:))
RETURN
END SUBROUTINE qtcrt
