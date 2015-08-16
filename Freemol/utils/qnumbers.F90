!
!H 
!H ---------------------------------------------------------------------
!H This File is part of the Frimol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H ---------------------------------------------------------------------
!H  MODULE qnumbers
!H ---------------------------------------------------------------------
!H  Routines to manage quantum numbers
!H ---------------------------------------------------------------------
!H ---------------------------------------------------------------------
!H
!
module qnumbers
  !
  use vartypes
  !
  implicit none
  private
  !
  !parameters
  real(FREAL), private, parameter :: zero = 0.D0, one = 1.D0
  !
  ! Public Procedures
  !
  public qn_ms
  !
contains

  subroutine qn_ms(m,msign,idmn)
    !
    integer(FINT), intent(in) :: idmn
    integer(FINT), dimension(:), intent(out) :: m
    integer(FINT), dimension(:), intent(out) :: msign
    !
    integer(FINT) :: k
    !
    k = (idmn/2)+1
#ifdef DEBUG
    if(ceiling(real(idmn,KIND=FREAL)/2.0_FREAL).ne.(int(idmn/2)+1)) then
       stop ' Internal error qn_ms k even '
    endif
    if((2*k-1).gt.size(m).or.(2*k-1).gt.size(msign)) then
       stop ' Internal error george '
    endif
#endif
    !
    m(1) = 0_FINT
    msign(1) = 0_FINT
#ifdef IFC
    print *,"NO IFC COMPATIBLE CODE: STOP"
    stop 3
! This code??
!    m(2:k) = (/1:k-1:1/)
#endif
!
#ifdef GFORTRAN
    print *,"NO IFC COMPATIBLE CODE: STOP"
    stop 3
! This code??
!    m(2:k) = (/1:k-1:1/)
#endif
    print *,"NO IFC COMPATIBLE CODE: STOP"
    stop 3
! This code??
!    m(2:k) = (/1:k-1:1/)
!
!    m(2:k) = (/1:k-1:1/)
!
    msign(2:k) = 0_FINT
    m(k+1:2*k-1) = m(2:k)
    msign(k+1:2*k-1) = 1_FINT
    !
  end subroutine qn_ms
  !
end module qnumbers

