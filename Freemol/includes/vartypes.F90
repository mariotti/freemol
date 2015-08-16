!
!H
!H  ---------------------------------------------------------------------------
!H  ---------------------------------------------------------------------------
!H  MODULE vartypes
!H  ---------------------------------------------------------------------------
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H  ---------------------------------------------------------------------------
!H  $Id: vartypes.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
!H  ---------------------------------------------------------------------------
!H  
!H  KeyWords: variables kind type character lenght real integer module
!H  
!H  These Modules define the variable types used in FREEMOL.
!H  Some of the definitions are equivalent To ADF1999 definitions in order to
!H  Interface some tools with ADF1999. It is supposed to work as well with
!H  ADF 2000.XX versions.
!H  
!H  Usage:
!H  use vartypes
!H
!H  Notes: define integer and real kind and default string lenght
!H  
!H  
!H  
!H  
!H  ---------------------------------------------------------------------------
!H  
!

module vartypes
  !
  implicit none
  public
  !
  ! The Biggest Integer
  ! TODO!!!
  integer, parameter :: FBIGGEST_INT  = kind(65560)
  !
  ! ADF compatibility
  !
  integer, parameter :: KINT  = kind(1)
  integer, parameter :: KREAL = kind(1.0d0)
  integer, parameter :: LCHARS = 160
  !
  ! Freemol Default
  !
  integer, parameter :: FINT  = kind(1)
  integer, parameter :: FREAL = kind(1.0d0)
  integer, parameter :: FLCHARS = 860
  !
  ! Freemol Small: for special use
  !
  integer, parameter :: SINT  = kind(1)
  integer, parameter :: SREAL = kind(1.0d0)
  integer, parameter :: SLCHARS = 10
  !
  ! Freemol Big: for special use
  !
  integer, parameter :: BINT  = kind(1)
  integer, parameter :: BREAL = kind(1.0d0)
  integer, parameter :: BLCHARS = 320
  !
  ! Some Machine parameters
  ! We state that machine precision is a limit!
  ! That's why we get one digit above: I won't play with machine precision!!
  !-------------------------------------------------------------------------
  !
  real(FREAL), parameter :: F_EPS = EPSILON(1.0_FREAL)*10.0_FREAL
  real(SREAL), parameter :: S_EPS = EPSILON(1.0_SREAL)*10.0_SREAL
  real(BREAL), parameter :: B_EPS = EPSILON(1.0_BREAL)*10.0_BREAL
  !
contains
!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H subroutine vartypes_prt_eps()
!H-----------------------------------------------------------------------------
!H
!
  subroutine vartypes_prt_eps()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine Prints internal EPS as it is obtained from F90 EPSILON value
!H times*10.0_FREAL with an  hash on top!
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    write(*,'("# FREAL EPS",1X,G24.14)') F_EPS
    write(*,'("# SREAL EPS",1X,G24.14)') S_EPS
    write(*,'("# BREAL EPS",1X,G24.14)') B_EPS
    !
  end subroutine vartypes_prt_eps
!
end module vartypes
