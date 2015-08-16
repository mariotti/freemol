!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H MODULE CHEMCONST 
!H-----------------------------------------------------------------------------
!H This File is part of the Frimol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H $Id: chemconst.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
!H-----------------------------------------------------------------------------
!H 
!H Update 20 Feb 2001: FM added comments
!H In this module we store some of foundamental constants
!H which are more related with a chemical environment.
!H TODO: Actually the constants are stored in a let say random precision!
!H       We have to set up a "constants precision" which will have to
!H       translate to the used precision. In order to achive this we
!H       will have to set up the higher experimental precision which will
!H       have to be hard code. Than we will have to find the higher machine
!H       precision and store the numbers rapresentation for them and
!H       translate it later on to the user program precision.
!H TODO: Define a default precision.
!H TODO: Introduce here as well some units conversions.
!H
!H HINTS: Units conversion should be performed at input level.
!H HINTS: Define an init function which select the required precision
!H        and return numbers in the requested precision.
!H 
!H-----------------------------------------------------------------------------
!H 
!H In practice the module should behave in to ways:
!H 1) propose functions to return constants (numbers) in the required precision
!H 2) propose a default for easy use.
!H 
!H-----------------------------------------------------------------------------
!H 
!H The problem is the interaction. How to get the best precision.
!H Or in other words How do I know the precision of the constant?
!H
!H-----------------------------------------------------------------------------
!H
!H WARN: So I warn you!! The numbers in here are just test level numbers!!
!H Check it on your own!!
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H 
!
module chemconst
  !
  use vartypes
  !
  implicit none
  private
  !
  ! meters to angstrom.
  !--------------------
  ! The best precision we set.
  integer, private, parameter :: CC_AU_M_REAL = kind(1.0000000000D-10)
  ! We set up the number.
  real(CC_AU_M_REAL), public, parameter :: cce_au_m = 1.0000000000D-10
  ! We translate the number to the used precision.
  real(FREAL), public, parameter :: cc_au2m = cce_au_m
  !
  ! bhor to angstroms: just up to the -8 digit!!! get care!!
  !---------------------------------------------------------
  integer, private, parameter :: CC_BHOR_ANGS_REAL = kind(0.52917706D0)
  ! We set up the number.
  real(CC_BHOR_ANGS_REAL), public, parameter :: cce_bhor_angs = 0.52917706D0
  ! We translate the number to the used precision.
  real(FREAL), public, parameter :: cc_au2ang = 0.52917706D0
  !
end module chemconst



