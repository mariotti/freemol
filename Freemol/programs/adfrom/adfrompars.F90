!
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!
!
module adfrompars
  !
  use vartypes
  !
  ! Define sto versions
  !
  integer(SINT), parameter :: STO_FMT_CART_SHORT = 1 !cart sto per atom type
  integer(SINT), parameter :: STO_FMT_CART_LONG = 2 !cart sto per atom
  integer(SINT), parameter :: STO_FMT_MOLDEN = STO_FMT_CART_LONG
  integer(SINT), parameter :: STO_FMT_SPHER_SHORT = 3 !spherical sto short
  integer(SINT), parameter :: STO_FMT_SPHER_LONG = 4 !spherical sto long
  integer(SINT), parameter :: STO_FMT_CART_ZCONTR = 5 !with zeta contraction
  integer(SINT), parameter :: STO_FMT_SPHER_ZCONTR = 6 !with zeta contraction
  integer(SINT), parameter :: STO_FMT_CART_SZCONTR = 7 !short zeta contraction
  integer(SINT), parameter :: STO_FMT_SPHER_SZCONTR = 8 !short zeta contraction
  !
  integer(SINT), parameter :: MO_FMT_LINES_CONTR = 1 !in lines exclude unused
  integer(SINT), parameter :: MO_FMT_MOLDEN = MO_FMT_LINES_CONTR
  integer(SINT), parameter :: MO_FMT_LIST_CONTR = 2 !multi lines exclude unused
  !
end module adfrompars
