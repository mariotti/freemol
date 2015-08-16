!
!H
!H-----------------------------------------------------------------------------
!H This File is part of the Frimol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H MODULE fit1Dpol_commons by F.Mariotti
!H-----------------------------------------------------------------------------
!H $Id: fit1Dpol_commons.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
!H-----------------------------------------------------------------------------
!H 
!H DESCRIPTION:
!H
!H A "modern" F90 common.
!H Contains all the data that should flows in all the modules!
!H 
!H TODO: 
!H 
!H 
!H WARN: 
!H 
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H 
!
module fit1Dpol_commons
!
  use vartypes
  !use messages
  implicit none
  !
  ! Just in case
  !-------------
  private
  !
  ! And now the public staff
  !-------------------------
  !
  ! fortran parameters
  !-------------------
  integer(FINT), public, parameter :: FLG_SIZE=7 ! f1dp_flags=Current Size
  !
  integer(FINT), public, parameter :: FLG_POLO=1 ! f1dp_flags=1 polyd dimension (for multidim)
  integer(FINT), public, parameter :: FLG_FUNC=2 ! f1dp_flags=2 type of fit function
  integer(FINT), public, parameter :: FLG_GPLT=3 ! f1dp_flags=3 produce gnuplot output
  integer(FINT), public, parameter :: FLG_IEPS=4 ! EPS is given in the input
  integer(FINT), public, parameter :: FLG_MINT=5 ! Minuit Type of run
  integer(FINT), public, parameter :: FLG_EINP=6 ! Echo of input: 1,in std out 2,infile(not yet)
  integer(FINT), public, parameter :: FLG_GNRT=7 ! GeNeRaTe: print the fitted function values
  !
  ! optional data for parameter: currently 2: 1,the initial value; 2,the step size(optional in input)
  integer(FINT), public, parameter :: PAR_SIZE=2 ! This is really code dependent check if you change it!
  !
  ! parameters for functions
  !-------------------------
  integer(FINT), public, parameter :: FUN_DUMMY=0 ! The test dummy function
  integer(FINT), public, parameter :: FUN_POLYN=1 ! A polynomial expansion
  integer(FINT), public, parameter :: FUN_MORSE=2 ! A simple morse function
  integer(FINT), public, parameter :: FUN_PMORS=3 ! Not yet implemented
  integer(FINT), public, parameter :: FUN_DUNH1=4 ! Not yet implemented
  integer(FINT), public, parameter :: FUN_DUNH2=5 ! Not yet implemented
  integer(FINT), public, parameter :: FUN_TEMPO=6 ! Not yet implemented
  integer(FINT), public, parameter :: FUN_POLYS=7 ! A polynomial expansion with shift of the X coordinates
  integer(FINT), public, parameter :: FUN_POLYD=8 ! A polynomial expansion of a given dimension
  !
  ! File units
  !-----------
  integer(FINT), public, save :: f1dp_iunin,f1dp_iunout
  integer(FINT), public, save :: f1dp_iungpltdata,f1dp_iungpltgplt
  integer(FINT), public, save :: f1dp_iungnrtdata,f1dp_iungnrtgnrt
  !
  ! and names
  !---------
  character(FLCHARS), public, save :: f1dp_gpltprefix
  character(FLCHARS), public, save :: f1dp_gnrtprefix
  !
  ! Local Data
  !-----------
  integer(FINT), public, save, dimension(FLG_SIZE) :: f1dp_flags
  integer(FINT), public, save :: f1dp_debug ! It is the verbosity level if 1 or 2 (debug if >2)
  character(FLCHARS), public, save :: f1dp_sessionname
  !
  ! Fit Data
  !---------
  integer(FINT), public, save :: f1dp_numpar
  real(FREAL),   public, save, dimension(:,:), allocatable :: f1dp_pars !note: currently the second dimension should be 2
  character(10), public, save, dimension(:), allocatable :: f1dp_parnames
  character(1),  public, save, dimension(:), allocatable :: f1dp_strfixed
  integer(FINT), public, save, dimension(:), allocatable :: f1dp_idxpars
  integer(FINT), public, save, dimension(:), allocatable :: f1dp_fixpars
  integer(FINT), public, save :: f1dp_numrec
  integer(FINT), public, save :: f1dp_numcol
  integer(FINT), public, save, dimension(:), allocatable :: f1dp_idxdata
  real(FREAL),   public, save, dimension(:,:), allocatable :: f1dp_data
  real(FREAL),   public, save, dimension(:,:), allocatable :: f1dp_cpydata !This is a copy of the original data
                                                                           ! if partial fit or combinations are set
  real(FREAL),   public, save, dimension(:), allocatable :: f1dp_ypol
  character(10), public, save :: f1dp_funame
  !
  ! Minuit Data
  !------------
  real(FREAL),   public, save :: f1dp_eps ! Input eps (accuracy)
  real(FREAL),   public, save :: f1dp_save_minuit, f1dp_error_minuit, f1dp_dummy_minuit, f1dp_svpar_minuit
  real(FREAL),   public, save, dimension(10) :: f1dp_temp_minuit
  real(FREAL),   public, save, dimension(6) :: f1dp_pars_minuit,f1dp_svbestpars
  character(FLCHARS), public, save :: f1dp_minuit_savename
  !
  !Make functions public (the ones that should be public!)
  !-------------------------------------------------------
  public fit1Dpol_commons_reset
  !
  contains
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function fit1Dpol_commons_reset()
!H-----------------------------------------------------------------------------
!H
!
    integer(FINT) function fit1Dpol_commons_reset()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine reset (put to zero or default) the common data
!H
!H-----------------------------------------------------------------------------
!H
!
      !
      fit1Dpol_commons_reset = 0
      !
      ! Set default debug/verbosity level now is 1 in production should be 0
      f1dp_debug=0
      !
      f1dp_flags(:)=0
      f1dp_numpar=0
      f1dp_numrec=0
      f1dp_numcol=0
      f1dp_funame="polynom"
      f1dp_sessionname=''
      !
      f1dp_iungpltdata=6
      f1dp_iungpltgplt=6
      !
    end function fit1Dpol_commons_reset

end module fit1Dpol_commons
