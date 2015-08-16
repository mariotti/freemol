!
!H
!H-----------------------------------------------------------------------------
!H This File is part of the Frimol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H MODULE csm_minuit_function by F.Mariotti: (c) F.Mariotti
!H-----------------------------------------------------------------------------
!H $Id: csm_minuit_function.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
!H-----------------------------------------------------------------------------
!H 
!H DESCRIPTION:
!H This is a the function called by minuit in order to minimise the CSM.
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
module csm_minuit_function
  !
  ! We use default varible types from Frimol Package.
  !---------------------------------------------------
  use vartypes
  use messages
  use modcsmmin
  use molecule
  implicit none
  private
  !
  public minuit_fcn
  !
contains
  !
!
!H
!H-----------------------------------------------------------------------------
!H subroutine minuit_fcn(npar,grad,fval,xval,iflag) 
!H-----------------------------------------------------------------------------
!H 
!  subroutine minuit_fcn(npar,gin,f,x,iflag) 
  subroutine minuit_fcn(npar,grad,fval,xval,iflag)
!
!H
!H-----------------------------------------------------------------------------
!H  
!H  From the manual of Minuit we read:
!H  INPUT parameters:
!H  npar	number of currently variable parameters.
!H  xval	vector of (constant orvariable) parameters.
!H  iflag	indicate what is to be calculated.
!H  (futil)	removed in this version.
!H  OUTPUT parameters:
!H  fval	The calculated function value.
!H  grad	The (optional) vector of first derivatives.
!H  
!H  You should declare the parameters as follow:
!H  npar	integer(FINT), intent(in) :: npar
!H  xval	real(FREAL), dimension(:), intent(in) ::  xval
!H  iflag	integer(FINT), intent(in) :: iflag
!H  fval	real(FREAL), intent(out) :: fval
!H  grad	real(FREAL), intent(out), optional :: grad
!H  
!H  
!H  
!H  
!H  
!H  
!H-----------------------------------------------------------------------------
!
    !
    ! FCN input and output parameters
    !--------------------------------
    integer(FINT), intent(in) :: npar
    real(FREAL), dimension(:), intent(in) ::  xval
    integer(FINT), intent(in) :: iflag
    real(FREAL), intent(out) :: fval
    real(FREAL), intent(out), dimension(:), optional :: grad
    !
    integer(FINT) :: irc
    !integer(FINT) :: nat
    !
    !
    ! We check how we are called: IFLAG
    !-------------------------------------------
    !write(*,'(1X,"#[csm_minuit_function]: Called FCN with IFLAG:",I2 )') iflag
    if (iflag.eq.1) then
       ! Init the function check consistency
       !------------------------------------
       if(npar.ne.6) then
          call message(MESERRO,"[csm_minuit_function] inconsistency.")
          call message(MESKILL,"[csm_minuit_function] require stop.")
       end if
       !
       ! Reset Coordinates
       !------------------
       irc = modcsmmin_restore()
       !
       call message(MESOUT,'[csm_minuit_function] Initialised.')
    else if (iflag.eq.2) then
       ! We cannot calculate GRAD now
       !-----------------------------------------------
       if(present(grad)) then
          grad(1)=0.0_FREAL
       endif
       call message(MESERRO,"[csm_minuit_function] Gradients not implemented.")
       call message(MESKILL,"[csm_minuit_function] require stop.")
    end if
    !
    ! Here we evaluate the function value
    !------------------------------------
    fval = modcsmmin_minuit_value(xval(1:6))
    if(iflag.ne.4) then
       ! print coords
       call molecule_print(6,'mldfromau')
    end if
    !
    ! The last iflag is evaluated here
    !---------------------------------
    if (iflag.eq.3) then
       ! Final calculations end/or print out
       !------------------------------------
       call message_value(MESOUT,'[csm_minuit_function]: Current CSM: ',fval)
       call message(MESOUT,'[csm_minuit_function]: Parameters:.')
       call message_value(MESOUT,'[csm_minuit_function]: x trans ',xval(1))
       call message_value(MESOUT,'[csm_minuit_function]: y trans ',xval(2))
       call message_value(MESOUT,'[csm_minuit_function]: z trans ',xval(3))
       call message_value(MESOUT,'[csm_minuit_function]: x   rot ',xval(4))
       call message_value(MESOUT,'[csm_minuit_function]: y   rot ',xval(5))
       call message_value(MESOUT,'[csm_minuit_function]: z   rot ',xval(6))
       call message(MESOUT,'[csm_minuit_function]: Coordinates:.')
       call molecule_print(6,'mldfromau')
    end if
    !
    return
    !
  end subroutine minuit_fcn
end module csm_minuit_function
