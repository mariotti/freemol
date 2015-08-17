!
!H
!H-----------------------------------------------------------------------------
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H MODULE dummy_minuit_function by F.Mariotti
!H-----------------------------------------------------------------------------
!H $Id: minuit_function.F90,v 1.1.1.1 2009/01/12 16:56:17 mariotti Exp $
!H-----------------------------------------------------------------------------
!H 
!H DESCRIPTION:
!H This is a simple example for MINUIT minimization program
!H and contains a simple example of a function to be minimised.
!H It contains a single subroutine which is the actual user function
!H to be minimised. I DO STRONGLY suggest to write in here only an
!H interface to the real function and not the real function itself.
!H This help in maintaining... hopefully..
!H By the way .. the example is the real function!! ;-)
!H
!H It does not contain any of the standard functions of Freemol
!H programming environment. So no "minuit_function_init"!!!
!H
!H In order to use it remove the dummy_ string from the name of the module.
!H It is here only as template.
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
module dummy_minuit_function
  !
  ! We use default varible types from Freemol Package.
  !---------------------------------------------------
  use vartypes
  use messages
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
    ! This example uses the z = (x^2+y^2) function
    ! or the z = ((x+4)^2+y^2)
    !---------------------------------------------
    !
    ! We check how we are called: IFLAG
    !-------------------------------------------
    if (iflag.eq.1) then
       ! Read input data: Here we do nothing
       !------------------------------------
       call message(MESOUT,'Dummy module: Request: Read Input')
    else if (iflag.eq.2) then
       ! Calculate GRAD, the first derivatives of FVAL.
       !-----------------------------------------------
       !grad(1) = 2.0_FREAL * xval(1)
       !grad(2) = 2.0_FREAL * xval(2)
       grad(1) = 2.0_FREAL * (xval(1) + 4.0_FREAL)
       grad(2) = 2.0_FREAL * xval(2)
    end if
    !
    ! Here we evaluate the function value
    !------------------------------------
    fval = ( (xval(1) + 4.0_FREAL)**2 + xval(2)**2)
    !
    ! The last iflag is evaluated here
    !---------------------------------
    if (iflag.eq.3) then
       ! Final calculations end print out
       !---------------------------------
       write(*,*) 'Calculation performed'
       write(*,*) 'Parameter 1:',xval(1)
       write(*,*) 'Parameter 2:',xval(2)
       write(*,*) 'Final Value:',fval
    end if
    !
    return
    !
  end subroutine minuit_fcn
end module dummy_minuit_function
