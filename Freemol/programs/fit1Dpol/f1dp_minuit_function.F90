!
!H
!H-----------------------------------------------------------------------------
!H This File is part of the Frimol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H MODULE f1dp_minuit_function by F.Mariotti
!H-----------------------------------------------------------------------------
!H $Id: f1dp_minuit_function.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
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
!H It does not contain any of the standard functions of Frimol
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
module f1dp_minuit_function
  !
  ! We use default varible types from Frimol Package.
  !---------------------------------------------------
  use vartypes
  use messages
  use fit1Dpol_commons! all the common data (local module)
  implicit none
  private
  !
  public minuit_fcn
  public f1dp_minuit_function_init
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
    real(FREAL), intent(out), dimension(:) :: grad
    !
    integer(FINT) :: irc
    !
    ! We siply bypass the problem to the real function
    select case (f1dp_flags(FLG_FUNC))
    case(FUN_DUMMY)
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
          write(f1dp_iunout,*) 'Calculation performed'
          write(f1dp_iunout,*) 'Parameter 1:',xval(1)
          write(f1dp_iunout,*) 'Parameter 2:',xval(2)
          write(f1dp_iunout,*) 'Final Value:',fval
       end if
       !
       return
    case(FUN_POLYN)
       irc=funct_polyn(npar,grad,fval,xval,iflag)
       return
    case(FUN_MORSE)
       irc=funct_morse(npar,grad,fval,xval,iflag)
       !call message(MESERRO,"[minuit_fcn] Some internal inconsistency: MORSE is NOT implemented")
       !call message(MESKILL,"[minuit_fcn] We stop now.")
       return
    case(FUN_POLYS)
       irc=funct_polys(npar,grad,fval,xval,iflag)
       return
    case(FUN_POLYD)
       irc=funct_polyd(npar,grad,fval,xval,iflag)
       return
    case default
       call message(MESERRO,"[minuit_fcn] Some internal inconsistency.")
       call message(MESKILL,"[minuit_fcn] We stop now.")
       return
    end select
    !
    return
    !
  end subroutine minuit_fcn
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function f1dp_minuit_function_init()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function f1dp_minuit_function_init(nfpars)
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize (put to zero or default) the data for minuit.
!H
!H-----------------------------------------------------------------------------
!H
!
    integer(FINT), intent(out) :: nfpars
!
! because we know commons we can eventually change the initial parameters
! like testing for sensible default ..
    !
    ! we default to no errors
    f1dp_minuit_function_init = 0
    !
    ! set npars for the dummy function
    nfpars = 2
    ! We should set here who we are
    !------------------------------
    f1dp_flags(FLG_FUNC)=-1
    if(index(f1dp_funame,'dummy').gt.0) then
       f1dp_flags(FLG_FUNC)=FUN_DUMMY
       nfpars = 2
    else if(index(f1dp_funame,'polynom').gt.0) then
       f1dp_flags(FLG_FUNC)=FUN_POLYN
       nfpars=f1dp_numpar
       !HERE WE SHOULD CHECK THE COMMAND LINE PARAMETER -p N
       !AND SET THE CORRECT STAFF. ALSO WE SHOULD CONSIDER REORDERING
       !FROM THE INPUT INDEX
    else if(index(f1dp_funame,'smorse').gt.0) then
       f1dp_flags(FLG_FUNC)=FUN_MORSE
       nfpars=f1dp_numpar
       !HERE WE SHOULD CHECK IF IT IS OK WITH NUM PARS!
       ! Actually this check can be done in the first function call
       ! with iflag=1 ...
    else if(index(f1dp_funame,'polys').gt.0) then
       f1dp_flags(FLG_FUNC)=FUN_POLYS
       nfpars=f1dp_numpar
    else if(index(f1dp_funame,'polyd').gt.0) then
       f1dp_flags(FLG_FUNC)=FUN_POLYD
       nfpars=f1dp_numpar
    else
       call message(MESERRO,"[f1dp_minuit_function] Internal inconsistency: why are you here?")
       call message(MESERRO,"[f1dp_minuit_function] We didn't set up a default function name?")
       call message(MESKILL,"[f1dp_minuit_function] We stop now")
    end if
    if(f1dp_flags(FLG_FUNC).eq.-1) then
       call message(MESERRO,"[f1dp_minuit_function] unkwon function type.")
       call message(MESKILL,"[f1dp_minuit_function] We Stop now.")
    end if
    !
!
!H
!H-----------------------------------------------------------------------------
!H
!
  end function f1dp_minuit_function_init
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function funct_polyn(npar,grad,fval,xval,iflag)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function funct_polyn(npar,grad,fval,xval,iflag)
!
!H
!H-----------------------------------------------------------------------------
!H
!H The acctual polynom function
!H
!H-----------------------------------------------------------------------------
!H
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
    integer(FINT) :: i,irc,iun,j,k
    character(FLCHARS) :: ctmp
    !
    funct_polyn = 0
    !
    if (iflag.eq.1) then
       ! Read input data: Here we do nothing
       !------------------------------------
       call message(MESOUT,'[funct_polyn] at present nothing to do at function call 1.')
    else if (iflag.eq.2) then
       !
       ! Calculate GRAD, the first derivatives of FVAL.
       !-----------------------------------------------
       !-----------------------------------------------
       grad(:)=0.0_FREAL
       !
       irc = f_poly(f1dp_numpar,xval(1:f1dp_numpar),f1dp_numrec,f1dp_data(1:f1dp_numrec,1),f1dp_ypol(1:f1dp_numrec))
       grad(1)=2.0_FREAL*sum(f1dp_ypol(1:f1dp_numrec)-f1dp_data(1:f1dp_numrec,2))
       grad(2)=2.0_FREAL*sum((f1dp_ypol(1:f1dp_numrec)-f1dp_data(1:f1dp_numrec,2))*f1dp_data(1:f1dp_numrec,1))
       do i=3,f1dp_numpar
          grad(i)=2.0_FREAL*&
               &sum((f1dp_ypol(1:f1dp_numrec)-f1dp_data(1:f1dp_numrec,2))*(f1dp_data(1:f1dp_numrec,1)**(i-1)))
       end do
    end if
    !
    ! Here we evaluate the function value
    !------------------------------------
    irc = f_poly(f1dp_numpar,xval(1:f1dp_numpar),f1dp_numrec,f1dp_data(1:f1dp_numrec,1),f1dp_ypol(1:f1dp_numrec))
    fval=0.0_FREAL
    fval=sum((f1dp_ypol(1:f1dp_numrec)-f1dp_data(1:f1dp_numrec,2))**2)
    !do i=1,f1dp_numrec
    !   fval=fval+(f1dp_ypol(i)-f1dp_data(i,2))**2
    !end do
    !
    ! The last iflag is evaluated here
    !---------------------------------
    if (iflag.eq.3) then
       ! Final calculations end print out
       !---------------------------------
       call message(MESOUT,"[funct_polyn] Calculation Done: fcn call 3 requested.")
       write(f1dp_iunout,'(1X,A)') ' Final values from the fcn call'
       do i=1,f1dp_numpar
          write(f1dp_iunout,'(1X,"Par:",1X,A10,2X,F18.10,2X,"fix: ",A1)') f1dp_parnames(i),xval(i),f1dp_strfixed(i)
       end do
       write(f1dp_iunout,'(1X,"Final FunctionValue:",1X,F18.10)') fval
       !
       ! Check if generate was a request
       !--------------------------------
       if (f1dp_flags(FLG_GNRT).gt.0) then
          if(f1dp_flags(FLG_GNRT).eq.2) then
             iun=f1dp_iungnrtdata
          else
             iun=f1dp_iunout
          endif
          irc = f_poly(f1dp_numpar,xval(1:f1dp_numpar),f1dp_numrec,f1dp_data(1:f1dp_numrec,1),f1dp_ypol(1:f1dp_numrec))
          write(iun,'(1X,"Index",1X,1X,"X value",11X,1X,"Y funct",11X,1X,"Delta**2")')
          do i=1,f1dp_numrec
             write(iun,'(1X,I5,1X,3(1X,F18.10))')&
                  &i,f1dp_data(i,1),f1dp_ypol(i),(f1dp_ypol(i)-f1dp_data(i,2))**2
          end do
          if(f1dp_flags(FLG_GPLT).eq.2) then
             iun=f1dp_iungpltgplt
             call message(MESERRO,"[funct_polyn] generate generate gnuplot file not yet implemented.")
          end if
       end if
       !
       ! Check if gnuplot was a request
       !-------------------------------
       if (f1dp_flags(FLG_GPLT).gt.0) then
          ! set the gnuplot dta unit
          if(f1dp_flags(FLG_GPLT).eq.2) then
             iun=f1dp_iungpltdata
          else
             iun=f1dp_iunout
          endif
          !
          write(iun,'(A1,1X,A)') '#','Gnuplot data file starts here'
          write(iun,'(A1,1X,"numrec=",I6,2X,"numcol=",I6)') '#',f1dp_numrec,f1dp_numcol
          do j=1,f1dp_numrec
             k=0
             do k=1,f1dp_numcol-1
                write(iun,'(F18.10,$)') f1dp_data(j,k)
             end do
             write(iun,'(F18.10)') f1dp_data(j,f1dp_numcol)
          end do
          ! set the gnuplot file unit
          if(f1dp_flags(FLG_GPLT).eq.2) then
             iun=f1dp_iungpltgplt
          else
             iun=f1dp_iunout
          endif
          !
          write(iun,'(A1,1X,A)') '#','Gnuplot autogenerated file starts here'
          write(iun,'("set style data linespoints")')
          write(iun,'(A1,1X,A)') '#','Parameters:'
          do i=1,f1dp_numpar
             write(iun,'(A10,"=",F18.10 )') f1dp_parnames(i),xval(i)
          end do
          write(iun,'(A1,1X,A)') '#','Function:'
          write(iun,'(A,$)') 'myfunct(x)='
          write(iun,'(A10,"+",$)') f1dp_parnames(1)
          do i=2,f1dp_numpar-1
             j=i-1
             write(iun,'(A10,"*(x**(",I6,")) + ",$)') f1dp_parnames(i),j
          end do
          j=f1dp_numpar-1
          write(iun,'(A10,"*(x**(",I6,"))")') f1dp_parnames(f1dp_numpar),j
          !
          j=len_trim(f1dp_gpltprefix)
          ctmp=f1dp_gpltprefix(1:j)//'.dat'
          ctmp=trim(ctmp)
          j=len_trim(ctmp)
          if(f1dp_flags(FLG_GPLT).eq.2) then
             write(iun,'("plot [",F12.6,":",F12.6,"] ",$)') f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             write(iun,'(A,$)') ctmp(1:j)
          else
             write(iun,'(A)') '# This would be the plot command'
             write(iun,'("# plot [",F12.6,":",F12.6,"] ",$)') f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             write(iun,'(A,$)') 'yourdatafile'
          end if
          write(iun,'(A1,$)') "'"
          write(iun,'(",myfunct(x)")')
          write(iun,'("pause -1")')
          write(iun,'("fit myfunct(x) ",A1,A,A1," via ",$)') "'",ctmp(1:j),"'"
          do i=1,f1dp_numpar-1
             j=len_trim(f1dp_parnames(i))
             write(iun,'(A,",",$)') f1dp_parnames(i)(1:j)
          end do
          j=len_trim(f1dp_parnames(f1dp_numpar))
          write(iun,'(A)') f1dp_parnames(f1dp_numpar)(1:j)
          !
          if(f1dp_flags(FLG_GPLT).eq.2) then
             write(iun,'("plot [",F12.6,":",F12.6,"] ",$)') f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             j=len_trim(ctmp)
             write(iun,'(A,$)') ctmp(1:j)
          else
             write(iun,'(A)') '# This would be the plot command'
             write(iun,'("# plot [",F12.6,":",F12.6,"] ",$)') f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             write(iun,'(A,$)') 'yourdatafile'
          end if
          write(iun,'(A1,$)') "'"
          write(iun,'(",myfunct(x)")')
          write(iun,'("pause -1")')
          !
       end if
       !
    end if
    !
    return
    !
  end function funct_polyn
!
!H
!H-----------------------------------------------------------------------------
!H
!


!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function f_poly(npc,pc,npt,vpt,ypt)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function f_poly(npc,pc,npt,vpt,ypt)
!
!H
!H-----------------------------------------------------------------------------
!H
!H A tools function to evaluate a polynimial
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    integer(FINT), intent(in) :: npc ! Number of coefficients
    integer(FINT), intent(in) :: npt ! number of (x,y) points
    real(FREAL), dimension(:), intent(in) :: pc ! polynom coefficients
    real(FREAL), dimension(:), intent(in) :: vpt ! the x points
    real(FREAL), dimension(:), intent(out) :: ypt ! The y returned points
    !
    integer(FINT) :: i,j
    !
    f_poly=0
    !
    do i=1,npt
       ypt(i)=pc(1)
       do j=2,npc
          ypt(i)=ypt(i)+pc(j)*(vpt(i)**(j-1))
       end do
    end do
    !
!
!H
!H-----------------------------------------------------------------------------
!H
!
  end function f_poly

!
!H
!H-----------------------------------------------------------------------------
!H
!

!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function funct_polys(npar,grad,fval,xval,iflag)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function funct_polys(npar,grad,fval,xval,iflag)
!
!H
!H-----------------------------------------------------------------------------
!H
!H The actual polys function which is like poly but shift the X as parameter!
!H
!H-----------------------------------------------------------------------------
!H
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
    integer(FINT) :: i,irc,iun,j,k
    character(FLCHARS) :: ctmp
    !
    funct_polys = 0
    !
    if (iflag.eq.1) then
       ! Read input data: Here we do nothing
       !------------------------------------
       call message(MESOUT,'[funct_polys] at present nothing to do at function call 1.')
    else if (iflag.eq.2) then
       !
       ! Calculate GRAD, the first derivatives of FVAL.
       !-----------------------------------------------
       ! We don't have it now
       !-----------------------------------------------
       grad(:)=0.0_FREAL
       !
       !irc = f_poly(f1dp_numpar,xval(:),f1dp_numrec,f1dp_data(:,1),f1dp_ypol(:))
       !grad(1)=2.0_FREAL*sum(f1dp_ypol(:)-f1dp_data(:,2))
       !grad(2)=2.0_FREAL*sum((f1dp_ypol(:)-f1dp_data(:,2))*f1dp_data(:,1))
       !do i=3,f1dp_numpar
       !   grad(i)=2.0_FREAL*sum((f1dp_ypol(:)-f1dp_data(:,2))*(f1dp_data(:,1)**(i-1)))
       !end do
    end if
    !
    ! Here we evaluate the function value
    ! first parameter is shift
    !------------------------------------
    irc = f_poly(f1dp_numpar-1,xval(2:),f1dp_numrec,f1dp_data(1:f1dp_numrec,1)+xval(1),f1dp_ypol(1:f1dp_numrec))
    fval=sum((f1dp_ypol(1:f1dp_numrec)-f1dp_data(1:f1dp_numrec,2))**2)
    !
    ! The last iflag is evaluated here
    !---------------------------------
    if (iflag.eq.3) then
       ! Final calculations end print out
       !---------------------------------
       call message(MESOUT,"[funct_polys] Calculation Done: fcn call 3 requested.")
       write(f1dp_iunout,'(1X,A)') ' Final values from the fcn call'
       do i=1,f1dp_numpar
          write(f1dp_iunout,'(1X,"Par:",1X,A10,2X,F18.10,2X,"fix: ",A1)') f1dp_parnames(i),xval(i),f1dp_strfixed(i)
       end do
       write(f1dp_iunout,'(1X,"Final FunctionValue:",1X,F18.10)') fval
       !
       ! Check if gnuplot was a request
       !-------------------------------
       if (f1dp_flags(FLG_GPLT).gt.0) then
          ! set the gnuplot dta unit
          if(f1dp_flags(FLG_GPLT).eq.2) then
             iun=f1dp_iungpltdata
          else
             iun=f1dp_iunout
          endif
          !
          write(iun,'(A1,1X,A)') '#','Gnuplot data file starts here'
          write(iun,'(A1,1X,"numrec=",I6,2X,"numcol=",I6)') '#',f1dp_numrec,f1dp_numcol
          do j=1,f1dp_numrec
             k=0
             do k=1,f1dp_numcol-1
                write(iun,'(F18.10,$)') f1dp_data(j,k)
             end do
             write(iun,'(F18.10)') f1dp_data(j,f1dp_numcol)
          end do
          ! set the gnuplot file unit
          if(f1dp_flags(FLG_GPLT).eq.2) then
             iun=f1dp_iungpltgplt
          else
             iun=f1dp_iunout
          endif
          !
          write(iun,'(A1,1X,A)') '#','Gnuplot autogenerated file starts here'
          write(iun,'("set style data linespoints")')
          write(iun,'(A1,1X,A)') '#','Parameters:'
          do i=1,f1dp_numpar
             write(iun,'(A10,"=",F18.10 )') f1dp_parnames(i),xval(i)
          end do
          write(iun,'(A1,1X,A)') '#','Function:'
          write(iun,'(A,$)') 'myfunct(x)='
          write(iun,'(A10,"+",$)') f1dp_parnames(2)
          do i=3,f1dp_numpar-1
             j=i-2
             write(iun,'(A10,"*((x+",A10,")**(",I6,")) + ",$)') f1dp_parnames(i),f1dp_parnames(1),j
          end do
          j=f1dp_numpar-2
          write(iun,'(A10,"*((x+",A10,")**(",I6,"))")') f1dp_parnames(f1dp_numpar),f1dp_parnames(1),j
          !
          j=len_trim(f1dp_gpltprefix)
          ctmp=f1dp_gpltprefix(1:j)//'.dat'
          ctmp=trim(ctmp)
          j=len_trim(ctmp)
          if(f1dp_flags(FLG_GPLT).eq.2) then
             write(iun,'("plot [",F12.6,":",F12.6,"] ",$)') f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             write(iun,'(A,$)') ctmp(1:j)
          else
             write(iun,'(A)') '# This would be the plot command'
             write(iun,'("# plot [",F12.6,":",F12.6,"] ",$)') f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             write(iun,'(A,$)') 'yourdatafile'
          end if
          write(iun,'(A1,$)') "'"
          write(iun,'(",myfunct(x)")')
          write(iun,'("pause -1")')
          !
       end if
       !
    end if
    !
    return
    !
  end function funct_polys
!
!H
!H-----------------------------------------------------------------------------
!H
!

!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function funct_morse(npar,grad,fval,xval,iflag)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function funct_morse(npar,grad,fval,xval,iflag)
!
!H
!H-----------------------------------------------------------------------------
!H
!H The acctual simple morse function
!H
!H-----------------------------------------------------------------------------
!H
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
    integer(FINT), parameter :: LOC_MAXPAR=3
    ! WARNING! the section iflag=3 depend on this too!
    !
    integer(FINT) :: i,irc,iun,j,k
    character(FLCHARS) :: ctmp
    !
    funct_morse = 0
    !
    if (iflag.eq.1) then
       ! Read input data: Here we do nothing
       !------------------------------------
       if(npar.gt.LOC_MAXPAR.or.f1dp_numpar.gt.LOC_MAXPAR) then
          call message(MESERRO,'[funct_morse] This simple morse fit accept only 3 parameters.')
          call message_value(MESERRO,'[funct_morse] Your parameters at input are:',f1dp_numpar)
          call message_value(MESERRO,'[funct_morse] Your parameters at minuit are:',npar)
          call message(MESKILL,'[funct_morse] We stop now.')
       end if
       call message(MESOUT,'[funct_morse] at present nothing to do at function call 1 beside a par check.')
    else if (iflag.eq.2) then
       !
       ! Calculate GRAD, the first derivatives of FVAL.
       !-----------------------------------------------
       !-----------------------------------------------
       grad(:)=0.0_FREAL
       !
    end if
    !
    ! Here we evaluate the function value
    !------------------------------------
    !morp(x,vd,va) =  vd * (1 - exp(-va*(x-re)))**2
    f1dp_ypol(1:f1dp_numrec)=xval(1) * ((1.0_FREAL - exp(-xval(2)*(f1dp_data(1:f1dp_numrec,1)-xval(3))))**2)
    fval=sum((f1dp_ypol(1:f1dp_numrec)-f1dp_data(1:f1dp_numrec,2))**2)
    !
    ! The last iflag is evaluated here
    !---------------------------------
    if (iflag.eq.3) then
       ! Final calculations end print out
       !---------------------------------
       call message(MESOUT,"[funct_morse] Calculation Done: fcn call 3 requested.")
       write(f1dp_iunout,'(1X,A)') ' Final values from the fcn call'
       do i=1,f1dp_numpar
          write(f1dp_iunout,'(1X,"Par:",1X,A10,2X,F18.10,2X,"fix: ",A1)') f1dp_parnames(i),xval(i),f1dp_strfixed(i)
       end do
       write(f1dp_iunout,'(1X,"Final FunctionValue:",1X,F18.10)') fval
       !
       ! Check if gnuplot was a request
       !-------------------------------
       if (f1dp_flags(FLG_GPLT).gt.0) then
          ! set the gnuplot dta unit
          if(f1dp_flags(FLG_GPLT).eq.2) then
             iun=f1dp_iungpltdata
          else
             iun=f1dp_iunout
          endif
          !
          write(iun,'(A1,1X,A)') '#','Gnuplot data file starts here'
          write(iun,'(A1,1X,"numrec=",I6,2X,"numcol=",I6)') '#',f1dp_numrec,f1dp_numcol
          do j=1,f1dp_numrec
             k=0
             do k=1,f1dp_numcol-1
                write(iun,'(F18.10,$)') f1dp_data(j,k)
             end do
             write(iun,'(F18.10)') f1dp_data(j,f1dp_numcol)
          end do
          ! set the gnuplot file unit
          if(f1dp_flags(FLG_GPLT).eq.2) then
             iun=f1dp_iungpltgplt
          else
             iun=f1dp_iunout
          endif
          !
          write(iun,'(A1,1X,A)') '#','Gnuplot autogenerated file starts here'
          write(iun,'("set style data linespoints")')
          write(iun,'(A1,1X,A)') '#','Parameters:'
          do i=1,f1dp_numpar
             write(iun,'(A10,"=",F18.10 )') f1dp_parnames(i),xval(i)
          end do
          write(iun,'(A1,1X,A)') '#','Function:'
          !morp(x,vd,va) =  vd * (1 - exp(-va*(x-re)))**2
          !
          ! We do it a bit less portable for aestetical matter
          !---------------------------------------------------
          write(iun,'(1X,"lvd=",F18.10)') xval(1)
          write(iun,'(1X,"lva=",F18.10)') xval(2)
          write(iun,'(1X,"lre=",F18.10)') xval(3)
          write(iun,'(A)') 'mymorse(x,lvd,lva)=  lvd * (1 - exp(-lva*(x-lre)))**2'
          !
          ! The data file name
          !-------------------
          j=len_trim(f1dp_gpltprefix)
          ctmp=f1dp_gpltprefix(1:j)//'.dat'
          ctmp=trim(ctmp)
          j=len_trim(ctmp)
          !
          ! The plot command
          !-----------------
          if(f1dp_flags(FLG_GPLT).eq.2) then
             write(iun,'("plot [",F12.6,":",F12.6,"] ",$)') f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             write(iun,'(A,$)') ctmp(1:j)
          else
             write(iun,'(A)') '# This would be the plot command'
             write(iun,'("# plot [",F12.6,":",F12.6,"] ",$)') f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             write(iun,'(A,$)') 'yourdatafile'
          end if
          write(iun,'(A1,$)') "'"
          write(iun,'(",mymorse(x,lvd,lva)")')
          write(iun,'("pause -1")')
          !
          ! The eventual fit check
          !-----------------------
          write(iun,'("fit mymorse(x,lvd,lva) ",A1,A,A1," via lvd,lva,lre")') "'",ctmp(1:j),"'"
          !
          ! The other fit plot
          !-------------------
          if(f1dp_flags(FLG_GPLT).eq.2) then
             write(iun,'("plot [",F12.6,":",F12.6,"] ",$)') f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             j=len_trim(ctmp)
             write(iun,'(A,$)') ctmp(1:j)
          else
             write(iun,'(A)') '# This would be the plot command'
             write(iun,'("# plot [",F12.6,":",F12.6,"] ",$)') f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             write(iun,'(A,$)') 'yourdatafile'
          end if
          write(iun,'(A1,$)') "'"
          write(iun,'(",mymorse(x,lvd,lva)")')
          write(iun,'("pause -1")')
          !
       end if
       !
    end if
    !
    return
    !
  end function funct_morse
!
!H
!H-----------------------------------------------------------------------------
!H
!


!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function funct_polyd(npar,grad,fval,xval,iflag)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function funct_polyd(npar,grad,fval,xval,iflag)
!
!H
!H-----------------------------------------------------------------------------
!H
!H The acctual polyd function: multidim fit
!H
!H-----------------------------------------------------------------------------
!H
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
    integer(FINT) :: i,irc,iun,j,k
    integer(FINT) :: mpar,mdim,dpar,ddim,matrdim
    character(FLCHARS) :: ctmp
    !
    ! The expansion data
    integer(FINT), dimension(:,:), allocatable,save :: coordexp
    integer(FINT) :: il, ib,j1,j2,j3,j4,j5,j6
    real(FREAL) :: rtmp
    !
    funct_polyd = 0
    !
    if (iflag.eq.1) then
       ! Read input data: Here we do nothing
       !------------------------------------
       call message(MESOUT,'[funct_polyd] at present nothing to do at function call 1.')
       mdim=f1dp_flags(FLG_POLO)
       ddim=f1dp_numcol-1
       mpar=(2**mdim)-1
       dpar=(2**ddim)-1
       call message_value(MESOUT,'[POLYD]          Problem Dimension: ',ddim)
       call message_value(MESOUT,'[POLYD] Max number of combinations: ',dpar)
       call message_value(MESOUT,'[POLYD] Dimension or Order of exp.: ',mdim)
       call message_value(MESOUT,'[POLYD]   Max number of parameters: ',mpar)
       !
       ! allocation data
       !Counts:           1           6
       !Counts:           2          27
       !Counts:           3          83
       !Counts:           4         209
       !Counts:           5         461
       !Counts:           6         923
       !Counts:           7        1715
       !Counts:           8        3002
       !Counts:           9        5004
       !Counts:          10        8007
       !Counts:          11       12375
       !Counts:          12       18563
       !
       ! ADD 1 to that number!
       !
       if(.not.allocated(coordexp)) then
          allocate(coordexp(6,924),STAT=irc)
          if(.not.irc.eq.0) then
             call message(MESERRO,'[fit1Dpol_polyd] Cannot allocate Exponents')
             call message(MESKILL,'[fit1Dpol_polyd] We stop now')
          end if
       endif
       !
       ! We warn because the method works only up to 6 at present
       !---------------------------------------------------------
       if(mdim.gt.6.or.ddim.ne.6) then
          call message(MESERRO,'[fit1Dpol_polyd] Max Expantion order is 6!')
          call message(MESKILL,'[fit1Dpol_polyd] We stop now.')
       end if
       !
       ! Get expansion
       !--------------
       coordexp(1,1) = 0
       coordexp(2,1) = 0
       coordexp(3,1) = 0
       coordexp(4,1) = 0
       coordexp(5,1) = 0
       coordexp(6,1) = 0
       ib=1
       do il = 1, mdim
          do j1 = 0, il
          do j2 = 0, j1
          do j3 = 0, j2
          do j4 = 0, j3
          do j5 = 0, j4
             ib=ib+1
             coordexp(1,ib) = il-j1
             coordexp(2,ib) = j1-j2
             coordexp(3,ib) = j2-j3
             coordexp(4,ib) = j3-j4
             coordexp(5,ib) = j4-j5
             coordexp(6,ib) = j5
          end do
          end do
          end do
          end do
          end do
       end do

       if(f1dp_debug.gt.0) then
          call message(MESOUT,"Here are the parameters names to be used")
          call message(MESOUT,"With this very same order.")
          write(f1dp_iunout,'(1X,A6,1X,A1,6(I1),2X,I6,2(2X,F18.8),1X,A1)')&
               &"CUTME0","f",0,0,0,0,0,0,1,0.0_FREAL,0.0_FREAL,"n"
          ib=1
          do il = 1, mdim
             do j1 = 0, il
                do j2 = 0, j1
                   do j3 = 0, j2
                      do j4 = 0, j3
                         do j5 = 0, j4
                            ib=ib+1
                            write(f1dp_iunout,'(1X,A6,1X,A1,6(I1),2X,I6,2(2X,F18.8),1X,A1)')&
                                 &"CUTME0","f",coordexp(:,ib),ib,0.0_FREAL,0.0_FREAL,"n"
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end if
       !
    else if (iflag.eq.2) then
       !
       ! Calculate GRAD, the first derivatives of FVAL.
       !-----------------------------------------------
       !-----------------------------------------------
       grad(:)=0.0_FREAL
       !
       !irc = f_poly(f1dp_numpar,xval(1:f1dp_numpar),f1dp_numrec,f1dp_data(1:f1dp_numrec,1),f1dp_ypol(1:f1dp_numrec))
       !grad(1)=2.0_FREAL*sum(f1dp_ypol(1:f1dp_numrec)-f1dp_data(1:f1dp_numrec,2))
       !grad(2)=2.0_FREAL*sum((f1dp_ypol(1:f1dp_numrec)-f1dp_data(1:f1dp_numrec,2))*f1dp_data(1:f1dp_numrec,1))
       !do i=3,f1dp_numpar
       !   grad(i)=2.0_FREAL*sum((f1dp_ypol(1:f1dp_numrec)-f1dp_data(1:f1dp_numrec,2))*(f1dp_data(1:f1dp_numrec,1)**(i-1)))
       !end do
    end if
    !
    ! Here we evaluate the function value
    !------------------------------------
    !irc = f_poly(f1dp_numpar,xval(1:f1dp_numpar),f1dp_numrec,f1dp_data(1:f1dp_numrec,1),f1dp_ypol(1:f1dp_numrec))
    f1dp_ypol(:)=0.0_FREAL
    do ib=1,f1dp_numrec
       f1dp_ypol(ib)=xval(1)
       do k=2,f1dp_numpar
          i=k-1
          rtmp = 1.0_FREAL
          do j=1,6
             if(coordexp(j,k).gt.0) then
                rtmp=rtmp*f1dp_data(ib,j)**coordexp(j,k)
             end if
          end do
          f1dp_ypol(ib)=f1dp_ypol(ib)+xval(k)*rtmp
       end do
       f1dp_ypol(ib)=f1dp_ypol(ib)+xval(1)
    end do
    fval=sum((f1dp_ypol(1:f1dp_numrec)-f1dp_data(1:f1dp_numrec,7))**2)
    !
    ! The last iflag is evaluated here
    !---------------------------------
    if (iflag.eq.3) then
       ! Final calculations end print out
       !---------------------------------
       call message(MESOUT,&
            &"[funct_polyd] Calculation Done: fcn call 3 requested.")
       write(f1dp_iunout,'(1X,A)') ' Final values from the fcn call'
       do i=1,f1dp_numpar
          write(f1dp_iunout,'(1X,"Par:",1X,A10,2X,F18.10,2X,"fix: ",A1)')&
               &f1dp_parnames(i),xval(i),f1dp_strfixed(i)
       end do
       write(f1dp_iunout,'(1X,"Final FunctionValue:",1X,F18.10)') fval
       !
       ! Check if gnuplot was a request
       !-------------------------------
       if (f1dp_flags(FLG_GPLT).gt.0) then
          ! set the gnuplot dta unit
          if(f1dp_flags(FLG_GPLT).eq.2) then
             iun=f1dp_iungpltdata
          else
             iun=f1dp_iunout
          endif
          !
          write(iun,'(A1,1X,A)') '#','Gnuplot data file starts here'
          write(iun,'(A1,1X,"numrec=",I6,2X,"numcol=",I6)')&
               &'#',f1dp_numrec,f1dp_numcol
          do j=1,f1dp_numrec
             k=0
             do k=1,f1dp_numcol-1
                write(iun,'(F18.10,$)') f1dp_data(j,k)
             end do
             write(iun,'(F18.10)') f1dp_data(j,f1dp_numcol)
          end do
          ! set the gnuplot file unit
          if(f1dp_flags(FLG_GPLT).eq.2) then
             iun=f1dp_iungpltgplt
          else
             iun=f1dp_iunout
          endif
          !
          write(iun,'(A1,1X,A)') '#','Gnuplot autogenerated file starts here'
          write(iun,'("set style data linespoints")')
          write(iun,'(A1,1X,A)') '#','Parameters:'
          do i=1,f1dp_numpar
             write(iun,'(A10,"=",F18.10 )') f1dp_parnames(i),xval(i)
          end do
          write(iun,'(A1,1X,A)') '#','Function:'
          write(iun,'(A,$)') 'myfunct(x)='
          write(iun,'(A10,"+",$)') f1dp_parnames(1)
          do i=2,f1dp_numpar-1
             j=i-1
             write(iun,'(A10,"*(x**(",I6,")) + ",$)') f1dp_parnames(i),j
          end do
          j=f1dp_numpar-1
          write(iun,'(A10,"*(x**(",I6,"))")') f1dp_parnames(f1dp_numpar),j
          !
          j=len_trim(f1dp_gpltprefix)
          ctmp=f1dp_gpltprefix(1:j)//'.dat'
          ctmp=trim(ctmp)
          j=len_trim(ctmp)
          if(f1dp_flags(FLG_GPLT).eq.2) then
             write(iun,'("plot [",F12.6,":",F12.6,"] ",$)')&
                  &f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             write(iun,'(A,$)') ctmp(1:j)
          else
             write(iun,'(A)') '# This would be the plot command'
             write(iun,'("# plot [",F12.6,":",F12.6,"] ",$)')&
                  &f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             write(iun,'(A,$)') 'yourdatafile'
          end if
          write(iun,'(A1,$)') "'"
          write(iun,'(",myfunct(x)")')
          write(iun,'("pause -1")')
          write(iun,'("fit myfunct(x) ",A1,A,A1," via ",$)') "'",ctmp(1:j),"'"
          do i=1,f1dp_numpar-1
             j=len_trim(f1dp_parnames(i))
             write(iun,'(A,",",$)') f1dp_parnames(i)(1:j)
          end do
          j=len_trim(f1dp_parnames(f1dp_numpar))
          write(iun,'(A)') f1dp_parnames(f1dp_numpar)(1:j)
          !
          if(f1dp_flags(FLG_GPLT).eq.2) then
             write(iun,'("plot [",F12.6,":",F12.6,"] ",$)')&
                  &f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             j=len_trim(ctmp)
             write(iun,'(A,$)') ctmp(1:j)
          else
             write(iun,'(A)') '# This would be the plot command'
             write(iun,'("# plot [",F12.6,":",F12.6,"] ",$)')&
                  &f1dp_data(1,1),f1dp_data(f1dp_numrec,1)
             write(iun,'(A1,$)') "'"
             write(iun,'(A,$)') 'yourdatafile'
          end if
          write(iun,'(A1,$)') "'"
          write(iun,'(",myfunct(x)")')
          write(iun,'("pause -1")')
          !
       end if
       !
    end if
    !
    return
    !
  end function funct_polyd
!
!H
!H-----------------------------------------------------------------------------
!H
!



end module f1dp_minuit_function
