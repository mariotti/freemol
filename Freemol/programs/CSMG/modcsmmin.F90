!
!H
!H-----------------------------------------------------------------------------
!H This File is part of the Frimol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H MODULE modcsmmin  Frimol by F.Mariotti: (c) F.Mariotti
!H-----------------------------------------------------------------------------
!H $Id: modcsmmin.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
!H-----------------------------------------------------------------------------
!H 
!H DESCRIPTION:
!H 
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
module modcsmmin
  !
  use vartypes
  use messages
  use molecule
  use csmmod
  !
  implicit none
  private
  !
  ! PUBLIC DECLARATIONS:
  !--------------------
  
  ! Public Procedures
  !------------------
  public modcsmmin_init
  public modcsmmin_isinit
  public modcsmmin_geterror
  public modcsmmin_restore
  public modcsm_nmopt
  !public mcsm_rotxyz
  public modcsmmin_setcoord
  public modcsmmin_minuit_value
  !public modcsmmin_
  !
  !
  ! Error Interface
  !----------------
  ! We define two error interfaces. Each one return the latest error
  ! code or string value. The String error function can be used to inqure
  ! about a given error code.
  !----------------------------------------------------------------------
  interface modcsmmin_geterror
     module procedure modcsmmin_geterror_i
     module procedure modcsmmin_geterror_c
  end interface !modcsmmin_geterror
  !
  real(FREAL), public, save :: modcsmmin_lastfunctionvalue
  !
  ! LOCAL DECLARATIONS:
  !--------------------
  !
  logical, save :: lcmin_isinit = .false.
  !
  ! error variables
  !----------------
  integer(FINT), save, private :: lcmin_error = 0
  character(FLCHARS), save, private :: lcmin_error_message = ""
  !
  ! Varia
  !------
  integer(FINT) :: irc, iunout
  logical, save :: verbose = .false.
  !
  ! Saved Coordinates
  !------------------
  real(FREAL), save, allocatable, public, dimension(:,:) :: saved_xyz
  !
  ! CSM stuff
  !----------
  integer(FINT), save :: nat, lnsymop
  real(FREAL), save :: srhorho, csmnorm
  real(FREAL), allocatable, dimension(:,:,:) :: impj
  !
  !
contains
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function modcsmmin_init(liunout,srr,cn,msymop,lver)
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function modcsmmin_init(liunout,srr,cn,msymop,lver)
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the modcsmmin module.
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    real(FREAL), intent(in) :: srr, cn
    integer(FINT), intent(in) :: msymop, liunout
    logical :: lver
    !
    lcmin_isinit = .true.
    iunout = liunout
    verbose = lver
    !
    ! We first save the original coordiantes
    !---------------------------------------
    lnsymop = msymop
    nat = molecule_getnat()
    allocate(saved_xyz(3,nat),STAT=irc)
    saved_xyz(1:3,1:nat) = molecule_xyz(1:3,1:nat)
    srhorho = srr
    csmnorm = cn
    allocate(impj(nat,lnsymop,nat),STAT=irc)
    !
    modcsmmin_init = 1
    !
  end function modcsmmin_init
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function modcsmmin_restore()
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function modcsmmin_restore()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the modcsmmin module.
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    ! We restore the original coordinates
    !------------------------------------
    nat = molecule_getnat()
    molecule_xyz(1:3,1:nat) = saved_xyz(1:3,1:nat)
    modcsmmin_restore = 1
    !
  end function modcsmmin_restore
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function modcsmmin_setcoord()
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function modcsmmin_setcoord()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the modcsmmin module.
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    ! We restore the original coordinates
    !------------------------------------
    nat = molecule_getnat()
    saved_xyz(1:3,1:nat) = molecule_xyz(1:3,1:nat)
    modcsmmin_setcoord = 1
    !
  end function modcsmmin_setcoord
!
!H
!H-----------------------------------------------------------------------------
!H logical function modcsmmin_isinit()
!H-----------------------------------------------------------------------------
!H
!
  logical function modcsmmin_isinit()
!
!H
!H-----------------------------------------------------------------------------
!H We merelyreturn the initialization flag.
!H-----------------------------------------------------------------------------
!H
!
    !
    modcsmmin_isinit = lcmin_isinit
    !
  end function modcsmmin_isinit
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function modcsmmin_geterror_i()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function modcsmmin_geterror_i()
!
!H
!H-----------------------------------------------------------------------------
!H 
!H Return the error code of the module.
!H 
!H-----------------------------------------------------------------------------
!H
!
    ! We merely return the integer current error code
    !------------------------------------------------
    modcsmmin_geterror_i = lcmin_error
    !
  end function modcsmmin_geterror_i
!
!H-----------------------------------------------------------------------------
!H
!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H character(FLCHARS) function modcsmmin_geterror_c()
!H-----------------------------------------------------------------------------
!H
!
  character(FLCHARS) function modcsmmin_geterror_c(code)
    !
    integer(FINT), intent(in) :: code
    !
!
!H
!H-----------------------------------------------------------------------------
!H 
!H Return the error code of the module in Object oriented
!H programming style. It is a part of a set of routine to return errors.
!H See the discussion in the interface.
!H 
!H-----------------------------------------------------------------------------
!H
!
    !
    ! We return the error message
    !----------------------------
    ! Note that we have not yet implement it
    ! so we report as example case 2_FINT as unwanted code
    ! but it has to disappear! And Default is a not reported
    ! error code.
    !-------------------------------------------------------
    select case(code)
    case(0_FINT)
       modcsmmin_geterror_c = "No Errors."
    case(1_FINT)
       modcsmmin_geterror_c = "Not documented error or generic error.."
    case(2_FINT)
       modcsmmin_geterror_c = "Double Error: Why did you get this message?."
    case default
       modcsmmin_geterror_c = "Not documented error or generic error.."
    end select
    !
    return
    !
  end function modcsmmin_geterror_c
!
!H
!H-----------------------------------------------------------------------------
!H
!H-----------------------------------------------------------------------------
!H
!
!  function 
!
!H
!H-----------------------------------------------------------------------------
!H
!
!
!H
!H-----------------------------------------------------------------------------
!H
!







!                                                                       
!-----------------------------------------------------------------------
!                                                                       
!
!H
!H-----------------------------------------------------------------------------
!H subroutine modcsm_nmopt(x,xm,ftol,fret,iter)
!H-----------------------------------------------------------------------------
!H
!H
!H-----------------------------------------------------------------------------
!H
!H This comes from a C.Daul Routine ... hope so.. at least ...
!H Should I do have a look to NR as well??
!H 
!H 
!H
!H-----------------------------------------------------------------------------
!H
!
!  subroutine modcsm_nmopt(funct,ndim,x,xm,ftol,fret,iter) 
  subroutine modcsm_nmopt(x,xm,ftol,fret,iter) 
    !
    integer(FINT), intent(out) :: iter
    real(FREAL), dimension(:), intent(inout) :: x
    real(FREAL), dimension(:), intent(in) :: xm
    real(FREAL), intent(in) :: ftol
    real(FREAL), intent(out) :: fret
    !dimension x(ndim),xm(ndim) 
    !
    ! Needed locals: seems to be!
    !----------------------------
    integer(FINT) :: ndim
    real(FREAL), allocatable, dimension(:) :: y,pr,prr,pbar
    real(FREAL), allocatable, dimension(:,:) :: p
    !dimension p(mmax,nmax),y(mmax),pr(nmax),prr(nmax),pbar(nmax) 
    !
    !parameter (nmax=6,mmax=nmax+1,alpha=1.0,beta=0.5,gamma=2.0,       &
    !     &           itmax=500,rltiny=1.e-7)                                  
    !parameter (nmax=6,mmax=nmax+1,alpha=1.0,beta=0.5,gamma=2.0,       &
    !     &           itmax=500,rltiny=1.e-7)                                  
    !
    ! Some pars
    !----------
    real(FREAL)  , parameter :: zero  = 0.0_FREAL
    real(FREAL)  , parameter :: half  = 0.5_FREAL
    real(FREAL)  , parameter :: one   = 1.0_FREAL
    real(FREAL)  , parameter :: two   = 2.0_FREAL
    real(FREAL)  , parameter :: rltiny  = 1.0E-8_FREAL
    integer(FINT), parameter :: itmax = 500
    !
    integer(FINT) :: nmax
    integer(FINT) :: mmax
    integer(FINT) :: mpts
    integer(FINT) :: i
    integer(FINT) :: j
    integer(FINT) :: ilo
    integer(FINT) :: ihi
    integer(FINT) :: inhi
    real(FREAL) :: rtol
    real(FREAL), parameter :: alpha = 1.0_FREAL
    real(FREAL) :: ypr
    real(FREAL), parameter :: gamma = 2.0_FREAL
    real(FREAL) :: yprr
    real(FREAL), parameter :: beta = 0.5_FREAL
    !
    !-----------------------------------------------------------------------
    !                                                                       
    !  minimization of a function using the downhill simplex of nelder and m
    ! ======================================================================
    !                                                                       
    !                                                                       
    !  funct: name of the routine used to calculate the function to be minim
    !  ndim : number of parameters                                          
    !  x(j) : parameter j (initial guesses on input, optimized values on out
    !  xm(j): characteristic step length of parameter j                     
    !  ftol : fractional convergence tolerance                              
    !  fret : optimal value of function at the local minimum                
    !  iter : number of iterations                                          
    !                                                                       
    !----------------------------------------------------------------------
    !
    ! We get the number of parameters
    !--------------------------------
    ndim = size(x,1)
    nmax = ndim
    mmax = ndim + 1
    !
    !dimension p(mmax,nmax),y(mmax),pr(nmax),prr(nmax),pbar(nmax) 
    !
    ! This stuff allocation
    !----------------------
    if(.not.allocated(p)) allocate(p(mmax,nmax),STAT=irc)
    if(.not.allocated(y)) allocate(y(mmax),STAT=irc)
    if(.not.allocated(pr)) allocate(pr(nmax),STAT=irc)
    if(.not.allocated(prr)) allocate(prr(nmax),STAT=irc)
    if(.not.allocated(pbar)) allocate(pbar(nmax),STAT=irc)
    !
    ! Reset these parameters
    !-----------------------
    p(:,:) = zero
    y(:) = zero
    pr(:) = zero
    prr(:) = zero
    pbar(:) = zero
    !
    !if(ndim.gt.nmax) pause 'too many parameters. increase dimension' 
    !
    mpts=ndim+1 
    !                                                                       
    !  initialize the simplex p(i,j)                                        
    !                                                                       
    do i=1,mpts 
       do j=1,ndim 
          p(i,j)=x(j)
          if(j.eq.i-1) p(i,j)=p(i,j)+xm(j) 
       enddo
    enddo
    !                                                                       
    !  get y(i) at vertices of simplex                                      
    !                                                                       
    do i=1,mpts 
       do j=1,ndim 
          x(j)=p(i,j) 
       enddo
       call funct(ndim,x,y(i)) 
    enddo
    !                                                                       
    iter=0 
1   ilo=1 
    if(y(1).gt.y(2))then 
       ihi=1 
       inhi=2 
    else 
       ihi=2 
       inhi=1 
    endif
    do i=1,mpts 
       if(y(i).lt.y(ilo)) ilo=i 
       if(y(i).gt.y(ihi))then 
          inhi=ihi 
          ihi=i 
       else if(y(i).gt.y(inhi))then 
          if(i.ne.ihi) inhi=i 
       endif
    enddo
    !                                                                       
    !
    rtol=two*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+rltiny) 
    if(verbose) then
       fret=zero 
       x(:)=zero
       do i=1,mpts 
          fret=fret+y(i) 
          do j=1,ndim 
             x(j)=x(j)+p(i,j) 
          end do
       end do
       fret=fret/mpts 
       do j=1,ndim 
          x(j)=x(j)/mpts 
       end do
       write(iunout,&
            &'("# Minimization at step:",1X,I5,1X,&
            &"with Value:",1X,F10.4,"Tol(Act/Req):",E10.4,"/",E10.4)')&
            &iter,fret,rtol,ftol
       write(iunout,'("# With parameters(tra.):",1X,3(2X,E12.4))') x(1:3)
       write(iunout,'("# With parameters(rot.):",1X,3(2X,E12.4))') x(4:6)
    end if
    if(rtol.lt.ftol) then 
       fret=zero 
       do j=1,ndim 
          x(j)=zero 
       end do
       do i=1,mpts 
          fret=fret+y(i) 
          do j=1,ndim 
             x(j)=x(j)+p(i,j) 
          end do
       end do
       fret=fret/mpts 
       do j=1,ndim 
          x(j)=x(j)/mpts 
       end do
       return 
    endif
    if(iter.eq.itmax) stop 'modcsm_nmopt exceeding maximum iterations.' 
    iter=iter+1 
    do j=1,ndim 
       pbar(j)=zero 
    end do
    do i=1,mpts 
       if(i.ne.ihi)then 
          do j=1,ndim 
             pbar(j)=pbar(j)+p(i,j) 
          end do
       end if
    end do
    do j=1,ndim 
       pbar(j)=pbar(j)/ndim 
       pr(j)=(one+alpha)*pbar(j)-alpha*p(ihi,j) 
    end do
    call funct(ndim,pr,ypr) 
    if(ypr.le.y(ilo))then 
       do j=1,ndim 
          prr(j)=gamma*pr(j)+(1.-gamma)*pbar(j) 
       end do
       call funct(ndim,prr,yprr) 
       if(yprr.lt.y(ilo))then 
          do j=1,ndim 
             p(ihi,j)=prr(j) 
          end do
          y(ihi)=yprr 
       else 
          do j=1,ndim 
             p(ihi,j)=pr(j) 
          end do
          y(ihi)=ypr 
       endif
    else if(ypr.ge.y(inhi))then 
       if(ypr.lt.y(ihi))then 
          do j=1,ndim 
             p(ihi,j)=pr(j) 
          end do
          y(ihi)=ypr 
       endif
       do j=1,ndim 
          prr(j)=beta*p(ihi,j)+(one-beta)*pbar(j) 
       end do
       call funct(ndim,prr,yprr) 
       if(yprr.lt.y(ihi))then 
          do j=1,ndim 
             p(ihi,j)=prr(j) 
          end do
          y(ihi)=yprr 
       else 
          do i=1,mpts 
             if(i.ne.ilo)then 
                do j=1,ndim 
                   pr(j)=half*(p(i,j)+p(ilo,j)) 
                   p(i,j)=pr(j) 
                end do
                call funct(ndim,pr,y(i)) 
             endif
          end do
       endif
    else 
       do j=1,ndim 
          p(ihi,j)=pr(j) 
       end do
       y(ihi)=ypr 
    endif
    go to 1 
  end subroutine modcsm_nmopt
  
  








!
!H
!H-----------------------------------------------------------------------------
!H subroutine mcsm_rotxyz(tx,ty,tz,rx,ry,rz)
!H-----------------------------------------------------------------------------
!H
!
  subroutine mcsm_rotxyz(tx,ty,tz,rx,ry,rz)
!
!H
!H-----------------------------------------------------------------------------
!H
!H tx,ty,tz is the translation vector
!H rx,ry,rz are the angles defining the amount the object is rotated
!H around x,y,z axes respectevly.
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    !
    real(FREAL), intent(in) :: tx,ty,tz,rx,ry,rz
    !
    integer(FINT) :: i,nat
    real(FREAL), dimension(3,3) :: rmat
    real(FREAL), dimension(3,3) :: amat,bmat,gmat
    real(FREAL), dimension(3) :: taxe
    real(FREAL) :: eps = 1.0E-12_FREAL
    real(FREAL) :: alpha,beta,gamma
    !
    ! Translate First
    !----------------
    molecule_xyz(1,:) = tx + saved_xyz(1,:)
    molecule_xyz(2,:) = ty + saved_xyz(2,:)
    molecule_xyz(3,:) = tz + saved_xyz(3,:)
    !
    ! Get rotational matrixes
    !------------------------
    !taxe(1:3) = (/1,0,0/)
    !call csmmod_euler(rx,taxe(1:3),.true.,alpha,beta,gamma,amat(1:3,1:3),eps)
    !taxe(1:3) = (/0,1,0/)
    !call csmmod_euler(ry,taxe(1:3),.true.,alpha,beta,gamma,bmat(1:3,1:3),eps)
    !taxe(1:3) = (/0,0,1/)
    !call csmmod_euler(rz,taxe(1:3),.true.,alpha,beta,gamma,gmat(1:3,1:3),eps)
    !
    taxe(1:3) = (/1,0,0/)
    call csmmod_euler(rx,taxe(1:3),.false.,alpha,beta,gamma,amat(1:3,1:3),eps)
    taxe(1:3) = (/0,1,0/)
    call csmmod_euler(ry,taxe(1:3),.false.,alpha,beta,gamma,bmat(1:3,1:3),eps)
    taxe(1:3) = (/0,0,1/)
    call csmmod_euler(rz,taxe(1:3),.false.,alpha,beta,gamma,gmat(1:3,1:3),eps)
    !
    rmat(1:3,1:3) = matmul(bmat,amat)
    rmat(1:3,1:3) = matmul(gmat,rmat)
    !
    ! Apply rotation
    !---------------
    nat = molecule_getnat()
    !call message(MESOUT,"INSIDE ROTATION")
    !call molecule_print(iunout,'mldfromau')
    do i=1,nat
       molecule_xyz(1:3,i) = matmul(rmat(1:3,1:3),molecule_xyz(1:3,i))
    end do
    !call molecule_print(iunout,'mldfromau')
    !
    !
    !
  end subroutine mcsm_rotxyz

!
!H
!H-----------------------------------------------------------------------------
!H subroutine funct(n,p,value)
!H-----------------------------------------------------------------------------
!H
!
  subroutine funct(n,p,value)
!
!H
!H-----------------------------------------------------------------------------
!H
!H 
!H 
!H 
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    integer(FINT), intent(in) :: n
    real(FREAL), intent(in), dimension(:) :: p
    real(FREAL), intent(out) :: value
    !
    if(n.ne.6) then
       call message(MESERRO,"[modcsm_funct] internal inconsistency.")
       call message(MESKILL,"[modcsm_funct] We stop now!")
    end if
    !
    ! Set the coordinates to the papameters
    !--------------------------------------
    call mcsm_rotxyz(p(1),p(2),p(3),p(4),p(5),p(6))
    !
    ! Call to get values for CSM
    !---------------------------
    call csmmod_getimpj(impj)
    value = ((real(lnsymop,FREAL)-1.0_FREAL)/(real(lnsymop,FREAL))) - ((sum(impj(:,:,:)))/(real(lnsymop,FREAL)*srhorho))
    !
  end subroutine funct
  !
  !

!
!H
!H-----------------------------------------------------------------------------
!H real(FREAL) function minuit_value(xvals)
!H-----------------------------------------------------------------------------
!H
!
  real(FREAL) function modcsmmin_minuit_value(xvals)
    !
    real(FREAL), intent(in), dimension(6) :: xvals
    real(FREAL) :: eps = 1.0E-12_FREAL
    !
    ! Set the coordinates to the papameters
    !--------------------------------------
    call mcsm_rotxyz(xvals(1),xvals(2),xvals(3),xvals(4),xvals(5),xvals(6))
    !
    ! Call to get values for CSM
    !---------------------------
    call csmmod_getimpj(impj)
    modcsmmin_minuit_value = &
         &((real(lnsymop,FREAL)-1.0_FREAL)/(real(lnsymop,FREAL))) - ((sum(impj(:,:,:)))/(real(lnsymop,FREAL)*srhorho))
    !
    modcsmmin_lastfunctionvalue = modcsmmin_minuit_value
    if (modcsmmin_minuit_value.lt.0.0_FREAL) then
       call message(MESWARN,"[munuit_value]: NEGATIVE FUNCTION VALUE.")
       call message(MESWARN,"[munuit_value]: Terms:")
       call message_value(MESWARN,"[munuit_value]:      CSM:",modcsmmin_minuit_value)
       call message_value(MESWARN,"[munuit_value]: srhosrho:",srhorho)
       call message_value(MESWARN,"[munuit_value]:    simpj:",sum(impj(:,:,:)))
       call message_value(MESWARN,"[munuit_value]:    mfact:",((real(lnsymop,FREAL)-1.0_FREAL)/(real(lnsymop,FREAL))))
       call message_value(MESWARN,"[munuit_value]:    ratio:",((sum(impj(:,:,:)))/(real(lnsymop,FREAL)*srhorho)))
       call csmmod_getimj(impj)
       call message_value(MESWARN,"[munuit_value]:   AltDef:",1 - (sum(impj(:,:,:))/real(lnsymop,FREAL))/(srhorho))
       call message_value(MESWARN,"[munuit_value]:     simj:",sum(impj(:,:,:)))
       call message_value(MESWARN,"[munuit_value]:    ratio:",((sum(impj(:,:,:)))/(real(lnsymop,FREAL)*srhorho)))
       call csmmod_getimpj(impj)
       if (abs(modcsmmin_minuit_value).lt.eps) then
          modcsmmin_minuit_value=0.0_FREAL
          call message(MESWARN,"[munuit_value]: Negative value smaller then EPS: Reset to ZERO.")
       end if
    end if
    !
  end function modcsmmin_minuit_value
  !
end module modcsmmin






