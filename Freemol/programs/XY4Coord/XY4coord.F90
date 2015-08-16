!! substitute:
!! PROGNAME with the PROGRAM NAME
!! DATE with the current date
!! NOTE with a NOTE
!!
!!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H PROGRAM XY4coord  Frimol by F.Mariotti:
!H-----------------------------------------------------------------------------
!H $Id: XY4coord.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
!H-----------------------------------------------------------------------------
!H This File is part of the Frimol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H 
!H DESCRIPTION:
!H  We implement some stuff from JCP 119,11513 plus some our tests and so on..
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
program XY4coord
  !
  use vartypes
  use messages
  use pcmdline
  use linetools
  use baseio
  use osec
  use molecule
  !
  implicit none
  !
  ! DECLARATIONS:
  !--------------------
  
  ! Procedures JUST a List dont remove comments
  !------------------
  !integer(FINT) :: XY4coord_init
  !logical ::  XY4coord_isinit
  !integer(FINT) :: XY4coord_geterror
  !integer(FINT), public :: XY4coord_
  !
  !
  ! LOCAL DECLARATIONS:
  !--------------------
  real(FREAL), parameter :: LPI = 3.1415926535897932384626433832795028841971694_FREAL
  !real(FREAL), parameter :: XY4LEPS = 1.0E-14_FREAL
  real(FREAL), parameter :: XY4LEPS = 1.0E-8_FREAL
  integer(FINT) :: i,j,k
  !
  ! return code from functions and the running one
  !-----------------------------------------------
  integer(FINT) :: irc,ircr
  character(FLCHARS) :: params
  character(FLCHARS) :: line
  real(FREAL) :: ftmp,frtmp,f1d144,fr2,fr3,fr6,f1d72,fa109
  character(FLCHARS) :: cscanopt
  real(FREAL),dimension(10) :: fscanpars
  integer(FINT),dimension(30) :: fscanvars
  integer(FINT) :: iscan
  real(FREAL), pointer :: fscancoord
  !
  !Coordinates
  !-----------
  ! Cartesian Coordinates
  ! Convention: X is atom number 1 and the 4 Y are atoms number 2-5
  !----------------------------------------------------------------
  real(FREAL), dimension(3,5), save :: mxyz
  !
  ! Initial Bonds and Angles
  real(FREAL), dimension(4), save :: xy,rxy
  real(FREAL), dimension(6), save :: yxy,ryxy
  ! Saved cartesian Coordinates: Same convention as above
  !----------------------------------------------------------
  real(FREAL), dimension(3,5), save :: smxyz
  ! Internal Coordinates
  ! Definition: dr(i)=X-Yi, da(1:6)=a12,a13,a14,a23,a24,a34
  !                                  1   2   3   4   5   6
  !--------------------------------------------------------
  real(FREAL), dimension(4), save, target :: dr,vr,svr
  real(FREAL), dimension(6), save, target :: da,va,sva
  ! Gamma Angles to store for calculation
  ! Gamma angles definition:
  ! G(1)=a12/3 G( 7)=a12/4
  ! G(2)=a13/2 G( 8)=a13/4
  ! G(3)=a14/2 G( 9)=a14/3
  ! G(4)=a23/1 G(10)=a23/4
  ! G(5)=a24/1 G(11)=a24/3
  ! G(6)=a34/1 G(12)=a34/2
  real(FREAL), dimension(12), save :: agam
  real(FREAL), dimension(12), save :: acagam
  integer(FINT) :: ig
  real(FREAL) :: fgamma
  real(FREAL), dimension(3) :: gatmp
  ! to use cosine SymCoord definition
  real(FREAL), dimension(6), save :: cda,vcda
  ! to store the actual alpha values from cosines
  real(FREAL), dimension(6), save :: acda,vacda
  ! The cosines and sines
  real(FREAL), dimension(6), save :: csda,snda
  ! Symmetry Coordinates
  ! Definition:
  ! Sdr(1)=S1, Sdr(2-4)=S2x,S2y,S2z, Sda(1-2)=S2a,S2b, Sda(3-5)=S4x,S4y,S4z, Sda(6)=Sr
  !-----------------------------------------------------------------------------------
  real(FREAL), dimension(4), save, target :: Sdr
  real(FREAL), dimension(6), save, target :: Sda
  !
  ! Rhos
  !-----
  ! Definition: rho1 = da12 + da13 + da14
  !             rho2 = da12 + da23 + da24
  !             rho3 = da13 + da23 + da34
  !             rho4 = da14 + da24 + da34
  !--------------------------------------
  real(FREAL), dimension(4), save :: rho, srho
  !
  ! Solutions for Sda(6)
  !---------------------
  complex(FREAL), dimension(4) :: ZSda
  complex(FREAL), dimension(4) :: ZOSda
  ! Dummy complex
  complex(FREAL) :: Ztmp
  ! index of max real Z solution
  integer(FINT) :: iZ
  ! Running index of max real Z solution
  integer(FINT) :: isZ
  !
  ! Files
  !------
  integer(FINT), save :: iuninput
  character(FLCHARS), save :: fileinput
  integer(FINT), save :: iunoutput
  character(FLCHARS), save :: fileoutput
  integer(FINT), save :: iunascan
  character(FLCHARS), save :: fileascan
  integer(FINT), save :: iunscan
  character(FLCHARS), save :: filescan
  integer(FINT), save :: iunxyz
  character(FLCHARS), save :: filexyz
  !
  ! Command Line
  !-------------
  integer(FINT) :: nargs
  character(FLCHARS),dimension(:), allocatable :: cargs
  !
  !
  ! error variables
  !----------------
  logical, save :: LOCXY4coord_isinit
  integer(FINT), save :: LOCXY4coord_error
  character(FLCHARS), save :: LOCXY4coord_error_message
  !
  ! Precompiute some values
  !------------------------
  f1d144 = 1.0_FREAL/144.0_FREAL
  f1d72  = 1.0_FREAL/72.0_FREAL
  fr2    = sqrt(2.0_FREAL)
  fr3    = sqrt(3.0_FREAL)
  fr6    = sqrt(6.0_FREAL)
  !
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Here we start
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !
  ! Initialization of routines
  !---------------------------
  irc = XY4coord_init()
  if(irc.lt.0) then
     call message(MESERRO,"[XY4C] Cannot initialize.")
     stop 1
  end if
  !
  ! Start reading
  !--------------
  !
  ! Molecule Section which will be our eckard frame
  !------------------------------------------------
  rewind(iuninput)
  irc = osec_set(iuninput,'molecule',params)
  if(irc.lt.0) then
     call message(MESERRO,"[XY4C] No molecule Section.")
     stop 1
  end if
  call molecule_read(iuninput,params)
  irc=molecule_getnat()
  if(irc.ne.5) then
     call message(MESERRO,"[XY4C] This is XY4 program: we need exactly 5 atoms!.")
     stop 1
  end if
  call message(MESOUT,"[XY4C] Eckard reference system readed in.")
  call molecule_print(iunoutput)
  !
  ! Read in xyz from input
  !-----------------------
  mxyz(1:3,1:5)=molecule_xyz(1:3,1:5)
  smxyz(1:3,1:5)=molecule_xyz(1:3,1:5)
  !
  ! Get bonds and angles: remember to skip atom at origin!
  !-------------------------------------------------------
  call getbonds(4,mxyz(1:3,2:5),xy(1:4))
  call getangles(4,mxyz(1:3,2:5),xy(1:4),yxy(1:6))
  write (*,'("#[x-out-comment] Initial Bonds and Angles (GRAD).")')
  write (*,'("#",1X,10(2X,A11))') 'r1','r2','r3','r4','a12','a13','a14','a23','a24','a34'
  write (*,'("#",1X,10(2X,F11.4))') xy(1:4),yxy(1:6)/LPI*180.0_FREAL
  !
  !
  ! We read the symmetric displacment coordinate system
  !----------------------------------------------------
  rewind(iuninput)
  irc = osec_set(iuninput,'x-xy4-symmcoord',params)
  if(irc.lt.0) then
     call message(MESERRO,"[XY4C] No x-xy4-symmcoord Section.")
     stop 1
  end if
  !
  !
  irc = line_getline(iuninput,line,3)
  read(line,*) Sdr(1:4),Sda(1:6)
  call message(MESOUT,"[XY4C] Symmetric coordinates data readed in.")
  write (*,'("[x-xy4-symmcoord] INPUT")')
  write (*,'("#",A11,9(1X,A11))') 'S1','S2x','S2y','S2z','S2a','S2b','S4x','S4y','S4z','Sr'
  write (*,'(10(1X,F11.4))') Sdr(1:4),Sda(1:6)
  !
  ! Run this input data
  !--------------------
  call message(MESOUT,"[XY4C] ")
  call message(MESOUT,"[XY4C] Start run on input values.")
  call message(MESOUT,"[XY4C] ")
  !
  ! Get internal coordinates: dr,da,vr,va
  !--------------------------------------
  call get_ra
  !
  ! Do Checkings
  !-------------
  irc=do_checks()
  if (irc.ne.0) then
     call message_value(MESWARN,"[XY4C] Check not Passed at Sr input value:",Sda(6))
     call message_value(MESOUT,"[XY4C] Check not Passed at Sr input value:",Sda(6))
  else
     !
     ! Get Cartesians
     !---------------
     irc=get_cart()
     if (irc.ne.0) then
        call message(MESOUT,"[XY4C] ")
        call message(MESERRO,"[XY4C] Error in Cartesian routine.")
        call message(MESOUT,"[XY4C] ")
     else
        call message(MESOUT,"[XY4C] ")
        call message_value(MESOUT,"[XY4C] Check OK at Sr input value:",Sda(6))
        call message(MESOUT,"[XY4C] ")
     end if
  end if
  !
  ! Evaluate possible Sr solutions
  !-------------------------------
  call eval_sr()
  call get_reZ()
  call message(MESOUT,"[XY4C] ")
  call message(MESOUT,"[XY4C] Generate Redundancies")
  call message(MESOUT,"[XY4C] ")
  write (*,'("#[x-out-comment] Generate Sr redundancies.")')
  if (iZ.gt.0) then
     write (*,'(A13,1X,4(2X,F12.4,1X,"i(",F12.4,")"))') '#Z Re Solut.:',ZSda(1:iZ)
  else
     write (*,'("#Z Re Solut.: NONE")')
     call message(MESERRO," No Available redundant solution. Keep going but there is an error!")
     call message(MESOUT," No Available redundant solution. Keep going but there is an error!")
  end if
  write (*,'(A13,1X,4(2X,F12.4,1X,"i(",F12.4,")"))') '#Z solutions:',ZOSda
  !
  ! Run on available solutions
  !---------------------------
  call message(MESOUT,"[XY4C] ")
  call message_value(MESOUT,"[XY4C] Start run on Available Sr solutions. N:",iZ)
  call message(MESOUT,"[XY4C] ")
  ircr=-1
  do isZ=1,iZ
     call message_value(MESOUT,"[XY4C] Checking Solution:",isZ)
     Sda(6)=ZSda(isZ)
     write (*,'("[x-xy4-symmcoord] Sol:",1X,I4,1X,"Sr Value:",1X,F11.4)') isZ,Sda(6)
     write (*,'("#",A11,9(1X,A11))') 'S1','S2x','S2y','S2z','S2a','S2b','S4x','S4y','S4z','Sr'
     write (*,'(10(1X,F11.4))') Sdr(1:4),Sda(1:6)
     !
     ! Get internal coordinates: dr,da,vr,va
     !--------------------------------------
     call get_ra
     !
     ! Do Checkings
     !-------------
     irc=do_checks()
     if (irc.ne.0) then
        call message_value(MESWARN,"[XY4C] Check not Passed at Sr solution:",isZ)
        call message_value(MESOUT,"[XY4C] Check not Passed at Sr solution:",isZ)
     else
        ircr=0
        !
        ! Get Cartesians
        !---------------
        irc=get_cart()
        if (irc.ne.0) then
           call message(MESOUT,"[XY4C] ")
           call message(MESERRO,"[XY4C] Error in Cartesian routine.")
           call message(MESOUT,"[XY4C] ")
        else
           call message_value(MESOUT,"[XY4C] Check OK at Sr solution:",isZ)
        end if
     end if
     !
  end do
  if (ircr.ne.0) then
     call message(MESERRO,"[XY4C]: No solution has been found at starting point. We DO stop.")
     call message(MESOUT,"[XY4C]: No solution has been found at starting point. We DO stop.")
     stop 3
  end if
  !
  ! Check if we have to perform a scan: i.e. x-XY4-scan section
  !------------------------------------------------------------
  rewind(iuninput)
  irc = osec_set(iuninput,'x-xy4-scan',params)
  if(irc.lt.0) then
     call message(MESOUT,"[XY4C] No scan section.")
  else
     call do_scan
  end if
  !
  !
  ! Check if we have to perform a revscan: i.e. x-XY4-scan section
  !------------------------------------------------------------
  rewind(iuninput)
  irc = osec_set(iuninput,'x-xy4-revscan',params)
  if(irc.lt.0) then
     call message(MESOUT,"[XY4C] No revscan section.")
  else
     call do_revscan
  end if
  !
  call message(MESOUT,"[XY4C] Normal Exit but check for messages.")
  !
  !
contains
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H subroutine get_ra()
!H-----------------------------------------------------------------------------
!H
!
  subroutine get_ra()
    implicit none

!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine 
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    !r1,r2,r3,r4
    !-----------
    dr(1)=0.125_FREAL*(Sdr(1)+Sdr(2)-Sdr(3)+Sdr(4))
    dr(2)=0.125_FREAL*(Sdr(1)-Sdr(2)+Sdr(3)+Sdr(4))
    dr(3)=0.125_FREAL*(Sdr(1)+Sdr(2)+Sdr(3)-Sdr(4))
    dr(4)=0.125_FREAL*(Sdr(1)-Sdr(2)-Sdr(3)-Sdr(4))
    !
    ! a12,a13,a14,a23,a24,a34
    !------------------------
    !
    da(1)= Sda(6)/fr6 + Sda(1)/fr3 - Sda(5)/fr2
    da(2)= Sda(6)/fr6 - Sda(1)/fr3 + Sda(2)/fr2 - Sda(3)/fr2
    da(3)= Sda(6)/fr6 - Sda(1)/fr3 - Sda(2)/fr2 - Sda(4)/fr2
    da(4)= Sda(6)/fr6 - Sda(1)/fr3 - Sda(2)/fr2 + Sda(4)/fr2
    da(5)= Sda(6)/fr6 - Sda(1)/fr3 + Sda(2)/fr2 + Sda(3)/fr2
    da(6)= Sda(6)/fr6 + Sda(1)/fr3 + Sda(5)/fr2
    !write (*,'(1X,A6,1X,6(2X,F14.6))') 'Angle:',da(1:6)
    !da(1:6)=acos(da(1:6))
    !write (*,'(1X,A6,1X,6(2X,F14.6))') 'Angle:',da(1:6)
    write (*,'("#[x-out-comment] Displacments from SymCoord: Bonds and Angles (GRAD).")')
    write (*,'("#",1X,10(2X,A11))') 'r1','r2','r3','r4','a12','a13','a14','a23','a24','a34'
    write (*,'("#",1X,10(2X,F11.4))' ) dr(1:4),da(1:6)/LPI*180.0_FREAL
    !
    ! Get full coordinates i.e. add displacments
    !--------------------------------------------
    vr(1:4)=xy(1:4)+dr(1:4)
    va(1:6)=yxy(1:6)+da(1:6)
    write (*,'("#[x-out-comment] Displaced Coordinates from SymCoord: Bonds and Angles (GRAD).")')
    write (*,'("#",1X,10(2X,A11))') 'r1','r2','r3','r4','a12','a13','a14','a23','a24','a34'
    write (*,'("#",1X,10(2X,F11.4))' ) vr(1:4),va(1:6)/LPI*180.0_FREAL
    !
    ! Evaluate Gammas
    !----------------
    agam(1)=(cos(va(1))-(cos(va(2))*cos(va(4))))/(sin(va(2))*sin(va(4)))
    agam(2)=(cos(va(2))-(cos(va(1))*cos(va(4))))/(sin(va(1))*sin(va(4)))
    agam(3)=(cos(va(3))-(cos(va(1))*cos(va(5))))/(sin(va(1))*sin(va(5)))
    agam(4)=(cos(va(4))-(cos(va(1))*cos(va(2))))/(sin(va(1))*sin(va(2)))
    agam(5)=(cos(va(5))-(cos(va(1))*cos(va(3))))/(sin(va(1))*sin(va(3)))
    agam(6)=(cos(va(6))-(cos(va(2))*cos(va(3))))/(sin(va(2))*sin(va(3)))
    !
    agam( 7)=(cos(va(1))-(cos(va(3))*cos(va(5))))/(sin(va(3))*sin(va(5)))
    agam( 8)=(cos(va(2))-(cos(va(3))*cos(va(6))))/(sin(va(3))*sin(va(6)))
    agam( 9)=(cos(va(3))-(cos(va(2))*cos(va(6))))/(sin(va(2))*sin(va(6)))
    agam(10)=(cos(va(4))-(cos(va(5))*cos(va(6))))/(sin(va(5))*sin(va(6)))
    agam(11)=(cos(va(5))-(cos(va(4))*cos(va(6))))/(sin(va(4))*sin(va(6)))
    agam(12)=(cos(va(6))-(cos(va(4))*cos(va(5))))/(sin(va(4))*sin(va(5)))
    !
    ! Evaluate ACOS of gammas
    !------------------------
    acagam(:)=acos(agam(:))
    !
    write (*,'("#[x-out-comment] Gamma Angles (GRAD).")')
    write (*,'("#",1X,6(2X,A11))') 'g12/3','g13/2','g14/2','g23/1','g24/1','g34/1'
    write (*,'("#",1X,6(2X,A11))') 'g12/4','g13/4','g14/3','g23/4','g24/3','g34/2'
    write (*,'("#",1X,6(2X,F11.4))' ) acagam(1:6)/LPI*180.0_FREAL
    write (*,'("#",1X,6(2X,F11.4))' ) acagam(7:12)/LPI*180.0_FREAL
    !
  !
  ! Get Rhos
  !---------
  rho(1)=da(1)+da(2)+da(3)
  rho(2)=da(1)+da(4)+da(5)
  rho(3)=da(2)+da(4)+da(6)
  rho(4)=da(3)+da(5)+da(6)
  write (*,'("#[x-out-comment] Rhos (RAD).")')
  write (*,'("#",1X,4(2X,F11.4))' ) rho(1:4)
  write (*,'("#[x-out-comment]  Sum of Rhos (RAD).",1X,2X,F11.4)') sum(rho(1:4))
  !
  !
  end subroutine get_ra
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function do_checks
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function do_checks()
    implicit none
    logical :: chkdebug
    !chkdebug=.TRUE.
    chkdebug=.FALSE.
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine 
!H
!H-----------------------------------------------------------------------------
!H
!
    do_checks=0
    !
    ! We require internal angles to be defined within
    if (chkdebug) call message(MESLOG,"[XY4C][CHECK] angles 0=<A=<180")
    do i=1,6
       if (chkdebug) call message_value(MESLOG,"[XY4C][CHECK] angles",va(i)/LPI*180.0_FREAL)
       if (va(i).lt.0.0_FREAL) then
          call message(MESERRO,"[XY4C] Inconsistent negative angle")
          call message_value(MESERRO,"[XY4C] term",i)
          do_checks=-1
       end if
       if (va(i).gt.LPI) then
          call message(MESWARN,"[XY4C] Too big angle: please checkit!")
          call message_value(MESERRO,"[XY4C] term",i)
          do_checks=-1
       end if
    end do
    !
    do i=1,4
       if(vr(i).lt.0.0001_FREAL) then
          call message(MESERRO,"[XY4C]: Very Small bond distance.")
          call message_value(MESERRO,"[XY4C] term",i)
          do_checks=1
       end if
    end do
    !
    ! Check XY3 type of sum
    !----------------------
    if (chkdebug) call message(MESLOG,"[XY4C][CHECK] angles Ammonia type < 360")
    if (chkdebug) call message_value(MESLOG,"[XY4C][CHECK] Y234",(va(4)+va(5)+va(6))/LPI*180.0_FREAL)
    if ((va(4)+va(5)+va(6)).gt.(2.0_FREAL*LPI)) then
       call message_value(MESERRO,"[XY4C]: Not realistic internal angles for Y234:",(va(4)+va(5)+va(6))/LPI*180.0_FREAL)
       do_checks=-1
    end if
    if (chkdebug) call message_value(MESLOG,"[XY4C][CHECK] Y134",(va(2)+va(3)+va(6))/LPI*180.0_FREAL)
    if ((va(2)+va(3)+va(6)).gt.(2.0_FREAL*LPI)) then
       call message_value(MESERRO,"[XY4C]: Not realistic internal angles for Y134:",(va(2)+va(3)+va(6))/LPI*180.0_FREAL)
       do_checks=-1
    end if
    if (chkdebug) call message_value(MESLOG,"[XY4C][CHECK] Y124",(va(1)+va(3)+va(5))/LPI*180.0_FREAL)
    if ((va(1)+va(3)+va(5)).gt.(2.0_FREAL*LPI)) then
       call message_value(MESERRO,"[XY4C]: Not realistic internal angles for Y124:",(va(1)+va(3)+va(5))/LPI*180.0_FREAL)
       do_checks=-1
    end if
    if (chkdebug) call message_value(MESLOG,"[XY4C][CHECK] Y123",(va(1)+va(2)+va(5))/LPI*180.0_FREAL)
    if ((va(1)+va(2)+va(4)).gt.(2.0_FREAL*LPI)) then
       call message_value(MESERRO,"[XY4C]: Not realistic internal angles for Y123:",(va(1)+va(2)+va(4))/LPI*180.0_FREAL)
       do_checks=-1
    end if
    !
    ! Check Gamma not to got infinity
    !--------------------------------
    if (chkdebug) call message(MESLOG,"[XY4C][CHECK] Gamma Defined")
    if ( (sin(va(3))*sin(va(5)) ).lt.XY4LEPS ) then
       call message(MESERRO,"[XY4C]: Error for Gamma: set 3 5")
       do_checks=-1
    end if
    if ( (sin(va(2))*sin(va(4)) ).lt.XY4LEPS ) then
       call message(MESERRO,"[XY4C]: Error for Gamma: set 2 4")
       do_checks=-1
    end if
    if ( (sin(va(3))*sin(va(6)) ).lt.XY4LEPS ) then
       call message(MESERRO,"[XY4C]: Error for Gamma: set 3 6")
       do_checks=-1
    end if
    if ( (sin(va(1))*sin(va(4)) ).lt.XY4LEPS ) then
       call message(MESERRO,"[XY4C]: Error for Gamma: set 1 4")
       do_checks=-1
    end if
    if ( (sin(va(2))*sin(va(6)) ).lt.XY4LEPS ) then
       call message(MESERRO,"[XY4C]: Error for Gamma: set 2 6")
       do_checks=-1
    end if
    if ( (sin(va(1))*sin(va(5)) ).lt.XY4LEPS ) then
       call message(MESERRO,"[XY4C]: Error for Gamma: set 1 5")
       do_checks=-1
    end if
    if ( (sin(va(5))*sin(va(6)) ).lt.XY4LEPS ) then
       call message(MESERRO,"[XY4C]: Error for Gamma: set 5 6")
       do_checks=-1
    end if
    if ( (sin(va(1))*sin(va(2)) ).lt.XY4LEPS ) then
       call message(MESERRO,"[XY4C]: Error for Gamma: set 1 2")
       do_checks=-1
    end if
    if ( (sin(va(4))*sin(va(6)) ).lt.XY4LEPS ) then
       call message(MESERRO,"[XY4C]: Error for Gamma: set 4 6")
       do_checks=-1
    end if
    if ( (sin(va(1))*sin(va(3)) ).lt.XY4LEPS ) then
       call message(MESERRO,"[XY4C]: Error for Gamma: set 1 3")
       do_checks=-1
    end if
    if ( (sin(va(4))*sin(va(5)) ).lt.XY4LEPS ) then
       call message(MESERRO,"[XY4C]: Error for Gamma: set 4 5")
       do_checks=-1
    end if
    if ( (sin(va(2))*sin(va(3)) ).lt.XY4LEPS ) then
       call message(MESERRO,"[XY4C]: Error for Gamma: set 2 3")
       do_checks=-1
    end if
    !
    ! Evaluate Gammas
    !----------------
    if (chkdebug) call message(MESLOG,"[XY4C][CHECK] Evaluate Gammas.")
    agam(1)=(cos(va(1))-(cos(va(2))*cos(va(4))))/(sin(va(2))*sin(va(4)))
    agam(2)=(cos(va(2))-(cos(va(1))*cos(va(4))))/(sin(va(1))*sin(va(4)))
    agam(3)=(cos(va(3))-(cos(va(1))*cos(va(5))))/(sin(va(1))*sin(va(5)))
    agam(4)=(cos(va(4))-(cos(va(1))*cos(va(2))))/(sin(va(1))*sin(va(2)))
    agam(5)=(cos(va(5))-(cos(va(1))*cos(va(3))))/(sin(va(1))*sin(va(3)))
    agam(6)=(cos(va(6))-(cos(va(2))*cos(va(3))))/(sin(va(2))*sin(va(3)))
    !
    agam( 7)=(cos(va(1))-(cos(va(3))*cos(va(5))))/(sin(va(3))*sin(va(5)))
    agam( 8)=(cos(va(2))-(cos(va(3))*cos(va(6))))/(sin(va(3))*sin(va(6)))
    agam( 9)=(cos(va(3))-(cos(va(2))*cos(va(6))))/(sin(va(2))*sin(va(6)))
    agam(10)=(cos(va(4))-(cos(va(5))*cos(va(6))))/(sin(va(5))*sin(va(6)))
    agam(11)=(cos(va(5))-(cos(va(4))*cos(va(6))))/(sin(va(4))*sin(va(6)))
    agam(12)=(cos(va(6))-(cos(va(4))*cos(va(5))))/(sin(va(4))*sin(va(5)))
    !
    ! Check Gamma Boudaries
    !----------------------
    if (chkdebug) call message(MESLOG,"[XY4C][CHECK] Gamma within bondaries.")
    do ig=1,12
       if (abs(agam(ig)).gt.1.0_FREAL) then
          call message_value(MESERRO,"[XY4C]: Gamma bigger then 1.0 for index:",ig)
          do_checks=-1
       end if
    end do
    !
    ! Evaluate ACOS of gammas
    !------------------------
    acagam(:)=acos(agam(:))
    !
    ! Check Gammas to add to 360
    !---------------------------
    if (chkdebug) call message(MESLOG,"[XY4C][CHECK] Gamma full rotation: GammaSum=180.")
    !
    ! H1 test
    if (chkdebug) call message_value(MESLOG,"[XY4C][CHECK] Gamma Sum for atom Y1:",acagam(4)+acagam(5)+acagam(6))
    if ((acagam(4)+acagam(5)+acagam(6)).gt.(2.0_FREAL*LPI)) then
          call message(MESERRO,"[XY4C]: Gamma Sum bigger then 180 for Y atom 1:")
          call message_value(MESLOG,"[XY4C][CHECK] Gamma Sum for atom Y1:",acagam(4)+acagam(5)+acagam(6))
          do_checks=-1
    end if
    !
    ! H2 test
    if (chkdebug) call message_value(MESLOG,"[XY4C][CHECK] Gamma Sum for atom Y2:",acagam(2)+acagam(3)+acagam(12))
    if ((acagam(2)+acagam(3)+acagam(12)).gt.(2.0_FREAL*LPI)) then
          call message(MESERRO,"[XY4C]: Gamma Sum bigger then 180 for Y atom 2:")
          call message_value(MESLOG,"[XY4C][CHECK] Gamma Sum for atom Y2:",acagam(2)+acagam(3)+acagam(12))
          do_checks=-1
    end if
    !
    ! H3 test
    if (chkdebug) call message_value(MESLOG,"[XY4C][CHECK] Gamma Sum for atom Y3:",acagam(1)+acagam(9)+acagam(11))
    if ((acagam(1)+acagam(9)+acagam(11)).gt.(2.0_FREAL*LPI)) then
          call message(MESERRO,"[XY4C]: Gamma Sum bigger then 180 for Y atom 3:")
          call message_value(MESLOG,"[XY4C][CHECK] Gamma Sum for atom Y3:",acagam(1)+acagam(9)+acagam(11))
          do_checks=-1
    end if
    !
    ! H4 test
    if (chkdebug) call message_value(MESLOG,"[XY4C][CHECK] Gamma Sum for atom Y4:",acagam(7)+acagam(8)+acagam(10))
    if ((acagam(2)+acagam(3)+acagam(12)).gt.(2.0_FREAL*LPI)) then
          call message(MESERRO,"[XY4C]: Gamma Sum bigger then 180 for Y atom 4:")
          call message_value(MESLOG,"[XY4C][CHECK] Gamma Sum for atom Y4:",acagam(7)+acagam(8)+acagam(10))
          do_checks=-1
    end if
    !
  end function do_checks
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function get_cart()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function get_cart()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine gets cartesians coordinates based on the scheme of:
!H JCP mirjana (2003) but we should check bounds for H4: 180 or 360?
!H
!H-----------------------------------------------------------------------------
!H
!
    get_cart=0
    write (*,'("#[x-out-comment] Call to Internal to cartesian coordinates. New Coordinates: X on origin.")')
    !
    ! H1
    !---
    mxyz(1,2)=vr(1)
    mxyz(2,2)=0.0_FREAL
    mxyz(3,2)=0.0_FREAL
    write(line,'(1X,A2,1X,3(2X,F12.4))') 'Y1',mxyz(1:3,2)
    call message(MESOUT,line)
    !
    ! H2
    !---
    mxyz(1,3)=vr(2)*cos(va(1))
    mxyz(2,3)=vr(2)*sin(va(1))
    mxyz(3,3)=0.0_FREAL
    write(line,'(1X,A2,1X,3(2X,F12.4))') 'Y2',mxyz(1:3,3)
    call message(MESOUT,line)
    !
    ! H3
    !---
    agam(4)=(cos(va(4))-(cos(va(1))*cos(va(2))))/(sin(va(1))*sin(va(2)))
    mxyz(1,4)=vr(3)*cos(va(2))
    mxyz(2,4)=vr(3)*sin(va(2))*agam(4)
    mxyz(3,4)=vr(3)*sin(va(2))*sin(acos(agam(4)))
    write(line,'(1X,A2,1X,3(2X,F12.4))') 'Y3',mxyz(1:3,4)
    call message(MESOUT,line)
    !
    ! H4
    !---
    agam(5)=(cos(va(5))-(cos(va(1))*cos(va(3))))/(sin(va(1))*sin(va(3)))
    mxyz(1,5)= vr(4)*cos(va(3))
    mxyz(2,5)= vr(4)*sin(va(3))*agam(5)
    mxyz(3,5)= vr(4)*sin(va(3))*sin(acos(agam(5)))
    agam(2)=(cos(va(2))-(cos(va(1))*cos(va(4))))/(sin(va(1))*sin(va(4)))
    agam(12)=(cos(va(6))-(cos(va(4))*cos(va(5))))/(sin(va(4))*sin(va(5)))
    !if ((acos(agam(2))+acos(agam(3))).gt.LPI) then
    if ((acos(agam(2))+acos(agam(3))+acos(agam(12)))-(2.0_FREAL*LPI).lt.XY4LEPS) then
       mxyz(3,5)=-mxyz(3,5)
    else
       ! All this shit down there is to get sorted agam....
       if (acos(agam(2)).gt.acos(agam(3))) then
          if (acos(agam(2)).gt.acos(agam(12))) then
             gatmp(1)=acos(agam(2))
             if (acos(agam(3)).gt.acos(agam(12))) then
                gatmp(2)=acos(agam(3))
                gatmp(3)=acos(agam(12))
             else
                gatmp(2)=acos(agam(12))
                gatmp(3)=acos(agam(3))
             end if
          else
             gatmp(1)=acos(agam(12))
             gatmp(2)=acos(agam(2))
             gatmp(3)=acos(agam(3))
          end if
       else
          if (acos(agam(3)).gt.acos(agam(12))) then
             gatmp(1)=acos(agam(3))
             if (acos(agam(2)).gt.acos(agam(12))) then
                gatmp(2)=acos(agam(2))
                gatmp(3)=acos(agam(12))
             else
                gatmp(2)=acos(agam(12))
                gatmp(3)=acos(agam(2))
             end if
          else
             gatmp(1)=acos(agam(12))
             gatmp(2)=acos(agam(3))
             gatmp(3)=acos(agam(2))
          end if
       endif
       if ((gatmp(1)-gatmp(2)-gatmp(3)).gt.XY4LEPS) then
          call message(MESERRO,"Not consistent gammas for H4 evaluation")
       end if
    endif
    !
    !call message_value(MESOUT,"H4 Gamma  2",acos(agam( 2))/LPI*180.0_FREAL)
    !call message_value(MESOUT,"H4 Gamma  3",acos(agam( 3))/LPI*180.0_FREAL)
    !call message_value(MESOUT,"H4 Gamma 12",acos(agam(12))/LPI*180.0_FREAL)
    !call message_value(MESOUT,"H4 FGamma 1",gatmp(1)/LPI*180.0_FREAL)
    !call message_value(MESOUT,"H4 FGamma 2",gatmp(2)/LPI*180.0_FREAL)
    !call message_value(MESOUT,"H4 FGamma 3",gatmp(3)/LPI*180.0_FREAL)
    !
    write(line,'(1X,A2,1X,3(2X,F12.4))') 'Y4',mxyz(1:3,5)
    call message(MESOUT,line)
    call message(MESWARN,"Assumed signed 'sin'. It can be inconsistent with input sym data.")
    !
    !
    ! Get New bonds and Angles
    !-------------------------
    rxy(1)=sqrt(sum(mxyz(1:3,2)*mxyz(1:3,2)))
    rxy(2)=sqrt(sum(mxyz(1:3,3)*mxyz(1:3,3)))
    rxy(3)=sqrt(sum(mxyz(1:3,4)*mxyz(1:3,4)))
    rxy(4)=sqrt(sum(mxyz(1:3,5)*mxyz(1:3,5)))
    ryxy(1)=acos(sum(mxyz(1:3,2)*mxyz(1:3,3))/(xy(1)*xy(2)))
    ryxy(2)=acos(sum(mxyz(1:3,2)*mxyz(1:3,4))/(xy(1)*xy(3)))
    ryxy(3)=acos(sum(mxyz(1:3,2)*mxyz(1:3,5))/(xy(1)*xy(4)))
    ryxy(4)=acos(sum(mxyz(1:3,3)*mxyz(1:3,4))/(xy(2)*xy(3)))
    ryxy(5)=acos(sum(mxyz(1:3,3)*mxyz(1:3,5))/(xy(2)*xy(4)))
    ryxy(6)=acos(sum(mxyz(1:3,4)*mxyz(1:3,5))/(xy(3)*xy(4)))
    !
    if (sum(abs(ryxy(1:6)-va(1:6))).gt.XY4LEPS) then
       call message(MESERRO,"Evaluation of Angles from cartesian doesn't match Alphas.")
       write(line,'(1X,6(2X,F11.4))') ryxy(1:6)
       call message(MESERRO,line)
       write(line,'(1X,6(2X,F11.4))') va(1:6)
       call message(MESERRO,line)
       get_cart=-1
       return
    end if
    !
    write(*,'("[x-xyz] 5")')
    write(*,'(1X,A1,1X,3(2X,F12.4))') 'X',mxyz(1:3,1)
    write(*,'(1X,A1,1X,3(2X,F12.4))') 'Y',mxyz(1:3,2)
    write(*,'(1X,A1,1X,3(2X,F12.4))') 'Y',mxyz(1:3,3)
    write(*,'(1X,A1,1X,3(2X,F12.4))') 'Y',mxyz(1:3,4)
    write(*,'(1X,A1,1X,3(2X,F12.4))') 'Y',mxyz(1:3,5)
    !
    !
  end function get_cart
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H subroutine get_symc
!H-----------------------------------------------------------------------------
!H
!
  subroutine get_symc()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine 
!H
!H-----------------------------------------------------------------------------
!H
!
    !Sdr
    Sdr(1)=( dr(1)+dr(2)+dr(3)+dr(4))/sqrt(4.0_FREAL)
    Sdr(2)=( dr(1)-dr(2)+dr(3)+dr(4))/sqrt(4.0_FREAL)
    Sdr(3)=(-dr(1)+dr(2)+dr(3)+dr(4))/sqrt(4.0_FREAL)
    Sdr(4)=( dr(1)+dr(2)-dr(3)-dr(4))/sqrt(4.0_FREAL)
    !
    !Sda
    Sda(1)=(2.0_FREAL*da(1)-da(2)-da(3)-da(4)-da(5)+2.0_FREAL*da(6))/sqrt(12.0_FREAL)
    Sda(2)=(da(2)-da(3)-da(4)+da(5))/2.0_FREAL
    Sda(3)=(da(5)-da(3))/sqrt(2.0_FREAL)
    Sda(4)=(da(4)-da(2))/sqrt(2.0_FREAL)
    Sda(5)=(da(6)-da(1))/sqrt(2.0_FREAL)
    Sda(6)=(da(1)+da(2)+da(3)+da(4)+da(5)+da(6))/sqrt(6.0_FREAL)
    ! We now get redundancy
    call eval_sr()
    call get_reZ()
    write (*,'("#[x-out-comment] Generate Sr redundancies.")')
    if (iZ.gt.0) then
       write (*,'(A13,1X,4(2X,F12.4,1X,"i(",F12.4,")"))') '#Z Re Solut.:',ZSda(1:iZ)
    else
       write (*,'("#Z Re Solut.: NONE")')
       call message(MESERRO," No Available redundant solution. Keep going but there is an error!")
       call message(MESOUT," No Available redundant solution. Keep going but there is an error!")
    end if
    write (*,'(1X,A12,1X,4(2X,F14.6,1X,"i(",F14.6,")"))') 'Z solutions:',ZOSda
!
!
    call message(MESERRO,"TODO: perform an internal check! redundant should be consistent!")
  end subroutine get_symc
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H subroutine do_scan
!H-----------------------------------------------------------------------------
!H
!
  subroutine do_scan()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine 
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    filescan='scan.dat'
    irc = baseio_open(iunscan,filescan,stat='UNKNOWN')
    if(irc.lt.0) then
       LOCXY4coord_error = 3
       call message(MESERRO,"[XY4C] Cannot Open SCAN File")
       stop 1
    end if
    fileascan='ascan.dat'
    irc = baseio_open(iunascan,fileascan,stat='UNKNOWN')
    if(irc.lt.0) then
       LOCXY4coord_error = 3
       call message(MESERRO,"[XY4C] Cannot Open ASCAN File")
       stop 1
    end if
    filexyz='xyz.dat'
    irc = baseio_open(iunxyz,filexyz,stat='UNKNOWN')
    if(irc.lt.0) then
       LOCXY4coord_error = 3
       call message(MESERRO,"[XY4C] Cannot Open XYZ File")
       stop 1
    end if
    !
    irc = line_getline(iuninput,line,3)
    if (line(1:1).eq.'E') then
       call message(MESOUT,"[XY4C] Stop request from scan.")
       stop 0
    end if
    if (line(1:1).eq.'R') then
       call message(MESOUT,"[XY4C] Running Reverse multi scan.")
    !==============================================================================================
       read(line,*) cscanopt,iscan
       call message_value(MESOUT,"Number of Points in Reverse MultiScan:",iscan)
       do k=1,iscan
          call message_value(MESOUT,"MultiScan Step:",k)
          irc = line_getline(iuninput,line,3)
          read(line,*) dr(1:4),da(1:6)
          da(1:6)=da(1:6)/180.0_FREAL*LPI
          call get_symc
          svr(1:4)=dr(1:4)
          sva(1:6)=da(1:6)
          ! Get full coordinates i.e. add displacments
          !--------------------------------------------
          vr(1:4)=xy(1:4)+dr(1:4)
          call message(MESOUT,"Bond lengths from internal displacments:")
          write (*,'(1X,A1,1X,4(2X,F14.6))') ' ',vr(1:4)
          !
          call message(MESWARN,"Assuming initial angles 109. NOT evaluated now...:")
          va(1:6)=((109.4712_FREAL/180.0_FREAL)*LPI)
          va(1:6)=va(1:6)+da(1:6)
          call message(MESOUT,"Angles (GRAD) from internal displacments:")
          write (*,'(1X,A1,1X,6(2X,F14.6))') ' ',va(1:6)/LPI*180.0_FREAL
          call get_ra
          do i=1,4
             if ((svr(i)-dr(i)).gt.XY4LEPS) then
                call message_value(MESERRO,"Double inversion doesnt Match. R Coord:",i)
             end if
             if ((sva(i)-da(i)).gt.XY4LEPS) then
                call message_value(MESERRO,"Double inversion doesnt Match. A Coord:",i)
                call message_value(MESERRO,"Value Origi:",sva(i))
                call message_value(MESERRO,"Value After:",da(i))
                call message_value(MESERRO,"Delta Value:",sva(i)-da(i))
             end if
          end do
          do i=5,6
             if ((sva(i)-da(i)).gt.XY4LEPS) then
                call message_value(MESERRO,"Double inversion doesnt Match. A Coord:",i)
                call message_value(MESERRO,"Value Origi:",sva(i))
                call message_value(MESERRO,"Value After:",da(i))
                call message_value(MESERRO,"Delta Value:",sva(i)-da(i))
             end if
          end do
          irc=do_checks()
          !
          ! Restore Original Values
          call message(MESOUT,"Resoring input coordinates for reverse scan")
          dr(1:4)=svr(1:4)
          da(1:6)=sva(1:6)
          !
          ! Get full coordinates i.e. add displacments
          !--------------------------------------------
          vr(1:4)=xy(1:4)+dr(1:4)
          call message(MESOUT,"Bond lengths from internal displacments:")
          write (*,'(1X,A1,1X,4(2X,F14.6))') ' ',vr(1:4)
          !
          call message(MESWARN,"Assuming initial angles 109. NOT evaluated now...:")
          va(1:6)=((109.4712_FREAL/180.0_FREAL)*LPI)
          va(1:6)=va(1:6)+da(1:6)
          call message(MESOUT,"Angles (GRAD) from internal displacments:")
          write (*,'(1X,A1,1X,6(2X,F14.6))') ' ',va(1:6)/LPI*180.0_FREAL
          !
          ! Check it back again
          call message(MESOUT,"Run Check values again.")
          irc=do_checks()
          if (irc.ne.0) then
             call message(MESERRO,"[XY4C] Check not Passed.")
             !stop 3
             write (iunscan,'(1X,I4,1X,20(2X,F14.6))') k,Sdr(1:4),Sda(1:6),-vr(1:4),va(1:6)
             write (iunascan,'(1X,I4,1X,20(2X,F14.6))') k,Sdr(1:4),Sda(1:6),-vr(1:4),va(1:6)/LPI*180.0_FREAL
          else
             write (iunscan,'(1X,I4,1X,20(2X,F14.6))') k,Sdr(1:4),Sda(1:6),vr(1:4),va(1:6)
             write (iunascan,'(1X,I4,1X,20(2X,F14.6))') k,Sdr(1:4),Sda(1:6),vr(1:4),va(1:6)/LPI*180.0_FREAL
             irc=get_cart()
             if (irc.ne.0) then
                call message(MESOUT,"[XY4C] ")
                call message(MESERRO,"[XY4C] Error in Cartesian routine.")
                call message(MESOUT,"[XY4C] ")
             else
                write (iunxyz,'(1X,I4)') 5
                write (iunxyz,'(1X,F14.6,1X)') real(k)
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'C',mxyz(1:3,1)*0.52917706E0_FREAL
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,2)*0.52917706E0_FREAL
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,3)*0.52917706E0_FREAL
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,4)*0.52917706E0_FREAL
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,5)*0.52917706E0_FREAL
             end if
          end if


          
       end do

    end if
    !==============================================================================================
    if (line(1:1).eq.'M') then
    !==============================================================================================
       call message(MESOUT,"[XY4C] ")
       call message(MESOUT,"[XY4C] Readed in Option M (Multi Scan)")
       call message(MESOUT,"[XY4C] ")
       read(line,*) cscanopt,iscan
       call message_value(MESOUT,"Number of Points in MultiScan:",iscan)
       do k=1,iscan
          call message_value(MESOUT,"MultiScan Step:",k)
          irc = line_getline(iuninput,line,3)
          read(line,*) Sdr(1:4),Sda(1:6)
          ! Eval redundancy
          call eval_sr()
          call get_reZ()
          write (*,'("#[x-out-comment] Generate Sr redundancies.")')
          if (iZ.gt.0) then
             write (*,'(A13,1X,4(2X,F12.4,1X,"i(",F12.4,")"))') '#Z Re Solut.:',ZSda(1:iZ)
          else
             write (*,'("#Z Re Solut.: NONE")')
             call message(MESERRO," No Available redundant solution. Keep going but there is an error!")
             call message(MESOUT," No Available redundant solution. Keep going but there is an error!")
          end if
          write (*,'(A13,1X,4(2X,F12.4,1X,"i(",F12.4,")"))') '#Z solutions:',ZOSda
          !Sda(6)=(fr3+3.0_FREAL*Sda(3)+2.0_FREAL*sqrt(2.0_FREAL+3.0_FREAL*sqrt(2.0_FREAL)*Sda(3)))/fr3
          !Sda(6)=(fr3+3.0_FREAL*Sda(3)-2.0_FREAL*sqrt(2.0_FREAL+3.0_FREAL*sqrt(2.0_FREAL)*Sda(3)))/fr3
          call message(MESOUT,"[XY4C] ")
          call message_value(MESOUT,"[XY4C] Start run on Available Sr solutions. N:",iZ)
          call message(MESOUT,"[XY4C] ")
          ircr=-1
          do isZ=1,iZ
             call message(MESOUT,"[XY4C] ")
             call message_value(MESOUT,"[XY4C] Checking Solution:",isZ)
             call message(MESOUT,"[XY4C] ")
             Sda(6)=ZSda(isZ)
             write (*,'("[x--step-xy4-symmcoord] Sol:",1X,I4,1X,"Sr Value:",1X,F11.4)') isZ,Sda(6)
             write (*,'("#",1X,10(1X,A11))') 'S1','S2x','S2y','S2z','S2a','S2b','S4x','S4y','S4z','Sr'
             write (*,'("#",1X,10(1X,F11.4))') Sdr(1:4),Sda(1:6)
             call get_ra
             irc=do_checks()
             if (irc.ne.0) then
                call message(MESERRO,"[XY4C] ")
                call message_value(MESERRO,"[XY4C] Check not Passed at Solution:",isZ)
                call message(MESERRO,"[XY4C] ")
                write (iunscan,'(1X,I4,1X,20(2X,F14.6))') k,Sdr(1:4),Sda(1:6),-vr(1:4),va(1:6)
                !write (iunascan,'(1X,I4,1X,20(2X,F14.6))') k,Sdr(1:4),Sda(1:6),-vr(1:4),va(1:6)/LPI*180.0_FREAL
             else
                call message(MESOUT,"[XY4C] ")
                call message_value(MESOUT,"[XY4C] Check OK at Solution:",isZ)
                call message(MESOUT,"[XY4C] ")
                write (iunscan,'(1X,I4,1X,20(2X,F14.6))') k,Sdr(1:4),Sda(1:6),vr(1:4),va(1:6)
                write (iunascan,'(1X,I4,1X,20(2X,F14.6))') k,Sdr(1:4),Sda(1:6),vr(1:4),va(1:6)/LPI*180.0_FREAL
                irc=get_cart()
                if (irc.ne.0) then
                   call message(MESOUT,"[XY4C] ")
                   call message(MESERRO,"[XY4C] Error in Cartesian routine.")
                   call message(MESOUT,"[XY4C] ")
                else
                   write (iunxyz,'(1X,I4)') 5
                   write (iunxyz,'(1X,F14.6,1X)') real(k)
                   write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'C',mxyz(1:3,1)*0.52917706E0_FREAL
                   write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,2)*0.52917706E0_FREAL
                   write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,3)*0.52917706E0_FREAL
                   write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,4)*0.52917706E0_FREAL
                   write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,5)*0.52917706E0_FREAL
                end if
             end if

          end do
          
       end do

    !==============================================================================================
    end if
    if (line(1:1).eq.'S') then
       call message(MESOUT,"[XY4C] ")
       call message(MESOUT,"[XY4C] Readed in Option S (Scan on one coordinate)")
       call message(MESOUT,"[XY4C] ")
       read(line,*) cscanopt,fscanvars(1),fscanpars(1:3)
       j=int((fscanpars(2)-fscanpars(1))/fscanpars(3))+1
       if (fscanvars(1).lt.5) then
          fscancoord => Sdr(fscanvars(1))
       else
          fscancoord => Sda(fscanvars(1)-4)
          if ((fscanvars(1)-4).eq.6) then
             call message(MESERRO,"Scan not allowed in redundand coordinate Sda(6)[Sr] at present.")
             stop 3
          end if
       end if
       call message_value(MESOUT,"Number of Points in Scan:",j)
       fscancoord=fscanpars(1)
       do k=1,j
          call message(MESOUT,"[XY4C] ")
          call message_value(MESOUT,"Scan Step:",k)
          call message(MESOUT,"[XY4C] ")
          !
          ! Eval redundancy
          call eval_sr()
          call get_reZ()
          write (*,'("#[x-out-comment] Generate Sr redundancies.")')
          if (iZ.gt.0) then
             write (*,'(A13,1X,4(2X,F12.4,1X,"i(",F12.4,")"))') '#Z Re Solut.:',ZSda(1:iZ)
          else
             write (*,'("#Z Re Solut.: NONE")')
             call message(MESERRO," No Available redundant solution. Keep going but there is an error!")
             call message(MESOUT," No Available redundant solution. Keep going but there is an error!")
          end if
          write (*,'(A13,1X,4(2X,F12.4,1X,"i(",F12.4,")"))') '#Z solutions:',ZOSda
          call message(MESOUT,"[XY4C] ")
          call message_value(MESOUT,"[XY4C] Start run on Available Sr solutions. N:",iZ)
          call message(MESOUT,"[XY4C] ")
          ircr=-1
          do isZ=1,iZ
             call message(MESOUT,"[XY4C] ")
             call message_value(MESOUT,"[XY4C] Checking Solution:",isZ)
             call message(MESOUT,"[XY4C] ")
             Sda(6)=ZSda(isZ)
             write (*,'("[x--step-xy4-symmcoord] Sol:",1X,I4,1X,"Sr Value:",1X,F11.4)') isZ,Sda(6)
             write (*,'("#",1X,10(1X,A11))') 'S1','S2x','S2y','S2z','S2a','S2b','S4x','S4y','S4z','Sr'
             write (*,'("#",1X,10(1X,F11.4))') Sdr(1:4),Sda(1:6)
             call get_ra
             irc=do_checks()
             if (irc.ne.0) then
                call message(MESERRO,"[XY4C] Check not Passed.")
                write (iunscan,'(1X,I4,1X,11(2X,F14.6))') k,fscancoord,-vr(1:4),va(1:6)
                !write (iunascan,'(1X,I4,1X,11(2X,F14.6))') k,fscancoord,-vr(1:4),va(1:6)/LPI*180.0_FREAL
             else
                write (iunscan,'(1X,I4,1X,11(2X,F14.6))') k,fscancoord,vr(1:4),va(1:6)
                write (iunascan,'(1X,I4,1X,11(2X,F14.6))') k,fscancoord,vr(1:4),va(1:6)/LPI*180.0_FREAL
                irc=get_cart()
                if (irc.ne.0) then
                   call message(MESOUT,"[XY4C] ")
                   call message(MESERRO,"[XY4C] Error in Cartesian routine.")
                   call message(MESOUT,"[XY4C] ")
                else
                   write (iunxyz,'(1X,I4)') 5
                   !write (iunxyz,'(1X,A4,1X,I4)') 'Step',k
                   write (iunxyz,'(1X,F14.6)') fscancoord
                   write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'C',mxyz(1:3,1)*0.52917706E0_FREAL
                   write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,2)*0.52917706E0_FREAL
                   write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,3)*0.52917706E0_FREAL
                   write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,4)*0.52917706E0_FREAL
                   write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,5)*0.52917706E0_FREAL
                end if
             end if
          end do
          fscancoord=fscancoord+fscanpars(3)
       end do
    end if
    !
    irc=baseio_close(iunscan)
    irc=baseio_close(iunascan)
    irc=baseio_close(iunxyz)
!
!
  end subroutine do_scan
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H subroutine do_mul_scan
!H-----------------------------------------------------------------------------
!H
!
  subroutine do_mul_scan()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine 
!H
!H-----------------------------------------------------------------------------
!H
!
    filescan='mul_scan.dat'
    irc = baseio_open(iunscan,filescan,stat='UNKNOWN')
    if(irc.lt.0) then
       LOCXY4coord_error = 3
       call message(MESERRO,"[XY4C] Cannot Open MUL_SCAN File")
       stop 1
    end if
    fileascan='ascan.dat'
    irc = baseio_open(iunascan,fileascan,stat='UNKNOWN')
    if(irc.lt.0) then
       LOCXY4coord_error = 3
       call message(MESERRO,"[XY4C] Cannot Open AMUL_SCAN File")
       stop 1
    end if
    filexyz='xyz.dat'
    irc = baseio_open(iunxyz,filexyz,stat='UNKNOWN')
    if(irc.lt.0) then
       LOCXY4coord_error = 3
       call message(MESERRO,"[XY4C] Cannot Open XYZ File")
       stop 1
    end if
    !
    irc = line_getline(iuninput,line,3)
    if (line(1:1).eq.'E') then
       call message(MESOUT,"[XY4C] Stop request from mul_scan.")
       stop 0
    end if
    if (line(1:1).eq.'S') then
       call message(MESOUT,"[XY4C] S not implemented.")
    end if
    if (line(1:1).eq.'M') then
    end if
    !==============================================================================================
    !==============================================================================================
    !==============================================================================================
    if (line(1:1).eq.'R') then
       read(line,*) cscanopt,fscanvars(1),fscanpars(1:3)
       j=int((fscanpars(2)-fscanpars(1))/fscanpars(3))+1
       if (fscanvars(1).lt.5) then
          fscancoord => Sdr(fscanvars(1))
       else
          fscancoord => Sda(fscanvars(1)-4)
       end if
       write (*,*) 'Num Steps=',j
       fscancoord=fscanpars(1)
       do k=1,j
          call get_ra
          irc=do_checks()
          if (irc.ne.0) then
             call message(MESERRO,"[XY4C] Check not Passed.")
             !stop 3
             write (iunscan,'(1X,I4,1X,11(2X,F14.6))') k,fscancoord,-vr(1:4),va(1:6)
             write (iunascan,'(1X,I4,1X,11(2X,F14.6))') k,fscancoord,-vr(1:4),va(1:6)/LPI*180.0_FREAL
          else
             write (iunscan,'(1X,I4,1X,11(2X,F14.6))') k,fscancoord,vr(1:4),va(1:6)
             write (iunascan,'(1X,I4,1X,11(2X,F14.6))') k,fscancoord,vr(1:4),va(1:6)/LPI*180.0_FREAL
             irc=get_cart()
             if (irc.ne.0) then
                call message(MESOUT,"[XY4C] ")
                call message(MESERRO,"[XY4C] Error in Cartesian routine.")
                call message(MESOUT,"[XY4C] ")
             else
                write (iunxyz,'(1X,I4)') 5
                !write (iunxyz,'(1X,A4,1X,I4)') 'Step',k
                write (iunxyz,'(1X,F14.6)') fscancoord
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'C',mxyz(1:3,1)*0.52917706E0_FREAL
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,2)*0.52917706E0_FREAL
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,3)*0.52917706E0_FREAL
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,4)*0.52917706E0_FREAL
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,5)*0.52917706E0_FREAL
             end if
          end if
          fscancoord=fscancoord+fscanpars(3)
       end do
    end if
    !
    irc=baseio_close(iunscan)
    irc=baseio_close(iunascan)
    irc=baseio_close(iunxyz)
!
!
  end subroutine do_mul_scan
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H subroutine do_revscan
!H-----------------------------------------------------------------------------
!H
!
  subroutine do_revscan()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine 
!H
!H-----------------------------------------------------------------------------
!H
!
    filescan='scan.dat'
    irc = baseio_open(iunscan,filescan,stat='UNKNOWN')
    if(irc.lt.0) then
       LOCXY4coord_error = 3
       call message(MESERRO,"[XY4C] Cannot Open SCAN File")
       stop 1
    end if
    fileascan='ascan.dat'
    irc = baseio_open(iunascan,fileascan,stat='UNKNOWN')
    if(irc.lt.0) then
       LOCXY4coord_error = 3
       call message(MESERRO,"[XY4C] Cannot Open ASCAN File")
       stop 1
    end if
    filexyz='xyz.dat'
    irc = baseio_open(iunxyz,filexyz,stat='UNKNOWN')
    if(irc.lt.0) then
       LOCXY4coord_error = 3
       call message(MESERRO,"[XY4C] Cannot Open XYZ File")
       stop 1
    end if
    !
    irc = line_getline(iuninput,line,3)
    if (line(1:1).eq.'E') then
       call message(MESOUT,"[XY4C] Stop request from scan.")
       stop 0
    end if
    if (line(1:1).eq.'S') then
       call message(MESOUT,"[XY4C] S not implemented.")
    end if
    if (line(1:1).eq.'M') then
       call message(MESOUT,"[XY4C] M not implemented.")
    end if
    if (line(1:1).eq.'R') then
       read(line,*) cscanopt,fscanvars(1),fscanpars(1:3)
       j=int((fscanpars(2)-fscanpars(1))/fscanpars(3))+1
       if (fscanvars(1).lt.5) then
          fscancoord => dr(fscanvars(1))
       else
          fscancoord => da(fscanvars(1)-4)
       end if
       write (*,*) 'Num Steps=',j
       fscancoord=fscanpars(1)
       do k=1,j
          call get_symc
          svr(1:4)=dr(1:4)
          sva(1:6)=da(1:6)
          ! Get full coordinates i.e. add displacments
          !--------------------------------------------
          vr(1:4)=xy(1:4)+dr(1:4)
          call message(MESOUT,"Bond lengths from internal displacments:")
          write (*,'(1X,A1,1X,4(2X,F14.6))') ' ',vr(1:4)
          !
          call message(MESWARN,"Assuming initial angles 109. NOT evaluated now...:")
          va(1:6)=((109.4712_FREAL/180.0_FREAL)*LPI)
          va(1:6)=va(1:6)+da(1:6)
          call message(MESOUT,"Angles (GRAD) from internal displacments:")
          write (*,'(1X,A1,1X,6(2X,F14.6))') ' ',va(1:6)/LPI*180.0_FREAL
          call get_ra
          do i=1,4
             if ((svr(i)-dr(i)).gt.XY4LEPS) then
                call message_value(MESERRO,"Double inversion doesnt Match. R Coord:",i)
             end if
             if ((sva(i)-da(i)).gt.XY4LEPS) then
                call message_value(MESERRO,"Double inversion doesnt Match. A Coord:",i)
                call message_value(MESERRO,"Value Origi:",sva(i))
                call message_value(MESERRO,"Value After:",da(i))
                call message_value(MESERRO,"Delta Value:",sva(i)-da(i))
             end if
          end do
          do i=5,6
             if ((sva(i)-da(i)).gt.XY4LEPS) then
                call message_value(MESERRO,"Double inversion doesnt Match. A Coord:",i)
                call message_value(MESERRO,"Value Origi:",sva(i))
                call message_value(MESERRO,"Value After:",da(i))
                call message_value(MESERRO,"Delta Value:",sva(i)-da(i))
             end if
          end do
          irc=do_checks()
          if (irc.ne.0) then
             call message(MESERRO,"[XY4C] Check not Passed.")
             !stop 3
             write (iunscan,'(1X,I4,1X,11(2X,F14.6))') k,fscancoord,-Sdr(1:4),Sda(1:6)
             write (iunascan,'(1X,I4,1X,11(2X,F14.6))') k,fscancoord,-Sdr(1:4),Sda(1:6)/LPI*180.0_FREAL
          else
             write (iunscan,'(1X,I4,1X,11(2X,F14.6))') k,fscancoord,Sdr(1:4),Sda(1:6)
             write (iunascan,'(1X,I4,1X,11(2X,F14.6))') k,fscancoord,Sdr(1:4),Sda(1:6)/LPI*180.0_FREAL
             irc=get_cart()
             if (irc.ne.0) then
                call message(MESOUT,"[XY4C] ")
                call message(MESERRO,"[XY4C] Error in Cartesian routine.")
                call message(MESOUT,"[XY4C] ")
             else
                write (iunxyz,'(1X,I4)') 5
                !write (iunxyz,'(1X,A4,1X,I4)') 'Step',k
                write (iunxyz,'(1X,F14.6)') fscancoord
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'C',mxyz(1:3,1)*0.52917706E0_FREAL
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,2)*0.52917706E0_FREAL
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,3)*0.52917706E0_FREAL
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,4)*0.52917706E0_FREAL
                write (iunxyz,'(1X,A3,2X,3(2X,F14.6))') 'H',mxyz(1:3,5)*0.52917706E0_FREAL
             end if
          end if
          fscancoord=fscancoord+fscanpars(3)
       end do
    end if
    !
    irc=baseio_close(iunscan)
    irc=baseio_close(iunascan)
    irc=baseio_close(iunxyz)
!
!
  end subroutine do_revscan
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H subroutine eval_sr()
!H-----------------------------------------------------------------------------
!H
!
  subroutine eval_sr()
    USE constants_NSWC
    IMPLICIT NONE

    !
    real(FREAL), dimension(0:4) :: A
    !
    INTERFACE
       SUBROUTINE qtcrt (b, z)
         USE constants_NSWC
         IMPLICIT NONE
         REAL (dp), INTENT(IN)     :: b(:)
         COMPLEX (dp), INTENT(OUT) :: z(:)
       END SUBROUTINE qtcrt
    END INTERFACE
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine 
!H
!H-----------------------------------------------------------------------------
!H
!
    ! A coefficients from JCP 118 (2003) 6260 - Wang, Carrington
    A(4)= -1.0_FREAL/12.0_FREAL

    A(3)=  (2.0_FREAL*sqrt(6.0_FREAL))/9.0_FREAL

    A(2)= -1.0_FREAL+ 0.5_FREAL * (Sda(1)*Sda(1)+Sda(2)*Sda(2)) + (1.0_FREAL/6.0_FREAL)*(Sda(3)*Sda(3)+Sda(4)*Sda(4)+Sda(5)*Sda(5))

    A(1)=      - (sqrt(6.0_FREAL)/3.0_FREAL)*Sda(2)*(Sda(3)*Sda(3)-Sda(4)*Sda(4))
    A(1)= A(1) + (sqrt(2.0_FREAL)/3.0_FREAL)*Sda(1)*(Sda(3)*Sda(3)+Sda(4)*Sda(4)-2.0_FREAL*Sda(5)*Sda(5))
    A(1)= A(1) + (sqrt(2.0_FREAL)/3.0_FREAL)*Sda(1)*(Sda(1)*Sda(1)-3.0_FREAL*Sda(2)*Sda(2))
    A(1)= A(1) - (sqrt(6.0_FREAL)/3.0_FREAL)*(Sda(1)*Sda(1)+Sda(2)*Sda(2))

    A(0)=        0.25_FREAL * (Sda(3)*Sda(3)+Sda(4)*Sda(4)+Sda(5)*Sda(5))
    A(0)= A(0) - 0.5_FREAL * (Sda(3)*Sda(3)*Sda(4)*Sda(4)+Sda(4)*Sda(4)*Sda(5)*Sda(5)+Sda(5)*Sda(5)*Sda(3)*Sda(3))
    A(0)= A(0) - (1.0_FREAL/6.0_FREAL) * (Sda(1)*Sda(1)-3.0_FREAL*Sda(2)*Sda(2)) * Sda(5)*Sda(5)
    A(0)= A(0) + (1.0_FREAL/3.0_FREAL) * Sda(1)*Sda(1)*(Sda(3)*Sda(3)+Sda(4)*Sda(4))
    A(0)= A(0) + (1.0_FREAL/sqrt(3.0_FREAL)) * Sda(1)*Sda(2)*(Sda(3)*Sda(3)-Sda(4)*Sda(4))
    A(0)= A(0) + 2.0_FREAL*sqrt(2.0_FREAL)*Sda(3)*Sda(4)*Sda(5)
    A(0)= A(0) + (2.0_FREAL/(3.0_FREAL*sqrt(3.0_FREAL))) * Sda(1) * (Sda(1)*Sda(1)-3.0_FREAL*Sda(2)*Sda(2))
    A(0)= A(0) - (Sda(3)*Sda(3)+Sda(4)*Sda(4)+Sda(5)*Sda(5))
    A(0)= A(0) - (Sda(1)*Sda(1)+Sda(2)*Sda(2))
    A(0)= A(0) + 1
    !
    call qtcrt(A, ZSda)
    ZOSda(:)=ZSda(:)
    !
!
!
  end subroutine eval_sr
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H subroutine get_reZ()
!H-----------------------------------------------------------------------------
!H
!
  subroutine get_reZ()
    implicit none
    integer(FINT) :: iexc,l,m,i2Z
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine 
!H
!H-----------------------------------------------------------------------------
!H
!
    iZ=0
    ! Put Imag values at end of vector
    do
       iexc=0
       do l=1,3
          if (abs(aimag(ZSda(l))).gt.abs(aimag(ZSda(l+1)))) then
             Ztmp=ZSda(l)
             ZSda(l)=ZSda(l+1)
             ZSda(l+1)=Ztmp
             iexc=1
          end if
       end do
       if(iexc.eq.0) then
          exit
       end if
    end do
    ! Find last real
    do l=1,4
       if (abs(aimag(ZSda(l))).lt.XY4LEPS) then
          iZ=l
       end if
    end do
    !
    ! On real solutions get degeneracy
    !
    ! first sort it.
    do
       iexc=0
       do l=1,iZ-1
          if (dble(ZSda(l)).gt.dble(ZSda(l+1))) then
             Ztmp=ZSda(l)
             ZSda(l)=ZSda(l+1)
             ZSda(l+1)=Ztmp
             iexc=1
          end if
       end do
       if(iexc.eq.0) then
          exit
       end if
    end do
    !return
    !
    ! now contraction
    do
       iexc=0
       !write (*,*) 'Z contraction loop over',1,iZ-1
       do l=1,iZ-1
          !write (*,*) 'Z contraction check',l,l+1
          if (abs(dble(ZSda(l))-dble(ZSda(l+1))).lt.XY4LEPS) then
             do m=l,iZ-1
                !write (*,*) 'Z contraction exchange',m,m+1
                ZSda(m)=ZSda(m+1)
             end do
             iexc=1
             i2Z=iZ-1
             exit
          end if
       end do
       if (iexc.eq.1) iZ=i2Z
       if(iexc.eq.0.or.iZ.eq.1) then
          exit
       end if
    end do

    !
    !

  end subroutine get_reZ
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H subroutine getbonds(n,lvx,lvr,debug)
!H-----------------------------------------------------------------------------
!H
!
  subroutine getbonds(n,lvx,lvr,debug)
    !
    implicit none
    integer(FINT), intent(in) :: n
    real(FREAL), dimension(:,:), intent(in) :: lvx
    real(FREAL), dimension(:), intent(out) :: lvr
    !
    integer(FINT) :: li
    !
    logical, optional :: debug
    logical :: ldebug
    logical :: lerror
    lerror=.FALSE.
    if (present(debug)) then
       ldebug=debug
    else
       ldebug=.TRUE.
    end if
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine returns bond lenghts from central atom at origin.
!H
!H-----------------------------------------------------------------------------
!H
!
    ! checks
    !-------
    if (ldebug) then
       if (size(lvx,1).ne.3) then
          call message(MESERRO,"[getbonds]: we need 3D cartesian input coordinates.")
          lerror=.TRUE.
       end if
       if (size(lvx,2).lt.n) then
          call message(MESERRO,"[getbonds]: Inconsistent number of atoms with respect cartesians coordinates.")
          lerror=.TRUE.
       end if
       if (size(lvx,2).gt.n) then
          call message(MESWARN,&
               &"[getbonds]: Inconsistent number of atoms with respect cartesians coordinates. But we can do the job.")
       end if
       if (size(lvr).lt.n) then
          call message(MESERRO,"[getbonds]: Inconsistent number of atoms with respect bonds vector.")
          lerror=.TRUE.
       end if
       if (size(lvr).gt.n) then
          call message(MESWARN,"[getbonds]: Inconsistent number of atoms with respect bonds vector. But we can do the job.")
       end if
       if (lerror) then
          call message(MESWARN,"[getbonds]: One or more errors in this routine. We do stop.")
          stop 3
       end if
    end if
    !
    ! END CHECKS
    !-----------
    !
    ! routine
    !--------
    do li=1,n
       lvr(li)=sqrt(sum(lvx(1:3,li)*lvx(1:3,li)))
    end do
!
!
  end subroutine getbonds
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H subroutine getangles(n,lvx,lvr,lva,debug)
!H-----------------------------------------------------------------------------
!H
!
  subroutine getangles(n,lvx,lvr,lva,debug)
    !
    implicit none
    integer(FINT), intent(in) :: n
    real(FREAL), dimension(:,:), intent(in) :: lvx
    real(FREAL), dimension(:), intent(in) :: lvr
    real(FREAL), dimension(:), intent(out) :: lva
    !
    integer(FINT) :: li,la,lb
    !
    logical, optional :: debug
    logical :: ldebug
    logical :: lerror
    lerror=.FALSE.
    if (present(debug)) then
       ldebug=debug
    else
       ldebug=.TRUE.
    end if
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine returns angles between vectors of catersian atoms from origin.
!H Order is defined by the loop check it or btw a12 a13 a14 a23 a24 a34 on 4 
!H
!H-----------------------------------------------------------------------------
!H
!
    ! checks
    !-------
    if (ldebug) then
       if (size(lvx,1).ne.3) then
          call message(MESERRO,"[getangles]: we need 3D cartesian input coordinates.")
          lerror=.TRUE.
       end if
       if (size(lvx,2).lt.n) then
          call message(MESERRO,"[getangles]: Inconsistent number of atoms with respect cartesians coordinates.")
          lerror=.TRUE.
       end if
       if (size(lvx,2).gt.n) then
          call message(MESWARN,&
               &"[getangles]: Inconsistent number of atoms with respect cartesians coordinates. But we can do the job.")
       end if
       if (size(lvr).lt.n) then
          call message(MESERRO,"[getangles]: Inconsistent number of atoms with respect bonds vector.")
          lerror=.TRUE.
       end if
       if (size(lvr).gt.n) then
          call message(MESWARN,"[getangles]: Inconsistent number of atoms with respect bonds vector. But we can do the job.")
       end if
       if (size(lva).lt.(n*n-n)/2) then
          call message(MESERRO,"[getangles]: Inconsistent number of atoms with respect angles vector.")
          lerror=.TRUE.
       end if
       if (size(lva).gt.(n*n-n)/2) then
          call message(MESERRO,"[getangles]: Inconsistent number of atoms with respect angles vector. But we can do the job.")
          lerror=.TRUE.
       end if
       if (lerror) then
          call message(MESWARN,"[getangles]: One or more errors in this routine. We do stop.")
          stop 3
       end if
    end if
    !
    ! END CHECKS
    !-----------
    !
    ! do it!
    !-------
    li=0
    do la=1,n-1
       do lb=la+1,n
          li=li+1
          lva(li)=acos(sum(lvx(1:3,la)*lvx(1:3,lb))/(lvr(la)*lvr(lb)))
       end do
    end do
!
!
  end subroutine getangles
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H subroutine XXXX
!H-----------------------------------------------------------------------------
!H
!
  subroutine XXXX
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine 
!H
!H-----------------------------------------------------------------------------
!H
!
  end subroutine XXXX
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H Function XY4coord_init()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function XY4coord_init()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the XY4coord module.
!H
!H-----------------------------------------------------------------------------
!H
!
    XY4coord_init = 0
    !
    ! Initialize messages and ourself
    !--------------------------------
    call messages_init()
    !
    ! File I/O stuff initialization
    !------------------------------
    irc = baseio_init()
    if(irc.lt.0) then
       call message(MESERRO,"[XY4C] Cannot initialize baseIO.")
       XY4coord_init = -1
       return
    end if
    !
    !
    ! Parse command line
    !-------------------
    call pcmd_iargc(nargs)
    allocate(cargs(nargs),STAT=irc)
    if(.not.irc.eq.0) then
       call message(MESERRO,"[XY4C] Cannot allocate for command line.")
       XY4coord_init = -1
       return
    end if
    call pcmd_getarg(nargs,cargs)
    call pcmd_getio(nargs,cargs,fileinput,fileoutput)
    !
    ! We open the input and output files
    !-----------------------------------
    fileinput = trim(fileinput)
    if(scan(fileinput," ").eq.0) then
       iuninput = 5
       fileinput = "stdin"
    else
       irc = baseio_open(iuninput,fileinput,stat='OLD')
       if(irc.lt.0) then
          LOCXY4coord_error = 3
          call message(MESERRO,"[XY4C] Cannot Open Input File")
          XY4coord_init = -1
          return
       end if
    end if
    !
    fileoutput = trim(fileoutput)
    if(scan(fileoutput," ").eq.0) then
       iunoutput = 6
       fileoutput = "stdout"
    else
       irc = baseio_open(iunoutput,fileoutput,stat='UNKNOWN')
       if(irc.lt.0) then
          LOCXY4coord_error = 3
          call message(MESERRO,"[XY4C] Cannot Open Output File")
          XY4coord_init = -1
          return
       end if
    end if
    !
    ! We need the molecule section
    !-----------------------------
    irc=molecule_init(.true.)
    if(irc.lt.0) then
       call message(MESERRO,"[XY4C] Cannot initialize molecule module..")
       XY4coord_init = -1
       return
    end if
    !
    return
    !
  end function XY4coord_init
!
!H
!H-----------------------------------------------------------------------------
!H logical function XY4coord_isinit()
!H-----------------------------------------------------------------------------
!H
!
  logical function XY4coord_isinit()
!
!H
!H-----------------------------------------------------------------------------
!H We merelyreturn the initialization flag.
!H-----------------------------------------------------------------------------
!H
!
    !
    XY4coord_isinit = LOCXY4coord_isinit
    !
  end function XY4coord_isinit
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H character(FLCHARS) function XY4coord_geterror()
!H-----------------------------------------------------------------------------
!H
!
  character(FLCHARS) function XY4coord_geterror(code)
    !
    integer(FINT), intent(in) :: code
    !
!
!H
!H-----------------------------------------------------------------------------
!H 
!H Return the error code of the program in Object oriented
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
       XY4coord_geterror = "No Errors."
    case(1_FINT)
       XY4coord_geterror = "Not documented error or generic error.."
    case(2_FINT)
       XY4coord_geterror = "Double Error: Why did you get this message?."
    case default
       XY4coord_geterror = "Not documented error or generic error.."
    end select
    !
    return
    !
  end function XY4coord_geterror
!
!H
!H-----------------------------------------------------------------------------
!H
!

  !
  !
end program XY4coord
