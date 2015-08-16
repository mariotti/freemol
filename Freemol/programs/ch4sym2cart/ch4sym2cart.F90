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
!H PROGRAM ch4sym2cart  Frimol by F.Mariotti:
!H-----------------------------------------------------------------------------
!H $Id: ch4sym2cart.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
!H-----------------------------------------------------------------------------
!H This File is part of the Frimol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
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
program ch4sym2cart
  !
  use vartypes
  use messages
  use pcmdline
  use linetools
  use baseio
  use osec
  use molecule
  !use ch4s2cmod
  !
  implicit none
  !
  ! DECLARATIONS:
  !--------------------
  
  ! Procedures JUST a List dont remove comments
  !------------------
  !integer(FINT) :: ch4sym2cart_init
  !logical ::  ch4sym2cart_isinit
  !integer(FINT) :: ch4sym2cart_geterror
  !integer(FINT), public :: ch4sym2cart_
  !
  !
  ! LOCAL DECLARATIONS:
  !--------------------
  real(FREAL), parameter :: LPI = 3.1415926535897932384626433832795028841971694_FREAL
  integer(FINT) :: i,j,k
  integer(FINT) :: irc
  character(FLCHARS) :: params
  character(FLCHARS) :: line
  real(FREAL) :: ftmp,frtmp,f1d144,fr2,fr3,f1d72,fgamma,fa109
  !
  !Coordinates
  !-----------
  real(FREAL), dimension(3,5), save :: mxyz
  real(FREAL), dimension(3,5), save :: rmxyz
  real(FREAL), dimension(10), save :: vsymcrd
  real(FREAL), dimension(10), save :: vzmatcrd
  !
  ! Files
  !------
  integer(FINT), save :: iuninput
  character(FLCHARS), save :: fileinput
  integer(FINT), save :: iunoutput
  character(FLCHARS), save :: fileoutput
  !
  ! Command Line
  !-------------
  integer(FINT) :: nargs
  character(FLCHARS),dimension(:), allocatable :: cargs
  !
  !
  ! error variables
  !----------------
  logical, save :: LOCch4sym2cart_isinit
  integer(FINT), save :: LOCch4sym2cart_error
  character(FLCHARS), save :: LOCch4sym2cart_error_message
  !
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Here we start
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !
  ! Initialize messages and ourself
  !--------------------------------
  call messages_init()
  irc=ch4sym2cart_init()
  !
  ! File I/O stuff initialization
  !------------------------------
  irc = baseio_init()
  if(irc.lt.0) then
     call message(MESERRO,"[CH4S2C] Cannot initialize baseIO.")
     stop 1
  end if
  !
  !
  ! Parse command line
  !-------------------
  call pcmd_iargc(nargs)
  allocate(cargs(nargs),STAT=irc)
  if(.not.irc.eq.0) then
     call message(MESERRO,"[CH4S2C] Cannot allocate for command line.")
     stop 1
     !return
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
        LOCch4sym2cart_error = 3
        call message(MESERRO,"[CH4S2C] Cannot Open Input File")
        stop 1
        !return
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
        LOCch4sym2cart_error = 3
        call message(MESERRO,"[CH4S2C] Cannot Open Output File")
        stop 1
        !return
     end if
  end if
  !
  ! Start reading
  !--------------
  !
  ! We need the molecule section
  !-----------------------------
  irc=molecule_init(.true.)
  if(irc.lt.0) then
     call message(MESERRO,"[CH4S2C] Cannot initialize molecule module..")
     stop 1
     !return
  end if
  rewind(iuninput)
  irc = osec_set(iuninput,'molecule',params)
  if(irc.lt.0) then
     call message(MESERRO,"[CH4S2C] No molecule Section.")
     stop 1
     !return
  end if
  call molecule_read(iuninput,params)
  call molecule_print(iunoutput)
  irc=molecule_getnat()
  if(irc.ne.5) then
     call message(MESERRO,"[CH4S2C] This is CH4 program: we need exactly 5 atoms!.")
     stop 1
  end if
  !
  call message(MESOUT,"[CH4S2C] Eckard reference system readed in..")
  !
  ! We read the symmetric coordinate system
  !----------------------------------------
  rewind(iuninput)
  irc = osec_set(iuninput,'x-ch4-symmcoord',params)
  if(irc.lt.0) then
     call message(MESERRO,"[CH4S2C] No x-ch4-symmcoord Section.")
     stop 1
     !return
  end if
  !
  !
  irc = line_getline(iuninput,line,3)
  read(line,*) vsymcrd(1:10)
  !write (*,*) vsymcrd(1:10)
  !
  !
  !r1,r2,r3,r4
  !-----------
  vzmatcrd(1)=0.125_FREAL*(vsymcrd(1)+vsymcrd(4)-vsymcrd(5)+vsymcrd(6))
  vzmatcrd(2)=0.125_FREAL*(vsymcrd(1)-vsymcrd(4)+vsymcrd(5)+vsymcrd(6))
  vzmatcrd(3)=0.125_FREAL*(vsymcrd(1)+vsymcrd(4)+vsymcrd(5)-vsymcrd(6))
  vzmatcrd(4)=0.125_FREAL*(vsymcrd(1)-vsymcrd(4)-vsymcrd(5)-vsymcrd(6))
  write (*,*) 'R', vzmatcrd(1:4)

  !
  ! a12,a13,a14,a23,a24,a34
  !------------------------
  f1d144 = 1.0_FREAL/144.0_FREAL
  f1d72  = 1.0_FREAL/72.0_FREAL
  fr2    = sqrt(2.0_FREAL)
  fr3    = sqrt(3.0_FREAL)
  !a12
  vzmatcrd(5)=f1d144*fr2 * ( fr2 * fr3 * vsymcrd(2) - 9_FREAL * fr2 * vsymcrd(3)&
       & + 36_FREAL * vsymcrd(7) - 72_FREAL * vsymcrd(9)&
       & - 4_FREAL * fr3 * vsymcrd(10))
  !a13
  vzmatcrd(6)= - f1d72 * fr3 * vsymcrd(2) + 0.125_FREAL * vsymcrd(3) &
       & + 0.25_FREAL * fr2 * vsymcrd(7) &
       & + (1.0_FREAL/36.0_FREAL) * fr2 * fr3 * vsymcrd(10)
  !a14
  vzmatcrd(7)=f1d144*fr2 * ( fr2 * fr3 * vsymcrd(2) &
       & - 36.0_FREAL * vsymcrd(8) &
       & - 4.0_FREAL * fr3 * vsymcrd(10))
  !a23
  vzmatcrd(8)= - f1d72 * fr3 * (vsymcrd(2) + 3.0_FREAL * fr3 * vsymcrd(3)&
       & - 6.0_FREAL * fr2 * fr3 * vsymcrd(8)&
       & - 2.0_FREAL * fr2 * vsymcrd(10))
  !a24
  vzmatcrd(9)=f1d144*fr2 * ( fr2 * fr3 * vsymcrd(2) - 9_FREAL * fr2 * vsymcrd(3)&
       & + 36_FREAL * vsymcrd(7)&
       & - 4.0_FREAL * fr3 * vsymcrd(10))
  !a34
  vzmatcrd(10)= f1d72 * fr3 * (5.0_FREAL * vsymcrd(2) - 3.0_FREAL * fr3 * vsymcrd(3)&
       & + 6.0_FREAL * fr2 * fr3 * vsymcrd(7) - 12.0_FREAL * fr2 * fr3 * vsymcrd(9)&
       & + 2.0_FREAL * fr2 * vsymcrd(10))
  !
  write (*,*) 'A', vzmatcrd(5:10)
  fa109=(109.4712_FREAL/180.0_FREAL)*LPI
  !
  !H1
  do i=5,10
     vzmatcrd(i)=vzmatcrd(i)+((109.4712_FREAL/180.0_FREAL)*LPI)
     if(vzmatcrd(i).lt.0.0_FREAL) then
        call message(MESERRO,"[CH4S2C] Inconsistent negative angle")
        call message_value(MESERRO,"[CH4S2C] term",i)
        stop 3
     endif
  end do
  write (*,*) 'A', vzmatcrd(5:10)
  ftmp=sqrt(sum(molecule_xyz(1:3,2)*molecule_xyz(1:3,2)))+vzmatcrd(1)
  frtmp=sqrt(sum(molecule_xyz(1:3,2)*molecule_xyz(1:3,2)))
  if(ftmp.lt.0.0001_FREAL) then
     call message(MESERRO,"[CH4S2C]: Very Small bond distance.")
  end if
  mxyz(1,2)=ftmp
  mxyz(2,2)=0.0_FREAL
  mxyz(3,2)=0.0_FREAL
  rmxyz(1,2)=frtmp
  rmxyz(2,2)=0.0_FREAL
  rmxyz(3,2)=0.0_FREAL
  write(line,'(1X,A18,1X,3(2X,F12.4))') 'NEW Coordinates H1',mxyz(1:3,2)
  call message(MESOUT,line)
  !H2
  ftmp=sqrt(sum(molecule_xyz(1:3,3)*molecule_xyz(1:3,3)))+vzmatcrd(2)
  frtmp=sqrt(sum(molecule_xyz(1:3,3)*molecule_xyz(1:3,3)))
  if(ftmp.lt.0.0001_FREAL) then
     call message(MESERRO,"[CH4S2C]: Very Small bond distance.")
  end if
  if(vzmatcrd(5).gt.LPI) then
     call message(MESWARN,"[CH4S2C]: BIG A12 Angle possible internal inversion.")
  end if
  mxyz(1,3)=ftmp*cos(vzmatcrd(5))
  mxyz(2,3)=abs(ftmp*sin(vzmatcrd(5)))
  mxyz(3,3)=0.0_FREAL
  rmxyz(1,3)=frtmp*cos(fa109)
  rmxyz(2,3)=abs(frtmp*sin(fa109))
  rmxyz(3,3)=0.0_FREAL
  write(line,'(1X,A18,1X,3(2X,F12.4))') 'NEW Coordinates H2',mxyz(1:3,3)
  call message(MESOUT,line)
  !H3
  ftmp=sqrt(sum(molecule_xyz(1:3,4)*molecule_xyz(1:3,4)))+vzmatcrd(3)
  frtmp=sqrt(sum(molecule_xyz(1:3,4)*molecule_xyz(1:3,4)))
  if((vzmatcrd(5)+vzmatcrd(6)+vzmatcrd(8)).gt.(2.0_FREAL*LPI)) then
     call message(MESERRO,"[CH4S2C]: Not realistic internal angles for H123.")
     stop 3
  end if
  if(vzmatcrd(6).gt.LPI) then
     call message(MESWARN,"[CH4S2C]: BIG A13 Angle possible internal inversion.")
  end if
  if(vzmatcrd(8).gt.LPI) then
     call message(MESWARN,"[CH4S2C]: BIG A23 Angle possible internal inversion.")
  end if
  fgamma=acos((cos(vzmatcrd(8))-cos(vzmatcrd(5))*cos(vzmatcrd(6)))/(sin(vzmatcrd(5))*sin(vzmatcrd(6))))
  call message_value(MESOUT,"H3 Gamma value:",fgamma/LPI*180.0_FREAL)
  if(fgamma.gt.LPI) then
     call message(MESWARN,"[CH4S2C]: BIG Gamma Angle at H3 possible internal inversion.")
  end if
  mxyz(1,4)=ftmp*cos(vzmatcrd(8))
  mxyz(2,4)=abs(ftmp*sin(vzmatcrd(8)))*cos(fgamma)
  mxyz(3,4)=abs(ftmp*sin(vzmatcrd(8)))*abs(sin(fgamma))
  rmxyz(1,4)=frtmp*cos(fa109)
  rmxyz(2,4)=abs(frtmp*sin(fa109))*cos(120.0_FREAL)
  rmxyz(3,4)=abs(frtmp*sin(fa109))*abs(sin(120.0_FREAL))
  write(line,'(1X,A18,1X,3(2X,F12.4))') 'NEW Coordinates H3',mxyz(1:3,4)
  call message(MESOUT,line)
  
  !H4
  ftmp=sqrt(sum(molecule_xyz(1:3,5)*molecule_xyz(1:3,5)))+vzmatcrd(4)
  frtmp=sqrt(sum(molecule_xyz(1:3,5)*molecule_xyz(1:3,5)))
  if((vzmatcrd(5)+vzmatcrd(7)+vzmatcrd(9)).gt.(2.0_FREAL*LPI)) then
     call message(MESERRO,"[CH4S2C]: Not realistic internal angles for H123.")
     stop 3
  end if
  if(vzmatcrd(7).gt.LPI) then
     call message(MESWARN,"[CH4S2C]: BIG A14 Angle possible internal inversion.")
  end if
  if(vzmatcrd(9).gt.LPI) then
     call message(MESWARN,"[CH4S2C]: BIG A24 Angle possible internal inversion.")
  end if
  fgamma=acos((cos(vzmatcrd(7))-cos(vzmatcrd(5))*cos(vzmatcrd(9)))/(sin(vzmatcrd(5))*sin(vzmatcrd(9))))
  call message_value(MESOUT,"H4 Gamma value:",fgamma/LPI*180.0_FREAL)
  if(fgamma.gt.LPI) then
     call message(MESWARN,"[CH4S2C]: BIG Gamma Angle at H4 possible internal inversion.")
  end if
  mxyz(1,5)=ftmp*cos(vzmatcrd(7))
  mxyz(2,5)=abs(ftmp*sin(vzmatcrd(7)))*cos(fgamma)
  mxyz(3,5)=abs(ftmp*sin(vzmatcrd(7)))*-abs(sin(fgamma))
  rmxyz(1,5)=frtmp*cos(fa109)
  rmxyz(2,5)=abs(frtmp*sin(fa109))*cos(120.0_FREAL)
  rmxyz(3,5)=abs(frtmp*sin(fa109))*-abs(sin(120.0_FREAL))
  write(line,'(1X,A18,1X,3(2X,F12.4))') 'NEW Coordinates H4',mxyz(1:3,5)
  call message(MESOUT,line)
  !
  ! Check consistency rotating back the molecule
  !---------------------------------------------
  write(line,'(1X,A18,1X,3(2X,F12.4))') 'Unchenged Coordinates H1',rmxyz(1:3,2)
  call message(MESOUT,line)
  write(line,'(1X,A18,1X,3(2X,F12.4))') 'Unchenged Coordinates H2',rmxyz(1:3,3)
  call message(MESOUT,line)
  write(line,'(1X,A18,1X,3(2X,F12.4))') 'Unchenged Coordinates H3',rmxyz(1:3,4)
  call message(MESOUT,line)
  write(line,'(1X,A18,1X,3(2X,F12.4))') 'Unchenged Coordinates H4',rmxyz(1:3,5)
  call message(MESOUT,line)

  !
  !
contains
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H Function ch4sym2cart_init()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function ch4sym2cart_init()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the ch4sym2cart module.
!H
!H-----------------------------------------------------------------------------
!H
!
    ch4sym2cart_init = 0
  end function ch4sym2cart_init
!
!H
!H-----------------------------------------------------------------------------
!H logical function ch4sym2cart_isinit()
!H-----------------------------------------------------------------------------
!H
!
  logical function ch4sym2cart_isinit()
!
!H
!H-----------------------------------------------------------------------------
!H We merelyreturn the initialization flag.
!H-----------------------------------------------------------------------------
!H
!
    !
    ch4sym2cart_isinit = LOCch4sym2cart_isinit
    !
  end function ch4sym2cart_isinit
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H character(FLCHARS) function ch4sym2cart_geterror()
!H-----------------------------------------------------------------------------
!H
!
  character(FLCHARS) function ch4sym2cart_geterror(code)
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
       ch4sym2cart_geterror = "No Errors."
    case(1_FINT)
       ch4sym2cart_geterror = "Not documented error or generic error.."
    case(2_FINT)
       ch4sym2cart_geterror = "Double Error: Why did you get this message?."
    case default
       ch4sym2cart_geterror = "Not documented error or generic error.."
    end select
    !
    return
    !
  end function ch4sym2cart_geterror
!
!H
!H-----------------------------------------------------------------------------
!H
!

  !
  !
end program ch4sym2cart
