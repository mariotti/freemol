!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H PROGRAM XY4PolySphere  Frimol by F.Mariotti:
!H-----------------------------------------------------------------------------
!H $Id: XY4PolySphere.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
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
program XY4PolySphere
  !
  use vartypes
  use messages
  use pcmdline
  use linetools
  use baseio
  use osec
  use molecule
  use mathtools
  use metpot4
  !
  implicit none
  !
  ! DECLARATIONS:
  !--------------------
  
  ! Procedures JUST a List dont remove comments
  !------------------
  !integer(FINT) :: XY4PolySphere_init
  !logical ::  XY4PolySphere_isinit
  !integer(FINT) :: XY4PolySphere_geterror
  !integer(FINT), public :: XY4PolySphere_
  !
  !
  ! LOCAL DECLARATIONS:
  !--------------------
  real(FREAL), parameter :: LPI = 3.1415926535897932384626433832795028841971694_FREAL
  real(FREAL), parameter :: LEPS = 1.0E-8_FREAL
  !
  ! Dummy indexes
  !--------------
  integer(FINT) :: i,j,k
  !
  ! return code from functions and the running one
  !-----------------------------------------------
  integer(FINT) :: irc
  !integer(FINT) :: irc,ircr
  character(FLCHARS) :: params
  character(FLCHARS) :: line
  !character(FLCHARS) :: cscanopt
  !real(FREAL),dimension(10) :: fscanpars
  !real(FREAL), pointer :: fscancoord
  !
  ! Num scans for polysphericals
  integer(FINT) :: iscan
  ! Num of generated coordinates: both gengrid and genrandom
  !---------------------------------------------------------
  integer(FINT) :: igenpt
  !
  !Coordinates
  !-----------
  ! Cartesian Coordinates
  ! Convention: X is atom number 1 and the 4 Y are atoms number 2-5
  !----------------------------------------------------------------
  real(FREAL), dimension(3,5), save :: mxyz
  !
  ! Catesian for Euler Routine
  !---------------------------
  real(FREAL), dimension(15), save :: vxyz,vexyz
  !
  ! Angle guess for Euler
  !---------------------
  real(FREAL), dimension(3) :: euguess
  !
  ! Ntriel for Euler
  !-----------------
  integer(FINT) :: eutrial
  !
  ! Saved cartesian Coordinates: Same convention as above
  !----------------------------------------------------------
  real(FREAL), dimension(3,5), save :: smxyz
  !
  ! Mass center and inertia matrixes
  !---------------------------------
  real(FREAL), dimension(5), save :: atmwe
  real(FREAL), dimension(3), save :: rmass
  real(FREAL), dimension(3,3), save :: rinertia
  !
  ! Bonds and Angles from cartesians
  !---------------------------------
  real(FREAL), dimension(4), save :: xy,rxy
  real(FREAL), dimension(6), save :: yxy,ryxy,cyxy,rcyxy
  !
  ! Metpot vector and outvalue
  !---------------------------
  real(FREAL), dimension(1:1,1:10), save :: vmpotin
  real(FREAL), dimension(1:1), save :: vmpotout
  real(FREAL), save :: maxseenvmpot = 0.0_FREAL
  !
  ! Potential Limit
  !----------------
  real(FREAL), save :: vlimit = -1.0_FREAL
  real(FREAL), dimension(4), save :: cvolume
  real(FREAL), dimension(4), save :: svolume
  !
  ! Done logical flag
  logical, save :: notdone
  !
  ! Internal Coordinates
  ! Definition: dr(i)=X-Yi, da(1:6)=a12,a13,a14,a23,a24,a34
  !                                  1   2   3   4   5   6
  !--------------------------------------------------------
  !real(FREAL), dimension(4), save, target :: dr,vr
  !real(FREAL), dimension(6), save, target :: da,va
  !
  ! Gamma Angles to store for calculation
  !--------------------------------------
  ! Gamma angles definition:
  ! G(1)=a12/3 G( 7)=a12/4
  ! G(2)=a13/2 G( 8)=a13/4
  ! G(3)=a14/2 G( 9)=a14/3
  ! G(4)=a23/1 G(10)=a23/4
  ! G(5)=a24/1 G(11)=a24/3
  ! G(6)=a34/1 G(12)=a34/2
  !real(FREAL), dimension(12), save :: agam
  !real(FREAL), dimension(12), save :: acagam
  !integer(FINT) :: ig
  !real(FREAL) :: fgamma
  !
  ! Polyspherical Coordinates
  ! Definition: plr(i)=X-Y(5-i), pla(1:5)=th3,th2,th1,phi2,phi1
  !                                     1   2   3    4    5   
  !--------------------------------------------------------
  real(FREAL), dimension(4), save, target :: plr,splr
  real(FREAL), dimension(5), save, target :: pla,spla
  !
  ! Random number generator variables
  real(FREAL), dimension(9,2) :: ranrange
  real(FREAL), dimension(9,2) :: sranrange
  integer(FINT), dimension(9,3) :: iranpars
  real(FREAL), dimension(9,3) :: ranpars
  real(FREAL), dimension(9) :: ranvalues
  !integer(FINT), dimension(9) :: iridx
  integer(FINT) :: ix1,ix2,ix3,ix4,ix5,ix6,ix7,ix8,ix9
  !integer(FINT) :: j1,j2,j3,j4,j5,j6,j7,j8,j9
  real(FREAL) :: xran
  !
  ! Symmetry Coordinates
  ! Definition:
  ! Sdr(1)=S1, Sdr(2-4)=S2x,S2y,S2z, Sda(1-2)=S2a,S2b, Sda(3-5)=S4x,S4y,S4z, Sda(6)=Sr
  !-----------------------------------------------------------------------------------
  !real(FREAL), dimension(4), save, target :: Sdr
  !real(FREAL), dimension(6), save, target :: Sda
  !
  ! Rhos
  !-----
  ! Definition: rho1 = da12 + da13 + da14
  !             rho2 = da12 + da23 + da24
  !             rho3 = da13 + da23 + da34
  !             rho4 = da14 + da24 + da34
  !--------------------------------------
  real(FREAL), dimension(4), save :: rho
  !real(FREAL), dimension(4), save :: rho, srho
  !
  ! Solutions for Sda(6)
  !---------------------
  !complex(FREAL), dimension(4) :: ZSda
  !complex(FREAL), dimension(4) :: ZOSda
  ! Dummy complex
  !complex(FREAL) :: Ztmp
  ! index of max real Z solution
  !integer(FINT) :: iZ
  ! Running index of max real Z solution
  !integer(FINT) :: isZ
  !
  ! Files
  !------
  integer(FINT), save :: iuninput
  character(FLCHARS), save :: fileinput
  integer(FINT), save :: iunoutput
  character(FLCHARS), save :: fileoutput
  !integer(FINT), save :: iunascan
  !character(FLCHARS), save :: fileascan
  !integer(FINT), save :: iunscan
  !character(FLCHARS), save :: filescan
  integer(FINT), save :: iunres
  character(FLCHARS), save :: fileres
  integer(FINT), save :: iunxyz
  character(FLCHARS), save :: filexyz
  !
  ! Command Line
  !-------------
  integer(FINT) :: nargs
  character(FLCHARS),dimension(:), allocatable :: cargs
  !
  ! Random Number information
  !--------------------------
  integer(FINT) :: rssize
  integer(FINT), dimension(1) :: rsget
  integer(FINT), dimension(1) :: rsput !=2147483562 !On pgf90 5.0-1
  !
  !
  ! debug and verbose stuff
  !------------------------
  logical, save :: XY4P_debug
  integer(FINT), save :: XY4P_verbose
  integer(FINT), dimension(1), save :: XY4P_rsset
  !
  ! error variables
  !----------------
  !logical, save :: LOCXY4PolySphere_isinit
  integer(FINT), save :: LOCXY4PolySphere_error
  !character(FLCHARS), save :: LOCXY4PolySphere_error_message
  !
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Here we start
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  XY4P_debug = .true.
  XY4P_verbose = 0
  !
  ! Initialization of routines
  !---------------------------
  irc = XY4PolySphere_init()
  if(irc.lt.0) then
     call message(MESERRO,"[XY4P] Cannot initialize.")
     stop 1
  end if
  !
  ! Start reading
  !--------------
  !
  ! Molecule Section which will be our eckard frame
  !------------------------------------------------
  rewind(iuninput)
  !
  !============================================================================
  !============================================================================
  ! SECTION moldecule
  !============================================================================
  !============================================================================
  !
  ! Set the input file to the molecule section...
  irc = osec_set(iuninput,'molecule',params)
  !
  ! Checks if we got it..
  if(irc.lt.0) then
     call message(MESERRO,"[XY4P] No molecule Section.")
     stop 1
  end if
  !
  ! Read Data using module molecule...
  call molecule_read(iuninput,params)
  !
  ! Do some checks...
  irc=molecule_getnat()
  if(irc.ne.5) then
     call message(MESERRO,"[XY4P] This is XY4 program: we need exactly 5 atoms!.")
     stop 1
  end if
  !
  ! Write out read input..
  write(iunoutput,'("[INPUT] # Extract from molecule input section.")')
  if (XY4P_verbose.gt.1) call message(MESOUT,"[XY4P] Eckard reference system readed in.")
  call molecule_print(iunoutput)
  write(iunoutput,'("[END-INPUT] # END extract from molecule input section.")')
  !
  ! Read in xyz from input
  !-----------------------
  mxyz(1:3,1:5)=molecule_xyz(1:3,1:5)
  ! save them...
  smxyz(1:3,1:5)=mxyz(1:3,1:5)
  !
  ! Check if atom 1 (the X or the C in CH4) is on coordinates 0.0 0.0 0.0
  ! NOTE: we only warn here.. we keep going!!!
  if (sum(abs(mxyz(1:3,1))).gt.LEPS) then
     call message(MESWARN,"[XY4P] Atom 1 is not in the origin: Please Check.")
  end if
  !
  ! Test for mass center
  !---------------------
  ! WARNING: We test on the original coordinates which we should use as eckart frame
  !          i.e. the correct way is: Use as eckard frame a molecule we center of
  !          mass on the central atom!
  !---------------------------------------------------------------------------------
  !
  ! collect data and evaluate mass center
  atmwe(1:5)=molecule_zchar(1:5)
  if (XY4P_verbose.gt.1) then
     write(iunoutput,'("#[x-out-comment] Mass charges used for Inertia and Center of Masses.")')
     write(iunoutput,'("#",1X,5(2X,F11.4))') atmwe(1:5)
  end if
  call mathtools_masscenter(molecule_xyz(:,1:5),atmwe,rmass,.false.)
  if (XY4P_verbose.gt.1) then
     write(iunoutput,'("#[x-out-comment] Initial Center of Masses.")')
     write(iunoutput,'("#",1X,3(2X,E14.6))') rmass(1:3)
  end if
  !
  ! If the result is not negligible we do it... and warn..
  if (sum(abs(rmass(1:3))).gt.LEPS) then
     call message(MESWARN,"[XY4P] Molecule not on center of masses. Performing Translation.")
     call mathtools_masscenter(molecule_xyz(:,1:5),atmwe,rmass,.true.)
     call message(MESWARN,"[XY4P] New Coordinates.")
     call molecule_print(iunoutput)
     if (XY4P_verbose.gt.1) then
        call mathtools_masscenter(molecule_xyz(:,1:5),atmwe,rmass,.false.)
        write(iunoutput,'("#[x-out-comment] Last Center of Masses.(CHECK to 0.0)")')
        write(iunoutput,'("#",1X,3(2X,E14.6))') rmass(1:3)
        call molecule_print(iunoutput)
     end if
  end if
  !
  ! Test for Inertia tensor
  !------------------------
  if (XY4P_verbose.gt.1) write(iunoutput,'("#[x-out-comment] Initial Inertia Tensor: Routine comments.")')
  call mathtools_inertia(molecule_xyz(:,1:5),atmwe,rinertia,.false.)
  if (XY4P_verbose.gt.1) then
     write(iunoutput,'("#[x-out-comment] Initial Inertia Tensor eigenvectors.")')
     write(iunoutput,'("#",1X,3(2X,F11.4))') rinertia(1,1:3)
     write(iunoutput,'("#",1X,3(2X,F11.4))') rinertia(2,1:3)
     write(iunoutput,'("#",1X,3(2X,F11.4))') rinertia(3,1:3)
  end if
  !
  if (XY4P_verbose.gt.1) write(iunoutput,'("#[x-out-comment] Inertia Rotate.")')
  call mathtools_inertia(molecule_xyz(:,1:5),atmwe,rinertia,.true.)
  if (XY4P_verbose.gt.1) then
     write(iunoutput,'("#[x-out-comment] Applied Inertia Tensor eigenvectors.")')
     write(iunoutput,'("#",1X,3(2X,F11.4))') rinertia(1,1:3)
     write(iunoutput,'("#",1X,3(2X,F11.4))') rinertia(2,1:3)
     write(iunoutput,'("#",1X,3(2X,F11.4))') rinertia(3,1:3)
  end if

  if (XY4P_verbose.gt.1) call molecule_print(iunoutput)
  if (XY4P_verbose.gt.2) then
     call mathtools_inertia(molecule_xyz(:,1:5),atmwe,rinertia,.false.)
     write(iunoutput,'("#[x-out-comment] Final Inertia Tensor eigenvectors.")')
     write(iunoutput,'("#",1X,3(2X,F11.4))') rinertia(1,1:3)
     write(iunoutput,'("#",1X,3(2X,F11.4))') rinertia(2,1:3)
     write(iunoutput,'("#",1X,3(2X,F11.4))') rinertia(3,1:3)
  end if
  !
  ! This is all about our reference frame...
  !
  !=======================================================================================
  ! The rest need X atom (C in CH4) at origin..
  ! So we do it..
  !=======================================================================================
  !
  mxyz(1,:) = mxyz(1,:) - mxyz(1,1)  
  mxyz(2,:) = mxyz(2,:) - mxyz(2,1)
  mxyz(3,:) = mxyz(3,:) - mxyz(3,1)
  !
  ! Get bonds and angles: remember to skip atom at origin!
  !-------------------------------------------------------
  call getbonds(4,mxyz(1:3,2:5),xy(1:4))
  call getangles(4,mxyz(1:3,2:5),xy(1:4),yxy(1:6))
  !
  ! Save the orinal coordinates bonds and angles
  !---------------------------------------------
  rxy(1:4) = xy(1:4)
  ryxy(1:6) = yxy(1:6)
  cyxy(1:6)=cos(yxy(1:6))
  rcyxy(1:6) = cyxy(1:6)
  !
  ! Evaluate Metpot4 as test
  !-------------------------
  vmpotin(1,1:4) = xy(1:4)
  vmpotin(1,5:10)= cyxy(1:6)
  call metpot4_vint(vmpotin,vmpotout)
  write(iunoutput,'("[x-out-result] Bonds, Angles DEGREE and Metpot4 data:")')
  write(iunoutput,'(1X,11(2X,F11.4))') xy(1:4),yxy(1:6)/LPI*180.0_FREAL,vmpotout(1)
  !
  if (XY4P_verbose.gt.0) then
     write(iunoutput,'("#[x-out-comment] Initial Bonds and Angles (GRAD).")')
     write(iunoutput,'("#",1X,10(2X,A14))') 'r1','r2','r3','r4','a12','a13','a14','a23','a24','a34'
     write(iunoutput,'("#",1X,10(2X,F14.8))') xy(1:4),yxy(1:6)/LPI*180.0_FREAL
  end if
  !
  !============================================================================
  !============================================================================
  ! SECTION x-xy4-polyspherical
  !============================================================================
  !============================================================================
  !
  ! We read the PolySpherical coordinate system
  !--------------------------------------------
  rewind(iuninput)
  iscan = 0
  irc = osec_set(iuninput,'x-xy4-polyspherical',params)
  if(irc.lt.0) then
     call message(MESERRO,"[XY4P] No x-xy4-polyspherical Section.")
     stop 1
  end if
  !
  !
  do while(.TRUE.)
     !
     irc = line_getline(iuninput,line,3)
     !
     ! if EOF exit loop
     if(irc.ne.0) then
        exit
     end if
     !
     ! if new section EXIT
     if(line(1:1).eq.'[') then
        exit
     end if
     iscan = iscan + 1
     read(line,*) plr(1:4),pla(1:5)
     splr(1:4)=plr(1:4)
     pla(1:5)=pla(1:5)/180.0_FREAL*LPI
     spla(1:5)=pla(1:5)
     if (XY4P_verbose.gt.1) then
        call message(MESOUT,"[XY4P] Polyspherical coordinates data readed in.")
        write(iunoutput,'("[x-xy4-polyspherical] INPUT")')
        write(iunoutput,'("#",A11,8(1X,A11))') 'r4','r3','r2','r1','th3','th2','th1','ph2','ph1'
        write(iunoutput,'(9(1X,F11.4))') plr(1:4),pla(1:5)/LPI*180.0_FREAL
        write(iunoutput,'("#",A11,8(1X,A11))') '(RAD) r4','r3','r2','r1','th3','th2','th1','ph2','ph1'
        write(iunoutput,'("#",9(1X,F11.4))') plr(1:4),pla(1:5)
     end if
     !
     !irc=poly2cart(plr(1:4),pla(1:5),mxyz(1:3,2:5),.true.)
     mxyz(:,:) = 0.0_FREAL
     irc=poly2cart(plr(1:4),pla(1:5),mxyz(1:3,2:5))
     if(irc.lt.0) then
        call message(MESERRO,"[XY4P] Cannot perform conversion.")
        stop 3
     end if
     smxyz(1:3,1:5)=mxyz(1:3,1:5)
     !
     ! Print partial result
     !---------------------
     if (XY4P_verbose.gt.1) then
        write(iunoutput,'("[x-out-result] 5 Cartesian From PolySpherical")')
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'X0',mxyz(1:3,1)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y1',mxyz(1:3,2)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y2',mxyz(1:3,3)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y3',mxyz(1:3,4)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y4',mxyz(1:3,5)
     end if
     !
     ! Get bonds and angles: remember to skip atom at origin!
     !-------------------------------------------------------
     call getbonds(4,mxyz(1:3,2:5),xy(1:4))
     call getangles(4,mxyz(1:3,2:5),xy(1:4),yxy(1:6))
     cyxy(1:6)=cos(yxy(1:6))
     !if (XY4P_verbose.gt.0) then
        vmpotin(1,1:4) = xy(1:4)
        vmpotin(1,5:10)= cyxy(1:6)
        !
        call metpot4_vint(vmpotin,vmpotout)
        write(iunoutput,'("[x-out-result] Bonds and Angles DEGREE")')
        write(iunoutput,'("#",1X,11(2X,A11))') 'r1','r2','r3','r4','a12','a13','a14','a23','a24','a34',&
             &'Vmpot4'
        write(iunoutput,'(1X,11(2X,F11.4))') xy(1:4),yxy(1:6)/LPI*180.0_FREAL,vmpotout(1)
        write(iunoutput,'("[x-out-result] Bonds and Angles RAD")')
        write(iunoutput,'("#",1X,10(2X,A11))') 'r1','r2','r3','r4','a12','a13','a14','a23','a24','a34'
        write(iunoutput,'(1X,10(2X,F11.4))') xy(1:4),yxy(1:6)
     !end if
     !
     ! Get Rhos
     !---------
     rho(1)=yxy(1)+yxy(2)+yxy(3)-((3.0_FREAL*109.4712_FREAL)/180.0_FREAL*LPI)
     rho(2)=yxy(1)+yxy(4)+yxy(5)-((3.0_FREAL*109.4712_FREAL)/180.0_FREAL*LPI)
     rho(3)=yxy(2)+yxy(4)+yxy(6)-((3.0_FREAL*109.4712_FREAL)/180.0_FREAL*LPI)
     rho(4)=yxy(3)+yxy(5)+yxy(6)-((3.0_FREAL*109.4712_FREAL)/180.0_FREAL*LPI)
     if (XY4P_verbose.gt.1) then
        write(iunoutput,'("#[x-out-comment] Rhos (RAD).")')
        write(iunoutput,'("#",1X,4(2X,F11.4))' ) rho(1:4)
        write(iunoutput,'("#[x-out-comment]  Sum of Rhos (RAD).",1X,2X,F11.4)') sum(rho(1:4))
     end if
     !
     ! Do center of mass translation and inertia tensor rotation
     !----------------------------------------------------------
     call mathtools_masscenter(mxyz,atmwe(1:5),rmass(1:3),.true.)
     if (XY4P_verbose.gt.1) then
        write(iunoutput,'("#[x-out-comment] PolyInput Translated for Center of Masses:")')
        write(iunoutput,'("#",1X,3(2X,F11.4))') rmass(1:3)
     end if
     call mathtools_inertia(mxyz(:,1:5),atmwe,rinertia,.true.)
     if (XY4P_verbose.gt.1) then
        write(iunoutput,'("#[x-out-comment] PolyInput Rotated for Inertia Tensor.")')
        write(iunoutput,'("#",1X,3(2X,F11.4))') rinertia(1,1:3)
        write(iunoutput,'("#",1X,3(2X,F11.4))') rinertia(2,1:3)
        write(iunoutput,'("#",1X,3(2X,F11.4))') rinertia(3,1:3)
     end if
     !
     if (XY4P_verbose.gt.1) then
        write(iunoutput,'("[x-out-result] 5 Cartesian Before Euler")')
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'X0',mxyz(1:3,1)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y1',mxyz(1:3,2)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y2',mxyz(1:3,3)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y3',mxyz(1:3,4)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y4',mxyz(1:3,5)
     end if
     !
     ! Call euler routine
     !-------------------
     do i=1,5
        vexyz(3*i-2)=molecule_xyz(1,i)
        vexyz(3*i-1)=molecule_xyz(2,i)
        vexyz(3*i-0)=molecule_xyz(3,i)
        !
        vxyz(3*i-2)=mxyz(1,i)
        vxyz(3*i-1)=mxyz(2,i)
        vxyz(3*i-0)=mxyz(3,i)
     end do
     if (XY4P_verbose.gt.1) then
        write(iunoutput,'("[x-comment] Euler Input Reference Frame: ")')
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'X0',molecule_xyz(1:3,1)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y1',molecule_xyz(1:3,2)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y2',molecule_xyz(1:3,3)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y3',molecule_xyz(1:3,4)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y4',molecule_xyz(1:3,5)
        call message(MESOUT,"[XY4P] Running Euler.")
     end if
     write(iunoutput,'("[x-out-euler] ")')
     euguess(1:3)=10.0_FREAL
     eutrial=100_FINT
     call euler(5,15,atmwe(1:5),vxyz,vexyz,euguess(1),euguess(2),euguess(3),eutrial,irc)
     !
     if(irc.lt.0) then
        call message(MESERRO,"[XY4P] Euler faild to find solutions.")
     else
        !
        ! Print partial result
        !---------------------
        do i=1,5
           mxyz(1,i)=vxyz(3*i-2)
           mxyz(2,i)=vxyz(3*i-1)
           mxyz(3,i)=vxyz(3*i-0)
        end do
        if (XY4P_verbose.gt.1) then
           write(iunoutput,'("[x-out-result] 5 Cartesian After Euler")')
           write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'X0',mxyz(1:3,1)
           write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y1',mxyz(1:3,2)
           write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y2',mxyz(1:3,3)
           write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y3',mxyz(1:3,4)
           write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y4',mxyz(1:3,5)
        end if
        write(iunoutput,'("[x-out-result] XYZ format")')
        write(iunoutput,'(1X,I4)') 5
        write(iunoutput,'(1X)')
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'C',mxyz(1:3,1)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'H',mxyz(1:3,2)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'H',mxyz(1:3,3)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'H',mxyz(1:3,4)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'H',mxyz(1:3,5)
        !
        ! At present we do not want initial input data in the res file
        !if (iunres.gt.0) then
        if (.false.) then
           write(iunres,'("# P  SCAN: ",I5)') iscan
           write(iunres,'("# P INPUT: ",9(2X,A18))') 'r1','r2','r3','r4','rho1','rho2','rho3','th1','th2'
           write(iunres,'("# P INPUT: ",10(2X,D18.10))') plr(1:4),pla(1:5),vmpotout(1)
           write(iunres,'("# P   INT: ",10(2X,A18))') 'r1','r2','r3','r4',&
                &                                     'a12','a13','a14','a23','a24','a34'
           write(iunres,'("# P   INT: ",11(2X,D18.10))')  xy(1:4),yxy(1:6)/LPI*180.0_FREAL,vmpotout(1)
           write(iunres,'("# P  CART: ",15(2X,A18))') 'x1','y1','z1','x2','y2','z2','x3','y3','z3',&
                &                                     'x4','y4','z4','x5','y5','z5'
           write(iunres,'("# P  CART: ",16(2X,D18.10))') mxyz(1:3,1:5),vmpotout(1)
           write(iunres,'("# RES FOR RegGrid")')
           write(iunres,'("P#GXC",2X,16(2X,D18.10))') mxyz(1:3,1:5),vmpotout(1)
           write(iunres,'("P#GXP",2X,10(2X,D18.10))') plr(1:4),pla(1:5),vmpotout(1)
           write(iunres,'("P#GXA",2X,10(2X,F18.10))') plr(1:4),pla(1:5)/LPI*180.0_FREAL,vmpotout(1)
        end if
        !
        if (iunxyz.gt.0) then
           write(iunxyz,'("5")')
           !write(iunxyz,'("[x-ingengrid]")')
           write(iunxyz,'(1X,F12.4)') vmpotout(1)
           write(iunxyz,'(1X,A2,1X,3(2X,F12.4))') 'C ',mxyz(1:3,1)
           write(iunxyz,'(1X,A2,1X,3(2X,F12.4))') 'H ',mxyz(1:3,2)
           write(iunxyz,'(1X,A2,1X,3(2X,F12.4))') 'H ',mxyz(1:3,3)
           write(iunxyz,'(1X,A2,1X,3(2X,F12.4))') 'H ',mxyz(1:3,4)
           write(iunxyz,'(1X,A2,1X,3(2X,F12.4))') 'H ',mxyz(1:3,5)
        end if
        !
     end if
  end do
  !
  !
  !============================================================================
  !============================================================================
  ! SECTION x-xy4-genrandom
  !============================================================================
  !============================================================================
  !
  rewind(iuninput)
  irc = osec_set(iuninput,'x-xy4-genrandom',params)
  if(irc.lt.0) then
     call message(MESWARN,"[XY4P] No Random Generation section.")
     call message(MESERRO,"[XY4P] We Use goto now...")
     !I know I shoudn't use it ut at preset here it is
     goto 555
     stop 1
  end if
  !
  ! Reset numer of generated points
  !--------------------------------
  igenpt=0
  !
  write(iunoutput,'(A)') "[x-xy4-genrandom] INPUT"
  do i=1,4
     irc = line_getline(iuninput,line,3)
     read(line,*) ranrange(i,1:2),iranpars(i,1:2)
     write(iunoutput,'(1X,2(1X,F12.4),1X,2(1X,I4))') ranrange(i,1:2),iranpars(i,1:2)
  end do
  !
  do i=5,9
     irc = line_getline(iuninput,line,3)
     read(line,*) ranrange(i,1:2),iranpars(i,1:2)
     write(iunoutput,'(1X,2(1X,F12.4),1X,2(1X,I4))') ranrange(i,1:2),iranpars(i,1:2)
     ranrange(i,1:2)=ranrange(i,1:2)/180.0_FREAL*LPI
  end do
  !
  if (XY4P_verbose.gt.1) then
     call random_seed(size=rssize)
     write(iunoutput,'("#[x-xy4-genrandom] INFO(SIZE):",2X,I20)') rssize
     call random_seed(get=rsget)
     write(iunoutput,'("#[x-xy4-genrandom] INFO(GET):",2X,I20)') rsget
  end if
  if (XY4P_verbose.gt.2) then
     call random_seed(size=rssize)
     write(iunoutput,'("#[x-xy4-genrandom] INFO(SIZE):",2X,I20)') rssize
     call random_seed(get=rsget)
     write(iunoutput,'("#[x-xy4-genrandom] INFO(GET):",2X,I20)') rsget
     !
     rsget=rsput
     !
     call random_seed(put=rsput)
     call random_number(xran)
     write(iunoutput,'("#[x-xy4-genrandom] TEST(0):",2X,F20.12)') xran
     call random_seed(put=rsput)
     call random_number(xran)
     write(iunoutput,'("#[x-xy4-genrandom] TEST(1):",2X,F20.12)') xran
     ! Reset before going on
     call random_seed(put=rsput)
  end if
  if (XY4P_rsset(1).eq.0) then
     call random_seed
  else
     call random_seed(put=XY4P_rsset)
  end if
  !
!  call random_seed(get=rsget)
!  write (iunoutput,'("# RNDSEED",I20)') rsget
  if (iunres.gt.0) then
!     call random_seed(get=rsget)
!     write (iunres,'("# RNDSEED",I20)') rsget
     write (iunres,'("# RNDSEED rs_get not implemented in gfortran")')
     write (iunres,'("# WARNING! Please check order in the code!!!!")')
     write(iunres,'("#",25(2X,A18))') 'r1','r2','r3','r4','rho1','rho2','rho3','th1','th2',&
          & 'x1','y1','z1','x2','y2','z2','x3','y3','z3','x4','y4','z4','x5','y5','z5',&
          & 'Vmetpot'!plr(1:4),pla(1:5),mxyz(1:3,2:5)
  end if
  !
  do ix1=1,iranpars(1,1)
     call random_number(xran)
     ranvalues(1)=xran*(ranrange(1,2)-ranrange(1,1))+ranrange(1,1)
  do ix2=1,iranpars(2,1)                                                     
     call random_number(xran)
     ranvalues(2)=xran*(ranrange(2,2)-ranrange(2,1))+ranrange(2,1)
  do ix3=1,iranpars(3,1)                                                     
     call random_number(xran)
     ranvalues(3)=xran*(ranrange(3,2)-ranrange(3,1))+ranrange(3,1)
  do ix4=1,iranpars(4,1)                                                     
     call random_number(xran)
     ranvalues(4)=xran*(ranrange(4,2)-ranrange(4,1))+ranrange(4,1)
  do ix5=1,iranpars(5,1)                                                     
     call random_number(xran)
     ranvalues(5)=xran*(ranrange(5,2)-ranrange(5,1))+ranrange(5,1)
  do ix6=1,iranpars(6,1)                                                     
     call random_number(xran)
     ranvalues(6)=xran*(ranrange(6,2)-ranrange(6,1))+ranrange(6,1)
  do ix7=1,iranpars(7,1)                                                     
     call random_number(xran)
     ranvalues(7)=xran*(ranrange(7,2)-ranrange(7,1))+ranrange(7,1)
  do ix8=1,iranpars(8,1)                                                     
     call random_number(xran)
     ranvalues(8)=xran*(ranrange(8,2)-ranrange(8,1))+ranrange(8,1)
  do ix9=1,iranpars(9,1)                                                     
     call random_number(xran)
     ranvalues(9)=xran*(ranrange(9,2)-ranrange(9,1))+ranrange(9,1)
     if (XY4P_verbose.gt.1) write(iunoutput,'("#",1X,A10,9(1X,F12.4))') 'RndPoint: ', ranvalues(:)
     !
     plr(1)=ranvalues(1)
     plr(2)=ranvalues(2)
     plr(3)=ranvalues(3)
     plr(4)=ranvalues(4)
     pla(1)=ranvalues(5)
     pla(2)=ranvalues(6)
     pla(3)=ranvalues(7)
     pla(4)=ranvalues(8)
     pla(5)=ranvalues(9)
     !irc=poly2cart(plr(1:4),pla(1:5),mxyz(1:3,2:5),.true.)
     mxyz(:,:) = 0.0_FREAL
     irc=poly2cart(plr(1:4),pla(1:5),mxyz(1:3,2:5))
     if(irc.lt.0) then
        call message(MESERRO,"[XY4P] Error in POLY2CART.")
     end if
     call getbonds(4,mxyz(1:3,2:5),xy(1:4))
     call getangles(4,mxyz(1:3,2:5),xy(1:4),yxy(1:6))
     cyxy(1:6)=cos(yxy(1:6))
     !
     ! Efore to write results we increese the numer of generated points
     !-----------------------------------------------------------------
     igenpt = igenpt + 1
     !
     write(iunoutput,'("[x-ran-result] Bonds and Angles DEGREE")')
     write(iunoutput,'(2X,10(2X,F11.4))') xy(1:4),yxy(1:6)/LPI*180.0_FREAL
     call mathtools_masscenter(mxyz,atmwe(1:5),rmass(1:3),.true.)
     euguess(1:3)=10.0_FREAL
     eutrial=100_FINT
     call euler(5,15,atmwe(1:5),vxyz,vexyz,euguess(1),euguess(2),euguess(3),eutrial,irc)
     if(irc.lt.0) then
        call message(MESERRO,"[XY4P] Euler faild to find solutions.")
     else
        write(iunoutput,'("[x-out-result] 5 Cartesian After Euler")')
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'X0',mxyz(1:3,1)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y1',mxyz(1:3,2)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y2',mxyz(1:3,3)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y3',mxyz(1:3,4)
        write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y4',mxyz(1:3,5)
     end if
     !
     write(iunoutput,'("[x-out-result] Poly Random in Degrees:")')
     write(iunoutput,'(1X,A2,1X,9(2X,F10.4))') 'PR',ranvalues(1:4),ranvalues(5:9)/LPI*180.0_FREAL
     vmpotin(1,1:4) = xy(1:4)
     vmpotin(1,5:10)= cyxy(1:6)
     !
     call metpot4_vint(vmpotin,vmpotout)
     if (vlimit.gt.LEPS) then
        if (vmpotout(1).lt.vlimit) then
           write(iunoutput,'("[x-out-result] METPOT4 potential")')
           write(iunoutput,'(1X,A2,1X,10(2X,F10.4))') 'VR',vmpotin(1,1:10)
           write(iunoutput,'(1X,A2,1X,1(2X,F22.12))') 'VI',vmpotout(1)
           if (iunres.gt.0) then
              write(iunres,'(1X,25(2X,D18.10))') plr(1:4),pla(1:5),mxyz(1:3,1:5),vmpotout(1)
           end if
        else
           write(iunoutput,'("[x-out-result] METPOT4 potential OUT OF LIMIT")')
           write(iunoutput,'(1X,A2,1X,10(2X,F10.4))') 'VR',vmpotin(1,1:10)
           write(iunoutput,'(1X,A2,1X,1(2X,F12.4))') 'VI',vmpotout(1)
           if (iunres.gt.0) then
              write(iunres,'(1X,24(2X,D18.10),1X,A3,2X,D18.10)') plr(1:4),pla(1:5),mxyz(1:3,1:5),'OUT',vmpotout(1)
           end if
        end if
     else
        if (iunres.gt.0) then
           write(iunres,'(1X,25(2X,D18.10))') plr(1:4),pla(1:5),mxyz(1:3,1:5),vmpotout(1)
           write(iunres,'(1X,"#RXC",2X,16(2X,D18.10))') mxyz(1:3,1:5),vmpotout(1)
           write(iunres,'(1X,"#RXP",2X,10(2X,D18.10))') plr(1:4),pla(1:5),vmpotout(1)
           write(iunres,'(1X,"#RXA",2X,10(2X,F18.10))') plr(1:4),pla(1:5)/LPI*180.0_FREAL,vmpotout(1)
        end if
        write(iunoutput,'("[x-out-result] METPOT4 potential")')
        write(iunoutput,'(1X,A2,1X,10(2X,F10.4))') 'VR',vmpotin(1,1:10)
        write(iunoutput,'(1X,A2,1X,1(2X,F12.4))') 'VI',vmpotout(1)
     end if
     !
  end do
  end do
  end do
  end do
  end do
  end do
  end do
  end do
  end do
  !
  !
  !============================================================================
  !============================================================================
  ! SECTION x-xy4-gengrid
  !============================================================================
  !============================================================================
  !
  ! Yep Again ... here we come from section random...
  555 continue
  rewind(iuninput)
  irc = osec_set(iuninput,'x-xy4-gengrid',params)
  if(irc.lt.0) then
     call message(MESERRO,"[XY4P] No GenGrid Generation section.")
     ! Here we do not use goto ecause we are the last at present...
     ! This part should be a it rewritten.. 
     goto 777
     stop 1
  end if
  !
  ! Reset Number of generated points
  !---------------------------------
  igenpt = 0
  !
  ! READ input and reproduce it as output for this INPUT section
  !-------------------------------------------------------------
  !-------------------------------------------------------------
  write(iunoutput,'(A)') "[x-xy4-gengrid] INPUT"
  !
  !Alternative code...
  do i=1,9
     irc = line_getline(iuninput,line,3)
     read(line,*) ranrange(i,1:2),iranpars(i,1)
     write(iunoutput,'(1X,2(1X,F12.4),1X,1(1X,I4))') ranrange(i,1:2),iranpars(i,1)
  end do
  !
  ! Convert angles in radiants
  do i=5,9
     ranrange(i,1:2)=ranrange(i,1:2)/180.0_FREAL*LPI
  end do
  !
  !End Alternative code...
  !
  ! Sustitute:
  !!do i=1,4
  !!   irc = line_getline(iuninput,line,3)
  !!   read(line,*) ranrange(i,1:2),iranpars(i,1)
  !!   write(iunoutput,'(1X,2(1X,F12.4),1X,1(1X,I4))') ranrange(i,1:2),iranpars(i,1)
  !!end do
  !!!
  !!do i=5,9
  !!   irc = line_getline(iuninput,line,3)
  !!   read(line,*) ranrange(i,1:2),iranpars(i,1)
  !!   write(iunoutput,'(1X,2(1X,F12.4),1X,1(1X,I4))') ranrange(i,1:2),iranpars(i,1)
  !!   ranrange(i,1:2)=ranrange(i,1:2)/180.0_FREAL*LPI
  !!end do
  !!
  ! End of Sustitute
  !!
  !
  write(iunoutput,'(A)') "[x-xy4-gengrid] INFO"
  write(iunoutput,'(A)') "[x-xy4-gengrid-info]"
  write(iunoutput,'(A)') "#We collected 9 hardcode coordinate"
  write(iunoutput,'(A)') "#1-4 are considered bonds"
  write(iunoutput,'(A)') "#5-9 are considered angles"
  write(iunoutput,'(A)') "#5-9 are directly conveted to radians for internal use."
  write(iunoutput,'(A)') "#The numer of points for each coords is (in order)"
  write(iunoutput,'(1X,9(1X,I3))') iranpars(:,1)
  write(iunoutput,'(A)') "#The total numer of points will be:"
  write(iunoutput,'(1X,I10)') product(iranpars(:,1))


  !
!  if (iunres.gt.0) then
!     write(iunres,'("# GENGRID running")')
!     write(iunres,'("# WARNING! Please check order in the code!!!!")')
!     write(iunres,'("#",25(2X,A18))') 'r1','r2','r3','r4','rho1','rho2','rho3','th1','th2','Vmetpot'
!     write(iunres,'("#",25(2X,A18))') 'x1','y1','z1','x2','y2','z2','x3','y3','z3','x4','y4','z4','x5','y5','z5','Vmetpot'
!     write(iunres,'("# WARNING! Most probaly is!!!!")')
!     write(iunres,'("#",25(2X,A18))') 'r4','r3','r2','r1','rho3','rho2','rho1','th2','th1','Vmetpot'
!     write(iunres,'("#",25(2X,A18))') 'x1','y1','z1','x2','y2','z2','x3','y3','z3','x4','y4','z4','x5','y5','z5','Vmetpot'
!  end if
  !
  write(iunoutput,'(A)') "[x-xy4-gengrid-info] Start of Loops..."
  write(iunoutput,'(A)') "# Results are given after:[x-gengrid-result-new-point] Point Numer:"
  write(iunoutput,'(A)') "# P is polyspherical"
  write(iunoutput,'(A)') "# P coordinates scheme is:"
  write(iunoutput,'("#",A1,1X,9(2X,A11))') 'P','r4','r3','r2','r1','t3','t2','t1','p2','p1'
  write(iunoutput,'(A)') "# I is internal"
  write(iunoutput,'(A)') "# I coordinate scheme is:"
  write(iunoutput,'("#",A1,1X,10(2X,A11))') 'I','r1','r2','r3','r4','a12','a13','a14','a23','a24','a34'
  write(iunoutput,'(A)') "# V is cosine internals with metpot results"
  write(iunoutput,'(A)') "# V scheme is cosines of I plus METPOT4 energy"
  write(iunoutput,'(A)') "# VI also repeat the current V value."
  write(iunoutput,'(A)') "# "
  write(iunoutput,'(A)') "# "
  !
  xran=1.0
  !
  ! Here we have a seriese of do loops with inner ifs.. I am sure that this can
  ! be done better for speed but at present we stay with a kinda clean code
  !
  do ix1=1,iranpars(1,1)
     if (iranpars(1,1).eq.1) then
        ranvalues(1)=ranrange(1,1)
     else
        ranvalues(1)=((ranrange(1,2)-ranrange(1,1))/(iranpars(1,1)-1))*(ix1-1)+ranrange(1,1)
     end if
  do ix2=1,iranpars(2,1)                                                     
     if (iranpars(2,1).eq.1) then
        ranvalues(2)=ranrange(2,1)
     else
        ranvalues(2)=((ranrange(2,2)-ranrange(2,1))/(iranpars(2,1)-1))*(ix2-1)+ranrange(2,1)
     end if
  do ix3=1,iranpars(3,1)                                                     
     if (iranpars(3,1).eq.1) then
        ranvalues(3)=ranrange(3,1)
     else
        ranvalues(3)=((ranrange(3,2)-ranrange(3,1))/(iranpars(3,1)-1))*(ix3-1)+ranrange(3,1)
     end if
  do ix4=1,iranpars(4,1)                                                     
     if (iranpars(4,1).eq.1) then
        ranvalues(4)=ranrange(4,1)
     else
        ranvalues(4)=((ranrange(4,2)-ranrange(4,1))/(iranpars(4,1)-1))*(ix4-1)+ranrange(4,1)
     end if
  do ix5=1,iranpars(5,1)                                                     
     if (iranpars(5,1).eq.1) then
        ranvalues(5)=ranrange(5,1)
     else
        ranvalues(5)=((ranrange(5,2)-ranrange(5,1))/(iranpars(5,1)-1))*(ix5-1)+ranrange(5,1)
     end if
  do ix6=1,iranpars(6,1)                                                     
     if (iranpars(6,1).eq.1) then
        ranvalues(6)=ranrange(6,1)
     else
        ranvalues(6)=((ranrange(6,2)-ranrange(6,1))/(iranpars(6,1)-1))*(ix6-1)+ranrange(6,1)
     end if
  do ix7=1,iranpars(7,1)                                                     
     if (iranpars(7,1).eq.1) then
        ranvalues(7)=ranrange(7,1)
     else
        ranvalues(7)=((ranrange(7,2)-ranrange(7,1))/(iranpars(7,1)-1))*(ix7-1)+ranrange(7,1)
     end if
  do ix8=1,iranpars(8,1)                                                     
     if (iranpars(8,1).eq.1) then
        ranvalues(8)=ranrange(8,1)
     else
        ranvalues(8)=((ranrange(8,2)-ranrange(8,1))/(iranpars(8,1)-1))*(ix8-1)+ranrange(8,1)
     end if
  do ix9=1,iranpars(9,1)                                                     
     if (iranpars(9,1).eq.1) then
        ranvalues(9)=ranrange(9,1)
     else
        ranvalues(9)=((ranrange(9,2)-ranrange(9,1))/(iranpars(9,1)-1))*(ix9-1)+ranrange(9,1)
     end if
     !if (XY4P_verbose.gt.1) write(iunoutput,'("#",1X,A10,9(1X,F12.4))') 'RndPoint: ', ranvalues(:)
     !
     !
     ! before to write the outpu we increase the generated data counter
     !-----------------------------------------------------------------
     igenpt = igenpt + 1
     !
     ! We also inform the out that this is a new point
     !------------------------------------------------
     write(iunoutput,'("[x-gengrid-result-new-point] Point Numer:",1X,I10)') igenpt
     !
     ! Write read in coordinates
     !write(iunoutput,'("[x-gengrid-result] Poly GenGrid in Degrees:")')
     !write(iunoutput,'("#",A2,1X,9(2X,A11))') 'PR','r4','r3','r2','r1','t3','t2','t1','p2','p1'
     write(iunoutput,'(1X,A1,1X,9(2X,F11.4))') 'P',ranvalues(1:4),ranvalues(5:9)/LPI*180.0_FREAL
     plr(1)=ranvalues(1)
     plr(2)=ranvalues(2)
     plr(3)=ranvalues(3)
     plr(4)=ranvalues(4)
     !
     pla(1)=ranvalues(5)
     pla(2)=ranvalues(6)
     pla(3)=ranvalues(7)
     pla(4)=ranvalues(8)
     pla(5)=ranvalues(9)
     !irc=poly2cart(plr(1:4),pla(1:5),mxyz(1:3,2:5),.true.)
     mxyz(:,:) = 0.0_FREAL
     irc=poly2cart(plr(1:4),pla(1:5),mxyz(1:3,2:5))
     if(irc.lt.0) then
        call message(MESERRO,"[XY4P] Error in POLY2CART.")
     end if
     call getbonds(4,mxyz(1:3,2:5),xy(1:4))
     call getangles(4,mxyz(1:3,2:5),xy(1:4),yxy(1:6))
     cyxy(1:6)=cos(yxy(1:6))
     !
     !write(iunoutput,'("[x-gengrid-result] Bonds and Angles DEGREE")')
     !write(iunoutput,'("#",3X,10(2X,A11))') 'r1','r2','r3','r4','a12','a13','a14','a23','a24','a34'
     write(iunoutput,'(1X,A1,1X,10(2X,F11.4))') 'I',xy(1:4),yxy(1:6)/LPI*180.0_FREAL
     !write(iunoutput,'("[x-gengrid-result] Poly GenGrid in Degrees:")')
     !write(iunoutput,'(1X,A2,1X,9(2X,F11.4))') 'PR',ranvalues(1:4),ranvalues(5:9)/LPI*180.0_FREAL

        !write(iunoutput,'("[x-out-result] 5 Cartesian asis")')
        !write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'X0',mxyz(1:3,1)
        !write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y1',mxyz(1:3,2)
        !write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y2',mxyz(1:3,3)
        !write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y3',mxyz(1:3,4)
        !write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y4',mxyz(1:3,5)

     call mathtools_masscenter(mxyz,atmwe(1:5),rmass(1:3),.true.)
     euguess(1:3)=10.0_FREAL
     eutrial=100_FINT
     call euler(5,15,atmwe(1:5),vxyz,vexyz,euguess(1),euguess(2),euguess(3),eutrial,irc)
     if(irc.lt.0) then
        call message(MESERRO,"[XY4P] Euler faild to find solutions.")
        stop 9
     else
        !write(iunoutput,'("[x-out-result] 5 Cartesian After Euler")')
        !write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'X0',mxyz(1:3,1)
        !write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y1',mxyz(1:3,2)
        !write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y2',mxyz(1:3,3)
        !write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y3',mxyz(1:3,4)
        !write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y4',mxyz(1:3,5)
     end if
     !
     !write(iunoutput,'("[x-out-result] Poly GenGrid in Degrees:")')
     !write(iunoutput,'(1X,A2,1X,9(2X,F10.4))') 'PR',ranvalues(1:4),ranvalues(5:9)/LPI*180.0_FREAL
     vmpotin(1,1:4) = xy(1:4)
     vmpotin(1,5:10)= cyxy(1:6)
     !
     call metpot4_vint(vmpotin,vmpotout)
     if(vmpotout(1).gt.maxseenvmpot) maxseenvmpot=vmpotout(1)
     !
     if (iunxyz.gt.0) then
        write(iunxyz,'("5")')
        !write(iunxyz,'("[x-ingengrid]")')
        write(iunxyz,'(1X,F12.4)') vmpotout(1)
        write(iunxyz,'(1X,A2,1X,3(2X,F12.4))') 'C ',mxyz(1:3,1)
        write(iunxyz,'(1X,A2,1X,3(2X,F12.4))') 'H ',mxyz(1:3,2)
        write(iunxyz,'(1X,A2,1X,3(2X,F12.4))') 'H ',mxyz(1:3,3)
        write(iunxyz,'(1X,A2,1X,3(2X,F12.4))') 'H ',mxyz(1:3,4)
        write(iunxyz,'(1X,A2,1X,3(2X,F12.4))') 'H ',mxyz(1:3,5)
     end if
     !
     if (vlimit.gt.LEPS) then
        if (vmpotout(1).lt.vlimit) then
           !write(iunoutput,'("[x-out-result] METPOT4 potential")')
           write(iunoutput,'(1X,A1,1X,10(2X,F11.4),1X,1(2X,F22.12))') 'V',vmpotin(1,1:10),vmpotout(1)
           write(iunoutput,'(1X,A2,1X,1(2X,F22.12))') 'VI',vmpotout(1)
           if (iunres.gt.0) then
              !write(iunres,'(1X,25(2X,D18.10))') plr(1:4),pla(1:5),mxyz(1:3,1:5),vmpotout(1)
              !
              ! Cartesina output
              write(iunres,'("#GXC",2X,16(2X,D18.10))') mxyz(1:3,1:5),vmpotout(1)
              ! polyspherical output in radians
              write(iunres,'("#GXP",2X,10(2X,D18.10))') plr(1:4),pla(1:5),vmpotout(1)
              ! polyspherical output in degrees
              write(iunres,'("#GXA",2X,10(2X,F18.10))') plr(1:4),pla(1:5)/LPI*180.0_FREAL,vmpotout(1)
              ! Internal output in radians
              write(iunres,'("#RFIR",4(1X,F8.4),1X,6(1X,F8.6),F14.4)') xy(1:4),yxy(1:6),vmpotout(1)
              ! Internal output in degrees
              write(iunres,'("#RFIA",4(1X,F8.4),1X,6(1X,F8.4),F14.4)') xy(1:4),yxy(1:6)/LPI*180.0_FREAL,vmpotout(1)
              write(iunres,'("#RFPA",4(1X,F8.4),1X,5(1X,F8.4),9X,F14.4)') plr(1:4),pla(1:5)/LPI*180.0_FREAL,vmpotout(1)
           end if
        else
           !On out of limits we print it as a negative number.
           !write(iunoutput,'("[x-out-result] METPOT4 potential OUT OF LIMIT")')
           write(iunoutput,'(1X,A1,1X,10(2X,F11.4),1X,1(2X,F22.12))') 'V',vmpotin(1,1:10),-vmpotout(1)
           !write(iunoutput,'(1X,A2,1X,10(2X,F11.4))') 'VR',vmpotin(1,1:10)
           write(iunoutput,'(1X,A2,1X,1(2X,F12.4))') 'VI',-vmpotout(1)
           if (iunres.gt.0) then
              !write(iunres,'(1X,24(2X,D18.10),1X,A3,2X,D18.10)') plr(1:4),pla(1:5),mxyz(1:3,1:5),'OUT',vmpotout(1)
              write(iunres,'("#GXC",2X,16(2X,D18.10))') mxyz(1:3,1:5),-vmpotout(1)
              write(iunres,'("#GXP",2X,10(2X,D18.10))') plr(1:4),pla(1:5),-vmpotout(1)
              write(iunres,'("#GXA",2X,10(2X,F18.10))') plr(1:4),pla(1:5)/LPI*180.0_FREAL,-vmpotout(1)
              write(iunres,'("#RFIR",4(1X,F8.4),1X,6(1X,F8.6),F14.4)') xy(1:4),yxy(1:6),-vmpotout(1)
              write(iunres,'("#RFPR",4(1X,F8.4),1X,5(1X,F8.6),9X,F14.4)') plr(1:4),pla(1:5),-vmpotout(1)
              write(iunres,'("#RFIA",4(1X,F8.4),1X,6(1X,F8.4),F14.4)') xy(1:4),yxy(1:6)/LPI*180.0_FREAL,-vmpotout(1)
              write(iunres,'("#RFPA",4(1X,F8.4),1X,5(1X,F8.4),9X,F14.4)') plr(1:4),pla(1:5)/LPI*180.0_FREAL,-vmpotout(1)
           end if
        end if
     else
        if (iunres.gt.0) then
           !write(iunres,'(1X,25(2X,D18.10))') plr(1:4),pla(1:5),mxyz(1:3,1:5),vmpotout(1)
           write(iunres,'("#GXC",2X,16(2X,D18.10))') mxyz(1:3,1:5),vmpotout(1)
           write(iunres,'("#GXP",2X,10(2X,D18.10))') plr(1:4),pla(1:5),vmpotout(1)
           write(iunres,'("#GXA",2X,10(2X,F18.10))') plr(1:4),pla(1:5)/LPI*180.0_FREAL,vmpotout(1)
           write(iunres,'("#GXI",2X,11(2X,F18.10))') xy(1:4),yxy(1:6)/LPI*180.0_FREAL,vmpotout(1)
           write(iunres,'("#RFIR",4(1X,F10.6),1X,6(1X,F10.8),F14.4)') xy(1:4),yxy(1:6),vmpotout(1)
           write(iunres,'("#RFPR",4(1X,F10.6),1X,5(1X,F10.8),11X,F14.4)') plr(1:4),pla(1:5),vmpotout(1)
           write(iunres,'("#RFIA",4(1X,F10.6),1X,6(1X,F10.6),F14.4)') xy(1:4),yxy(1:6)/LPI*180.0_FREAL,vmpotout(1)
           write(iunres,'("#RFPA",4(1X,F10.6),1X,5(1X,F10.6),11X,F14.4)') plr(1:4),pla(1:5)/LPI*180.0_FREAL,vmpotout(1)
        end if
        !write(iunoutput,'("[x-out-result] METPOT4 potential")')
        write(iunoutput,'(1X,A1,1X,10(2X,F11.4),1X,1(2X,F22.12))') 'V',vmpotin(1,1:10),vmpotout(1)
        !write(iunoutput,'(1X,A2,1X,10(2X,F10.4))') 'VR',vmpotin(1,1:10)
        write(iunoutput,'(1X,A2,1X,1(2X,F22.12))') 'VI',vmpotout(1)
     end if
     !
  end do
  end do
  end do
  end do
  end do
  end do
  end do
  end do
  end do
  !
  write(iunoutput,'("[x-out-result] Statistics")')
  write(iunoutput,'("Largest METPOT seen so far:",1X,F12.4)') maxseenvmpot
  !
  !


  !
  !
  !============================================================================
  !============================================================================
  ! SECTION x-xy4-srcgrid
  !============================================================================
  !============================================================================
  !
  ! Yep Again ... here we come from section random...
  777 continue
  rewind(iuninput)
  irc = osec_set(iuninput,'x-xy4-srcgrid',params)
  if(irc.lt.0) then
     call message(MESERRO,"[XY4P] No SearchGrid Generation section.")
     ! Here we do not use goto ecause we are the last at present...
     ! This part should be a it rewritten.. 
     goto 999
     stop 1
  end if
  !
  ! Reset Number of generated points
  !---------------------------------
  igenpt = 0
  !
  ! READ input and reproduce it as output for this INPUT section
  !-------------------------------------------------------------
  !-------------------------------------------------------------
  !
  write(iunoutput,'(A)') "[x-xy4-srcgrid] INPUT"
  !
  ! Read the vlim potential
  irc = line_getline(iuninput,line,3)
  read(line,*) vlimit
  if(vlimit.lt.LEPS) then
     write(iunoutput,'(A)') "# Error! Input vlimit too small!"
     stop 45
  end if
  !
  write(iunoutput,'(1X,F12.4)') vlimit
  !
  !Alternative code...
  do i=1,9
     irc = line_getline(iuninput,line,3)
     read(line,*) ranrange(i,1:2),ranpars(i,1)
     write(iunoutput,'(1X,2(1X,F12.4),1X,1(1X,F12.4))') ranrange(i,1:2),ranpars(i,1)
  end do
  !Save Ranges
  sranrange(1:9,1:2) = ranrange(1:9,1:2)
! We do not convert it yet!
!  !
!  ! Convert angles in radiants
!  do i=5,9
!     ranrange(i,1:2)=ranrange(i,1:2)/180.0_FREAL*LPI
!  end do
  !
  !End Alternative code...
  !
  ! Sustitute:
  !!do i=1,4
  !!   irc = line_getline(iuninput,line,3)
  !!   read(line,*) ranrange(i,1:2),iranpars(i,1)
  !!   write(iunoutput,'(1X,2(1X,F12.4),1X,1(1X,I4))') ranrange(i,1:2),iranpars(i,1)
  !!end do
  !!!
  !!do i=5,9
  !!   irc = line_getline(iuninput,line,3)
  !!   read(line,*) ranrange(i,1:2),iranpars(i,1)
  !!   write(iunoutput,'(1X,2(1X,F12.4),1X,1(1X,I4))') ranrange(i,1:2),iranpars(i,1)
  !!   ranrange(i,1:2)=ranrange(i,1:2)/180.0_FREAL*LPI
  !!end do
  !!
  ! End of Sustitute
  !!
  !
  write(iunoutput,'(A)') "[x-xy4-srcgrid] INFO"
  write(iunoutput,'(A)') "[x-xy4-srcgrid-info]"
  write(iunoutput,'(A)') "#We collected 9 hardcode coordinate"
  write(iunoutput,'(A)') "#1-4 are considered bonds"
  write(iunoutput,'(A)') "#5-9 are considered angles"
  write(iunoutput,'(A)') "#5-9 will be conveted to radians for internal use."
  write(iunoutput,'(A)') "#The step size for each coords is (in order: au, grad)"
  write(iunoutput,'(1X,9(1X,F12.4))') ranpars(:,1)


!  !
!  if (iunres.gt.0) then
!     write (iunres,'("# GENGRID running")')
!     write (iunres,'("# WARNING! Please check order in the code!!!!")')
!     write(iunres,'("#",25(2X,A18))') 'r1','r2','r3','r4','rho1','rho2','rho3','th1','th2','Vmetpot'
!     write(iunres,'("#",25(2X,A18))') 'x1','y1','z1','x2','y2','z2','x3','y3','z3','x4','y4','z4','x5','y5','z5','Vmetpot'
!     write (iunres,'("# WARNING! Most probaly is!!!!")')
!     write(iunres,'("#",25(2X,A18))') 'r4','r3','r2','r1','rho3','rho2','rho1','th2','th1','Vmetpot'
!     write(iunres,'("#",25(2X,A18))') 'x1','y1','z1','x2','y2','z2','x3','y3','z3','x4','y4','z4','x5','y5','z5','Vmetpot'
!  end if
!  !
!  write(iunoutput,'(A)') "[x-xy4-srcgrid-info] Start of Loops..."
!  !
  xran=1.0
  !
  !================
  ! SRC algorithms
  !================
  !
  ! REMOVE THIS PART: NOT CONSISTENT WITH THE CODE AND
  ! ALSO NOT MUCH USE OF IT!!!
  !
!  ! Read starting coordinates: same as molecule section
!  !----------------------------------------------------
!  ! THIS IS REDUNDANT CODE AT PRESENT!
!  !
!  xy(1:4) = rxy(1:4)
!  cyxy(1:6) = rcyxy(1:6)
!  !
!  ! The rxy and rcyxy are coming from the [molecule] section right on top!
!  ! Yep Technically speaking we need to give some order to this code.
!  ! Right now we check that the [molecule] section coordinates are not off...
!  !
!  ! Check starting point in internals
!  !----------------------------------
!  vmpotin(1,1:4) = xy(1:4)
!  vmpotin(1,5:10)= cyxy(1:6)
!  call metpot4_vint(vmpotin,vmpotout)
!  !
!  if(vmpotout(1).lt.LEPS) then
!     write(iunoutput,'(A)') "# Error! Input start coord give a too small vlimit vs LEPS!"
!     stop 46
!  end if  
!  !
!  if(vmpotout(1).gt.vlimit) then
!     write(iunoutput,'(A)') "# Error! Input start coord give a too small vlimit vs start coord potential"
!     stop 47
!  end if
!  !
!

  !
  ! Check starting point in Poluspherical
  !--------------------------------------
  !
  ! determine our best guess for Polyspherical coordinates
  ! THe starting point is now the medium geometry! Out of the given
  ! grid limits.
  ! It is given
  ! in the next two lines: i.e. medium euclidean point.
  ! Note convertion to radian is done here!
  splr(1:4) = ((ranrange(1:4,2)-ranrange(1:4,1))/2.0_FREAL)+ranrange(1:4,1)
  spla(1:5) = ((ranrange(5:9,2)-ranrange(5:9,1))/2.0_FREAL)+ranrange(5:9,1)
  if( XY4P_debug) then
     write(iunoutput,'("[x-xy4-polyspherical] Debug info central starting point")')
     write(iunoutput,'("#",9(1X,F11.4))') splr(1:4),spla(1:5)
  end if
  ! we get the metpot here: pla goes to radians
  spla(1:5)=spla(1:5)/180.0_FREAL*LPI
  mxyz(:,:) = 0.0_FREAL
  irc=poly2cart(splr(1:4),spla(1:5),mxyz(1:3,2:5))
  if(irc.lt.0) then
     call message(MESERRO,"[XY4P] Cannot perform conversion.")
     stop 3
  end if
  call getbonds(4,mxyz(1:3,2:5),xy(1:4))
  call getangles(4,mxyz(1:3,2:5),xy(1:4),yxy(1:6))
  cyxy(1:6)=cos(yxy(1:6))
  vmpotin(1,1:4) = xy(1:4)
  vmpotin(1,5:10)= cyxy(1:6)
  call metpot4_vint(vmpotin,vmpotout)
  if( XY4P_debug) then
     write(iunoutput,'("[x-xy4-polyspherical] Debug info central starting point (RAD)")')
     write(iunoutput,'("#",10(1X,F11.4))') splr(1:4),spla(1:5),vmpotout(1)
  end if

  if(vmpotout(1).lt.LEPS) then
     write(iunoutput,'(A)') "# Error! Input start coord give a too small vlimit vs LEPS!"
     stop 56
  end if
  !
  if(vmpotout(1).gt.vlimit) then
     write(iunoutput,'(A)') "# Error! Input start coord give a too small vlimit vs start coord potential"
     write(iunoutput,*) vmpotout(1), vlimit
     stop 57
  end if
 
  !
  ! Set starting limits
  !--------------------
  ! i.e. Given the grid limits we set a starting hypercubes from the
  ! guessed medium point as +- the step size.
  !
  ! Here we use spl[ra].. to be consistent we should use some ranrenge value..
  ! i.e. tha same as it is done to evaluate spl[ra] just above...
  write(iunoutput,'(A)') "[x-xy4-srcgrid-info] Start limits"
  do i=1,4
     ! The radial part should kinda stay as is
     ranrange(i,1)=splr(i) - ranpars(i,1)
     ranrange(i,2)=splr(i) + ranpars(i,1)
     !print *,ranrange(i,1),splr(i) , ranpars(i,1)
     ! this should set it as not done....
     ranpars(i,2) = 0.0_FREAL
     write(iunoutput,'(1X,2(1X,F12.4),1X,F12.4)') ranrange(i,1), ranrange(i,2), ranpars(i,1)
  end do
!
! FXT:ME:THere is something odd in here! Check! TODO!
  do i=5,9
     ranrange(i,1)=spla(i-4)*180.0_FREAL/LPI - ranpars(i,1)
     ranrange(i,2)=spla(i-4)*180.0_FREAL/LPI + ranpars(i,1)
     ! this should set it as not done....
     ranpars(i,2) = 0.0_FREAL
     write(iunoutput,'(1X,2(1X,F12.4),1X,F12.4)') ranrange(i,1), ranrange(i,2), ranpars(i,1)
  end do
  ! convert it back to radians
  ranrange(5:9,1:2) = ranrange(5:9,1:2)/180.0_FREAL*LPI
  ! might invert the problem and give convertion of the print print out...
  !print *, ranrange(1:4,1),ranrange(5:9,1)*180.0_FREAL/LPI
  !print *, ranrange(1:4,2),ranrange(5:9,2)*180.0_FREAL/LPI
  !
  ! Check starting limits: The initial values should be within the vlim
  !--------------------------------------------------------------------
  plr(1:4)=ranrange(1:4,1)
  pla(1:5)=ranrange(5:9,1)
  !
  ! First pointL The left corner!
  mxyz(:,:) = 0.0_FREAL
  irc=poly2cart(plr(1:4),pla(1:5),mxyz(1:3,2:5))
  if(irc.lt.0) then
     call message(MESERRO,"[XY4P] Cannot perform conversion 884723.")
     stop 3
  end if
  call getbonds(4,mxyz(1:3,2:5),xy(1:4))
  call getangles(4,mxyz(1:3,2:5),xy(1:4),yxy(1:6))
  cyxy(1:6)=cos(yxy(1:6))
  vmpotin(1,1:4) = xy(1:4)
  vmpotin(1,5:10)= cyxy(1:6)
  call metpot4_vint(vmpotin,vmpotout)
  if(vmpotout(1).gt.maxseenvmpot) maxseenvmpot=vmpotout(1)
  if(vmpotout(1).lt.LEPS) then
     write(iunoutput,'(A)') "#Error! Input start coord give a too small vlimit vs LEPS!"
     stop 66
  end if
  !
  if(vmpotout(1).gt.vlimit) then
     write(iunoutput,'(A)') "#Error! Input start coord give a too small vlimit vs start coord potential"
     write(iunoutput,'(A)') "We do not stop here but be aware!!"
     write(iunoutput,*) vmpotout(1), vlimit
     !stop 67
  end if
  !
  ! The code below is suppressed in favor of a function - CHECK TODO

!  do j1=1,2
!  do j2=1,2
!  do j3=1,2
!  do j4=1,2
!  do j5=1,2
!  do j6=1,2
!  do j7=1,2
!  do j8=1,2
!  do j9=1,2
!     plr(1) = ranrange(1,j1)
!     plr(2) = ranrange(2,j2)
!     plr(3) = ranrange(3,j3)
!     plr(4) = ranrange(4,j4)
!     pla(1) = ranrange(5,j5)
!     pla(2) = ranrange(6,j6)
!     pla(3) = ranrange(7,j7)
!     pla(4) = ranrange(8,j8)
!     pla(5) = ranrange(9,j9)
!     irc=poly2cart(plr(1:4),pla(1:5),mxyz(1:3,2:5))
!     if(irc.lt.0) then
!        call message(MESERRO,"[XY4P] Cannot perform conversion 884723.")
!        stop 3
!     end if
!     call getbonds(4,mxyz(1:3,2:5),xy(1:4))
!     call getangles(4,mxyz(1:3,2:5),xy(1:4),yxy(1:6))
!     cyxy(1:6)=cos(yxy(1:6))
!     vmpotin(1,1:4) = xy(1:4)
!     vmpotin(1,5:10)= cyxy(1:6)
!     call metpot4_vint(vmpotin,vmpotout)
!     if(vmpotout(1).gt.maxseenvmpot) maxseenvmpot=vmpotout(1)
!     if(vmpotout(1).lt.LEPS) then
!        write(iunoutput,'(A)') "#Error! Input start coord give a too small vlimitvs LEPS!"
!        stop 76
!     end if
!     !
!     if(vmpotout(1).gt.vlimit) then
!        write(iunoutput,'(A)') "#Error! Input start coord give a too small vlimit vs start coord potential"
!        write(iunoutput,*) vmpotout(1), vlimit
!        stop 77
!     end if
!     if( XY4P_debug) then
!        write(iunoutput,'("#",10(1X,F11.4))') plr(1:4),pla(1:5)/LPI*180.0_FREAL,vmpotout(1)
!     end if
!  end do
!  end do
!  end do
!  end do
!  end do
!  end do
!  end do
!  end do
!  end do
!
! The function goes here!
  vmpotout(1) = chk_rangevol(ranrange(1:9,1:2),cvolume(1:4))
  !
  ! For the sake of documentation which will never appear (yep maybe):
  ! This function will check:
  ! - Poly 2 cart convertion: Will stop.
  ! - Will indeed update the MaxSeenPotential
  ! - Will also stop if the potential is below the Local Machine Precision.
  ! - will give a negative potential if one of the points are over the gloal vlimit.
  !
  if( XY4P_debug) then
     write(iunoutput,'("[x-xy4-srcgrid] Debug  Vols:",1X,4(1X,E14.6),1X,F11.4)') cvolume(1:4),vmpotout(1)
  end if
  !
  ! Check potential
  if(vmpotout(1).lt.0.0_FREAL) then
     write(iunoutput,'("[x-xy4-srcgrid] WARNING  Vols:",1X,4(1X,E14.6),1X,F11.4)') cvolume(1:4),vmpotout(1)
  else
     write(iunoutput,'("[x-xy4-srcgrid] info Volumes and largest pot:",1X,4(1X,E14.6),1X,F11.4)') cvolume(1:4),vmpotout(1)
  end if
  !

  !
  ! Search symmetrically around the center
  !---------------------------------------
  ranpars(1:9,2)=0.0_FREAL
  ranpars(1:9,3)=0.0_FREAL
  !
  notdone = .TRUE.
  j=0
  do while (notdone)
     notdone = .FALSE.
     j=j+1
     do i=1,4
       ! Positive test
       ranrange(i,2)=ranrange(i,2)+ranpars(i,1)  
! TODO WARINING !!!
! ! Here we use .gt. ut we might want to use a REAL comparison
! ! More over: we do not need to get metpot4 if we are out of range!!
! !
! ! Same for all following coordinates!!!!
! !
!
 !      if(ranrange(i,2).gt.sranrange(i,2)) then
 !        write(iunoutput,'("Coord",1X,I3,1X,"Hit the limit +",1X,F14.6)') i,ranrange(i,2)
 !        ranrange(i,2)=ranrange(i,2)-ranpars(i,1)
 !        ranpars(i,3)=1.0_FREAL
 !      else
          vmpotout(1) = chk_rangevol(ranrange(1:9,1:2),cvolume(1:4))
          if(vmpotout(1).lt.0.0_FREAL) then
             ranrange(i,2)=ranrange(i,2)-ranpars(i,1)
             ranpars(i,3)=1.0_FREAL 
             if( XY4P_debug) then
                write(iunoutput,'(A,1X,I2)') "#NOTE: point discarded! R+",i
             end if
          end if
 !      end if
       ! Negative test
       ranrange(i,1)=ranrange(i,1)-ranpars(i,1)
 !      if(ranrange(i,1).lt.sranrange(i,1)) then
 !        write(iunoutput,'("Coord",1X,I3,1X,"Hit the limit -",1X,F14.6)') i,ranrange(i,1)
 !        ranrange(i,1)=ranrange(i,1)+ranpars(i,1)
 !        ranpars(i,2)=1.0_FREAL
 !      else
          vmpotout(1) = chk_rangevol(ranrange(1:9,1:2),cvolume(1:4))
          if(vmpotout(1).lt.0.0_FREAL) then
             ranrange(i,1)=ranrange(i,1)+ranpars(i,1)
             ranpars(i,2)=1.0_FREAL
             if( XY4P_debug) then
                write(iunoutput,'(A,1X,I2)') "#NOTE: point discarded! R-",i
             end if
          end if
 !      end if
     end do
     !
     !
     do i=5,9
       ! Positive test
       ranrange(i,2)=ranrange(i,2)+ranpars(i,1)*LPI/180.0_FREAL
       !
!       if(ranrange(i,2).ge.sranrange(i,2)) then
!         write(iunoutput,'("Coord",1X,I3,1X,"Hit the limit +",1X,F14.6)') i,ranrange(i,2)
!         ranrange(i,2)=ranrange(i,2)-ranpars(i,1)*LPI/180.0_FREAL
!         ranpars(i,3)=1.0_FREAL
!       end if
       !
       vmpotout(1) = chk_rangevol(ranrange(1:9,1:2),cvolume(1:4))
       if(vmpotout(1).lt.0.0_FREAL) then
         ranrange(i,2)=ranrange(i,2)-ranpars(i,1)*LPI/180.0_FREAL
         ranpars(i,3)=1.0_FREAL 
         if( XY4P_debug) then
            write(iunoutput,'(A,1X,I2)') "#NOTE: point discarded! A+",i
         end if
       end if
       ! Negative test
       ranrange(i,1)=ranrange(i,1)-ranpars(i,1)*LPI/180.0_FREAL
 !      if(ranrange(i,2).le.sranrange(i,2)) then
 !        write(iunoutput,'("Coord",1X,I3,1X,"Hit the limit +",1X,F14.6)') i,ranrange(i,2)
 !        ranrange(i,2)=ranrange(i,2)+ranpars(i,1)*LPI/180.0_FREAL
 !        ranpars(i,2)=1.0_FREAL
 !      end if
       vmpotout(1) = chk_rangevol(ranrange(1:9,1:2),cvolume(1:4))
       if(vmpotout(1).lt.0.0_FREAL) then
         ranrange(i,1)=ranrange(i,1)+ranpars(i,1)*LPI/180.0_FREAL
         ranpars(i,2)=1.0_FREAL
         if( XY4P_debug) then
            write(iunoutput,'(A,1X,I2)') "#NOTE: point discarded! A-",i
         end if
       end if
     end do
     !
     ! test if all are within limits
     !------------------------------
     if ((sum(ranpars(1:9,2)).lt.8.5_FREAL).or.(sum(ranpars(1:9,3)).lt.8.5_FREAL)) then
        notdone = .TRUE.
        !write(iunoutput,'("# NO Check is: Low bound",1X,F12.4,1X,"High bound",1X,F12.4)') sum(ranpars(1:9,2)),sum(ranpars(1:9,3))
        !write(iunoutput,'("# NO Low  Vec is:",1X,9(1X,F10.2))') ranpars(1:9,2)
        !write(iunoutput,'("# NO High Vec is:",1X,9(1X,F10.2))') ranpars(1:9,3)
        !print *,(sum(ranpars(1:9,2)).lt.8.5_FREAL),(sum(ranpars(1:9,3)).lt.8.5_FREAL)
     else
        notdone = .FALSE.
        write(iunoutput,'("# We Are at the limits")')
        write(iunoutput,'("# Check is: Low bound",1X,F12.4,1X,"High bound",1X,F12.4)') sum(ranpars(1:9,2)),sum(ranpars(1:9,3))
        write(iunoutput,'("# Low  Vec is:",1X,9(1X,F10.2))') ranpars(1:9,2)
        write(iunoutput,'("# High Vec is:",1X,9(1X,F10.2))') ranpars(1:9,3)
        !print *,(sum(ranpars(1:9,2)).lt.8.5_FREAL),(sum(ranpars(1:9,3)).lt.8.5_FREAL)
     end if
     !
     !
     ! Limit the search anyway!
     if (j.gt.100000) then
        notdone = .FALSE.
        write(iunoutput,'("Search will stop due to J")')
     end if
  end do
  !
  !
  !
  maxseenvmpot = 0.0_FREAL
  vmpotout(1) = chk_rangevol(ranrange(1:9,1:2),cvolume(1:4))
  write(iunoutput,'("Latest Lower:")')
  write(iunoutput,'(1X,9(1X,F11.4))') ranrange(1:4,1),ranrange(5:9,1)/LPI*180.0_FREAL
  write(iunoutput,'("Latest Upper:")') 
  write(iunoutput,'(1X,9(1X,F11.4))') ranrange(1:4,2),ranrange(5:9,2)/LPI*180.0_FREAL
  write(iunoutput,'("Largest potential:",1X,I8,1X,F11.4)') j,vmpotout(1)
  write(iunoutput,'("Volumes (Tot, Rad, Ang, Fact):",1X,4(1X,E14.4))') cvolume(1:4)
  !
  ! Save last volume info
  svolume(:)=cvolume(:)
  !
  ! Check that all points are inside the potential
  !-----------------------------------------------
  ! TODO!

  !
  !
  write(iunoutput,'("[x-out-result] Statistics")')
  write(iunoutput,'("Largest METPOT seen so far:",1X,F12.4)') maxseenvmpot
  !
  !

  999 continue

contains
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function poly2cart
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function poly2cart(lr,la,lx,debug)
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
    implicit none
    real(FREAL), dimension(:), intent(in) :: lr,la
    real(FREAL), dimension(:,:), intent(out) :: lx
    !
    logical, optional :: debug
    logical :: ldebug
    !
    ! Start with no errors
    !---------------------
    poly2cart=0
    !
    ! Set default debug level
    !------------------------
    if (present(debug)) then
       ldebug=debug
    else
       ldebug=.false.
    end if
    !
    ! checks
    !-------
    if (ldebug) then
       ! Sizes
       !------
       if(size(lr).gt.4) then
          call message(MESWARN,"[poly2cart]: Poly R vector bigger than 4: is that ok?.")
       end if
       if(size(lr).lt.4) then
          call message(MESERRO,"[poly2cart]: Poly R vector smaller than 4. We STOP.")
          poly2cart=-1
       end if
       if(size(la).gt.5) then
          call message(MESWARN,"[poly2cart]: Poly A vector bigger than 5: is that ok?.")
       end if
       if(size(la).lt.5) then
          call message(MESERRO,"[poly2cart]: Poly A vector smaller than 5. We STOP.")
          poly2cart=-1
          return
       end if
       if (size(lx,1).ne.3) then
          call message(MESERRO,"[poly2cart]: we need 3D cartesian input coordinates. We STOP.")
          poly2cart=-1
          return
       end if
       if (size(lx,2).gt.4) then
          call message(MESWARN,"[poly2cart]: Cartesian atoms vector bigger then 4.")
       end if
       if (size(lx,2).lt.4) then
          call message(MESERRO,"[poly2cart]: Cartesian atoms vector smaller then 4. We STOP.")
          poly2cart=-1
          return
       end if
       ! ranges
       !-------
       if(maxval(la(1:3)).gt.LPI) then
          call message(MESERRO,"[poly2cart]: One of the Theta is out of range (180). We STOP.")
          poly2cart=-1
          return
       end if
       !if(maxval(la(4:5)).gt.2.0_FREAL*LPI) then
       if(maxval(la(4:5)).gt.LPI) then
          call message(MESERRO,"[poly2cart]: One of the Phi is out of range (360). We STOP.")
          poly2cart=-1
          return
       end if
       if(minval(la(1:5)).lt.LEPS) then
          call message(MESERRO,"[poly2cart]: One of the Theta or Phi is out of range (0). We STOP.")
          poly2cart=-1
          return
       end if
    end if
    !
    ! H4
    lx(1,4) = 0.0_FREAL
    lx(2,4) = 0.0_FREAL
    lx(3,4) = lr(1)
    !
    ! H3
    lx(1,3) = -lr(2)*sin(la(1))
    lx(2,3) = 0.0_FREAL
    lx(3,3) = lr(2)*cos(la(1))
    !
    ! H2
    lx(1,2) =  lr(3)*sin(la(2))*cos(la(4))
    lx(2,2) =  lr(3)*sin(la(2))*sin(la(4))
    lx(3,2) =  lr(3)*cos(la(2))
    !
    ! H1
    lx(1,1) =  lr(4)*sin(la(3))*cos(la(5))
    lx(2,1) = -lr(4)*sin(la(3))*sin(la(5))
    lx(3,1) =  lr(4)*cos(la(3))
    !
    ! C0
    !lx()
    !
    ! Write Partial Data
    !-------------------
    if (ldebug) then
       write(iunoutput,'("[x-poly2cart] 4 Y4 atoms")')
       write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y4',lx(1:3,4)
       write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y3',lx(1:3,3)
       write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y2',lx(1:3,2)
       write(iunoutput,'(1X,A2,1X,3(2X,F12.4))') 'Y1',lx(1:3,1)
    end if
    !
  end function poly2cart
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
!H Function XY4PolySphere_init()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function XY4PolySphere_init()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the XY4PolySphere module.
!H
!H-----------------------------------------------------------------------------
!H
!
    XY4PolySphere_init = 0
    !
    ! Initialize messages and ourself
    !--------------------------------
    call messages_init()
    !
    ! File I/O stuff initialization
    !------------------------------
    irc = baseio_init()
    if(irc.lt.0) then
       call message(MESERRO,"[XY4P] Cannot initialize baseIO.")
       XY4PolySphere_init = -1
       return
    end if
    !
    !
    ! Parse command line
    !-------------------
    call pcmd_iargc(nargs)
    allocate(cargs(nargs),STAT=irc)
    if(.not.irc.eq.0) then
       call message(MESERRO,"[XY4P] Cannot allocate for command line.")
       XY4PolySphere_init = -1
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
          LOCXY4PolySphere_error = 3
          call message(MESERRO,"[XY4P] Cannot Open Input File")
          XY4PolySphere_init = -1
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
          LOCXY4PolySphere_error = 3
          call message(MESERRO,"[XY4P] Cannot Open Output File")
          XY4PolySphere_init = -1
          return
       end if
    end if
    !
    ! Other input stuff
    !------------------
    ! The HELP
    !---------
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-h',params)
    if(irc.gt.0) then
       irc = index(params,'full')
       if(irc.gt.0) then
          write (*,*) 'Full Help not yet Available...'
          stop 0
       else
          write (*,*) 
          write (*,*) 'XY4P Simple Help: Long version to be implemented'
          write (*,*) 
          write (*,*) ' Options:'
          write (*,*) 
          write (*,*) '-i [<inputfilename>]'
          write (*,*) '-o [<outputfilename>]'
          write (*,*) 
          write (*,*) '-d <N>   Set Debug Level to N'
          write (*,*) '-v <N>   Set Verbosity Level to N'
          write (*,*) 
          write (*,*) '-V <F>   Set top limit of the potential to exclude points'
          write (*,*) 
          write (*,*) '-r <resultfilename>'
          write (*,*) '-x <xyzfilename> seen coordinated saved in XYZ format'
          write (*,*) '-R <random seed>  Integer. Set random seed.'
          write (*,*) 
          stop 0
       end if
    end if
    ! Debug
    !------
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-d',params)
    if(irc.lt.0) then
       call message(MESERRO,"[rdinput] -d Require a logical value parameter (T/F).")
       stop 1
    else if(irc.gt.0) then
       params = trim(params)
       if(scan(params," ").eq.0) then
          call message(MESERRO,"[rdinput] -d Require a logical value parameter (T/F).")
          stop 1
       end if
       read(params,'(L1)',IOSTAT=irc) XY4P_debug
       if(irc.ne.0) then
          call message(MESERRO,"[rdinput] Option error -d Require a logical value parameter (T/F).")
          stop 1
       end if
    end if
    ! Verbosity
    !----------
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-v',params)
    if(irc.lt.0) then
       call message(MESERRO,"[rdinput] -v Require an integer value parameter.")
       stop 1
    else if(irc.gt.0) then
       params = trim(params)
       if(scan(params," ").eq.0) then
          call message(MESERRO,"[rdinput] -v Require an integer value parameter.")
          stop 1
       end if
       read(params,'(I8)',IOSTAT=irc) XY4P_verbose
       if(irc.ne.0) then
          call message(MESERRO,"[rdinput] Option error -v Require an integer value parameter.")
          stop 1
       end if
    end if

    ! Random Seed
    !------------
    XY4P_rsset(1) = 0
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-R',params)
    if(irc.lt.0) then
       call message(MESERRO,"[rdinput] -R Require an integer value parameter.")
       stop 1
    else if(irc.gt.0) then
       params = trim(params)
       if(scan(params," ").eq.0) then
          call message(MESERRO,"[rdinput] -R Require an integer value parameter.")
          stop 1
       end if
       read(params,'(I8)',IOSTAT=irc) XY4P_rsset(1)
       if(irc.ne.0) then
          call message(MESERRO,"[rdinput] Option error -R Require an integer value parameter.")
          stop 1
       end if
    end if

    ! Potential Limit
    !----------------
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-V',params)
    if(irc.lt.0) then
       call message(MESERRO,"[rdinput] -V Require an real value parameter.")
       stop 1
    else if(irc.gt.0) then
       params = trim(params)
       if(scan(params," ").eq.0) then
          call message(MESERRO,"[rdinput] -V Require an real value parameter.")
          stop 1
       end if
       !read(params,'(G)',IOSTAT=irc) vlimit
       !read(params,'(G24.12)',IOSTAT=irc) vlimit
       read(params,*,IOSTAT=irc) vlimit
       if(irc.ne.0) then
          call message_value(MESERRO,"[rdinput] Option error -V Require an real value parameter.",vlimit)
          stop 1
       end if
       if (vlimit.lt.0.0_FREAL) then
          call message_value(MESERRO,"[rdinput] -V Require an real POSITIVE value parameter.",vlimit)
          stop 1
       end if
       if (abs(vlimit).lt.LEPS) then
          call message_value(MESERRO,"[rdinput] -V Require an real POSITIVE/NOT-ZERO value parameter.",vlimit)
          stop 1
       end if
    end if
    !
    ! RESULTS file name
    iunres =0
    fileres=""
    irc = pcmd_checkarg(nargs,cargs,'-r',params)
    if(irc.lt.0) then
       call message(MESERRO,"[rdinput] -r Require a valid file name.")
       stop 1
    else if(irc.gt.0) then
       params = trim(params)
       if(scan(params," ").eq.0) then
          call message(MESERRO,"[rdinput] -r Require a valid file name.")
          stop 1
       end if
       fileres = adjustl(params)
       fileres = trim(fileres)
          irc = baseio_open(iunres,fileres,stat='UNKNOWN')
          if(irc.lt.0) then
             call message(MESERRO,"[rdinput] -r Require a valid file name. Faild Open Statment.")
             stop 1
          end if
    end if
    !
    ! XYZ file name
    iunxyz =0
    filexyz=""
    irc = pcmd_checkarg(nargs,cargs,'-x',params)
    if(irc.lt.0) then
       call message(MESERRO,"[rdinput] -x Require a valid file name.")
       stop 1
    else if(irc.gt.0) then
       params = trim(params)
       if(scan(params," ").eq.0) then
          call message(MESERRO,"[rdinput] -x Require a valid file name.")
          stop 1
       end if
       filexyz = adjustl(params)
       filexyz = trim(filexyz)
          irc = baseio_open(iunxyz,filexyz,stat='UNKNOWN')
          if(irc.lt.0) then
             call message(MESERRO,"[rdinput] -x Require a valid file name. Faild Open Statment.")
             stop 1
          end if
    end if
    !
    ! (RE)Initialize messages and ourself
    !------------------------------------
    call messages_init(iunoutput,iunoutput,0,iunoutput)
    !
    ! We need the molecule section
    !-----------------------------
    irc=molecule_init(.true.)
    if(irc.lt.0) then
       call message(MESERRO,"[XY4P] Cannot initialize molecule module..")
       XY4PolySphere_init = -1
       return
    end if
    !
    ! METPOT4_initialization
    !-----------------------
    call metpot4_initpot()
    !
    return
    !
  end function XY4PolySphere_init
!
!H
!H-----------------------------------------------------------------------------
!H logical function XY4PolySphere_isinit()
!H-----------------------------------------------------------------------------
!H
!
  logical function XY4PolySphere_isinit()
!
!H
!H-----------------------------------------------------------------------------
!H We merelyreturn the initialization flag.
!H-----------------------------------------------------------------------------
!H
!
    !
    XY4PolySphere_isinit = .TRUE.
    !
  end function XY4PolySphere_isinit
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H character(FLCHARS) function XY4PolySphere_geterror()
!H-----------------------------------------------------------------------------
!H
!
  character(FLCHARS) function XY4PolySphere_geterror(code)
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
       XY4PolySphere_geterror = "No Errors."
    case(1_FINT)
       XY4PolySphere_geterror = "Not documented error or generic error.."
    case(2_FINT)
       XY4PolySphere_geterror = "Double Error: Why did you get this message?."
    case default
       XY4PolySphere_geterror = "Not documented error or generic error.."
    end select
    !
    return
    !
  end function XY4PolySphere_geterror
!
!H
!H-----------------------------------------------------------------------------
!H
!


!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H real(FREAL) function chk_rangevol(ranrange,vol)
!H-----------------------------------------------------------------------------
!H
!H In simple words.. This function return the potential: with a negative
!H sign if it is over the given limit (which is read from the external vlimit)
!H It does also return a volume for the given coordinates. It does distinguish
!H the r(1:4) and the a(5:9) via an hardcoded scale factor which is now=1.0.
!H Volume is given back as a vector=(total,rvolume,avolume,scalefactor)
!H
!
  real(FREAL) function chk_rangevol(ranrange,vol)
    !
    ! Scale factor for volumes if needed.
    ! Please note that it DOES scale volumes not
    ! coordinates!
    integer(FINT), parameter :: vscale = 1.0_FREAL
    real(FREAL), dimension(9,2), intent(in) :: ranrange
    real(FREAL), intent(out), dimension(4) :: vol
    !
    real(FREAL), dimension(9) ::  cpoint
    real(FREAL), dimension(3,5) :: lxyz
    real(FREAL), dimension(4) :: xy
    real(FREAL), dimension(6) :: yxy
    !real(FREAL) :: lmaxp
    !
    integer(FINT) :: j1,j2,j3,j4,j5,j6,j7,j8,j9
    real(FREAL) :: tvolume, rvolume, avolume, lvmax
    real(FREAL), dimension(9) :: tmpvolume
    !
    !lmaxp = 0.0_FREAL
    lvmax = 0.0_FREAL
    !
    !
    ! Evaluate volume of the cood hypercube
    !--------------------------------------
    !rvolume = (ranrange(1,2)-ranrange(1,1))*&
    !        & (ranrange(2,2)-ranrange(2,1))*&
    !        & (ranrange(3,2)-ranrange(3,1))*&
    !        & (ranrange(4,2)-ranrange(4,1))
    ! DEUG MY F90 knowledge and performances:
    ! yep.. here not too much of a performance
    ! will e tested... ut we give us a chance..
    ! [Whe you address a couple of languages:
    !  How do you write the here sound of ch?
    !  Chance - schanze - scianse - what the french now?]
    tmpvolume(1:4) = ranrange(1:4,2) - ranrange(1:4,1)
    rvolume = product(tmpvolume(1:4))
    tmpvolume(1:5) = ranrange(5:9,2) - ranrange(5:9,1)
    avolume = product(tmpvolume(1:5))
    !
    tvolume = rvolume * vscale * avolume
    !
    vol(1) = tvolume
    vol(2) = rvolume
    vol(3) = avolume
    vol(4) = vscale
    !
    !
    ! The code here is not efficient but uses
    ! only one variables... prolem shited to the function!
    !
    if( XY4P_debug) then
       write(iunoutput,'("# Debug GRID Low  point: ",1X,9(1X,E12.4))') ranrange(1:4,1), ranrange(5:9,1)
       write(iunoutput,'("# Debug GRID High point: ",1X,9(1X,E12.4))') ranrange(1:4,2), ranrange(5:9,2)
    end if
    !
    do j1=1,2
    do j2=1,2
    do j3=1,2
    do j4=1,2
    do j5=1,2
    do j6=1,2
    do j7=1,2
    do j8=1,2
    do j9=1,2
       cpoint(1) = ranrange(1,j1)
       cpoint(2) = ranrange(2,j2)
       cpoint(3) = ranrange(3,j3)
       cpoint(4) = ranrange(4,j4)
       cpoint(5) = ranrange(5,j5)
       cpoint(6) = ranrange(6,j6)
       cpoint(7) = ranrange(7,j7)
       cpoint(8) = ranrange(8,j8)
       cpoint(9) = ranrange(9,j9)
       lxyz(:,:) = 0.0_FREAL
       irc=poly2cart(cpoint(1:4),cpoint(5:9),lxyz(1:3,2:5))
       if(irc.lt.0) then
          call message(MESERRO,"[XY4P] Cannot perform conversion 84583723.")
          stop 3
       end if
       call getbonds(4,lxyz(1:3,2:5),xy(1:4))
       call getangles(4,lxyz(1:3,2:5),xy(1:4),yxy(1:6))
       cyxy(1:6)=cos(yxy(1:6))
       vmpotin(1,1:4) = xy(1:4)
       vmpotin(1,5:10)= cyxy(1:6)
       call metpot4_vint(vmpotin,vmpotout)
       if(vmpotout(1).gt.maxseenvmpot) maxseenvmpot=vmpotout(1)
       if(vmpotout(1).lt.LEPS) then
          write(iunoutput,'(A)') "#Error! Input start coord gives a too small vlimitvs LEPS!"
          stop 76
       end if
       !
       if(vmpotout(1).gt.lvmax) lvmax = vmpotout(1)
       !
       !if( XY4P_debug) then
       !   write(iunoutput,'("#FD: ",10(1X,F11.4))') cpoint(1:4),cpoint(5:9)/LPI*180.0_FREAL,vmpotout(1)
       !end if
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    !
    if(lvmax.gt.vlimit) then
       chk_rangevol = -lvmax
    else
       chk_rangevol = lvmax
    end if
    !
    !if( XY4P_debug) then
    !   write(iunoutput,'("# Debug GRID Vols:",1X,5(1X,E14.6))') vol(1:4),chk_rangevol
    !end if
    !
    !
    !
end function chk_rangevol
  !
  !
end program XY4PolySphere
