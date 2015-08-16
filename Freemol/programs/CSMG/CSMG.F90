!
!H
!H-----------------------------------------------------------------------------
!H This File is part of the Frimol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H PROGRAM CSMG  Frimol by F.Mariotti: (c) F.Mariotti
!H-----------------------------------------------------------------------------
!H $Id: CSMG.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
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
program CSMG
  !
  use vartypes
  use messages
  use csmmod
  use baseio
  use molecule
  use mathtools
  use modcsmmin
  use csm_minuit
  !
  implicit none
  !
!------------------------------------------------------------------------------
  ! DECLARATIONS:
!------------------------------------------------------------------------------
  !--------------------
  
  ! Procedures
  !------------------
  !integer(FINT) CSMG_init
  !logical CSMG_isinit
  !integer(FINT) CSMG_
  !
  ! LOCAL DECLARATIONS:
  !--------------------
  !
  logical, save :: LocCSMG_isinit = .false.
  !
  ! error variables
  !----------------
  integer(FINT), save :: LocCSMG_error
  character(FLCHARS), save :: LocCSMG_error_message
  !
  ! I/O Files
  !----------
  character(FLCHARS), save :: finput, foutput, filelog, fileerr
  integer(FINT), save :: iunin, iunout, iunlog, iunerr, iunsave, mininp
  !
  ! dimension stuff
  !----------------
  integer(FINT), save :: nat, lnsymop
  integer(FINT) :: ierflg,i
  !
  ! Others
  !-------
  integer(FINT) :: irc
  character(FLCHARS) :: errmsg
  ! a local PI value
  real(FREAL), parameter :: LPI = 3.1415926535897932384626433832795028841971694_FREAL
  !
  ! The Integrals
  !--------------
  real(FREAL) :: nelec
  real(FREAL) :: csmnorm
  real(FREAL) :: srrii, srrij, srhorho
  real(FREAL), dimension(3) :: masscenter
  real(FREAL), dimension(3,3) :: rotinertia
  real(FREAL), allocatable, dimension(:) :: rhoint
  real(FREAL), allocatable, dimension(:) :: rrii
  real(FREAL), allocatable, dimension(:) :: rrij
  real(FREAL), allocatable, dimension(:,:,:) :: rrimj
  real(FREAL), allocatable, dimension(:,:,:) :: rrimpj
  !
  ! CSM
  !----
  real(FREAL) :: CSMI,CSM, CSMmain, CSMx, CSMy, CSMz, CSMmin
  !
  ! pars
  !-----
  real(FREAL), dimension(6) :: curparmin, parsteps
  real(FREAL) :: ftol, resid
  integer(FINT) :: iter
  integer(FINT), dimension(6) :: parvect
  character(FLCHARS), dimension(6) :: parlabel
  !
  ! pars for minuit
  !----------------
  integer(FINT) :: csm_ipar_minuit
  real(FREAL) :: csm_save_minuit, csm_error_minuit, csm_dummy_minuit, csm_svpar_minuit
  real(FREAL), dimension(10) :: csm_temp_minuit
  real(FREAL), dimension(6) :: csm_pars_minuit,csm_svbestpars
  !
  ! Dummies
  !--------
  real(FREAL) :: ra,rb,rg
  character(FLCHARS) :: chardummy
  !
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Here we start
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !
!------------------------------------------------------------------------------
  ! Initialization
!------------------------------------------------------------------------------
  ! Get messages initialization as default
  ! it will be changed later in the code
  !-------------------------------------
  call messages_init()
  !
  ! Just follow the rules and call this program init
  !-------------------------------------------------
  irc = CSMG_init()
  !
  ! File I/O stuff initialization
  !------------------------------
  irc = baseio_init()
  !
  ! We get initialized the CSM module as well
  ! Pay attention that we will get a csm DATA module as well
  !---------------------------------------------------------
  irc = csmmod_init()
  !
  ! We get input data and parameters
  !---------------------------------
  irc = csmmod_rdinput(iunin,finput,iunout,foutput,iunlog,filelog,iunerr,fileerr)
  if (irc.lt.0) then
     irc = csmmod_geterror() ! Fetch error number
     errmsg = csmmod_geterror(irc) ! Fetch error code
     call message(MESERRO,errmsg)
  end if
  !stop "here"
  !
  ! We get back some input pars like:
  ! Number of atoms and symmetry operations
  !----------------------------------------
  irc = csmmod_getnat()
  nat = irc
  if(irc.le.0) then
     call message(MESERRO,"[CSMG] what is nat than?")
     stop 1
  end if
  irc = csmmod_getnsymop()
  lnsymop = irc
  if(irc.le.0) then
     call message(MESERRO,"[CSMG] what is lnsymop than?")
     stop 1
  end if
  !
  ! Print the actual coordinates at input
  !--------------------------------------
  call molecule_print(iunout,'mldfromau')
  !
  ! We adjust the molecule a bit
  !-----------------------------
  ! We check first that is not required a fixed position
  !-----------------------------------------------------
  if(.not.csm_noshift) then
     !
     ! Move to the Center of Mass
     !---------------------------
     molecule_atmchar(1:nat) = real(molecule_zchar(1:nat),FREAL)
     call mathtools_masscenter(molecule_xyz(1:3,1:nat),&
          &molecule_atmchar(1:nat),masscenter,.true.)
     write(iunout,'("# Molecule moved to the center of mass")') 
     write(iunout,'("# Shift Vector is:",1X,3(1X,F16.8))') masscenter
     !
     ! Rotate to Inertia Axis
     !-----------------------
     if(.not.csm_noinertia) then
        call mathtools_inertia(molecule_xyz(1:3,1:nat),molecule_atmchar(1:nat),rotinertia,.true.)
        write(iunout,'("# Molecule Rotated to match inertia axes")') 
     end if
     !
     call molecule_print(iunout,'mldfromau')
     !
  end if
  write(iunout,'("# Initial Integrals follow ...")') 
  !
  ! We compute The Rho
  !-------------------
  allocate(rhoint(nat),STAT=irc)
  call csmmod_getrho(rhoint)
  nelec = sum(rhoint(:))
  write(iunout,'("# INT: Number of electrons from integrated RHO:",1X,F18.11)') nelec
  !
  !
  ! We compute first g_i*g_i in a vector(1:n)
  !------------------------------------------
  allocate(rrii(nat),STAT=irc)
  call csmmod_getrrii(rrii)
  srrii = sum(rrii(:))
  !
  ! We compute g_i*g_j in a vector as well
  !---------------------------------------
  allocate(rrij((nat*(nat-1_FINT))/2_FINT),STAT=irc)
  call csmmod_getrrij(rrij)
  srrij = sum(rrij(:))
  !
  srhorho = srrii +  2.0_FREAL * srrij
  write(iunout,'("# INT: \int \rho \rho dr: ",1X,F18.11)') srhorho
  !
  ! then the other integrals(g_i S_m g_j) on a matrix (i,m,j)
  !-------------------------------------------------------------
  allocate(rrimj(nat,lnsymop,nat),STAT=irc)
  call csmmod_getimj(rrimj)
  write(iunout,'("# INT: \int \rho \tilde{\rho} dr/ng: ",1X,F18.11)') sum(rrimj(:,:,:))/real(lnsymop,FREAL)
  !
  ! then the other integrals(g_i S_m' g_j) on a matrix (i,m',j)
  !-------------------------------------------------------------
  allocate(rrimpj(nat,lnsymop,nat),STAT=irc)
  call csmmod_getimpj(rrimpj)
  write(iunout,'("# INT: \int \rho \tilde{\rho}m dr/ng: ",1X,F18.11)') sum(rrimpj(:,:,:))/real(lnsymop,FREAL)
  !
  ! DEBUGGING
  !-------------------------------------------
  ! New Equation: See paper draft 07 July 2005
  !-------------------------------------------
  ! Pseudo code:
  ! D = (ng-1)/(ng) * (sum(rrimj(:,:,:)))/(srhorho)
  !-----------------------------------------------
  ! terms
  write(iunout,'("#  DEBUG     \int rr  =",1X,F18.11)') srhorho
  !
  CSMI = sum(rrimj(:,:,:))/real(lnsymop,FREAL)
  write(iunout,'("#  DEBUG     \int rr~ =",1X,F18.11)') CSMI
  !
  CSMI = (srhorho+sum(rrimpj(:,:,:)))/real(lnsymop,FREAL)
  write(iunout,'("#  DEBUG     \int rr~ =",1X,F18.11)') CSMI
  !
  CSMI = 1 - (sum(rrimj(:,:,:))/real(lnsymop,FREAL))/(srhorho)
  write(iunout,'("#  DEBUG     CSM value =",1X,F18.11)') CSMI
  !
  CSMI = ((real(lnsymop,FREAL)-1.0_FREAL)/(real(lnsymop,FREAL))) - ((sum(rrimpj(:,:,:)))/(real(lnsymop,FREAL)*srhorho))
  write(iunout,'("#  DEBUG     CSM value =",1X,F18.11)') CSMI
  !
  ! the first initial CSM Value and summary
  !----------------------------------------
  csmnorm = srhorho
  CSMI = ((real(lnsymop,FREAL)-1.0_FREAL)/(real(lnsymop,FREAL))) - ((sum(rrimpj(:,:,:)))/(real(lnsymop,FREAL)*srhorho))
  !
  write(iunout,'("# Number of electrons =",1X,F18.11)') nelec
  write(iunout,'("# Initial rho^2 value =",1X,F18.11)') srhorho
  write(iunout,'("#  Initial rrii value =",1X,F18.11)') srrii
  write(iunout,'("#  Initial rrij value =",1X,F18.11)') 2.0_FREAL * srrij
  write(iunout,'("#  Initial imj  value =",1X,F18.11)') (1.0_FREAL/real(lnsymop,FREAL)) * sum(rrimj(:,:,:))
  write(iunout,'("#  Initial impj value =",1X,F18.11)') (1.0_FREAL/real(lnsymop,FREAL)) * sum(rrimpj(:,:,:))
  write(iunout,'("#  Initial CSM  value =",1X,F18.11)') CSMI
  write(iunout,'("#  Initial CSM% value =",1X,F18.11)') CSMI*100.0_FREAL
  !
  ! Print the actual coordinates after eventual shift or rotation
  !--------------------------------------------------------------
  call molecule_print(iunout,'mldfromau')
  !
  ! We do check wether we did ask for single step
  !----------------------------------------------
  if(csm_onestep) then
     call message(MESOUT,"[CSMG] Required onestep: no optimization.")
     call message(MESOUT,"[CSMG] We stop right now.")
     stop 0
  end if
  !
  ! We start to run the minimization procedure
  !-------------------------------------------
  irc = modcsmmin_init(iunout,srhorho,csmnorm,lnsymop,.true.)
  !
  !
  ! Option Read Minuit
  !-----------------------
  if(csm_readminuit) then
     write(iunout,'("# Request for Minuit: Read from input file.")') 
     irc = baseio_open(iunsave,'minuit.save',stat='UNKNOWN')
     ! this sets minuit to read from input file
     !-----------------------------------------
     irc = minuit_init(iunin,iunout,iunsave)
     mininp = csmmod_setrdminuit(iunin)
     if(mininp.lt.0) then
        call message(MESERRO,'[CSMG] unable to open minuit section.')
        call message(MESKILL,'[CSMG] kill request.')
     end if
     ! here we start minuit
     !---------------------
     call minuit_minuit
  !
  ! Option Start Interctive Minuit Session
  !---------------------------------------
  else if(csm_interminuit) then
     write(iunout,'("# Request for Minuit: Read from stdin.")') 
     ! this sets minuit to read from stdin
     !------------------------------------
     mininp = csmmod_setrdminuit(iunin)
     if(mininp.lt.0) then
        call message(MESERRO,'[CSMG] unable to open minuit section.')
        call message(MESKILL,'[CSMG] kill request.')
     end if
     call minuit_minuit
     write(iunout,'("# MINUIT CSM VALUE in % is:",1X,F18.11)')&
          &modcsmmin_lastfunctionvalue*100.0_FREAL
  !
  !
  ! Option NOT Minuit: Run simplex from C.Daul code
  !------------------------------------------------
  else if(.not.csm_useminuit) then
     write(iunout,'("# Request Simplex: Run it with internal procedure. NOT IMPLEMENTED!")') 
     curparmin(:) = 0.0_FREAL
     parsteps(1:3) = 0.1_FREAL
     parsteps(4:6) = 0.1_FREAL
     ftol = 0.00001_FREAL
     resid = CSMI
     !
     write(iunout,'("# Request for Simplex: Initial Coordinates.")') 
     call molecule_print(iunout,'mldfromau')
     if(.not.csm_onestep) then
        call modcsm_nmopt(curparmin(1:6),parsteps(1:6),ftol,resid,iter)
        write(iunout,'("# CSM Simplex Values from input minimization.")')
        write(iunout,'("# Final Number of Simplex optim. Steps:",1X,I5)') iter
     end if
     CSMmain = resid
     write(iunout,'("# Final CSM Simplex Value is:",1X,F18.11)') CSMmain
     call molecule_print(iunout,'mldfromau')
     !
     ! restore old coordinates: we start over
     !---------------------------------------
     irc = modcsmmin_restore()
     !
     !
     ! We rotate on x
     !---------------
     curparmin(:) = 0.0_FREAL
     !
     ra = 90.0_FREAL
     rb =  0.0_FREAL
     call csmmod_rotxyz(0.0_FREAL,0.0_FREAL,0.0_FREAL,ra,rb,rb,molecule_xyz)
     !
     write(iunout,'("# Coordinates rotated along X. New coordinates:")')
     call molecule_print(iunout,'mldfromau')
     call csmmod_getimpj(rrimpj)
     CSMx = ((real(lnsymop,FREAL)-1.0_FREAL)/(real(lnsymop,FREAL))) - ((sum(rrimpj(:,:,:)))/(real(lnsymop,FREAL)*srhorho))
     if(.not.csm_onestep) then
        call modcsm_nmopt(curparmin(1:6),parsteps(1:6),ftol,resid,iter)
        write(iunout,'("# CSM Values after x minimization.")')
        write(iunout,'("# Final Number of optim. Steps:",1X,I5)') iter
        write(iunout,'("#      Before min is:",1X,F18.11)') CSMx
        CSMx = resid
     else
        write(iunout,'("# CSM Values after x rotation.")')
     end if
     write(iunout,'("# Final CSM Value is:",1X,F18.11)') CSMx
     irc = modcsmmin_restore()
     !
     ! We rotate on y
     !---------------
     curparmin(:) = 0.0_FREAL
     ra = 90.0_FREAL
     rb =  0.0_FREAL
     call csmmod_rotxyz(0.0_FREAL,0.0_FREAL,0.0_FREAL,rb,ra,rb,molecule_xyz)
     !
     write(iunout,'("# Coordinates rotated along Y.")')
     call molecule_print(iunout,'mldfromau')
     call csmmod_getimpj(rrimpj)
     CSMy = ((real(lnsymop,FREAL)-1.0_FREAL)/(real(lnsymop,FREAL))) - ((sum(rrimpj(:,:,:)))/(real(lnsymop,FREAL)*srhorho))
     if(.not.csm_onestep) then
        call modcsm_nmopt(curparmin(1:6),parsteps(1:6),ftol,resid,iter)
        write(iunout,'("# CSM Values after y minimization.")')
        write(iunout,'("# Final Number of optim. Steps:",1X,I5)') iter
        write(iunout,'("#      Before min is:",1X,F18.11)') CSMy
        CSMy = resid
     else
        write(iunout,'("# CSM Values after y rotation.")')
     end if
     write(iunout,'("# Final CSM Value is:",1X,F18.11)') CSMy
     irc = modcsmmin_restore()
     !call molecule_print(77,'pxyz')
     !
     ! We rotate on z
     !---------------
     curparmin(:) = 0.0_FREAL
     ra = 90.0_FREAL
     rb =  0.0_FREAL
     call csmmod_rotxyz(&
          &0.0_FREAL,0.0_FREAL,0.0_FREAL,&
          &rb,rb,ra,molecule_xyz)
     !irc = modcsmmin_setcoord()
     !call molecule_print(77,'pxyz')
     !
     write(iunout,'("# Coordinates rotated along Z.")')
     call molecule_print(iunout,'mldfromau')
     call csmmod_getimpj(rrimpj)
     CSMz = ((real(lnsymop,FREAL)-1.0_FREAL)/(real(lnsymop,FREAL))) - ((sum(rrimpj(:,:,:)))/(real(lnsymop,FREAL)*srhorho))
     if(.not.csm_onestep) then
        call modcsm_nmopt(curparmin(1:6),parsteps(1:6),ftol,resid,iter)
        write(iunout,'("# CSM Values after z minimization.")')
        write(iunout,'("# Final Number of optim. Steps:",1X,I5)') iter
        write(iunout,'("#      Before min is:",1X,F18.11)') CSMz
        CSMz = resid
     else
        write(iunout,'("# CSM Values after z rotation.")')
     end if
     write(iunout,'("# Final CSM Simplex Value is:",1X,F18.11)') CSMz
     irc = modcsmmin_restore()
     !
     ! Write out summary
     !------------------
     write(iunout,'("# CSM Summary.")')
     write(iunout,'("#    Inizial CSM Value is:",1X,F18.11)') CSMI
     write(iunout,'("#   After optimization is:",1X,F18.11)') CSMmain
     write(iunout,'("#          Along x axe is:",1X,F18.11)') CSMx
     write(iunout,'("#          Along y axe is:",1X,F18.11)') CSMy
     write(iunout,'("#          Along z axe is:",1X,F18.11)') CSMz
     CSMmin = min(CSMI,CSMmain,CSMx,CSMy,CSMz)
     write(iunout,'("# The Lowest CSM Value is:",1X,F18.11)') CSMmin
     write(iunout,'("#    Obtained CSM in % is:",1X,F18.11)') CSMmin*100.0_FREAL
     !
     !
     !
     !
  else
     write(iunout,'("# Request for Minuit: Run it with internal procedure.")') 
     !
     ! Here start minuit
     !------------------
     write(iunout,'("# WE RUN NOW MINUIT")')
     !
     ! Reset of parameters and Init
     !-----------------------------
     curparmin(:) = 0.0_FREAL
     parsteps(1:3) = 0.1_FREAL
     parsteps(4:6) = 0.1_FREAL
     parvect(:) = (/1,2,3,4,5,6/)
     parlabel(:) =  (/'disp x','disp y','disp z','rot x ','rot y ','rot z '/)
     !
     irc = bio_getfreefh(iunsave)
     irc = minuit_init(iunin,iunout,iunsave)
     !
     ! Initialize minuit parameters
     !-----------------------------
     do i= 1,6
        call minuit_mnparm(parvect(i),parlabel(i),curparmin(i),parsteps(i),&
             &0.0_FREAL,0.0_FREAL,ierflg) 
        if (ierflg .ne. 0)  then 
           call message_value(MESERRO,'[CSMG]: unable to define parameter no.',i)
           call message(MESKILL,'[CSMG]: we stop now')
        endif
     end do
     call minuit_mnseti('CSM Minuit') 
     !
     ! Be carefull this parameter is important
     !----------------------------------------
     call minuit_mnexcm('set eps  ',(/1.0E-6_FREAL/),1,ierflg)
     !
     call minuit_mnexcm('call fcn ',(/1.0_FREAL/),1,ierflg)
     call minuit_mnexcm('set print',(/3.0_FREAL/),1,ierflg)
     !
     ! We fix all the parameter first
     !-------------------------------
     call minuit_mnexcm('fix      ',(/1.0_FREAL/),1,ierflg)
     call minuit_mnexcm('fix      ',(/2.0_FREAL/),1,ierflg)
     call minuit_mnexcm('fix      ',(/3.0_FREAL/),1,ierflg)
     call minuit_mnexcm('fix      ',(/4.0_FREAL/),1,ierflg)
     call minuit_mnexcm('fix      ',(/5.0_FREAL/),1,ierflg)
     call minuit_mnexcm('fix      ',(/6.0_FREAL/),1,ierflg)
     !
     ! We now release the rotations
     !-----------------------------
     call minuit_mnexcm('release  ',(/4.0_FREAL/),1,ierflg)
     call minuit_mnexcm('migrad   ',(/0.0_FREAL/),0,ierflg)
     call minuit_mnexcm('release  ',(/5.0_FREAL/),1,ierflg)
     call minuit_mnexcm('migrad   ',(/0.0_FREAL/),0,ierflg)
     call minuit_mnexcm('release  ',(/6.0_FREAL/),1,ierflg)
     call minuit_mnexcm('migrad   ',(/0.0_FREAL/),0,ierflg)
     call minuit_mnexcm('minos    ',(/0.0_FREAL/),0,ierflg)
     call minuit_mnexcm('call fcn ',(/3.0_FREAL/),1,ierflg)
     !
     ! We release one by one the displacements
     !----------------------------------------
     call minuit_mnexcm('release  ',(/1.0_FREAL/),1,ierflg)
     call minuit_mnexcm('migrad   ',(/0.0_FREAL/),0,ierflg)
     call minuit_mnexcm('release  ',(/2.0_FREAL/),1,ierflg)
     call minuit_mnexcm('migrad   ',(/0.0_FREAL/),0,ierflg)
     call minuit_mnexcm('release  ',(/3.0_FREAL/),1,ierflg)
     call minuit_mnexcm('migrad   ',(/0.0_FREAL/),0,ierflg)
     call minuit_mnexcm('minos    ',(/0.0_FREAL/),0,ierflg)
     !
     call minuit_getamin(modcsmmin_lastfunctionvalue)
     write(iunout,'("# Minuit CSM initial in % is:",1X,F18.11)')&
          &modcsmmin_lastfunctionvalue*100.0_FREAL
     !
     csm_save_minuit = modcsmmin_lastfunctionvalue
     !
     ! We save the best parameters so far
     !-----------------------------------
     !
     ! We now start it again from the other directions
     ! We assume that displacment is quite ok now.
     ! And we just rotate by 90 every angle....
     !-----------------------------------------------
     ! Save current params
     ! Each iteration start from the original position
     !------------------------------------------------
     do i= 1,6
        call minuit_mnpout(parvect(i),parlabel(i),csm_pars_minuit(i),csm_error_minuit,&
             &csm_dummy_minuit,csm_dummy_minuit,ierflg) 
        if (ierflg .lt. 0)  then 
           call message_value(MESERRO,'[CSMG]: unable to define parameter no.',i)
           call message(MESKILL,'[CSMG]: we stop now')
        endif
     end do
     ! X axis first
     !-------------
     !call minuit_mnpout(4,chardummy,csm_svpar_minuit,&
     !     &csm_error_minuit,csm_dummy_minuit,csm_dummy_minuit,csm_ipar_minuit)
     !
     ! We set it (Xrot=4)
     !------------------
     csm_temp_minuit(1) = 4.0_FREAL
     csm_temp_minuit(2) = csm_pars_minuit(4) + (LPI/2.0_FREAL)
     call minuit_mnexcm('set par  ',csm_temp_minuit(1:2),2,ierflg)
     call minuit_mnexcm('migrad   ',(/0.0_FREAL/),0,ierflg)
     !
     call minuit_getamin(modcsmmin_lastfunctionvalue)
     if(csm_save_minuit.gt.modcsmmin_lastfunctionvalue) then
        !
        ! It is better... so we save the value
        !-------------------------------------
        csm_save_minuit = modcsmmin_lastfunctionvalue
        !
        ! We save the best parameters so far
        !-----------------------------------
        do i= 1,6
           call minuit_mnpout(parvect(i),parlabel(i),csm_svbestpars(i),csm_error_minuit,&
                &csm_dummy_minuit,csm_dummy_minuit,ierflg) 
           if (ierflg .lt. 0)  then 
              call message_value(MESERRO,'[CSMG]: unable to define parameter no.',i)
              call message(MESKILL,'[CSMG]: we stop now')
           endif
        end do
     end if
     write(iunout,'("# Minuit CSM on X rot. in % is:",1X,F18.11)') csm_save_minuit*100.0_FREAL
     ! We restore the parameters
     !------------------------------------------------
     do i= 1,6
        call minuit_mnparm(parvect(i),parlabel(i),csm_pars_minuit(i),parsteps(i),&
             &0.0_FREAL,0.0_FREAL,ierflg) 
        if (ierflg .ne. 0)  then 
           call message_value(MESERRO,'[CSMG]: unable to define parameter no.',i)
           call message(MESKILL,'[CSMG]: we stop now')
        endif
     end do
     !
     ! Y axis second
     !--------------
     !
     ! We set it (Yrot=5)
     !------------------
     csm_temp_minuit(1) = 5.0_FREAL
     csm_temp_minuit(2) = csm_pars_minuit(5) + (LPI/2.0_FREAL)
     call minuit_mnexcm('set par  ',csm_temp_minuit(1:2),2,ierflg)
     call minuit_mnexcm('migrad   ',(/0.0_FREAL/),0,ierflg)
     !
     call minuit_getamin(modcsmmin_lastfunctionvalue)
     if(csm_save_minuit.gt.modcsmmin_lastfunctionvalue) then
        !
        ! It is better so....we save the value
        !-------------------------------------
        csm_save_minuit = modcsmmin_lastfunctionvalue
        !
        ! We save the best parameters so far
        !-----------------------------------
        do i= 1,6
           call minuit_mnpout(parvect(i),parlabel(i),csm_svbestpars(i),csm_error_minuit,&
                &csm_dummy_minuit,csm_dummy_minuit,ierflg) 
           if (ierflg .lt. 0)  then 
              call message_value(MESERRO,'[CSMG]: unable to define parameter no.',i)
              call message(MESKILL,'[CSMG]: we stop now')
           endif
        end do
     end if
     write(iunout,'("# Minuit CSM on Y rot. in % is:",1X,F18.11)') csm_save_minuit*100.0_FREAL
     ! We restore the parameters
     !------------------------------------------------
     do i= 1,6
        call minuit_mnparm(parvect(i),parlabel(i),csm_pars_minuit(i),parsteps(i),&
             &0.0_FREAL,0.0_FREAL,ierflg) 
        if (ierflg .ne. 0)  then 
           call message_value(MESERRO,'[CSMG]: unable to define parameter no.',i)
           call message(MESKILL,'[CSMG]: we stop now')
        endif
     end do
     !
     ! Z axis 3rd
     !-----------
     !
     ! We set it (Zrot=6)
     !------------------
     csm_temp_minuit(1) = 6.0_FREAL
     csm_temp_minuit(2) = csm_pars_minuit(6) + (LPI/2.0_FREAL)
     call minuit_mnexcm('set par  ',csm_temp_minuit(1:2),2,ierflg)
     call minuit_mnexcm('migrad   ',(/0.0_FREAL/),0,ierflg)
     !
     call minuit_getamin(modcsmmin_lastfunctionvalue)
     if(csm_save_minuit.gt.modcsmmin_lastfunctionvalue) then
        !
        ! It is better so....we save the value
        !-------------------------------------
        csm_save_minuit = modcsmmin_lastfunctionvalue
        !
        ! We save the best parameters so far
        !-----------------------------------
        do i= 1,6
           call minuit_mnpout(parvect(i),parlabel(i),csm_svbestpars(i),csm_error_minuit,&
                &csm_dummy_minuit,csm_dummy_minuit,ierflg) 
           if (ierflg .lt. 0)  then 
              call message_value(MESERRO,'[CSMG]: unable to define parameter no.',i)
              call message(MESKILL,'[CSMG]: we stop now')
           endif
        end do
     end if
     write(iunout,'("# Minuit CSM on Z rot. in % is:",1X,F18.11)') csm_save_minuit*100.0_FREAL
     ! We restore the parameters
     !------------------------------------------------
     do i= 1,6
        call minuit_mnparm(parvect(i),parlabel(i),csm_pars_minuit(i),parsteps(i),&
             &0.0_FREAL,0.0_FREAL,ierflg) 
        if (ierflg .ne. 0)  then 
           call message_value(MESERRO,'[CSMG]: unable to define parameter no.',i)
           call message(MESKILL,'[CSMG]: we stop now')
        endif
     end do
     !
     ! We should have got the best so far
     !-----------------------------------
     ! We set the best parameters
     !---------------------------
     do i= 1,6
        call minuit_mnparm(parvect(i),parlabel(i),csm_svbestpars(i),parsteps(i),&
             &0.0_FREAL,0.0_FREAL,ierflg) 
        if (ierflg .ne. 0)  then 
           call message_value(MESERRO,'[CSMG]: unable to define parameter no.',i)
           call message(MESKILL,'[CSMG]: we stop now')
        endif
     end do
     !
     !We run migrad on strategy 2
     !---------------------------
     call minuit_mnexcm('set str  ',(/2.0_FREAL/),1,ierflg)
     call minuit_mnexcm('migrad   ',(/0.0_FREAL/),0,ierflg)
     call minuit_getamin(modcsmmin_lastfunctionvalue)
     if(csm_save_minuit.gt.modcsmmin_lastfunctionvalue) then
        !
        ! It is better so....we save the value
        !-------------------------------------
        csm_save_minuit = modcsmmin_lastfunctionvalue
     end if
     !
     ! And print out
     !--------------
     write(iunout,'("# MINUIT BEST CSM VALUE in % is:",1X,F18.11)')&
          &csm_save_minuit*100.0_FREAL

     !
  end if
  !
  !
contains
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H Function CSMG_init()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function CSMG_init()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the CSMG module.
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    !
    LocCSMG_isinit = .true.
    CSMG_init = 1
    !
    !
  end function CSMG_init
!
!H
!H-----------------------------------------------------------------------------
!H logical function CSMG_isinit()
!H-----------------------------------------------------------------------------
!H
!
  logical function CSMG_isinit()
!
!H
!H-----------------------------------------------------------------------------
!H We merely return the initialization flag.
!H-----------------------------------------------------------------------------
!H
!
    !
    CSMG_isinit = LocCSMG_isinit
    !
  end function CSMG_isinit
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function CSMG_geterror_i()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function CSMG_geterror_i()
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
    CSMG_geterror_i = LocCSMG_error
    !
  end function CSMG_geterror_i
!
!H-----------------------------------------------------------------------------
!H
!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H character(FLCHARS) function CSMG_geterror_c()
!H-----------------------------------------------------------------------------
!H
!
  character(FLCHARS) function CSMG_geterror_c(code)
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
       CSMG_geterror_c = "No Errors."
    case(1_FINT)
       CSMG_geterror_c = "Not documented error or generic error.."
    case(2_FINT)
       CSMG_geterror_c = "Double Error: Why did you get this message?."
    case default
       CSMG_geterror_c = "Not documented error or generic error.."
    end select
    !
    return
    !
  end function CSMG_geterror_c
!
!H
!H-----------------------------------------------------------------------------
!H
!

  !
  !
end program CSMG
