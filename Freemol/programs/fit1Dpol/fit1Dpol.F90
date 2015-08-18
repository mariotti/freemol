!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H PROGRAM fit1Dpol  Frimol by F.Mariotti:
!H-----------------------------------------------------------------------------
!H $Id: fit1Dpol.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
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
!H ADDING A NEW FUNCTION HELP: 
!H 
!H 1) Write the code
!H    
!H    - open the file fit1Dpol_commons.F90
!H    - add you function as integer identifier of the type
!H      FUN_XXXXX=Y where Y is a free handle
!H      just after the other definitions: search for "FUN_DUMMY=0".
!H    
!H    - open the file f1dp_minuit_function.F90
!H      Search for the subroutine minuit_fcn ("subroutine minuit_fcn")
!H      Add a case in the select tree of this function: Last case is "case default"
!H      just before it add a case statment like this one:
!H    
!H    case(FUN_XXXXX)
!H       irc=funct_XXXXX(npar,grad,fval,xval,iflag)
!H       return
!H    
!H    The actual arguments to the function can change with versions please check it!
!H    These have to be the same as for "subroutine minuit_fcn".
!H    
!H    - Add your actual function code.
!H      The best idea is to copy an existing one like for example "funct_polyn"
!H      and even a best idea is to call your function "funct_XXXXX"!
!H    
!H    - You have to make sure that your function gets readed in the input file
!H      the function f1dp_minuit_function_init sets it. So search for
!H      this function
!H      ("function f1dp_minuit_function_init"). At present the only thing this
!H      init
!H      does is to set the correct parameter once it has the function name
!H      so don't
!H      get lost and add something like this just before the last (lonly) else:
!H    
!H    else if(index(f1dp_funame,'XXXXX').gt.0) then
!H       f1dp_flags(FLG_FUNC)=FUN_XXXXX
!H       nfpars=f1dp_numpar
!H    
!H    Note that we set the actual number of parameters just in case
!H    in principle you should WARN about that!
!H    
!H    - You are done!
!H    
!H    - PLEASE TEST IT!
!H    
!H    - THINGS TO PAY ATTENTION AND REMINDER
!H    
!H    Data are in f1dp_data(f1dp_numrec,f1dp_numcol)
!H    with dimensions: f1dp_numrec,f1dp_numcol
!H    
!H    Parameters are in xval(f1dp_numpar) if you used very same variable names 
!H    
!H    WARN: The npar input variable is the actual number of free parameters! WARN
!H    
!H    For iflag read MINUIT manual but in general: iflag=
!H    1 called once at the beginning for any setup/init (note data are already stored)
!H    2 called to obtain derivatives: if unsure reset them and minuit will check
!H      versus numerical ones.
!H    3 called at end and should print the results
!H    
!H    The function should always return it's own value (in fval if very same vars)
!H    
!H    The function is an integer function and should return 0 if no error or
!H    return a negative number if error (eventual value can set the error code)
!H    
!H    Order of parameters is your own buisness i.e. it is the very same as imput
!H    but you can eventually use the name of the parameters stored in 
!H    f1dp_parnames(f1dp_numpar). But you will slow down the code!
!H    
!H    You can use all the fit1Dpol_commons.F90 data but remember that it is Minuit
!H    that calls you!
!H    
!H    
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H 
!
!
! f1dp_flags ALLOCATED FLAGS
!---------------------------
!
! f1dp_flags(1) = order of the polinomial. Given by -p N
! f1dp_flags(2) = Fit function to be used: they have to be implemented
!
!---------------------------
!
program fit1Dpol
  !
  ! When including modules pay attention to lib ordering
  !
  use vartypes        ! require "-lincludes" in USEDLIBS var of Makelocal
  !
  ! Used modules
  use messages        ! require "-lincludes" in USEDLIBS var of Makelocal
  use baseio          ! require "-lmodules" in USEDLIBS var of Makelocal
  use pcmdline        ! require "-lincludes" in USEDLIBS var of Makelocal
  use osec            ! require "-lmodules" in USEDLIBS var of Makelocal
  use linetools       ! require "-lincludes" in USEDLIBS var of Makelocal
  use fit1Dpol_commons! all the common data (local module)
  use fit1dpol_minuit ! it will be present if the normal configure procedure if followed
  use f1dp_minuit_function ! We need it for some initialization staff
  !
  implicit none
  !
  !
  ! Procedures present in this program
  !-----------------------------------
  ! fit1Dpol_init
  ! fit1Dpol_isinit
  ! fit1Dpol_geterror
  ! fit1Dpol_
  !
  !
  ! DECLARATIONS:
  !--------------------
  !
  !
  ! init and error variables
  !-------------------------
  logical, save :: LOCfit1Dpol_isinit
  integer(FINT), save :: LOCfit1Dpol_error
  character(FLCHARS), save :: LOCfit1Dpol_error_message
  !
  ! Local Dummy Vars
  !-----------------
  integer(FINT) :: i
  integer(FINT) :: ierflg
  integer(FINT) :: irc
  integer(FINT) :: nfpars
  !
  ! Local Vars
  !-----------------
  ! for eventual minuit use
  integer(FINT) :: iunsave
  ! NOTE: we choosed in this case to put eventual global flags in a flags vector
  !
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Here we start
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  !
  ! Reset common data
  irc=fit1Dpol_commons_reset()
  if(irc.lt.0) then
     LOCfit1Dpol_error = 1
     call message(MESERRO,"[fit1Dpol] Cannot reset common data.")
     call message(MESKILL,"[fit1Dpol] We exit now.")
  end if
  !
  ! Get messages initialization with default parameters first
  ! it will be changed later in the code if needed
  !----------------------------------------------------------
  call messages_init()
  !
  ! File I/O stuff initialization
  !------------------------------
  irc = baseio_init()
  if(irc.lt.0) then
     LOCfit1Dpol_error = 1
     call message(MESERRO,"[fit1Dpol] Cannot Initialize BaseIO Module.")
     call message(MESKILL,"[fit1Dpol] We exit now.")
  end if
  !
  ! Just follow the rules and call this program init
  !-------------------------------------------------
  irc = fit1Dpol_init(f1dp_flags)
  if(irc.lt.0) then
     LOCfit1Dpol_error = 1
     call message(MESERRO,"[fit1Dpol] Cannot Initialize.")
     call message(MESKILL,"[fit1Dpol] We exit now.")
  end if
  !
  ! We reinitialize messages: we have new units
  !--------------------------------------------
  ! We actually use defaults: everything to output except error duplicated on stderr
  call messages_init(f1dp_iunout,f1dp_iunout,0,f1dp_iunout)
  !
  ! Now we read the input file
  !---------------------------
  !irc = fit1Dpol_rdinput(f1dp_iunin,f1dp_iunout,f1dp_flags)
  irc = fit1Dpol_rdinput(f1dp_iunin,f1dp_iunout)
  if(irc.lt.0) then
     LOCfit1Dpol_error = 1
     call message(MESERRO,"[fit1Dpol] Cannot read the input file.")
     call message(MESKILL,"[fit1Dpol] We exit now.")
  end if
  !
  ! We now initialize minuit
  !-------------------------
  !
  ! we first setup the eventual name for minuit session to eventually anable restart we use temporarly 'i'
  if(len_trim(f1dp_sessionname).gt.0) then
     f1dp_sessionname=trim(f1dp_sessionname) ! just if we didn't do it yet!
     i = len_trim(f1dp_sessionname)
     f1dp_minuit_savename = f1dp_minuit_savename(1:i)//'_minuit.save'
     f1dp_minuit_savename = trim(f1dp_minuit_savename)
  else
     f1dp_minuit_savename = 'minuit.save'
     f1dp_minuit_savename = trim(f1dp_minuit_savename)
  end if
  i = len_trim(f1dp_minuit_savename)
  ! we are forced here to give an unknown status
  i = baseio_open(iunsave,f1dp_minuit_savename(1:i),stat='UNKNOWN')
  !
  ! We initialize other function specific staff
  !--------------------------------------------
  ! In this case we have a modified version of the function module
  ! which has an extra fuction to actually perform the reset and fill in of the
  ! data because it could depend from the function type and specific options.
  ! This call should modify data and return the dimesion in such a way that we
  ! can call the minuit_mnparm in a do loop of this type
  !
  ! do i=1,nfpars
  !    call minuit_mnparm(parvect(i),parlabel(i),curparmin(i),parsteps(i),0.0_FREAL,0.0_FREAL,ierflg)
  ! end do
  !
  ! which will fill in the whole minimization process (see minuit manual)
  !
  !----------------------------------------------------------------------------
  irc = f1dp_minuit_function_init(nfpars)
  !
  irc = minuit_init(f1dp_iunin,f1dp_iunout,iunsave)
  !
  ! Fill in minuit with data
  !-------------------------
  do i=1,nfpars
     call minuit_mnparm(f1dp_idxpars(i),f1dp_parnames(i),f1dp_pars(i,1),f1dp_pars(i,2),0.0_FREAL,0.0_FREAL,ierflg)
     if (ierflg.ne.0)  then
        call message_value(MESERRO,'[fit1Dpol]: unable to define parameter no.',i)
        call message(MESKILL,'[fit1Dpol]: we stop now: problem in minuit subsystem.')
     endif
  end do
  !
  ! The -m 2 option has been set: minuit is interactive and it stop just after
  !---------------------------------------------------------------------------
  if (f1dp_flags(FLG_MINT).eq.2) then
     irc = minuit_init(5,6,iunsave)
     call minuit_minuit
  end if
  !
  ! Set the title of the session
  !-----------------------------
  if(len_trim(f1dp_sessionname).gt.0) then
     call minuit_mnseti(f1dp_sessionname)
  else
     call minuit_mnseti('Minuit Session XX')
  end if
  !
  ! We call the function initialization
  !------------------------------------
  call minuit_mnexcm('call fcn ',(/1.0_FREAL/),1,ierflg)
  !
  ! Fix what has to be fixed
  !-------------------------
  do i=1,nfpars
     if(f1dp_fixpars(i).eq.1) then
        if(f1dp_debug.gt.0) then
           call message_value(MESLOG,"[fit1Dpol] Fixing parameter no: ",i)
        end if
        call minuit_mnexcm('fix      ',(/real(i,FREAL)/),1,ierflg)
     end if
  end do
  !
  ! We assume we have gradients at present WARNING
  call minuit_mnexcm('set grad ',(/0.0_FREAL/),1,ierflg)
  !
  ! set at present an intermediate print
  call minuit_mnexcm('set pri  ',(/2.0_FREAL/),1,ierflg)
  ! 
  ! set accurate strategy
  call minuit_mnexcm('set str  ',(/2.0_FREAL/),1,ierflg)
  !
  ! We setup the EPS for the fit if present in the input
  !-----------------------------------------------------------------------------
  if(f1dp_flags(FLG_IEPS).eq.1) then
     call minuit_mnexcm('set eps  ',(/f1dp_eps/),1,ierflg)
  end if
  if (ierflg.ne.0)  then
     call message(MESERRO,'[fit1Dpol]: unable to setup the EPS parameter')
     call message(MESKILL,'[fit1Dpol]: we stop now: problem in minuit subsystem.')
  endif
  !
  ! be sure that debug level gt.1 get output
  !-----------------------------------------
  if(f1dp_debug.gt.1) then
     call flush(f1dp_iunout)
  end if
  !
  !
  ! Call first time function to initialize it: NOTE data are there already
  !-----------------------------------------------------------------------
  call minuit_mnexcm('call fcn ',(/1.0_FREAL/),1,ierflg)
  !
  ! be sure that debug level gt.1 get output
  !-----------------------------------------
  if(f1dp_debug.gt.1) then
     call flush(f1dp_iunout)
  end if
  !
  ! Minuit run
  !-------------------
  !call minuit_mnexcm('simplex  ',(/0.0_FREAL/),0,ierflg)
  call minuit_mnexcm('migrad   ',(/0.0_FREAL/),0,ierflg)
  !
  ! be sure that debug level gt.1 get output
  !-----------------------------------------
  if(f1dp_debug.gt.1) then
     call flush(f1dp_iunout)
  end if
  !
  call minuit_mnexcm('minos    ',(/0.0_FREAL/),0,ierflg)
  !
  ! be sure that debug level gt.1 get output
  !-----------------------------------------
  if(f1dp_debug.gt.1) then
     call flush(f1dp_iunout)
  end if
  !
  !call minuit_mnexcm('imp      ',(/0.0_FREAL/),0,ierflg)
  call minuit_mnexcm('migrad   ',(/1000.0_FREAL,1.0D-08/),2,ierflg)
  !
  ! be sure that debug level gt.1 get output
  !-----------------------------------------
  if(f1dp_debug.gt.1) then
     call flush(f1dp_iunout)
  end if
  !
  !
  ! Minuit ask function to terminate
  !---------------------------------
  call minuit_mnexcm('call fcn ',(/3.0_FREAL/),1,ierflg)
  !
  ! Warn for debugging stage
  call message(MESOUT,"Program not completed: mesout test")
  call message(MESERRO,"Program not completed: meserro test")
  call message(MESLOG,"Program not completed: meslog test")

contains
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H Function fit1Dpol_init(f1dp_flags)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function fit1Dpol_init(f1dp_flags)
    !
    !integer(FINT), intent(out) :: f1dp_iunin,f1dp_iunout
    integer(FINT), intent(out), dimension(:) :: f1dp_flags
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the fit1Dpol module.
!H
!H-----------------------------------------------------------------------------
!H
!
    ! Local declarations
    !-------------------
    integer(FINT) :: nargs
    character(FLCHARS),dimension(:), allocatable :: cargs
    character(FLCHARS) :: finput
    character(FLCHARS) :: foutput
    character(FLCHARS) :: params
    !
    integer(FINT) :: iostat
    integer(FINT) :: itmp
    character(FLCHARS) :: rfrmt
    !
    LOCfit1Dpol_isinit = .true.
    fit1Dpol_init = 0
    !
    ! Some defaults
    f1dp_iunin=5
    f1dp_iunout=6
    !
    ! Use the command line facility
    !------------------------------
    call pcmd_iargc(nargs)
    allocate(cargs(nargs),STAT=irc)
    if(.not.irc.eq.0) then
       LOCfit1Dpol_error = 2
       LOCfit1Dpol_isinit = .false.
       fit1Dpol_init = -1
       call message(MESDEBG,'[fit1Dpol_init] Cannot allocate cmd line parameters')
       return
    end if
    call pcmd_getarg(nargs,cargs)
    !
    ! Get the input and out put files if specified with -i and -o
    !------------------------------------------------------------
    call pcmd_getio(nargs,cargs,finput,foutput)
    !
    ! We open the input and output files
    !-----------------------------------
    finput = trim(finput)
    if(scan(finput," ").eq.0) then
       f1dp_iunin = 5
       finput = "stdin"
    else
       irc = baseio_open(f1dp_iunin,finput,stat='OLD')
       if(irc.lt.0) then
          LOCfit1Dpol_error = 3
          LOCfit1Dpol_isinit = .false.
          fit1Dpol_init = -1
          call message(MESERRO,"[fit1Dpol_init] Cannot Open Input File")
          call message_value(MESLOG,"[fit1Dpol_init] Input file name: ",finput)
          return
       end if
    end if
    if(f1dp_debug.gt.1) then
       call message_value(MESLOG,"[fit1Dpol_init] Input file name: ",finput)
       call message_value(MESLOG,"[fit1Dpol_init] Input file unit: ",f1dp_iunin)
    end if
    !
    foutput = trim(foutput)
    if(scan(foutput," ").eq.0) then
       f1dp_iunout = 6
       foutput = "stdout"
    else
       irc = baseio_open(f1dp_iunout,foutput,stat='UNKNOWN')
       if(irc.lt.0) then
          LOCfit1Dpol_error = 3
          LOCfit1Dpol_isinit = .false.
          fit1Dpol_init = -1
          call message(MESERRO,"[fit1Dpol_init] Cannot Open Output File")
          call message_value(MESLOG,"[fit1Dpol_init] Output file name: ",foutput)
          return
       end if
    end if
    if(f1dp_debug.gt.1) then
       call message_value(MESLOG,"[fit1Dpol_init] Output file name: ",foutput)
       call message_value(MESLOG,"[fit1Dpol_init] Output file unit: ",f1dp_iunout)
    end if
    !
    ! Other options
    !--------------
    !
    ! -p N gets the order of polynom if you need to reduce it
    !--------------------------------------------------------
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-p',params)
    params=trim(params)
    if(irc.gt.0) then
       call message(MESERRO,"[fit1Dpol_init] -p option is not implemented!")
       call message(MESERRO,"[fit1Dpol_init] Check if new versions do!")
       call message(MESKILL,"[fit1Dpol_init] We stop now.")
       if(len_trim(params).lt.1) then
          LOCfit1Dpol_error = 1
          LOCfit1Dpol_isinit = .false.
          fit1Dpol_init = -1
          call message(MESERRO,"[fit1Dpol_init] -p option require an integer argument")
          return
       end if
       read(params,'(I8)', IOSTAT=iostat) itmp
       if(iostat.ne.0) then
          LOCfit1Dpol_error = 1
          LOCfit1Dpol_isinit = .false.
          fit1Dpol_init = -1
          call message(MESERRO,"[fit1Dpol_init] -p option require an integer argument")
          return
       end if
       f1dp_flags(FLG_POLO)=itmp
    end if
    !
    ! -gnuplot <file>  specify to produce gnuplot file <file> is optional
    !--------------------------------------------------------------------
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-gnuplot',params)
    if(irc.gt.0) then
       f1dp_flags(FLG_GPLT)=1 ! 1 is produce gnuplot, 2 will be with file name
       params=adjustl(params)
       itmp=len_trim(params)
       if(itmp.gt.1) then
          f1dp_gpltprefix=adjustl(params(1:itmp))
          irc = baseio_open(f1dp_iungpltdata,f1dp_gpltprefix(1:itmp)//'.dat',stat='UNKNOWN')
          if(irc.lt.0) then
             call message(MESERRO,"[fit1Dpol_init] cannot open gnuplot data file")
             call message(MESKILL,"[fit1Dpol_init] We stop now")
          end if
          irc = baseio_open(f1dp_iungpltgplt,f1dp_gpltprefix(1:itmp)//'.gplt',stat='UNKNOWN')
          if(irc.lt.0) then
             call message(MESERRO,"[fit1Dpol_init] cannot open gnuplot file")
             call message(MESKILL,"[fit1Dpol_init] We stop now")
          end if
          f1dp_flags(FLG_GPLT)=2 ! gnuplot called with file names
       end if
       !
       if (f1dp_debug.gt.1) then
          call message(MESOUT,"[fit1Dpol_init] -gnuplot request will be set")
       end if
       !
    end if
    !
    !
    ! -d N debug option
    !--------------------------------------------------------------------------------------------
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-d',params)
    params=trim(params)
    if(irc.gt.0) then
       ! we have benn colled so we default at least to level 1
       call message(MESLOG,"[fit1Dpol_init] -d request: setting debug level")
       f1dp_debug=1
       if(len_trim(params).gt.1) then
          call line_getform(params,'I',rfrmt,i)
          rfrmt='('//rfrmt(1:i)//')'
          read (params,rfrmt(1:i+2)) irc
          if(irc.gt.0) then
             f1dp_debug=irc
             call message_value(MESLOG,"[fit1Dpol_init] debug level set to: ",irc)
          else
             call message(MESERRO,"[fit1Dpol_init] -d request with invalid level.")
          end if
       end if
       !
    end if
    !
    !
    !
    ! -m N 0 normal, 1 read (not implmented), 2 interactive minuit
    !--------------------------------------------------------------------------------------------
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-m',params)
    params=trim(params)
    if(irc.gt.0) then
       ! we have been called but we are still buggy so warn
       call message(MESLOG,"[fit1Dpol_init] -m request: Still development")
       f1dp_flags(FLG_MINT)=0
       if(len_trim(params).gt.1) then
          call line_getform(params,'I',rfrmt,i)
          rfrmt='('//rfrmt(1:i)//')'
          read (params,rfrmt(1:i+2)) irc
          if(irc.gt.-1.and.irc.lt.3) then
             f1dp_flags(FLG_MINT)=irc
             if(f1dp_debug.gt.0) then
                call message_value(MESLOG,"[fit1Dpol_init] Minuit goes with mode: ",irc)
             end if
          else
             call message(MESERRO,"[fit1Dpol_init] -m request with invalid level.")
             call message(MESKILL,"[fit1Dpol_init] We stop now.")
          end if
       end if
       !
    end if
    !
    ! -generate <file>  specify to produce generate file. <file> is optional and not implemented now
    !-----------------------------------------------------------------------------------------------
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-generate',params)
    if(irc.gt.0) then
       f1dp_flags(FLG_GNRT)=1 ! 1 is produce data in out, 2 will be with file name
       params=adjustl(params)
       itmp=len_trim(params)
       if(itmp.gt.1) then
          f1dp_gnrtprefix=adjustl(params(1:itmp))
          irc = baseio_open(f1dp_iungnrtdata,f1dp_gnrtprefix(1:itmp)//'.dat',stat='UNKNOWN')
          if(irc.lt.0) then
             call message(MESERRO,"[fit1Dpol_init] cannot open generate data file")
             call message(MESKILL,"[fit1Dpol_init] We stop now")
          end if
          irc = baseio_open(f1dp_iungnrtgnrt,f1dp_gnrtprefix(1:itmp)//'.gplt',stat='UNKNOWN')
          if(irc.lt.0) then
             call message(MESERRO,"[fit1Dpol_init] cannot open generate file")
             call message(MESKILL,"[fit1Dpol_init] We stop now")
          end if
          f1dp_flags(FLG_GNRT)=2 ! generate called with file names
       end if
       !
       if (f1dp_debug.gt.1) then
          call message(MESOUT,"[fit1Dpol_init] -generate request will be set")
       end if
       !
    end if
    !
  end function fit1Dpol_init
!
!H
!H-----------------------------------------------------------------------------
!H logical function fit1Dpol_isinit()
!H-----------------------------------------------------------------------------
!H
!
  logical function fit1Dpol_isinit()
!
!H
!H-----------------------------------------------------------------------------
!H We merelyreturn the initialization flag.
!H-----------------------------------------------------------------------------
!H
!
    !
    fit1Dpol_isinit = LOCfit1Dpol_isinit
    !
  end function fit1Dpol_isinit
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H character(FLCHARS) function fit1Dpol_geterror()
!H-----------------------------------------------------------------------------
!H
!
  character(FLCHARS) function fit1Dpol_geterror(code)
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
       fit1Dpol_geterror = "No Errors."
    case(1_FINT)
       fit1Dpol_geterror = "Not documented error or generic error.."
    case(2_FINT)
       fit1Dpol_geterror = "Unable to allocate for command line parameters."
    case(3_FINT)
       fit1Dpol_geterror = "Unable to open the input strim."
    case(4_FINT)
       fit1Dpol_geterror = "Double Error: Why did you get this message?."
    case(5_FINT)
       fit1Dpol_geterror = "Double Error: Why did you get this message?."
    case(6_FINT)
       fit1Dpol_geterror = "Double Error: Why did you get this message?."
    case(7_FINT)
       fit1Dpol_geterror = "Double Error: Why did you get this message?."
    case(8_FINT)
       fit1Dpol_geterror = "Double Error: Why did you get this message?."
    case default
       fit1Dpol_geterror = "Not documented error or generic error.."
    end select
    !
    return
    !
  end function fit1Dpol_geterror
!
!H
!H-----------------------------------------------------------------------------
!H
!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H Function fit1Dpol_rdinput(f1dp_iunin,f1dp_iunout,f1dp_flags)
!H-----------------------------------------------------------------------------
!H
!
!  integer(FINT) function fit1Dpol_rdinput(f1dp_iunin,f1dp_iunout,f1dp_flags)
  integer(FINT) function fit1Dpol_rdinput(f1dp_iunin,f1dp_iunout)
!
  integer(FINT), intent(in) :: f1dp_iunin,f1dp_iunout
!  integer(FINT), intent(in), dimension(:) :: f1dp_flags
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine read the input file
!H
!H-----------------------------------------------------------------------------
!H
!
  ! line_read required variables
  !-----------------------------
  ! we compute the required dimension
  ! in the following way:
  ! the maxim line lenght we will read is FLCHARS
  ! than if we separate words by blanks we have:
  ! for integer and characters the maximum number of words
  ! is LCHARS/2 while for reals is LCHARS/3
  ! and it is the worst case! (+1 just to be sure)
  !
  integer(FINT) :: niwrd
  integer(FINT) :: nrwrd
  integer(FINT) :: ncwrd
  integer(FINT), dimension(FLCHARS/2) :: iwrds
  real(FREAL), dimension((FLCHARS/3)+1) :: rwrds
  character(FLCHARS), dimension(FLCHARS/2) :: cwrds
  character(FLCHARS) :: lrform
  !
  !
  logical :: done
  character(FLCHARS) :: params, line
  character(FLCHARS) :: cfun, cmdim
  character(FLCHARS) :: cleps
  character(FLCHARS) :: cnumrec
  character(FLCHARS) :: cnumcol
  character(FLCHARS) :: cnumpars
  character(FLCHARS) :: rfrmt
  real(FREAL) :: leps
  integer(FINT) :: numrec
  integer(FINT) :: numcol
  integer(FINT) :: numpars
  integer(FINT) :: i,j,k
  !
  fit1Dpol_rdinput = 0
  done = .false.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! READING THE SESSION NAME SECTION
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  rewind(f1dp_iunin)
  params=''
  irc = osec_set(f1dp_iunin,'x-name',params)
  if(irc.gt.0) then
     irc = line_getline(f1dp_iunin,line,3)
     read(line,'(A)') f1dp_sessionname
     if(f1dp_debug.gt.0) then
        call message_value(MESOUT,"The session name given in the input is: ",f1dp_sessionname)
     end if
  end if
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! READING THE DATA SECTION
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  rewind(f1dp_iunin)
  params=''
  irc = osec_set(f1dp_iunin,'x-fit1dpol-data',params)
  if(irc.lt.0) then
     fit1Dpol_rdinput = -1
     LOCfit1Dpol_error = 1
     call message(MESERRO,'[fit1Dpol_rdinput] No x-fit1dpol-data section.')
     return
  end if
  !
  irc=line_getval(params,"numrec",cnumrec)
  if(irc.lt.0) then
     fit1Dpol_rdinput = -1
     LOCfit1Dpol_error = 1
     call message(MESERRO,'[fit1Dpol_rdinput] This version require numrec=N as parameter.')
     return
  end if
  call line_getform(cnumrec,'I',rfrmt,i)
  rfrmt='('//rfrmt(1:i)//')'
  read (cnumrec,rfrmt(1:i+2)) numrec
  !
  irc=line_getval(params,"numcol",cnumcol)
  if(irc.lt.0) then
     fit1Dpol_rdinput = -1
     LOCfit1Dpol_error = 1
     call message(MESERRO,'[fit1Dpol_rdinput] This version require numcol=N as parameter.')
     return
  end if
  call line_getform(cnumcol,'I',rfrmt,i)
  rfrmt='('//rfrmt(1:i)//')'
  read (cnumcol,rfrmt) numcol
  !
  !
  !irc=line_getval(params,"echo",ctmp)
  !if(irc.lt.0) then
  if(index(params,'echo').gt.0) then
     f1dp_flags(FLG_EINP)=1
     if(f1dp_debug.gt.1) then
        call message(MESLOG,'[fit1Dpol_rdinput] Echo of input request is set.')
     end if
     write(f1dp_iunout,*) "[x-fit1dpol-data] numrec=",numrec," numcol=",numcol," echo"
  end if
  call line_getform(cnumcol,'I',rfrmt,i)
  rfrmt='('//rfrmt(1:i)//')'
  read (cnumcol,rfrmt) numcol
  !
  f1dp_numrec = numrec
  f1dp_numcol = numcol
  if (f1dp_debug.gt.1) then
     call message_value(MESLOG,"[fit1Dpol_rdinput] Number of Records: ",numrec)
     call message_value(MESLOG,"[fit1Dpol_rdinput] Number of Columns: ",numcol)
  end if
  !
  allocate(f1dp_data(numrec,numcol),STAT=irc)
  if(irc.lt.0) then
     fit1Dpol_rdinput = -1
     LOCfit1Dpol_error = 1
     call message(MESERRO,'[fit1Dpol_rdinput] Cannot allocate data matrix.')
     return
  end if
  !
  allocate(f1dp_idxdata(numrec),STAT=irc)
  if(irc.lt.0) then
     fit1Dpol_rdinput = -1
     LOCfit1Dpol_error = 1
     call message(MESERRO,'[fit1Dpol_rdinput] Cannot allocate data index matrix.')
     return
  end if
  !
  allocate(f1dp_ypol(numrec),STAT=irc)
  if(irc.lt.0) then
     fit1Dpol_rdinput = -1
     LOCfit1Dpol_error = 1
     call message(MESERRO,'[fit1Dpol_rdinput] Cannot allocate ypol matrix.')
     return
  end if
  !
  i=0
  done = .false.
  do while(.not.done)
     irc = line_getline(f1dp_iunin,line,3)
     if (irc.lt.0) then
        done = .true.
        cycle
     end if
     !
     if(line(1:1).eq.'[') then   !New Section
        if(i.ne.numrec) then
           call message(MESERRO,"[fit1dpol_rdinput] less record that specified: decrease numrec.")
           call message_value(MESWARN,"[fit1dpol_rdinput] Actual number of record readed: ",i)
           call message(MESWARN,"[fit1dpol_rdinput] We ignore the error and set a new numrec.")
           numrec=i
           f1dp_numrec=numrec
        end if
        done = .true.
        cycle
     end if
     !
     i=i+1
     if(i.gt.numrec) then
        call message(MESERRO,"[fit1dpol_rdinput] Exceeded maximum number of data.")
        call message(MESWARN,"[fit1dpol_rdinput] We ignore the extra data.")
        i=i-1
        done = .true.
        cycle
     end if
     !
     call line_read(line,iwrds,rwrds,cwrds,niwrd,nrwrd,ncwrd,lrform)
     if(nrwrd.lt.numcol) then
        call message(MESERRO,"[fit1dpol_rdinput] Incompatible number of columns.")
        call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
        call message_value(MESERRO,"[fit1dpol_rdinput] Number of Reals is: ",nrwrd)
        call message_value(MESERRO,"[fit1dpol_rdinput] Number of  Ints is: ",niwrd)
        call message_value(MESERRO,"[fit1dpol_rdinput] Number of Chars is: ",ncwrd)
        call message(MESERRO,line)
        call message(MESKILL,"[fit1dpol_rdinput] We exit now.")
     end if
     if(nrwrd.gt.numcol) then
        call message(MESWARN,"[fit1dpol_rdinput] More columns than actual readed.")
        call message_value(MESWARN,"[fit1dpol_rdinput] Input section line: ",i)
        call message(MESWARN,"[fit1dpol_rdinput] We ignore and continue.")
     end if
     if(niwrd.gt.1.or.ncwrd.ne.0) then
        call message(MESERRO,"[fit1dpol_rdinput] Detected spurius data: please check the input.")
        call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
        call message(MESWARN,"[fit1dpol_rdinput] We ignore the extra data.")
     end if
     !
     if(niwrd.gt.0) then
        f1dp_idxdata(i)=iwrds(1)
     else
        f1dp_idxdata(i)=i
     end if
     do j=1,numcol
        f1dp_data(i,j)=rwrds(j)
     end do
     !
  end do!END while(.not.done)
  !
  ! Check input echo
  !-----------------
  if(f1dp_flags(FLG_EINP).gt.0) then
     do i=1,numrec
        write(f1dp_iunout,'(1X,I6,$)') f1dp_idxdata(i)
        do j=1,numcol
           write(f1dp_iunout,'(1X,F18.12,$)') f1dp_data(i,j)
        end do
        write(f1dp_iunout,'(1X)')
     end do
  end if
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! READING THE PARAMETER SECTION
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  rewind(f1dp_iunin)
  params=''
  irc = osec_set(f1dp_iunin,'x-fit1dpol-pars',params)
  if(irc.lt.0) then
     fit1Dpol_rdinput = -1
     LOCfit1Dpol_error = 1
     call message(MESERRO,'[fit1Dpol_rdinput] No x-fit1dpol-pars section.')
     return
  end if
  !
  irc=line_getval(params,"numpars",cnumpars)
  if(irc.lt.0) then
     fit1Dpol_rdinput = -1
     LOCfit1Dpol_error = 1
     call message(MESERRO,'[fit1Dpol_rdinput] This version require numpars=N as parameter.')
     return
  end if
  call line_getform(cnumpars,'I',rfrmt,i)
  rfrmt='('//rfrmt(1:i)//')'
  read (cnumpars,rfrmt(1:i+2)) numpars
  f1dp_numpar = numpars
  if (f1dp_debug.gt.0) then
     call message_value(MESLOG,"[fit1Dpol_rdinput] Number of Parameters: ",numpars)
  end if
  !
  irc=line_getval(params,"eps",cleps)
  if(irc.lt.0) then
     f1dp_flags(FLG_IEPS)=0
  else
     f1dp_flags(FLG_IEPS)=1
     call line_getform(cleps,'D',rfrmt,i)
     rfrmt='('//rfrmt(1:i)//')'
     read (cleps,rfrmt) leps
     f1dp_eps = leps
     if (f1dp_debug.gt.0) then
        call message_value(MESLOG,"[fit1Dpol_rdinput] EPS precision set by input to: ",f1dp_eps)
     end if
  end if
  !
  allocate(f1dp_pars(numpars,PAR_SIZE),STAT=irc)
  if(irc.lt.0) then
     fit1Dpol_rdinput = -1
     LOCfit1Dpol_error = 1
     call message(MESERRO,'[fit1Dpol_rdinput] Cannot allocate parameters matrix.')
     return
  end if
  !
  allocate(f1dp_idxpars(numpars),STAT=irc)
  if(irc.lt.0) then
     fit1Dpol_rdinput = -1
     LOCfit1Dpol_error = 1
     call message(MESERRO,'[fit1Dpol_rdinput] Cannot allocate parameters index matrix.')
     return
  end if
  !
  allocate(f1dp_parnames(numpars),STAT=irc)
  if(irc.lt.0) then
     fit1Dpol_rdinput = -1
     LOCfit1Dpol_error = 1
     call message(MESERRO,'[fit1Dpol_rdinput] Cannot allocate parameter names matrix.')
     return
  end if
  !
  allocate(f1dp_strfixed(numpars),STAT=irc)
  if(irc.lt.0) then
     fit1Dpol_rdinput = -1
     LOCfit1Dpol_error = 1
     call message(MESERRO,'[fit1Dpol_rdinput] Cannot allocate parameter names matrix.')
     return
  end if
  !
  allocate(f1dp_fixpars(numpars),STAT=irc)
  if(irc.lt.0) then
     fit1Dpol_rdinput = -1
     LOCfit1Dpol_error = 1
     call message(MESERRO,'[fit1Dpol_rdinput] Cannot allocate parameter names matrix.')
     return
  end if
  !
  irc=line_getval(params,"dimension",cmdim)
  if(irc.eq.0) then
     call line_getform(cmdim,'I',rfrmt,i)
     rfrmt='('//rfrmt(1:i)//')'
     read (cmdim,rfrmt(1:i+2)) f1dp_flags(FLG_POLO)
     if (f1dp_debug.gt.1) then
        call message_value(MESLOG,"[fit1Dpol_rdinput] MultiDim Dim set at input: ",f1dp_flags(FLG_POLO))
     end if
     if (f1dp_debug.gt.0.and.f1dp_flags(FLG_POLO).ne.(f1dp_numcol-1)) then
        call message_value(MESWARN,"[fit1Dpol_rdinput] MultiDim Dim set at input: ",f1dp_flags(FLG_POLO))
        call message_value(MESWARN,"[fit1Dpol_rdinput]  Differ from data columns: ",(f1dp_numcol-1))
     end if
  end if
  !
  irc=line_getval(params,"funct",cfun)
  if(irc.lt.0) then
     call message(MESWARN,'[fit1Dpol_rdinput] No funct=NAME in the input: assuming polynom.')
     f1dp_funame='polynom'
  else
     cfun=trim(cfun)
     i=len_trim(cfun)
     if(i.gt.10) i=10
     f1dp_funame=cfun(1:i)
  end if
  if (f1dp_debug.gt.0) then
     call message(MESLOG,"[fit1Dpol_rdinput] Function set at input is "//cfun(1:i))
  end if
  !
  i=0
  done = .false.
  do while(.not.done)
     irc = line_getline(f1dp_iunin,line,3)
     if (irc.lt.0) then
        done = .true.
        cycle
     end if
     !
     if(line(1:1).eq.'[') then   !New Section
        done = .true.
        cycle
     end if
     !
     i=i+1
     if(i.gt.numpars) then
        call message(MESERRO,"[fit1dpol_rdinput] Exceeded maximum number of parameters.")
        call message(MESWARN,"[fit1dpol_rdinput] We ignore the extra parameters.")
        i=i-1
        done = .true.
        cycle
     end if
     !
     call line_read(line,iwrds,rwrds,cwrds,niwrd,nrwrd,ncwrd,lrform)
     !
     ! note these warning are PAR_SIZE dependent!
     if(ncwrd.lt.1.or.nrwrd.lt.1) then
        call message(MESERRO,"[fit1dpol_rdinput] Cannot read parameter.")
        call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
        call message(MESWARN,"[fit1dpol_rdinput] We skip this line: Trying to continue.")
     end if
     if(ncwrd.gt.2.or.nrwrd.gt.2.or.niwrd.gt.1) then
        call message(MESERRO,"[fit1dpol_rdinput] Detected spurius data: please check the input.")
        call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
        call message(MESWARN,"[fit1dpol_rdinput] We skip this line: Trying to continue.")
     end if
     !
     ! get the initial value
     f1dp_pars(i,1)=rwrds(1)
     ! get the step size if there otherwise see the next comment
     if(nrwrd.gt.1) then
        f1dp_pars(i,2)=rwrds(2)
     else
        ! We have to guess the step size so we assume 1/1000 should be ok
        f1dp_pars(i,2)=abs(f1dp_pars(i,1))/1000.0_FREAL
     end if
     ! get the symbolic name ! WARN We truncate the string
     f1dp_parnames(i)=cwrds(1)(1:10)
     ! get if we need it to be fixed
     if(ncwrd.gt.1) then
        f1dp_strfixed(i)=cwrds(2)(1:1)
     else
        f1dp_strfixed(i)='n'
     end if
     ! we fix only if you put YES!
     if(f1dp_strfixed(i).eq.'Y'.or.f1dp_strfixed(i).eq.'y') then
        f1dp_fixpars(i)=1
     else
        f1dp_fixpars(i)=0
     end if
     ! get the new index if required
     if(niwrd.gt.0) then
        f1dp_idxpars(i)=iwrds(1)
     else
        f1dp_idxpars(i)=i
     end if
     !
  end do!END while(.not.done)
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! READING THE dataops SECTION
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  rewind(f1dp_iunin)
  params=''
  irc = osec_set(f1dp_iunin,'x-fit1dpol-dataops',params)
  if(irc.gt.0) then
     !
     if (f1dp_debug.gt.1) then
        call message(MESOUT,"[fit1dpol_rdinput] got dataops section.")
     end if
     !
     i=0
     done = .false.
     do while(.not.done)
        !
        i=i+1
        !
        irc = line_getline(f1dp_iunin,line,3)
        if (irc.lt.0) then
           done = .true.
           cycle
        end if
        !
        if(line(1:1).eq.'[') then   !New Section
           done = .true.
           cycle
        end if
        !
        call line_read(line,iwrds,rwrds,cwrds,niwrd,nrwrd,ncwrd,lrform)
        ! check if enough data
        if(ncwrd.lt.1.or.nrwrd.lt.1.or.niwrd.lt.1.or.niwrd.eq.2.or.nrwrd.eq.2) then
           call message(MESERRO,"[fit1dpol_rdinput] Cannot read operations.")
           call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
           call message(MESKILL,"[fit1dpol_rdinput] We stop now.")
        end if
        if(ncwrd.ne.1.or.nrwrd.gt.3.or.niwrd.gt.3) then
           call message(MESERRO,"[fit1dpol_rdinput] Spurius data in dataops section.")
           call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
           call message(MESWARN,"[fit1dpol_rdinput] We keep only valid data: Trying to continue.")
        end if
        if(index(cwrds(1),'shift').gt.0) then
           if (niwrd.gt.2) then
              if(iwrds(3).gt.f1dp_numrec.or.iwrds(2).lt.1) then
                 call message(MESERRO,"[fit1dpol_rdinput] Invalid operation definition: index out of bound.")
                 call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
                 call message(MESKILL,"[fit1dpol_rdinput] We stop now.")
              end if
              do j=iwrds(2),iwrds(3)
                 f1dp_data(j,iwrds(1))=f1dp_data(j,iwrds(1))+rwrds(1)
              end do
           else
              if (nrwrd.gt.2) then
                 do j=1,f1dp_numrec
                    if(f1dp_data(j,iwrds(1)).ge.rwrds(2).and.f1dp_data(j,iwrds(1)).le.rwrds(3)) then
                       f1dp_data(j,iwrds(1))=f1dp_data(j,iwrds(1))+rwrds(1)
                    end if
                 end do
              else
                 do j=1,f1dp_numrec
                    f1dp_data(j,iwrds(1))=f1dp_data(j,iwrds(1))+rwrds(1)
                 end do
              end if
           end if
        else if(index(cwrds(1),'scale').gt.0) then
           do j=1,f1dp_numrec
              f1dp_data(j,iwrds(1))=f1dp_data(j,iwrds(1))*rwrds(1)
           end do
        else
           call message(MESERRO,"[fit1dpol_rdinput] Cannot read operations: unknown operation.")
           call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
           call message(MESKILL,"[fit1dpol_rdinput] We stop now.")
        end if
        !
     end do
     !
     if (f1dp_debug.gt.0) then
        write (f1dp_iunout,'(1X,"Set of new transformed points is:")')
        !
        !write (f1dp_iunout,'(1X,"[x-fid1dpol-data] numrec=",I6,2X,"numcol=",I6)') f1dp_numrec, f1dp_numcol
        write (f1dp_iunout,'(1X,A,I6,2X,A,I6)')&
             &"[x-fid1dpol-data] numrec=",f1dp_numrec,"numcol=",f1dp_numcol
        do j=1,numrec
           k=0
           write(f1dp_iunout,'(1X,I6,3X,$)') j
           do k=1,numcol-1
              write(f1dp_iunout,'(F18.10,1X,$)') f1dp_data(j,k)
           end do
           write(f1dp_iunout,'(F18.10)') f1dp_data(j,numcol)
        end do
     end if
     !
  end if
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! READING THE what SECTION
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  rewind(f1dp_iunin)
  params=''
  irc = osec_set(f1dp_iunin,'x-fit1dpol-what',params)
  if(irc.gt.0) then
     !
     if (f1dp_debug.gt.1) then
        call message(MESOUT,"[fit1dpol_rdinput] got what section.")
     end if
     !
     ! We allocate and save a copy of the data if required
     ! WARNING These are already transformed data by dataops because dataops comes first
     !----------------------------------------------------------------------------------
     if(index(params,'save').gt.0) then
        allocate(f1dp_cpydata(f1dp_numrec,f1dp_numcol),STAT=irc)
        if(irc.lt.0) then
           fit1Dpol_rdinput = -1
           LOCfit1Dpol_error = 1
           call message(MESERRO,'[fit1Dpol_rdinput] Cannot allocate copy data matrix.')
           return
        end if
        if(f1dp_debug.gt.1) then
           call message(MESLOG,"[fit1Dpol_rdinput] We have allocated a copy of the data.")
        end if
        !
        f1dp_cpydata(:,:)=f1dp_data(:,:)
     end if
     !
     i=0
     done = .false.
     do while(.not.done)
        !
        i=i+1
        !
        irc = line_getline(f1dp_iunin,line,3)
        if (irc.lt.0) then
           done = .true.
           cycle
        end if
        !
        if(line(1:1).eq.'[') then   !New Section
           done = .true.
           cycle
        end if
        !
        call line_read(line,iwrds,rwrds,cwrds,niwrd,nrwrd,ncwrd,lrform)
        ! check if enough data
        if(ncwrd.lt.1) then
           call message(MESERRO,"[fit1dpol_rdinput] Cannot read what command.")
           call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
           call message(MESKILL,"[fit1dpol_rdinput] We stop now.")
        end if
        if(ncwrd.gt.1.or.nrwrd.gt.0.or..not.(niwrd.eq.6.or.niwrd.eq.2)) then
           call message(MESERRO,"[fit1dpol_rdinput] Spurius data in what section.")
           call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
           call message(MESWARN,"[fit1dpol_rdinput] We keep only valid data: Trying to continue.")
        end if
        !
        ! WARNING: We check input consistency now assuming all operations get
        !          6 or 2 integers and in case to we resize them.
        !          If new operations do not follow the rule you should recode this part!
        !-------------------------------------------------------------------------------
        if (.not.(niwrd.eq.6.or.niwrd.eq.2)) then
           call message(MESERRO,"[fit1dpol_rdinput] Invalid what operation definition: We require 6 integers.")
           call message_value(MESERRO,"[fit1dpol_rdinput] We got: ",niwrd)
           call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
           call message(MESKILL,"[fit1dpol_rdinput] We stop now.")
        end if
        if (nrwrd.gt.0) then
           call message(MESERRO,"[fit1dpol_rdinput] Spurius data (real number found): We assume error.")
           call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
           call message(MESKILL,"[fit1dpol_rdinput] We stop now.")
        end if
        !
        ! if only 2 integers means all column so we reset the data at proper values
        !--------------------------------------------------------------------------
        if (niwrd.eq.2) then
           if (iwrds(1).eq.iwrds(2)) then
              call message(MESERRO,"[fit1dpol_rdinput] You surely did  something wrong in what.")
              call message_value(MESERRO,"[fit1dpol_rdinput] copy on itself column in input section line: ",i)
              call message(MESKILL,"[fit1dpol_rdinput] We stop now.")
           end if
           niwrd=6
           !iwrds(1)=iwrds(1) ! just for comptibility ;)
           iwrds(4)=iwrds(2)
           iwrds(2)=1
           iwrds(3)=f1dp_numrec
           iwrds(5)=1
           iwrds(6)=f1dp_numrec
        end if
        !
        ! Check consistency
        !------------------
        if((.not.(iwrds(2).le.iwrds(3).or.iwrds(5).le.iwrds(6))).or.((iwrds(3)-iwrds(2)).ne.(iwrds(6)-iwrds(5)))) then
           call message(MESERRO,"[fit1dpol_rdinput] Conformation Error in what.")
           call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
           call message(MESKILL,"[fit1dpol_rdinput] We stop now.")
        end if
        !
        if(iwrds(1).eq.iwrds(4)) then !Same column
           if((iwrds(5).ge.iwrds(2).and.iwrds(5).le.iwrds(3)).or.(iwrds(6).ge.iwrds(2).and.iwrds(6).le.iwrds(3))) then
              call message(MESWARN,"[fit1dpol_rdinput] Overlapping regions in what operation.")
           end if
        end if
        !
        ! Starts operations
        !------------------
        if(index(cwrds(1),'copy').gt.0) then
           ! vector copy
           !-------------------
           f1dp_data(iwrds(5):iwrds(6),iwrds(4))=f1dp_data(iwrds(2):iwrds(3),iwrds(1))
        else if(index(cwrds(1),'addr').gt.0) then
           f1dp_data(iwrds(5):iwrds(6),iwrds(4))=f1dp_data(iwrds(2):iwrds(3),iwrds(1))+f1dp_data(iwrds(5):iwrds(6),iwrds(4))
        else if(index(cwrds(1),'addl').gt.0) then
           f1dp_data(iwrds(2):iwrds(3),iwrds(1))=f1dp_data(iwrds(2):iwrds(3),iwrds(1))+f1dp_data(iwrds(5):iwrds(6),iwrds(4))
        else if(index(cwrds(1),'subr').gt.0) then
           f1dp_data(iwrds(5):iwrds(6),iwrds(4))=f1dp_data(iwrds(2):iwrds(3),iwrds(1))-f1dp_data(iwrds(5):iwrds(6),iwrds(4))
        else if(index(cwrds(1),'subl').gt.0) then
           f1dp_data(iwrds(2):iwrds(3),iwrds(1))=f1dp_data(iwrds(2):iwrds(3),iwrds(1))-f1dp_data(iwrds(5):iwrds(6),iwrds(4))
        else if(index(cwrds(1),'mulr').gt.0) then
           f1dp_data(iwrds(5):iwrds(6),iwrds(4))=f1dp_data(iwrds(2):iwrds(3),iwrds(1))*f1dp_data(iwrds(5):iwrds(6),iwrds(4))
        else if(index(cwrds(1),'mull').gt.0) then
           f1dp_data(iwrds(2):iwrds(3),iwrds(1))=f1dp_data(iwrds(2):iwrds(3),iwrds(1))*f1dp_data(iwrds(5):iwrds(6),iwrds(4))
        else if(index(cwrds(1),'divr').gt.0) then
           f1dp_data(iwrds(5):iwrds(6),iwrds(4))=f1dp_data(iwrds(2):iwrds(3),iwrds(1))/f1dp_data(iwrds(5):iwrds(6),iwrds(4))
        else if(index(cwrds(1),'divl').gt.0) then
           f1dp_data(iwrds(2):iwrds(3),iwrds(1))=f1dp_data(iwrds(2):iwrds(3),iwrds(1))/f1dp_data(iwrds(5):iwrds(6),iwrds(4))
        else if(index(cwrds(1),'use12s').gt.0) then
           f1dp_numrec=iwrds(3)-iwrds(2)+1
           numrec=f1dp_numrec
           f1dp_data(1:f1dp_numrec,1)=f1dp_data(iwrds(2):iwrds(3),iwrds(1))
           f1dp_data(1:f1dp_numrec,2)=f1dp_data(iwrds(5):iwrds(6),iwrds(4))
        !else if(index(cwrds(1),'add').gt.0) then
        !   do j=1,f1dp_numrec
        !      f1dp_data(j,iwrds(1))=f1dp_data(j,iwrds(1))*rwrds(1)
        !   end do
        !else if(index(cwrds(1),'add').gt.0) then
        !   do j=1,f1dp_numrec
        !      f1dp_data(j,iwrds(1))=f1dp_data(j,iwrds(1))*rwrds(1)
        !   end do
        else
           call message(MESERRO,"[fit1dpol_rdinput] Cannot read operations: unknown operation.")
           call message_value(MESERRO,"[fit1dpol_rdinput] Input section line: ",i)
           call message(MESKILL,"[fit1dpol_rdinput] We stop now.")
        end if
        !
     end do
     !
     if (f1dp_debug.gt.0) then
        write (f1dp_iunout,'(1X,"Set of new transformed points is:")')
        !write (f1dp_iunout,'(1X,"[x-fid1dpol-data] numrec=",I6,2X,"numcol=",I6)') f1dp_numrec, f1dp_numcol
        write (f1dp_iunout,'(1X,A,I6,2X,A,I6)')&
             &"[x-fid1dpol-data] numrec=",f1dp_numrec,"numcol=",f1dp_numcol
        do j=1,numrec
           k=0
           write(f1dp_iunout,'(1X,I6,3X,$)') j
           do k=1,numcol-1
              write(f1dp_iunout,'(F18.10,1X,$)') f1dp_data(j,k)
           end do
           write(f1dp_iunout,'(F18.10)') f1dp_data(j,numcol)
        end do
     end if
     !
  end if
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! END READING SECTIONS
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !
end function fit1Dpol_rdinput
!
!H
!H-----------------------------------------------------------------------------
!H
!
  !
  !
end program fit1Dpol
