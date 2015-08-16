!
!H
!H-----------------------------------------------------------------------------
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H MODULE csmmod  Freemol by F.Mariotti: (c) F.Mariotti
!H-----------------------------------------------------------------------------
!H $Id: csmmod.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
!H-----------------------------------------------------------------------------
!H 
!H DESCRIPTION:
!H Name convention used within this module:
!H csmmod_XXX   are Public Procedures/Variables
!H csm_XXX      are Private Procedures/Variables
!H
!H Each public routine should set the module error code to 0.
!H So each routine contains a line as:
!H    csm_error = 0
!H 
!H PUBLIC PROCEDURES:
!H 
!H   public csmmod_init
!H          This routine initialize the csmmod module.
!H          At present sets the initialization and error flags and obtains a coputed value
!H          for PI (180 degree).
!H
!H          RETURNS[ineger(FINT)]: 0 always.
!H
!H          INPUT: NONE
!H
!H          OUTPUT: NONE
!H 
!H   public csmmod_isinit
!H          Returns a logical stating true if the module has been initialized.
!H
!H          RETURNS[logical]: .true. (for initialised) or .false. (for not initialised).
!H
!H          INPUT: NONE
!H
!H          OUTPUT: NONE
!H 
!H   public csmmod_geterror
!H          It is an interface and it has this two masks:
!H          integer(FINT) function csmmod_geterror_i()
!H          character(FLCHARS) function csmmod_geterror_c(icode)
!H                  integer(FINT), optional, intent(in) :: icode
!H          It is a standard modules routine to report errors.
!H          
!H 
!H   public csmmod_rdinput
!H          This routine read input data from a molden/freemol/Frimol file.
!H          
!H          
!H 
!H   public csmmod_getnat
!H          
!H          
!H          
!H 
!H   public csmmod_getnsymop
!H          
!H          
!H          
!H 
!H   public csmmod_getrho
!H          
!H          
!H          
!H 
!H   public csmmod_getrrii
!H          
!H          
!H          
!H 
!H   public csmmod_getrrij
!H          
!H          
!H          
!H 
!H   public csmmod_getimj
!H          
!H          
!H          
!H 
!H   public csmmod_euler
!H          
!H          
!H          
!H 
!H   public csmmod_rotxyz
!H          
!H          
!H          
!H 
!H 
!H 
!H PUBLIC VARIABLES:
!H 
!H atg_index      
!H   integer(FINT), allocatable, save, public, dimension(:) :: atg_index
!H   real(FREAL), allocatable, save, public, dimension(:,:) :: atg_pars
!H   ! SymOp Data
!H   !
!H   logical, save, public :: symopt
!H   logical, save, public, allocatable, dimension(:) :: syminv
!H   integer(KINT), save, public :: nsymop = DEF_nsymop
!H   integer(KINT), save, public :: symident = 0
!H   character(LCHARS), save, public,allocatable, dimension(:) :: symopstr
!H   real(KREAL), save, public, allocatable, dimension(:) :: symalpha
!H   real(KREAL), save, public, allocatable, dimension(:,:) :: symaxe
!H   real(KREAL), save, public, allocatable, dimension(:,:) :: symtrans
!H   real(KREAL), save, public, allocatable, dimension(:,:,:) :: symrmat
!H   real(KREAL), save, public, allocatable, dimension(:,:) :: symeuler
!H   real(KREAL), save, public, allocatable, dimension(:) :: symcsmval
!H   !
!H   ! More options
!H   !-------------
!H   logical, save, public :: csm_noshift = .false.
!H   logical, save, public :: csm_onestep = .false.
!H   logical, save, public :: csm_noinertia = .false.
!H   logical, save, public :: csm_useminuit = .false.
!H   logical, save, public :: csm_readminuit = .false.
!H   logical, save, public :: csm_interminuit = .false.
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
module csmmod
  !
  use vartypes
  use messages
  use baseio
  use osec
  use linetools
  use molecule
  !
  implicit none
  private
  !
  ! PUBLIC DECLARATIONS:
  !--------------------
  !
  ! Public Procedures
  !------------------
  public csmmod_init
  public csmmod_isinit
  public csmmod_geterror
  public csmmod_rdinput
  public csmmod_getnat
  public csmmod_getnsymop
  public csmmod_getrho
  public csmmod_getrrii
  public csmmod_getrrij
  public csmmod_getimj
  public csmmod_getimpj
  public csmmod_euler
  public csmmod_rotxyz
  public csmmod_setrdminuit
  ! public csmmod_
  !
  !
  ! Error Interface
  !----------------
  ! We define two error interfaces. Each one return the latest error
  ! code or string value. The string version can inquire about an error code.
  !--------------------------------------------------------------------------
  interface csmmod_geterror
     module procedure csmmod_geterror_i
     module procedure csmmod_geterror_c
  end interface !csmmod_geterror

  ! LOCAL DECLARATIONS:
  !--------------------
  !
  logical, save :: csm_isinit = .false.
  !
  ! The gaussian functions
  !-----------------------
  integer(FINT), allocatable, save, public, dimension(:) :: atg_index
  real(FREAL), allocatable, save, public, dimension(:,:) :: atg_pars
  !
  !Common data
  !-----------
  integer(FINT), save, private :: nat
  real(FREAL), save, private :: pi
  !
  ! SymOp Data
  !
  logical, save, public :: symopt
  logical, save, public, allocatable, dimension(:) :: syminv
  integer(KINT), parameter :: DEF_nsymop = 180
  integer(KINT), save, public :: nsymop = DEF_nsymop
  integer(KINT), save, public :: symident = 0
  integer(KINT), save, public :: syminvidx = 0
  character(LCHARS), save, public,allocatable, dimension(:) :: symopstr
  real(KREAL), save, public, allocatable, dimension(:) :: symalpha
  real(KREAL), save, public, allocatable, dimension(:,:) :: symaxe
  real(KREAL), save, public, allocatable, dimension(:,:) :: symtrans
  real(KREAL), save, public, allocatable, dimension(:,:,:) :: symrmat
  !real(KREAL), save, public, allocatable, dimension(:,:) :: symeuler
  real(KREAL), allocatable, dimension(:,:) :: symeuler
  real(KREAL), save, public, allocatable, dimension(:) :: symcsmval
  !
  ! More options
  !-------------
  logical, save, public :: csm_noshift = .false.
  logical, save, public :: csm_onestep = .false.
  logical, save, public :: csm_noinertia = .false.
  logical, save, public :: csm_useminuit = .false.
  logical, save, public :: csm_readminuit = .false.
  logical, save, public :: csm_interminuit = .false.
  !
  ! error variables
  !----------------
  integer(FINT), save, private :: csm_error
  character(FLCHARS), save, private :: csm_error_message
  !
contains
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function csmmod_init()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function csmmod_init()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the csmmod module.
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    !
    csm_isinit = .true.
    csm_error = 0
    csmmod_init = 0
    !
    ! Evaluate pi
    pi  = atan2(1.0_FREAL,1.0_FREAL) * 4.0_FREAL
    !
    !
  end function csmmod_init
!
!H
!H-----------------------------------------------------------------------------
!H logical function csmmod_isinit()
!H-----------------------------------------------------------------------------
!H
!
  logical function csmmod_isinit()
!
!H
!H-----------------------------------------------------------------------------
!H We merelyreturn the initialization flag.
!H-----------------------------------------------------------------------------
!H
!
    csm_error = 0
    csmmod_isinit = csm_isinit
    !
  end function csmmod_isinit
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function csmmod_geterror_i()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function csmmod_geterror_i()
!
!H
!H-----------------------------------------------------------------------------
!H Return the last integer error code of the module.
!H-----------------------------------------------------------------------------
!H
    ! We merely return the integer current error code
    !------------------------------------------------
    csmmod_geterror_i = csm_error
    !
  end function csmmod_geterror_i
!
!H-----------------------------------------------------------------------------
!H
!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H character(FLCHARS) function csmmod_geterror_c(icode)
!H-----------------------------------------------------------------------------
!H
!
  character(FLCHARS) function csmmod_geterror_c(icode)
    !
!    integer(FINT), optional, intent(in) :: icode
    integer(FINT), intent(in) :: icode
    integer(FINT) :: code
    !
!
!H
!H-----------------------------------------------------------------------------
!H 
!H Return the latest error message of the module if any or the error string
!H if a code is given as input.
!H 
!H-----------------------------------------------------------------------------
!H
!
    ! If a code is given we return that message.
    ! If a code is note given we return the current status.
    !------------------------------------------------------
!    if(present(icode)) then
       code = icode
!    else
!       code = csm_error
!    end if
    !
    ! We return the error message
    !----------------------------
    ! Note that we have not yet implement it completly.
    ! Some error codes are just fake messagges.
    !-------------------------------------------------------
    select case(code)
    case(0_FINT)
       csmmod_geterror_c = "No Errors."
    case(1_FINT)
       csmmod_geterror_c = "Not documented error or generic error.."
    case(2_FINT)
       csmmod_geterror_c = "[csmmod_rdinput] Cannot allocate cmd line parameters"
    case(3_FINT)
       csmmod_geterror_c = "[csmmod_rdinput] Cannot Open Input/Output File"
    case(4_FINT)
       csmmod_geterror_c = "[csmmod_rdinput]&
            & -l/-e Require a file name. (- is a valid file name."
    case(5_FINT)
       csmmod_geterror_c = "Double Error: Why did you get this message?."
    case(6_FINT)
       csmmod_geterror_c = "Double Error: Why did you get this message?."
    case(7_FINT)
       csmmod_geterror_c = "Double Error: Why did you get this message?."
    case(8_FINT)
       csmmod_geterror_c = "Double Error: Why did you get this message?."
    case default
       csmmod_geterror_c = "Not documented error."
    end select
    !
    return
    !
  end function csmmod_geterror_c
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function csmmod_rdinput&
!H               &(iunin,finput,iunout,foutput,iunlog,filelog,iunerr,fileerr)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function csmmod_rdinput&
       &(iunin,finput,iunout,foutput,iunlog,filelog,iunerr,fileerr)
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine read input data from a molden/freemol/Freemol file.
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    !
    use pcmdline
    !
    character(FLCHARS), intent(out) :: finput, foutput, filelog, fileerr
    integer(FINT), intent(out) :: iunin, iunout, iunlog, iunerr
    !
    integer(FINT) :: nargs
    character(FLCHARS),dimension(:), allocatable :: cargs
    character(FLCHARS) :: params
    integer(FINT) :: irc
    !character(FLCHARS) :: cdummy
    !
    ! The buffer input file
    !----------------------
    integer(FINT) :: iunbin
    character(FLCHARS) :: fbufin = "CSMD.bufferInput"
    !
    ! Defaults
    !---------
    csm_error = 0
    csmmod_rdinput = 0
    !
    ! We read the command line for I/O file names
    !--------------------------------------------
    call pcmd_iargc(nargs)
    allocate(cargs(nargs),STAT=irc)
    if(.not.irc.eq.0) then
       csm_error = 2
       ! csm_error_message = 'Too many cmd line parameters'
       csmmod_rdinput = -1
       call message(MESDEBG,&
            &'[csmmod_rdinput] Cannot allocate cmd line parameters')
       return
    end if
    call pcmd_getarg(nargs,cargs)
    call pcmd_getio(nargs,cargs,finput,foutput)
    !
    ! We open the input and output files
    !-----------------------------------
    finput = trim(finput)
    if(scan(finput," ").eq.0) then
       iunin = 5
       finput = "stdin"
    else
       irc = baseio_open(iunin,finput,stat='OLD')
       if(irc.lt.0) then
          csm_error = 3
          csmmod_rdinput = -1
          call message(MESERRO,"[csmmod_rdinput] Cannot Open Input File")
          return
       end if
    end if
    foutput = trim(foutput)
    if(scan(foutput," ").eq.0) then
       iunout = 6
       foutput = "stdout"
    else
       irc = baseio_open(iunout,foutput,stat='UNKNOWN')
       if(irc.lt.0) then
          csm_error = 3
          csmmod_rdinput = -1
          call message(MESERRO,"[csmmod_rdinput] Cannot Open Output File")
          return
       end if
    end if
    !
    ! We check other options
    !-----------------------
    ! LOGFILE
    filelog = "stdlog"
    irc = pcmd_checkarg(nargs,cargs,'-l',params)
    if(irc.lt.0) then
       csm_error = 1
       csmmod_rdinput = -1
       call message(MESERRO,"[csmmod_rdinput]&
            & -l Require a file name. (- is a valid file name.")
       return
    else if (irc.gt.0) then
       params = trim(params)
       if(scan(params," ").eq.0) then
          csm_error = 4
          csmmod_rdinput = -1
          call message(MESERRO,"[csmmod_rdinput]&
               & -l Require a file name. (- is a valid file name.")
          return
       end if
       filelog = trim(params)
       if ((len_trim(filelog).eq.1).and.(filelog(1:1).eq."-")) then
          iunlog = iunout
          filelog = "stdlog"
       else
          irc = baseio_open(iunlog,filelog,stat='UNKNOWN')
          if(irc.lt.0) then
             csm_error = 3
             csmmod_rdinput = -1
             call message(MESERRO,"[csmmod_rdinput] Cannot Open Log File")
             return
          end if
       end if
    else
       iunlog = iunout
    end if
    ! ERRFILE
    fileerr="stderr"
    irc = pcmd_checkarg(nargs,cargs,'-e',params)
    if(irc.lt.0) then
       csm_error = 1
       csmmod_rdinput = -1
       call message(MESERRO,"[csmmod_rdinput]&
            & -e Require a file name. (- is a valid file name.")
       return
    else if(irc.gt.0) then
       params = trim(params)
       if(scan(params," ").eq.0) then
          csm_error = 4
          csmmod_rdinput = -1
          call message(MESERRO,"[csmmod_rdinput]&
               & -e Require a file name. (- is a valid file name.")
          return
       end if
       fileerr = trim(params)
       if ((len_trim(fileerr).eq.1).and.(fileerr(1:1).eq."-")) then
          iunerr = iunout
          fileerr = "stderr"
       else
          irc = baseio_open(iunerr,fileerr,stat='UNKNOWN')
          if(irc.lt.0) then
             csm_error = 3
             csmmod_rdinput = -1
             call message(MESERRO,"[csmmod_rdinput] Cannot Open Err File")
             return
          end if
       end if
    else
       iunerr = iunout
    end if
    !
    ! We check for shift
    !-------------------
    csm_noshift = .false.
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-p',params)
    irc = index(params,'noshift')
    if(irc.gt.0) then
       csm_noshift = .true.
       write(iunout,'("# NOShift option selected.")')
    end if
    !
    ! We check for no inertia rotation
    !---------------------------------
    csm_noinertia = .false.
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-p',params)
    irc = index(params,'noinertia')
    if(irc.gt.0) then
       csm_noinertia = .true.
       write(iunout,'("# NOInertia option selected.")')
    end if
    !
    ! We check for onestep
    !---------------------
    csm_onestep = .false.
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-p',params)
    irc = index(params,'onestep')
    if(irc.gt.0) then
       csm_onestep = .true.
       write(iunout,'("# OneStep option selected.")')
    end if
    !
    ! We check for Minuit
    !---------------------
    csm_useminuit = .false.
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-p',params)
    irc = index(params,'useminuit')
    if(irc.gt.0) then
       csm_useminuit = .true.
       write(iunout,'("# UseMinuit option selected.")')
    end if
    !
    ! We check for ReadMinuit
    !---------------------
    csm_readminuit = .false.
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-p',params)
    irc = index(params,'readminuit')
    if(irc.gt.0) then
       csm_readminuit = .true.
       write(iunout,'("# ReadMinuit option selected.")')
    end if
    !
    ! We check for Interactive Minuit
    !--------------------------------
    csm_interminuit = .false.
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-p',params)
    irc = index(params,'interminuit')
    if(irc.gt.0) then
       csm_interminuit = .true.
       write(iunout,'("# Interactive Minuit option selected.")')
    end if
    !
    ! The Help
    !---------
    params =''
    irc = pcmd_checkarg(nargs,cargs,'-h',params)
    if(irc.gt.0) then
       irc = index(params,'full')
       if(irc.gt.0) then
          call csmmod_printhelp(6,'E')
          stop 0
       else
          call csmmod_printhelp()
          stop 0
       end if
    end if
    !
    ! We setup message function again
    !--------------------------------
    ! We send messages to the output file
    ! and all the rest in the proper one
    ! messages_init should be:
    !                 (messgs,output,errors,logs  ,optional mode)
    call messages_init(iunout,iunout,iunerr,iunlog)
    ! We annouce it
    call message(MESOUT,"[csmmod_rdinput] messages reinitialized.")
    call message(MESLOG,"[csmmod_rdinput] messages reinitialized.")
    ! We print out some infos about files
    !------------------------------------
    call message_value(MESOUT,"Input  file unit is:",iunin)
    call message_value(MESOUT,"Input  file name is:",finput)
    call message_value(MESOUT,"Output file unit is:",iunout)
    call message_value(MESOUT,"Output file name is:",foutput)
    call message_value(MESOUT,"Error  file unit is:",iunerr)
    call message_value(MESOUT,"Error  file name is:",fileerr(1:len_trim(fileerr)))
    call message_value(MESOUT,"Log    file unit is:",iunlog)
    call message_value(MESOUT,"Log    file name is:",filelog)
    !
    ! Print out infos about the readed files
    ! we use a buffer file to store input information
    !------------------------------------------------
    irc = baseio_open(iunbin,fbufin,stat='UNKNOWN')
    if(irc.lt.0) then
       csm_error = 1
       csmmod_rdinput = -1
       call message(MESERRO,"[csmmod_rdinput] Cannot Open Buffer File")
       return
    end if
    !
    ! Files are done we start to parse the input
    !-------------------------------------------
    write(iunbin,'(A14)') '[input-parsed]'
    irc = csm_rdinput(iunin,iunbin)
    !
    !
  end function csmmod_rdinput
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function csm_rdinput(iunin,iunbin)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function csm_rdinput(iunin,iunbin)
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine parse input data from a molden/freemol file (Real)
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    use chemconst
    !
    !
    integer(FINT), intent(in) :: iunin,iunbin
    !
    integer(FINT) :: irc,ig
    character(FLCHARS) :: params
    character(FLCHARS) :: line
    logical :: done
    !
    csm_rdinput=0
    !
    ! We need the molecule section
    !-----------------------------
    rewind(iunin)
    irc = osec_set(iunin,'molecule',params)
    if(irc.lt.0) then
       csm_rdinput = -1
       csm_error = 1
       call message(MESERRO,'[csmmod_rdinput] No molecule section.')
       return
    end if
    !
    call message_comment(iunbin,'Got molecule section. With Parameters:')
    call message_comment(iunbin,adjustl(params))
    !
    !
    call message(MESLOG,"[csm_rdinput] molecule_init called with debug.")
    irc=molecule_init(.true.)
    if(irc.lt.0) then
       csm_rdinput = -1
       csm_error = 1
       call message(MESERRO,'[csm_rdinput] Cannot Initialize the molecule.')
       return
    end if
    !
    call molecule_read(iunin,params)
    call molecule_print(iunbin)
    !
    ! we perform a conversion to au
    !------------------------------
    call molecule_convert(cc_au2ang)
    !call molecule_print(iunbin)
    !
    ! Now we read the gaussian functions in a dummy way
    !--------------------------------------------------
    rewind(iunin)
    irc = osec_set(iunin,'x-csmg-cgauss',params)
    if(irc.lt.0) then
       csm_rdinput = -1
       csm_error = 1
       call message(MESERRO,'[csmmod_rdinput] No gaussian functions section.')
       return
    end if
    !
    nat = molecule_getnat()
    allocate(atg_index(nat),STAT=irc)
    if(irc.ne.0) then
       call message(MESERRO,"[csmmod_rdinput] cannot allocate atg_index.")
    end if
    allocate(atg_pars(nat,2),STAT=irc)
    if(irc.ne.0) then
       call message(MESERRO,"[csmmod_rdinput] cannot allocate atg_pars.")
    end if
    write(iunbin,'("[x-csmg-cgauss] Nrec=",I4)') nat
    do ig=1,nat
       irc = line_getline(iunin,line,3)
       read(line,*) atg_index(ig),atg_pars(ig,1),atg_pars(ig,2)
       write(iunbin,'(1X,I4,1x,2(2X,F18.12))')&
            &atg_index(ig),atg_pars(ig,1),atg_pars(ig,2)
    end do
    !
    ! Now The SymOps
    !---------------
    rewind(iunin)
    irc = osec_set(iunin,'x-csmg-symop',params)
    if(irc.lt.0) then
       csm_rdinput = -1
       csm_error = 1
       call message(MESERRO,'[csmmod_rdinput] No symmetry operations section.')
       return
    end if
    !
    allocate(symcsmval(nsymop),STAT=irc)
    if(irc.ne.0) then
       stop 'cannot allocate for symop csmval.'
    endif
    allocate(syminv(nsymop),STAT=irc)
    if(irc.ne.0) then
       stop 'cannot allocate for symop logicals.'
    endif
    syminv(:) = .false.
    !
    allocate(symopstr(nsymop),STAT=irc)
    if(irc.ne.0) then
       stop 'cannot allocate for symop strings.'
    endif
    allocate(symalpha(nsymop),symaxe(3,nsymop),&
         &symtrans(3,nsymop),symrmat(3,3,nsymop),STAT=irc)
    if(irc.ne.0) then
       stop 'cannot allocate for symop reals.'
    endif
    allocate(symeuler(3,nsymop),STAT=irc)
    if(irc.ne.0) then
       stop 'cannot allocate for symop eulers.'
    endif
    symtrans(:,:) = 0.0_KREAL
    !
    ig = 0
    done = .false.
    do while(.not.done)
       irc = line_getline(iunin,line,3)
       if (irc.lt.0) then
          done = .true.
       else
          !
          if(line(1:1).eq.'[') then   !New Section
             done = .true.
          else
             ig = ig + 1
             if(ig.gt.nsymop) then
                call message(MESERRO,"[csm_rdinput] Exceeded maximum number of operations.")
                call message(MESKILL,"[csm_rdinput] We stop now.")
             end if
             read(line,*) symopstr(ig),symaxe(1,ig),symaxe(2,ig),symaxe(3,ig),&
                  &symalpha(ig),symtrans(1,ig),symtrans(2,ig),symtrans(3,ig)
             !WARNING CHECK HOW WE READ NOW THE INVERSION
             params = 'F'
             if(index(line(10:),'i').gt.0) then
                syminv(ig) = .true.
                params = 'T'
             end if
             if(index(line(10:),'I').gt.0) then
                syminv(ig) = .true.
                params = 'T'
             end if
             symalpha(ig)=symalpha(ig) *&
                  &(4.0_KREAL*atan(1.0_KREAL)) / 180.0_KREAL
             !write(*,*) line
             !write(*,'(1X,A8,1X,7(1X,F12.4),1X,A1)') &
             !     &symopstr(ig),symaxe(1,ig),symaxe(2,ig),symaxe(3,ig),&
             !     &symalpha(ig),symtrans(1,ig),symtrans(2,ig),symtrans(3,ig),&
             !     &params
          endif
          !
       end if!irc get line
    end do
    nsymop = ig
    !
    write(iunbin,'("[x-csmg-cgauss] Nrec=",I4)') nsymop
    do ig=1,nsymop
       params = 'F'
       if(syminv(ig)) then
          params = 'T'
       end if
       write(iunbin,'(1X,A8,1X,7(1X,F12.4),1X,A1)') &
            &symopstr(ig),symaxe(1,ig),symaxe(2,ig),symaxe(3,ig),&
            &symalpha(ig),symtrans(1,ig),symtrans(2,ig),symtrans(3,ig),&
            &params

    end do
    !
    ! Build the symmetry matrix
    call csmgetrmtv()
    !
  end function csm_rdinput
!
!H
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function csmmod_getnat()
    !
    csmmod_getnat = nat
    !
  end function csmmod_getnat
!
!H
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function csmmod_getnsymop()
    !
    csmmod_getnsymop = nsymop
    !
  end function csmmod_getnsymop
!
!H
!H-----------------------------------------------------------------------------
!H
!
  subroutine csmgetrmtv
    !
    integer(FINT) :: iop
    real(FREAL) :: thr = 1.0E-8_FREAL
    !
    logical :: debug = .false.
    !
    symident = 0_FINT
    syminvidx = 0_FINT
    do iop=1,nsymop
       if(abs(symalpha(iop)).lt.thr) then
          if(symident.ne.0.and..not.syminv(iop)) then
             call message(MESWARN,'[csmgetrmtv] WARNING MORE THAN ONE IDENTITY')
          endif
          if(syminvidx.ne.0.and.syminv(iop)) then
             call message(MESWARN,'[csmgetrmtv] WARNING MORE THAN ONE INVERSION')
          endif
          !
          ! Build Identity or Inversion Matrix
          if(syminv(iop)) then
             syminvidx = iop
             symrmat(1:3,1:3,iop) = 0.0_FREAL
             symrmat(1,1,iop) = -1.0_FREAL
             symrmat(2,2,iop) = -1.0_FREAL
             symrmat(3,3,iop) = -1.0_FREAL
          else
             symident = iop
             symrmat(1:3,1:3,iop) = 0.0_FREAL
             symrmat(1,1,iop) = 1.0_FREAL
             symrmat(2,2,iop) = 1.0_FREAL
             symrmat(3,3,iop) = 1.0_FREAL
          endif
          if(debug) then
             if(symident.ne.1.and..not.syminv(iop)) then
                call message(MESWARN,&
                     &'[csmgetrmtv] WARNING IDENTITY NOT FIRST SYMOP!')
             endif
          end if
       else
          !print *,iop, symalpha(iop),symaxe(1:3,iop),syminv(iop),&
          !&symeuler(1,iop),symeuler(2,iop),symeuler(3,iop),&
          !&symrmat(1:3,1:3,iop),thr
          !stop
          call csmmod_euler(symalpha(iop),symaxe(1:3,iop),syminv(iop),&
               &symeuler(1,iop),symeuler(2,iop),symeuler(3,iop),&
               &symrmat(1:3,1:3,iop),thr)
       endif
       if(debug) then
          write(0,'(A10)') symopstr(iop)(1:10)
          write(0,'("Transla. Vector")')
          write(0,'(3(1X,F12.4))') symtrans(:,iop)
!          write(0,'(1(1X,F12.4))') symeuler(1,iop)
!          write(0,'(1(1X,F12.4))') symeuler(2,iop)
!          write(0,'(1(1X,F12.4))') symeuler(3,iop)
!          write(0,'(3(1X,F12.4))') symeuler(:,iop)
          write(0,'("Rotation Vector")')
          write(0,'(3(1X,F12.4))') symaxe(:,iop)
          write(0,'("Rotation Matrix")')
          write(0,'(3(1X,F12.4))') symrmat(:,:,iop)
       endif
    enddo
    !
    !
  end subroutine csmgetrmtv
!
!
!H
!H----------------------------------------------------------------------
!H subroutine csmmod_euler(angle,rotax,inv,alpha,beta,gamma,rmat,leps)
!H----------------------------------------------------------------------
!H
!
  subroutine csmmod_euler(rangle, rotax, inv, alpha, beta, gamma, rmat,leps)
    !
    !H----------------------------------------------------------------------
    !H
    !H    calculates the euler angles 'alpha,beta,gamma' and the 3-d
    !H    rotation matrix 'rmat' corresponding to a rotation by 'angle'
    !H    radians about the rotation axis through point 'x,y,z' and
    !H    the origin. if 'inv'=.true. , an inversion is performed
    !H    following the pure rotation.
    !H
    !H    alpha = 1st rotation: about the z-axis,    -pi.lt.alpha.le.pi
    !H    beta  = 2nd rotation: about the y'-axis,    0 .le.beta .le.pi
    !H    gamma = 3rd rotation: about the z"-axis,   -pi.lt.gamma.le.pi
    !H
    !H    alpha = 1st rotation: about the z-axis,    -pi.< .alpha.<=.pi
    !H    beta  = 2nd rotation: about the y'-axis,    0 .<=.beta .<=.pi
    !H    gamma = 3rd rotation: about the z"-axis,   -pi.< .gamma.<=.pi
    !H
    !H    all rotations are passive (rotations of the coord system
    !H    rather than the object vectors). the coord system is assumed to
    !H    be right-handed, and rotations are positive if counterclockwise.
    !H    the axis vector 'rotax' is returned normalized and oriented so
    !H    that z is positive. (if z=0, so that y is pos; if z & y =0, so
    !H    that x is pos.) angle is returned with the appropriate sign,
    !H    within the range -pi.lt.angle.le.pi
    !H
    !H    the euler angles are extracted from the rotation matrix which
    !H    transforms z-axis into the rot'n axis, rotates by angle, and
    !H    then inverts the first transformation.
    !H
    !H    cf. rose,  p. 50ff, p.65
    !H
    !H    called by: symops,closur,molrot (Orignal CALL/ not valid here)
    !H
    !H----------------------------------------------------------------------
    !
    !W    F.Mariotti LEPS parameter in this case can be the machine precision
    !W    and not the symmetry precision!!
    !
    !
    !
    !
    logical, intent(in) :: inv
    real(FREAL), intent(in) :: rangle
    real(FREAL), intent(in) :: leps
    real(FREAL), intent(out) :: alpha
    real(FREAL), intent(out) :: beta
    real(FREAL), intent(out) :: gamma
    real(FREAL), dimension(3) :: rotax
    real(FREAL), dimension(3,3) :: rmat
    !
    !
    real(FREAL) :: angle
    !
    !     FROM INLINE CODE
    !
    !integer(FINT) :: i
    !integer(FINT) :: j
    integer(FINT) :: nmod
    real(FREAL) :: a
    real(FREAL) :: b
    real(FREAL) :: ca
    real(FREAL) :: cb
    real(FREAL) :: cosa
    real(FREAL) :: cosb
    real(FREAL) :: cosg
    real(FREAL) :: ct
    real(FREAL) :: rax
    real(FREAL) :: sa
    real(FREAL) :: sb
    real(FREAL) :: sina
    real(FREAL) :: sinb
    real(FREAL) :: sing
    real(FREAL) :: st
    real(FREAL) :: t
    real(FREAL) :: xd1
    real(FREAL) :: xd2
    real(FREAL) :: xn1
    real(FREAL) :: xn2
    real(FREAL) :: xyproj
    real(FREAL) :: x
    real(FREAL) :: y
    real(FREAL) :: z
    real(FREAL) :: pi
    real(FREAL) :: twopi
    real(FREAL) :: zero
    real(FREAL) :: one
    !
    !
    !
    angle=rangle
    pi = 4.0D0 * atan(1.0D0)
    twopi = 2.D0 * pi
    zero = 0.0_FREAL
    one = 1.0_FREAL
    !
    a = zero
    b = zero
    t = zero
    x = rotax(1)
    y = rotax(2)
    z = rotax(3)
    !
    !-----normalize the rotation axis vector ------------------------------
    !
    rax = sqrt(x * x + y * y + z * z)
    !
    !-----if no axis of rotation was specified, error stop ----------------
    !
    if (rax.lt.leps) then
       print *,'DEBUG Euler: undefined rotational axis.'
       stop 3
    endif
    !
    x = x / rax
    y = y / rax
    z = z / rax
    if (abs(x).lt.leps) x = zero
    if(abs(y).lt.leps) y = zero
    if (abs(z).lt.leps) z = zero
    !
    !-----bring angle within the range -pi to pi --------------------------
    !
    !print *,angle,twopi
    !call flush(6)
    nmod = int(angle/twopi)
    angle = angle-real(nmod,KIND=FREAL) * twopi
    if (angle.gt.pi) angle = -(twopi - angle)
    if (angle.lt. - pi) angle = twopi + angle
    !
    !-----direct the rotation axis vector into the upper
    !     hemisphere if it points toward negative z ---------------------
    !
    !
    !OLDCODE
    !      if (z.gt.zero)                             go to 10
    !      if (z.eq.zero.and.y.gt.zero)               go to 10
    !      if (z.eq.zero.and.y.eq.zero.and.x.gt.zero) go to 10
    !
    !      x=-x
    !      y=-y
    !      z=-z
    !      angle=-angle
    ! 10   if (abs(abs(angle)-pi).lt.leps) angle=pi
    !NEWLINES
    if (z.le.zero) then
       if (.not.(z.eq.zero.and.y.gt.zero) ) then
          if (.not.(z.eq.zero.and.y.eq.zero.and.x.gt.zero) ) then
             x = - x
             y = - y
             z = - z
             angle = - angle
          endif
       endif
    endif
    !
    if (abs(abs(angle)-pi).lt.leps) angle = pi
    !NEWLINES
    !
    !
    if (abs(angle).lt.leps) angle = zero
    xyproj = sqrt(x * x + y * y)
    rotax(1) = x
    rotax(2) = y
    rotax(3) = z
    !
    !-----find rot angles to transform z-axis into rotax,
    !     rotate by angle, & transform back: (ab)**-1 * t * (ab) --------
    !
    if (xyproj.ne.zero) then
       if (y.eq.zero.and.x.lt.zero) a = pi
       if (y.ne.zero) a = - atan2( - y, x)
       b = atan2(xyproj, z)
    endif
    !
    t = angle
    !
    !-----extract euler angles for the net rotation -----------------------
    !
    sa = sin(a)
    ca = cos(a)
    sb = sin(b)
    cb = cos(b)
    st = sin(t)
    ct = cos(t)
    beta = acos(one-(one-ct) * sb * sb)
    !
    !-----special cases: beta=zero or beta=pi ( sin(beta)=zero ) ----------
    !
    if (beta.lt.leps) then
       alpha = zero
       beta = zero
       gamma = angle
    elseif ((pi - beta) .lt.leps) then
       alpha = - atan2(x, y)
       beta = pi
       gamma = - alpha
    else
       xn1 = sa * cb *(one-ct)
       xn2 = ca * st
       xd1 = ca * cb *(one-ct)
       xd2 = sa * st
       alpha = atan2((xn1 - xn2),(xd1 + xd2))
       gamma = atan2((xn1 + xn2),(xd2 - xd1))
    endif
    !
    !-----------------------------------------------------------------------
    !     construct the 3-d rotation matrix for this set of euler angles
    !-----------------------------------------------------------------------
    !
    sina = sin(alpha)
    cosa = cos(alpha)
    sinb = sin(beta)
    cosb = cos(beta)
    sing = sin(gamma)
    cosg = cos(gamma)
    rmat(1, 1) = cosa * cosb * cosg - sina * sing
    rmat(1, 2) = sina * cosb * cosg + cosa * sing
    rmat(1, 3) = - sinb * cosg
    rmat(2, 1) = - cosa * cosb * sing - sina * cosg
    rmat(2, 2) = - sina * cosb * sing + cosa * cosg
    rmat(2, 3) = sinb * sing
    rmat(3, 1) = cosa * sinb
    rmat(3, 2) = sina * sinb
    rmat(3, 3) = cosb
    !
    !-----if this is an improper rotation, invert through the origin ------
    !
    if (inv) then
       rmat = rmat * (-1.0D0)
    endif
    !
    !
    !
    return
  end subroutine csmmod_euler
  !
  !

!
!H
!H-----------------------------------------------------------------------------
!H
!
  subroutine csmmod_getrho(vint)
    !
    real(FREAL), dimension(:), intent(out) :: vint
    !
    integer(FINT) :: i
    !
    if(size(vint).lt.nat) then
       call message(MESERRO,'[csmmod_getrho] Internal inconsistency.')
       call message(MESKILL,'[csmmod_getrho] Kill Request.')
    end if
    !
    do i=1,nat
       !
       vint(i) = atg_pars(i,1)*(sqrt(pi/(atg_pars(i,2)))**3)
       !
    end do
    !
  end subroutine csmmod_getrho
!
!H
!H-----------------------------------------------------------------------------
!H
!
  subroutine csmmod_getrrii(vint)
    !
    real(FREAL), dimension(:), intent(out) :: vint
    !
    integer(FINT) :: i
    !
    if(size(vint).lt.nat) then
       call message(MESERRO,'[csmmod_getrrii] Internal inconsistency.')
       call message(MESKILL,'[csmmod_getrrii] Kill Request.')
    end if
    !
    do i=1,nat
       !
       vint(i) = atg_pars(i,1)*atg_pars(i,1)&
            &*(sqrt(pi/(atg_pars(i,2)+atg_pars(i,2)))**3)
       !
    end do
    !
  end subroutine csmmod_getrrii
!
!H
!H-----------------------------------------------------------------------------
!H
!
  subroutine csmmod_getrrij(vint)
    !
    real(FREAL), dimension(:), intent(out) :: vint
    !
    integer(FINT) :: i,j,k
    real(FREAL) :: r
    !
    k = 0
    do i=1,nat
       do j=(i+1),nat
          k=k+1
          r = sqrt(&
               &(molecule_xyz(1,i)-molecule_xyz(1,j))**2 +&
               &(molecule_xyz(2,i)-molecule_xyz(2,j))**2 +&
               &(molecule_xyz(3,i)-molecule_xyz(3,j))**2)
          !
          vint(k) = atg_pars(i,1)*atg_pars(j,1)&
               &*(sqrt(pi/(atg_pars(i,2)+atg_pars(j,2)))**3)
          vint(k) = vint(k) * exp(&
               &-((atg_pars(i,2)*atg_pars(j,2))/(atg_pars(i,2)+atg_pars(j,2)))&
               & * (abs(r)**2))
          !vint(k) = atg_pars(i,1)*atg_pars(j,1)&
          !     &*sqrt(pi/(atg_pars(i,2)+atg_pars(j,2)))**3&
          !     &*exp(&
          !     &-((atg_pars(i,2)*atg_pars(j,2))/(atg_pars(i,2)+atg_pars(j,2)))&
          !     &*abs(r)**2)
          !write(*,'("GIJ",1X,2(1X,I4),2(1X,F16.9))') i,j,vint(k),r
       end do
    end do
    
    !
  end subroutine csmmod_getrrij
!
!H
!H-----------------------------------------------------------------------------
!H
!
  subroutine csmmod_getimj(mint)
    !
    real(FREAL), dimension(:,:,:), intent(out) :: mint
    !
    integer(FINT) :: i,j,m
    real(FREAL), dimension(3) :: lxyz
    real(FREAL) :: r
    !
    do m=1,nsymop
       !print *,'SymOp ',m
       do j=1,nat
          lxyz = matmul(symrmat(1:3,1:3,m),molecule_xyz(1:3,j))
          lxyz(1:3) = lxyz(1:3) + symtrans(1:3,m)
          do i=1,nat
             r = sqrt(&
                  &(molecule_xyz(1,i)-lxyz(1))**2 +&
                  &(molecule_xyz(2,i)-lxyz(2))**2 +&
                  &(molecule_xyz(3,i)-lxyz(3))**2)
             !print *,'Running ',m,i,j
             if (abs(r).lt.1.0D-6)then
                !print *,'Closed ',i,j
                mint(i,m,j) = atg_pars(i,1)*atg_pars(j,1)&
                     &*sqrt(pi/(atg_pars(i,2)+atg_pars(j,2)))**3
             else
                mint(i,m,j) = atg_pars(i,1)*atg_pars(j,1)&
                     &*sqrt(pi/(atg_pars(i,2)+atg_pars(j,2)))**3&
                     &*exp(&
                     &-((atg_pars(i,2)*atg_pars(j,2))/(atg_pars(i,2)+atg_pars(j,2)))&
                     &*abs(&
                     &sqrt(&
                     &(molecule_xyz(1,i)-lxyz(1))**2 +&
                     &(molecule_xyz(2,i)-lxyz(2))**2 +&
                     &(molecule_xyz(3,i)-lxyz(3))**2)&
                     &)**2)
             end if
             !write (*,'("IMJ",1X,3(1x,I4),1X,4(1X,F16.9))')&
             !     & i,j,m,mint(i,m,j),lxyz
          end do
       end do
    end do
    !
  end subroutine csmmod_getimj
!
!H
!H-----------------------------------------------------------------------------
!H
!
  subroutine csmmod_getimpj(mint)
    !
    real(FREAL), dimension(:,:,:), intent(out) :: mint
    !
    integer(FINT) :: i,j,m
    real(FREAL), dimension(3) :: lxyz
    real(FREAL) :: r
    !
    mint(:,:,:)=0.0_FREAL
    do m=2,nsymop
       !print *,'SymOp ',m
       do j=1,nat
          lxyz = matmul(symrmat(1:3,1:3,m),molecule_xyz(1:3,j))
          lxyz(1:3) = lxyz(1:3) + symtrans(1:3,m)
          do i=1,nat
             r = sqrt(&
                  &(molecule_xyz(1,i)-lxyz(1))**2 +&
                  &(molecule_xyz(2,i)-lxyz(2))**2 +&
                  &(molecule_xyz(3,i)-lxyz(3))**2)
             !print *,'Running ',m,i,j
             if (abs(r).lt.1.0D-6)then
                !print *,'Closed ',i,j
                mint(i,m,j) = atg_pars(i,1)*atg_pars(j,1)&
                     &*sqrt(pi/(atg_pars(i,2)+atg_pars(j,2)))**3
             else
                mint(i,m,j) = atg_pars(i,1)*atg_pars(j,1)&
                     &*sqrt(pi/(atg_pars(i,2)+atg_pars(j,2)))**3&
                     &*exp(&
                     &-((atg_pars(i,2)*atg_pars(j,2))/(atg_pars(i,2)+atg_pars(j,2)))&
                     &*abs(&
                     &sqrt(&
                     &(molecule_xyz(1,i)-lxyz(1))**2 +&
                     &(molecule_xyz(2,i)-lxyz(2))**2 +&
                     &(molecule_xyz(3,i)-lxyz(3))**2)&
                     &)**2)
             end if
             !write (*,'("IMJ",1X,3(1x,I4),1X,4(1X,F16.9))')&
             !     & i,j,m,mint(i,m,j),lxyz
          end do
       end do
    end do
    !
  end subroutine csmmod_getimpj
!
!H
!H-----------------------------------------------------------------------------
!H
!

!
!H
!H-----------------------------------------------------------------------------
!H subroutine csmmod_rotxyz(tx,ty,tz,rx,ry,rz,xyz)
!H-----------------------------------------------------------------------------
!H
!
  subroutine csmmod_rotxyz(tx,ty,tz,rx,ry,rz,xyz)
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
    real(FREAL), intent(in) :: tx,ty,tz
    real(FREAL), intent(inout) :: rx,ry,rz
    real(FREAL), dimension(:,:) :: xyz
    !
    integer(FINT) :: i,nat
    real(FREAL), dimension(3,3) :: rmat
    real(FREAL), dimension(3,3) :: amat,bmat,gmat
    real(FREAL), dimension(3) :: taxe
    real(FREAL) :: eps = 1.0E-12_FREAL
    real(FREAL) :: alpha,beta,gamma
    !
    ! Transform in radians
    !---------------------
    rx = rx * (4.0_KREAL*atan(1.0_KREAL)) / 180.0_KREAL
    ry = ry * (4.0_KREAL*atan(1.0_KREAL)) / 180.0_KREAL
    rz = rz * (4.0_KREAL*atan(1.0_KREAL)) / 180.0_KREAL
    !
    ! Translate First
    !----------------
    xyz(1,:) = tx + xyz(1,:)
    xyz(2,:) = ty + xyz(2,:)
    xyz(3,:) = tz + xyz(3,:)
    !
    ! Get rotational matrixes
    !------------------------
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
    ! Some debugging
    !---------------
    if(.false.) then
       write(0,'("Required Rotation")')
       write(0,'("Translation Vector:")')
       write(0,'(3(1X,F12.4))') tx,ty,tz
       write(0,'("Euler Vector:")')
       write(0,'(3(1X,F12.4))') rx,ry,rz
       write(0,'("X Matrix:")')
       write(0,'(3(1X,F12.4))') amat(:,:)
       write(0,'("Y Matrix:")')
       write(0,'(3(1X,F12.4))') bmat(:,:)
       write(0,'("Z Matrix:")')
       write(0,'(3(1X,F12.4))') gmat(:,:)
       write(0,'("Final Rotation Matrix")')
       write(0,'(3(1X,F12.4))') rmat(:,:)
       call molecule_print(0,'mldfromau')
    endif
    !
    ! Apply rotation
    !---------------
    nat = molecule_getnat()
    do i=1,nat
       xyz(1:3,i) = matmul(rmat(1:3,1:3),xyz(1:3,i))
    end do
    !
    !
    !
  end subroutine csmmod_rotxyz

!
!H
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function csmmod_setrdminuit(iunin)
    !
    integer(FINT), intent(in) :: iunin
    character(FLCHARS) :: params
    integer(FINT) :: irc
    !
    csmmod_setrdminuit = iunin
    rewind(iunin)
    irc = osec_set(iunin,'x-minuit',params)
    if(irc.lt.0) then
       csmmod_setrdminuit = -1
       csm_error = 1
       call message(MESERRO,'[csmmod_setrdminuit] No minuit section.')
       return
    end if
    !
  end function csmmod_setrdminuit
!
!H
!H-----------------------------------------------------------------------------
!H subroutine csmmod_printhelp(iun,what)
!H-----------------------------------------------------------------------------
!H
!
  subroutine csmmod_printhelp(iun,what)
!
!H
!H-----------------------------------------------------------------------------
!H
!H Print Help on iun. There is an optional 'what'.
!H what at present can have to possibilities: 'S' and 'E' for short and
!H extended help. Possibly in the future there will be an 'M' for manual.
!H It defaults to 'S'
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    integer(FINT), intent(in), optional :: iun
    character(1), intent(in), optional :: what
    !
    integer(FINT) :: liun
    character(FLCHARS)  :: lwhat
    !
    if(present(iun)) then
       liun = iun
    else
       liun = 6
    end if
    !
    if(present(what)) then
       lwhat = what
    else
       lwhat = 'S'
    end if
    !
    write(*,*) ""
    write(*,*) " CSMD - Continuous Symmetry Measure with electron Density"
    write(*,*) " Usage:"
    write(*,*) " CSMD [-i input] [-o output] [-e errorfile] [-l logfile] [-p noshift]"
    write(*,*) "      [-p onestep] [-p useminuit] [-p readminuit] [-p interminuit]"
    write(*,*) " CSMD -h"
    write(*,*) " CSMD -h full"
    if (lwhat.eq.'E') then
    write(*,*) ""
    write(*,*) " -i F        Set input file name to F"
    write(*,*) " -o F        Set output file name to F"
    write(*,*) " -e F        Set error file name to F"
    write(*,*) " -l F        Set log file name to F"
    write(*,*) " -p S        Set option to string S. It can (and should) be present more then"
    write(*,*) "             once in the command line."
    write(*,*) "             Available [S] options are:"
    write(*,*) "             noinertia"
    write(*,*) "                        Do not rotate the initial coordinates to diagonalise"
    write(*,*) "                        inertia tensor. Meaningful only if noshift is not set."
    write(*,*) "             noshift"
    write(*,*) "                        Do not translate the molecule to have the center of"
    write(*,*) "                        mass at the origin of coordinates. Imply noinertia."
    write(*,*) "             onestep"
    write(*,*) "                        Do not perform minimization."
    write(*,*) "                        Only one initial step of CSM."
    write(*,*) "             useminuit"
    write(*,*) "                        Use minuit as optimized."
    write(*,*) "                        There are some bugs around so use it and test."
    write(*,*) "             readminuit"
    write(*,*) "                        Read minuit commands from the input file in the section"
    write(*,*) "                        [x-minuit]. Refer to minuit (version f77) for available"
    write(*,*) "                        commands."
    write(*,*) "             interminuit"
    write(*,*) "                        Start an interactive minuit session to fine tune your"
    write(*,*) "                        optimization."
    write(*,*) ""
    else if (lwhat.eq.'M') then
    write(*,*) ""
    write(*,*) " Print Manual not yet implemented. Please run: CSMD -h full"
    write(*,*) ""
    end if
    !
    return
    !
  end subroutine csmmod_printhelp
  !
  !
end module csmmod

