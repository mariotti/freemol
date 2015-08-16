!
!H
!H-----------------------------------------------------------------------------
!H This File is part of the Frimol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H MODULE extio: Frimol by F.Mariotti
!H-----------------------------------------------------------------------------
!H $Id: extio.F90,v 1.1.1.1 2009/01/12 16:56:17 mariotti Exp $
!H-----------------------------------------------------------------------------
!H 
!H This module perform controlled io.
!H It is a part of the file manager modules system.
!H File manager moldule system is composed by:
!H basic I/O: baseio, Extended File manager: extio, Section Manager: sections.
!H 
!H This module define some routines to perform i/o using a new file type.
!H 
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H In this module I would like to implement:
!H OPEN statment that returns extfile type:
!H open()               Assume scratch binary file for local storage
!H open(name)           Assume ascii file for reading/overwriting = UNKNOWN
!H                      TODO: in the future: implement autodetection of
!H                            the url string format.
!H open(url)            Same as before but with url.
!H                      At present we should use an open(URL=STR)
!H open(XX,FORMAT=YY,STATUS=ZZ,ACCESS=QQ)
!H           Where:
!H           XX define the three different way to open the file
!H           as above. (note URL)
!H 
!H           YY is a string in the following list:
!H           ASCII      Open an STD ascii file.
!H           BINARY     Open a STD binary file.
!H           COMPRESSED Open a compressed file (bzip2,gzip).
!H           BYTECODE   Open a byte coded file
!H           BASE64     Open a base64 file.
!H           PORTABLE   Open with default definition of portab.
!H           DEFAULT    Open default selected.
!H                      TODO: Which one will be default?
!H           
!H           ZZ is a string defining the status as:
!H           INPUT      File is the input (same as read).
!H           OUPUT      File os the output (same as write but it is local
!H                      and we try to flush it in order to garantee
!H                      as much as possible a match between the status
!H                      of the program and the status of the output).
!H           READ       File must exist and is opened for reading.
!H           WRITE      File must NOT exist and is opened for writing.
!H           APPEND     File is open for append it can exist or not.
!H           UNKNOWN    File can be readed or written in (NOT APPEND!)
!H           SCRATCH    Used for local storage it does not exist
!H                      a new unit is always created.
!H           REQUEST    It is quite new as a concept: We perform
!H                      a request in order to get the current status
!H                      of the file. Right now the request goes
!H                      thrue the identification by name or url or
!H                      (see later) via the current status of an
!H                      existing unit of extfile type.
!H                      The possibility to define a protocol of comunication
!H                      with a "file server" is anyway open.
!H                      It can be used for flow control.
!H           
!H           QQ is a string defining the ACCESS (Today STD FORTRAN):
!H           DIRECT
!H           SEQUENTIAL
!H           
!H
!H same as above but return the iun (anable to get out of the extfile
!H                                   struct: deprecated).
!H                                   IUN is returned as procedure parameter
!H                                   This routine still return an extfile type.
!H open(iun)            Assume scratch binary file for local storage
!H open(iun,name)       Assume ascii file for reading/overwriting = UNKNOWN
!H                      TODO: in the future: implement autodetection of
!H                            the url string format.
!H open(iun,url)        Same as before but with url.
!H
!H open(iun,XX,FORMAT=YY,STATUS=ZZ,ACCESS=QQ)
!H                      Same as the description above but with one main
!H                      difference: If REQUEST is given we check the table
!H                      of the opened units against the iun value. In this
!H                      case iun is used as input and not as output parameter.
!H                      To force iun to be ignored and to use the routine
!H                      as in the previous examples set iun to -1 before the
!H                      call (deprecated).
!H
!H The same set of routines can be called as well with an extfile type
!H instead of the iun parameters i.e.:
!H open(extfile)
!H In this case the procedure return an integer error code.
!H
!H strings given as name or url are always processed to check trailing spaces.
!H
!H To summarize:
!H the open procedure can be called as:
!H
!H extfile = open() Open Scratch
!H extfile = open(iun) Open Scratch
!H extfile = open(extfile) Open extfile dependent
!H extfile = open(params) Open params dependent
!H extfile = open(iun,params) Open params dependent
!H extfile = open(extfile,params) Open params dependent
!H WARNING!! Changed.. the extfile version require a pointer to extfile
!H           type to be given
!H iun = open() Open Scratch
!H iun = open(iun) Open Scratch
!H iun = open(extfile) Open extfile dependent
!H iun = open(params) Open params dependent
!H iun = open(iun,params) Open params dependent
!H iun = open(extfile,params) Open params dependent
!H
!H With this reminding:
!H The routines with integer return code will return errors via the return
!H code: or in other words check if the iun return code is .lt. 0 to catch
!H errors.
!H The error code in the routines with an extfile as return code is set
!H as property: extfile%error.
!H
!H TODO: follow.....
!H It is much more what we really miss!!!!! ...
!H
!H I miss:
!H A connection by the meaning of a unique link between fortran units
!H and the new extfile... functions like: extfile = enquire(unit).
!H These kind of functions are part of more extended INQUIRE project,
!H or are part of "less coded" extension of extfile control which
!H will have to be introduced.
!H Starting from the "simple" case of extfiles we would like to inquire
!H about the current state of the file.
!H This means also that we will have to introduce a file lock facility
!H at this level.
!H
!H We want to shift to a much more OOP environment, still keeping fortan
!H standards.
!H so in the fututre we will move to this conceptand procedure:
!H
!H 1) I ask for a new file handle
!H              The (ext)file hanlde is stored for updates(and we will keep
!H              keep track of what's going on).
!H 2) I get a file handle with its defaults.
!H 3) I modify the default as needed.
!H 4) I ask to open the file(handle) with the current properties.
!H 5) I get the file handle.
!H
!H Because we are not java and we do not need to mess up all the concept
!H We need simple cases to be typed at once and if possible
!H Specifing the use we are going to do of the file.
!H It means: we will have to hide the 1-5 procedure.
!H
!H our OOP concept will reduce to:
!H I need a file for input/output
!H I need a file for fast access of working data
!H I need a file for temporary storage which will be usefull for restart.
!H
!H and our OOP concept will reduce further:
!H
!H whenever I need a fast access file (SCRATCH) I will try to 
!H be strict to the fortran file concept and get the best out of it.
!H
!H here come the next abstraction we will do: the section concept.
!H
!H so generally speacking: Why all of this?
!H
!H We built a baseio module which has to deal with
!H open and close statment and furnish us of a fortran valid unit.
!H The baseio module has not to interfer with other problems.
!H The baseio module deals at machine level and will implement
!H "methods" at this level "only". It deal with pwrformance, so it
!H just handle fortran units.
!H
!H-----------------------------------------------------------------------------
!H We try do do the same as baseio:
!H extio_xxxx public procedures
!H eio_xxxx   private procedures
!H 
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H 
!
module extio
  !
  use vartypes
  use linetools
  use messages
  use baseio
  !
  implicit none
  private
  !
  ! The Public extfile TYPE
  public extfile
  !
  type extfile
     ! The STD Fortran Unit
     !---------------------
     integer(FINT) :: unit
     character(FLCHARS) :: name          ! The local STD Name
     character(FLCHARS) :: url           ! If Available the name as URL
     logical :: issection                ! Is A Frimol Section File
     logical :: isurl                    ! Do we have an URL
     logical :: isascii                  ! Is Binary or Ascii
     character(FLCHARS) :: format        ! The STD Fortran format
     character(FLCHARS) :: status        ! The STD Fortran status
     character(FLCHARS) :: access        ! The STD Fortran access
     integer(FINT) :: error              ! The error code if any
     character(FLCHARS) :: error_message ! The error message if any
  end type extfile
  !
  ! Internal Error Handler
  integer(FINT), save, private :: eio_error = 0
  !
  !Is the module initialized? up to now always initialized
  logical,save, private :: eio_isinit = .true.
  !
  public extio_init
  public extio_open
  public extio_geterror
  !
  interface extio_open
     !
     module procedure extio_open_null !extfile = open()
     module procedure extio_open_null_i !iun = open()
     !
     module procedure extio_open_iun !extfile = open(iun)
     module procedure extio_open_iun_i !iun = open(iun)
     !
!     module procedure extio_open_efh !True open routine. extfile = open(extfile)
!     module procedure extio_open_efh_i !iun = open(extfile)
     !
     module procedure extio_open_nullp !extfile = open(params)
     module procedure extio_open_null_ip !iun = open(params)
     !
     module procedure extio_open_iunp !extfile = open(iun,params)
     module procedure extio_open_iun_ip !iun = open(iun,params)
     !
!     module procedure extio_open_efhp !extrfile = open(extfile,params)
!     module procedure extio_open_efh_ip !iun = open(extfile,params)
     !
  end interface !extio_open
  !
  !
  !
  ! Error Interface
  !----------------
  ! We define two error interfaces. Each one return the latest error
  ! code or string value. The String error function can be used to inqure
  ! about a given error code.
  !----------------------------------------------------------------------
  interface extio_geterror
     module procedure extio_geterror_i
     module procedure extio_geterror_c
  end interface !extio_geterror
  !
contains
  !
!
!H
!H-----------------------------------------------------------------------------
!H SOME RETURN VALUE FUNCTIONS
!H-----------------------------------------------------------------------------
!H 
!
!
!H
!H-----------------------------------------------------------------------------
!H logical function extio_isinit()
!H-----------------------------------------------------------------------------
!H 
!
  logical function extio_isinit()
    extio_isinit = eio_isinit
  end function extio_isinit
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function extio_geterror_i()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function extio_geterror_i()
!
!H
!H-----------------------------------------------------------------------------
!H 
!H Return the error code of the module.
!H It is a part of a set of routine to return errors.
!H See the discussion in the interface.
!H 
!H-----------------------------------------------------------------------------
!H
!
    ! We merely return the integer current error code
    !------------------------------------------------
    extio_geterror_i = eio_error
    !
  end function extio_geterror_i
!
!H-----------------------------------------------------------------------------
!H
!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H character(FLCHARS) function extio_geterror_c(code)
!H-----------------------------------------------------------------------------
!H
!
  character(FLCHARS) function extio_geterror_c(code)
    !
!    integer(FINT), optional, intent(in) :: code
! This has been changed for IFC compiler which was getting confused...
! to this:
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
    integer(FINT) :: lcode !the local used error code
    !
    lcode = eio_error
!    if(present(code)) then
       lcode = code
!    endif
    !
    ! We return the error code
    !-------------------------
    ! Note that we have not yet implement it
    ! so we report as example case 2_FINT as unwanted code
    ! but it has to disappear! And Default is a not reported
    ! error code.
    !-------------------------------------------------------
    select case(code)
    case(0_FINT)
       extio_geterror_c = "No Errors."
    case(1_FINT)
       extio_geterror_c = "Not documented error or generic error.."
    case(2_FINT)
       extio_geterror_c = "Double Error: Why did you get this message?."
    case default
       extio_geterror_c = "Not documented error or generic error.."
    end select
    !
    return
    !
  end function extio_geterror_c

!
!H
!H-----------------------------------------------------------------------------
!H END OF: SOME RETURN VALUE FUNCTIONS
!H-----------------------------------------------------------------------------
!H 
!
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function extio_init()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function extio_init()
    !
    integer(FINT) :: irc
!
!H
!H-----------------------------------------------------------------------------
!H Initialization function:
!H Set the eio_isinit flag and call the baseio initialization.
!H Perform as well some resets.
!H-----------------------------------------------------------------------------
!H
!
    irc = baseio_init()
    if(irc.lt.0) then
       eio_error = 1
       eio_isinit = .false.
       extio_init = -1
       call message(MESDEBG,"[extio_init]: cannot initialize baseio.")
       return
    endif
    !
    ! There are not errors in baseio initialization
    ! so we return the true...
    !-------------------------
    eio_error = 0
    eio_isinit = .true.
    extio_init = 0
    return
    !
  end function extio_init
  !
  !
  !
!    type(extfile) function extio_openefh(file,status,frmt,acc) 
!      !
!      character(FLCHARS), intent(in), optional :: file
!      character(FLCHARS), intent(in), optional :: status
!      character(FLCHARS), intent(in), optional :: frmt
!      character(FLCHARS), intent(in), optional :: acc
!      !
!      type(extfile) :: lefh
!      integer(FINT) :: irc
!      !
!      irc = extio_open(lefh,file,status,frmt,acc)
!      if(irc.lt.0) then
!         call message(MESERRO,'[extio_open]: cannot open the file.')
!      end if
!      extio_openefh = lefh
!      !
!   end function extio_openefh
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function extio_open_null_i()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function extio_open_null_i()
!
!H
!H-----------------------------------------------------------------------------
!H So in this routine we assume a call like:
!H irc/iun = extio_open()
!H which intrinsically mean to open e temporary scratch file
!H where we ask for the fortran unit.
!H irc/iun is the fortran unit unless .lt. 0 in which case is an error code.
!H As before there is no meaning to use this set of function to speed up
!H so we use the extfile type to construct a file type and than we will
!H destroy it. This function is present for compatibility. The advantage
!H of using this module comes whem we will be able to allocate the spece in the
!H "best" position: Fast, enought big.
!H
!H-----------------------------------------------------------------------------
!H
!
    !Locals
    !------
    type(extfile) :: efh
    !
    ! We set our defaults
    !--------------------
    ! Note we leave the other the other parameters unchanged
    ! right now for future use.
    !--------------------------
    efh%isascii = .false.
    efh%format = 'UNFORMATTED'
    efh%status = 'SCRATCH'
    efh%access = 'SEQUENTIAL'
    ! We reset the error label
    !-------------------------
    efh%error = 0
    efh%error_message = ''
    extio_open_null_i = 0
    !
    efh = extio_open_efh(efh)
    !
    ! Have we got any error?
    !-----------------------
    if(efh%error.gt.0) then
       extio_open_null_i = -1
       return
    endif
    !
    extio_open_null_i = efh%unit
    return
    !
  end function extio_open_null_i
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function extio_open_iun_i(iun)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function extio_open_iun_i(iun)
    !
    integer(FINT), intent(out) :: iun
!
!H
!H-----------------------------------------------------------------------------
!H So in this routine we assume a call like:
!H irc = extio_open(iun)
!H which intrinsically mean to open e temporary scratch file
!H where we ask for the fortran unit.
!H irc/iun is the fortran unit unless .lt. 0 in which case is an error code.
!H As before there is no meaning to use this set of function to speed up
!H so we use the extfile type to construct a file type and than we will
!H destroy it. This function is present for compatibility. The advantage
!H of using this module comes whem we will be able to allocate the spece in the
!H "best" position: Fast, enought big.
!H The input/output iun parameter is the iun/irc as well.
!H
!H-----------------------------------------------------------------------------
!H
!
    !Locals
    !------
    type(extfile) :: efh
    !
    ! We set our defaults
    !--------------------
    ! Note we leave the other the other parameters unchanged
    ! right now for future use.
    !--------------------------
    efh%isascii = .false.
    efh%format = 'UNFORMATTED'
    efh%status = 'SCRATCH'
    efh%access = 'SEQUENTIAL'
    ! We reset the error label
    !-------------------------
    efh%error = 0
    efh%error_message = ''
    extio_open_iun_i = 0
    !
    efh = extio_open_efh(efh)
    !
    ! Have we got any error?
    !-----------------------
    if(efh%error.gt.0) then
       extio_open_iun_i = -1
       iun = -1
       return
    endif
    !
    extio_open_iun_i = efh%unit
    iun = efh%unit
    return
    !
    !
  end function extio_open_iun_i
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function extio_open_efh_i(efh)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function extio_open_efh_i(efh)
    !
    type(extfile), intent(inout) :: efh
    !XXXXXXXXXXXX
    extio_open_efh_i = 0
    return
    !
  end function extio_open_efh_i
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function extio_open_null_ip(status,format,access)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function extio_open_null_ip(status,format,access)
    !
    !character(FLCHARS), optional :: status
    ! Force the status to be given or the we cannot distinguish
    ! from the empty version..........without _i(p)
    character(FLCHARS) :: status
    character(FLCHARS), optional :: format
    character(FLCHARS), optional :: access
    !
    !XXXXXXXXXXXX
    extio_open_null_ip = 0
    !
    return
    !
  end function extio_open_null_ip
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function extio_open_iun_ip(status,format,access)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function extio_open_iun_ip(iun,status,format,access)
    !
    integer(FINT), intent(out) :: iun
    character(FLCHARS) :: status
    character(FLCHARS), optional :: format
    character(FLCHARS), optional :: access
    !
    !XXXXXXXXXXXX
    !
    extio_open_iun_ip = 0
    !
    return
    !
  end function extio_open_iun_ip
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function extio_open_efh_ip(efh,status,format,access)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function extio_open_efh_ip(efh,status,format,access)
    !
    type(extfile), intent(inout) :: efh
    character(FLCHARS), optional :: status
    character(FLCHARS), optional :: format
    character(FLCHARS), optional :: access

    !XXXXXXXXXXXX
    !
    extio_open_efh_ip = 0
    !
    return
    !
  end function extio_open_efh_ip
!
!H
!H-----------------------------------------------------------------------------
!H  type(extfile) function extio_open_null()
!H-----------------------------------------------------------------------------
!H
!
  type(extfile) function extio_open_null(efh)
    !
!
!H
!H-----------------------------------------------------------------------------
!H So in this routine we assume a call like:
!H extfile = extio_open()
!H which intrinsically mean to open e temporary scratch file
!H where we ask for the fortran unit.
!H As before there is no meaning to use this set of function to speed up
!H so we use the extfile type to construct a file type and than we will
!H destroy it. This function is present for compatibility. The advantage
!H of using this module comes whem we will be able to allocate the spece in the
!H "best" position: Fast, enought big.
!H We return a correct extfile handle.
!H
!H-----------------------------------------------------------------------------
!H
!
    !Locals
    !------
    type(extfile) :: efh
    !
    ! We set our defaults
    !--------------------
    ! Note we leave the other the other parameters unchanged
    ! right now for future use.
    !--------------------------
    efh%isascii = .false.
    efh%format = 'UNFORMATTED'
    efh%status = 'SCRATCH'
    efh%access = 'SEQUENTIAL'
    ! We reset the error label
    !-------------------------
    efh%error = 0
    efh%error_message = ''
    !
    efh = extio_open_efh(efh)
    extio_open_null = efh
    !
    ! We don't need to check for error
    ! we suppose extio_open_efh does it for us
    ! so basically we just return the error or the extfile
    !-----------------------------------------------------
    !
    return
    !
  end function extio_open_null
!
!H
!H-----------------------------------------------------------------------------
!H  type(extfile) function extio_open_iun(iun)
!H-----------------------------------------------------------------------------
!H
!
  type(extfile) function extio_open_iun(efh,iun)
    !
    integer(FINT), intent(out) :: iun
    type(extfile) :: efh
    !
!
!H
!H-----------------------------------------------------------------------------
!H So in this routine we assume a call like:
!H extfile = extio_open(iun)
!H which intrinsically mean to open e temporary scratch file
!H where we ask for the fortran unit.
!H As before there is no meaning to use this set of function to speed up
!H so we use the extfile type to construct a file type and than we will
!H destroy it. This function is present for compatibility. The advantage
!H of using this module comes whem we will be able to allocate the spece in the
!H "best" position: Fast, enought big.
!H We return a correct extfile handle.
!H iun is retuner as the correct fortran unit.
!H
!H-----------------------------------------------------------------------------
!H
!
    !Locals
    !------
    !
    ! We set our defaults
    !--------------------
    ! Note we leave the other the other parameters unchanged
    ! right now for future use.
    !--------------------------
    efh%isascii = .false.
    efh%format = 'UNFORMATTED'
    efh%status = 'SCRATCH'
    efh%access = 'SEQUENTIAL'
    ! We reset the error label
    !-------------------------
    efh%error = 0
    efh%error_message = ''
    !
    efh = extio_open_efh(efh)
    extio_open_iun = efh
    !
    ! We don't need to check for error
    ! we suppose extio_open_efh does it for us
    ! so basically we just return the error or the extfile.
    ! We just set the fortran unit.
    !-----------------------------------------------------
    iun = efh%unit
    !
    return
    !
  end function extio_open_iun
!
!H
!H-----------------------------------------------------------------------------
!H  type(extfile) function extio_open_nullp(name,status,format,access)
!H-----------------------------------------------------------------------------
!H
!
  type(extfile) function extio_open_nullp(efh,name,status,format,access)
!
!H
!H-----------------------------------------------------------------------------
!H 
!H We are likely to go to open a fully described file.
!H We have the name so it not (in principle) for temporary storage and
!H anyway we get the details via parameters. We just have to check if the
!H parameters are there.
!H The dafault case as stated before is an ASCII file.
!H 
!H 
!H-----------------------------------------------------------------------------
!H
!
    !
    character(FLCHARS) :: name
    type(extfile) :: efh
    character(FLCHARS), optional :: status
    character(FLCHARS), optional :: format
    character(FLCHARS), optional :: access
    !
    !Locals
    !------
    !
    ! We set our defaults
    !--------------------
    ! Note we leave the other the other parameters unchanged
    ! right now for future use.
    !--------------------------
    efh%isascii = .true.
    efh%format = 'FORMATTED'
    efh%status = 'UNKNOWN'
    efh%access = 'SEQUENTIAL'
    ! We reset the error label
    !-------------------------
    efh%error = 0
    efh%error_message = ''
    ! We set the name
    !----------------
    efh%name = trim(name)
    !TODO: check the URL version!!
    !-----------------------------
    !
    ! We have to set the ither parameters.
    !-------------------------------------
    if(present(status)) then
       efh%status = trim(status)
    endif
    if(present(format)) then
       efh%format = trim(format)
    endif
    if(present(access)) then
       efh%access = trim(access)
    endif

    efh = extio_open_efh(efh)
    !
    ! Have we got any error?
    !-----------------------
    ! Let the return code to set it!
    !-------------------------------
    extio_open_nullp = efh
    return

    !
  end function extio_open_nullp
!
!H
!H-----------------------------------------------------------------------------
!H  type(extfile) function extio_open_iunp(status,format,access)
!H-----------------------------------------------------------------------------
!H
!
  type(extfile) function extio_open_iunp(efh,iun,status,format,access)
    !
    integer(FINT), intent(out) :: iun
    type(extfile) :: efh
    character(FLCHARS) :: name
    character(FLCHARS) :: status
    character(FLCHARS), optional :: format
    character(FLCHARS), optional :: access
!
!H
!H-----------------------------------------------------------------------------
!H 
!H We are likely to go to open a fully described file.
!H We have the name so it not (in principle) for temporary storage and
!H anyway we get the details via parameters. We just have to check if the
!H parameters are there.
!H The dafault case as stated before is an ASCII file.
!H 
!H 
!H-----------------------------------------------------------------------------
!H
!
    !
    !Locals
    !------
    !
    ! We set our defaults
    !--------------------
    ! Note we leave the other the other parameters unchanged
    ! right now for future use.
    !--------------------------
    efh%isascii = .true.
    efh%format = 'FORMATTED'
    efh%status = 'UNKNOWN'
    efh%access = 'SEQUENTIAL'
    ! We reset the error label
    !-------------------------
    efh%error = 0
    efh%error_message = ''
    ! We set the name
    !----------------
    name="x.unknown"
    efh%name = trim(name)
    stop "still developping. extio:864"
    !TODO: check the URL version!!
    !-----------------------------
    !
    ! We have to set the ither parameters.
    !-------------------------------------
!    if(present(status)) then
       efh%status = trim(status)
!    endif
    if(present(format)) then
       efh%format = trim(format)
    endif
    if(present(access)) then
       efh%access = trim(access)
    endif
    !
    efh = extio_open_efh(efh)
    !
    ! Have we got any error?
    !-----------------------
    ! Let the return code to set it!
    !-------------------------------
    extio_open_iunp = efh
    iun = efh%unit
    return
    !
  end function extio_open_iunp
!
!H
!H-----------------------------------------------------------------------------
!H  type(extfile) function extio_open_efhp(efh,status,format,access)
!H-----------------------------------------------------------------------------
!H
!
  type(extfile) function extio_open_efhp(efh,status,format,access)
    !
    type(extfile), intent(inout) :: efh
    character(FLCHARS), optional,intent(in) :: status
    character(FLCHARS), optional,intent(in) :: format
    character(FLCHARS), optional,intent(in) :: access
    !
!
!H
!H-----------------------------------------------------------------------------
!H 
!H We are likely to go to open a fully described file.
!H We have the extfile so it not (in principle) for temporary storage and
!H anyway we get the details via parameters. We just have to check if the
!H parameters are there.
!H The dafault case as stated before is an ASCII file.
!H 
!H 
!H-----------------------------------------------------------------------------
!H
!
    !
    !Locals
    !------
    !type(extfile) :: efh
    !
    ! We set our defaults
    !--------------------
    ! Note we leave the other the other parameters unchanged
    ! right now for future use.
    !--------------------------
    efh%unit = -1
    efh%name = ""
    efh%url = ""
    efh%isurl = .false.
    efh%isascii = .true.
    efh%format = 'FORMATTED'
    efh%status = 'UNKNOWN'
    efh%access = 'SEQUENTIAL'
    efh%error = 0
    efh%error_message = "No Errors"
    !
    ! We have to check calling parameters.
    !-------------------------------------
    if(present(status)) then
       efh%status = trim(status)
    endif
    if(present(format)) then
       efh%format = trim(format)
    endif
    if(present(access)) then
       efh%access = trim(access)
    endif
    !
    efh = extio_open_efh(efh)
    !
    ! return the error code of this function
    !-------------------------------
    extio_open_efhp = efh
    return
    !
  end function extio_open_efhp
!
!H
!H-----------------------------------------------------------------------------
!H type(extfile) function extio_open_efh(efh)
!H-----------------------------------------------------------------------------
!H
!
  type(extfile) function extio_open_efh(efh)
!
!H 
!H-----------------------------------------------------------------------------
!H This is the basic function which opens the file. It suppose to get as
!H input parameter an efh handle filled in with all required data.
!H There is a little check about these data (NOT TRUE ANYMORE).
!H It return the efh structure
!H with the coddect unit in it or an error code.
!H On error it set the variable efh%error and the global variable.
!H UPDATE: we do check of the input data ... it means we do this routine
!H the main one which will do all the checks!!!...
!H Actually a check on the keywords like status etc. is done by
!H baseio_procedures and we will not perform them again... except some
!H extension that we will implement in the future ...
!H Like: MEMORY(to force to memory), DUMP(to create a static dump) etc..
!H 
!H-----------------------------------------------------------------------------
!H 
!
    !
    type(extfile), intent(inout) :: efh
    !
    character(FLCHARS) :: status
    character(FLCHARS) :: format
    character(FLCHARS) :: access
    integer(FINT) :: i
    integer(FINT) :: irc
    !
    !
    ! We TRIM all the strings!
    ! and reset some values
    !-------------------------
    efh%name = trim(efh%name)
    efh%url = trim(efh%url)
    efh%status = trim(efh%status)
    efh%format = trim(efh%format)
    efh%access = trim(efh%access)
    efh%error = 0_FINT
    efh%error_message = trim(efh%error_message)
    !if(.not.(efh%isbinary)) efh%isascii = .true.
    eio_error = 0
    !
    ! are we ready?
    !--------------
    if(.not.(eio_isinit)) then
       eio_error = 1
       efh%error = 1
       extio_open_efh = efh
       call message(MESDEBG,"[extio_open]: module not initialized.")
       return
    end if
    !
    ! Do we have correct status format and access strings?
    i = len_trim(efh%status)
    if(windex(baseio_statuskeys,efh%status(1:i)).lt.0) then
       call message(MESERRO,"[extio_open]: Invalid Status.")
       eio_error = 1
       efh%error = 1
       extio_open_efh = efh
       return
    endif
    !
    i = len_trim(efh%format)
    if(windex(baseio_formatkeys,efh%format(1:i)).lt.0) then
       call message(MESERRO,"[extio_open]: Invalid Format.")
       eio_error = 1
       efh%error = 1
       extio_open_efh = efh
       return
    endif
    !
    i = len_trim(efh%access)
    if(windex(baseio_accesskeys,efh%access(1:i)).lt.0) then
       call message(MESERRO,"[extio_open]: Invalid Access.")
       eio_error = 1
       efh%error = 1
       extio_open_efh = efh
       return
    endif
    ! TODO we should check that status format and access are
    !      not in competition with isascii!! ...
    !
    ! Finally we are asked for a file!
    !---------------------------------
    irc = baseio_open(efh%unit,efh%name,efh%status,efh%format,efh%access)
    if(irc.lt.0) then
       call message(MESERRO,'[extio_open]: Cannot open the file.')
       eio_error = 1
       efh%error = 1
       extio_open_efh = efh
       return
    end if
    !
    extio_open_efh = efh
    !
  end function extio_open_efh
  !
  !
  !
!   type(extfile) function extio_openp(status,format,access)
!     !
!     character(FLCHARS), optional :: status
!     character(FLCHARS), optional :: format
!     character(FLCHARS), optional :: access
!     !
!     !
!   end function extio_openp
!   !
!   !
!   !
!   integer(FINT) function extio_open_efh_p(efh,file,status,frmt,acc)
!     !
!     type(extfile), intent(inout) :: efh
!     character(FLCHARS), intent(in), optional :: file
!     character(FLCHARS), intent(in), optional :: status
!     character(FLCHARS), intent(in), optional :: frmt
!     character(FLCHARS), intent(in), optional :: acc
!     !
!     character(LCHARS) :: lfile
!     character(LCHARS) :: lstatus
!     character(LCHARS) :: lfrmt
!     character(LCHARS) :: lacc
!     !
!     integer(FINT) :: liun
!     integer(FINT) :: irc
!     integer(FINT) :: i
!     !
!     ! Init
!     !-----
!     eio_error = 0
!     lfile = ''
!     lstatus = 'OLD'
!     lfrmt = 'FORMATTED'
!     lacc = 'SEQUENTIAL'
!     !
!     ! We get the file name first
!     !---------------------------
!     !---------------------------
!     if(present(file)) then
!        lfile = file
!        lfile = trim(lfile)
!        if(len_trim(efh%name).ne.0) then
!           ! We use the routine with double names
!           ! we trust the given name in the call
!           ! but we produce as well a warning message.
!           !------------------------------------------
!           call message(MESWARN,'[extio_open]: Called with double names.')
!        endif
!        efh%name = lfile
!     else
!        efh%name = trim(efh%name)
!        lfile = efh%name
!     end if
!     !
!     ! Have we got it?
!     if(len_trim(lfile).eq.0) then
!        ! There is not a file name:
!        ! at least there is not any string there.
!        !----------------------------------------
!           eio_error = 1
!           extio_open_efh_p = -1
!           call message(MESERRO,'[extio_open]: No file name given.')
!           return
!     end if
!     !
!     ! We check for the status string.
!     !--------------------------------
!     !--------------------------------
!     if(present(status)) then
!        lstatus = status
!        lstatus = trim(lstatus)
!        i = len_trim(lstatus)
!        if(windex(baseio_statuskeys,lstatus(1:i)).gt.0) then
!           efh%status = lstatus(1:i)
!        else
!           ! Invalid status
!           !---------------
!           eio_error = 1
!           extio_open_efh_p = -1
!           call message(MESERRO,'[extio_open]: Not valid status code.')
!           return
!        end if
!     endif
!     !
!     ! We check for the format string.
!     !--------------------------------
!     !--------------------------------
!     if(present(frmt)) then
!        lfrmt = frmt
!        lfrmt = trim(lfrmt)
!        i = len_trim(lfrmt)
!        if(windex(baseio_formatkeys,lfrmt(1:i)).gt.0) then
!           efh%format = lfrmt(1:i)
!        else
!           ! Invalid format
!           !---------------
!           eio_error = 1
!           extio_open_efh_p = -1
!           call message(MESERRO,'[extio_open]: Not valid format code.')
!           return
!        end if
!     endif
!     !
!     ! We check for the access string.
!     !--------------------------------
!     !--------------------------------
!     if(present(acc)) then
!        lacc = acc
!        lacc = trim(lacc)
!        i = len_trim(lacc)
!        if(windex(baseio_accesskeys,lacc(1:i)).gt.0) then
!           efh%access = lacc(1:i)
!        else
!           ! Invalid access
!           !---------------
!           eio_error = 1
!           extio_open_efh_p = -1
!           call message(MESERRO,'[extio_open]: Not valid access code.')
!           return
!        end if
!     endif
!     !
!     ! We use baseio so we initialize it
!     !----------------------------------
!     if (.not.baseio_isinit) then
!        irc = baseio_init()
!        if (irc.lt.0) then
!           eio_error = 1
!           extio_open_efh_p = -1
!           call message(MESERRO,'[extio_open]: cannot initialize baseio module.')
!           return
!        end if
!     end if
!     !
!     ! We are asked for a file
!     !------------------------
!     irc = baseio_open(liun,lfile,lstatus,lfrmt,lacc)
!     if(irc.lt.0) then
!           call message(MESERRO,'[extio_open]: Cannot open the file.')
!        eio_error = 1
!        extio_open_efh_p = -1
!        return
!     end if
!     !
!     ! We  built the efh type completely
!     !----------------------------------
!     efh%unit = liun
!     if(index('UNFORMATTED',lfrmt).gt.0) then
!        efh%isbinary = .true.
!        efh%isascii = .false.
!     else
!        efh%isbinary = .false.
!        efh%isascii = .true.
!     endif
!     !
!   end function extio_open_efh_p
  !
  !
  !
end module extio
