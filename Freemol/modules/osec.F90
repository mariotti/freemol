!
!H
!H-----------------------------------------------------------------------------
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H MODULE osec  Freemol by F.Mariotti: (c) F.Mariotti
!H-----------------------------------------------------------------------------
!H $Id: osec.F90,v 1.1.1.1 2009/01/12 16:56:17 mariotti Exp $
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
module osec
  !
  use vartypes
  use strtools
  use linetools
  !
  implicit none
  private
  !
  ! PUBLIC DECLARATIONS:
  !--------------------
  
  ! Public Procedures
  !------------------
  public osec_init
  public osec_isinit
  public osec_geterror
  public osec_set
  !public osec_
  !
  !
  ! Error Interface
  !----------------
  ! We define two error interfaces. Each one return the latest error
  ! code or string value. The String error function can be used to inqure
  ! about a given error code.
  !----------------------------------------------------------------------
  interface osec_geterror
     module procedure osec_geterror_i
     module procedure osec_geterror_c
  end interface !osec_geterror

  ! LOCAL DECLARATIONS:
  !--------------------
  !
  logical, save :: losec_isinit = .false.
  !
  ! error variables
  !----------------
  integer(FINT), save, private :: losec_error
  character(FLCHARS), save, private :: losec_error_message
  
contains
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H Function osec_init()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function osec_init()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the osec module.
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    osec_init = 1
    !
  end function osec_init
!
!H
!H-----------------------------------------------------------------------------
!H logical function osec_isinit()
!H-----------------------------------------------------------------------------
!H
!
  logical function osec_isinit()
!
!H
!H-----------------------------------------------------------------------------
!H We merelyreturn the initialization flag.
!H-----------------------------------------------------------------------------
!H
!
    !
    osec_isinit = losec_isinit
    !
  end function osec_isinit
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function osec_geterror_i()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function osec_geterror_i()
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
    osec_geterror_i = losec_error
    !
  end function osec_geterror_i
!
!H-----------------------------------------------------------------------------
!H
!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H character(FLCHARS) function osec_geterror_c()
!H-----------------------------------------------------------------------------
!H
!
  character(FLCHARS) function osec_geterror_c(code)
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
       osec_geterror_c = "No Errors."
    case(1_FINT)
       osec_geterror_c = "Not documented error or generic error.."
    case(2_FINT)
       osec_geterror_c = "Double Error: Why did you get this message?."
    case default
       osec_geterror_c = "Not documented error or generic error.."
    end select
    !
    return
    !
  end function osec_geterror_c
!H
!H
!H------------------------------------------------------------------------------
!H integer(FINT) function osec_set(iun,ssec,scpars,ilev)
!H------------------------------------------------------------------------------
!H
!
  integer(FINT) function osec_set(iun,ssec,scpars,ilev)
!
!H------------------------------------------------------------------------------
!H    Freemol Format Set Section
!H
!H    Search a section and set current file record to the first
!H    section line. The file is searched from the current position.
!H
!H
!H    input
!H          IUN          integer file unit
!H
!H          SSEC         string section title with no sqare braket
!H
!H    outout
!H          mr_setsec   integer recturn code:
!H                      -n Undefined error (Not implemented)
!H                      -3 File is not opened for FreemolFormat
!H                      -2 Received a void sec name
!H                      -1 EOF Error
!H                       0 Not found
!H                       1 found
!H                       n To be defined (Not implemented)
!H
!H          scpars       String with inlined parameters of the section
!H
!H    Variables
!H
!H          line         string*(LCHARS) readed line from mld file
!H
!
!
    !
    !
    integer(FINT), intent(in) :: iun
    character(*),intent(in) :: ssec
    character(*),intent(out) :: scpars
    integer(FINT), intent(in), optional :: ilev ! optional search level
    !
    integer(FINT) :: isec
    integer(FINT) :: irc
    integer(FINT) :: lilev
    character(FLCHARS) :: lsec
    character(FLCHARS) :: line
    logical :: done
    logical :: llev
    !
    llev = .false.
    lilev = 0
    if(present(ilev)) then
       lilev = ilev
       llev = .true.
    end if
    !
    osec_set = 0
    losec_error = 0
    losec_error_message = " "
    scpars = " "
    !
    if (len_trim(ssec).gt.FLCHARS) then
       osec_set = -1
       losec_error = 1
       losec_error_message = "[osec_set] Too long section name."
       return
    end if
    lsec = trim(ssec)
    call str_lowcase(lsec)
    if(scan(lsec," ").eq.0) then
       osec_set = -1
       losec_error = 1
       losec_error_message = "[osec_set] Not a valid section."
       return
    end if
    isec = len_trim(lsec)
    !
    done = .false.
    do while (.not.done)
       irc = line_getline(iun,line)
       if (irc.lt.0) then
          osec_set = -1
          losec_error = 1
          losec_error_message = "[osec_set] Error reading line."
          return
       end if
       !
       !--------Search section line
       if (line(1:1).eq.'[') then
          irc = index(line,'[',.true.)
          if (llev) then
             if (lilev.ne.irc) cycle
          endif
          line = line(irc+1:)
          if(line(1:len_trim(lsec)).eq.lsec(1:len_trim(lsec))) then
             irc = index(line,']',.true.)
             if (irc.gt.1.and.irc.le.len(line)) then
                scpars=trim(line(irc+1:))
             end if
             osec_set = 1
             return
          end if
       end if
       !
    end do
    !
    !
  end function osec_set
!
!H
!H-----------------------------------------------------------------------------
!H
!

  !
  !
end module osec






