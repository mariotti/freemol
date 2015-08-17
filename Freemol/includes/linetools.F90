!
!H
!H  
!H----------------------------------------------------------------------
!H----------------------------------------------------------------------
!H  Module linetools
!H----------------------------------------------------------------------
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H----------------------------------------------------------------------
!H $Id: linetools.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
!H----------------------------------------------------------------------
!H
!
module linetools
!
!H
!H----------------------------------------------------------------------
!H      Some Routine are adapted from FREEREAD package taken from CCL.
!H      CCL: Computational Chemistry List http://www.ccl.org/
!H      A statment from Jan Labanowski:
!H      This software was taken from anon. ftp on rani.chem.yale.edu
!H      It is a set of routines which allow format free input in FORTRAN
!H
!H      Jan Labanowski
!H      jkl@osc.edu
!H----------------------------------------------------------------------
!H
!H This module contain some string manipulations dedicated routines.
!H Specifically we deal with strings readed from an ascii file.
!H
!H----------------------------------------------------------------------
!H
!H  the variable line_error keep track of errors:
!H  It is initialized to 0 and on error is set to >0.
!H
!H----------------------------------------------------------------------
!H We use this convenction:
!H line_xxxx    are public routines
!H l_ln_xxxx    are private routines
!H----------------------------------------------------------------------
!H
!
  !
  use vartypes
  use strtools
  !
  implicit none
  private
  !
  ! Public Declaration
  !
  public line_read
  public line_gettype
  public line_getform
  public line_allrea
  public line_getline
  public line_getval
  public line_get_oneval
  public line_get_varval
  !
  public windex
  !
  ! Interfaces
  !
  interface line_get_oneval
     module procedure line_get_oneval_i
     module procedure line_get_oneval_r
  end interface !line_get_oneval
  !
  interface line_get_varval
     module procedure line_get_varval_i
     module procedure line_get_varval_r
  end interface !line_get_varval
  !
  !
  !
  ! Error Interface
  !----------------
  ! We define two error interfaces. Each one return the latest error
  ! code or string value. The String error function can be used to inquire
  ! about a given error code.
  !-----------------------------------------------------------------------
  interface line_geterror
     module procedure line_geterror_i
     module procedure line_geterror_c
  end interface !line_geterror
  !
  !
  ! Integer to keep track of errors
  integer(FINT), public, save :: l_ln_error = 0
  ! Logical to define the initialization state
  !  logical, private, save :: l_ln_isinit = .false.
  ! We set that (up to date) this module is always
  ! initialized
  logical, private, save :: l_ln_isinit = .true.
  !
contains
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
!H integer(FINT) function line_init()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function line_init()
    !
    line_init = 0
    l_ln_isinit = .true.
    !
  end function line_init
!
!H
!H-----------------------------------------------------------------------------
!H logical function line_isinit()
!H-----------------------------------------------------------------------------
!H 
!
  logical function line_isinit()
    line_isinit = l_ln_isinit
  end function line_isinit
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function line_geterror_i()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function line_geterror_i()
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
    line_geterror_i = l_ln_error
    !
  end function line_geterror_i
!
!H-----------------------------------------------------------------------------
!H
!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H character(FLCHARS) function line_geterror_c()
!H-----------------------------------------------------------------------------
!H
!
  character(FLCHARS) function line_geterror_c(code)
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
    ! We return the error code
    !-------------------------
    ! Note that we have not yet implement it
    ! so we report as example case 2_FINT as unwanted code
    ! but it has to disappear! And Default is a not reported
    ! error code.
    !-------------------------------------------------------
    select case(code)
    case(0_FINT)
       line_geterror_c = "No Errors."
    case(1_FINT)
       line_geterror_c = "Not documented error or generic error.."
    case(2_FINT)
       line_geterror_c = "Double Error: Why did you get this message?."
    case default
       line_geterror_c = "Not documented error or generic error.."
    end select
    !
    return
    !
  end function line_geterror_c

!
!H
!H-----------------------------------------------------------------------------
!H END OF: SOME RETURN VALUE FUNCTIONS
!H-----------------------------------------------------------------------------
!H 
!
!
!H
!H----------------------------------------------------------------------
!H   SUBROUTINE LINE_READ(LINE,IWRD,RWRD,SWRD,NIWRD,NRWRD,NSWRD,FORM)
!H----------------------------------------------------------------------
!H   
!H   This routine read a line and return integers reals and string
!H   in the given vectors. A string FORM gives back the format
!H   as in the freeread routine.
!H   
!H
!H         This subroutine does take an input line, separates it into
!H    integers, reals, and characters as guided by blanks, and load it
!H    into the corresponding arrays.
!H
!H    INPUT :  LINE
!H    OUTPUT : NITGR (the number of integers found)
!H             NREA (the number of reals found)
!H             NCHA (the number of character strings found)
!H             INVA (a *-Element array containing all the integers)
!H             REA (a *-Element array containing all the reals)
!H             CHA (a *-Element array containing the strings)
!H
!H
!H----------------------------------------------------------------------
!H   
!
  subroutine line_read(line,iwrd,rwrd,swrd,niwrd,nrwrd,nswrd,form)
    implicit none
    !
    character(*), intent(inout) :: line
    character(*), intent(out) :: form
    integer(FINT), intent(out) :: niwrd
    integer(FINT), intent(out) :: nrwrd
    integer(FINT), intent(out) :: nswrd
    integer(FINT), intent(out), dimension(:) :: iwrd
    real(FREAL), intent(out), dimension(:) :: rwrd
    character(*), intent(out), dimension(:) :: swrd
    !
    logical :: notdone
    integer(FINT) :: fdmn
    integer(FINT) :: idmn
    integer(FINT) :: rdmn
    integer(FINT) :: sdmn
    integer(FINT) :: wdmn
    integer(FINT) :: i
    integer(FINT) :: wcount
    character(1) :: blk
    character(FLCHARS) :: word
    !
    blk = ' '
    notdone = .true.
    i=1
    wcount = 0
    niwrd = 0
    nrwrd = 0
    nswrd = 0
    fdmn = len(form)
    idmn = size(iwrd)
    rdmn = size(rwrd)
    sdmn = size(swrd)
    wdmn = len(swrd(1))
    !
    do while (notdone)
       line = adjustl(line(i:))
       i = index(line,blk)
       wcount = wcount + 1
       iswordif: if(i.gt.1) then
          word = line(1:i-1)
          if(wcount.gt.fdmn) then
             print*,'MS manage this error!'
             stop 3
          endif
          form(wcount:wcount) = line_gettype(word)
          select case(form(wcount:wcount))
          case('I')
             niwrd = niwrd + 1
             if (niwrd.gt.idmn) then
                print*,'MS too much!!'
                stop 3
             else
                read(word,*) iwrd(niwrd)
             endif
          case('F','D')        !case('D')
             nrwrd = nrwrd + 1
             if (nrwrd.gt.rdmn) then
                print*,'MS too much!!'
                stop 3
             else
                read(word,*) rwrd(nrwrd)
             endif
          case('A')
             nswrd = nswrd + 1
             if(i-1.gt.wdmn) then
                print*,'MS: word maybe truncated'
             endif
             if (nswrd.gt.sdmn) then
                print*,'MS too much!!'
                stop 3
             else
                swrd(nswrd) = word
             endif
          case default
          end select
       else
          notdone = .false.
       endif iswordif
    end do
    !
    return
    !
  end subroutine line_read
!
!H
!H----------------------------------------------------------------------
!H
!
!
!H
!H----------------------------------------------------------------------
!H   CHARACTER(1) FUNCTION LINE_GETTYPE(WORD)
!H----------------------------------------------------------------------
!H   
!H   This function return the type of a word. type is 'I','F','D' or 'A'.
!H   
!H----------------------------------------------------------------------
!H   
!
  character(1) function line_gettype(word)

    character(*), intent(in) :: word

    logical :: notdone
    logical :: isdecimal
    logical :: isexpo
    logical :: issign
    integer(FINT) :: i,iwdmn
    
    line_gettype = 'I'
    notdone = .true.
    isdecimal = .false.
    isexpo = .false.
    issign = .false. !sign is considered only in exponentials
    iwdmn = len_trim(word)
    i = 1
    
    if(ichar(word(1:1)).ge.ichar('0').and.ichar(word(1:1)).le.ichar('9')) then
       line_gettype = 'I'
    else if (word(1:1).eq.'-'.or.word(1:1).eq.'+') then
       line_gettype = 'I'
       if(iwdmn.eq.1) then
          line_gettype = 'A'
          notdone = .false.
       endif
    else if (word(1:1).eq.'.') then
       line_gettype = 'F'
    else
       line_gettype = 'A'
       notdone = .false.
    endif

    do while (notdone)
       i = i + 1
       noendwrdif: if(i.gt.iwdmn) then
          notdone = .false.
       else
          if(word(i:i).eq.'-'.or.word(i:i).eq.'+') then
             if(.not.isexpo.or.(isexpo.and.issign)) then
                line_gettype = 'A'
                notdone = .false.
             else
                issign = .true.
             endif
          else if(word(i:i).eq.'.') then
             if(line_gettype.eq.'I') then
                line_gettype = 'F'
             else
                line_gettype = 'A'
                notdone = .false.
             endif
          else if(word(i:i).eq.'E'.or.word(i:i).eq.'e') then
             if(line_gettype.eq.'F') then
                if(isexpo) then
                   line_gettype = 'A'
                   notdone = .false.
                else if ((i.eq.2).and.(word(1:1).eq.'-'.or.word(1:1).eq.'+')) then
                   line_gettype = 'A'
                   notdone = .false.
                else
                   line_gettype = 'F'
                   isexpo = .true.
                endif
             else
                line_gettype = 'A'
                notdone = .false.
             endif
          else if(word(i:i).eq.'D'.or.word(i:i).eq.'d') then
             if(line_gettype.eq.'F') then
                if(isexpo) then
                   line_gettype = 'A'
                   notdone = .false.
                else if ((i.eq.2).and.(word(1:1).eq.'-'.or.word(1:1).eq.'+')) then
                   line_gettype = 'A'
                   notdone = .false.
                else
                   line_gettype = 'D'
                   isexpo = .true.
                endif
             else
                line_gettype = 'A'
                notdone = .false.
             endif
          else if((ichar(word(i:i)).lt.ichar('0'))&
               &.or.(ichar(word(i:i)).gt.ichar('9'))) then
             line_gettype = 'A'
             notdone = .false.
          endif
       endif noendwrdif
    end do
    !
    return
  end function line_gettype
!
!H
!H----------------------------------------------------------------------
!H
!
!
!H
!H----------------------------------------------------------------------
!H   SUBROUTINE LINE_GETFORM(WORD,TYP,EFRMT[,IDMN,RKIND])
!H----------------------------------------------------------------------
!H   
!H   This subroutine return the format to read the word of type typ.
!H   The format is complete and optionally it returns the string dimension
!H   of the format string and in the future the needed kind.
!H   
!H----------------------------------------------------------------------
!H   
!
  subroutine line_getform(word,typ,efrmt,idmn,rkind)
    !
    !
    character(*), intent(in) :: word
    character(1), intent(in) :: typ
    character(*), intent(out) :: efrmt
    integer(FINT), intent(out), optional :: idmn
    integer(FINT), intent(out), optional :: rkind
    !
    character(2) :: lstr
    character(2) :: sdec
    integer(FINT) :: dec
    integer(FINT) :: lwrd
    integer(FINT) :: istr
    integer(FINT) :: idec
    !
    lwrd = len_trim(word)
    if(present(rkind)) then
       rkind = -1
       print *,'MS: function not yet supported!'
    endif
    !
    if(lwrd.gt.99) then
       print *,'MS: internal limit: too big number?'
    else if(lwrd.eq.0) then
       print*,'MS: Is it a word?'
    endif
    !
    write(lstr,'(I2)') lwrd
    lstr = adjustl(lstr)
    istr = int(log10(real(lwrd))) + 1
    !
    idec = 0
    efrmt = typ
    efrmt = efrmt(1:1)//lstr
    istr = len_trim(efrmt)
    !
    dec = index(word,'.')
    if(typ.eq.'F'.and.dec.gt.0) then
       dec = lwrd - dec
       idec = max(index(word,'E'),index(word,'e'))
       if(idec.gt.0) dec = dec - lwrd + idec - 1
       if(dec.gt.0) then
          idec = int(log10(real(dec)))+1
          write(sdec,'(I2)') dec
          sdec=adjustl(sdec)
          efrmt = efrmt(1:istr)//'.'//sdec(1:idec)
          istr = len_trim(efrmt)
       endif
    else if(typ.eq.'D'.and.dec.gt.0) then
       dec = lwrd - index(word,'.')
       idec = max(index(word,'D'),index(word,'d'))
       if(idec.gt.0) dec = dec - lwrd + idec - 1
       if(dec.gt.0) then
          idec = int(log10(real(dec)))+1
          write(sdec,'(I2)') dec
          sdec=adjustl(sdec)
          efrmt = efrmt(1:istr)//'.'//sdec(1:idec)
          istr = len_trim(efrmt)
       endif
    end if
    !
    if(present(idmn)) idmn = istr
    !
    return
    !
  end subroutine line_getform
!
!H
!H----------------------------------------------------------------------
!H    SUBROUTINE LINE_ALLREA(FORM, INVA, REA, NINT, NREA)
!H----------------------------------------------------------------------
!H
!-----------------------------------------------------------------------
  subroutine line_allrea(form, inva, rea, nitgr, nrea)
!-----------------------------------------------------------------------
!
!H
!H    This subroutine takes the INVA and REA arrays obtained from
!H    FREEREAD, converts all integers to real numbers and stores them in
!H    the REA array while restoring the same order in which they were in
!H    the original line.
!H
!H    INPUT :  FORM, a string containing the types of the different
!H             substrings found in the LINE by FREEREAD (I, F, or A).
!H             INVA, the array containing all the integers
!H             REA, the array containing the real numbers
!H             NITGR and NREA, the number of integers and reals,
!H             respectively.
!H
!H    OUTPUT : NREA and REA
!H
!W    USES IMPLICIT DIMENSIONS
!
    character(*) :: form
    integer(FINT) :: nitgr
    integer(FINT) :: nrea
    real(FREAL), dimension(:) :: rea
    integer(FINT), dimension(:) :: inva
    !
    integer(FINT) :: iitgr
    integer(FINT) :: irea
    real(FREAL) :: oldrea (80)
    !
    integer(FINT) :: i
    integer(FINT) :: is
    integer(FINT) :: ie
    !
    if ( (nitgr.eq.0) .and. (nrea.eq.0) ) return
    do i = 1, nrea
       oldrea (i) = rea (i)
    enddo
    iitgr = 0
    nrea = 0
    irea = 0
    call str_limits (form, is, ie)
    do i = is, ie
       if (form (i:i) .ne.'A') then
          nrea = nrea + 1
          if (form (i:i) .eq.'F') then
             irea = irea + 1
             rea (nrea) = oldrea (irea)
          else
             iitgr = iitgr + 1
             rea (nrea) = inva (iitgr)
          endif
       endif
    enddo
    !
    !
    !
    return
  end subroutine line_allrea
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function line_getline (iun, line, mode)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function line_getline (iun, line, mode)
!
!
!H
!H-----------------------------------------------------------------------------
!H    Read a LINE from unit IUN. If the first line is blank read it again.
!H    INPUT     IUN    file unit
!H              IOMODE    set behaviour (optional, default = 0)
!H                     0 Default: Read the line,
!H                                strip initial blanks,
!H                                convert TABS to 1 (one) blank char.
!H                     1 MoreCHR: As 0 but
!H                                strip =.
!H                     3 Escape commented lines: this mode
!H                                handle commented lines ( first
!H                                character is # as blank lines.)
!H                                It behave as mode 0 for the other cases.
!H
!H                     9          Just read the line.
!H
!H
!H    OUTPUT    LINE   the readed and processed line.
!H
!H              LINE_GETLINE IOSTAT READ VALUE.
!H                     0  no errors.
!H                     <0 on ERROR.
!H-----------------------------------------------------------------------------
!H
!
    character(*) :: line
    integer(FINT) :: iun
    integer(FINT), optional :: mode
    !
    logical :: notdone
    integer(FINT) :: i
    integer(FINT) :: ios
    integer(FINT) :: iomode
    !
    notdone = .true.
    l_ln_error = 0
    line_getline = 0
    iomode = 0
    if(present(mode)) then
       iomode = mode
    endif

    do while(notdone)
        ! We default to read only one line
        ! the only exception is a blank line:
        ! we read it again.
       notdone = .false.
       read (iun,'(a)',iostat=ios) line
       if (ios.eq.0) then
          if (.not.iomode.eq.9) then
             !--------Now every other iomode are defaults
             ios = len(line)
             !
             if (iomode.eq.1) then
                do i = 1,ios
                   if (ichar(line(i:i)).eq.61.or.ichar(line(i:i)).eq.9) then
                      line(i:i) = ' '
                   endif
                enddo
             else ! this else include all modes but 1
                do i = 1,ios
                   if (ichar(line(i:i)).eq.9) line(i:i) = ' '
                enddo
             endif
             line = adjustl(line)
             if(line(1:1).eq.' ') then
                !Blank Line
                notdone = .true.
             endif
             if(iomode.eq.3.and.(line(1:1).eq.'#')) then
                ! Comment Line and mode is 3
                notdone = .true.
             end if
          endif
       else
          l_ln_error = 1
          line_getline = -1
       endif
    enddo
    !
    !
    !
    return
  end function line_getline


!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function line_get_varval_r(line,var,val)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function line_get_varval_r(line,var,val)
!
!H Get if present the first valued of couple i.e. line as "val=3.14159" return
!H the value. Lines as "val1 val2=3.14159" return val2 value (i.e. the first)!
!H var is destroied anyway!
!H
!

    character(FLCHARS), intent(in) :: line
    character(FLCHARS), intent(out) :: var
    real(FREAL), intent(out) :: val

    character(FLCHARS) :: sdummy
    integer(FINT) :: i,j

    line_get_varval_r = 0 ! No data in the line

    sdummy = line
    i = index(sdummy,'=')
    if(i.gt.0) then
       do while(.true.)
          j = index(sdummy,' ')
          if(j.lt.1) then
             exit
          endif
          if(j.lt.i) then
             j = j + 1
             sdummy = sdummy(j:)
          else
             line_get_varval_r = line_get_oneval_r(line,var,val)
             exit
          endif
          i = index(sdummy,'=')
       enddo
    endif

  end function line_get_varval_r




!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function line_get_varval_i(line,var,val)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function line_get_varval_i(line,var,val)
!
!H Get if present the first valued of couple i.e. line as "val=3.14159" return
!H the value. Lines as "val1 val2=3.14159" return val2 value (i.e. the first)!
!H var is destroied anyway!
!H
!

    character(FLCHARS), intent(in) :: line
    character(FLCHARS), intent(out) :: var
    integer(FINT), intent(out) :: val

    character(FLCHARS) :: sdummy
    integer(FINT) :: i,j

    line_get_varval_i = 0 ! No data in the line

    sdummy = line
    i = index(sdummy,'=')
    if(i.gt.0) then
       do while(.true.)
          j = index(sdummy,' ')
          if(j.lt.i) then
             j = j + 1
             sdummy = sdummy(j:)
          else
             line_get_varval_i = line_get_oneval_i(line,var,val)
             exit
          endif
          i = index(sdummy,'=')
       enddo
    endif

  end function line_get_varval_i




!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function line_get_oneval_r(line,var,val)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function line_get_oneval_r(line,var,val)
!
!H
!H Get if present the first valued of couple i.e. line as "val=3.14159" return
!H the value. Lines as "val1 val2=3.14159" return nothing (i.e. only first)!
!H var is destroied anyway!
!H
!
    character(*), intent(in) :: line
    character(*), intent(out) :: var
    real(FREAL), intent(out) :: val

    integer(FINT) :: i,j
    character(FLCHARS) :: sdummy
    character(FLCHARS) :: efrmt
    character(1) :: styp

    line_get_oneval_r = 0 ! No data in the line

    sdummy = adjustl(line)
    call str_lowcase(sdummy)
    if(.not.(len_trim(sdummy).lt.1)) then
       i = index(sdummy,'=')
       if(i.gt.0) then
          var = sdummy(1:i)
          j = index(sdummy,' ')
          if(j.gt.i) then
             i = i -1
             var = sdummy(1:i)
             i = i + 2
             sdummy = adjustl(sdummy(i:j))
             styp = line_gettype(sdummy)
             if(styp.eq.'F'.or.styp.eq.'D') then
                call line_getform(sdummy,styp,efrmt)
                read(sdummy,'('//efrmt//')') val
                line_get_oneval_r = 1
             endif
          else
             sdummy = adjustl(sdummy(j:))
             if(sdummy(1:1).eq.'=') then
                sdummy = adjustl(sdummy(2:))
                styp = line_gettype(sdummy)
                if(styp.eq.'F'.or.styp.eq.'D') then
                   call line_getform(sdummy,styp,efrmt)
                   read(sdummy,'('//efrmt//')') val
                   line_get_oneval_r = 1
                endif
             endif
          endif
       endif
    endif

  end function line_get_oneval_r



!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function line_get_oneval_i(line,var,val)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function line_get_oneval_i(line,var,val)
!
!H Get if present the first valued of couple i.e. line as "val=3.14159" return
!H the value. Lines as "val1 val2=3.14159" return nothing (i.e. only first)!
!H var is destroied anyway!
!H
!
    character(*), intent(in) :: line
    character(*), intent(out) :: var
    integer(FINT), intent(out) :: val

    integer(FINT) :: i,j
    character(FLCHARS) :: sdummy
    character(FLCHARS) :: efrmt
    character(1) :: styp

    line_get_oneval_i = 0 ! No data in the line

    sdummy = adjustl(line)
    call str_lowcase(sdummy)
    if(.not.(len_trim(sdummy).lt.1)) then
       i = index(sdummy,'=')
       if(i.gt.0) then
          j = index(sdummy,' ')
          if(j.gt.i) then
             i = i -1
             var = sdummy(1:i)
             i = i + 2
             sdummy = sdummy(i:j)
             styp = line_gettype(sdummy)
             if(styp.eq.'I') then
                call line_getform(sdummy,styp,efrmt)
                read(sdummy,'('//efrmt//')') val
                line_get_oneval_i = 1
             endif
          else
             sdummy = adjustl(sdummy(j:))
             if(sdummy(1:1).eq.'=') then
                sdummy = adjustl(sdummy(2:))
                styp = line_gettype(sdummy)
                if(styp.eq.'I') then
                   call line_getform(sdummy,styp,efrmt)
                   read(sdummy,'('//efrmt//')') val
                   line_get_oneval_i = 1
                endif
             endif
          endif
       endif
    endif
    !
  end function line_get_oneval_i
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function line_getval(line,var,val)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function line_getval(line,var,val)
!
!H Read line to search for the couple var=val were var match the input
!H string var. Returns val as string. It gets the first instance of var.
!H VAR is checked case insensitive.
!H
!
    character(*), intent(in) :: line
    character(*), intent(in) :: var
    character(*), intent(out) :: val
    !
    character(FLCHARS) :: lline
    character(FLCHARS) :: lvar
    integer(FINT) :: i,ipos
    character(FLCHARS) :: sdummy
    !character(FLCHARS) :: efrmt
    !character(1) :: styp
    !
    line_getval = 0
    val=''
    if(len(line).gt.FLCHARS) then
       line_getval = -1
       return
    end if
    lline = adjustl(line)
    call str_lowcase(lline)
    lvar=var
    call str_lowcase(lvar)
    i=index(lline,lvar(1:len_trim(lvar)))
    if(i.gt.0) then
       ipos=i+len_trim(lvar)
       sdummy=lline(ipos:)
       sdummy=adjustl(sdummy)
       if(sdummy(1:1).eq.'=') then
          if(len_trim(sdummy).eq.1) then
             return
          end if
          sdummy=sdummy(2:)
          ipos=index(sdummy,' ')
          if(ipos.gt.0) then
             sdummy=sdummy(1:ipos)
          end if
          val=sdummy
          !That's the normal exit point
          return
       end if
    end if
    ! The true Exit point is inside the ifs stuff!
    line_getval = -1
    !
  end function line_getval
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) recursive function windex(line,key,back)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) recursive function windex(line,key,back) result(windexr)
!
!H
!H-----------------------------------------------------------------------------
!H This routine is en equivalent of the routine index but
!H it return the index only if the key match an exact word in line.
!H In this preliminar version a word is separated by blanks.
!H TODO: set up something to define the word separators.
!H-----------------------------------------------------------------------------
!H
!
    !
    character(*), intent(in) :: line
    character(*), intent(in) :: key
    logical, optional :: back
    !
    character(FLCHARS) :: lkey
    integer(FINT) :: lilen
    integer(FINT) :: i
    logical :: lback
    ! Now we use wsep as parameter
    ! but in the future can be a variable to assign the value
    ! and is we need more separators we will need to manage
    ! with the code.
    character(1), parameter :: wsep = ' '
    !
    l_ln_error = 0
    windexr = 0
    !
    lback = .false.
    if(present(back)) then
       lback = back
    end if
    if(len_trim(key).gt.FLCHARS) then
       ! WARNING: actually we have to take
       !          into account the blanks at the top
       !          But right now we don't manage them.
       l_ln_error = 1
       windexr = -1
       return
    end if
    lkey = trim(key)
    if(len_trim(lkey).gt.0) then
       lilen = len_trim(lkey)
    else
       ! WARNING: we produce an error if the
       !          key is not a word!
       ! This is different from the default of the
       ! function index. But in this case
       ! we search for words.. so we need it at least.
       ! to bypass this and get the same behaviour
       ! of the index function check the return
       ! code in this way to catch not in line:
       ! if(windex(....).lt.1) print *, ' error '
       !
       windexr = -1
       l_ln_error = 1
       return
    end if
    !
    if(lback) then
       i = index(line,lkey(1:lilen))
    else
       i = index(line,lkey(1:lilen),lback)
    end if
    ! We check the previous character
    !--------------------------------
    if(i.gt.1) then
       if(line(i-1:i-1).ne.' ') then
          windexr = -1
       end if
    endif
    !
    ! We check the following character
    !---------------------------------
    if((i+lilen-1).lt.len(line)) then
       if(line(i+lilen:i+lilen).ne.' ') then
          windexr = -1
       end if
    endif
    !
    ! We check windexr to decide to search over
    !-----------------------------------------
    if (windexr.eq.-1) then
       if ((i+lilen).le.len(line)) then
          i = windex(line(i+lilen:len(line)),lkey,lback)
       end if
    end if
    !
    windexr = i
    return
    !
  end function windex
  !
  !
  !
end module linetools
