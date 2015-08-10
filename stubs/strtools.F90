!
!H  
!H----------------------------------------------------------------------
!H----------------------------------------------------------------------
!H  Module strtools
!H----------------------------------------------------------------------
!H----------------------------------------------------------------------
!H $Id: strtools.F90,v 1.4 2002/02/25 05:32:35 fmariot Exp $
!H----------------------------------------------------------------------
!H  
!H This module contain some string manipulation routines.
!H
!H  the variable strtools_error keep track of errors:
!H  It is initialized to 0 and on error is set to >0.
!H  The variable is not reset! unless a function that return
!H  strtools_error goes OK!!
!H
!H Public Procedures:
!H
!H   str_limits    Get string bounds
!H   str_upcase    Transform to UpCase (accept string and vector of strings)
!H   str_lowcase   Transform to LowCase (accept string and vector of strings)
!H   str_tocap     Copy to UpCase
!H   str_tolow     Copy to LowCase
!H   str_strcmp    Compare two strings
!H   str_getfrmt   Take an input and return the fortran format
!H
!H----------------------------------------------------------------------
!H
!
module strtools
  !
  use vartypes
  !
  implicit none
  private
  !
  ! Public Declaration
  !
  public str_limits
  public str_upcase
  public str_lowcase
  public str_tocap
  public str_tolow
  public str_strcmp
  public str_getfrmt
  !
  ! Interfaces declaration
  !
  interface str_upcase
     module procedure str_upcases
     module procedure str_upcasev
  end interface
  !
  interface str_lowcase
     module procedure str_lowcases
     module procedure str_lowcasev
  end interface
  !
  interface str_strcmp
     module procedure str_strcmpss
     module procedure str_strcmpssn
     module procedure str_strcmpssm
     module procedure str_strcmpssnm
     module procedure str_strcmpssnn
  end interface
  !
  interface str_getfrmt
     module procedure str_getfrmti
     module procedure str_getfrmtr
     module procedure str_getfrmtc !for compatibility
  end interface !str_getfrmt
  !
  ! Integer to keep track of errors
  integer(FINT), save :: strtools_error = 0
contains
!
!H
!H----------------------------------------------------------------------
!H    SUBROUTINE STR_LIMITS(STR, FIRST, LAST)
!H----------------------------------------------------------------------
!H
!-----------------------------------------------------------------------
  subroutine str_limits(str, first, last)
!-----------------------------------------------------------------------
!
!
!H
!H----------------------------------------------------------------------
!H
!H    This subroutine finds the "FIRST" and the "LAST" non-blank
!H    characters in the string "STR". The length of the string is not
!H    numerically  limited, but its length is determined when called.
!H
!H
!
    character(*), intent(in) :: str
    integer(FINT),intent(out) :: first
    integer(FINT),intent(out) :: last
    !
    integer(FINT) :: i
    !
    first = 0
    last = len_trim(str)
    do i = 1,last
       if (first.eq.0) then
          if (str(i:i).ne.' ') then
             first = i
             exit
          endif
       endif
    enddo
    !
    !
    return
  end subroutine str_limits
!
!H
!H----------------------------------------------------------------------
!H    SUBROUTINE STR_UPCASES(STRING)
!H----------------------------------------------------------------------
!H
!-----------------------------------------------------------------------
  subroutine str_upcases(string)
!-----------------------------------------------------------------------
!
!H
!H Conversion of string  to uppercase letters.
!H The string is returned in upper case.
!H
!
    character(*), intent(inout) :: string

    integer(FINT) :: i
    integer(FINT) :: iasc

    do i = 1, len_trim(string)
       iasc = ichar(string(i:i))
       if ((iasc.ge.97).and.(iasc.le.122)) then
          string(i:i)= char(iasc - 32)
       endif
    enddo

    return
  end subroutine str_upcases
!
!H
!H----------------------------------------------------------------------
!H    SUBROUTINE UPCASEV(STRING(:))
!H----------------------------------------------------------------------
!H
!-----------------------------------------------------------------------
  subroutine str_upcasev(string)
!-----------------------------------------------------------------------
!
!H
!H Conversion of lowcase letters to up case.
!H It operate on a vector of strings.
!H The string vector is returned in upcase.
!H TODO: now it just call n times upcases
!H
!
    character(*), intent(inout), dimension(:) :: string
    !
    integer(FINT) :: i,idmn
    !
    idmn = size(string)
    !
    do i=1,idmn
       call str_upcases(string(i))
    end do
    !
  end subroutine str_upcasev
!
!H
!H----------------------------------------------------------------------
!H    SUBROUTINE LOWCASES(STRING)
!H----------------------------------------------------------------------
!H
!-----------------------------------------------------------------------
  subroutine str_lowcases(string)
!-----------------------------------------------------------------------
!
!H
!H Conversion of uppercase letters to lower case.
!H The string is returned in lowcase.
!H
!
    character(*), intent(inout) :: string
    !
    integer(FINT) :: i
    integer(FINT) :: iasc
    !
    do i = 1, len_trim(string)
       iasc = ichar(string(i:i))
       if ((iasc.ge.65).and.(iasc.le.90)) then
          string (i:i) = char(iasc + 32)
       endif
    enddo
    !
    !
    return
  end subroutine str_lowcases
!
!H
!H----------------------------------------------------------------------
!H    SUBROUTINE LOWCASEV(STRING(:))
!H----------------------------------------------------------------------
!H
!-----------------------------------------------------------------------
  subroutine str_lowcasev(string)
!-----------------------------------------------------------------------
!
!H
!H Conversion of uppercase letters to lower case.
!H It operate on a vector of strings.
!H The string vector is returned in lowcase.
!H TODO now just call n times lowcases
!H
!
    character(*), intent(inout), dimension(:) :: string
    !
    integer(FINT) :: i,idmn
    !
    idmn = size(string)
    !
    do i=1,idmn
       call str_lowcases(string(i))
    end do
    !
  end subroutine str_lowcasev
!
!H
!H----------------------------------------------------------------------
!H integer(FINT) function str_tocap(instr, line, nchars)
!H----------------------------------------------------------------------
!H
!
  integer(FINT) function str_tocap(instr, line, nchars)
!
!H
!H transform the first NCARS left line (character*(*)) to Upper Case
!H and copy in line.  instr(1:ncars) => lowcase => line(1:ncars)
!H    nchars is optional (behave as UPCASE)
!H    if(nchars.gt.len(line)) the line is returned upcase and
!H    strtools_error is set to nchars.
!H
!
    character(*), intent(out) :: line
    character(*),intent(in) :: instr
    integer(FINT), intent(in), optional :: nchars
    !
    integer(FINT) :: ncars
    integer(FINT) :: i, j
    !
    strtools_error = 0
    str_tocap = 0
    ncars = len(instr)
    if(present(nchars)) then
       if(nchars.lt.ncars) then
          ncars = nchars
       else
          strtools_error = nchars
          str_tocap = -1
          return
       endif
    endif
    if(ncars.gt.len(line)) then
       ncars = len(line)
       strtools_error = ncars
    endif
    do i = 1,ncars
       j = ichar(instr(i:i))
       if (j.ge.97.and.j.le.122) then
          line(i:i) = char(j - 32)
       else
          line(i:i) = instr(i:i)
       endif
    enddo
    !
    str_tocap = ncars
    !
    return
  end function str_tocap
!
!H
!H----------------------------------------------------------------------
!H integer(FINT) function str_tolow(instr, line, nchars)
!H----------------------------------------------------------------------
!H
!H     SubRoutine tolow Modified from old Molden code
!H     .. quite different
!H    nchars is optional (behave as UPCASE)
!H    if(nchars.gt.len(line)) the line is returned upcase and
!H    strtools_error is set to nchars.
!H
!
  integer(FINT) function str_tolow(instr, line, nchars)
!
!H
!H    transform the first NCARS left line (character*(*)) to Lower Case
!H
!
    character(*), intent(in) :: instr
    character(*), intent(out) :: line
    integer(FINT), intent(in), optional :: nchars
    !
    integer(FINT) :: ncars
    integer(FINT) :: i, j
    !
    strtools_error = 0
    str_tolow = 0
    ncars = len(instr)
    if(present(nchars)) then
       if(nchars.lt.ncars) then
          ncars = nchars
       else
          strtools_error = nchars
       endif
    endif
    if(ncars.gt.len(line)) then
       ncars = len(line)
       strtools_error = ncars
       str_tolow = -1
       return
    endif
    do i = 1,ncars
       j = ichar(instr(i:i))
       if (j.ge.65.and.j.le.90) then
          line(i:i) = char(j + 32)
       else
          line(i:i) = instr(i:i)
       endif
    enddo
    !
    str_tolow = ncars
    !
    return
  end function str_tolow
!
!H
!H----------------------------------------------------------------------
!H  INTEGER FUNCTION STR_STRCMPSS(STR1,STR2)
!H----------------------------------------------------------------------
!H
!
  integer(FINT) function str_strcmpss(str1,str2)
!
!H    This routine compare two string and return ilen
!H    if the strings are equal. The Strings are compared
!H    without take into account leading blanks.
!H    Only starting blanks are considered.
!H
!H    return
!H             ilen  if strings are equal
!H             0     if not
!H             -1    on error
!H
!
!
    character(*), intent(in) :: str1
    character(*), intent(in) :: str2
    integer(FINT) :: ilen
    !
    ilen = len_trim(str1)
    if(ilen.eq.len_trim(str2)) then
       if(str1(1:ilen).eq.str2(1:ilen)) then
          str_strcmpss = ilen
       else
          str_strcmpss = 0
       endif
    else
       str_strcmpss = 0
    endif
    return
  end function str_strcmpss
!
!H
!H----------------------------------------------------------------------
!H  INTEGER FUNCTION STR_STRCMPSSN(STR1,STR2,N)
!H----------------------------------------------------------------------
!H
!
  integer(FINT) function str_strcmpssn(str1,str2,n)
!
!H
!H    This routine compare two string.
!H    The Strings are compared up to n chars.
!H
!H    return
!H             >1  if strings are equal
!H             0     if not
!H             -1    on error
!H    strtools_error is set if n.gt.min(len(str1),len(str2))
!H
!
    character(*), intent(in) :: str1
    character(*), intent(in) :: str2
    integer(FINT), intent(in) :: n
    !
    strtools_error = 0
    if(n.gt.min(len_trim(str1),len_trim(str2))) then
       str_strcmpssn = 0
    else
       if(str1(1:n).eq.str2(1:n)) then
          str_strcmpssn = n
       else
          str_strcmpssn = 0
       endif
    endif
    return
  end function str_strcmpssn
!
!
!
!H
!H    ------------------------------------------------------------------
!H    *    INTEGER FUNCTION STR_STRCMPSS(STR1,STR2,MODE)
!H    ------------------------------------------------------------------
!H
!
  integer(FINT) function str_strcmpssm(str1,str2,mode)
!
!H    This routine compare two string and return 1
!H    if the strings are equal. The Strings are compared
!H    to len_trim and
!H    case insensitive if mode='I'or'i' .
!H    case insensitive comarison is made up to FLCHARS (see vartypes)
!H
!H    mode determine: (optional, default = 0)
!H    'I'        Default: Strings are compared case insensitive.
!H
!H    return
!H              1     if strings are equal
!H              0     if not
!H             -1     on error
!H
!H  strtools_error is set on errors.
!H
!
    character(*), intent(in) :: str1
    character(*), intent(in) :: str2
    character(1), intent(in) :: mode
    !
    character(FLCHARS) :: lstr1
    character(FLCHARS) :: lstr2
    integer(FINT) :: il1,il2
    !
    str_strcmpssm = 1
    strtools_error = 0
    !
    il1 = len_trim(str1)
    il2 = len_trim(str2)
    if(il1.eq.il2) then
       if(mode.eq.'I'.or.mode.eq.'i') then
          if(il1.gt.FLCHARS) then
             str_strcmpssm = -1
             strtools_error = 1
          else
             lstr1 = str1
             lstr2 = str2
             call str_lowcase(lstr1)
             call str_lowcase(lstr2)
             if(lstr1(1:il1).ne.lstr2(1:il1)) then
                str_strcmpssm = 0
             endif
          endif
       else
          str_strcmpssm = str_strcmpss(str1,str2)
       endif
    else
       str_strcmpssm = 0
    endif
    !
    return
  end function str_strcmpssm
!
!
!
!H
!H    ------------------------------------------------------------------
!H    *    INTEGER FUNCTION STR_STRCMPSS(STR1,STR2,N,MODE)
!H    ------------------------------------------------------------------
!H
!
  integer(FINT) function str_strcmpssnn(str1,str2,n,mode)
!
!H    This routine compare two string and return ilen
!H    if the strings are equal. The Strings are compared
!H    to the lenght n. if inlen = 0 the two strings
!H    are compared to the whole dimension.
!H
!H    iMode determine: (optional, default = 0)
!H    0        Default - Strings are compared case sensitive.
!H    1        Strings are compared case insensitive.
!H
!H    return
!H              1  if strings are equal
!H              0     if not
!H             -1    on error
!H
!
    character(*), intent(in) :: str1
    character(*), intent(in) :: str2
    integer(FINT), intent(in) :: n
    integer(FINT), intent(in) :: mode
    !
    character(FLCHARS) :: lstr1
    character(FLCHARS) :: lstr2
    integer(FINT) :: il1,il2,il

    str_strcmpssnn = 1

    if(mode.eq.1) then
       il1 = len_trim(str1)
       il2 = len_trim(str2)
       if(il1.eq.il2) then
          if(n.gt.il1) then
             il = il1
          else
             il = n
          endif
       else
          if(n.gt.min(il1,il2)) then
             str_strcmpssnn = 0
             return
          else
             il = n
          endif
       endif
       !
       if(il.gt.FLCHARS) then
          str_strcmpssnn = -1
       else
          lstr1 = str1
          lstr2 = str2
          call str_lowcase(lstr1)
          call str_lowcase(lstr2)
          str_strcmpssnn = 0
          if(lstr1(1:il).eq.lstr2(1:il)) then
             str_strcmpssnn = 1
          endif
       endif
    else
       str_strcmpssnn = str_strcmpssn(str1,str2,n)
    endif
    return
  end function str_strcmpssnn
!
!H
!H    ------------------------------------------------------------------
!H    *    INTEGER FUNCTION STR_STRCMPSS(STR1,STR2)
!H    ------------------------------------------------------------------
!H
!
  integer(FINT) function str_strcmpssnm(str1,str2,n,mode)
!
!H    This routine compare two string and return ilen
!H    if the strings are equal. The Strings are compared
!H    to the lenght inlen. if inlen = 0 the two strings
!H    are compared to the whole dimension.
!H
!H    iMode determine: (optional, default = 0)
!H    0        Default - Strings are compared case sensitive.
!H    1        Strings are compared case insensitive.
!H
!H    return
!H             ilen  if strings are equal
!H             0     if not
!H             -1    on error
!H
!
    character(*), intent(in) :: str1
    character(*), intent(in) :: str2
    integer(FINT), intent(in) :: n
    character(1), intent(in) :: mode
    !
    character(FLCHARS) :: lstr1
    character(FLCHARS) :: lstr2
    integer(FINT) :: il1,il2
    !
    str_strcmpssnm = 1
    !
    il1 = len_trim(str1)
    il2 = len_trim(str2)
    if(n.gt.min(il1,il2)) then
       strtools_error = 1
       str_strcmpssnm = 0
    else
       if(mode.eq.'I'.or.mode.eq.'i') then
          if(il1.gt.FLCHARS.or.il2.gt.FLCHARS) then
             str_strcmpssnm = -1
          else
             lstr1 = str1
             lstr2 = str2
             call str_lowcase(lstr1)
             call str_lowcase(lstr2)
             str_strcmpssnm = str_strcmpssn(lstr1,lstr2,n)
          endif
       else
          str_strcmpssnm = str_strcmpssn(str1,str2,n)
       endif
    endif
    return
  end function str_strcmpssnm
!
!H
!H------------------------------------------------------------------
!H    INTEGER FUNCTION STR_GETFRMTI(STR,IVAL)
!H------------------------------------------------------------------
!H
!
  integer(FINT) function str_getfrmti(str,ival)
!
!H
!H WARNING: This format allows up to I9 Fortran .. which shoud be enought
!H but keep it in mind...........
!H
    !
    character(FLCHARS), intent(out) :: str
    integer(FINT), intent(in) :: ival
    !
    integer(FINT) :: ifield
    integer(FINT) :: ivaltenth
    !
    ifield = 1
    ivaltenth = int( ival / 10_FINT )
    do while (ivaltenth.gt.0)
       ifield = ifield + 1
       ivaltenth = int( ivaltenth / 10_FINT )
    end do
    write(str,'(A1,I1)') 'I',ifield
    !
    str_getfrmti = ifield
    return
    !
  end function str_getfrmti
  !
!
!H
!H------------------------------------------------------------------
!H    INTEGER FUNCTION STR_GETFRMTC(STR,CVAL)
!H------------------------------------------------------------------
!H
!H WARNING: here we return simply 'A'//len(cval)
!
  integer(FINT) function str_getfrmtc(str,cval)
    !
    character(FLCHARS), intent(out) :: str
    character(FLCHARS), intent(in) :: cval
    !
    integer(FINT) :: ifield
    !
    ifield = len(cval)
    str = 'A'
    write(str(2:),*) ifield
    !
    str_getfrmtc = 1
    return
    !
  end function str_getfrmtc
  !
!
!H
!H------------------------------------------------------------------
!H    INTEGER FUNCTION STR_GETFRMTR(STR,RVAL)
!H------------------------------------------------------------------
!H
!H WARNING: we return the E format actually allways E15.8
!H
!
  integer(FINT) function str_getfrmtr(str,rval)
    !
    character(FLCHARS), intent(out) :: str
    real(FREAL), intent(in) :: rval
    !
    str = 'E15.8'
    !
    str_getfrmtr = 1
    return
    !
  end function str_getfrmtr
  !
end module strtools

