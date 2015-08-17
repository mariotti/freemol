!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H  MODULE pcmdline  Freemol by F.Mariotti: (c) F.Mariotti
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H $Id: pcmdline.F90,v 1.4 2002/03/06 12:15:27 fmariot Exp $
!H-----------------------------------------------------------------------------
!H  
!H  KeyWords: standard std f90 iargc iargv getarg
!H  
!H  This module define some routines for cmdline mamagement.
!H  
!H  Usage:
!H  Use pcmdline
!H      Define subroutines:
!H      md_iargc(integer(KINT)::N)
!H          Return the Number of cmdline arguments.
!H      md_getarg(integer(KINT)::N,character(*)::Str)
!H          Return the N-th command line argument.
!H  
!H  NOT COMPLETED: TODO ERROR HANDLING
!H  
!H  
!H  
!H  
!H  ---------------------------------------------------------------------------
!H  
!
module pcmdline

  use vartypes

  implicit none
  private

  integer(FINT), public, save :: pcmd_error = 0

  public pcmd_iargc
  public pcmd_getarg
  public pcmd_getio
  public pcmd_checkarg

  interface pcmd_checkarg
     module procedure pcmd_checkargn
     module procedure pcmd_checkargp
  end interface !pcmd_checkarg

contains
  !
!H 
!H-----------------------------------------------------------------------------
!H   subroutine pcmd_getio(N,cmdstr,fnin,fnout)
!H-----------------------------------------------------------------------------
!H 
!H   
!H   This routine search for -i and -o command line options
!H   and return the argument if any.
!H   
!H   
!H   
!H   
  subroutine pcmd_getio(N,cmdstr,fnin,fnout)
    !
    integer(FINT), intent(in) :: N
    character(FLCHARS), dimension(:), intent(in) :: cmdstr
    character(FLCHARS), intent(inout) :: fnin,fnout
    !
    integer(FINT) :: i
    !
    pcmd_error = 0
    !
    ! Search for -i -o options
    !
    do i=1,N
       if(len_trim(cmdstr(i)).eq.2.and.cmdstr(i)(1:2).eq.'-i') then
          if(i.lt.N) then
             select case (cmdstr(i+1)(1:1))
             case ('-') ! check if std input required
                if(len_trim(cmdstr(i+1)).eq.1) then
                   fnin = '-'
                else
                   pcmd_error = -1
                   write(0,'("Error: change this code to use messages")')
                   write(0,'("Error: invalid file name given with -i option")')
                   stop 1
                endif
             case (' ','$') ! it not a file name
                pcmd_error = -1
                write(0,'("Error: change this code to use messages")')
                write(0,'("Error: invalid file name given with -i option")')
                stop 1
             case default
                fnin = cmdstr(i+1)
             end select
          else
             pcmd_error = -1
             write(0,'("Error: change this code to use messages")')
             write(0,'("Error: no input file name given with -i option")')
             stop 1
          endif
       elseif(len_trim(cmdstr(i)).eq.2.and.cmdstr(i)(1:2).eq.'-o') then
          if(i.lt.N) then
             select case (cmdstr(i+1)(1:1))
                case ('-') ! check if std input required
                   if(len_trim(cmdstr(i+1)).eq.1) then
                      fnout = '-'
                   else
                      pcmd_error = -1
                      write(0,'("Error: change this code to use messages")')
                      write(0,'("Error: invalid file name given with -o option")')
                      stop 1
                   endif
                case (' ','$') ! it not a file name
                   pcmd_error = -1
                   write(0,'("Error: change this code to use messages")')
                   write(0,'("Error: invalid file name given with -o option")')
                   stop 1
                case default
                   fnout = cmdstr(i+1)
                end select
          else
             pcmd_error = -1
             write(0,'("Error: change this code to use messages")')
             write(0,'("Error: no output file name given with -o option")')
             stop 1
          endif
       endif
    enddo
    !
  end subroutine pcmd_getio
  !
  !
  !
!H 
!H-----------------------------------------------------------------------------
!H   subroutine pcmd_checkargn(N,cmdstr,arg)
!H-----------------------------------------------------------------------------
!H 
!H   
!H   This routine search for the given argument arg
!H   and return the the number of presence on the command line.
!H   Does not check for parameters to arg.
!H   
!H   
!H   
  integer(FINT) function pcmd_checkargn(N,cmdstr,arg)
    !
    integer(FINT), intent(in) :: N
    character(FLCHARS), dimension(:), intent(in) :: cmdstr
    character(*), intent(in) :: arg
    !
    integer(FINT) :: i
    !
    pcmd_error = 0
    pcmd_checkargn = 0
    !
    ! Search for arg
    !
    do i=1,N
       if(arg.eq.cmdstr(i)) then
          pcmd_checkargn = pcmd_checkargn + 1
       endif
    enddo

  end function pcmd_checkargn

!H 
!H   subroutine pcmd_checkargp(N,cmdstr,arg,pars)
!H 
!H   
!H   This routine search for the given argument arg
!H   and return the the number of presence on the command line.
!H   It checks for parameters and return ' ' empty string if not
!H   
!H   
!H   
  integer(FINT) function pcmd_checkargp(N,cmdstr,arg,pars)
    !
    integer(FINT), intent(in) :: N
    character(FLCHARS), dimension(:), intent(in) :: cmdstr
    character(*), intent(in) :: arg
    character(FLCHARS), intent(out) :: pars
    !
    integer(FINT) :: i
    !
    pcmd_error = 0
    pcmd_checkargp = 0
    pars = ' '
    !
    ! Search for arg
    !
    do i=1,N
       if(arg.eq.cmdstr(i)) then
          pcmd_checkargp = pcmd_checkargp + 1
          if(i.lt.N) then
             select case (cmdstr(i+1)(1:1))
                case ('-') ! check if std input required
                   if(len_trim(cmdstr(i+1)).eq.1) then
                      pars = adjustl(pars)
                      pars = pars(1:len_trim(pars))//' -'
                   endif
                !case (' ','$') ! it not a file name
                case default
                   pars = adjustl(pars)
                   pars = pars(1:len_trim(pars))//' '//cmdstr(i+1)
                end select
             
          endif
       endif
    enddo

  end function pcmd_checkargp

  subroutine pcmd_iargc(N)
    !
    integer(FINT) , intent(out) :: N
    integer :: iargc
    !
    pcmd_error = 0
    N = iargc()
    return
    !
  end subroutine pcmd_iargc
  !
  !
  !
  subroutine pcmd_getarg(N,StrD)
    integer(FINT), intent(in) :: N
    character(FLCHARS),dimension(:) ,intent(out) :: StrD
    !
    integer(FINT) :: i
    !
    pcmd_error = 0
    !
    if(size(StrD,1).lt.N) then
       pcmd_error = -1
       return
    endif
    !
    do i=1,N
       call getarg(i,StrD(i))
    enddo
    !
    return
  end subroutine pcmd_getarg
!
end module pcmdline
