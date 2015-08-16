!
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!
!
module adfrommodule
  !
  use vartypes
  use pcmdline
  use adfrompars
  !
  implicit none
  private
  !
  ! controls
  integer(FINT), public, save :: adfrom_errors = 0
  integer(FINT), public, save :: sto_fmt = 2 !define sto putput format
  integer(FINT), public, save :: mo_fmt = 1 !define sto putput format
  !
  ! Files Data
  integer(FINT), public, save :: unit_in
  integer(FINT), public, save :: unit_out
  character(FLCHARS), public, save :: file_in
  character(FLCHARS), public, save :: file_out
  !
  ! command line stuff
  integer(FINT), save :: cmdnum
  character(FLCHARS), save, dimension(:), allocatable :: cmdpars
  !
  integer(FINT) :: irc

  !
  public adfrom_init
  public adfrom_end
  public adfstop
  public adfstopon
  public adfmessage
  !
  interface adfmessage
     module procedure adfmessagel
     module procedure adfmessagec
     module procedure adfmessagei
     module procedure adfmessager
     module procedure adfmessagern
  end interface !adfmessage
  !
contains
  !
  subroutine adfrom_init
    !
    integer(FINT) :: irc
    character(FLCHARS) :: args
    !
    ! process the command line
    !
    call pcmd_iargc(cmdnum)
    !
    allocate(cmdpars(cmdnum),STAT=irc)
    if(irc.ne.0) then
       adfrom_errors = -1
       return
    endif
    !
    call pcmd_getarg(cmdnum,cmdpars)
    !
    ! Check for help
    irc = pcmd_checkarg(cmdnum,cmdpars,'-h')
    if(irc.gt.0) then
       write (*,*) ' Write the help!!'
       stop 1
    endif
    !
    ! set default file names and read command line
    unit_in = 11
    unit_out = 12
    file_in = 'TAPE21'
    file_out = 'orbital.mld'
    call pcmd_getio(cmdnum,cmdpars,file_in,file_out)
    !catch stdin and stdout and open the files
    if(file_in.eq.'-') then
       unit_in = 5
       write (0,*) 'Error: use messages in this code.'
       write (0,*) 'Error: Input from stdin not yet implemented'
       stop 1
    endif
    if(file_out.eq.'-') then
       unit_out = 6
    else
       open(unit_out,FILE=file_out,STATUS='UNKNOWN',IOSTAT=irc)
       if(irc.ne.0) then
          write (*,*) ' Error: opening the output file '
       endif
    endif
    !
    ! Check for sto outpout format
    irc = pcmd_checkarg(cmdnum,cmdpars,'-sto',args)
    if(irc.gt.0) then
       if(len_trim(args) < 2 .or. args(1:1).eq.' ') then
          call AdfStop('Error: no arguments to -sto option')
       endif
       select case (args)
       case ('molden','cart','sto','xsto')
          sto_fmt = STO_FMT_CART_LONG !STO_FMT_MOLDEN
       case ('cartshort')
          sto_fmt = STO_FMT_CART_SHORT
       case ('cartzc','cartzcontraction')
          sto_fmt = STO_FMT_CART_ZCONTR
       case ('cartshortz','cartszc')
          sto_fmt = STO_FMT_CART_SZCONTR
       case ('spher','ysto','sphersto')
          sto_fmt = STO_FMT_SPHER_LONG
       case ('sphershort')
          sto_fmt = STO_FMT_SPHER_SHORT
       case ('spherzc','spherzcontraction')
          sto_fmt = STO_FMT_SPHER_ZCONTR
       case ('sphershortz','spherszc')
          sto_fmt = STO_FMT_SPHER_SZCONTR
       case default
          call AdfStop('Error: not a valid argument to -sto option')
       end select
    endif
    !
    irc = pcmd_checkarg(cmdnum,cmdpars,'-mo',args)
    if(irc.gt.0) then
       if(len_trim(args) < 2 .or. args(1:1).eq.' ') then
          call AdfStop('Error: no arguments to -mo option')
       endif
       select case (args)
       case ('molden','mo','molines')
          mo_fmt = MO_FMT_LINES_CONTR !MO_FMT_MOLDEN
       case ('listc','listcontr')
          mo_fmt = MO_FMT_LIST_CONTR
       case default
          call AdfStop('Error: not a valid argument to -mo option')
       end select
    endif
    !
  end subroutine adfrom_init
  !
  !
  !
  subroutine adfrom_end()
    !
    if(allocated(cmdpars)) deallocate(cmdpars)
    !
  end subroutine adfrom_end
  !
  subroutine adfmessagel(mess)
    !
    logical, intent(in) :: mess
    !
    write(*,*) 'MES: ',mess
    !
  end subroutine adfmessagel
  !
  subroutine adfmessagec(mess)
    !
    character(*), intent(in) :: mess
    !
    write(*,*) 'MES: '//mess
    !
  end subroutine adfmessagec
  !
  subroutine adfmessagei(mess)
    !
    integer(FINT), intent(in) :: mess
    !
    write(*,*) 'MES: ',mess
    !
  end subroutine adfmessagei
  !
  subroutine adfmessager(mess)
    !
    real(FREAL), intent(in) :: mess
    !
    write(*,*) 'MES: ',mess
    !
  end subroutine adfmessager
  !
  subroutine adfmessagern(mess,width)
    !
    real(FREAL), intent(in), dimension(:) :: mess
    integer(FINT), intent(in), optional :: width
    !
    integer(FINT) :: iwt
    integer(FINT) :: isz,iln,i
    character(FLCHARS) :: frmt
    !
    if(present(width)) then
       iwt = width
    else
       iwt = 6
    endif
    !
    isz = size(mess,1)
    iln = isz / iwt
    if(mod(isz,iwt).gt.0) iln = iln + 1
    !
    if(iwt.gt.9) then
       write(frmt,'("1X,A4,1X,",I1,"(1X,F15.8)")') iwt
    elseif(iwt.gt.99) then
       write(0,*) 'error unsupported format'
       stop 1
    else
       write(frmt,'("1X,A4,1X,"I1,"(1X,F15.8)")') iwt
    endif
    !
    do i=1,iln
       write(*,frmt) 'MES:',mess
    enddo
    !
  end subroutine adfmessagern
  !
  subroutine adfstop(mes)
    !
    character(*) :: mes
    !
    write (0,*) mes
    stop 1
    !
  end subroutine adfstop
  !
  subroutine adfstopon(on,mes)
    !
    logical, intent(in) :: on
    character(*), intent(in) :: mes
    !
    if (on) then
       write (0,*) mes
       stop 1
    endif
    !
  end subroutine adfstopon
  !
end module adfrommodule
