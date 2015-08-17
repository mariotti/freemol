!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H MODULE MESSAGES Copyright(c) Freemol2000 F.Mariotti
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H $Id: messages.F90,v 1.10 2002/07/21 18:21:36 fmariot Exp $
!H-----------------------------------------------------------------------------
!H 
!H 
!
module messages
  !
  use vartypes
  !
  implicit none
  private
!
!H
!H-----------------------------------------------------------------------------
!H      Defined Functions
!H-----------------------------------------------------------------------------
!H
!H      messages_init
!H                    This function initialize the messages module:
!H                    In input expect four integers identifing the
!H                    fortran units where to write messages. (see units)
!H                    and an extra optional argument to be defined.
!H      message
!H                    It is the actual messages facility.
!H                    In input get the type of the message and
!H                    the message string itself.
!H
!H      message_value
!H                    same as the message routine but
!H                    accept a value to be printed as well.
!H                    A tipical example can be:
!H                    call message_value(MESWARN,'current too low:',current)
!H
!H      message_setmessage
!H                    Define a storable message. It is supposed to be used
!H                    to define a general message wich will describe the real
!H                    problem and not the dedicated problem of the procedure.
!H                    I mean .. if you get it . store it there .. tobe retrived.
!H
!H-----------------------------------------------------------------------------
!H      Public variables ( ..... is it an incomplete list?) anyway.
!H-----------------------------------------------------------------------------
!H
!H     messages_lastcode
!H                   It contains the last colle code in the messages module.
!H     messages_lastmessage
!H                   It contains the last message give from the calling routine.
!H NOTE: Actually the message_setmessage has been thoguht to clarify this problem.
!H
!H-----------------------------------------------------------------------------
!H      Fortran Units
!H-----------------------------------------------------------------------------
!H
!H      fhout
!H             It is supposed to be the output file handler.
!H             write here what you want to show in the output.
!H             See also default behaviour.
!H      fhmes
!H             It is an axtra handler used for printing messages.
!H             It is the plain messages utility, somethig like:
!H             I did it!
!H      fhlog
!H             It must be used to record in a log file status infos.
!H             It is supposed that you log states of the program.
!H             You can have your log file handler.
!H      
!H      fherr
!H             File handler where to report errors.
!H      
!H
!H-----------------------------------------------------------------------------
!H      how it works and defaults
!H-----------------------------------------------------------------------------
!H
!H      The standard use will be to call message as:
!H      
!H      message(<mescode>,truemessage)
!H      
!H      The mescode will define where you what to print the message
!H      as defined by the init procedure. the mescode is and integer
!H      value and it is mapped to a set of parameters.
!H      
!H      The general behaviour will be defined by the integer mescode
!H      and the imode.
!H      
!H      The messages are printed in the given file handler
!H      in the form:
!H      MESXXX: <message>
!H      I suggest a 71 max length for the message in order to have 80 column
!H      result string.
!H      If it is printed in the output the char '#' is added before the line!
!H      A plain message to the output does not contain the string 'MESXXXX:'
!H      but only the first 'comment out' character '#'.
!H      
!H      An undefined message type is anyway printed to out err and mes fh.
!H      
!H-----------------------------------------------------------------------------
!H      Parameters defining different message code
!H-----------------------------------------------------------------------------
!H
!H      These are integer codes to give to message routine.
!H
!H      MESPLAI
!H               Plain Message. Print this message in the message file.
!H
!H      MESLOG
!H               Log Message. Print this message in the logfile.
!H
!H      MESOUT
!H               Output Message. Print this message in the output.
!H
!H      MESWARN
!H               Warning Message. Print this message in the output.
!H
!H      MESERRO
!H               Error Message. Print this message in the output and in
!H               the error files. but do not kill the procedure.
!H      
!H      MESKILL
!H               Print the message and kill the program.
!H               This mesage is printed in the output and in the error
!H               file. DO NOT GIVE BACK THE CONTROL. WE WOULD LIKE TO
!H               KILL THE PROGRAM!
!H
!H      MESINTR
!H               This error must be handled throght interactive mode.
!H               Interactive mode:
!H               if the program is running in a interactive mode
!H               the program will produce a prompt asking
!H               for the next action. Not Yet implemented.
!H               If the program is running in NOT-interactive mode
!H               a message is produced and default routine will be called.
!H               NOT IMPLEMENTED. This Should be used by graphical procedures!!
!H               ...and we have to handle it!!
!H
!H     MESDEBG   Messages for debugging porposes: activated only if
!H               debug lgical is true in this module.
!H
!H
!H-----------------------------------------------------------------------------
!H TODO
!H-----------------------------------------------------------------------------
!H
!H This module never produce an error! We will have to introduce a check
!H to determine if actually it is possible to write in the specified
!H unit.
!H-----------------------------------------------------------------------------
!H
!
  !Parameters
  !----------
  ! We set-up messages with a prestring of 9 chars
  ! The total lenght of a message cannot go over of FLCHARS-9
  ! So: PreString
  !--------------
  integer(FINT), parameter :: MESPRESTR = 9
  !Public parameters
  !-----------------
  !
  ! are defined as bit mask because maybe in the future
  ! we will need it....
  integer(FINT), parameter, public :: MESERRO = 0
  integer(FINT), parameter, public :: MESPLAI = 1
  integer(FINT), parameter, public :: MESWARN = 2
  integer(FINT), parameter, public :: MESINTR = 4
  integer(FINT), parameter, public :: MESOUT  = 8
  integer(FINT), parameter, public :: MESLOG  = 16
  integer(FINT), parameter, public :: MESKILL = 32
  integer(FINT), parameter, public :: MESDEBG = 64
  !
  ! the debug flag: please change this if you feel like.
  !---------------
  logical, save, public :: messages_debug = .false.
  !
  ! The mode flag: not used now
  !--------------
  integer(FINT) :: msg_mode
  !
  ! Public Vars
  !------------
  integer(FINT), save, public :: messages_lastcode
  character(FLCHARS),save, public :: messages_lastmessage
  character(FLCHARS),save, public :: messages_lastsetmessage
  !
  !Private Vars
  !------------
  ! check if it is initialized
  logical, save, private :: messages_isinit = .false.
  ! units defaults
  integer(FINT), save, private :: fhlog = 6
  integer(FINT), save, private :: fhout = 6
  integer(FINT), save, private :: fhmes = 6
  integer(FINT), save, private :: fherr = 0
  integer(FINT), save, private :: lastcode = 0

  !Subroutines
  public messages_init
  public message
  public message_value
  public message_setmessage
  !
  public message_comment
!
!H
!H-----------------------------------------------------------------------------
!H      messages_init interface
!H-----------------------------------------------------------------------------
!H      
!H      This interface resolves into two subroutine procedures:
!H      messages_initn()
!H      messages_initw(<pars>)
!H      The first one is called with no argument and sets file handler
!H      defaults if are not given.
!H      The second one sets the file handlers.
!H
!H      TODO:
!H      For future extension we should set up a new routine
!H      which is called with a vector or better a hash table.
!H      The concept is to write what and where.
!H      
!H-----------------------------------------------------------------------------
!H      
!H
!
  interface messages_init
     module procedure messages_initw
     module procedure messages_initn
  end interface
  !
  !
  !
  interface message_value
     module procedure message_valuei
     module procedure message_valuer
     module procedure message_valuec
  end interface
  !
contains
!
!H
!H-----------------------------------------------------------------------------
!H Subroutine messages_initn
!H-----------------------------------------------------------------------------
!H
!
  subroutine messages_initn(vdebug)
    !
    logical, optional, intent(in) :: vdebug
    !
    if(present(vdebug)) then
       messages_debug = vdebug
    end if
    ! we simply call the true routine
    call messages_initw(6,6,0,6,1,messages_debug)
    !
  end subroutine messages_initn
!
!H
!H-----------------------------------------------------------------------------
!H Subroutine messages_initw
!H-----------------------------------------------------------------------------
!H
!
  subroutine messages_initw(iumes,iuout,iuerr,iulog,imode,vdebug)
!
!H
!H    Set Up init value for Messages
!H
!H    INPUT: iuout    file unit for output
!H           iuerr    file unit for errors
!H           iulog    file unit for log
!H           iumes    file unit for messages
!H           imode    Not used now (optional)
!H
!H
!H
!H
!W    Not used mode
!W    Do not test initial units and there is not a default
!H-----------------------------------------------------------------------------
!H
!
    integer(FINT), intent(in) :: iuout
    integer(FINT), intent(in) :: iuerr
    integer(FINT), intent(in) :: iumes
    integer(FINT), intent(in) :: iulog
    integer(FINT), optional, intent(in) :: imode
    logical, optional, intent(in) :: vdebug
    !
    ! Set files handlers
    !-------------------
    fhmes = iumes
    fhout = iuout
    fherr = iuerr
    fhlog = iulog
    !
    messages_isinit = .true.
    !
    ! We nullify some data
    !---------------------
    messages_lastcode = -1
    messages_lastmessage = ''
    messages_lastsetmessage = ''
    !
    ! We check the debug flag
    !------------------------
    if(present(vdebug)) then
       messages_debug = vdebug
    end if
    if(messages_debug) then
       call message(MESLOG,'Set debug to true.')
    endif
    !
    ! We check the imode flag
    !------------------------
    msg_mode = 1
    if(present(imode)) then
       msg_mode = imode
    end if
    !
    return
    !
  end subroutine messages_initw
!
!H
!H-----------------------------------------------------------------------------
!H Subroutine message
!H-----------------------------------------------------------------------------
!H
!
  subroutine message(icode,cmes,imode)
    !
    !H
    !H    This subroutine will write the message in the appropriate
    !H    handler
    !H
    !H    icode   code of the message see parameter card
    !H    cmes    message string (no implicit dimension but it is better
    !H            if cmes(72).
    !H    imode   Not used now ... for future porpose. (optional)
    !H
    !H
    character(*), intent(in) :: cmes
    integer(FINT), intent(in) :: icode
    integer(FINT), intent(in), optional :: imode
    !The Local Var
    character(FLCHARS) :: lcmes
    integer(FINT) :: meslen
    !
    ! What ever is going on we set the last call
    !-------------------------------------------
    messages_lastcode = icode
    !
    ! we check out of boundaries stuff: should be checked before
    !                                   but we are the message routine.
    !---------------------------------
    if(len(cmes).gt.FLCHARS-MESPRESTR) then
       lcmes=cmes(1:FLCHARS-MESPRESTR)
    else
       lcmes=cmes(1:)
    end if
    !
    if(len(cmes).gt.FLCHARS) then
       messages_lastmessage = lcmes(1:LCHARS)
    else
       messages_lastmessage = lcmes
    endif
    meslen = len_trim(lcmes)
    !
    !
    ! We want to be sure we have the correct I/O units
    !-------------------------------------------------
    if(.not.messages_isinit) call messages_init()
    !
    ! We save how we are called
    !--------------------------
    lastcode = icode
    !
    ! Select the call
    !----------------
    if (icode.eq.MESPLAI) then
       write (fhmes, '(A9,(A))') '#MESSAGE:', lcmes(1:meslen)
    elseif (icode.eq.MESWARN) then
       write (fhout, '(A9,(A))') '#MESWARN:', lcmes(1:meslen)
    elseif (icode.eq.MESERRO) then
       write (fhout, '(A9,(A))') '#MESERRO:', lcmes(1:meslen)
       write (fherr, '(A9,(A))') '#MESERRO:', lcmes(1:meslen)
    elseif (icode.eq.MESKILL) then
       write (fhout, '(A9,(A))') '#MESKILL:', lcmes(1:meslen)
       write (fherr, '(A9,(A))') '#MESKILL:', lcmes(1:meslen)
       !the kill procedure is not yet implemented so we just stop it!
       write (fhout, '(A)') ' KILL GOES WITH STOP NOW '
       write (fherr, '(A)') ' KILL GOES WITH STOP NOW '
       stop 1
    elseif (icode.eq.MESINTR) then
       write (fhout, '(A9,(A))') '#MESINTR:', lcmes(1:meslen)
       write (fhmes, '(A9,(A))') '#MESINTR:', lcmes(1:meslen)
       write (fherr, '(A9,(A))') '#MESINTR:', lcmes(1:meslen)
       write (fhout, '(A)') ' MESSAGE TYPE NOT YET IMPLEMENTED '
       write (fhmes, '(A)') ' MESSAGE TYPE NOT YET IMPLEMENTED '
       write (fherr, '(A)') ' MESSAGE TYPE NOT YET IMPLEMENTED '
    elseif (icode.eq.MESLOG) then
       write (fhlog, '(A9,(A))') '# MESLOG:', lcmes(1:meslen)
    elseif (icode.eq.MESOUT) then
       write (fhout, '(A1,1X,(A))') '#', lcmes(1:meslen)
    elseif (icode.eq.MESDEBG) then
       if (messages_debug) then
          write (fhout, '(A9,(A))') '#MESDEBG:', lcmes(1:meslen)
          write (fherr, '(A9,(A))') '#MESDEBG:', lcmes(1:meslen)
       endif
    else
       write (fhout, '(A)') '# UNKNOWN MESSAGE TYPE'
       write (fhmes, '(A)') ' UNKNOWN MESSAGE TYPE'
       write (fherr, '(A)') ' UNKNOWN MESSAGE TYPE'
       write (fhout, '(A9,(A))') '#  ERROR:', lcmes(1:meslen)
       write (fhmes, '(A9,(A))') '#  ERROR:', lcmes(1:meslen)
       write (fherr, '(A9,(A))') '#  ERROR:', lcmes(1:meslen)
    endif
    !
    return
  end subroutine message
!
!H
!H-----------------------------------------------------------------------------
!H Subroutine message_valuei
!H-----------------------------------------------------------------------------
!H
!
  subroutine message_valuei(icode,cmes,ival,imode)
    !
    use strtools
    !
    !H
    !H    This subroutine call the message facility
    !H    after builted the appropriate string.
    !H WARNING: This routine cat too long messages.
    !H The longest supported message is LLENM characters!
    !H
    !
    integer(FINT), intent(in) :: icode
    character(*), intent(in) :: cmes
    integer(FINT), intent(in) :: ival
    integer(FINT), intent(in), optional :: imode
    !
    integer(FINT) :: i
    character(FLCHARS) :: newmessage
    character(FLCHARS) :: val_frmt
    character(FLCHARS) :: val_data
    !
    i = str_getfrmt(val_frmt,ival)
    write(val_data,val_frmt(1:len_trim(val_frmt))) ival
    newmessage = cmes(1:len_trim(cmes))//' '//val_data
    !
    ! than we call the normal message routine
    if (present(imode)) then
       call message(icode,newmessage,imode)
    else
       call message(icode,newmessage)
    end if
    !
    return
  end subroutine message_valuei
  !
!
!H
!H-----------------------------------------------------------------------------
!H Subroutine message_valuer
!H-----------------------------------------------------------------------------
!H
!
  subroutine message_valuer(icode,cmes,rval,imode)
    !
    use strtools
    !
    !H
    !H    This subroutine call the message facility
    !H    after builted the appropriate string.
    !H WARNING: This routine cat too long messages.
    !H The longest supported message is LLENM characters!
    !H
    !
    character(*), intent(in) :: cmes
    integer(FINT), intent(in) :: icode
    real(FREAL), intent(in) :: rval
    integer(FINT), intent(in), optional :: imode
    !
    integer(FINT) :: i
    character(FLCHARS) :: newmessage
    character(FLCHARS) :: val_frmt
    character(FLCHARS) :: val_data
    !
    i = str_getfrmt(val_frmt,rval)
    i = len(val_frmt)
    val_frmt = '('//val_frmt(1:i)//')'
    write(val_data,val_frmt) rval
    newmessage = cmes(1:len_trim(cmes))//' '//val_data
    !
    ! than we call the normal message routine
    if (present(imode)) then
       call message(icode,newmessage,imode)
    else
       call message(icode,newmessage)
    end if
    !
    return
  end subroutine message_valuer
  !
!
!H
!H-----------------------------------------------------------------------------
!H Subroutine message_valuec
!H-----------------------------------------------------------------------------
!H
!
  subroutine message_valuec(icode,cmes,cval,imode)
    !
    use strtools
    !
    !H
    !H This subroutine call the message facility
    !H after builted the appropriate string.
    !H WARNING: This routine cat too long messages.
    !H The longest supported message is LLENM characters!
    !H
    !
    character(*), intent(in) :: cmes
    integer(FINT), intent(in) :: icode
    character(*), intent(in) :: cval
    integer(FINT), intent(in), optional :: imode
    !
    integer(FINT) :: i
    character(FLCHARS) :: newmessage
    !
    i = len(cval)
    newmessage = cmes(1:len_trim(cmes))//cval(1:i)
    !
    ! than we call the normal message routine
    if (present(imode)) then
       call message(icode,newmessage,imode)
    else
       call message(icode,newmessage)
    end if
    !
    return
  end subroutine message_valuec
  !
!
!H
!H-----------------------------------------------------------------------------
!H Subroutine message_setmessage
!H-----------------------------------------------------------------------------
!H
!
  subroutine message_setmessage(str)
    !
    character(*), intent(in) :: str
    !
    if(len(str).gt.FLCHARS) then
       messages_lastmessage = str(1:LCHARS)
    else
       messages_lastmessage = str
    endif
    !
    !
  end subroutine message_setmessage
!
!H
!H-----------------------------------------------------------------------------
!H Subroutine message_comment(iun,cmes)
!H-----------------------------------------------------------------------------
!H
!
  subroutine message_comment(iun,cmes)
!
!H
!H Write a comment in unit iun. Comments are '#' lines
!H
!
    integer(FINT), intent(in) :: iun
    character(*), intent(in) :: cmes
    !
    character(FLCHARS) :: lcmes
    !
    if(len(cmes).gt.FLCHARS-2) then
       lcmes = cmes(1:FLCHARS-2)
       lcmes = adjustl(lcmes)
    else
       lcmes = cmes
    endif
    !
    write(iun,'(A2,A)') "# ",lcmes(1:len_trim(lcmes))
    !
  end subroutine message_comment
  !
end module messages
!
!
!
