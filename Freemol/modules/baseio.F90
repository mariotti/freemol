!
!H
!H-----------------------------------------------------------------------------
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H MODULE baseio: Freemol by F.Mariotti: (c) F.Mariotti 2002
!H-----------------------------------------------------------------------------
!H $Id: baseio.F90,v 1.1.1.1 2009/01/12 16:56:17 mariotti Exp $
!H-----------------------------------------------------------------------------
!H 
!H This module perform the basic io process.
!H It is a part of the file manager modules system.
!H File manager moldule system is composed by:
!H basic I/O: baseio, Extended File manager: extio, Section Manager: sections.
!H 
!H This module define a manager for the file handlers and redefine the
!H standard routines open close and filestat.
!H NOT yet everything is implemented.
!H
!H This module is used directly by extio (EXTended file manager) and section
!H manager. The target of the module is to give an abstraction over fortran
!H file handlers.
!H 
!H 
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H 
!
module baseio
  !
  use vartypes
  use messages
  !
  implicit none
  private
  !
!
!H
!H-----------------------------------------------------------------------------
!H 
!H As a convenction we adopt here two differents names:
!H baseio_xxxx
!H   and
!H bio_xxxx
!H
!H Where baseio_xxxx identify public routines and
!H bio_xxxx private routines with some exceptions.
!H 
!H 
!H In the module we have:
!H Available Data are:
!H Private
!H   bio_fh                     Matrix with file handler numbers
!H                              and state value.
!H                              bio_fh(i,1) return the handle number
!H                              bio_fh(i,2) return:
!H                                          0   if handler is not used
!H                                          >0  if handler is used
!H                              Higher values for bio_fh(i,2) are reserver
!H                              for future use like file lock facility.
!H BIO_FHMIN                    default minimum value for the file handle(35)
!H                              we skip default fortran units like 0,5,6.
!H                              The idea is to reserve more units for example:
!H                              0    error device
!H                              1
!H                              2
!H                              3
!H                              4
!H                              5    input device
!H                              6    output device
!H                              7    messages device
!H                              8    log facility device
!H                              9    printer device
!H                              10
!H                              11   graphical device (remind X11)
!H                              12
!H                              13
!H                              14
!H                              15
!H                              16
!H                              17
!H                              18
!H                              19
!H                              20   frimol private atom module
!H                              21   frimol private molecule module
!H                              22   frimol private molecules module
!H                              23   frimol private 
!H                              24   frimol private 
!H                              25   frimol private 
!H                              26   frimol private 
!H                              27   frimol private 
!H                              28   frimol private 
!H                              29   frimol private 
!H                              30
!H                              31
!H                              32
!H                              33
!H                              34
!H                              These devices must be defined elsewhere!!!
!H                              Inside the Freemol Package:
!H                              Do we need to extend this??
!H                              In principle forever in practice we want
!H                              to keep an easy use!!
!H BIO_FHMAX                    default maxmum value for the file handle(90)
!H                              It is system dependent!! Is 90 enough?
!H bio_fhbeg                    Actual used variable for fh min
!H bio_fhend                    Actual used variable for fh max
!H bio_fhtot                    Total number of available units
!H bio_fhtop                    The latest allocated fh
!H bio_fhfst                    The first lower available fh
!H
!H bio_error                    integer error flag ( >0 if errors)
!H                              It is declared private you have to use
!H                              The new function: baseio_geterror().
!H 
!H CHECK THESE COMMENTS!!!
!H Available procedures are:
!H Private (sub or fun)
!H 
!H 
!H 
!H Public (fun or sub)
!H baseio_init                  Initialize the baseio module
!H baseio_open                  Open a file
!H baseio_close                 Close the file unit
!H baseio_kill                  Shutdown the module
!H bio_getfreefh                We get the next file handler
!H bio_removefh                 We remove a file handler from the table
!H 
!H Following the Computationl chemistry ideas we have:
!H 
!H baseio_open_scratch          We want a temporary storage as a work space
!H baseio_open_temp             We want a temporary storage to be used
!H                              for future reading even on a restart procedure.
!H                              Actually a file on a disk.
!H These two last function require a comment:
!H We implement these funtion here because we guess that are required
!H at low level implemetation in order to achive speed up.
!H It would be possible to have equivalent funtion on the
!H extio module but maybe a faster implemetation is usefull.
!H We have two ways to implement it and hopefully we would like to
!H be implemetation indipendent.
!H The first trivial solution is to call baseio_open(...) with
!H appropriate parameters. The second solution is to call
!H the function on its own: it anable for better optimization.
!H With the second solution we will be able to decide where our scratch
!H (temporary) storage has to go: memory vs disk for scrath files.
!H Maybe it gets work to implement it andto optimize it but
!H I would like to keep this possibility open.
!H So the actual implementation uses baseio_open!!!
!H 
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H 
!H TODO
!H A lot!
!H 
!H   - Implement correctly the errors codes.
!H   - Implements these routines
!H     baseio_bwrite() Write whatere in binary format
!H     baseio_bread() Write whatere in ascii format
!H     baseio_awrite() Read whatere in binary format
!H     baseio_aread() Read whatere in ascii format
!H     
!H   
!H   
!H   
!H   
!H   
!H 
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H   
!
  !
  ! Parameters
  !-----------
  ! Maybe we need a different approach to this.
  ! But right now seems to work!
  !--------------------------------------------
  character(LCHARS), public, parameter :: baseio_statuskeys = 'APPEND OLD NEW UNKNOWN REPLACE SCRATCH'
  character(LCHARS), public, parameter :: baseio_formatkeys = 'FORMATTED UNFORMATTED'
  character(LCHARS), public, parameter :: baseio_accesskeys = 'SEQUENTIAL DIRECT'
!
  ! private parameters:
  !--------------------
  ! We would like to reseve some units (up to 35)
  ! And we do not want to go over some OS/COMP limits
  !-------------------------------------------------------
  integer(FINT), private, parameter :: BIO_FHMIN = 35
  integer(FINT), private, parameter :: BIO_FHMAX = 99
  !
  ! private variables:
  !-----------------------
  ! DEBUG: check this initialization
  !---------------------------------
  logical, save, private :: bio_isinit = .false.
  integer(FINT), save, private :: bio_fhbeg
  integer(FINT), save, private :: bio_fhend
  integer(FINT), save, private :: bio_fhtot
  integer(FINT), save, private :: bio_fhfst
  integer(FINT), save, private :: bio_fhtop
  !
  ! The file handle matrix: The Table!
  !-----------------------------------
  integer(FINT), save, private, dimension(:,:), allocatable :: bio_fh
  !
  ! error variables
  !----------------
  integer(FINT), save, private :: bio_error
  character(FLCHARS), save, private :: bio_error_message
  !
  ! public parameters:
  !-----------------------
  ! integer(FINT), save, public :: baseio_error
  !
  ! public procedures
  !------------------
  public baseio_init
  public baseio_open
  public baseio_close
  public baseio_isinit
  public baseio_geterror
  public baseio_open_scratch
  public baseio_open_temp
  public baseio_bwrite
  public baseio_bread
  public baseio_reclen
  !
  ! Interfaces
  !-----------
  interface baseio_bwrite
     module procedure baseio_bwrite_in ! Write an integer number
     module procedure baseio_bwrite_iv ! Write an integer vector
     module procedure baseio_bwrite_im ! Write an integer matrix
     !module procedure baseio_bwrite_iu ! Write an integer multid
     module procedure baseio_bwrite_rn
     module procedure baseio_bwrite_rv
     module procedure baseio_bwrite_rm
     !module procedure baseio_bwrite_ru
     module procedure baseio_bwrite_cn
     module procedure baseio_bwrite_cv
     module procedure baseio_bwrite_cm
     !module procedure baseio_bwrite_cu
  end interface !baseio_bwrite
  !
  !
  interface baseio_bread
     module procedure baseio_bread_in ! Read an integer number
     !module procedure baseio_bread_iv ! Read an integer vector
     !module procedure baseio_bread_im ! Read an integer matrix
     !module procedure baseio_bread_iu ! Read an integer multid
     !module procedure baseio_bread_rn
     !module procedure baseio_bread_rv
     !module procedure baseio_bread_rm
     !module procedure baseio_bread_ru
     !module procedure baseio_bread_cn
     !module procedure baseio_bread_cv
     !module procedure baseio_bread_cm
     !module procedure baseio_bread_cu
  end interface !baseio_bread
  !
  !
  interface baseio_reclen
     !module procedure baseio_reclen_in ! Get reclen for an integer number
     module procedure baseio_reclen_iv ! Get reclen for an integer vector
     module procedure baseio_reclen_im ! Get reclen for an integer matrix
     module procedure baseio_reclen_iu ! Get reclen for an integer multid
     !module procedure baseio_reclen_rn
     module procedure baseio_reclen_rv
     module procedure baseio_reclen_rm
     module procedure baseio_reclen_ru
     !module procedure baseio_reclen_cn
     module procedure baseio_reclen_cv
     module procedure baseio_reclen_cm
     module procedure baseio_reclen_cu

  end interface !baseio_reclen
  !
  !
  ! Error Interface
  !----------------
  ! We define two error interfaces. Each one return the latest error
  ! code or string value. The String error function can be used to inqure
  ! about a given error code.
  !----------------------------------------------------------------------
  interface baseio_geterror
     module procedure baseio_geterror_i
     module procedure baseio_geterror_c
  end interface !baseio_geterror
  !
  ! public procedures but deprecated
  !---------------------------------
  ! These procedures are public but are deprecated.
  ! These procedures can provide a tool to handle
  ! file handlers when you do not want to use the provided
  ! baseio_open, baseio_close and other baseio_procedures.
  ! It is suppose you know what you are doing but
  ! it is not a garanty to still find these procedure in the future.
  ! Most probably we will provide a new module: filehandler module?
  ! TODO: crate the file handler module ;)
  !-----------------------------------------------------------------
  public bio_getfreefh
  public bio_removefh
  !
contains
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_geterror_i()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_geterror_i()
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
    baseio_geterror_i = bio_error
    !
  end function baseio_geterror_i
!
!H-----------------------------------------------------------------------------
!H
!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H character(FLCHARS) function baseio_geterror_c()
!H-----------------------------------------------------------------------------
!H
!
  character(FLCHARS) function baseio_geterror_c(code)
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
       baseio_geterror_c = "No Errors."
    case(1_FINT)
       baseio_geterror_c = "Not documented error or generic error.."
    case(2_FINT)
       baseio_geterror_c = "Double Error: Why did you get this message?."
    case default
       baseio_geterror_c = "Not documented error or generic error.."
    end select
    !
    return
    !
  end function baseio_geterror_c
!
!H
!H-----------------------------------------------------------------------------
!H
!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H Function baseio_init(fhbeg,fhend)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_init(fhbeg,fhend)
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the baseio module and get some optional parameters:
!H Up to date we have:
!H fhbeg                        define the first handle available
!H fhend                        define the last handle available
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    integer(FINT), optional, intent(in) :: fhbeg
    integer(FINT), optional, intent(in) :: fhend
    !
    integer(FINT) :: irc
    integer(FINT) :: i
    !
    ! We don't want to be initialized twice!
    ! DEBUG: NOTE: do we need to call this initialisation
    ! in more part of the code?
    ! At least we give the function: baseio_isinit()
    ! to be inqired! So go on and produce an error!
    ! As being Init we return error code -1 on error
    ! inquire with baseio_error() funtion
    !-----------------------------------------------
    if(bio_isinit) then
       bio_error = 1
       baseio_init = -1
       ! 
       ! We use messages to report the double call
       !-------------------------------------------
       call message(MESDEBG,"[baseio_init]: Called twice, Internal Error.")
       return
    end if
    !
    ! Some reset
    !-----------
    baseio_init = 0
    bio_error = 0
    bio_isinit = .false.
    bio_fhbeg = BIO_FHMIN
    bio_fhend = BIO_FHMAX
    bio_fhtop = BIO_FHMIN
    bio_fhfst = BIO_FHMIN
    !
    ! We require different bounding.
    !-------------------------------
    if(present(fhbeg)) then
       if(fhbeg.gt.BIO_FHMIN) then
          bio_error = 1
          baseio_init = -1
          call message(MESDEBG,&
               &"[baseio_init]: Too low file handle, Internal Error.")
          return
       endif
       bio_fhbeg = fhbeg
    endif
    !
    if(present(fhend)) then
       if(fhbeg.gt.BIO_FHMAX) then
          bio_error = 1
          baseio_init = -1
          call message(MESDEBG,&
               &"[baseio_init]: Too High file handle, Internal Error.")
          return
       endif
       bio_fhend = fhend
    endif
    !
    ! We do not axcept negative values
    ! as well some incocitency!
    !---------------------------------
    if(bio_fhbeg.lt.0) then
       bio_error = 1
       baseio_init = -1
       call message(MESDEBG,"[baseio_init]: negative fh, Internal Error.")
       return
    else if (bio_fhend.lt.bio_fhbeg) then
       bio_error = 1
       baseio_init = -1
       call message(MESDEBG,"[baseio_init]: No fh window, Internal Error.")
       return
    end if
    !
    ! We get the nuber of handles we can manage.
    !-------------------------------------------
    bio_fhtot = bio_fhend - bio_fhbeg + 1
    !
    ! Something is wrong: we have allocated vectors.
    !-----------------------------------------------
    ! DEBUG: NOTE: Actually in a proper use we will never fall
    ! in this if statment. But we check it anyway:
    ! On errors and file initialization I do not care
    ! about time performance! Anyway we force deallocation!
    !---------------------------------------------------------
    if(allocated(bio_fh)) then
       call message(MESDEBG,"[baseio_init]: fh deallocation forced.")
       deallocate(bio_fh,STAT=irc)
       if(irc.ne.0) then
          call message(MESDEBG,"[baseio_init]: unable to deallocate fh.")
          bio_error = 1
          baseio_init = -1
          return
       endif
    end if
    !
    ! NOTE: In this init we assume that exist a corrispondence
    !       between the fh unit id and the index of the vector
    !       which contains the information.
    !       Actually this is the way how it is implemented in
    !       the init function today, but do not trust it.
    !       A possible future version of this function
    !       or of a file handlers module will be able
    !       to select a set of file handlers not necessarly
    !       sequential.
    !       SEE: the implementation of bio_getfreefh and bio_freefh.
    !       
    !       
    !---------------------------------------------------------------
    !
    ! We actually allocate the handlers vector.
    !------------------------------------------
    ! This allocation is a key point!
    ! Right now we allocate a vector defined from
    ! the starting handler(it can be >0) to the final one.
    ! This can be changed in the future where we do
    ! not assume a corrispondence between file handle
    ! the index of the matrix.
    !-----------------------------------------------------
    allocate(bio_fh(bio_fhbeg:bio_fhend,2),STAT=irc)
    if(irc.ne.0) then
       call message(MESERRO,"[baseio_init]: unable to allocate fh.")
       bio_error = 1
       baseio_init = -1
       return
    end if
    !
    ! We set up the handler as specifyied before:
    ! with corrispondence file handle<->units.
    !--------------------------------------------
    bio_fh(:,1) = (/ (i,i=bio_fhbeg,bio_fhend) /)
    bio_fh(:,2) = 0
    !
    ! We are initialized.
    !--------------------
    bio_isinit = .true.
    !
    ! Done!
    !------
    !
  end function baseio_init
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H Function baseio_close(iun)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_close(iun)
!
!H
!H    input IUN = file unit to close
!H
!
    ! We declare it as inout because we want to nullify it as well
    ! or in other words we want to set it to -1 ..
    integer(FINT), intent(inout) :: iun
    !
    integer(FINT) :: irc
    !
    ! We hope to get no errors.
    !--------------------------
    baseio_close = 0
    bio_error = 0
    !
    ! We actually close the unit and we check for errors.
    !----------------------------------------------------
    close(iun,IOSTAT=irc)
    if(irc.eq.0) then
       ! DEBUG: Shall we check for locking??
       !        Locking avaiability will change
       !        thi code.
       !---------------------------------------
       baseio_close = -1
       bio_error = 1
    end if
    !
    ! We remove the handler from our table.
    !--------------------------------------
    irc = bio_removefh(iun)
    if (irc.lt.0) then
       ! What's going on?
       ! We are not able to close our file ...
       ! We just need to remove from table!
       ! DEBUG: do we need this piece of code?
       ! We warn it up but we go on with less resorces "maybe"
       ! and we return error code.
       !------------------------------------------------------
       baseio_close = -1
       bio_error = 1
       call message(MESWARN,'[baseio_close]: unable to free fh on error.')
    end if
!
! Done!
!
  end function baseio_close
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function bio_removefh(iun)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function bio_removefh(iun)
    !
    integer(FINT), intent(inout) :: iun
    !
    integer(FINT) :: fhindex
    integer(FINT) :: i
    !
    ! NOTE: We do not suppose to have the iun stored
    !       in the iun position of the fh matrix.
    ! TODO: in the next picture of the program
    !       we would like to have an extra procedure
    !       a kind of locking control ...
    !       .. a ne procedure will be bio_freefh() ..
    ! WARN: We check the existence of a sequential storing
    !       of the iun IDs just to speed up the procedure ...
    !
    !---------------------------------------------------------
    !
    ! We require a null index at least.
    !----------------------------------
    fhindex = -1
    !
    bio_removefh = iun
    bio_error = 0
    !
    ! In the following procedure we trust the first
    ! matching index!
    ! WANRING: we assume an increse of the handlers
    ! values.
    !----------------------------------------------
    !
    ! If we are inside the table index.
    !----------------------------------
    if((iun.ge.lbound(bio_fh,1)).and.(iun.le.ubound(bio_fh,1))) then
       ! We check first if the corrispondece file unit<->index is correct.
       !------------------------------------------------------------------
       if(iun.eq.bio_fh(iun,1)) then
          ! We got the index.
          !------------------
          fhindex = iun
       else
          ! We have to check in cases where the corrispondece
          ! file unit<->index is not correct.
          do i=lbound(bio_fh,1),ubound(bio_fh,1)
             if (iun.eq.bio_fh(i,1)) then
                fhindex = i
                ! we get the first one!
                !----------------------
                exit
             end if
          end do
       end if
    else
       ! The given unit was outside the table index
       ! so we check it unit by unit.
       !-------------------------------------------
       do i=lbound(bio_fh,1),ubound(bio_fh,1)
          if (iun.eq.bio_fh(i,1)) then
             fhindex = i
             ! we get the first one!
             !----------------------
             exit
          end if
       end do
    end if
    !
    if(fhindex.eq.-1) then
       bio_removefh = -1
       bio_error = 1
       call message(MESDEBG,"[bio_removefh]: iun index not found.")
       return
    end if
    !
    ! We got fhindex!!
    !-----------------
    !
    ! We free an index and if it is a low
    ! values we set the new lowest available index.
    !----------------------------------------------
    if(fhindex.lt.bio_fhfst) bio_fhfst = fhindex
    !
    ! We check if we have it already free!! It must not happen!
    !----------------------------------------------------------
    if(bio_fh(fhindex,2).le.0) then
       !
       ! Now in this case we suppose it is a problem
       ! because we are freeing an actually free handler!
       ! at least for our implementation!
       ! but it could not affect the application
       ! unless you did something wrog!!
       !
       ! Then we return beacause actually the iun is free at least for us
       ! TODO: check if actually the unit is used or not inquring the file!!
       !--------------------------------------------------------------------
       !
       call message(MESDEBG,&
            &'[bio_removefh]: Attempt to remove an empty handler.')
       bio_removefh = -1
       bio_error = -1
       return
    else
       ! Here we destroy the handler from the table!
       !--------------------------------------------
       bio_fh(fhindex,2) = 0_FINT
       ! In a file lock environment maybe we need something like:
       ! bio_fh(fhindex,2) = bio_fh(fhindex,2) - 1_FINT
       ! and check if it is actually free (bio_fh(fhindex,2) < 1).
       !----------------------------------------------------------
    endif
    !
  end function bio_removefh
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function bio_getfreefh(iun)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function bio_getfreefh(iun)
!
!H
!H-----------------------------------------------------------------------------
!H We want to give back a valid iun (file handler) and
!H set it up in the table.
!H 
!H-----------------------------------------------------------------------------
!H
!
    !
    integer(FINT), intent(out) :: iun
    !
    integer(FINT) :: i
    logical :: done
    !
    ! Now no errors
    !--------------
    bio_getfreefh = 0
    bio_error = 0
    done = .false.
    !
    ! We get the first free unit.
    !----------------------------
    ! This code use bio_fhfst but I am not sure it is
    ! the best solution.
    !------------------------------------------------
    if(bio_fhfst.gt.ubound(bio_fh,1)) then
       ! Too many file opened
       iun = -1
       bio_getfreefh = -1
       bio_error = 1
       call message(MESERRO,'[bio_getfreefh]: Too many files open.')
       return
    end if
    !
    ! We have it free.
    !-----------------
    iun = bio_fh(bio_fhfst,1)
    ! We set the lock flag.
    !----------------------
    ! But first we ensure it is free.
    ! This will be usefull maybe when we will
    ! implement a lock file facility.
    !----------------------------------------
    if(bio_fh(iun,2).eq.0) then
       iun = bio_fh(iun,1)
       bio_fh(iun,2) = 1
       ! This can be this = this + one
       ! It will exclude the next else if ...
       ! and we will have to check in which mode
       ! we open the file to include lock facility.
       !-------------------------------------------
       bio_getfreefh = iun
    else
       ! We got an internal error:
       ! Lock facility is not implemented so
       ! we cannot ask for two handler on the same file.
       !------------------------------------------------
       ! I removed the nullification of the iun:
       ! This because I trust the main code and
       ! Maybe the file is still used or better the
       ! variable containing the file is on use.
       ! But we produce errors code anyway.
       !-------------------------------------------
       !iun = -1
       !--------
       bio_getfreefh = -1
       bio_error = 1
       ! TODO: insert a check on the status of the file iun
       !       right now we simply lock it.
       ! We lock the unit to damaged.
       bio_fh(bio_fhfst,2) = -1
       call message(MESDEBG,'[bio_getfreefh]: File already in use.')
    end if
    ! We set the top value.
    !----------------------
    if(bio_fhfst.gt.bio_fhtop) bio_fhtop = bio_fhfst
    !
    ! We Search for the next one.
    !----------------------------
    ! We need to do it to keep bio_fhfst up to date.
    ! The search start from the next iun. This algoritm
    ! using this index (bio_fhfst) maybe is not faster
    ! but anyway the procdure will depend on the application.
    ! We check it after the definition of the top index
    ! just because it is easier.
    !--------------------------------------------------------
    bio_fhfst = bio_fhfst + 1
    do i=bio_fhfst,ubound(bio_fh,1)
       if(bio_fh(i,2).eq.0) then
          bio_fhfst = i
          exit
       end if
    enddo
    !
    return
!
! ! This code is obsolete
! ! Or better is a different algoritm.
!
!     do i=lbound(bio_fh,1),ubound(bio_fh,1)
!        ! Only 0 is a valid free iun
!        ! TODO: we will heve to implement state = state + 1
!        !       in order to get a kind of file lock ..
!        !       Maybe .. or we do it some whereelse ...
!        if(bio_fh(i,2).eq.0) then
!           iun = bio_fh(i,1)
!           bio_fh(i,2) = 1
!           ! This can be this = this + one
!           ! It will exclude the next else if in the next loop ...
!           !
!           ! We founf it .. so return
!           return
!        end if
!     end do
!     do i=lbound(bio_fh,1),ubound(bio_fh,1)
!        !We handle as well this method but we warn it up ...
!        if(bio_fh(i,2).gt.0) then
!           iun = bio_fh(i,1)
!           bio_fh(i,2) = bio_fh(i,2) + 1
!           call message(MESWARN,'[bio_getfreefh]: unstable procedure used.')
!           return
!        else
!           ! Here we return errors ...
!           iun = -1
!           bio_getfreefh = -1
!           bio_error = -1
!           return
!        end if
!     end do
!
    !
  end function bio_getfreefh
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_kill()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_kill()
!
!H
!H-----------------------------------------------------------------------------
!H This subroutine just try to close all open files
!H and to remove all allocated vectors.
!H It can be called before to kill the program
!H In case of errors.
!H-----------------------------------------------------------------------------
!H
!
    !
    integer(FINT) :: irc
    integer(FINT) :: i
    ! The usual function initialization.
    !-----------------------------------
    bio_error = 0
    baseio_kill = 0
    !
    ! We will need to run on over the table to free files.
    !-----------------------------------------------------
    !    do i=fhbeg,fhend !Old code?? is it better??
    ! Actually we check bounds.
    !--------------------------
    do i=lbound(bio_fh,1),ubound(bio_fh,1)
       if(bio_fh(i,2).gt.0) then
          !
          ! If we have any kind of allocated unit we free it.
          !--------------------------------------------------
          ! And actually we try as well to close the file.
          !-----------------------------------------------
          close(i,IOSTAT=irc)
          if(irc.ne.0) then
             !
             ! When ever we get an error we set it.
             !-------------------------------------
             ! Actually this subroutine has to work always!
             !---------------------------------------------
             bio_error = 1
             baseio_kill = -1
          endif
       endif
    enddo
    !
    ! Here we reset the module.
    ! We can think about a call of the type:
    ! baseio_init('FORCE')
    ! to reforce an initialization and to be used here
    ! but I guess better we keep these separated.
    !-------------------------------------------------
    bio_isinit = .false.
    bio_fhbeg = BIO_FHMIN
    bio_fhend = BIO_FHMAX
    bio_fhtop = BIO_FHMIN
    bio_fhfst = BIO_FHMIN
    !
    ! We deallocate the module.
    !--------------------------
    ! This save also time to check the table:
    ! we just free it.
    !----------------------------------------
    if(allocated(bio_fh)) deallocate(bio_fh)
    !
  end function baseio_kill
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H logical function baseio_isinit()
!H-----------------------------------------------------------------------------
!H
!
  logical function baseio_isinit()
!
!H
!H-----------------------------------------------------------------------------
!H We merelyreturn the initialization flag.
!H-----------------------------------------------------------------------------
!H
!
    baseio_isinit = bio_isinit
  end function baseio_isinit
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H Function baseio_open(iun,file,stat,frmt,acc)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_open(iun,file,stat,frmt,acc)
    !
    use linetools
    !
    integer(KINT), intent(out) :: iun
    character(*), intent(in) :: file
    character(*), intent(in), optional :: stat
    character(*), intent(in), optional :: frmt
    character(*), intent(in), optional :: acc
    !
!
!H
!H-----------------------------------------------------------------------------
!H         This routine open a file by file name
!H
!H            OUTPUT
!H            iun        File Handler
!H            RETVALUE   File Handler or error (<0)
!H
!H            INPUT
!H            file       Name of the file (string)
!H            stat       Open For Read Write and mode (character,optional)
!H                       Valid format: Are the same as the fortran open
!H                       STATUS.
!H
!H
!H            frmt       ASCII or BINARY file (character, optional)
!H                       Valid format: Are the same as the fortran open
!H                       FORMAT but we allow also ASCII and BINARY.
!H
!H            acc        The ACCESS method as in fortran open intrinsic.
!H
!H
!H        DEFAULT for stat and frmt the same as for fortran open intrinsic.
!H
!H
!H TODO: add append method
!H TODO: add operative system indipendece from file names with dir structure ..
!H TODO: implement a kind of locking when reading from an existing file
!H       but maybe we will do on an extramodule with use this one.
!H       Say that up to now we enable this procedure to store information
!H       on the second column of fh matrix and we don't handle them.
!H-----------------------------------------------------------------------------
!H
!
    !
    ! TODO: .. actually we don't need such a long strings ..
    character(FLCHARS) :: lfile
    character(FLCHARS) :: lstat
    character(FLCHARS) :: lfrmt
    character(FLCHARS) :: lacc
    !
    logical :: isscratch
    integer(FINT) :: ilfile
    character(LCHARS) :: keywords
    !
    integer(FINT) :: irc
    integer(FINT) :: i
    !
    ! Set Defaults and Init
    !----------------------
    ! Today defaults are: existing file for ascii read.
    !--------------------------------------------------
    isscratch = .false.
    lfile = ''
    lstat = 'OLD'
    lfrmt = 'FORMATTED'
    lacc = 'SEQUENTIAL'
    !
    baseio_open = 0
    bio_error = 0
    !
    if(len_trim(file).gt.FLCHARS) then
       bio_error = 1
       baseio_open = -1
       call message(MESERRO,'[baseio_open]: File name too long.')
       return
    endif
    !
    ! We check for the status string.
    !--------------------------------
    if(present(stat)) then
       lstat = stat
       lstat = trim(lstat)
       i = len_trim(lstat)
       if(windex(baseio_statuskeys,lstat(1:i)).le.0) then
          ! Invalid status
          !---------------
          bio_error = 1
          baseio_open = -1
          call message(MESERRO,'[baseio_open]: Not valid status code.')
          return
       end if
       ! We need to know if it is a scrathc file: no name required!
       !-----------------------------------------------------------
       !if(windex(lstat,'SCRATCH')) then
       if(windex(lstat,'SCRATCH').gt.0) then
          isscratch = .true.
       end if
       ! catch the append status because is not
       ! implemented.
       !---------------------------------------
       ! WARNING: we suppose that it is the first!!!
       !--------------------------------------------
       if(windex(baseio_statuskeys,lstat(1:i)).eq.1) then
          bio_error = 1
          baseio_open = -1
          call message(MESERRO,'[baseio_open]: function not implemented: APPEND')
       end if
    endif
    !
    ! We check for the format string.
    !--------------------------------
    if(present(frmt)) then
       lfrmt = frmt
       lfrmt = trim(lfrmt)
       i = len_trim(lfrmt)
       if(windex(baseio_formatkeys,lfrmt(1:i)).le.0) then
          ! Invalid format
          !---------------
          bio_error = 1
          baseio_open = -1
          call message(MESERRO,'[baseio_open]: Not valid format code.')
          return
       end if
    endif
    !
    ! We check for the access string.
    !--------------------------------
    !--------------------------------
    if(present(acc)) then
       lacc = acc
       lacc = trim(lacc)
       i = len_trim(lacc)
       if(windex(baseio_accesskeys,lacc(1:i)).le.0) then
          ! Invalid access
          !---------------
          bio_error = 1
          baseio_open = -1
          call message(MESERRO,'[baseio_open]: Not valid access code.')
          return
       end if
    endif
    !
    ! We get the file name.
    !----------------------
    lfile = file
    lfile = trim(lfile)
    !
    ! Have we got it?
    if((len_trim(lfile).eq.0).and..not.(isscratch)) then
       ! There is not a file name:
       ! at least there is not any string there.
       !----------------------------------------
       bio_error = 1
       baseio_open = -1
       call message(MESERRO,'[baseio_open]: No file name given.')
       return
    end if
    !
    ilfile = len_trim(lfile)
    !
    ! We are done with init
    ! Now we ask for FH
    irc = bio_getfreefh(iun)
    ! we catch errors
    if (irc.lt.0) then
       bio_error = 1
       baseio_open = -1
       call message(MESERRO,'[baseio_open]: Cannot Get Free FH.')
       return
    endif
    !We default the return value to the file handle
    baseio_open = iun
    !
    if(isscratch) then
       open(iun,status=lstat)
    else
       open(iun,file=lfile(:ilfile),&
            &status=lstat,form=lfrmt,ACCESS=lacc,IOSTAT=irc)
    end if
    if(irc.gt.0) then
       ! On error we remove the lock on the file handler
       irc = bio_removefh(iun)
       if (irc.lt.0) then
          ! What's going on?
          ! We are not able to close our file ...
          ! We warn it up but we go on with less resorces
          call message(MESWARN,'[baseio_open]: unable to free fh on error.')
       end if
       baseio_open = -1
       bio_error = 1
    end if
    !
  end function baseio_open
!
!H
!H-----------------------------------------------------------------------------
!H integer function baseio_open_scratch(iun)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_open_scratch(iun)
!
!H
!H-----------------------------------------------------------------------------
!H We open a scratch temporary file
!H-----------------------------------------------------------------------------
!H
!
    !
    integer(FINT), intent(inout) :: iun
    !
!
!H
!H-----------------------------------------------------------------------------
!H This procedure has to be used for opening temporary storage
!H files.
!H In priciple the file has not to be usefull for a restart
!H but has to be used for procedure temporary storage.
!H Even in the program module we can use this kind of file but
!H the main reason to use the opened file  is that such a file
!H has to be accessed fast at least locally (local in the procedure)!
!H
!H The actual implementation lack abit of fantasy!
!H We are going to use the standard baseio procedue but I guess I can implement
!H it in a better way. Shall we do it machine dedicated?
!H 
!H-----------------------------------------------------------------------------
!H
!
    !
    integer(FINT) :: irc
    !
    bio_error = 0
    !
    irc = baseio_open(iun,file='',stat='SCRATCH')
    baseio_open_scratch = irc
    !
  end function baseio_open_scratch
!
!H
!H-----------------------------------------------------------------------------
!H integer function baseio_open_temp(iun)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_open_temp(iun)
!
!H
!H-----------------------------------------------------------------------------
!H We open a temporary file: see discussion in the module.
!H-----------------------------------------------------------------------------
!H
!
    !
    integer(FINT), intent(inout) :: iun
    !
!
!H
!H-----------------------------------------------------------------------------
!H This procedure has to be used for opening temporary storage
!H files. Actually we differ from the scratch in the 
!H that we require it to be really fast and we
!H in principle) localozed it to a given procedure.
!H The target will be to have a scratch file in memory
!H but we let the routine to decide whether we have enogh
!H memoery or not.
!H In priciple the file has not to be usefull for a restart
!H but has to be used for procedure temporary storage.
!H
!H The actual implementation lack a bit of fantasy!
!H We are going to use the standard baseio procedue but I guess I can implement
!H it in a better way. Shall we do it machine dedicated?
!H Atually this procedure is not implemented.
!H 
!H-----------------------------------------------------------------------------
!H
!
    !
    integer(FINT) :: irc
    !
    bio_error = 0
    !
    irc = baseio_open(iun,file='',stat='SCRATCH')
    baseio_open_temp = irc
    !
  end function baseio_open_temp
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_reclen_iv(data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_reclen_iv(data)
    !
    integer(FINT), intent(in), dimension(:) :: data
    !
    !we return just the lenght
    baseio_reclen_iv=size(data)
    !
    !
  end function baseio_reclen_iv
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_reclen_im(data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_reclen_im(data)
    !
    integer(FINT), intent(in), dimension(:,:) :: data
    !
    !we return just the lenght
    baseio_reclen_im=size(data,1)*size(data,2)
    !
    !
  end function baseio_reclen_im
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_reclen_iu(data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_reclen_iu(data)
    !
    integer(FINT), pointer :: data
    !
    !we return just the lenght
    baseio_reclen_iu=0
    call message(MESKILL,"[baseio_reclen]: function not implemented")
    !
    !
  end function baseio_reclen_iu
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_reclen_rv(data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_reclen_rv(data)
    !
    real(FREAL), intent(in), dimension(:) :: data
    !
    !we return just the lenght
    baseio_reclen_rv=size(data)
    !
    !
  end function baseio_reclen_rv
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_reclen_rm(data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_reclen_rm(data)
    !
    real(FREAL), intent(in), dimension(:,:) :: data
    !
    !we return just the lenght
    baseio_reclen_rm=size(data,1)*size(data,2)
    !
    !
  end function baseio_reclen_rm
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_reclen_ru(data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_reclen_ru(data)
    !
    real(FREAL), pointer :: data
    !
    !we return just the lenght
    baseio_reclen_ru=0
    call message(MESKILL,"[baseio_reclen]: function not implemented")
    !
    !
  end function baseio_reclen_ru
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_reclen_cv(data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_reclen_cv(data)
    !
    character(FLCHARS), intent(in), dimension(:) :: data
    !
    !we return just the lenght
    baseio_reclen_cv=size(data)
    !
    !
  end function baseio_reclen_cv
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_reclen_cm(data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_reclen_cm(data)
    !
    character(FLCHARS), intent(in), dimension(:,:) :: data
    !
    !we return just the lenght
    baseio_reclen_cm=size(data,1)*size(data,2)
    !
    !
  end function baseio_reclen_cm
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_reclen_cu(data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_reclen_cu(data)
    !
    character(FLCHARS), pointer :: data
    !
    !we return just the lenght
    baseio_reclen_cu=0
    call message(MESKILL,"[baseio_reclen]: function not implemented")
    !
    !
  end function baseio_reclen_cu
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_bwrite_in(iun,data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_bwrite_in(iun,data)
    !
    integer(FINT), intent(in) :: iun
    integer(FINT), intent(in) :: data
    !
    integer(FINT) :: reclen
    !
    reclen = 1
    write(iun) 'i',kind(FBIGGEST_INT),reclen
    write(iun) data
    baseio_bwrite_in = reclen
    !
  end function baseio_bwrite_in
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_bwrite_iv(iun,data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_bwrite_iv(iun,data)
    !
    integer(FINT), intent(in) :: iun
    integer(FINT), intent(in), dimension(:) :: data
    !
    integer(FINT) :: reclen
    integer(FBIGGEST_INT) :: lsize
    !
    lsize = size(data)
    reclen = lsize
    if (.true.)then
       write(iun) 'i',kind(data(1)),lsize
       write(iun) data
    else
       !reclen = size(data) !baseio_reclen(data)
       reclen = baseio_reclen(data)
       ! we should do something different
       write(iun) 'i',kind(FBIGGEST_INT),lsize
       write(iun) data
    end if
    baseio_bwrite_iv = reclen
    !
  end function baseio_bwrite_iv
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_bwrite_im(iun,data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_bwrite_im(iun,data)
    !
    integer(FINT), intent(in) :: iun
    integer(FINT), intent(in), dimension(:,:) :: data
    !
    integer(FINT) :: reclen
    integer(FBIGGEST_INT) :: lsize
    !
    lsize = size(data,1)*size(data,2)
    reclen = lsize
    if (.true.)then
       write(iun) 'i',kind(data),lsize
       write(iun) data
    else
       !reclen = size(data) !baseio_reclen(data)
       reclen = baseio_reclen(data)
       ! we should do something different
       write(iun) 'i',kind(FBIGGEST_INT),lsize
       write(iun) data
    end if
    baseio_bwrite_im = reclen
    !
  end function baseio_bwrite_im
!!$!
!!$!H
!!$!H-----------------------------------------------------------------------------
!!$!H integer(FINT) function baseio_bwrite_iu(iun,data)
!!$!H-----------------------------------------------------------------------------
!!$!H
!!$!
!!$  integer(FINT) function baseio_bwrite_iu(iun,data)
!!$    !
!!$    integer(FINT), intent(in) :: iun
!!$    integer(FINT), pointer :: data
!!$    !
!!$    integer(FINT) :: reclen
!!$    integer(FBIGGEST_INT) :: lsize
!!$    !
!!$    !lsize = size(data)
!!$    lsize = 0
!!$    call message(MESWARN,"[baseio_bwrite]: Warning unsupported feature.")
!!$    if (.true.)then
!!$       write(iun) 'i',kind(data),lsize
!!$       write(iun) data
!!$    else
!!$       !reclen = size(data) !baseio_reclen(data)
!!$       reclen = baseio_reclen(data)
!!$       ! we should do something different
!!$       write(iun) 'i',kind(FBIGGEST_INT),lsize
!!$       write(iun) data
!!$    end if
!!$    !
!!$  end function baseio_bwrite_iu
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_bwrite_rn(iun,data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_bwrite_rn(iun,data)
    !
    integer(FINT), intent(in) :: iun
    real(FREAL), intent(in) :: data
    !
    integer(FINT) :: reclen
    !
    reclen = 1
    write(iun) 'i',kind(FBIGGEST_INT),reclen
    write(iun) data
    baseio_bwrite_rn = reclen
    !
  end function baseio_bwrite_rn
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_bwrite_rv(iun,data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_bwrite_rv(iun,data)
    !
    integer(FINT), intent(in) :: iun
    real(FREAL), intent(in), dimension(:) :: data
    !
    integer(FINT) :: reclen
    integer(FBIGGEST_INT) :: lsize
    !
    lsize = size(data)
    reclen = lsize
    if (.true.)then
       write(iun) 'i',kind(data),lsize
       write(iun) data
    else
       !reclen = size(data) !baseio_reclen(data)
       reclen = baseio_reclen(data)
       ! we should do something different
       write(iun) 'i',kind(FBIGGEST_INT),lsize
       write(iun) data
    end if
    baseio_bwrite_rv = reclen
    !
  end function baseio_bwrite_rv
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_bwrite_rm(iun,data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_bwrite_rm(iun,data)
    !
    integer(FINT), intent(in) :: iun
    real(FREAL), intent(in), dimension(:,:) :: data
    !
    integer(FINT) :: reclen
    integer(FBIGGEST_INT) :: lsize
    !
    lsize = size(data,1)*size(data,2)
    reclen = lsize
    if (.true.)then
       write(iun) 'i',kind(data),lsize
       write(iun) data
    else
       !reclen = size(data) !baseio_reclen(data)
       reclen = baseio_reclen(data)
       ! we should do something different
       write(iun) 'i',kind(FBIGGEST_INT),lsize
       write(iun) data
    end if
    baseio_bwrite_rm = reclen
    !
  end function baseio_bwrite_rm
!!$!
!!$!H
!!$!H-----------------------------------------------------------------------------
!!$!H integer(FINT) function baseio_bwrite_ru(iun,data)
!!$!H-----------------------------------------------------------------------------
!!$!H
!!$!
!!$  integer(FINT) function baseio_bwrite_ru(iun,data)
!!$    !
!!$    integer(FINT), intent(in) :: iun
!!$    real(FREAL), pointer :: data
!!$    !
!!$    integer(FINT) :: reclen
!!$    integer(FBIGGEST_INT) :: lsize
!!$    !
!!$    !lsize = size(data)
!!$    lsize = 0
!!$    call message(MESWARN,"[baseio_bwrite]: Warning unsupported feature.")
!!$    if (.true.)then
!!$       write(iun) 'i',kind(data),lsize
!!$       write(iun) data
!!$    else
!!$       !reclen = size(data) !baseio_reclen(data)
!!$       reclen = baseio_reclen(data)
!!$       ! we should do something different
!!$       write(iun) 'i',kind(FBIGGEST_INT),lsize
!!$       write(iun) data
!!$    end if
!!$    !
!!$  end function baseio_bwrite_ru
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_bwrite_cn(iun,data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_bwrite_cn(iun,data)
    !
    integer(FINT), intent(in) :: iun
    character(FLCHARS), intent(in) :: data
    !
    integer(FINT) :: reclen
    !
    reclen = 1
    write(iun) 'i',kind(FBIGGEST_INT),reclen
    write(iun) data
    baseio_bwrite_cn = reclen
    !
  end function baseio_bwrite_cn
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_bwrite_cv(iun,data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_bwrite_cv(iun,data)
    !
    integer(FINT), intent(in) :: iun
    character(FLCHARS), intent(in), dimension(:) :: data
    !
    integer(FINT) :: reclen
    integer(FBIGGEST_INT) :: lsize
    !
    lsize = size(data)
    reclen = lsize
    if (.true.)then
       write(iun) 'i',kind(data),lsize
       write(iun) data
    else
       !reclen = size(data) !baseio_reclen(data)
       reclen = baseio_reclen(data)
       ! we should do something different
       write(iun) 'i',kind(FBIGGEST_INT),lsize
       write(iun) data
    end if
    baseio_bwrite_cv = reclen
    !
  end function baseio_bwrite_cv
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_bwrite_cm(iun,data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_bwrite_cm(iun,data)
    !
    integer(FINT), intent(in) :: iun
    character(FLCHARS), intent(in), dimension(:,:) :: data
    !
    integer(FINT) :: reclen
    integer(FBIGGEST_INT) :: lsize
    !
    lsize = size(data,1)*size(data,2)
    reclen = lsize
    if (.true.)then
       write(iun) 'i',kind(data),lsize
       write(iun) data
    else
       !reclen = size(data) !baseio_reclen(data)
       reclen = baseio_reclen(data)
       ! we should do something different
       write(iun) 'i',kind(FBIGGEST_INT),lsize
       write(iun) data
    end if
    baseio_bwrite_cm = reclen
    !
  end function baseio_bwrite_cm
!!$!
!!$!H
!!$!H-----------------------------------------------------------------------------
!!$!H integer(FINT) function baseio_bwrite_cu(iun,data)
!!$!H-----------------------------------------------------------------------------
!!$!H
!!$!
!!$  integer(FINT) function baseio_bwrite_cu(iun,data)
!!$    !
!!$    integer(FINT), intent(in) :: iun
!!$    character(FLCHARS), pointer :: data
!!$    !
!!$    integer(FINT) :: reclen
!!$    integer(FBIGGEST_INT) :: lsize
!!$    !
!!$    !lsize = size(data)
!!$    lsize = 0
!!$    call message(MESWARN,"[baseio_bwrite]: Warning unsupported feature.")
!!$    if (.true.)then
!!$       write(iun) 'i',kind(data),lsize
!!$       write(iun) data
!!$    else
!!$       !reclen = size(data) !baseio_reclen(data)
!!$       reclen = baseio_reclen(data)
!!$       ! we should do something different
!!$       write(iun) 'i',kind(FBIGGEST_INT),lsize
!!$       write(iun) data
!!$    end if
!!$    !
!!$  end function baseio_bwrite_cu
  !
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function baseio_bread_in(iun,data)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function baseio_bread_in(iun,data)
    !
    integer(FINT), intent(in) :: iun
    integer(FINT), intent(out) :: data
    !
    character(1) :: dtype
    integer(FINT) :: dkind
    integer(FINT) :: reclen
    !
    integer(FINT) :: irc
    !
    reclen = 1
    !write(iun) 'i',kind(FBIGGEST_INT),reclen
    read(iun) dtype,dkind,reclen
    if (dtype.ne.'i') then
       call message(MESERRO,"[baseio_bread]: attempt to read wrong type.")
       bio_error = 1
       baseio_bread_in = -1
       return
    end if
    read(iun,IOSTAT=irc) data
    if (irc.gt.0) then
       call message(MESERRO,"[baseio_bread]: could not read data.")
       bio_error = 1
       baseio_bread_in = -1
    end if
    return
    !
  end function baseio_bread_in
  !
  !
  !



end module baseio
