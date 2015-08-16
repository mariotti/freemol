!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H MODULE molecule  Freemol by F.Mariotti
!H-----------------------------------------------------------------------------
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H $Id: molecule.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
!H-----------------------------------------------------------------------------
!H 
!H This module set the required variables to load molecule data.
!H It prepare data for a single molecule.
!H 
!H [Atoms], [Molecule]
!H 
!H-----------------------------------------------------------------------------
!H 
!
module molecule
  !
  ! modules
  use vartypes
  use messages
  use linetools
  use mathtools
  use chemconst
  !
  implicit none
  private
!
!H
!H This module contains molecule Data.
!H Available Data are:
!H Private
!H   MDATOMSMAXDIM              Default maximun number of allowed atoms
!H   mol_isdist                 logical flag to check on the fly if
!H                              the distance matric has been computed
!H   mol_force_atoms            logical flag to determine if we need
!H                              to force naming convention i.e.
!H                              change user defined name to fit the atomic
!H                              charge.
!H   mol_nat                     Current Number of atoms
!H   mol_natmax                  Current Maximum dimension or
!H                              allocated dimension
!H   mol_numichar                Number of different atoms type by charge
!H   mol_numtypes                Number of different atoms type by name
!H   mol_numids                  Number of different atoms type by IDs
!H   mol_sts                    Current module Status
!H                              Up to date we have only two status:
!H                              0 not initialized
!H                              1 intialized
!H                              But maybe we can have more sttus in the future.
!H   mol_frmt                    Current format: defaulted to 'XYZ'
!H   mol_afrmt                   List of available Available formats
!H   
!H 
!H Public
!H   molecule_iscart            flag to recognize the internal format
!H   molecule_iszmat            flag to recognize the internal format
!H   molecule_isintr            flag to recognize the internal format
!H   molecule_debug             logical debug flag for the module
!H   molecule_error             integer error code
!H                              TODO: Now it return 1 on error but
!H                              in the future will return an integer
!H                              identifying the type of error.
!H                              1 will be reserved for generic error
!H                              and we suppose the check will be of the type:
!H                              if (molecule_error.gt.0) stop 1
!H   molecule_zchar             Vector of integer atomic Z values
!H   molecule_atmID             Vector of Atoms ID (Identification code)
!H   molecule_name              Vector of user defined atomic names
!H   molecule_xyz               Matrix of atomic cartesian cordinates
!H   molecule_atmchar           Vector of atomic nuclear charges
!H   molecule_distm             Distance matrix
!H   molecule_atmseq            Sequential number given in the input
!H   TODO:
!H   mol_zmat_int               Matrix of integers for Zmatrix
!H   mol_zmat_real              Matrix of reals for Zmatrix
!H   mol_cint_int               Matrix of integers for Internal coordinates
!H   mol_cint_real              Matrix of reals for Internal coordinates
!H
!H   molecule_numsym            Maximum number of different atoms by Z (= 105)
!H   mol_symtbl                 Atoms Z to symbols table
!H   
!H   
!H 
!H Available procedures are:
!H Private (sub or fun)
!H   sub pmol_readxyz           local routine to read effectively
!H                              the cartesian coordinates.
!H 
!H Public (fun or sub)
!H   fun molecule_init          initialize the module
!H   fun molecule_terminate     terminate the module
!H   fun molecule_allocate      allocate all the workspace
!H   fun molecule_deallocate    deallocate all the workspace
!H   fun molecule_setnat        set the current number of atoms
!H   fun molecule_getnat        get the current number of atoms
!H   fun molecule_getnumtypes   get number of different atom types by Z
!H   fun molecule_getnumids     get the number of different atoms by ID
!H   fun molecule_setnumtypes   set the number of different atoms type by Z
!H                              WARNING: This function will be removed
!H   fun molecule_setnumids     set the number of different atoms type by ID
!H                              WARNING: This function will be removed
!H   fun molecule_getdim        Return the maximum numer of records storable
!H   fun molecule_isinit        Query if the module is initialized
!H   fun molecule_isallocated   Query if the records are allocated
!H   fun molecule_setpar        Modify some default parameters (on fly)
!H   fun molecule_compdist      compute the distance matrix
!H   sub molecule_print         Print ascii version of the section
!H   sub molecule_readascii     Read ascii version of the section
!H   sub molecule_save          Save the data in a temporary device
!H   sub molecule_store         Save the data in static device
!H   sub molecule_restore       Restore the data from a temporary/static device
!H   sub molecule_convert       Perform some data conversion
!H   sub molecule_getfrmt       Get the current format
!H   sub molecule_check         Perform an integrity check
!H   sub molecule_sort          Sort atoms by Z values
!H   sub molecule_sortby        Sort atoms by integer indexed list
!H   
!H-----------------------------------------------------------------------------
!H   
!H FOLLOW COMMENTED CODE AND PROCEDURE DESCRIPTION
!H   
!H 

  ! Default Values
  !----------------
  integer(FINT), private, parameter :: MDATOMSMAXDIM = 500

  ! Public variables
  !------------------
  logical, public, save :: molecule_debug = .false.
  integer(FINT), public, save :: molecule_error = 0
  ! Atom Z
  integer(FINT), public, save, allocatable, dimension(:) :: molecule_zchar
  ! Atom ID
  integer(FINT), public, save, allocatable, dimension(:) :: molecule_atmID
  ! Atom Names
  character(FLCHARS), public, save, allocatable, dimension(:) :: molecule_name
  ! Atom Sequential number from the input
  integer(FINT), public, save, allocatable, dimension(:) :: molecule_atmseq
  ! Atom Coordinates
  real(FREAL), public, save, allocatable, dimension(:,:) :: molecule_xyz
  ! Atom Effective Nuclear Charge
  real(FREAL), public, save, allocatable, dimension(:) :: molecule_atmchar
  ! Distance matrix
  real(FREAL), public, save, allocatable, dimension(:,:) :: molecule_distm
  !
  logical, public, save :: molecule_iscart
  logical, public, save :: molecule_iszmat
  logical, public, save :: molecule_isintr

  ! Private General use Vars
  !--------------------------
  ! Is the dist computed?
  logical, private, save :: mol_isdist = .false.
  !
  ! Now used to force a change in name based on charge
  logical, private, save :: mol_force_atoms = .false.
  !
  ! check if we use angstrom or bohr
  logical, private, save :: molecule_angs = .false.
  !
  ! check if we use an alternative unit
  logical, private, save :: molecule_xunit = .false.
  real(FREAL), private, save :: molecule_vxunit = 0.0_FREAL
  !
  ! Number of atoms
  integer(FINT), private, save :: mol_nat = 0
  !
  ! Max dimension or allocated dimension
  integer(FINT), private, save :: mol_natmax = 0
  !
  ! Number of different atoms type by charge
  integer(FINT), private, save :: mol_numichar = 0
  !
  ! Number of different atoms type by name
  integer(FINT), private, save :: mol_numtypes = 0
  !
  ! Number of different atoms type by IDs
  integer(FINT), private, save :: mol_numids = 0
  !
  ! Current Status
  integer(FINT), private, save :: mol_sts = 0
  !
  ! Current format
  character(FLCHARS), private, save :: mol_frmt = 'XYZ'
  !
  ! Available formats
  character(FLCHARS), private, save, dimension(2) :: mol_afrmt
  !
  ! Functions
  !-----------
  private pmol_readxyz
  !
  public molecule_init
  public molecule_terminate
  public molecule_allocate
  public molecule_deallocate
  public molecule_setnat
  public molecule_getnat
  public molecule_getnumtypes
  public molecule_setnumtypes
  public molecule_getnumids
  public molecule_setnumids
  public molecule_getdim
  public molecule_isinit
  public molecule_isallocated
  public molecule_setpar
  public molecule_compdist
  public molecule_rotxyz
  !
  ! Subroutines
  !-------------
  public molecule_print
  public molecule_readascii
  public molecule_read
  public molecule_save
  public molecule_store
  public molecule_restore
  public molecule_convert
  public molecule_getfrmt
  public molecule_check
  public molecule_sort
  public molecule_sortby

  ! Possible next functions?
  !--------------------------
  ! _isvalued  => to check if there are stored values
  !
  !H
  !H We store an internal table for transaltion from
  !H atomic Z value and lower case atomic symbols.
  !H Shall we in the future define a symbols module?
  !H We made it public here anyway!
  !H
  !
  ! data
  !------
  integer(FINT), public, parameter :: molecule_numsym = 105
  character(2), public, dimension(molecule_numsym):: mol_symtbl
  data mol_symtbl /'h ','he','li','be','b ','c ','n ','o ','f ',&
       &'ne','na','mg','al','si','p ','s ','cl','ar','k ',&
       &'ca','sc','ti','v ','cr','mn','fe','co','ni','cu',&
       &'zn','ga','ge','as','se','br','kr','rb','sr','y ',&
       &'zr','nb','mo','tc','ru','rh','pd','ag','cd','in',&
       &'sn','sb','te','i ','xe','cs','ba','la','ce','pr',&
       &'nd','pm','sm','eu','gd','tb','dy','ho','er','tm',&
       &'yb','lu','hf','ta','w ','re','os','ir','pt','au',&
       &'hg','tl','pb','bi','po','at','rn','fr','ra','ac',&
       &'th','pa','u ','np','pu','am','cm','bk','cf','es',&
       &'fm','md','no','lr','rf','ha'/
!
!
contains
!
!H
!H-----------------------------------------------------------------------------
!H Function molecule_init
!H-----------------------------------------------------------------------------
!H
  integer(FINT) function molecule_init(debugflag)
    !
    logical, intent(in), optional :: debugflag
    !
    integer(FINT) :: irc
    !
    ! Set the debug flag
    !
    ! we will produce at least no errors
    molecule_init = 0
    molecule_error = 0
    !
    if(present(debugflag)) then
       molecule_debug = debugflag
    else
       molecule_debug = .false.
    end if
    !
    ! We call molecule_setpar in order to set all
    ! the dafault parameters.
    irc = molecule_setpar(usename=.true.,angs=.false.,xunit=.false.)
    if(irc.ne.0) then
       molecule_error = 1
       molecule_init = -1
       return
    end if
    !
    ! change of the following parameters
    ! is not allowed at run time, so we set them here.
    ! The default is to have cartesian coordinates.
    molecule_iscart = .true.
    molecule_iszmat = .false.
    molecule_isintr = .false.
    !
    ! We switch to status 1: the module has been initialezed
    mol_sts = 1
  end function molecule_init
!
!H
!H-----------------------------------------------------------------------------
!H Function molecule_terminate
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function molecule_terminate()
    !
    ! We just call the deallocating function
    ! and we set some variables
    !
    integer(FINT) :: irc
    !
    ! Reset error flag
    molecule_error = 0
    !
    irc = molecule_deallocate()
    !
    ! The function molecule_deallocate will set
    ! for us the number of atoms to zero.
    !
    mol_sts = 0 !UnInitialized
    molecule_terminate = 0
    !
  end function molecule_terminate
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function molecule_allocate(idmn,forceit)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function molecule_allocate(idmn,forceit)
!
!H
!H This function allocate all the vectors in order to store
!H molecule data.
!H WANRING! it does set mol_natmax but NOT mol_nat.
!H
!
    !
    ! We try to allocate all the data variables to the
    ! maximum/requestet value. We even reset the values!!
    !
    ! The dimension of the module: nomber of atoms
    integer(FINT), optional, intent(in) :: idmn
    ! Shal we force allocation (if it is already allocated?)
    logical, optional, intent(in) :: forceit
    !
    integer(FINT) :: ival
    integer(FINT) :: irc
    !
    !
    ! Resetdefaults
    molecule_error = 0
    molecule_allocate = 0
    !Set de default dimension if not given in the parameter
    ival = MDATOMSMAXDIM

    if(present(idmn)) then
       if(idmn.le.0) then
          ! We cannot allocate .lt.0
          molecule_allocate = -1
          molecule_error = 1
          !
          ! We call deallocate anyway to prduce runtime errors
          ! and to try an interpretation of this call.
          irc = molecule_deallocate()
          return
       else
          ! We set the dimension from the call
          ival = idmn
       endif
    endif
    !
    if(present(forceit)) then
       if(forceit) then
          ! We asked to force allocation so we deallocate first.
          irc = molecule_deallocate()
       end if
       !
    end if
    !
    !STEP: Allocate characters
    if(allocated(molecule_name)) then
       molecule_allocate = -1
       molecule_error = 1
       call message(MESDEBG,'[molecule_allocate] Already allocated.')
       return
    end if
    allocate(molecule_name(ival),STAT=molecule_error)
    !
    !on no error we proceed to the others values
    if(molecule_error.eq.0) then
       !
       molecule_name(1:ival) = ' ' !Reset
       !STEP: Allocate Integers
       if   (allocated(molecule_atmID).or.&
            &allocated(molecule_zchar).or.&
            &allocated(molecule_atmseq)) then
          !
          molecule_allocate = -1
          molecule_error = 1
          call message(MESDEBG,'[molecule_allocate] Already allocated.')
          return
          !
       end if
       !
       allocate(molecule_atmID(ival),&
               &molecule_zchar(ival),&
               &molecule_atmseq(ival),STAT=molecule_error)
       !
       if(molecule_error.eq.0) then
          molecule_atmID(1:ival) = 0_FINT !Reset Dafault
          molecule_zchar(1:ival) = 0_FINT !Reset Dafault
          molecule_atmseq(1:ival) = 0_FINT
          !STEP: allocate reals
          if(allocated(molecule_xyz).or.allocated(molecule_atmchar)) then
             molecule_allocate = -1
             molecule_error = 1
             call message(MESDEBG,'[molecule_allocate] Already allocated.')
             return
          end if
          !
          allocate(molecule_xyz(3,ival),&
                  &molecule_atmchar(ival),STAT=molecule_error)
          !
          if(molecule_error.eq.0) then
             molecule_atmchar(1:ival) = 0.0_FREAL
             molecule_xyz(1:3,1:ival) = 0.0_FREAL
             !set the maximun storable number of atoms
             mol_natmax = ival
          else
             ! we are not able to allocate reals
             ! so we download previous allocated stuff
             deallocate(molecule_name)
             deallocate(molecule_atmID,molecule_zchar,molecule_atmseq)
             molecule_error = 1
             molecule_allocate = -1
             return
          endif
       else
          ! we are not able to allocate molecule_atmID and molecule_zchar
          ! so we download previous allocated stuff
          deallocate(molecule_name)
          molecule_error = 1
          molecule_allocate = -1
       endif
    else
       !error on alloc of molecule_name
       molecule_error = 1
       molecule_allocate = -1
       return
    endif
    !
    ! We did it!
    !
    return
    !
  end function molecule_allocate
!H
!H-----------------------------------------------------------------------------
!H Function molecule_deallocate
!H-----------------------------------------------------------------------------
!H
  integer(FINT) function molecule_deallocate()
    !
    if(allocated(molecule_atmID)) deallocate(molecule_atmID)
    if(allocated(molecule_zchar)) deallocate(molecule_zchar)
    if(allocated(molecule_atmseq)) deallocate(molecule_atmseq)
    if(allocated(molecule_name)) deallocate(molecule_name)
    if(allocated(molecule_atmchar)) deallocate(molecule_atmchar)
    if(allocated(molecule_xyz)) deallocate(molecule_xyz)
    !
    !at this point we are quite sure there are no atoms!!
    mol_nat = 0
    mol_natmax = 0
    !
    !We suppose allway no errors
    molecule_deallocate = 0
    !
  end function molecule_deallocate
!
!H
!H-----------------------------------------------------------------------------
!H Function molecule_getnat
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function molecule_getnat()
!
!H
!H  Return the current number of atoms
!H
!
    if(mol_sts.eq.0) then
       molecule_getnat = -1
    else
       molecule_getnat = mol_nat
    endif
  end function molecule_getnat
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function molecule_setnat(inat)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function molecule_setnat(inat)
!
!H
!H  Set the current number of atoms
!H  Use this function carefully: no way to get the old number of atoms
!H
!
    integer(FINT), intent(in) :: inat
    !
    ! reset to no errors
    molecule_error = 0
    !
    if(mol_sts.eq.0) then
       molecule_setnat = -1
       mol_nat = 0
       molecule_error = 1
    else
       molecule_setnat = 0
       ! TODO: Add a check to the allocated dimension...
       mol_nat = inat
    endif
    !
  end function molecule_setnat
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function molecule_getdim()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function molecule_getdim()
!
!H
!H  Return the maximun storable number of atoms
!H  Return -1 if not allocated! NOT IMPLEMENTED
!H
!H TODO:
!H Plase check this function: we trust mol_natmax
!H but we dont check if it is correct!!
!H
!
    molecule_getdim = mol_natmax
    if(mol_natmax.lt.0) then
       molecule_getdim = -1
    end if
    !
  end function molecule_getdim
!
!H
!H-----------------------------------------------------------------------------
!H logical function molecule_isinit()
!H-----------------------------------------------------------------------------
!H
!
  logical function molecule_isinit()
!
!H
!H  Check if the module has been initialized
!H We just query for internal variables....
!H
!
    if(mol_sts.eq.0) then
       molecule_isinit = .false.
    else
       molecule_isinit = .true.
    endif
  end function molecule_isinit
!
!H
!H-----------------------------------------------------------------------------
!H logical function molecule_isallocated()
!H-----------------------------------------------------------------------------
!H

  logical function molecule_isallocated()
!
!H
!H  Check if all the vectore are allocated
!H We just step over allocated quesries ....
!H
!
    molecule_isallocated = .true.
    molecule_isallocated = allocated(molecule_atmID).and.molecule_isallocated
    molecule_isallocated = allocated(molecule_zchar).and.molecule_isallocated
    molecule_isallocated = allocated(molecule_atmseq).and.molecule_isallocated
    molecule_isallocated = allocated(molecule_name).and.molecule_isallocated
    molecule_isallocated = allocated(molecule_atmchar).and.molecule_isallocated
    molecule_isallocated = allocated(molecule_xyz).and.molecule_isallocated
!
  end function molecule_isallocated
!
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function molecule_setpar(usename,angs,xunit,vxunit)
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function molecule_setpar(usename,angs,xunit,vxunit)
!
!H
!H This function set up some parameters for the molecule module
!H Here we trust the default value so we do not change
!H any value if the corresponding parameter is not given.
!H REMEMBER this function does not change the data value.
!H in order to transform from one unit to an other
!H use the function molecule_convert.
!H Use it always with optional parameters by name
!H i.e.: molecule_setpar(usename=.false.)
!H
!H
!H INPUT
!H 
!H      USENAME      Usename is used in atoms_check in order
!H                   to generate new names if the charge and name
!H                   does not match. If .false. (Default) only
!H                   a warning is generated. This set the
!H                   mol_force_atoms variable.
!H      ANGS         Determine if we are using angstroms
!H                   Within the module bohrs are the default units.
!H                   This set the molecule_angs variable
!H      XUNIT        We propose an alternative measure unit.
!H      VXUNIT       This contain the corresponding conversion
!H                   factor defined by NEWUNIT = VXUNIT * bohrs
!H OUTPUT
!H     Return an error code (-1) in case of errors.
!H
!
    logical, intent(in), optional :: usename
    logical, intent(in), optional :: angs
    logical, intent(in), optional :: xunit
    real(FREAL), intent(in), optional :: vxunit
    !
    ! As usual we reset the error code
    molecule_error = 0
    ! Check usename
    if(present(usename)) then
       mol_force_atoms = usename
    end if
    ! Check angs
    if(present(angs)) then
       molecule_angs = angs
    endif
    ! Check xunit
    if(present(xunit)) then
       if(xunit) then
          if(present(vxunit)) then
             molecule_xunit = .true.
             molecule_vxunit = vxunit
          else
             ! This is an internal error
             ! we cannot set up an alternative unit measure
             ! without defining it
             molecule_xunit = .false.
             molecule_setpar = -1
             molecule_error = 1
             return
          endif
       end if
    endif
    !
    molecule_setpar = 0
    !
  end function molecule_setpar
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function molecule_getnumtypes()
!H-----------------------------------------------------------------------------
!H
  integer(FINT) function molecule_getnumtypes()
!
!H
!H  Return the number of different atom types
!H
!
    if(mol_numtypes.gt.0) then
       molecule_getnumtypes = mol_numtypes
    endif

  end function molecule_getnumtypes
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function molecule_getnumids()
!H-----------------------------------------------------------------------------
!H
  integer(FINT) function molecule_getnumids()
!
!H
!H  Return the number of different atom IDs !NOT IMPLEMENTED IN THE CODE
!H
!
    if(mol_numtypes.gt.0) then
       molecule_getnumids = mol_numids
    endif

  end function molecule_getnumids
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function molecule_setnumtypes(ntypes)
!H-----------------------------------------------------------------------------
!H
  integer(FINT) function molecule_setnumtypes(ntypes)
!
!H
!H  Set the number of different atoms types !! TO BE REMOVED
!H
!H We will not allowed to hack the code ... we will kill this function
!H
!
    integer(FINT), intent(in) :: ntypes

    molecule_setnumtypes = 0
    mol_numtypes = ntypes

  end function molecule_setnumtypes
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function molecule_setnumids(nids)
!H-----------------------------------------------------------------------------
!H
  integer(FINT) function molecule_setnumids(nids)
!
!H
!H  Set the number of different atoms IDs !! TO BE REMOVED
!H
!H We will not allowed to hack the code ... we will kill this function
!H
!
    integer(FINT), intent(in) :: nids

    molecule_setnumids = 0
    mol_numids = nids

  end function molecule_setnumids
!H
!H-----------------------------------------------------------------------------
!H integer(FINT) function molecule_compdist(diag)
!H-----------------------------------------------------------------------------
!H
  integer(FINT) function molecule_compdist(diag)
!
!H
!H  Compute the distance matrix
!H  If the optional logic parameted 'diag' is true
!H  the result is stored in a vector else in a whole
!H  complete matrix. VECTOR FORM NOT YET IMPLEMENTED!
!H
!
    logical, optional :: diag
    integer(FINT) :: i,j
    integer(FINT) :: irc

    if(present(diag)) then
       if(.not.diag) then
          allocate(molecule_distm(mol_nat,mol_nat),STAT=irc)
          if(irc.ne.0) then
             molecule_compdist = -1
             return
          endif
       else
          molecule_compdist = -1 !Not available linear distance matrix
       endif
    else
       allocate(molecule_distm(mol_nat,mol_nat),STAT=irc)
       if(irc.ne.0) then
          molecule_compdist = -1
          return
       endif
    endif
    !
    ! WARNING
    ! BE CAREFUL: this code assume that you are
    ! not allowed to get the vector form of distance matrix ...
    !
    !DIAG DISTANCE MATRIX
!
!   Compute the distance matrix
!-----------------------------------------------------------------------
!
    do i=1,mol_nat
       do j=i+1,mol_nat
          molecule_distm(i,j) = sqrt(&
               (molecule_xyz(1,i)-molecule_xyz(1,j))**2 +&
               (molecule_xyz(2,i)-molecule_xyz(2,j))**2 +&
               (molecule_xyz(3,i)-molecule_xyz(3,j))**2 )
       enddo
!------diagonal side
       molecule_distm(i,i) = 0.0_FREAL
    enddo
!---copy second side
    do i=1,mol_nat
       do j=i+1,mol_nat
          molecule_distm(j,i) = molecule_distm(i,j)
       enddo
    enddo
    !
    mol_isdist = .true.
    !
  end function molecule_compdist
!H
!H-----------------------------------------------------------------------------
!H Function 
!H-----------------------------------------------------------------------------
!H

!H
!H-----------------------------------------------------------------------------
!H Function 
!H-----------------------------------------------------------------------------
!H


!H
!H-----------------------------------------------------------------------------
!H subroutine molecule_print(iun)
!H-----------------------------------------------------------------------------
!H
!
  subroutine molecule_print(iun, frmt)
    !
    integer(FINT), intent(in) :: iun
    character(*), optional, intent(in) :: frmt
    !
    character(FLCHARS) :: lfrmt
    integer(FINT) :: i
    !
    lfrmt = adjustl(mol_frmt)
    if(present(frmt)) then
       i=len(frmt)
       if(i.gt.FLCHARS) i=FLCHARS
       lfrmt=adjustl(frmt(1:i))
    endif
!
!H
!H  WARNING: Warning there is a limit on max atoms due to format
!H
!
    if(index(lfrmt,'pxyz').gt.0) then
       if (lfrmt(i:i+4).eq.'pxyz') then
          write(iun,'(1X,I5)') mol_nat
          write(iun,'(1X,"CSM pxyz print out")')
          do i=1,mol_nat
             !          write(iun,'(2X,A12,2X,3(1X,D16.8))')&
             !               &molecule_name(i),molecule_xyz(1:3,i)
             write(iun,'(2X,A12,2X,3(1X,F12.8))')&
                  &molecule_name(i),molecule_xyz(1:3,i)
          end do
       end if
       !
    else if(index(lfrmt,'mldfromau').gt.0) then
       !
       write(iun,'(A16,I4,1X,A7,A)')&
            &'[molecule] &
            &Nrec=',mol_nat,'Format=',mol_frmt(1:len_trim(mol_frmt))
       !
       do i=1,mol_nat
          !write(iun,'(2X,A12,2X,I3,2X,F10.4,1X,3(1X,D16.8))')&
          !     &molecule_name(i),molecule_zchar(i),&
          !     &molecule_atmchar(i),molecule_xyz(1:3,i)
          !write(iun,'(2X,A12,2X,I3,2X,3(1X,D16.8))')&
          !     &molecule_name(i),molecule_zchar(i),molecule_xyz(1:3,i)
          write(iun,'(2X,A12,1X,2(2X,I3),1X,3(1X,D16.8))')&
               &molecule_name(i),i,&
               &molecule_zchar(i),molecule_xyz(1:3,i)/cc_au2ang
       end do
       !
    else
!
!This is the MOLDEN format compatibility .. old code ..
!
!    write(iun,'(A13,I4,1X,A5,A)')&
!         &'[Atoms] Nrec=',mol_nat,'FRMT=',mol_frmt(1:len_trim(mol_frmt))
!
       write(iun,'(A16,I4,1X,A7,A)')&
            &'[molecule] &
            &Nrec=',mol_nat,'Format=',mol_frmt(1:len_trim(mol_frmt))

       if(molecule_debug) then
          call message(MESDEBG,&
               &'[molecule_print] Remember Format string not supported')
       endif
       
       do i=1,mol_nat
          !write(iun,'(2X,A12,2X,I3,2X,F10.4,1X,3(1X,D16.8))')&
          !     &molecule_name(i),molecule_zchar(i),&
          !     &molecule_atmchar(i),molecule_xyz(1:3,i)
          !write(iun,'(2X,A12,2X,I3,2X,3(1X,D16.8))')&
          !     &molecule_name(i),molecule_zchar(i),molecule_xyz(1:3,i)
          write(iun,'(2X,A12,2X,I3,2X,F10.4,1X,3(1X,D16.8))')&
               &molecule_name(i),i,&
               &molecule_atmchar(i),molecule_xyz(1:3,i)
       end do

    end if

  end subroutine molecule_print
!H
!H-----------------------------------------------------------------------------
!H Subroutine molecule_readascii
!H-----------------------------------------------------------------------------
!H
!
  subroutine molecule_readascii(iun,frmt)
    !
    ! We use some constants
    use chemconst
!
!H
!H    This routine reads atom card from a Freemol Format file.
!H
!H    input
!H          IUN          Unit file for read.
!H          FRMT         Will Define a "force" key for the format.
!H                       NOT YET IMPLEMENTED. (optional)
!H
!H    Otput
!H          Try to fill the molecule module with all
!H          the available data.
!H          Call the routine molecule_check to
!H          perform a check of the readed data.
!H
!H NOTES about the format (see also inline code).
!H
!H The atom alias string must be the first string in the
!H record if present.
!H
!
    !
    integer(FINT) :: iun
    character(FLCHARS), optional :: frmt
    !
    ! Local variables
    !
    ! wether or not perform a conversion and its value
    logical :: doxcon
    real(FREAL) :: vconv
    !
    ! do we have a given format?
    logical :: havfrm
    character(FLCHARS) :: lfrmt
    !
    ! INITIALIZATION
    !
    ! reset the error flag
    molecule_error = 0
    !
    ! we have not yet read anything: so reset.
    ! The conversion is set: We read Angs and need bohrs
    doxcon = .true.
    vconv = cc_au2ang
    !
    ! we will check for the format later
    havfrm = .false.
    lfrmt = ''
    !
    if (.not.molecule_angs) then
       if (molecule_xunit) then
          vconv = molecule_vxunit
       else
          vconv = 1.0_FREAL ! Actually do not perform the conversion
       endif
    endif
    !
    if(present(frmt)) then
       havfrm = .true.
       lfrmt = frmt
       call message(MESWARN,'RDMOL: frmt not yet implemented.')
       molecule_error = 1
       return
    endif
    !
    ! In this routine we do not check for the main format
    if(molecule_iscart) then
       call pmol_readxyz(iun,doxcon,vconv,lfrmt)
    else if (molecule_iszmat) then
       call message(MESERRO,'RDMOL: ZMATRIX Not Yet Implemented.')
       molecule_error = 1
       return
    else if (molecule_isintr) then
       call message(MESERRO,'RDMOL: INTERNAL Not Yet Implemented.')
       molecule_error = 1
       return
    else
       call message(MESERRO,'RDMOL: Unknown format!!.')
       molecule_error = 1
       return
    end if
    !
    !
    !
  end subroutine molecule_readascii
!
!H
!H-----------------------------------------------------------------------------
!H Subroutine pmol_readxyz
!H-----------------------------------------------------------------------------
!H
!
  subroutine pmol_readxyz(iun,doxcon,vcon,frmt)
    !
    use linetools
    !
    integer(FINT), intent(in) :: iun
    logical, intent(in) :: doxcon
    real(FREAL), intent(in) :: vcon
    character(FLCHARS), intent(in), optional :: frmt
    !
    integer(FINT) :: irc
    logical ::done
    character(FLCHARS) :: line
    integer(FINT) :: iline
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
    integer(FINT) :: maxlines
    !
    ! Init
    !------
    done = .false.
    molecule_error = 0
    iline = 0
    !
    ! we save the maximum number of allocated data
    maxlines = mol_natmax
    !
    do while(.not.done)
       ! we use line_getline with mode 3
       ! because we want to strip out comment lines as well.
       irc = line_getline(iun,line,3)
       if(irc.ge.0) then
          ! are we done?
          ! TODO: put on getline the check about new section ...
          if(.not.(line(1:1).eq.'[')) then
             call line_read(line,iwrds,rwrds,cwrds,niwrd,nrwrd,ncwrd,lrform)
             !
             ! We require at least the nuclei charge or the atom label
             !--------------------------------------------------------
             if (niwrd.eq.0.and.ncwrd.eq.0) then
                call message(MESERRO,"[rdxyz]: Cannot identify the atom.")
                call message_value(MESERRO,"[rdxyz]: line = ",iline)
                molecule_error = 1
                return
             endif
             !
             ! We require the coordinates
             !---------------------------
             if (nrwrd.lt.3) then
                call message(MESERRO,"[rdxyz]: Cannot find coordinates.")
                call message_value(MESERRO,"[rdxyz] At line= ",iline)
                molecule_error = 1
                return
             endif
             iline = iline + 1
             if(iline.gt.maxlines) then
                done = .true.
                molecule_error = 1
                call message(MESERRO,"[rdxyz]: Too many lines of data.")
                call message_value(MESERRO,"[rdxyz]: max lines = ",maxlines)
                ! before return we set the number of readed atoms
                !------------------------------------------------
                mol_nat = maxlines
                return
             endif
             ! The atom name is the first and
             ! we are sure we have at least one cword
             !---------------------------------------
             molecule_name(iline) = cwrds(1)
             ! We read the charge if given and the
             ! cartesian coordinates CXYZ or XYZ
             !----------------------------------
             if (nrwrd.gt.3) then
                molecule_atmchar(iline) = rwrds(1)
                molecule_xyz(1,iline) = rwrds(2)
                molecule_xyz(2,iline) = rwrds(3)
                molecule_xyz(3,iline) = rwrds(4)
             else
                molecule_xyz(1,iline) = rwrds(1)
                molecule_xyz(2,iline) = rwrds(2)
                molecule_xyz(3,iline) = rwrds(3)
             endif
             ! Shall we have to perform the unit conversion
             if(doxcon) then
                molecule_xyz(1,iline) = molecule_xyz(1,iline)/vcon
                molecule_xyz(2,iline) = molecule_xyz(2,iline)/vcon
                molecule_xyz(3,iline) = molecule_xyz(3,iline)/vcon
             end if
             ! If we have an integer and it comes before the name
             ! it is a sequential number.
             !---------------------------
             if ((lrform(1:1).eq.'I').and.(ncwrd.gt.0)) then
                molecule_atmseq(iline) = iwrds(1)
                ! if we have more than 1 integer and there is not
                ! the atom name then the first number
                ! is the sequential number.
                !--------------------------
             !else if ((niwrd.gt.1).and.(ncwrd.eq.0)) then
                ! We change it in:
                ! Two integers means that the second is the charge
                ! and the first is the sequential number.
                !-------------------------------------------------
             else if ((niwrd.gt.1)) then
                molecule_atmseq(iline) = iwrds(1)
                molecule_zchar(iline) = iwrds(2)
                ! Every other case we suppose to have the
                ! integer charge identifying the atom
                !------------------------------------
             else
                molecule_zchar(iline) = iwrds(1)
             end if
             ! We are done: so set the number of atoms
             !----------------------------------------
             mol_nat = iline
          else
             done = .true.
          end if
       else
          ! We do not matter: we are done if we get errors on reading a line
          ! TODO: maybe we can check which kind of error we get....
          !--------------------------------------------------------
          done = .true.
       end if
    end do
    !
  end subroutine pmol_readxyz
!
!H
!H-----------------------------------------------------------------------------
!H Subroutine molecule_save
!H-----------------------------------------------------------------------------
!H
!
  subroutine molecule_save

  end subroutine molecule_save

!H
!H-----------------------------------------------------------------------------
!H Subroutine molecule_store
!H-----------------------------------------------------------------------------
!H

  subroutine molecule_store

  end subroutine molecule_store

!H
!H-----------------------------------------------------------------------------
!H Subroutine molecule_restore
!H-----------------------------------------------------------------------------
!H

  subroutine molecule_restore

  end subroutine molecule_restore

!H
!H-----------------------------------------------------------------------------
!H Subroutine molecule_convert
!H-----------------------------------------------------------------------------
!H

  subroutine molecule_convert(convf)
    !
    real(FREAL), intent(in) :: convf
    !
    integer(FINT) :: lnat
    !
    lnat = molecule_getnat()
    !
    molecule_xyz(1:3,1:lnat) = molecule_xyz(1:3,1:lnat) * convf
    !
  end subroutine molecule_convert

!H
!H-----------------------------------------------------------------------------
!H Subroutine molecule_getfrmt
!H-----------------------------------------------------------------------------
!H

  subroutine molecule_getfrmt

  end subroutine molecule_getfrmt


!H
!H-----------------------------------------------------------------------------
!H Subroutine molecule_check
!H-----------------------------------------------------------------------------
!H
  
  subroutine molecule_check

    !CHECK IT
    use strtools

  !
  !
  !     Check Names and Atom types
  !
  !
  !     Local
  !
  logical :: notfound
  logical :: notdone
  character(FLCHARS) :: cdummy
  integer(FINT) :: i
  integer(FINT) :: j
  integer(FINT) :: ic
  integer(FINT) :: irc
  !
  integer(FINT) :: ityp
  !
  integer(FINT) :: idx1
  integer(FINT) :: idx2
  integer(FINT) :: idx3
  !
  !     Check Atom Types
  !
  notfound = .true.
  notdone = .true.
  ityp = 1
  molecule_atmID (1) = 1
  !
  ! Translate each atom name in lower case (I like it!)
  !----------------------------------------------------
  call str_lowcase(molecule_name)
  !
  ! Build Atoms type based on atom name
  do i = 2, mol_nat
     if(molecule_name(i)(1:1).eq.' ') then !not yet named
        if(molecule_zchar(i).ne.0) then
           molecule_name(i) = mol_symtbl(molecule_zchar(i))
        endif
     endif
     do j = 1, i - 1
        if(molecule_name(j)(1:1).eq.' ') then !not yet named
           if(molecule_zchar(j).ne.0) then
              molecule_name(j) = mol_symtbl(molecule_zchar(j))
           endif
        endif
        irc = str_strcmp(molecule_name(i),molecule_name(j))
        if (irc.gt.0) then
           ! Atom Type exist. Already stored
           if(molecule_zchar(i).eq.molecule_zchar(j)) then
              molecule_atmID(i) = molecule_atmID(j)
              notfound = .false.
           endif
        endif
     enddo
     if (notfound) then
        ! New Atom Type
        ityp = ityp + 1
        molecule_atmID(i) = ityp
     endif
     notfound = .true.
  enddo
  !
  do i=1,mol_nat
     !
     ic = molecule_zchar(i) ! get Z and check if not out of bounds
     if (ic.gt.molecule_numsym) then
        call message(MESERRO,'SYM: Nuclear Charge Bigger than Max!')
        call message(MESKILL,'SYM: Kill request.')
     endif
     !
     ! Find molecule_name AtomName
     idx1 = index(molecule_name(i),'.')
     idx2 = index(molecule_name(i),'-')
     idx3 = index(molecule_name(i),'_')
     if(idx1.eq.0) then
        if(idx2.eq.0) then
           idx1 = idx3
           if(idx3.eq.0) then
              idx1 = len_trim(molecule_name(i)) + 1
           endif
        else
           if(idx3.eq.0) then
              idx1 = idx2
           else
              idx1 = min(idx2,idx3)
           endif
        endif
     else
        if(idx2.eq.0) then
           if(idx3.ne.0) then
              idx1 = min(idx1,idx3)
           endif
        else
           if(idx3.eq.0) then
              idx1 = min(idx1,idx2)
           else
              idx1 = min(idx1,idx2,idx3)
           endif
        endif
     endif
     idx1 = idx1 - 1
     if (ic.gt.0) then ! I have the nuclear charge!!!
        !Compare up to idx1
        irc = str_strcmp(molecule_name(i),mol_symtbl(ic),idx1)
        if (irc.eq.0) then
           call message(MESWARN,'SYM: Atom charge and atom name differ.')
           if(mol_force_atoms) then
              call message(MESWARN,'SYM: use charge definition.')
              j = len_trim(mol_symtbl(ic))
              cdummy = molecule_name(i)
              molecule_name(i) = mol_symtbl(ic)(1:j)//molecule_name(i)(idx1+1:)
              cdummy = 'SYM: '//cdummy(1:len_trim(cdummy))//&
                   &' Changed in '//molecule_name(i)
              call message(MESWARN,cdummy)
           endif
        endif
     else  ! No charge ..oopss..
        notfound = .true.
        !
        j = 0
        do while(notdone)
           j = j + 1
           if(j.le.molecule_numsym) then
              irc = str_strcmp(molecule_name(i),mol_symtbl(j),idx1)
              if (irc.gt.0) then
                 molecule_zchar(i) = j
                 notfound = .false.
                 notdone = .false.
              endif
           else
              notdone = .false.
              !write (cdummy,'(1X,12A,1X,I4)')&
              !     &'SYM: Unable to find Atom Number:', i
              !call message(MESERRO, cdummy)
              call message_value(MESERRO,&
                   &'[molecule_check] Unable to find Atom Number:',i)
              call message(MESERRO,'SYM: AtomType Card not yet supported.')
              call message(MESKILL,'SYM: Kill request.')
           endif
        enddo
     endif
  enddo
  !
  !
  !
  return
    

  end subroutine molecule_check
!
!H
!H-----------------------------------------------------------------------------
!H Subroutine molecule_sort
!H-----------------------------------------------------------------------------
!H
!
  recursive subroutine molecule_sort(order)
    !
    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C.
    ! (1990) "Programmer's Guide to Fortran 90",
    ! McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated
    ! integer array which gives
    ! the positions of the elements in the original order.
    ! 
    ! Adapted to sort Atoms in FREEMOL by F.Mariotti.
    !
    integer(FINT), dimension (:), intent(OUT) :: order
    !
    ! Local variable
    integer(FINT) :: i
    !
    ! For md_atoms .. the distance matric must be compuoted again
    mol_isdist = .false.
    !
    do i = 1, size(order)
       order(i) = i
    end do
    !
    call quick_sort_1(1, size(order))
    !
  contains
    !
    recursive subroutine quick_sort_1(left_end, right_end)
      !
      integer(FINT), intent(IN) :: left_end, right_end
      !
      !     Local variables
      integer(FINT)             :: reference
      integer(FINT)             :: i, j, itemp
      real(FREAL)               :: rtemp
      character(FLCHARS)         :: stemp
      integer(FINT), parameter  :: max_simple_sort_size = 6
      !
      if (right_end < left_end + max_simple_sort_size) then
         ! Use interchange sort for small molecule_zchars
         call interchange_sort(left_end, right_end)
         !
      else
         ! Use partition ("quick") sort
         itemp = (left_end + right_end)/2 !Get rid of f90/Sun compiler bug
         reference = molecule_zchar(itemp)
         i = left_end - 1; j = right_end + 1
         !
         do
            ! Scan molecule_zchar from left end until element >= reference is found
            do
               i = i + 1
               if (molecule_zchar(i) >= reference) exit
            end do
            ! Scan molecule_zchar from right end until element <= reference is found
            do
               j = j - 1
               if (molecule_zchar(j) <= reference) exit
            end do
            !
            if (i < j) then
               ! Swap two out-of-order elements
               !
               ! Swap Integer charges
               itemp = molecule_zchar(i)
               molecule_zchar(i) = molecule_zchar(j)
               molecule_zchar(j) = itemp
               ! Swap IDs
               itemp = molecule_atmID(i); molecule_atmID(i) = molecule_atmID(j); molecule_atmID(j) = itemp
               ! Swap Names
               stemp = molecule_name(i); molecule_name(i) = molecule_name(j); molecule_name(j) = stemp
               ! Swap Real charges
               rtemp = molecule_atmchar(i); molecule_atmchar(i) = molecule_atmchar(j); molecule_atmchar(j) = rtemp
               ! Swap Coordinates
               rtemp = molecule_xyz(1,i)
               molecule_xyz(1,i) = molecule_xyz(1,j)
               molecule_xyz(1,j) = rtemp
               rtemp = molecule_xyz(2,i)
               molecule_xyz(2,i) = molecule_xyz(2,j)
               molecule_xyz(2,j) = rtemp
               rtemp = molecule_xyz(3,i)
               molecule_xyz(3,i) = molecule_xyz(3,j)
               molecule_xyz(3,j) = rtemp
               ! Swap Index
               itemp = order(i); order(i) = order(j); order(j) = itemp
            else if (i == j) then
               i = i + 1
               exit
            else
               exit
            end if
         end do

         if (left_end < j) call quick_sort_1(left_end, j)
         if (i < right_end) call quick_sort_1(i, right_end)
      end if

    end subroutine quick_sort_1


    subroutine interchange_sort(left_end, right_end)

      integer(FINT), intent(IN) :: left_end, right_end

      !     Local variables
      integer(FINT)             :: i, j, itemp
      character(FLCHARS)         :: stemp
      real(FREAL)               :: rtemp

      do i = left_end, right_end - 1
         do j = i+1, right_end
            if (molecule_zchar(i) > molecule_zchar(j)) then
               ! Swap Integer charges
               itemp = molecule_zchar(i)
               molecule_zchar(i) = molecule_zchar(j)
               molecule_zchar(j) = itemp
               ! Swap IDs
               itemp = molecule_atmID(i); molecule_atmID(i) = molecule_atmID(j); molecule_atmID(j) = itemp
               ! Swap Names
               stemp = molecule_name(i); molecule_name(i) = molecule_name(j); molecule_name(j) = stemp
               ! Swap Real charges
               rtemp = molecule_atmchar(i); molecule_atmchar(i) = molecule_atmchar(j); molecule_atmchar(j) = rtemp
               ! Swap Coordinates
               rtemp = molecule_xyz(1,i)
               molecule_xyz(1,i) = molecule_xyz(1,j)
               molecule_xyz(1,j) = rtemp
               rtemp = molecule_xyz(2,i)
               molecule_xyz(2,i) = molecule_xyz(2,j)
               molecule_xyz(2,j) = rtemp
               rtemp = molecule_xyz(3,i)
               molecule_xyz(3,i) = molecule_xyz(3,j)
               molecule_xyz(3,j) = rtemp
               ! Swap Index
               itemp = order(i); order(i) = order(j); order(j) = itemp
            end if
         end do
      end do

    end subroutine interchange_sort

  end subroutine molecule_sort
!H
!H-----------------------------------------------------------------------------
!H Subroutine 
!H-----------------------------------------------------------------------------
!H
  recursive subroutine molecule_sortby(list,order,from,to)
    !
    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C.
    ! (1990) "Programmer's Guide to Fortran 90",
    ! McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an
    ! associated integer array which gives
    ! the positions of the elements in the original order.
    ! 
    ! Adapted to sort Atoms in FREEMOL by F.Mariotti.
    !
    integer(FINT), dimension (:), intent(OUT) :: order
    integer(FINT), dimension (:), intent(OUT) :: list
    integer(FINT), intent(in), optional :: from,to
    integer(FINT) :: ib,ie
    !
    ! Local variable
    integer(FINT) :: i
    !
    ib = 1
    ie = size(list)
    if(present(from)) then
       ib = from
    endif
    if(present(to)) then
       ie = to
    endif
    !
    !For md_atoms .. the distance matric must be compuoted again
    mol_isdist = .false.
    !
    do i = 1, size(list)
       order(i) = i
    end do
    !
    call quick_sort_1(1, size(list))
    !
  contains
    !
    recursive subroutine quick_sort_1(left_end, right_end)
      !
      integer(FINT), intent(IN) :: left_end, right_end
      !
      !     Local variables
      integer(FINT)             :: reference
      integer(FINT)             :: i, j, itemp, it, jt
      real(FREAL)               :: rtemp
      character(FLCHARS)         :: stemp
      integer(FINT), parameter  :: max_simple_sort_size = 6
      !
      if (right_end < left_end + max_simple_sort_size) then
         ! Use interchange sort for small lists
         call interchange_sort(left_end, right_end)
         !
      else
         ! Use partition ("quick") sort
         itemp = (left_end + right_end)/2 !Get rid of f90/Sun compiler bug
         reference = list(itemp)
         i = left_end - 1; j = right_end + 1
         !
         do
            ! Scan list from left end until element >= reference is found
            do
               i = i + 1
               if (list(i) >= reference) exit
            end do
            ! Scan list from right end until element <= reference is found
            do
               j = j - 1
               if (list(j) <= reference) exit
            end do
            !
            if (i < j) then
               ! Swap two out-of-order elements
               !
               ! Swap list
               itemp = list(i); list(i) = list(j); list(j) = itemp
               it = ib + i - 1
               jt = ib + j - 1
               ! Swap Integer charges
               itemp = molecule_zchar(it)
               molecule_zchar(it) = molecule_zchar(jt)
               molecule_zchar(jt) = itemp
               ! Swap IDs
               itemp = molecule_atmID(it)
               molecule_atmID(it) = molecule_atmID(jt)
               molecule_atmID(jt) = itemp
               ! Swap Names
               stemp = molecule_name(it)
               molecule_name(it) = molecule_name(jt)
               molecule_name(jt) = stemp
               ! Swap Real charges
               rtemp = molecule_atmchar(it)
               molecule_atmchar(it) = molecule_atmchar(jt)
               molecule_atmchar(jt) = rtemp
               ! Swap Coordinates
               rtemp = molecule_xyz(1,it)
               molecule_xyz(1,it) = molecule_xyz(1,jt)
               molecule_xyz(1,jt) = rtemp
               rtemp = molecule_xyz(2,it)
               molecule_xyz(2,it) = molecule_xyz(2,jt)
               molecule_xyz(2,jt) = rtemp
               rtemp = molecule_xyz(3,it)
               molecule_xyz(3,it) = molecule_xyz(3,jt)
               molecule_xyz(3,jt) = rtemp
               ! Swap Index
               itemp = order(i); order(i) = order(j); order(j) = itemp
            else if (i == j) then
               i = i + 1
               exit
            else
               exit
            end if
         end do

         if (left_end < j) call quick_sort_1(left_end, j)
         if (i < right_end) call quick_sort_1(i, right_end)
      end if

    end subroutine quick_sort_1


    subroutine interchange_sort(left_end, right_end)

      integer(FINT), intent(IN) :: left_end, right_end

      !     Local variables
      integer(FINT)             :: i, j, itemp, it, jt
      character(FLCHARS)         :: stemp
      real(FREAL)               :: rtemp

      do i = left_end, right_end - 1
         do j = i+1, right_end
            if (list(i) > list(j)) then
               ! Swap list
               itemp = list(i); list(i) = list(j); list(j) = itemp
               it = ib + i - 1
               jt = ib + j - 1
               ! Swap Integer charges
               itemp = molecule_zchar(it)
               molecule_zchar(it) = molecule_zchar(jt)
               molecule_zchar(jt) = itemp
               ! Swap IDs
               itemp = molecule_atmID(it); molecule_atmID(it) = molecule_atmID(jt); molecule_atmID(jt) = itemp
               ! Swap Names
               stemp = molecule_name(it)
               molecule_name(it) = molecule_name(jt)
               molecule_name(jt) = stemp
               ! Swap Real charges
               rtemp = molecule_atmchar(it)
               molecule_atmchar(it) = molecule_atmchar(jt)
               molecule_atmchar(jt) = rtemp
               ! Swap Coordinates
               rtemp = molecule_xyz(1,it)
               molecule_xyz(1,it) = molecule_xyz(1,jt)
               molecule_xyz(1,jt) = rtemp
               rtemp = molecule_xyz(2,it)
               molecule_xyz(2,it) = molecule_xyz(2,jt)
               molecule_xyz(2,jt) = rtemp
               rtemp = molecule_xyz(3,it)
               molecule_xyz(3,it) = molecule_xyz(3,jt)
               molecule_xyz(3,jt) = rtemp
               ! Swap Index
               itemp = order(i); order(i) = order(j); order(j) = itemp
            end if
         end do
      end do

    end subroutine interchange_sort

  end subroutine molecule_sortby
!
!H
!H-----------------------------------------------------------------------------
!H subroutine molecule_read(iunin,params)
!H-----------------------------------------------------------------------------
!H
!
  subroutine molecule_read(iunin,params)
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine parse input data from a molden/freemol/Freemol file (Real)
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    integer(FINT), parameter :: MAXREC = MDATOMSMAXDIM
    !
    integer(FINT), intent(in) :: iunin
    character(FLCHARS), intent(in) :: params
    !
    integer(FINT) :: irc
    !
    ! Process Parameters
    !-------------------
    call molecule_parspar(params)
    !
    ! Parspar allocate only if it has nrec in the parameters.
    ! So we check if we are allocated.
    !---------------------------------
    if(.not.molecule_isallocated()) then
       irc = molecule_allocate()
       if (irc.lt.0) then
          call message(MESERRO,"[molecule_read] cannot allocate molecule.")
          molecule_error = 1
          return
       end if
    end if
    !
    ! Format not implemented so:
    !call molecule_readascii(iunin,mol_frmt)
    call molecule_readascii(iunin)
    !
  end subroutine molecule_read

!H
!H-----------------------------------------------------------------------------
!H Subroutine molecule_parspar(params)
!H-----------------------------------------------------------------------------
!H
  subroutine molecule_parspar(params)
    !
    !
    integer(FINT), parameter :: MAXREC = MDATOMSMAXDIM
    !
    character(FLCHARS), intent(in) :: params
    !
    integer(FINT) :: irc
    ! NREC
    integer(FINT) :: nrec
    character(FLCHARS) :: cnrec
    ! FRMT
    character(FLCHARS) :: frmt
    ! Others
    !character(FLCHARS) :: cdummy
    ! NREC
    nrec=MAXREC
    irc = line_getval(params,'nrec',cnrec)
    if(irc.ge.0) then
       read(cnrec,*) nrec
    end if
    call message_value(MESDEBG,'[molecule_parspar] Number of Records=',nrec)
    !
    !We spread it out to molecule module and allocate it
    irc = molecule_setnat(nrec)
    irc = molecule_allocate(nrec)
    ! FRMT
    frmt='xyz'
    irc = line_getval(params,'frmt',frmt)
    if(irc.lt.0) then
       irc = line_getval(params,'format',frmt)
    end if
    frmt=adjustl(frmt)
    call message_value(MESDEBG,'[molecule_parspar] Format=',frmt)
    !
    ! Format in not yet implemented so we keep going.
    !------------------------------------------------
    mol_frmt = frmt
    return
    !
  end subroutine molecule_parspar
!
!H
!H-----------------------------------------------------------------------------
!H subroutine molecule_rotxyz(tx,ty,tz,rx,ry,rz)
!H-----------------------------------------------------------------------------
!H
!
  subroutine molecule_rotxyz(tx,ty,tz,rx,ry,rz,vangle)
    !
    real(FREAL), intent(in) :: tx,ty,tz
    real(FREAL), intent(in) :: rx,ry,rz
    real(FREAL), intent(in) :: vangle
    !
    integer(FINT) :: i,nat
    real(FREAL), dimension(3,3) :: rmat
    real(FREAL), dimension(3) :: taxe
    real(FREAL) :: eps = 1.0E-8_FREAL
    real(FREAL) :: alpha,beta,gamma,angle
    !
    ! Transform in radians
    !---------------------
    angle = vangle * (4.0_KREAL*atan(1.0_KREAL)) / 180.0_KREAL
    !
    ! Translate First
    !----------------
    molecule_xyz(1,:) = tx + molecule_xyz(1,:)
    molecule_xyz(2,:) = ty + molecule_xyz(2,:)
    molecule_xyz(3,:) = tz + molecule_xyz(3,:)
    !
    ! Get rotational matrixes
    !------------------------
    taxe(1) = rx
    taxe(2) = ry
    taxe(3) = rz
    call mathtools_euler(angle,taxe(1:3),&
         &.false.,alpha,beta,gamma,rmat(1:3,1:3),eps)
    !
    ! Some debugging
    !---------------
    if(.false.) then
       write(0,'("Required Rotation")')
       write(0,'(3(1X,F12.4))') tx,ty,tz
       write(0,'(3(1X,F12.4))') rx,ry,rz
       write(0,'("Final Rotation Matrix")')
       write(0,'(3(1X,F12.4))') rmat(:,:)
    endif
    !
    ! Apply rotation
    !---------------
    nat = molecule_getnat()
    do i=1,nat
       molecule_xyz(1:3,i) = matmul(rmat(1:3,1:3),molecule_xyz(1:3,i))
    end do
    !
    !

  end subroutine molecule_rotxyz
!
!H
!H-----------------------------------------------------------------------------
!H Subroutine 
!H-----------------------------------------------------------------------------
!H
!


end module molecule

