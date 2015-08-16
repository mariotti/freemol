!
!H
!H  
!H----------------------------------------------------------------------
!H----------------------------------------------------------------------
!H  MODULE adf_module for adfrom program
!H----------------------------------------------------------------------
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H----------------------------------------------------------------------
!H
!H  Routines related to ADF TAPE21 keyed file
!H
!H----------------------------------------------------------------------
!H
!
module adfmodule
  !
  use kf
  use vartypes
  use adfrommodule
  !
  implicit none
  private
  !
  ! controls
  integer(FINT), public, save :: adf_errors = 0
  !
  ! ADF TAPE21 DATA
  !
  ! General Section: NOTE ADF LCHARS, KINT and KREAL are used
  !
  character(LCHARS), public, save :: title
  character(LCHARS), public, save :: runtype
  character(LCHARS), public, save :: jobid
  integer(KINT), public, save :: nspin
  real(KREAL), public, save :: electrons
  !
  ! Geometry Section: NOTE ADF LCHARS, KINT and KREAL are used
  !
  integer(KINT), public, save :: nnuc
  integer(KINT), public, save :: ntyp
  integer(KINT), public, save, allocatable,dimension(:) :: nqptr
  integer(KINT), public, save, allocatable,dimension(:) :: natyp
  character(FLCHARS), public, save, allocatable, dimension(:) :: atomtype
  real(KREAL), public, save, allocatable, dimension(:,:) :: xyz
  real(KREAL), public, save, allocatable, dimension(:) :: qtch
  !
  ! Basis Section: NOTE ADF LCHARS, KINT and KREAL are used
  !
  integer(KINT), public, save :: naos
  integer(KINT), public, save :: nbos
  integer(KINT), public, save, allocatable,dimension(:) :: nbptr
  integer(KINT), public, save, allocatable,dimension(:,:) :: stoTint
  real(KREAL), public, save, allocatable,dimension(:,:) :: stoTreal
  !
  ! Symmetry Section: NOTE ADF LCHARS, KINT and KREAL are used
  !
  integer(KINT), public, save :: nsym
  integer(KINT), public, save :: npeq
  integer(KINT), public, save :: imos
  integer(KINT), public, save, allocatable,dimension(:) :: nfcn
  integer(KINT), public, save, allocatable,dimension(:,:) :: npart
  integer(KINT), public, save, allocatable,dimension(:) :: norb
  integer(KINT), public, save, allocatable,dimension(:,:) :: nmos
  character(FLCHARS), public, save, allocatable, dimension(:) :: symlab
  real(KREAL), public, save, allocatable,dimension(:,:,:) :: froc
  real(KREAL), public, save, allocatable,dimension(:,:,:) :: eigval
  real(KREAL), public, save, allocatable,dimension(:,:,:,:) :: eigvec
  !Symmetry Indexes
  integer(KINT), public, save, allocatable,dimension(:,:,:) :: neneidx
  integer(KINT), public, save, allocatable,dimension(:,:) :: nenesort
  real(KREAL), public, save, allocatable,dimension(:,:) :: orbenergy
  !
  
  !
  !
  public adf99_init
  public adf99_deallocate
  public read_general
  public read_geometry
  public read_basis
  public read_symmetry
  public write_general
  public write_geometry
  public write_basis
  public write_symmetry
  !
contains
  !
  subroutine adf99_init(unit_in,file_in)
    !
    integer(FINT), intent(in) :: unit_in
    character(FLCHARS), intent(in) :: file_in
    !
    if (.not.kfexfl(file_in)) then
       write(0,*) 'Error: use messages for this code.'
       write(0,*) 'Error: Input file does not exist.'
       stop 1
    endif
    !
    write(0,*) 'Error: fix it: kfopfl does not skip std units.!!'
    call kfopfl(unit_in,file_in)
    !
  end subroutine adf99_init
  !
  !
  subroutine adf99_deallocate()
    !
    write (0,*) ' Insert adf99_deallocate '
    !
  end subroutine adf99_deallocate
  !
  !
  subroutine read_general(unit_in)
    !
    integer(FINT), intent(in) :: unit_in
    !
    ! Code
    !
    call adfmessage('General Section Begin ... ')
    !
    ! Open General Section
    !
    call kfopsc(unit_in,'General')
    !
    ! Read Data
    !
    call kfread(unit_in,'jobid',jobid)
    call kfread(unit_in,'title',title)
    call kfread(unit_in,'runtype',runtype)
    call kfread(unit_in,'nspin',nspin)
    call kfread(unit_in,'electrons',electrons)
    !
    ! Close Section
    !
    call kfclsc(unit_in)
    !
    ! Write Out Infos
    !
    call adfmessage(' Title: '//title)
    call adfmessage(' Job ID: '//jobid)
    call adfmessage(' RUN TYPE: '//runtype)
    !
    if (nspin.eq.1) then
       call adfmessage(' Run Mode: '//'RESTRICTED CALCULATION')
    else
       call adfmessage(' Run Mode: '//'UNRESTRICTED CALCULATION')
    end if
    call adfmessage(' Valence Electrons: ')
    call adfmessage(electrons)
    !
    !
    call adfmessage(' ... General Section end.')
    !
  end subroutine read_general
  !
  !
  subroutine read_geometry(unit_in)
    !
    integer(FINT), intent(in) :: unit_in
    !
    integer(FINT) :: i
    integer(FINT) :: irc
    real(KREAL), dimension(:), allocatable :: rtemp
    !
    ! Code
    !
    call adfmessage('Geometry Section Begin ... ')
    !
    ! Open Geometry Section
    !
    call kfopsc(unit_in,'Geometry')
    !
    ! Read Data
    !
    call kfread(unit_in,'nnuc',nnuc)
    !
    allocate(rtemp(3*nnuc),STAT=irc)
    call adfstopon(irc.ne.0,'ERROR: in read_geometry: cannot allocate rtemp')
    call kfread(unit_in,'xyz',rtemp(:))
    allocate(xyz(3,nnuc),STAT=irc)
    call adfstopon(irc.ne.0,'ERROR: in read_geometry: cannot allocate xyz')
    xyz = reshape(rtemp,(/3,int(nnuc)/))
    deallocate(rtemp,STAT=irc)
    call adfstopon(irc.ne.0,'ERROR: in read_geometry: cannot deallocate rtemp')
    !
    call kfread(unit_in,'ntyp',ntyp)
    !
    allocate(atomtype(ntyp),STAT=irc)
    call adfstopon(irc.ne.0,' cannot allocate atomtype ')
    !
    allocate(qtch(ntyp),STAT=irc)
    call adfstopon(irc.ne.0,' cannot allocate qtch ')
    !
    allocate(natyp(ntyp),STAT=irc)
    call adfstopon(irc.ne.0,' cannot allocate natyp ')
    !
    allocate(nqptr(ntyp+1),STAT=irc)
    call adfstopon(irc.ne.0,' cannot allocate nqptr ')
    !
    call kfread(unit_in,'atomtype',atomtype)
    call kfread(unit_in,'qtch',qtch)
    call kfread(unit_in,'nqptr',nqptr)

    !
    ! Close Section
    !
    call kfclsc(unit_in)
    !
    ! Write Out Infos
    !
    call adfmessage('Number of nuclei')
    call adfmessage(nnuc)
    call adfmessage('Coordinates')
    do i=1,nnuc
       call adfmessage(xyz(1:3,i))
    enddo
    !
    call adfmessage(' ... Geometry Section end.')
    !
  end subroutine read_geometry
  !
  !
  subroutine read_basis(unit_in)
    !
    integer(FINT), intent(in) :: unit_in
    !
    integer(FINT) :: irc
    !
    ! Code
    !
    call adfmessage('Basis Section begin ... ')
    !
    ! Open Basis Section
    !
    call kfopsc(unit_in,'Basis')
    !
    allocate(nbptr(ntyp+1),STAT=irc)
    call AdfStopOn(irc.ne.0,'Cannot allocate nbptr')
    !allocate(,STAT=irc)
    !call AdfStopOn(irc.ne.0,'')
    !
    call kfread(unit_in,'nbptr',nbptr)
    call kfread(unit_in,'naos',naos)
    call kfread(unit_in,'nbos',nbos)
    !
    allocate(stoTint(4,nbos),STAT=irc)
    call AdfStopOn(irc.ne.0,'Cannot allocate stoTint')
    !
    call kfread(unit_in,'kx',stoTint(1,:))
    call kfread(unit_in,'ky',stoTint(2,:))
    call kfread(unit_in,'kz',stoTint(3,:))
    call kfread(unit_in,'kr',stoTint(4,:))
    !
    allocate(stoTreal(2,nbos),STAT=irc)
    call AdfStopOn(irc.ne.0,'Cannot allocate stoTreal')
    !
    call kfread(unit_in,'alf',stoTreal(1,:))
    call kfread(unit_in,'bnorm',stoTreal(2,:))
    !
    ! Close Section
    !
    call kfclsc(unit_in)
    !
    ! Write Out Infos
    !
    call adfmessage('Number of atomic orbitals')
    call adfmessage(naos)
    call adfmessage('Number different of atomic orbitals')
    call adfmessage(nbos)
    !
    call adfmessage(' ... Basis Section end.')
    !
  end subroutine read_basis
  !
  !
  subroutine read_symmetry(unit_in)
    !
    use adfmathtools
    !
    integer(FINT), intent(in) :: unit_in
    !
    integer(FINT) :: ind
    integer(FINT) :: k
    integer(FINT) :: ls
    integer(FINT) :: irc
    integer(FINT) :: ispin
    integer(FINT) :: isym
    integer(FINT) :: nsfcn
    integer(FINT) :: nsos
    integer(FINT) :: isos
    integer(KINT) :: ntotfcn
    integer(KINT) :: nmaxfcn
    character(LCHARS) :: keysec
    character(LCHARS) :: keynpart
    character(LCHARS) :: keynmo
    character(LCHARS) :: keyfroc
    character(LCHARS) :: keyeps
    character(LCHARS) :: keyeigen
    !
    real(KREAL), allocatable, dimension(:) :: rwrk
    !
    ! Code
    !
    call adfmessage('Symmetry Section begin ... ')
    !
    call kfopsc(unit_in,'Symmetry')
    !
    call kfread(unit_in,'nsym',nsym)
    call kfread(unit_in,'npeq',npeq)
    !
    allocate(symlab(nsym),STAT=irc)
    call AdfStopOn(irc.ne.0,' Cannot allocate symlab.')
    allocate(nfcn(nsym),STAT=irc)
    call AdfStopOn(irc.ne.0,' Cannot allocate nfcn.')
    allocate(norb(nsym),STAT=irc)
    call AdfStopOn(irc.ne.0,' Cannot allocate norb.')
    !
    !call kfrdns(unit_in,'symlab',rep,nsym,1)
    call kfread(unit_in,'symlab',symlab,nsym)
    call kfread(unit_in,'nfcn',nfcn)
    call kfread(unit_in,'norb',norb)
    !
    ntotfcn = sum(nfcn(1:nsym))
    allocate(npart(ntotfcn,nsym),STAT=irc)
    call AdfStopOn(irc.ne.0,' Cannot allocate npart.')
    !
    allocate(nmos(nsym,nspin),STAT=irc)
    call AdfStopOn(irc.ne.0,' Cannot allocate nmos.')
    !
    ! we suppose (as it must be) that the MOs are .le. NAOS!!
    nmaxfcn = maxval(nfcn(1:nsym))
    allocate(froc(nsym,nmaxfcn,nspin),STAT=irc)
    call AdfStopOn(irc.ne.0,' Cannot allocate froc.')
    !
    allocate(eigval(nsym,nmaxfcn,nspin),STAT=irc)
    call AdfStopOn(irc.ne.0,' Cannot allocate eigval.')
    !
    allocate(eigvec(nsym,nmaxfcn,nmaxfcn,nspin),STAT=irc)
    call AdfStopOn(irc.ne.0,' Cannot allocate eigvec.')
    !
    allocate(rwrk(nmaxfcn*nmaxfcn),STAT=irc)
    call AdfStopOn(irc.ne.0,' Cannot allocate rwrk.')
    !
    keynpart = 'npart'
    do ispin = 1,nspin
       if (ispin.eq.1) then
          keynmo = 'nmo_A'
          keyfroc = 'froc_A'
          keyeps = 'eps_A'
          keyeigen = 'Eigen-Bas_A'
       else
          keynmo = 'nmo_B'
          keyfroc = 'froc_B'
          keyeps = 'eps_B'
          keyeigen = 'Eigen-Bas_B'
       end if
       do isym=1,nsym
          nsfcn = nfcn(isym)
          ls = len_trim(symlab(isym))
          keysec = symlab(isym)(1:ls)
          call kfopsc(unit_in,keysec(1:ls))
          if (ispin.eq.1) then
             call kfread(unit_in,'npart',npart(1:nsfcn,isym),nsfcn,1)
          end if
          call kfread(unit_in,keynmo,nmos(isym,ispin))
          nsos=nmos(isym,ispin)
          call kfread(unit_in,keyfroc,froc(isym,1:nsos,ispin),nsos)
          call kfread(unit_in,keyeps,eigval(isym,1:nsos,ispin),nsos)
          call kfread(unit_in,keyeigen,rwrk(1:nsos*nsfcn),nsos*nsfcn)
          ind=0
          do isos = 1, nsos
             do k=1,nsfcn
                ind = ind+1
                eigvec(isym,k,isos,ispin)=rwrk(ind)
             enddo
          enddo
       enddo
    enddo
    !
    !     Sorting the eigenvalues in ascending energy ordering.
    !
    allocate(orbenergy(nsym*nmaxfcn,nspin),STAT=irc)
    call AdfStopOn(irc.ne.0,' Cannot allocate orbenergy.')
    !
    allocate(neneidx(nsym*nmaxfcn,2,nspin),STAT=irc)
    call AdfStopOn(irc.ne.0,' Cannot allocate neneidx.')
    !
    allocate(nenesort(nsym*nmaxfcn,nspin),STAT=irc)
    call AdfStopOn(irc.ne.0,' Cannot allocate nenesort.')
    !
    imos = 0
    do ispin=1,nspin
       do isym=1,nsym
          nsos = nmos(isym,ispin)
          do isos=1,nsos
             imos = imos + 1
             orbenergy(imos,ispin) = eigval(isym,isos,ispin)
             neneidx(imos,1,ispin) = isym
             neneidx(imos,2,ispin) = isos
          enddo
       enddo
       call quicksortidx(imos,orbenergy(1:imos,ispin),nenesort(1:imos,ispin))
    end do
    !
    call kfclsc(unit_in)
    !
    call adfmessage(' ... Symmetry Section end.')
    !
  end subroutine read_symmetry
  !
  !
  subroutine write_general(unit_out)
    !
    integer(FINT), intent(in) :: unit_out
    !
    write(unit_out,'(A8)') '[ADFROM]'
    write(unit_out,'(A7)') '[Title]'
    write(unit_out,'(A80)') title
    write(unit_out,'(A80)') jobid
    write(unit_out,'(A80)') runtype
    if(nspin.eq.1) then
       write(unit_out,'(A10)') 'Restricted'
    else
       write(unit_out,'(A12)') 'Unrestricted'
    endif
    write(unit_out,'(A10,F15.8)') 'Electrons=',electrons
    !
  end subroutine write_general
  !
  !
  subroutine write_geometry(unit_out)
    !
    use chemconst
    !
    integer(FINT), intent(in) :: unit_out
    !
    integer(FINT) :: inuc
    integer(FINT) :: ityp
    integer(FINT) :: iatyp
    !
    write(unit_out,'(A13,I5)') '[Atoms] NREC=',nnuc
    !
    inuc = 0
    do ityp=1,ntyp
       do iatyp=nqptr(ityp),nqptr(ityp+1)-1
          !write(*,*) ityp,ntyp,iatyp,
          inuc = inuc + 1
          write (unit_out,'(1X,A10,1X,I4,1X,I4,1X,3(1X,E15.8))')&
               &atomtype(ityp)(1:10),&
               &inuc,int(qtch(ityp)),xyz(1:3,inuc)*cc_au2ang
          call flush(unit_out)
       enddo
    enddo
    !
  end subroutine write_geometry
  !
  !
  subroutine write_basis(unit_out,iopt)
    !
    use adfrompars
    !
    integer(FINT), intent(in) :: unit_out
    integer(FINT), intent(in), optional :: iopt
    integer(FINT) :: lopt
    !
    if (present(iopt)) then
       lopt = iopt
    else
       lopt = STO_FMT_MOLDEN !default format for molden
    endif
    !
    select case ((lopt))
    case(STO_FMT_CART_SHORT)     !cart sto per atom type
       write (0,*) 'Not yet implemented'
       stop 1
    case (STO_FMT_CART_LONG)     !STO_FMT_MOLDEN
       call OutBasisMolden(unit_out)
    case (STO_FMT_SPHER_SHORT)   !spherical sto short
       write (0,*) 'Not yet implemented'
       stop 1
    case (STO_FMT_SPHER_LONG)    !spherical sto long
       write (0,*) 'Not yet implemented'
       stop 1
    case (STO_FMT_CART_ZCONTR)   !with zeta contraction
       write (0,*) 'Not yet implemented'
       stop 1
    case (STO_FMT_SPHER_ZCONTR)  !with zeta contraction
       write (0,*) 'Not yet implemented'
       stop 1
    case (STO_FMT_CART_SZCONTR)  !short zeta contraction
       write (0,*) 'Not yet implemented'
       stop 1
    case (STO_FMT_SPHER_SZCONTR) !short zeta contraction
       call AdfStop('sto format not yet implemented')
    case default
       call AdfStop('Internal Error in sto formats')
    end select
    !

    !
  end subroutine write_basis
  !
  !
  subroutine OutBasisMolden(unit_out)
    !
    integer(FINT), intent(in) :: unit_out
    !
    integer(FINT) :: inuc
    integer(FINT) :: ityp
    integer(FINT) :: iatyp
    integer(FINT) :: iaos
    integer(FINT) :: ibtyp
    !
    write(unit_out,'(A11,1X,I7)') '[STO] NREC=',naos
    write(unit_out,'(A71)') '# atom    kx    ky    kz    kr alpha           bnorm         # from ADF'
    !
    iaos = 0
    inuc = 0
    do ityp=1,ntyp
       do iatyp=nqptr(ityp),nqptr(ityp+1)-1
          inuc = inuc + 1
          do ibtyp=nbptr(ityp),nbptr(ityp+1)-1
             iaos=iaos+1
             write(unit_out,'(5(1X,I5),1X,2(1X,E15.8),1X,I5)')&
                  &inuc,stoTint(1:4,ibtyp),stoTreal(1:2,ibtyp),iaos
          enddo
       enddo
    enddo
    !
  end subroutine OutBasisMolden
  !
  !
  subroutine write_symmetry(unit_out,iopt)
    !
    use adfrompars
    !
    integer(FINT), intent(in) :: unit_out
    integer(FINT), intent(in), optional :: iopt
    integer(FINT) :: lopt
    !
    integer(FINT) :: i
    !
    if (present(iopt)) then
       lopt = iopt
    else
       lopt = MO_FMT_MOLDEN !default format for molden
    endif
    !
    write(unit_out,'(A16,I3)') '[SYMMETRY] NREC=',nsym
    !
    if(nspin.eq.1) then
       write(unit_out,'("#",A3,1X,A10,2(1X,A6))')&
            &'num','SymLbl    ','# AOs','# MOs'
       do i=1,nsym
          write(unit_out,'(1X,I3,1X,A10,2(1X,I6))')&
               &i,symlab(i),nfcn(i),nmos(i,nspin)
       enddo
    else
       write(unit_out,'(1X,A3,1X,A10,3(1X,A6))')&
            &'num','SymLbl    ','# AOs','# MOsA','# MOsB'
       do i=1,nsym
          write(unit_out,'(1X,I3,1X,A10,3(1X,I6))')&
               &i,symlab(i),nfcn(i),nmos(i,1),nmos(i,2)
       enddo
    endif
    !
    select case ((lopt))
    case(MO_FMT_LINES_CONTR) !MO_FMT_MOLDEN
       call OutMOMolden(unit_out)
    case (MO_FMT_LIST_CONTR)
       call AdfStop('mo format not yet implemented')
    case default
       call AdfStop('Internal Error in mo formats')
    end select
    !
  end subroutine write_symmetry
  !
  subroutine OutMOMolden(unit_out)
    !
    integer(FINT), intent(in) :: unit_out
    !
    character(5) :: cspin
    integer(KINT) :: i
    integer(KINT) :: ispin
    integer(KINT) :: ifcn
    integer(KINT) :: isym
    integer(KINT) :: kmos
    !
    write(unit_out,'(A10,I5)') '[MO] NREC=',imos
    do ispin=1,nspin
       if(ispin.eq.1) then
          cspin = 'Alpha'
       else
          cspin = 'Beta '
       endif
       do i=1,imos
          isym = neneidx(nenesort(i,ispin),1,ispin)
          kmos = neneidx(nenesort(i,ispin),2,ispin)
          write(unit_out,'(1X,"Sym=",I3,A10,2X,I3)')&
               &kmos,symlab(isym),isym
          write(unit_out,'(1X,"Ene=",E15.8)')&
               &eigval(isym,kmos,ispin)
          write(unit_out,'(1X,"Spin=",A5,1X,I1)')&
               &cspin,ispin
          write(unit_out,'(1X,"Occup=",E15.8)')&
               &froc(isym,kmos,ispin)
          do ifcn=1,nfcn(isym)
             write(unit_out,'(1X,I5,2X,E15.8)')&
                  &npart(ifcn,isym),eigvec(isym,ifcn,kmos,ispin)
          enddo
       enddo
    enddo
    !
  end subroutine OutMOMolden
  !
end module adfmodule
