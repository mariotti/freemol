!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H MODULE metpot4
!H-----------------------------------------------------------------------------
!H This File is part of the Frimol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H $Id: metpot4.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
!H-----------------------------------------------------------------------------
!H
!H
!
module metpot4
!
!H
!H Original copyright statment
!H
!H ============================================================================
!H
!H  ! potential model METPOT4 for methane and subgroups
!H  ! version 31AUG97
!H  !  authors
!H  !         ----------------------------------------
!H  !        | Roberto Marquardt and Martin Quack     |
!H  !        | Laboratorium fuer Physikalische Chemie |
!H  !        | ETH-Zuerich(Zentrum)                   |
!H  !        | 8092 Zuerich                           |
!H  !        |                                        |
!H  !        | roma@phys.chem.ethz.ch                 |
!H  !         ----------------------------------------
!H  !-----------------------------------------------------------------------
!H
!H ============================================================================
!H
  !
  use vartypes
  use messages
  !
  implicit none
  private
  !
  ! PUBLIC DECLARATIONS:
  !--------------------

  ! Public Procedures
  !------------------
  public metpot4_init
  public metpot4_isinit
  public metpot4_initpot
  public metpot4_vcar
  public metpot4_vint
  !
  ! Public Vars
  !------------
  character(LCHARS), public :: metpot4_version
  !
  ! Private Variables
  !------------------
  !
  logical, save :: LOCmetpot4_isinit = .false.
  !
  ! old parameters definition
  !--------------------------
  real(FREAL), parameter :: one = 1.0_FREAL
  real(FREAL), parameter :: small = 1.0E-06_FREAL
  real(FREAL), parameter :: swthre = 200.0_FREAL
  real(FREAL), parameter :: ajoule = 50341.12503_FREAL
  real(FREAL), parameter :: qjmol = 83.59346114_FREAL
  real(FREAL), parameter :: abohr = 0.529177249_FREAL
  !
  ! old commons decls changed to save modules datatype
  !---------------------------------------------------
  !
  ! PARDEF old commons and dimensions.
  ! i.e. parameters definitions for the potential
  ! I will try to comment them on code base....
  ! These are hard coded at present in metpot4_initpot()
  !-----------------------------------------------------
  integer(FINT), parameter :: n = 1
  integer(FINT), parameter :: ni = 4
  integer(FINT), parameter :: mi = ni - 1
  integer(FINT), parameter :: nni = ni*(ni-1)/2
  integer(FINT), parameter :: nahb = 40
  integer(FINT), parameter :: nad = 3
  ! V0 used only in calling SWIT4 needs to be domentioned!
  !-----------------------------------------------------
  real(FREAL), save, allocatable, dimension(:) :: pars_v0
  real(FREAL), save, allocatable, dimension(:) :: pars_ar
  real(FREAL), save, allocatable, dimension(:) :: pars_re
  real(FREAL), save, allocatable, dimension(:) :: pars_fs
  real(FREAL), save, allocatable, dimension(:) :: pars_as
  real(FREAL), save, allocatable, dimension(:) :: pars_e6
  real(FREAL), save, allocatable, dimension(:) :: pars_e8
  real(FREAL), save, allocatable, dimension(:) :: pars_rs
  real(FREAL), save, allocatable, dimension(:) :: pars_ac
  real(FREAL), save, allocatable, dimension(:) :: pars_ce
  real(FREAL), save, allocatable, dimension(:) :: pars_fb0
  real(FREAL), save, allocatable, dimension(:) :: pars_fb1
  real(FREAL), save, allocatable, dimension(:) :: pars_fb2
  real(FREAL), save, allocatable, dimension(:) :: pars_fb3
  real(FREAL), save, allocatable, dimension(:) :: pars_dii
  real(FREAL), save, allocatable, dimension(:) :: pars_aii
  real(FREAL), save, allocatable, dimension(:,:) :: pars_ad0
  real(FREAL), save, allocatable, dimension(:,:) :: pars_ad1
  real(FREAL), save, allocatable, dimension(:,:) :: pars_ad2
  real(FREAL), save, allocatable, dimension(:,:) :: pars_ad3
  real(FREAL), save, allocatable, dimension(:,:) :: pars_ad4
  real(FREAL), save, allocatable, dimension(:,:) :: pars_ab
  !
  integer(FINT), save :: pars_nsw
  real(FREAL), save :: pars_ssw
  integer(FINT), save :: pars_nrjj
  real(FREAL), save :: pars_srjj
  integer(FINT), save :: pars_ndjj
  real(FREAL), save :: pars_sdjj
  integer(FINT), save :: pars_najj
  real(FREAL), save :: pars_sajj
  real(FREAL), save :: pars_riii
  real(FREAL), save :: pars_aiii
  real(FREAL), save :: pars_diii
  !
  ! dummies
  integer(FINT) :: irc
contains
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H Function metpot4_init()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function metpot4_init()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the metpot4 module.
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    ! default to return with no errors
    !---------------------------------
    metpot4_init = 1
    !
    ! set the current version of metpot4: old stile from marquardt and quack
    !-----------------------------------------------------------------------
    metpot4_version = 'VERS. 53.970831'
    !
    ! Initialize the potential: Internally stored
    !--------------------------------------------
    call metpot4_initpot() !Old pardef oddly called by old VINI
    !
    LOCmetpot4_isinit = .true.
    !
  end function metpot4_init
!
!
!H
!H-----------------------------------------------------------------------------
!H logical function metpot4_isinit()
!H-----------------------------------------------------------------------------
!H
!
  logical function metpot4_isinit()
!
!H
!H-----------------------------------------------------------------------------
!H We merelyreturn the initialization flag.
!H-----------------------------------------------------------------------------
!H
!
    !
    metpot4_isinit = LOCmetpot4_isinit
    !
  end function metpot4_isinit
!-----------------------------------------------------------------------
!
!H
!H-----------------------------------------------------------------------------
!H subroutine metpot4_vint(RINT,V)
!H-----------------------------------------------------------------------------
!H
!
  subroutine metpot4_vint(RINT,V)
    !
!! ! nominal parameters
!! ! units are aJ=10**(-18)J for energies and A=100pm for lengths
!!       common /PAR0/   V0(1:NI)
!!       common /PARS1/  AR(1:NI),RE(1:NI)
!!       common /PARS2/  FS(1:NI)
!!       common /PARS3/  AS(1:NI)
!!       common /PARS4/  E6(1:NI)
!!       common /PARS5/  E8(1:NI)
!!       common /PARS6/  RS(1:NI)
!!       common /PARB1/  AC(1:MI),CE(1:MI)
!!       common /PARB2/  FB0(1:MI)
!!       common /PARB3/  FB1(1:MI)
!!       common /PARB4/  FB2(1:MI)
!!       common /PARB5/  FB3(1:MI)
!!       common /PARB6/  AD0(1:NAD,1:MI)
!!       common /PARB8/  AD1(1:NAD,1:MI)
!!       common /PARB9/  AD2(1:NAD,1:MI)
!!       common /PARB10/ AD3(1:NAD,1:MI)
!!       common /PARB11/ AD4(1:NAD,1:MI)
!!       common /PARB12/ AB(1:NAHB,1:MI)
!!       common /PARI1/  DII(1:MI)
!!       common /PARI2/  AII(1:MI)
!!       common /PARI3/  NRJJ,SRJJ,RIII
!!       common /PARI4/  NAJJ,SAJJ,AIII
!!       common /PARI5/  NDJJ,SDJJ,DIII
!!       common /PARSW/  NSW,SSW
! local parameters
    real(FREAL), dimension(1:n), intent(out) :: v
    real(FREAL), dimension(1:n,1:10), intent(in) :: rint
    !
    real(FREAL), dimension(n,ni) :: r
    real(FREAL), dimension(n,ni) :: r2
    real(FREAL), dimension(n,ni) :: dr
    real(FREAL), dimension(n,nni) :: c
    real(FREAL), dimension(n,nni) :: dc
    real(FREAL), dimension(n,nni) :: d
    real(FREAL), dimension(n,nni) :: dd
    real(FREAL), dimension(n,ni) :: spsw
    real(FREAL), dimension(n,nahb,ni) :: spab
    real(FREAL), dimension(n,nad,ni) :: spad
    real(FREAL), dimension(n,ni) :: sprjj
    real(FREAL), dimension(n,ni) :: spdjj
    real(FREAL), dimension(n,ni) :: spajj
    real(FREAL), dimension(n,nni) :: spd3
    real(FREAL), dimension(n,ni) :: spre
    real(FREAL), dimension(n,ni) :: spce
    real(FREAL), dimension(n,ni) :: sp6
    real(FREAL), dimension(n,ni) :: sp8
    real(FREAL), dimension(n) :: v0j
    real(FREAL), dimension(n) :: req
    real(FREAL), dimension(n) :: ceq
    real(FREAL), dimension(n,ni) :: rej
    real(FREAL), dimension(n,ni) :: fsj
    real(FREAL), dimension(n,ni) :: asj
    real(FREAL), dimension(n,ni) :: e6j
    real(FREAL), dimension(n,ni) :: e8j
    real(FREAL), dimension(n,ni) :: rsj
    real(FREAL), dimension(n,nni) :: cej
    real(FREAL), dimension(n) :: fb0j
    real(FREAL), dimension(n) :: fb1j
    real(FREAL), dimension(n) :: fb2j
    real(FREAL), dimension(n) :: fb3j
    real(FREAL), dimension(n,nad) :: ad0j
    real(FREAL), dimension(n,nad) :: ad1j
    real(FREAL), dimension(n,nad) :: ad2j
    real(FREAL), dimension(n,nad) :: ad3j
    real(FREAL), dimension(n,nad) :: ad4j
    real(FREAL), dimension(n,nahb) :: abj
    real(FREAL), dimension(n,nni) :: diij
    real(FREAL), dimension(n) :: aiij
    real(FREAL), dimension(n,nni) :: d3f
    real(FREAL), dimension(n,ni) :: y
    real(FREAL), dimension(n,ni) :: yd0
    real(FREAL), dimension(n,nni) :: x0
    real(FREAL), dimension(n,ni) :: yd1
    real(FREAL), dimension(n,nni) :: x1
    real(FREAL), dimension(n,ni) :: yd2
    real(FREAL), dimension(n,nni) :: x2
    real(FREAL), dimension(n,ni) :: yd3
    real(FREAL), dimension(n,nni) :: x3
    real(FREAL), dimension(n,ni) :: yd4
    real(FREAL), dimension(n,nni) :: x4
    real(FREAL), dimension(n,4) :: s0
    real(FREAL), dimension(n,4) :: sa
    real(FREAL), dimension(n,4) :: sb
    real(FREAL), dimension(n,4) :: s1
    real(FREAL), dimension(n,4) :: s2
    real(FREAL), dimension(n,4) :: s3
    real(FREAL), dimension(n) :: sb0
    real(FREAL), dimension(n) :: sb1
    real(FREAL), dimension(n,2) :: sb2
    real(FREAL), dimension(n,3) :: sb3
    real(FREAL), dimension(n,nni) :: rjj
    real(FREAL), dimension(n,nni) :: ajj
    real(FREAL), dimension(n,nni) :: djj
    real(FREAL), dimension(n) :: arj
    real(FREAL), dimension(n) :: acj
    !
    !
    integer(FINT) :: j,k,kk,l
    !
    !
    real(FREAL) :: pi,t1,csq3
    !
    ! the number pi
    pi=acos(-one)
    !
    ! sqrt(3)
    csq3=sqrt(3.0_FREAL)
    ! map coordinates
    do j=1,n
       r(j,1)=rint(j,1)
       r(j,2)=rint(j,2)
       r(j,3)=rint(j,3)
       r(j,4)=rint(j,4)
       c(j,1)=rint(j,5)
       c(j,2)=rint(j,6)
       c(j,3)=rint(j,7)
       c(j,4)=rint(j,8)
       c(j,5)=rint(j,9)
       c(j,6)=rint(j,10)
    end do
    !
    do k=1,ni
       do j=1,n
          r2(j,k)=r(j,k)**2
       end do
    end do
    !
    kk=0
    do k=1,ni
       do l=(k+1),ni
          kk=kk+1
          do j=1,n
             d(j,kk) = sqrt(r2(j,k)+r2(j,l)-2.0_FREAL*r(j,k)*r(j,l)*c(j,kk))
          end do
       end do
    end do
    !
    ! evaluate generic one dimensional switching functions
    do k=1,ni
       do j=1,n
          spsw(j,k)=sp(pars_nsw,pars_ssw,r(j,k))
       end do
       !
       do j=1,n
          sprjj(j,k)=sp(pars_nrjj,pars_srjj,r(j,k))
       end do
       do j=1,n
          spdjj(j,k)=sp(pars_ndjj,pars_sdjj,r(j,k))
       end do
       do j=1,n
          spajj(j,k)=sp(pars_najj,pars_sajj,r(j,k))
       end do
    end do
    do k=1,ni
       do l=1,nahb
          do j=1,n
             spab(j,l,k)=spsw(j,k)
          end do
       end do
       do l=1,nad
          do j=1,n
             spad(j,l,k)=spsw(j,k)
          end do
       end do
    end do
    kk=0
    do k=1,ni
       do l=(k+1),ni
          kk=kk+1
          do j=1,n
             spd3(j,kk)=sp(pars_ndjj,pars_sdjj,d(j,kk))
          end do
       end do
    end do
! switch parameters
    call swit4(pars_v0,  spsw,v0j  )
    call swit4(pars_ar,  spsw,arj  )
    call swit4(pars_re,  spsw,req  )
    call swit3(pars_fs,  spsw,fsj  )
    call swit3(pars_as,  spsw,asj  )
    call swit3(pars_e6,  spsw,e6j  )
    call swit3(pars_e8,  spsw,e8j  )
    call swit3(pars_rs,  spsw,rsj  )
    call swip4(pars_ac,  spsw,acj  )
    call swip4(pars_fb0, spsw,fb0j )
    call swip4(pars_fb1, spsw,fb1j )
    call swip4(pars_fb2, spsw,fb2j )
    call swip4(pars_fb3, spsw,fb3j )
    call swip4(pars_aii, spsw,aiij )
    call swip4(pars_ce , spsw,ceq  )
    call swit2(pars_dii, spsw,diij )
    call swip4p(pars_ad0,spad,ad0j ,nad)
    call swip4p(pars_ad1,spad,ad1j ,nad)
    call swip4p(pars_ad2,spad,ad2j ,nad)
    call swip4p(pars_ad3,spad,ad3j ,nad)
    call swip4p(pars_ad4,spad,ad4j ,nad)
    call swip4p(pars_ab ,spab,abj  ,nahb)
    !
    do k=1,ni
       do j=1,n
          sp6(j,k)=sp(6,rsj(j,k),r(j,k))
       end do
       do j=1,n
          sp8(j,k)=sp(8,rsj(j,k),r(j,k))
       end do
    end do
    ! calculate equilibrium structure
    do k=1,ni
       do j=1,n
          spre(j,k)=tanh(arj(j)*(r(j,k)-req(j)))
       end do
    end do
    call swit3(pars_re,spre,rej)
    do k=1,ni
       do j=1,n
          spce(j,k)=tanh(acj(j)*(r(j,k)-req(j)))
       end do
    end do
    call swit2(pars_ce,spce,cej)
    ! prepare coordinates
    do k=1,ni
       do j=1,n
          dr(j,k)=r(j,k)-rej(j,k)
          y(j,k)=(1.0_FREAL-exp(-asj(j,k)*dr(j,k)))/asj(j,k)*(1.0_FREAL+e6j(j,k)*sp6(j,k)+e8j(j,k)*sp8(j,k))
          yd0(j,k)=exp(-ad0j(j,1)*dr(j,k)-ad0j(j,2)*dr(j,k)**2-ad0j(j,3)*dr(j,k)**3)
          yd1(j,k)=exp(-ad1j(j,1)*dr(j,k)-ad1j(j,2)*dr(j,k)**2-ad1j(j,3)*dr(j,k)**3)
          yd2(j,k)=exp(-ad2j(j,1)*dr(j,k)-ad2j(j,2)*dr(j,k)**2-ad2j(j,3)*dr(j,k)**3)
          yd3(j,k)=exp(-ad3j(j,1)*dr(j,k)-ad3j(j,2)*dr(j,k)**2-ad3j(j,3)*dr(j,k)**3)
          yd4(j,k)=exp(-ad4j(j,1)*dr(j,k)-ad4j(j,2)*dr(j,k)**2-ad4j(j,3)*dr(j,k)**3)
       end do
    end do
    kk=0
    do k=1,ni
       do l=(k+1),ni
          kk=kk+1
          do j=1,n
             dc(j,kk)= c(j,kk)-cej(j,kk)
             x0(j,kk)= c(j,kk)*yd0(j,k)*yd0(j,l)
             x1(j,kk)=dc(j,kk)*yd1(j,k)*yd1(j,l)
             x2(j,kk)=dc(j,kk)*yd2(j,k)*yd2(j,l)
             x3(j,kk)=dc(j,kk)*yd3(j,k)*yd3(j,l)
             x4(j,kk)=dc(j,kk)*yd4(j,k)*yd4(j,l)
          end do
       end do
    end do
! symmetrized anharmonic bending coordinates
    do j=1,n
       s0(j,1)=x1(j,1)+x1(j,2)+x1(j,3)+x1(j,4)+x1(j,5)+x1(j,6)
       sa(j,1)=(2*(x1(j,1)+x1(j,6))-(x1(j,2)+x1(j,3)+x1(j,4)+x1(j,5)))/csq3
       sb(j,1)=x1(j,2)-x1(j,3)-x1(j,4)+x1(j,5)
       s1(j,1)=x1(j,1)-x1(j,6)
       s2(j,1)=x1(j,2)-x1(j,5)
       s3(j,1)=x1(j,3)-x1(j,4)
       s0(j,2)=x2(j,1)+x2(j,2)+x2(j,3)+x2(j,4)+x2(j,5)+x2(j,6)
       sa(j,2)=(2*(x2(j,1)+x2(j,6))-(x2(j,2)+x2(j,3)+x2(j,4)+x2(j,5)))/csq3
       sb(j,2)=x2(j,2)-x2(j,3)-x2(j,4)+x2(j,5)
       s1(j,2)=x2(j,1)-x2(j,6)
       s2(j,2)=x2(j,2)-x2(j,5)
       s3(j,2)=x2(j,3)-x2(j,4)
       s0(j,3)=x3(j,1)+x3(j,2)+x3(j,3)+x3(j,4)+x3(j,5)+x3(j,6)
       sa(j,3)=(2*(x3(j,1)+x3(j,6))-(x3(j,2)+x3(j,3)+x3(j,4)+x3(j,5)))/csq3
       sb(j,3)=x3(j,2)-x3(j,3)-x3(j,4)+x3(j,5)
       s1(j,3)=x3(j,1)-x3(j,6)
       s2(j,3)=x3(j,2)-x3(j,5)
       s3(j,3)=x3(j,3)-x3(j,4)
       s0(j,4)=x4(j,1)+x4(j,2)+x4(j,3)+x4(j,4)+x4(j,5)+x4(j,6)
       sa(j,4)=(2*(x4(j,1)+x4(j,6))-(x4(j,2)+x4(j,3)+x4(j,4)+x4(j,5)))/csq3
       sb(j,4)=x4(j,2)-x4(j,3)-x4(j,4)+x4(j,5)
       s1(j,4)=x4(j,1)-x4(j,6)
       s2(j,4)=x4(j,2)-x4(j,5)
       s3(j,4)=x4(j,3)-x4(j,4)
       sb0(j)=x0(j,1)+x0(j,2)+x0(j,3)+x0(j,4)+x0(j,5)+x0(j,6)
       t1 = s0(j,1)+abj(j,1)*s0(j,2)**2+abj(j,2)*(sa(j,2)**2+sb(j,2)**2)+&
            &abj(j,3)*(s1(j,2)**2+s2(j,2)**2+s3(j,2)**2)+abj(j,4)*s0(j,3)**3+ab&
            &j(j,5)*(sa(j,3)**3-3*sa(j,3)*sb(j,3)**2)+abj(j,6)*(sa(j,3)*(2*s1(j&
            &,3)**2-s2(j,3)**2-s3(j,3)**2)+sb(j,3)*csq3*(s2(j,3)**2-s3(j,&
            &3)**2))
       sb1(j) = t1+abj(j,7)*s1(j,3)*s2(j,3)*s3(j,3)+abj(j,8)*s0(j,4)**4+a&
            &bj(j,9)*(sa(j,4)**2+sb(j,4)**2)**2+abj(j,10)*(s1(j,4)**2+s2(j,4)**&
            &2+s3(j,4)**2)*(sa(j,4)**2+sb(j,4)**2)+abj(j,11)*((s1(j,4)**2-s2(j,&
            &4)**2/2-s3(j,4)**2/2)*(sa(j,4)**2-sb(j,4)**2)-csq3*sa(j,4)*s&
            &b(j,4)*(s2(j,4)**2-s3(j,4)**2))+abj(j,12)*(s1(j,4)**4+s2(j,4)**4+s&
            &3(j,4)**4)+abj(j,13)*(s1(j,4)**2*s2(j,4)**2+s1(j,4)**2*s3(j,4)**2+&
            &s2(j,4)**2*s3(j,4)**2)
       t1 = sa(j,1)+abj(j,14)*(sa(j,2)**2-sb(j,2)**2)+abj(j,15)*(s1(j,2)*&
            &*2-s2(j,2)**2/2-s3(j,2)**2/2)+abj(j,16)*(sa(j,3)**2+sb(j,3)**2)*sa&
            &(j,3)+abj(j,17)*(s1(j,3)**2+s2(j,3)**2+s3(j,3)**2)*sa(j,3)+abj(j,1&
            &8)*(sa(j,4)**2-sb(j,4)**2)*(sa(j,4)**2+sb(j,4)**2)
       sb2(j,1) = t1+abj(j,19)*((sa(j,4)**2-sb(j,4)**2)**2-4*sa(j,4)**2*s&
            &b(j,4)**2)+abj(j,20)*(sa(j,4)**2-sb(j,4)**2)*(s1(j,4)**2+s2(j,4)**&
            &2+s3(j,4)**2)+abj(j,21)*(sa(j,4)**2+sb(j,4)**2)*(s1(j,4)**2-s2(j,4&
            &)**2/2-s3(j,4)**2/2)+abj(j,22)*sa(j,4)*s1(j,4)*s2(j,4)*s3(j,4)+abj&
            &(j,23)*(s2(j,4)**2*s3(j,4)**2-s1(j,4)**2*(s2(j,4)**2+s3(j,4)**2)/2&
            &)+abj(j,24)*(s1(j,4)**2-s2(j,4)**2/2-s3(j,4)**2/2)*(s1(j,4)**2+s2(&
            &j,4)**2+s3(j,4)**2)
       sb2(j,2) = sb(j,1)-2*abj(j,14)*sa(j,2)*sb(j,2)+abj(j,15)*csq3&
            &*(s2(j,2)**2-s3(j,2)**2)/2+abj(j,16)*(sa(j,3)**2+sb(j,3)**2)*sb(j&
            &,3)+abj(j,17)*(s1(j,3)**2+s2(j,3)**2+s3(j,3)**2)*sb(j,3)-2*abj(j,1&
            &8)*sa(j,4)*sb(j,4)*(sa(j,4)**2+sb(j,4)**2)+4*abj(j,19)*sa(j,4)*sb(&
            &j,4)*(sa(j,4)**2-sb(j,4)**2)-2*abj(j,20)*sa(j,4)*sb(j,4)*(s1(j,4)*&
            &*2+s2(j,4)**2+s3(j,4)**2)+abj(j,21)*(sa(j,4)**2+sb(j,4)**2)*&
            &csq3*(s2(j,4)**2-s3(j,4)**2)/2+abj(j,22)*sb(j,4)*s1(j,4)*s2(j,4)*s&
            &3(j,4)+abj(j,23)*csq3*s1(j,4)**2*(s3(j,4)**2-s2(j,4)**2)/2+a&
            &bj(j,24)*csq3*(s2(j,4)**2-s3(j,4)**2)*(s1(j,4)**2+s2(j,4)**2&
            &+s3(j,4)**2)/2
       t1 = s1(j,1)+abj(j,25)*s2(j,2)*s3(j,2)+abj(j,26)*sa(j,2)*s1(j,2)+a&
            &bj(j,27)*s1(j,3)**3+abj(j,28)*s1(j,3)*(s2(j,3)**2+s3(j,3)**2)+abj(&
            &j,29)*s1(j,3)*(sa(j,3)**2+sb(j,3)**2)+abj(j,30)*s1(j,3)*(sa(j,3)**&
            &2-sb(j,3)**2)+abj(j,31)*s2(j,3)*s3(j,3)*sa(j,3)
       sb3(j,1) = t1+abj(j,32)*sa(j,4)*s1(j,4)*(sa(j,4)**2+sb(j,4)**2)+ab&
            &j(j,33)*s1(j,4)*(sa(j,4)**3-3*sa(j,4)*sb(j,4)**2)+abj(j,34)*s2(j,4&
            &)*s3(j,4)*(sa(j,4)**2-sb(j,4)**2)+abj(j,35)*s2(j,4)*s3(j,4)*(sa(j,&
            &4)**2+sb(j,4)**2)+abj(j,36)*sa(j,4)*s1(j,4)*(s1(j,4)**2+s2(j,4)**2&
            &+s3(j,4)**2)+abj(j,37)*sa(j,4)*s1(j,4)**3+abj(j,38)*sb(j,4)*s1(j,4&
            &)*(s2(j,4)**2-s3(j,4)**2)+abj(j,39)*s2(j,4)*s3(j,4)*(s1(j,4)**2+s2&
            &(j,4)**2+s3(j,4)**2)+abj(j,40)*s1(j,4)**2*s2(j,4)*s3(j,4)
       t1 = s2(j,1)+abj(j,25)*s1(j,2)*s3(j,2)-abj(j,26)*s2(j,2)*(sa(j,2)/&
            &2-csq3*sb(j,2)/2)+abj(j,27)*s2(j,3)**3+abj(j,28)*s2(j,3)*(s1&
            &(j,3)**2+s3(j,3)**2)+abj(j,29)*s2(j,3)*(sa(j,3)**2+sb(j,3)**2)+abj&
            &(j,30)*s2(j,3)*(sb(j,3)**2/2-sa(j,3)**2/2-csq3*sa(j,3)*sb(j,&
            &3))-abj(j,31)*s3(j,3)*s1(j,3)*(sa(j,3)-csq3*sb(j,3))/2
       sb3(j,2) = t1-abj(j,32)*s2(j,4)*(sa(j,4)/2-csq3*sb(j,4)/2)*(&
            &sa(j,4)**2+sb(j,4)**2)+abj(j,33)*s2(j,4)*(sa(j,4)**3-3*sa(j,4)*sb(&
            &j,4)**2)-abj(j,34)*s1(j,4)*s3(j,4)*(sa(j,4)**2-sb(j,4)**2+2*csq3&
            &*sa(j,4)*sb(j,4))/2+abj(j,35)*s1(j,4)*s3(j,4)*(sa(j,4)**2+sb(j&
            &,4)**2)-abj(j,36)*s2(j,4)*(sa(j,4)/2-csq3*sb(j,4)/2)*(s1(j,4&
            &)**2+s2(j,4)**2+s3(j,4)**2)-abj(j,37)*(sa(j,4)/2-csq3*sb(j,4&
            &)/2)*s2(j,4)**3+abj(j,38)*(csq3*sa(j,4)+sb(j,4))*s2(j,4)*(s1&
            &(j,4)**2-s3(j,4)**2)/2+abj(j,39)*s1(j,4)*s3(j,4)*(s1(j,4)**2+s2(j,&
            &4)**2+s3(j,4)**2)+abj(j,40)*s2(j,4)**2*s1(j,4)*s3(j,4)
       t1 = s3(j,1)+abj(j,25)*s1(j,2)*s2(j,2)-abj(j,26)*s3(j,2)*(sa(j,2)/&
            &2+csq3*sb(j,2)/2)+abj(j,27)*s3(j,3)**3+abj(j,28)*s3(j,3)*(s1&
            &(j,3)**2+s2(j,3)**2)+abj(j,29)*s3(j,3)*(sa(j,3)**2+sb(j,3)**2)+abj&
            &(j,30)*s3(j,3)*(sb(j,3)**2/2-sa(j,3)**2/2+csq3*sa(j,3)*sb(j,&
            &3))-abj(j,31)*s1(j,3)*s2(j,3)*(sa(j,3)+csq3*sb(j,3))/2
       sb3(j,3) = t1-abj(j,32)*s3(j,4)*(sa(j,4)/2+csq3*sb(j,4)/2)*(&
            &sa(j,4)**2+sb(j,4)**2)+abj(j,33)*s3(j,4)*(sa(j,4)**3-3*sa(j,4)*sb(&
            &j,4)**2)-abj(j,34)*s1(j,4)*s2(j,4)*(sa(j,4)**2-sb(j,4)**2-2*csq3&
            &*sa(j,4)*sb(j,4))/2+abj(j,35)*s1(j,4)*s2(j,4)*(sa(j,4)**2+sb(j&
            &,4)**2)-abj(j,36)*s3(j,4)*(sa(j,4)/2+csq3*sb(j,4)/2)*(s1(j,4&
            &)**2+s2(j,4)**2+s3(j,4)**2)-abj(j,37)*(sa(j,4)/2+csq3*sb(j,4&
            &)/2)*s3(j,4)**3+abj(j,38)*(csq3*sa(j,4)-sb(j,4))*s3(j,4)*(s1&
            &(j,4)**2-s2(j,4)**2)/2+abj(j,39)*s1(j,4)*s2(j,4)*(s1(j,4)**2+s2(j,&
            &4)**2+s3(j,4)**2)+abj(j,40)*s3(j,4)**2*s1(j,4)*s2(j,4)
    end do
    ! calculate potential
    do j=1,n
       v(j) = v0j(j)-fb0j(j)*ceq(j)                                     &
            &*(yd0(j,1)*yd0(j,2)+yd0(j,1)*yd0(j,3)+yd0(j,1)*yd0(j,4)+yd0(j,2)*yd0(j,3)+yd0(j,2)*yd0(j,4)+yd0(j,3)*yd0(j,4))
    end do
    ! stretching potential
    do k=1,ni
       do j=1,n
          v(j) = v(j)+0.5_FREAL*fsj(j,k)*y(j,k)**2
       end do
    end do
    ! bending potential
    do j=1,n
       v(j) = v(j)+ fb0j(j)*sb0(j)
       v(j) = v(j)+ fb1j(j)*sb1(j)**2
       v(j) = v(j)+ fb2j(j)*(sb2(j,1)**2+sb2(j,2)**2)
       v(j) = v(j)+ fb3j(j)*(sb3(j,1)**2+sb3(j,2)**2+sb3(j,3)**2)
    end do
    ! asymptotic values for the pair potential
    kk=0
    do k=1,ni
       do l=(k+1),ni
          kk=kk+1
          do j=1,n
             rjj(j,kk)=sqrt(rej(j,k)**2+rej(j,l)**2-2.0_FREAL*rej(j,k)*rej(j,l)*cej(j,kk))*(1.0_FREAL-sprjj(j,k)*sprjj(j,l))&
                  & +pars_riii*(sprjj(j,k)*sprjj(j,l))
             ajj(j,kk)=aiij(j)*(1.0_FREAL-spajj(j,k)*spajj(j,l))+pars_aiii*(spajj(j,k)*spajj(j,l))
             djj(j,kk)=diij(j,kk)*(1.0_FREAL-spdjj(j,k))*(1.0_FREAL-spdjj(j,l))+pars_diii*(spdjj(j,k)*spdjj(j,l))
             dd(j,kk)=d(j,kk)-rjj(j,kk)
          end do
       end do
    end do
    ! pair potential
    do kk=1,nni
       do j=1,n
          v(j)=v(j)+djj(j,kk)*(1.0_FREAL-exp(-ajj(j,kk)*dd(j,kk)))**2
       end do
    end do
! substract diii binding energy at the rate of increase of djj
! to avoid 3-fold clustering, substraction is inhibited whenever
! the distances of more than two yy bonds become too small
    do j=1,n
       d3f(j,1)=spd3(j,2)*spd3(j,3)*spd3(j,4)*spd3(j,5)
       d3f(j,2)=spd3(j,1)*spd3(j,3)*spd3(j,4)*spd3(j,6)
       d3f(j,3)=spd3(j,2)*spd3(j,1)*spd3(j,6)*spd3(j,5)
       v(j)=v(j)-pars_diii*spdjj(j,1)*spdjj(j,2)*d3f(j,1)
       v(j)=v(j)-pars_diii*spdjj(j,1)*spdjj(j,3)*d3f(j,2)
       v(j)=v(j)-pars_diii*spdjj(j,1)*spdjj(j,4)*d3f(j,3)
       v(j)=v(j)-pars_diii*spdjj(j,2)*spdjj(j,3)*d3f(j,3)
       v(j)=v(j)-pars_diii*spdjj(j,2)*spdjj(j,4)*d3f(j,2)
       v(j)=v(j)-pars_diii*spdjj(j,3)*spdjj(j,4)*d3f(j,1)
    end do
    ! energy conversion
    do j=1,n
       v(j)=v(j)*ajoule
    end do
    return
  end subroutine metpot4_vint
!-----------------------------------------------------------------------
  subroutine metpot4_vcar(X,V)
    !
    real(FREAL), dimension(1:N,1:3,1:5), intent(in) :: X
    real(FREAL), dimension(1:N), intent(out) :: V
    real(FREAL), dimension(1:N,1:10) :: R
    real(FREAL), dimension(1:N,1:4) :: R2
    !
    integer(FINT) :: j
    !
    do J=1,N
       R2(J,1)= ((X(J,1,1)-X(J,1,5))**2+(X(J,2,1)-X(J,2,5))**2+(X(J,3,1)-X(J,3,5))**2)
       R2(J,2)= ((X(J,1,2)-X(J,1,5))**2+(X(J,2,2)-X(J,2,5))**2+(X(J,3,2)-X(J,3,5))**2)
       R2(J,3)= ((X(J,1,3)-X(J,1,5))**2+(X(J,2,3)-X(J,2,5))**2+(X(J,3,3)-X(J,3,5))**2)
       R2(J,4)= ((X(J,1,4)-X(J,1,5))**2+(X(J,2,4)-X(J,2,5))**2+(X(J,3,4)-X(J,3,5))**2)
       R(J,1)=sqrt(R2(J,1))
       R(J,2)=sqrt(R2(J,2))
       R(J,3)=sqrt(R2(J,3))
       R(J,4)=sqrt(R2(J,4))
       R(J,5)= (R2(J,1) + R2(J,2)-(X(J,1,1)-X(J,1,2))**2-(X(J,2,1)-X(J,2,2))**2-(X(J,3,1)-X(J,3,2))**2)/(2.0_FREAL*R(J,1)*R(J,2))
       R(J,6)= (R2(J,1) + R2(J,3)-(X(J,1,1)-X(J,1,3))**2-(X(J,2,1)-X(J,2,3))**2-(X(J,3,1)-X(J,3,3))**2)/(2.0_FREAL*R(J,1)*R(J,3))
       R(J,7)= (R2(J,1) + R2(J,4)-(X(J,1,1)-X(J,1,4))**2-(X(J,2,1)-X(J,2,4))**2-(X(J,3,1)-X(J,3,4))**2)/(2.0_FREAL*R(J,1)*R(J,4))
       R(J,8)= (R2(J,3) + R2(J,2)-(X(J,1,3)-X(J,1,2))**2-(X(J,2,3)-X(J,2,2))**2-(X(J,3,3)-X(J,3,2))**2)/(2.0_FREAL*R(J,3)*R(J,2))
       R(J,9)= (R2(J,4) + R2(J,2)-(X(J,1,4)-X(J,1,2))**2-(X(J,2,4)-X(J,2,2))**2-(X(J,3,4)-X(J,3,2))**2)/(2.0_FREAL*R(J,4)*R(J,2))
       R(J,10)= (R2(J,3) + R2(J,4)-(X(J,1,3)-X(J,1,4))**2-(X(J,2,3)-X(J,2,4))**2-(X(J,3,3)-X(J,3,4))**2)/(2.0_FREAL*R(J,3)*R(J,4))
!      R(J,1)=ABOHR*R(J,1)
!      R(J,2)=ABOHR*R(J,2)
!      R(J,3)=ABOHR*R(J,3)
!      R(J,4)=ABOHR*R(J,4)
    end do
    write (*,*) r
    call metpot4_vint(R,V)
    return
  end subroutine metpot4_vcar
!-----------------------------------------------------------------------
  real(FREAL) function sp(ne,s,r)
    !
    integer(FINT), intent(in) :: ne
    real(FREAL), intent(in) :: s,r
    !
    integer(FINT) :: j
    !
    sp=0.0_FREAL
    if(ne.le.0_FINT) stop ' Nsw.LE.0'
    if(R.le.0.0_FREAL) return
    if((s/r)**ne.gt.SWTHRE) then
       sp=0.0_FREAL
    elseif((s/r)**ne.lt.SMALL) then
       sp=1.0_FREAL
    else
       sp=exp(-(s/r)**ne)
    endif
    return
  end function sp
  !-----------------------------------------------------------------------
  subroutine SWIT4(P,SP,PJ)
    !
    real(FREAL), dimension(1:4), intent(in) :: p
    real(FREAL), dimension(1:N,1:4), intent(in) :: sp
    real(FREAL), dimension(1:N), intent(out) :: pj
    real(FREAL), dimension(1:N,1:4) :: sq
    !dimension P(1:4),SP(N,4),PJ(N),SQ(N,4)
    !
    integer(FINT) :: j
    !
    do 1 J=1,N
       SQ(J,1)=1.0_FREAL-SP(J,1)
       SQ(J,2)=1.0_FREAL-SP(J,2)
       SQ(J,3)=1.0_FREAL-SP(J,3)
       SQ(J,4)=1.0_FREAL-SP(J,4)
       PJ(J)=P(1)*(SQ(J,1)*SQ(J,2)*SQ(J,3)*SQ(J,4))                      &
            &     +P(2)*(SP(J,1)*SQ(J,2)*SQ(J,3)*SQ(J,4)+                      &
            &            SQ(J,1)*SP(J,2)*SQ(J,3)*SQ(J,4)+                      &
            &            SQ(J,1)*SQ(J,2)*SP(J,3)*SQ(J,4)+                      &
            &            SQ(J,1)*SQ(J,2)*SQ(J,3)*SP(J,4))                      &
            &     +P(3)*(SP(J,1)*SP(J,2)*SQ(J,3)*SQ(J,4)+                      &
            &            SP(J,1)*SQ(J,2)*SP(J,3)*SQ(J,4)+                      &
            &            SQ(J,1)*SP(J,2)*SP(J,3)*SQ(J,4)+                      &
            &            SP(J,1)*SQ(J,2)*SQ(J,3)*SP(J,4)+                      &
            &            SQ(J,1)*SP(J,2)*SQ(J,3)*SP(J,4)+                      &
            &            SQ(J,1)*SQ(J,2)*SP(J,3)*SP(J,4))                      &
            &     +P(4)*(SQ(J,1)*SP(J,2)*SP(J,3)*SP(J,4)+                      &
            &            SP(J,1)*SQ(J,2)*SP(J,3)*SP(J,4)+                      &
            &            SP(J,1)*SP(J,2)*SQ(J,3)*SP(J,4)+                      &
            &            SP(J,1)*SP(J,2)*SP(J,3)*SQ(J,4))
1   end do
    return
  end subroutine SWIT4
!-----------------------------------------------------------------------
  subroutine SWIP4(P,SP,PJ)
    !implicit double precision (A-H,O-Z)
    real(FREAL), dimension(1:3), intent(in) :: p
    real(FREAL), dimension(1:N,1:4), intent(in) :: sp
    real(FREAL), dimension(1:N), intent(out) :: pj
    real(FREAL), dimension(1:N,1:4) :: sq
    !dimension P(1:3),SP(N,4),PJ(N),SQ(N,4)
    !
    integer(FINT) :: j
    !
    do 1 J=1,N
       SQ(J,1)=1.0_FREAL-SP(J,1)
       SQ(J,2)=1.0_FREAL-SP(J,2)
       SQ(J,3)=1.0_FREAL-SP(J,3)
       SQ(J,4)=1.0_FREAL-SP(J,4)
       PJ(J)=P(1)*(SQ(J,1)*SQ(J,2)*SQ(J,3)*SQ(J,4))                      &
            &     +P(2)*(SP(J,1)*SQ(J,2)*SQ(J,3)*SQ(J,4)+                      &
            &            SQ(J,1)*SP(J,2)*SQ(J,3)*SQ(J,4)+                      &
            &            SQ(J,1)*SQ(J,2)*SP(J,3)*SQ(J,4)+                      &
            &            SQ(J,1)*SQ(J,2)*SQ(J,3)*SP(J,4))                      &
            &     +P(3)*(SP(J,1)*SP(J,2)*SQ(J,3)*SQ(J,4)+                      &
            &            SP(J,1)*SQ(J,2)*SP(J,3)*SQ(J,4)+                      &
            &            SQ(J,1)*SP(J,2)*SP(J,3)*SQ(J,4)+                      &
            &            SP(J,1)*SQ(J,2)*SQ(J,3)*SP(J,4)+                      &
            &            SQ(J,1)*SP(J,2)*SQ(J,3)*SP(J,4)+                      &
            &            SQ(J,1)*SQ(J,2)*SP(J,3)*SP(J,4))
1   end do
    return
  end subroutine SWIP4
!-----------------------------------------------------------------------
  subroutine SWIP4P(p,sp,pj,nl)
    !implicit double precision (A-H,O-Z)
    integer(FINT), parameter :: NNL=100
    integer(FINT), intent(in) :: nl
    real(FREAL), dimension(1:nl,1:3), intent(in) :: p
    real(FREAL), dimension(1:N,1:nl,1:4), intent(in) :: sp
    real(FREAL), dimension(1:N,1:nl), intent(out) :: pj
    real(FREAL), dimension(1:N,1:nl,1:4) :: sq
    !dimension P(1:NL,1:3),SP(N,1:NL,4),PJ(N,1:NL),SQ(N,1:NNL,4)
    !
    integer(FINT) :: j,k
    !
    if(NL.gt.NNL) stop ' RANGE ERROR IN SWIP4P. INTERRUPTED.'
    do K=1,NL
       do J=1,N
          SQ(J,K,1)=1.0_FREAL-SP(J,K,1)
          SQ(J,K,2)=1.0_FREAL-SP(J,K,2)
          SQ(J,K,3)=1.0_FREAL-SP(J,K,3)
          SQ(J,K,4)=1.0_FREAL-SP(J,K,4)
       end do
    end do
    !
    do K=1,NL
       do J=1,N
          PJ(J,K)=P(K,1)*(SQ(J,K,1)*SQ(J,K,2)*SQ(J,K,3)*SQ(J,K,4))          &
               &       +P(K,2)*(SP(J,K,1)*SQ(J,K,2)*SQ(J,K,3)*SQ(J,K,4)+          &
               &                SQ(J,K,1)*SP(J,K,2)*SQ(J,K,3)*SQ(J,K,4)+          &
               &                SQ(J,K,1)*SQ(J,K,2)*SP(J,K,3)*SQ(J,K,4)+          &
               &                SQ(J,K,1)*SQ(J,K,2)*SQ(J,K,3)*SP(J,K,4))          &
               &       +P(K,3)*(SP(J,K,1)*SP(J,K,2)*SQ(J,K,3)*SQ(J,K,4)+          &
               &                SP(J,K,1)*SQ(J,K,2)*SP(J,K,3)*SQ(J,K,4)+          &
               &                SQ(J,K,1)*SP(J,K,2)*SP(J,K,3)*SQ(J,K,4)+          &
               &                SP(J,K,1)*SQ(J,K,2)*SQ(J,K,3)*SP(J,K,4)+          &
               &                SQ(J,K,1)*SP(J,K,2)*SQ(J,K,3)*SP(J,K,4)+          &
               &                SQ(J,K,1)*SQ(J,K,2)*SP(J,K,3)*SP(J,K,4))
       end do
    end do
    !
    return
  end subroutine SWIP4P
!-----------------------------------------------------------------------
  subroutine SWIT3(P,SP,PJ)
    implicit double precision (A-H,O-Z)
    real(FREAL), dimension(1:4), intent(in) :: p
    real(FREAL), dimension(1:N,1:4), intent(in) :: sp
    real(FREAL), dimension(1:N,1:4), intent(out) :: pj
    real(FREAL), dimension(1:N,1:4) :: sq
    !dimension P(1:4),SP(N,4),PJ(N,4),SQ(N,4)
    !
    integer(FINT) :: j,i,k,l,m
    !
    do 1 J=1,N
       SQ(J,1)=1.0_FREAL-SP(J,1)
       SQ(J,2)=1.0_FREAL-SP(J,2)
       SQ(J,3)=1.0_FREAL-SP(J,3)
       SQ(J,4)=1.0_FREAL-SP(J,4)
1   end do
    do 6 I=1,4
       do 5 K=1,4
          if(K.eq.I) goto 5
          do 4 L=(K+1),4
             if(L.eq.I) goto 4
             do 3 M=(L+1),4
                if(M.eq.I) goto 3
                do 2 J=1,N
                   PJ(J,I)=P(1)*(SQ(J,L)*SQ(J,M)*SQ(J,K))                            &
                        &       +P(2)*(SP(J,L)*SQ(J,M)*SQ(J,K)+                            &
                        &              SQ(J,L)*SP(J,M)*SQ(J,K)+                            &
                        &              SQ(J,L)*SQ(J,M)*SP(J,K))                            &
                        &       +P(3)*(SP(J,L)*SP(J,M)*SQ(J,K)+                            &
                        &              SP(J,L)*SQ(J,M)*SP(J,K)+                            &
                        &              SQ(J,L)*SP(J,M)*SP(J,K))                            &
                        &       +P(4)*(SP(J,L)*SP(J,M)*SP(J,K))
2               end do
3            end do
4         end do
5      end do
6   end do
    return
  end subroutine SWIT3
!-----------------------------------------------------------------------
  subroutine SWIT2(P,SP,PJ)
    implicit double precision (A-H,O-Z)
    real(FREAL), dimension(1:3), intent(in) :: p
    real(FREAL), dimension(1:N,1:4), intent(in) :: sp
    real(FREAL), dimension(1:N,1:6), intent(out) :: pj
    real(FREAL), dimension(1:N,1:4) :: sq
    !dimension P(1:3),SP(N,4),PJ(N,6),SQ(N,4)
    !
    integer(FINT) :: j,kk,i,k,l,m
    !
    do 1 J=1,N
       SQ(J,1)=1.0_FREAL-SP(J,1)
       SQ(J,2)=1.0_FREAL-SP(J,2)
       SQ(J,3)=1.0_FREAL-SP(J,3)
       SQ(J,4)=1.0_FREAL-SP(J,4)
1   end do
    KK=0
    do 6 I=1,4
       do 5 K=(I+1),4
          KK=KK+1
          do 4 L=1,4
             if(L.eq.I.or.L.eq.K) goto 4
             do 3 M=(L+1),4
                if(M.eq.I.or.M.eq.K) goto 3
                do 2 J=1,N
                   PJ(J,KK)=P(1)*(SQ(J,L)*SQ(J,M))                                   &
                        &        +P(2)*(SP(J,L)*SQ(J,M)+                                   &
                        &               SQ(J,L)*SP(J,M))                                   &
                        &        +P(3)*(SP(J,L)*SP(J,M))
2               end do
3            end do
4         end do
5      end do
6   end do
    return
  end subroutine SWIT2
!-----------------------------------------------------------------------
  real(FREAL) function tanh(x)
    !implicit double precision(A-H,O-Z)
    real(FREAL), parameter :: TOL=100.0_FREAL
    real(FREAL), intent(in) :: x
    !
    real(FREAL) :: ep
    !
    if(x.ge.0.0_FREAL) then
       tanh=1.0_FREAL
       if(x.lt.tol) then
          ep=exp(-x)
          tanh=(1.0_FREAL-ep**2)/(1.0_FREAL+ep**2)
       endif
    else
       tanh=-1.0_FREAL
       if(x.gt.-tol) then
          ep=exp(x)
          tanh=-(1.0_FREAL-ep**2)/(1.0_FREAL+ep**2)
       endif
    endif
    return
  end function tanh
!-----------------------------------------------------------------------
!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H subroutine metpot4_initpot()
!H-----------------------------------------------------------------------------
!H
!
  subroutine metpot4_initpot()
    !
    !
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the metpot4 potential: default values.
!H
!H Old comment: ! definition of molecule specific potential parameters
!H
!H-----------------------------------------------------------------------------
!H
!
! definition of molecule specific potential parameters
    !parameter(NI=4,MI=NI-1,NAHB=40,NAD=3)
! nominal parameters
! units are aJ=10**(-18)J for energies and A=100pm for lengths
!!      common /PAR0/   V0(1:NI)
!!      common /PARS1/  AR(1:NI),RE(1:NI)
!!      common /PARS2/  FS(1:NI)
!!      common /PARS3/  AS(1:NI)
!!      common /PARS4/  E6(1:NI)
!!      common /PARS5/  E8(1:NI)
!!      common /PARS6/  RS(1:NI)
!!      common /PARB1/  AC(1:MI),CE(1:MI)
!!      common /PARB2/  FB0(1:MI)
!!      common /PARB3/  FB1(1:MI)
!!      common /PARB4/  FB2(1:MI)
!!      common /PARB5/  FB3(1:MI)
!!      common /PARB6/  AD0(1:NAD,1:MI)
!!      common /PARB8/  AD1(1:NAD,1:MI)
!!      common /PARB9/  AD2(1:NAD,1:MI)
!!      common /PARB10/ AD3(1:NAD,1:MI)
!!      common /PARB11/ AD4(1:NAD,1:MI)
!!      common /PARB12/ AB(1:NAHB,1:MI)
!!      common /PARI1/  DII(1:MI)
!!      common /PARI2/  AII(1:MI)
!!      common /PARI3/  NRJJ,SRJJ,RIII
!!      common /PARI4/  NAJJ,SAJJ,AIII
!!      common /PARI5/  NDJJ,SDJJ,DIII
!!      common /PARSW/  NSW,SSW
    !
    ! Allocate variables
    !-------------------
    ! we check for init so that we can simply call _initpot
    ! again to reset the potential
    !------------------------------------------------------
    if (.not.LOCmetpot4_isinit) then
       if(allocated(pars_v0)) then
          call message(MESWARN,"[initpot]: v0 already alllocated.")
          deallocate(pars_v0,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate v0.")
             stop 3
          endif
       end if
       if(allocated(pars_ar)) then
          call message(MESWARN,"[initpot]: ar already alllocated.")
          deallocate(pars_ar,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate ar.")
             stop 3
          endif
       end if
       if(allocated(pars_re)) then
          call message(MESWARN,"[initpot]: re already alllocated.")
          deallocate(pars_re,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate re.")
             stop 3
          endif
       end if
       if(allocated(pars_fs)) then
          call message(MESWARN,"[initpot]: fs already alllocated.")
          deallocate(pars_fs,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate fs.")
             stop 3
          endif
       end if
       if(allocated(pars_as)) then
          call message(MESWARN,"[initpot]: as already alllocated.")
          deallocate(pars_as,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate as.")
             stop 3
          endif
       end if
       if(allocated(pars_e6)) then
          call message(MESWARN,"[initpot]: e6 already alllocated.")
          deallocate(pars_e6,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate e6.")
             stop 3
          endif
       end if
       if(allocated(pars_e8)) then
          call message(MESWARN,"[initpot]: e8 already alllocated.")
          deallocate(pars_e8,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate e8.")
             stop 3
          endif
       end if
       if(allocated(pars_rs)) then
          call message(MESWARN,"[initpot]: rs already alllocated.")
          deallocate(pars_rs,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate rs.")
             stop 3
          endif
       end if
       if(allocated(pars_ac)) then
          call message(MESWARN,"[initpot]: ac already alllocated.")
          deallocate(pars_ac,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate ac.")
             stop 3
          endif
       end if
       if(allocated(pars_ce)) then
          call message(MESWARN,"[initpot]: ce already alllocated.")
          deallocate(pars_ce,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate ce.")
             stop 3
          endif
       end if
       if(allocated(pars_fb0)) then
          call message(MESWARN,"[initpot]: fb0 already alllocated.")
          deallocate(pars_fb0,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate fb0.")
             stop 3
          endif
       end if
       if(allocated(pars_fb1)) then
          call message(MESWARN,"[initpot]: fb1 already alllocated.")
          deallocate(pars_fb1,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate fb1.")
             stop 3
          endif
       end if
       if(allocated(pars_fb2)) then
          call message(MESWARN,"[initpot]: fb2 already alllocated.")
          deallocate(pars_fb2,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate fb2.")
             stop 3
          endif
       end if
       if(allocated(pars_fb3)) then
          call message(MESWARN,"[initpot]: fb3 already alllocated.")
          deallocate(pars_fb3,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate fb3.")
             stop 3
          endif
       end if
       if(allocated(pars_ad0)) then
          call message(MESWARN,"[initpot]: ad0 already alllocated.")
          deallocate(pars_ad0,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate ad0.")
             stop 3
          endif
       end if
       if(allocated(pars_ad1)) then
          call message(MESWARN,"[initpot]: ad1 already alllocated.")
          deallocate(pars_ad1,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate ad1.")
             stop 3
          endif
       end if
       if(allocated(pars_ad2)) then
          call message(MESWARN,"[initpot]: ad2 already alllocated.")
          deallocate(pars_ad2,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate ad2.")
             stop 3
          endif
       end if
       if(allocated(pars_ad3)) then
          call message(MESWARN,"[initpot]: ad3 already alllocated.")
          deallocate(pars_ad3,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate ad3.")
             stop 3
          endif
       end if
       if(allocated(pars_ad4)) then
          call message(MESWARN,"[initpot]: ad4 already alllocated.")
          deallocate(pars_ad4,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate ad4.")
             stop 3
          endif
       end if
       if(allocated(pars_ab)) then
          call message(MESWARN,"[initpot]: ab already alllocated.")
          deallocate(pars_ab,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate ab.")
             stop 3
          endif
       end if
       if(allocated(pars_dii)) then
          call message(MESWARN,"[initpot]: dii already alllocated.")
          deallocate(pars_dii,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate dii.")
             stop 3
          endif
       end if
       if(allocated(pars_aii)) then
          call message(MESWARN,"[initpot]: aii already alllocated.")
          deallocate(pars_aii,STAT=irc)
          if(irc.ne.0) then
             call message(MESERRO,"[initpot]: unable to deallocate aii.")
             stop 3
          endif
       end if
       !
       allocate(pars_v0(1:ni),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate v0.")
          stop 3
       endif

       allocate(pars_ar(1:ni),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate ar.")
          stop 3
       endif

       allocate(pars_re(1:ni),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate re.")
          stop 3
       endif

       allocate(pars_fs(1:ni),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate fs.")
          stop 3
       endif

       allocate(pars_as(1:ni),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate as.")
          stop 3
       endif

       allocate(pars_e6(1:ni),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate e6.")
          stop 3
       endif

       allocate(pars_e8(1:ni),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate e8.")
          stop 3
       endif

       allocate(pars_rs(1:ni),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate rs.")
          stop 3
       endif

       allocate(pars_ac(1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate ac.")
          stop 3
       endif

       allocate(pars_ce(1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate ce.")
          stop 3
       endif

       allocate(pars_fb0(1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate fb0.")
          stop 3
       endif

       allocate(pars_fb1(1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate fb1.")
          stop 3
       endif

       allocate(pars_fb2(1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate fb2.")
          stop 3
       endif

       allocate(pars_fb3(1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate fb3.")
          stop 3
       endif

       allocate(pars_ad0(1:nad,1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate ad0.")
          stop 3
       endif

       allocate(pars_ad1(1:nad,1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate ad1.")
          stop 3
       endif

       allocate(pars_ad2(1:nad,1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate ad2.")
          stop 3
       endif

       allocate(pars_ad3(1:nad,1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate ad3.")
          stop 3
       endif

       allocate(pars_ad4(1:nad,1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate ad4.")
          stop 3
       endif

       allocate(pars_ab(1:nahb,1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate ab.")
          stop 3
       endif

       allocate(pars_dii(1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate dii.")
          stop 3
       endif

       allocate(pars_aii(1:mi),STAT=irc)
       if(irc.ne.0) then
          call message(MESERRO,"[initpot]: unable to allocate aii.")
          stop 3
       endif
    end if
    !
    ! parameters for the potential
    ! optionally parameters may be read from an external file
    ! by assigning with
    !     READ(*,*) V0(1)
    ! etc, which was not done here
    pars_v0(1)        =  0.000000000e+00_FREAL
    pars_ar(1)        =  0.138066204_FREAL
    pars_re(1)        =   1.08580000_FREAL
    pars_fs(1)        =   5.36719376_FREAL
    pars_as(1)        =   1.76960545_FREAL
    pars_e6(1)        =  0.576357501e-03_FREAL
    pars_e8(1)        = -0.437869154e-01_FREAL
    pars_rs(1)        =   3.00000000_FREAL
    pars_ac(1)        =  0.912461600_FREAL
    pars_ce(1)        = -0.333333300_FREAL
    pars_fb0(1)       =  0.146835789_FREAL
    pars_fb1(1)       =  0.196780787e-01_FREAL
    pars_fb2(1)       =  0.820392218e-01_FREAL
    pars_fb3(1)       =  0.988915918e-01_FREAL
    pars_ad0(1,1)     =  0.712116053_FREAL
    pars_ad1(1,1)     =  0.250514171_FREAL
    pars_ad2(1,1)     =  0.439195737_FREAL
    pars_ad3(1,1)     = -0.288854799_FREAL
    pars_ad4(1,1)     =  0.000000000e+00_FREAL
    pars_ad0(2,1)     =   1.00101245_FREAL
    pars_ad1(2,1)     = -0.484864657e-01_FREAL
    pars_ad2(2,1)     =  0.315480471_FREAL
    pars_ad3(2,1)     =  0.315480471_FREAL
    pars_ad4(2,1)     =  0.000000000e+00_FREAL
    pars_ad0(3,1)     =  0.000000000e+00_FREAL
    pars_ad1(3,1)     =  0.161650526_FREAL
    pars_ad2(3,1)     =  0.161650526_FREAL
    pars_ad3(3,1)     =  0.161650526_FREAL
    pars_ad4(3,1)     =  0.000000000e+00_FREAL
    pars_ab( 1,1)     =  0.000000000e+00_FREAL
    pars_ab( 2,1)     =  0.227772946_FREAL
    pars_ab( 3,1)     = -0.590449983_FREAL
    pars_ab( 4,1)     =  0.000000000e+00_FREAL
    pars_ab( 5,1)     = -0.749425140e-01_FREAL
    pars_ab( 6,1)     =  0.158542242_FREAL
    pars_ab( 7,1)     = -0.899683975e-01_FREAL
    pars_ab( 8,1)     =  0.000000000e+00_FREAL
    pars_ab( 9,1)     =  0.000000000e+00_FREAL
    pars_ab(10,1)     =  0.000000000e+00_FREAL
    pars_ab(11,1)     =  0.000000000e+00_FREAL
    pars_ab(12,1)     =  0.000000000e+00_FREAL
    pars_ab(13,1)     =  0.000000000e+00_FREAL
    pars_ab(14,1)     = -0.947988981e-01_FREAL
    pars_ab(15,1)     =  0.379195592_FREAL
    pars_ab(16,1)     =  0.193537976e-01_FREAL
    pars_ab(17,1)     = -0.243975134_FREAL
    pars_ab(18,1)     =  0.000000000e+00_FREAL
    pars_ab(19,1)     =  0.000000000e+00_FREAL
    pars_ab(20,1)     =  0.000000000e+00_FREAL
    pars_ab(21,1)     =  0.000000000e+00_FREAL
    pars_ab(22,1)     =  0.000000000e+00_FREAL
    pars_ab(23,1)     =  0.000000000e+00_FREAL
    pars_ab(24,1)     =  0.000000000e+00_FREAL
    pars_ab(25,1)     =  0.635332915e-01_FREAL
    pars_ab(26,1)     = -0.366809636e-01_FREAL
    pars_ab(27,1)     =  0.336515626e-01_FREAL
    pars_ab(28,1)     =  0.728605440e-01_FREAL
    pars_ab(29,1)     =  0.219828219_FREAL
    pars_ab(30,1)     =  0.130696605e-01_FREAL
    pars_ab(31,1)     = -0.226373155e-01_FREAL
    pars_ab(32,1)     =  0.000000000e+00_FREAL
    pars_ab(33,1)     =  0.000000000e+00_FREAL
    pars_ab(34,1)     =  0.000000000e+00_FREAL
    pars_ab(35,1)     =  0.000000000e+00_FREAL
    pars_ab(36,1)     =  0.000000000e+00_FREAL
    pars_ab(37,1)     =  0.000000000e+00_FREAL
    pars_ab(38,1)     =  0.000000000e+00_FREAL
    pars_ab(39,1)     =  0.000000000e+00_FREAL
    pars_ab(40,1)     =  0.000000000e+00_FREAL
    pars_dii(1)       =  0.594584559e-04_FREAL
    pars_aii(1)       =   3.60000000_FREAL
    pars_v0(2)        =  0.000000000e+00_FREAL
    pars_ar(2)        =   1.42153110_FREAL
    pars_re(2)        =   1.07900000_FREAL
    pars_fs(2)        =   5.62068862_FREAL
    pars_as(2)        =   1.76960545_FREAL
    pars_e6(2)        =  0.576357501e-03_FREAL
    pars_e8(2)        = -0.587475070e-01_FREAL
    pars_rs(2)        =   3.00000000_FREAL
    pars_ac(2)        =  0.912461600_FREAL
    pars_ce(2)        = -0.500000000_FREAL
    pars_fb0(2)       =  0.545632984e-01_FREAL
    pars_fb1(2)       =  0.546400258e-01_FREAL
    pars_fb2(2)       =  0.116832553_FREAL
    pars_fb3(2)       =  0.000000000e+00_FREAL
    pars_ad0(1,2)     =  0.000000000e+00_FREAL
    pars_ad1(1,2)     =  0.250514171_FREAL
    pars_ad2(1,2)     =  0.439195737_FREAL
    pars_ad3(1,2)     = -0.288854799_FREAL
    pars_ad4(1,2)     =  0.000000000e+00_FREAL
    pars_ad0(2,2)     =  0.800000000_FREAL
    pars_ad1(2,2)     = -0.484864657e-01_FREAL
    pars_ad2(2,2)     =  0.315480471_FREAL
    pars_ad3(2,2)     =  0.315480471_FREAL
    pars_ad4(2,2)     =  0.000000000e+00_FREAL
    pars_ad0(3,2)     =  0.000000000e+00_FREAL
    pars_ad1(3,2)     =  0.161650526_FREAL
    pars_ad2(3,2)     =  0.161650526_FREAL
    pars_ad3(3,2)     =  0.161650526_FREAL
    pars_ad4(3,2)     =  0.000000000e+00_FREAL
    pars_ab( 1,2)     =  0.000000000e+00_FREAL
    pars_ab( 2,2)     =  0.000000000e+00_FREAL
    pars_ab( 3,2)     =  0.000000000e+00_FREAL
    pars_ab( 4,2)     =  0.000000000e+00_FREAL
    pars_ab( 5,2)     =  0.000000000e+00_FREAL
    pars_ab( 6,2)     =  0.000000000e+00_FREAL
    pars_ab( 7,2)     =  0.000000000e+00_FREAL
    pars_ab( 8,2)     =  0.000000000e+00_FREAL
    pars_ab( 9,2)     =  0.000000000e+00_FREAL
    pars_ab(10,2)     =  0.000000000e+00_FREAL
    pars_ab(11,2)     =  0.000000000e+00_FREAL
    pars_ab(12,2)     =  0.000000000e+00_FREAL
    pars_ab(13,2)     =  0.000000000e+00_FREAL
    pars_ab(14,2)     =  0.000000000e+00_FREAL
    pars_ab(15,2)     =  0.000000000e+00_FREAL
    pars_ab(16,2)     =  0.000000000e+00_FREAL
    pars_ab(17,2)     =  0.000000000e+00_FREAL
    pars_ab(18,2)     =  0.000000000e+00_FREAL
    pars_ab(19,2)     =  0.000000000e+00_FREAL
    pars_ab(20,2)     =  0.000000000e+00_FREAL
    pars_ab(21,2)     =  0.000000000e+00_FREAL
    pars_ab(22,2)     =  0.000000000e+00_FREAL
    pars_ab(23,2)     =  0.000000000e+00_FREAL
    pars_ab(24,2)     =  0.000000000e+00_FREAL
    pars_ab(25,2)     =  0.000000000e+00_FREAL
    pars_ab(26,2)     =  0.000000000e+00_FREAL
    pars_ab(27,2)     =  0.000000000e+00_FREAL
    pars_ab(28,2)     =  0.000000000e+00_FREAL
    pars_ab(29,2)     =  0.000000000e+00_FREAL
    pars_ab(30,2)     =  0.000000000e+00_FREAL
    pars_ab(31,2)     =  0.000000000e+00_FREAL
    pars_ab(32,2)     =  0.000000000e+00_FREAL
    pars_ab(33,2)     =  0.000000000e+00_FREAL
    pars_ab(34,2)     =  0.000000000e+00_FREAL
    pars_ab(35,2)     =  0.000000000e+00_FREAL
    pars_ab(36,2)     =  0.000000000e+00_FREAL
    pars_ab(37,2)     =  0.000000000e+00_FREAL
    pars_ab(38,2)     =  0.000000000e+00_FREAL
    pars_ab(39,2)     =  0.000000000e+00_FREAL
    pars_ab(40,2)     =  0.000000000e+00_FREAL
    pars_dii(2)       =  0.108273496e-01_FREAL
    pars_aii(2)       =   2.00000000_FREAL
    pars_v0(3)        =  0.000000000e+00_FREAL
    pars_ar(3)        =  0.568493183_FREAL
    pars_re(3)        =   1.07500000_FREAL
    pars_fs(3)        =   5.11607882_FREAL
    pars_as(3)        =   1.76960545_FREAL
    pars_e6(3)        =  0.576357501e-03_FREAL
    pars_e8(3)        = -0.300000000e-01_FREAL
    pars_rs(3)        =   3.00000000_FREAL
    pars_ac(3)        =  0.912461600_FREAL
    pars_ce(3)        = -0.692646889_FREAL
    pars_fb0(3)       =  0.000000000e+00_FREAL
    pars_fb1(3)       =  0.281148286_FREAL
    pars_fb2(3)       =  0.000000000e+00_FREAL
    pars_fb3(3)       =  0.000000000e+00_FREAL
    pars_ad0(1,3)     =  0.000000000e+00_FREAL
    pars_ad1(1,3)     =  0.250514171_FREAL
    pars_ad2(1,3)     =  0.439195737_FREAL
    pars_ad3(1,3)     = -0.288854799_FREAL
    pars_ad4(1,3)     =  0.000000000e+00_FREAL
    pars_ad0(2,3)     =  0.800000000_FREAL
    pars_ad1(2,3)     = -0.484864657e-01_FREAL
    pars_ad2(2,3)     =  0.315480471_FREAL
    pars_ad3(2,3)     =  0.315480471_FREAL
    pars_ad4(2,3)     =  0.000000000e+00_FREAL
    pars_ad0(3,3)     =  0.000000000e+00_FREAL
    pars_ad1(3,3)     =  0.161650526_FREAL
    pars_ad2(3,3)     =  0.161650526_FREAL
    pars_ad3(3,3)     =  0.161650526_FREAL
    pars_ad4(3,3)     =  0.000000000e+00_FREAL
    pars_ab( 1,3)     = -0.600000000_FREAL
    pars_ab( 2,3)     =  0.000000000e+00_FREAL
    pars_ab( 3,3)     =  0.000000000e+00_FREAL
    pars_ab( 4,3)     = -0.290000000_FREAL
    pars_ab( 5,3)     =  0.000000000e+00_FREAL
    pars_ab( 6,3)     =  0.000000000e+00_FREAL
    pars_ab( 7,3)     =  0.000000000e+00_FREAL
    pars_ab( 8,3)     =  0.000000000e+00_FREAL
    pars_ab( 9,3)     =  0.000000000e+00_FREAL
    pars_ab(10,3)     =  0.000000000e+00_FREAL
    pars_ab(11,3)     =  0.000000000e+00_FREAL
    pars_ab(12,3)     =  0.000000000e+00_FREAL
    pars_ab(13,3)     =  0.000000000e+00_FREAL
    pars_ab(14,3)     =  0.000000000e+00_FREAL
    pars_ab(15,3)     =  0.000000000e+00_FREAL
    pars_ab(16,3)     =  0.000000000e+00_FREAL
    pars_ab(17,3)     =  0.000000000e+00_FREAL
    pars_ab(18,3)     =  0.000000000e+00_FREAL
    pars_ab(19,3)     =  0.000000000e+00_FREAL
    pars_ab(20,3)     =  0.000000000e+00_FREAL
    pars_ab(21,3)     =  0.000000000e+00_FREAL
    pars_ab(22,3)     =  0.000000000e+00_FREAL
    pars_ab(23,3)     =  0.000000000e+00_FREAL
    pars_ab(24,3)     =  0.000000000e+00_FREAL
    pars_ab(25,3)     =  0.000000000e+00_FREAL
    pars_ab(26,3)     =  0.000000000e+00_FREAL
    pars_ab(27,3)     =  0.000000000e+00_FREAL
    pars_ab(28,3)     =  0.000000000e+00_FREAL
    pars_ab(29,3)     =  0.000000000e+00_FREAL
    pars_ab(30,3)     =  0.000000000e+00_FREAL
    pars_ab(31,3)     =  0.000000000e+00_FREAL
    pars_ab(32,3)     =  0.000000000e+00_FREAL
    pars_ab(33,3)     =  0.000000000e+00_FREAL
    pars_ab(34,3)     =  0.000000000e+00_FREAL
    pars_ab(35,3)     =  0.000000000e+00_FREAL
    pars_ab(36,3)     =  0.000000000e+00_FREAL
    pars_ab(37,3)     =  0.000000000e+00_FREAL
    pars_ab(38,3)     =  0.000000000e+00_FREAL
    pars_ab(39,3)     =  0.000000000e+00_FREAL
    pars_ab(40,3)     =  0.000000000e+00_FREAL
    pars_dii(3)       =  0.483034861e-01_FREAL
    pars_aii(3)       =   2.00000000_FREAL
    pars_v0(4)        =  0.000000000e+00_FREAL
    pars_ar(4)        =  0.568493183_FREAL
    pars_re(4)        =   1.11980000_FREAL
    pars_fs(4)        =   4.48563674_FREAL
    pars_as(4)        =  1.800000000e+00_FREAL
    pars_e6(4)        =  0.576357501e-03_FREAL
    pars_e8(4)        = -0.300000000e-01_FREAL
    pars_rs(4)        =   3.00000000_FREAL
    pars_nsw          =   8.00000000_FREAL
    pars_ssw          =   3.00000000_FREAL
    pars_nrjj         =   12.0000000_FREAL
    pars_srjj         =   1.25742563_FREAL
    pars_riii         =  0.741400000_FREAL
    pars_najj         =   4.00000000_FREAL
    pars_sajj         =   1.80000000_FREAL
    pars_aiii         =   1.94534000_FREAL
    pars_ndjj         =   6.00000000_FREAL
    pars_sdjj         =   1.66154170_FREAL
    pars_diii         =  0.760596828_FREAL
    return
  end subroutine metpot4_initpot
    !
    !
    !
  end module metpot4
!
