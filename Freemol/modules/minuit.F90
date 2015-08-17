!H
!H-----------------------------------------------------------------------------
!H This File is part of the Freemol package
!H This is a fortran 90 port of the original (1999) F77 minuit downloaded
!H from http://www.cern.ch/. See there the current copyright statment.
!H-----------------------------------------------------------------------------
!H MODULE MINUIT for Freemol by F.Mariotti
!H-----------------------------------------------------------------------------
!H $Id: minuit.F90,v 1.1.1.1 2009/01/12 16:56:17 mariotti Exp $
!H-----------------------------------------------------------------------------
!H
!H DESCRIPTION:
!H This is the cern package minuit ported to fortran90 for the Freemol package
!H and it is quite in beta release!!!! Please check for minuit at the
!H CERN web page: http://www.cern.ch/
!H
!H TODO:
!H
!H
!H WARN:
!H
!H HERE SOME LINES FROM THE ORIGINAL VERSION:
!H     *
!H     * Revision 1.1  2002/07/19 01:44:03  fmariot
!H     * Just added minuit files but to be configured
!H     *
!H     * Revision 1.1.1.1  1996/03/07 14:31:28  mclareni
!H     * Minuit
!H     *
!H     *
!H
!H HERE THE OLD DOCS:
!H
!H     CPNAM   Parameter name (10 characters)
!H     U       External (visible to user in FCN) value of parameter
!H     ALIM, BLIM Lower and upper parameter limits. If both zero, no limi
!H     ERP,ERN Positive and negative MINOS errors, if calculated.
!H     WERR    External parameter error (standard deviation, defined by U
!H     GLOBCC  Global Correlation Coefficient
!H     NVARL   =-1 if parameter undefined,      =0 if constant,
!H     = 1 if variable without limits,  =4 if variable with limits
!H     (Note that if parameter has been fixed, NVARL=1 or =4, and NIOFEX=
!H     NIOFEX  Internal parameter number, or zero if not currently variab
!H     NEXOFI  External parameter number for currently variable parameter
!H     X, XT   Internal parameter values (X are sometimes saved in XT)
!H     DIRIN   (Internal) step sizes for current step
!H     variables with names ending in ..S are saved values for fixed para
!H     VHMAT   (Internal) error matrix stored as Half MATrix, since
!H     it is symmetric
!H     VTHMAT  VHMAT is sometimes saved in VTHMAT, especially in MNMNOT
!H
!H     ISW definitions:
!H      ISW(1) =0 normally, =1 means CALL LIMIT EXCEEDED
!H      ISW(2) =0 means no error matrix
!H             =1 means only approximate error matrix
!H             =2 means full error matrix, but forced pos-def.
!H             =3 means good normal full error matrix exists
!H      ISW(3) =0 if Minuit is calculating the first derivatives
!H             =1 if first derivatives calculated inside FCN
!H      ISW(4) =-1 if most recent minimization did not converge.
!H             = 0 if problem redefined since most recent minimization.
!H             =+1 if most recent minimization did converge.
!H      ISW(5) is the PRInt level.  See SHO PRIntlevel
!H      ISW(6) = 0 for batch mode, =1 for interactive mode
!H                      =-1 for originally interactive temporarily batch
!H
!H     LWARN is true if warning messges are to be put out (default=true)
!H            SET WARN turns it on, set NOWarn turns it off
!H     LREPOR is true if exceptional conditions are put out (default=fals
!H            SET DEBUG turns it on, SET NODebug turns it off
!H     LIMSET is true if a parameter is up against limits (for MINOS)
!H     LNOLIM is true if there are no limits on any parameters (not yet u
!H     LNEWMN is true if the previous process has unexpectedly improved F
!H     LPHEAD is true if a heading should be put out for the next paramet
!H        definition, false if a parameter has just been defined
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H
!
module minuit
  !
  use vartypes
  use messages
  use dummy_minuit_function
  !
  implicit none
  private
  !
  ! PUBLIC DECLARATIONS:
  !--------------------

  ! Public Procedures
  !------------------
  public minuit_init  ! Ex MNINIT
  public minuit_isinit
  public minuit_geterror
  !public minuit_
  public minuit_mnparm
  public minuit_mnseti
  public minuit_mnexcm
  public minuit_mnemat
  public minuit_mncont
  public minuit_mnpout
  public minuit_minuit
  !
  ! Some Minuit Subroutines are not here .... like:
  !------------------------------------------------
  ! mncomd
  ! mnstat
  ! mnerrs
  ! mnintr
  ! mninpu
  !
  !
  ! Error Interface
  !----------------
  ! We define two error interfaces. Each one return the latest error
  ! code or string value. The String error function can be used to inqure
  ! about a given error code.
  !----------------------------------------------------------------------
  interface minuit_geterror
     module procedure minuit_geterror_i
     module procedure minuit_geterror_c
  end interface !minuit_geterror

  ! LOCAL DECLARATIONS:
  !--------------------
  !
  logical, save :: mnit_isinit = .false.
  !
  ! error variables
  !----------------
  integer(FINT), save, private :: mnit_error = 0
  character(FLCHARS), save, private :: mnit_error_message = ''
  !
  !Minuit I/O units
  !----------------
  integer(FINT), save :: minuit_filein !OLD:isysrd
  integer(FINT), save :: minuit_fileout !OLD:isyswr
  integer(FINT), save :: minuit_filesave !OLD:isyssa
  !
  ! Others
  !-------
  integer(FINT), parameter :: mne=100
  integer(FINT), parameter :: mni=50
  integer(FINT), parameter :: mnihl=mni*(mni+1)/2
  integer(FINT), parameter :: maxdbg=10
  integer(FINT), parameter :: maxstk=10
  integer(FINT), parameter :: maxcwd=20
  integer(FINT), parameter :: maxp=30
  integer(FINT), parameter :: maxcpt=101
  real(FREAL), parameter :: zero=0.0, one=1.0, half=0.5
  !
  integer(FINT), save :: icomnd, istrat, itaur
  !integer(FINT), save :: ke1cr, ke2cr, limset, lnewmn, lnolim, lphead
  !integer(FINT), save :: lrepor, lwarn, maxext, maxint, nblock, newpag
  integer(FINT), save :: ke1cr, ke2cr
  integer(FINT), save :: maxext, maxint, nblock, newpag
  logical, save :: lrepor, lwarn, limset, lnewmn, lnolim, lphead
  integer(FINT), save :: nfcn, nfcnfr, nfcnlc, nfcnmx, npagln, npagwd
  integer(FINT), save :: npar, npfix, nstkrd, nstkwr, nu
  !
  integer(FINT), save, dimension(0:maxdbg) :: idbg !(0:maxdbg)
  integer(FINT), save, dimension(mni) :: ipfix !mni
  integer(FINT), save, dimension(maxstk) :: istkrd !maxstk
  integer(FINT), save, dimension(maxstk) :: istkwr !maxstk
  integer(FINT), save, dimension(7) :: isw !7
  integer(FINT), save, dimension(1:mni) :: nexofi !mni
  integer(FINT), save, dimension(mne) :: niofex !mne
  integer(FINT), save, dimension(mne) :: nvarl !mne
  integer(FINT), save, dimension(2) :: nwrmes !2
  !
  real(FREAL), save :: amin, apsi, bigedm, dcovar, edm, epsi, epsma2
  real(FREAL), save :: epsmac, fval3, undefi, up, updflt, vlimhi
  real(FREAL), save :: vlimlo, xdircr, xmidcr, ydircr, ymidcr
  !
  real(FREAL), save, dimension(mne) :: alim !(mne)
  real(FREAL), save, dimension(mne) :: blim !mne
  real(FREAL), save, dimension(mni) :: dgrd !mni
  real(FREAL), save, dimension(mni) :: dirin !mni
  real(FREAL), save, dimension(mni) :: dirins !mni
  real(FREAL), save, dimension(mni) :: ern !mni
  real(FREAL), save, dimension(mni) :: erp !mni
  real(FREAL), save, dimension(mni) :: g2 !mni
  real(FREAL), save, dimension(mni) :: g2s !mni
  real(FREAL), save, dimension(mne) :: gin !mne
  real(FREAL), save, dimension(mni) :: globcc !mni
  real(FREAL), save, dimension(mni) :: grd !mni
  real(FREAL), save, dimension(mni) :: grds !mni
  real(FREAL), save, dimension(mni) :: gstep !mni
  real(FREAL), save, dimension(mni) :: gsteps !mni
  real(FREAL), save, dimension(mni,mni+1) :: p !mni,mni+1
  real(FREAL), save, dimension(mni) :: pbar !mni
  real(FREAL), save, dimension(mni) :: prho !mni
  real(FREAL), save, dimension(mni) :: pstar !mni
  real(FREAL), save, dimension(mni) :: pstst !mni
  real(FREAL), save, dimension(mne) :: u !mne
  real(FREAL), save, dimension(mnihl) :: vhmat !mnihl
  real(FREAL), save, dimension(mnihl) :: vthmat !mnihl
  real(FREAL), save, dimension(mni) :: werr !mni
  real(FREAL), save, dimension(maxp) :: word7 !maxp
  real(FREAL), save, dimension(mni) :: x !mni
  real(FREAL), save, dimension(maxcpt) :: xpt !maxcpt
  real(FREAL), save, dimension(mni) :: xs !mni
  real(FREAL), save, dimension(mni) :: xt !mni
  real(FREAL), save, dimension(mni) :: xts !mni
  real(FREAL), save, dimension(maxcpt) :: ypt !maxcpt
  !
  character(FLCHARS), save :: cfrom
  character(FLCHARS), dimension(0:3) :: covmes
  character(FLCHARS), save :: cstatu
  character(FLCHARS), save :: ctitl
  character(FLCHARS), save :: cundef
  character(FLCHARS), save :: cvrsn
  character(FLCHARS), save :: cword
  !
  character(FLCHARS), save, dimension(maxcpt) :: chpt !maxcpt
  character(FLCHARS), save, dimension(mne) :: cpnam !mne
  !
contains
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function minuit_init(ird,iwr,isv)
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function minuit_init(ird,iwr,isv)
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the minuit module.
!H
!H-----------------------------------------------------------------------------
!H
!
    integer(FINT), intent(in) :: ird,iwr,isv
    !
    logical :: ltmp
    integer(FINT) :: i
    integer(FINT) :: idb
    real(FREAL) :: epstry
    real(FREAL) :: epsp1
    real(FREAL) :: epsbak
    real(FREAL) :: piby2
    real(FREAL) :: distnn
    !
    !integer(FINT), dimension(1:1) :: istkwr
    !integer(FINT), dimension(1:6) :: isw
    !integer(FINT), dimension(0:maxdbg) :: idbg
    !character(FLCHARS) :: cvrsn, cundef, cfrom, cstatu, ctitl
    !
    mnit_isinit = .true.
    minuit_init = 1
    !
    minuit_filein = ird
    minuit_fileout = iwr
    istkwr(1) = minuit_fileout
    nstkwr = 1
    minuit_filesave = isv
    nstkrd = 0
!               version identifier
    cvrsn = '96.03 '
!               some constant constants in common
    maxint=mni
    maxext=mne
    undefi = -54321.
    bigedm = 123456.
    cundef = ')UNDEFINED'
    covmes(0) = 'NO ERROR MATRIX       '
    covmes(1) = 'ERR MATRIX APPROXIMATE'
    covmes(2) = 'ERR MATRIX NOT POS-DEF'
    covmes(3) = 'ERROR MATRIX ACCURATE '
!                some starting values in COMMON
    nblock = 0
    icomnd = 0
    ctitl = cundef
    cfrom = 'INPUT   '
    nfcnfr = nfcn
    cstatu= 'INITIALIZE'
    isw(3) = 0
    isw(4) = 0
    isw(5) = 1
!         isw(6)=0 for batch jobs,  =1 for interactive jobs
!                      =-1 for originally interactive temporarily batch
! REMOVED NOW
    isw(6) = 0
!!     #ifndef cernlib_msstdcall
!      if (intrac(dummy))  isw(6) = 1
!!     #else
!      if (intrac())  isw(6) = 1
!!     #endif
!     debug options set to default values
    do idb= 0, maxdbg
       idbg(idb) = 0
    end do
    lrepor = .false.
    lwarn  = .true.
    limset = .false.
    lnewmn = .false.
    istrat = 1
    itaur = 0
!        default page dimensions and 'new page' carriage control integer
    npagwd = 120
    npagln = 56
    newpag = 1
    if (isw(6).gt.0) then
       npagwd = 80
       npagln = 30
       newpag = 0
    end if
    up = 1.0
    updflt = up
!                   determine machine accuracy epsmac
    epstry = 0.5
    ltmp = .true.
    do i= 1, 100
       epstry = epstry * 0.5
       epsp1 = one + epstry
       call mntiny(epsp1, epsbak)
       if (epsbak.lt.epstry) then
          ltmp = .false.
          exit
       end if
    end do
    if(ltmp) then
       epstry = 1.0e-7
       epsmac = 4.0*epstry
       write (minuit_fileout,'(A,A,E10.2)') ' MNINIT UNABLE TO DETERMINE',&
            &' ARITHMETIC PRECISION. WILL ASSUME:',epsmac
    end if
35  epsmac = 8.0 * epstry
    epsma2 = 2.0 * sqrt(epsmac)
!                 the vlims are a non-negligible distance from pi/2
!         used by mnpint to set variables "near" the physical limits
    piby2 = 2.0*atan(1.0)
    distnn = 8.0*sqrt(epsma2)
    vlimhi =  piby2 - distnn
    vlimlo = -piby2 + distnn
    call mncler
    write (minuit_fileout,'(3A,I3,A,I3,A,E10.2)')  '  MINUIT RELEASE ',cvrsn,&
         &' INITIALIZED.   DIMENSIONS ',mne,'/',mni,'  EPSMAC=',epsmac
    mnit_isinit = .true.
    return
    !
  end function minuit_init
!
!H
!H-----------------------------------------------------------------------------
!H logical function minuit_isinit()
!H-----------------------------------------------------------------------------
!H
!
  logical function minuit_isinit()
!
!H
!H-----------------------------------------------------------------------------
!H We merelyreturn the initialization flag.
!H-----------------------------------------------------------------------------
!H
!
    !
    minuit_isinit = mnit_isinit
    !
  end function minuit_isinit
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function minuit_geterror_i()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function minuit_geterror_i()
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
    minuit_geterror_i = mnit_error
    !
  end function minuit_geterror_i
!
!H-----------------------------------------------------------------------------
!H
!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H character(FLCHARS) function minuit_geterror_c()
!H-----------------------------------------------------------------------------
!H
!
  character(FLCHARS) function minuit_geterror_c(code)
    !
    !integer(FINT), optional, intent(in) :: code
    !
    ! removed optional due to IFC confusion
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
    integer(FINT) :: lcode !the local used error code
    !
    lcode = mnit_error
    ! It is not optional anymore!!
    !if(present(code)) then
       lcode = code
    !endif
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
       minuit_geterror_c = "No Errors."
    case(1_FINT)
       minuit_geterror_c = "Not documented error or generic error.."
    case(2_FINT)
       minuit_geterror_c = "Double Error: Why did you get this message?."
    case default
       minuit_geterror_c = "Not documented error or generic error.."
    end select
    !
    return
    !
  end function minuit_geterror_c
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H HERE STARTS THE PORTED CODE!
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine minuit_mnparm(k,cnamj,uk,wk,a,b,ierflg)
        !
!c        called from mnpars and user-callable
!c    implements one parameter definition, that is:
!c          k     (external) parameter number
!c          cnamk parameter name
!c          uk    starting value
!c          wk    starting step size or uncertainty
!c          a, b  lower and upper physical parameter limits
!c    and sets up (updates) the parameter lists.
!c    output: ierflg=0 if no problems
!c                  >0 if mnparm unable to implement definition
!c
    character*(*) cnamj
    character  cnamk*10, chbufi*4
    !
    integer(FINT) :: kint
    integer(FINT) :: k, ierflg, ktofix, ix, nvl, lastin, in, kinfix, ierr
    real(FREAL) :: uk, wk, a, b, sav, danger, pinti, sav2, vplu
    real(FREAL) :: vminu, gsmin, small
    !
    cnamk = cnamj
    kint = npar
    !
    if (k.lt.1.or.k.gt.maxext) then
!                     parameter number exceeds allowed maximum value
       write (minuit_fileout,9)  k,maxext
9      format (/' MINUIT USER ERROR.  PARAMETER NUMBER IS',I11/&
            &',  ALLOWED RANGE IS ONE TO',I4/)
       go to 800
    endif
!                     normal parameter request
    ktofix = 0
    if (nvarl(k).lt.0) go to 50
!         previously defined parameter is being redefined
!                                     find if parameter was fixed
    do ix= 1, npfix
       if (ipfix(ix).eq.k)  ktofix = k
    end do
    if (ktofix.gt.0)  then
       call mnwarn('W','PARAM DEF','REDEFINING A FIXED PARAMETER.')
       if (kint.ge.maxint)  then
          write (minuit_fileout,'(A)') ' CANNOT RELEASE. MAX NPAR EXCEEDED.'
          go to 800
       endif
       call mnfree(-k)
    endif
!                       if redefining previously variable parameter
    if(niofex(k).gt.0) kint = npar-1
50  continue
!
!                                      . . .print heading
    if (lphead.and.isw(5).ge.0) then
       write (minuit_fileout,61)
       lphead = .false.
    endif
61  format(/' PARAMETER DEFINITIONS:'/&
         &        '    NO.   NAME         VALUE      STEP SIZE      LIMITS')
    if (wk.gt.zero)  go to 122
!                                        . . .constant parameter . . . .
    if (isw(5).ge.0)  write (minuit_fileout, 82)  k,cnamk,uk
82  format (1X,I5,1X,1H',A10,1H',1X,G13.5, '  CONSTANT')
    nvl = 0
    go to 200
122 if (a.eq.zero.and.b.eq.zero) then
!                                      variable parameter without limits
       nvl = 1
       if (isw(5).ge.0)  write (minuit_fileout, 127)  k,cnamk,uk,wk
127    format (1X,I5,1X,1H',A10,1H',1X,2G13.5, '     NO LIMITS')
    else
!                                         variable parameter with limits
       nvl = 4
       lnolim = .false.
       if (isw(5).ge.0)  write (minuit_fileout, 132)  k,cnamk,uk,wk,a,b
132    format(1X,I5,1X,1H',A10,1H',1X,2G13.5,2X,2G13.5)
    endif
!                             . . request for another variable parameter
    kint = kint + 1
    if (kint.gt.maxint)  then
       write (minuit_fileout,135)  maxint
135    format (/' MINUIT USER ERROR.   TOO MANY VARIABLE PARAMETERS.'/&
            &   ' THIS VERSION OF MINUIT DIMENSIONED FOR',I4//)
       go to 800
    endif
    if (nvl.eq.1)  go to 200
    if (a.eq.b)  then
       write (minuit_fileout,'(/A,A/A/)') ' USER ERROR IN MINUIT PARAMETER',&
            &   ' DEFINITION',' UPPER AND LOWER LIMITS EQUAL.'
       go to 800
    endif
    if (b.lt.a) then
       sav = b
       b = a
       a = sav
       call mnwarn('W','PARAM DEF','PARAMETER LIMITS WERE REVERSED.')
       if (lwarn) lphead=.true.
    endif
    if ((b-a).gt.1.0e7)  then
       write (chbufi,'(I4)') k
       call mnwarn('W','PARAM DEF',&
            &               'LIMITS ON PARAM'//chbufi//' TOO FAR APART.')
       if (lwarn) lphead=.true.
    endif
    danger = (b-uk)*(uk-a)
    if (danger.lt.0.)&
         &     call mnwarn('W','PARAM DEF','STARTING VALUE OUTSIDE LIMITS.')
    if (danger.eq.0.)&
         &     call mnwarn('W','PARAM DEF','STARTING VALUE IS AT LIMIT.')
200 continue
!                           . . . input ok, set values, arrange lists,
!                                    calculate step sizes gstep, dirin
    cfrom = 'PARAMETR'
    nfcnfr = nfcn
    cstatu= 'NEW VALUES'
    nu = max(nu,k)
    cpnam(k) = cnamk
    u(k) = uk
    alim(k) = a
    blim(k) = b
    nvarl(k) = nvl
!                             k is external number of new parameter
!           lastin is the number of var. params with ext. param. no.< k
    lastin = 0
    do ix= 1, k-1
       if (niofex(ix).gt.0)  lastin=lastin+1
    end do
!                 kint is new number of variable params, npar is old
    if (kint.eq.npar)  go to 280
    if (kint.gt.npar) then
!                          insert new variable parameter in list
       do in= npar,lastin+1,-1
          ix = nexofi(in)
          niofex(ix) = in+1
          nexofi(in+1)= ix
          x    (in+1) = x    (in)
          xt   (in+1) = xt   (in)
          dirin(in+1) = dirin(in)
          g2   (in+1) = g2   (in)
          gstep(in+1) = gstep(in)
       end do
    else
!                          remove variable parameter from list
       do in= lastin+1,kint
          ix = nexofi(in+1)
          niofex(ix) = in
          nexofi(in)= ix
          x     (in)= x    (in+1)
          xt    (in)= xt   (in+1)
          dirin (in)= dirin(in+1)
          g2    (in)= g2   (in+1)
          gstep (in)= gstep(in+1)
       end do
    endif
280 continue
    ix = k
    niofex(ix) = 0
    npar = kint
    call mnrset(1)
!                                       lists are now arranged . . . .
    if (nvl.gt.0)  then
       in = lastin+1
       nexofi(in) = ix
       niofex(ix) = in
       sav = u(ix)
       call mnpint(sav,ix,pinti)
       x(in) = pinti
       xt(in) = x(in)
       werr(in) = wk
       sav2 = sav + wk
       call mnpint(sav2,ix,pinti)
       vplu = pinti - x(in)
       sav2 = sav - wk
       call mnpint(sav2,ix,pinti)
       vminu = pinti - x(in)
       dirin(in) = 0.5 * (abs(vplu) +abs(vminu))
       g2(in) = 2.0*up / dirin(in)**2
       gsmin = 8.*epsma2*abs(x(in))
       gstep(in) = max (gsmin, 0.1_FREAL*dirin(in))
       if (amin.ne.undefi) then
          small = sqrt(epsma2*(amin+up)/up)
          gstep(in) = max(gsmin, small*dirin(in))
       endif
       grd  (in) = g2(in)*dirin(in)
!                   if parameter has limits
       if (nvarl(k).gt.1) then
          if (gstep(in).gt.0.5)  gstep(in)=0.5
          gstep(in) = -gstep(in)
       endif
    endif
    if (ktofix.gt.0)  then
       kinfix = niofex(ktofix)
       if (kinfix.gt.0)  call mnfixp(kinfix,ierr)
       if (ierr.gt.0)  go to 800
    endif
    ierflg = 0
    return
!                   error on input, unable to implement request  . . . .
800 continue
    ierflg = 1
    return
  end subroutine minuit_mnparm
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
      subroutine minuit_mnseti(tit)

        character(*) :: tit
        ctitl = tit
        return
      end subroutine minuit_mnseti
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
      subroutine mntiny(epsp1,epsbak)
!
!cc        compares its argument with the value 1.0, and returns
!cc        the value .true. if they are equal.  to find epsmac
!cc        safely by foiling the fortran optimizer
!cc
!
        real(FREAL) :: epsp1,epsbak
        real(FREAL), parameter :: one=1.0
        !
        epsbak =  epsp1  - one
        return
      end subroutine mntiny
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
      subroutine mncler
!c        called from minuit and by option from mnexcm
!c        resets the parameter list to undefined
        !
        integer(FINT) :: i
        !
        npfix = 0
        nu = 0
        npar = 0
        nfcn = 0
        nwrmes(1) = 0
        nwrmes(2) = 0
        do i= 1, maxext
           u(i) = 0.0
           cpnam(i) = cundef
           nvarl(i) = -1
           niofex(i) = 0
        end do
        call mnrset(1)
        cfrom = 'CLEAR   '
        nfcnfr = nfcn
        cstatu ='UNDEFINED '
        lnolim = .true.
        lphead = .true.
        return
      end subroutine mncler
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
      subroutine mnwarn(copt,corg,cmes)
!     if copt='w', cmes is a warning message from corg.
!     if copt='d', cmes is a debug message from corg.
!         if set warnings is in effect (the default), this routine
!             prints the warning message cmes coming from corg.
!         if set nowarnings is in effect, the warning message is
!             stored in a circular buffer of length maxmes.
!         if called with corg=cmes='sho', it prints the messages in
!             the circular buffer, fifo, and empties the buffer.

      character copt*1, corg*(*), cmes*(*), ctyp*7
      integer(FINT), parameter :: maxmes=10
      integer(FINT) :: i
      character     origin(maxmes,2)*10, warmes(maxmes,2)*60
      common/mn7wrc/origin,              warmes
      real(FREAL), save, dimension(maxmes,2) :: nfcwar
      real(FREAL), save, dimension(2) :: icirc
      !common/mn7wri/nfcwar(maxmes,2),icirc(2)
      character englsh*20
      !
      integer(FINT) :: ityp, ic, nm
      !
!
      if (corg(1:3).eq.'SHO'.and.cmes(1:3).eq.'SHO')  go to 200
!             either print warning or put in buffer
      if (copt.eq.'W')  then
        ityp = 1
        if (lwarn) then
          write (minuit_fileout,'(A,A/A,A)') ' MINUIT WARNING IN ',corg,&
     &              ' ============== ',cmes
          return
        endif
      else
        ityp = 2
        if (lrepor) then
          write (minuit_fileout,'(A,A/A,A)') ' MINUIT DEBUG FOR  ',corg,&
     &              ' ============== ',cmes
          return
        endif
      endif
!                 if appropriate flag is off, fill circular buffer
      if (nwrmes(ityp).eq.0)  icirc(ityp) = 0
      nwrmes(ityp) = nwrmes(ityp) + 1
      icirc(ityp) = icirc(ityp) + 1
      if (icirc(ityp).gt.maxmes) icirc(ityp) = 1
      ic = icirc(ityp)
      origin(ic,ityp) = corg
      warmes(ic,ityp) = cmes
      nfcwar(ic,ityp) = nfcn
      return
!
!             'sho warnings', ask if any suppressed mess in buffer
  200 continue
      if (copt.eq.'W') then
        ityp = 1
        ctyp = 'WARNING'
      else
        ityp = 2
        ctyp = '*DEBUG*'
      endif
      if (nwrmes(ityp).gt.0) then
         englsh = ' WAS SUPPRESSED.  '
         if (nwrmes(ityp).gt.1) englsh = 'S WERE SUPPRESSED.'
         write (minuit_fileout,'(/1X,I5,A,A,A,A/)') nwrmes(ityp),&
     &    ' MINUIT ',ctyp,' MESSAGE', englsh
         nm = nwrmes(ityp)
         ic = 0
         if (nm.gt.maxmes) then
              write (minuit_fileout,'(A,I2,A)')  ' ONLY THE MOST RECENT ',&
     &          maxmes,' WILL BE LISTED BELOW.'
              nm = maxmes
              ic = icirc(ityp)
         endif
         write (minuit_fileout,'(A)') '  CALLS  ORIGIN         MESSAGE'
           do i= 1, nm
              ic = ic + 1
              if (ic.gt.maxmes)  ic = 1
              write (minuit_fileout,'(1X,I6,1X,A,1X,A)')&
                   & nfcwar(ic,ityp),origin(ic,ityp),warmes(ic,ityp)
           end do
           nwrmes(ityp) = 0
           write (minuit_fileout,'(1H )')
        endif
        return
      end subroutine mnwarn

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
      subroutine mnrset(iopt)
!c        called from mncler and whenever problem changes, for example
!c        after set limits, set param, call minuit_fcn 6
!c    if iopt=1,
!c        resets function value and errors to undefined
!c    if iopt=0, sets only minos errors to undefined
        !
        integer(FINT) :: iopt, iext
        integer(FINT) :: i
        !
        cstatu = 'RESET     '
        if (iopt.ge.1)  then
           amin = undefi
           fval3 = 2.0*abs(amin) + 1.
           edm = bigedm
           isw(4) = 0
           isw(2) = 0
           dcovar = 1.
           isw(1) = 0
        endif
        lnolim = .true.
        do i= 1, npar
           iext = nexofi(i)
           !F. Mariotti
           if (iext.gt.0) then
              if (nvarl(iext).ge.4) lnolim=.false.
           end if
           erp(i) = zero
           ern(i) = zero
           globcc(i) = zero
        end do
        if (isw(2).ge.1)  then
           isw(2) = 1
           dcovar = max(dcovar,half)
        endif
        return
      end subroutine mnrset

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnfree(k)
!c        restores one or more fixed parameter(s) to variable status
!c        by inserting it into the internal parameter list at the
!c        appropriate place.
!c
!--       k = 0 means restore all parameters
!--       k = 1 means restore the last parameter fixed
!--       k = -i means restore external parameter i (if possible)
!--       iq = fix-location where internal parameters were stored
!--       ir = external number of parameter being restored
!--       is = internal number of parameter being restored


    integer(FINT) :: k, ka, ik, ipsav, ir, is, lc, iq, i

    real(FREAL) :: xv, xtv, dirinv, grdv, g2v, gstepv


    if (k.gt.1)  write (minuit_fileout,510)
    if (npfix.lt.1)  write (minuit_fileout,500)
    if (k.eq.1.or.k.eq.0)  go to 40
!                   release parameter with specified external number
    ka = iabs(k)
    if (niofex(ka).eq.0)  go to 15
    write (minuit_fileout,540)
540 format (' IGNORED.  PARAMETER SPECIFIED IS ALREADY VARIABLE.')
    return
15  if (npfix.lt.1)  go to 21
    do ik= 1, npfix
       if (ipfix(ik).eq.ka)  go to 24
    end do
21  write (minuit_fileout,530) ka
530 format (' PARAMETER',I4,' NOT FIXED.  CANNOT BE RELEASED.')
    return
24  if (ik.eq.npfix)  go to 40
!                   move specified parameter to end of list
    ipsav = ka
    xv = xs(ik)
    xtv = xts(ik)
    dirinv = dirins(ik)
    grdv = grds(ik)
    g2v = g2s(ik)
    gstepv = gsteps(ik)
    do i= ik+1,npfix
       ipfix(i-1) = ipfix(i)
       xs(i-1) = xs(i)
       xts(i-1) = xts(i)
       dirins(i-1) = dirins(i)
       grds(i-1) = grds(i)
       g2s(i-1) = g2s(i)
       gsteps(i-1) = gsteps(i)
    end do
    ipfix(npfix) = ipsav
    xs(npfix) = xv
    xts(npfix) = xtv
    dirins(npfix) = dirinv
    grds(npfix) = grdv
    g2s(npfix) = g2v
    gsteps(npfix) = gstepv
    !                restore last parameter in fixed list  -- ipfix(npfix)
40  continue
    if (npfix.lt.1)  go to 300
    ir = ipfix(npfix)
    is = 0
    do ik= nu, ir, -1
       if (niofex(ik).gt.0) then
          lc = niofex(ik) + 1
          is = lc - 1
          niofex(ik) = lc
          nexofi(lc) = ik
          x(lc)     = x(lc-1)
          xt(lc)    = xt(lc-1)
          dirin(lc) = dirin(lc-1)
          werr(lc)  = werr(lc-1)
          grd(lc)   = grd(lc-1)
          g2(lc)    = g2(lc-1)
          gstep(lc) = gstep(lc-1)
       endif
    end do
    npar = npar + 1
    if (is.eq.0)   is = npar
    niofex(ir) = is
    nexofi(is) = ir
    iq = npfix
    x(is) = xs(iq)
    xt(is) = xts(iq)
    dirin(is) = dirins(iq)
    werr(is)  = dirins(iq)
    grd(is) = grds(iq)
    g2(is) = g2s(iq)
    gstep(is) = gsteps(iq)
    npfix = npfix - 1
    isw(2) = 0
    dcovar = 1.
    if (isw(5)-itaur.ge.1)  write(minuit_fileout,520) ir,cpnam(ir)
    if (k.eq.0)  go to 40
300 continue
    !         if different from internal, external values are taken
    call mnexin(x)
400 return
500 format (' CALL TO MNFREE IGNORED.  THERE ARE NO FIXED PA',&
         & 'RAMETERS'/)
510 format (' CALL TO MNFREE IGNORED.  ARGUMENT GREATER THAN ONE'/)
520 format (20X, 9HPARAMETER,I4,2H, ,A10,' RESTORED TO VARIABLE.')
  end subroutine mnfree

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
      subroutine mnexin(pint)
!
!cc        transforms the external parameter values u to internal
!cc        values in the dense array pint. subroutine mnpint is used.
!cc
      real(FREAL), dimension(:) :: pint !(*)
      !
      integer(FINT) :: liint, iext
      REAL(FREAL) :: pinti
      !
      limset = .false.
      do liint= 1, npar
         iext = nexofi(liint)
         call mnpint(u(iext),iext,pinti)
         pint(liint) = pinti
      end do
      return
    end subroutine mnexin

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
      subroutine mnpint(pexti,i,pinti)
!c        calculates the internal parameter value pinti corresponding
!c        to the external value pexti for parameter i.
!c
        integer(FINT) :: i
        real(FREAL) :: pexti,pinti
      character chbufi*4, chbuf2*30

      integer(FINT) :: igo
      real(FREAL) :: alimi, blimi, yy, yy2, a

      pinti = pexti
      igo = nvarl(i)
      if (igo.eq.4)  then
!--                          there are two limits
        alimi = alim(i)
        blimi = blim(i)
        yy=2.0*(pexti-alimi)/(blimi-alimi) - 1.0
        yy2 = yy**2
        if (yy2.ge.(1.0- epsma2))  then
           if (yy.lt.0.) then
               a = vlimlo
               chbuf2 = 'IS AT ITS LOWER ALLOWED LIMIT.'
           else
               a = vlimhi
               chbuf2 = 'IS AT ITS UPPER ALLOWED LIMIT.'
           endif
           pinti = a
           pexti = alimi + 0.5* (blimi-alimi) *(sin(a) +1.0)
           limset = .true.
           write (chbufi,'(I4)') i
           if (yy2.gt.1.0) chbuf2 = ' BROUGHT BACK INSIDE LIMITS.'
           call mnwarn('W',cfrom,'VARIABLE'//chbufi//chbuf2)
         else
           pinti = asin(yy)
         endif
      endif
      return
    end subroutine mnpint

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
      subroutine mnfixp(liint,ierr)
!c        removes parameter liint from the internal (variable) parameter
!c        list, and arranges the rest of the list to fill the hole.
!c
        integer(FINT) :: liint, ierr
        integer(FINT) :: i
        real(FREAL), dimension(mni) :: yy !(mni)
        !
        integer(FINT) :: iext, nold, lc, ik, m, n, ndex, knew, kold, j
        real(FREAL) :: yyover
        !
!                           first see if it can be done
      ierr = 0
      if (liint.gt.npar.or.liint.le.0)  then
         ierr = 1
         write (minuit_fileout,'(A,I4)')&
     &       ' MINUIT ERROR.  ARGUMENT TO MNFIXP=',liint
         go to 300
      endif
      iext = nexofi(liint)
      if (npfix.ge.mni) then
         ierr = 1
         write (minuit_fileout,'(A,I4,A,I4)') ' MINUIT CANNOT FIX PARAMETER',&
     &   IEXT,' MAXIMUM NUMBER THAT CAN BE FIXED IS',mni
         go to 300
      endif
!                           reduce number of variable parameters by one
      niofex(iext) = 0
      nold = npar
      npar = npar - 1
!                       save values in case parameter is later restored
      npfix = npfix + 1
      ipfix(npfix) = iext
      lc = liint
      xs(npfix) = x(lc)
      xts(npfix) = xt(lc)
      dirins(npfix) = werr(lc)
      grds(npfix) = grd(lc)
      g2s(npfix) = g2(lc)
      gsteps(npfix) = gstep(lc)
!                        shift values for other parameters to fill hole
      do ik= iext+1, nu
         if  (niofex(ik).gt.0)  then
            lc = niofex(ik) - 1
            niofex(ik) = lc
            nexofi(lc) = ik
            x(lc)     = x(lc+1)
            xt(lc)    = xt(lc+1)
            dirin(lc) = dirin(lc+1)
            werr(lc)  = werr(lc+1)
            grd(lc)   = grd(lc+1)
            g2(lc)    = g2(lc+1)
            gstep(lc) = gstep(lc+1)
         endif
      end do
      if (isw(2).le.0)  go to 300
!                    remove one row and one column from variance matrix
      if (npar.le.0)  go to 300
      do i= 1, nold
         m = max(i,liint)
         n = min(i,liint)
         ndex = m*(m-1)/2 + n
         yy(i)=vhmat(ndex)
      end do
      yyover = 1.0/yy(liint)
      knew = 0
      kold = 0
      do i= 1, nold
         do j= 1, i
            kold = kold + 1
            if (.not.(j.eq.liint.or.i.eq.liint)) then
               knew = knew + 1
               vhmat(knew) = vhmat(kold) - yy(j)*yy(i)*yyover
            end if
         end do
      end do
300   return
    end subroutine mnfixp

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
   subroutine minuit_mnexcm(comand,plist,llist,ierflg)
!c        interprets a command and takes appropriate action,
!c        either directly by skipping to the corresponding code in
!c        mnexcm, or by setting up a call to a subroutine
!c
     !
     character(*) :: comand
!   cannot say dimension plist(llist) since llist can be =0.
     real(FREAL), dimension(:) :: plist !(*)
     integer(FINT) :: llist, ierflg
     integer(FINT) :: i
     !
     integer(FINT), parameter :: mxpt=101
     real(FREAL), dimension(mxpt) ::  xptu, yptu !(mxpt)
     integer(FINT) :: nntot
!  alphabetical order of command names!
     !character*10 cname(40), cneway, chwhy*18, c26*30, cvblnk*2
     character(10), dimension(40) :: cname !(40)
     character(10) :: cneway
     character(18) :: chwhy !*18
     character(30) :: c26 !*30
     character(2) :: cvblnk !*2
     logical ltofix, lfixed, lfreed
     !
     !
     integer(FINT) :: lk, icol, let, iw, lnow, inonde, kll, nf
     integer(FINT) :: nsuper, ilist, iext, liint, ierr, krl, it
     integer(FINT) :: it2, ke1, ke2, ierrf, iflag, nparx, nowprt
     integer(FINT) :: kcol, nptu, izero
     real(FREAL) :: f, step, rno
     !
     character(4) :: comd !*4
     character(26) :: clower, cupper !*26
     data clower/'abcdefghijklmnopqrstuvwxyz'/
     data cupper/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
!
!  recognized minuit commands:
     data cname( 1) / 'MINIMIZE  ' /
     data cname( 2) / 'SEEK      ' /
     data cname( 3) / 'SIMPLEX   ' /
     data cname( 4) / 'MIGRAD    ' /
     data cname( 5) / 'MINOS     ' /
     data cname( 6) / 'SET XXX   ' /
     data cname( 7) / 'SHOW XXX  ' /
     data cname( 8) / 'TOP OF PAG' /
     data cname( 9) / 'FIX       ' /
     data cname(10) / 'RESTORE   ' /
     data cname(11) / 'RELEASE   ' /
     data cname(12) / 'SCAN      ' /
     data cname(13) / 'CONTOUR   ' /
     data cname(14) / 'HESSE     ' /
     data cname(15) / 'SAVE      ' /
     data cname(16) / 'IMPROVE   ' /
     data cname(17) / 'CALL FCN  ' /
     data cname(18) / 'STANDARD  ' /
     data cname(19) / 'END       ' /
     data cname(20) / 'EXIT      ' /
     data cname(21) / 'RETURN    ' /
     data cname(22) / 'CLEAR     ' /
     data cname(23) / 'HELP      ' /
     data cname(24) / 'MNCONTOUR ' /
     data cname(25) / 'STOP      ' /
     data cname(26) / 'JUMP      ' /
     data cname(27) / '          ' /
     data cname(28) / '          ' /
     data cname(29) / '          ' /
     data cname(30) / '          ' /
     data cname(31) / '          ' /
     data cname(32) / '          ' /
     data cname(33) / '          ' /
!  OBSOLETE COMMANDS:
     data cname(34) / 'COVARIANCE' /
     data cname(35) / 'PRINTOUT  ' /
     data cname(36) / 'GRADIENT  ' /
     data cname(37) / 'MATOUT    ' /
     data cname(38) / 'ERROR DEF ' /
     data cname(39) / 'LIMITS    ' /
     data cname(40) / 'PUNCH     ' /
     data nntot/40/
     !
!      ierflg is now (94.5) defined the same as icondn in mncomd
!c            = 0: command executed normally
!c              1: command is blank, ignored
!c              2: command line unreadable, ignored
!c              3: unknown command, ignored
!c              4: abnormal termination (e.g., migrad not converged)
!c              9: reserved
!c             10: end command
!c             11: exit or stop command
!c             12: return command
     lk = len(comand)
     if (lk.gt.maxcwd) lk=maxcwd
     cword = comand(1:lk)
!              get upper case
     do icol= 1, lk
        do let= 1, 26
           if (cword(icol:icol).eq.clower(let:let))&
                &cword(icol:icol) = cupper(let:let)
        end do
     end do
!           copy the first maxp arguments into common (word7), making
!           sure that word7(1)=0. if llist=0
     do iw= 1, maxp
        word7(iw) = zero
        if (iw.le.llist) word7(iw) = plist(iw)
     end do
     icomnd = icomnd + 1
     nfcnlc = nfcn
     if (cword(1:7).ne.'SET PRI'.or.word7(1).ge.0.)  then
        if (isw(5).ge.0) then
           lnow = llist
           if (lnow.gt.4)  lnow=4
           write (minuit_fileout,"(1H ,10(1H*)/' **',I5,' **',A,4G12.4)")&
                &icomnd, cword(1:lk), (plist(i),i=1,lnow)
           inonde = 0
           if (llist.gt.lnow) then
              kll = llist
              if (llist.gt.maxp) then
                 inonde = 1
                 kll = maxp
              endif
              write (cvblnk,'(I2)') lk
              c26 = '(11H **********,'//cvblnk//'X,4G12.4)'
              write (minuit_fileout,c26) (plist(i),i=lnow+1,kll)
           endif
           write (minuit_fileout, '(1H ,10(1H*))' )
           if (inonde.gt.0)&
                &write (minuit_fileout, '(1H ,10(1H*),A,I3,A)')&
                &'  ERROR: ABOVE CALL TO MNEXCM TRIED TO PASS MORE THAN ',&
                &maxp,' PARAMETERS.'
        endif
     endif
     nfcnmx = word7(1)
     if (nfcnmx.le.0)  nfcnmx = 200 + 100*npar + 5*npar**2
     epsi = word7(2)
     if (epsi.le.zero)  epsi = 0.1 * up
     lnewmn = .false.
     lphead = .true.
     isw(1) = 0
     ierflg = 0
!                look for command in list cname . . . . . . . . . .
     do i= 1, nntot
        if (cword(1:3).eq.cname(i)(1:3))  go to 90
     end do
     write (minuit_fileout,'(11X,''UNKNOWN COMMAND IGNORED:'',A)') comand
     ierflg = 3
     go to 5000
!                normal case: recognized minuit command . . . . . . .
90   continue
     if (cword(1:4).eq.'MINO') i = 5
     if (i.ne.6.and.i.ne.7.and.i.ne.8.and.i.ne.23)  then
        cfrom = cname(i)
        nfcnfr = nfcn
     endif
!              1    2    3    4    5    6    7    8    9   10
      go to ( 400, 200, 300, 400, 500, 700, 700, 800, 900,1000,&
     &       1100,1200,1300,1400,1500,1600,1700,1800,1900,1900,&
     &       1900,2200,2300,2400,1900,2600,3300,3300,3300,3300,&
     &       3300,3300,3300,3400,3500,3600,3700,3800,3900,4000) , i
!                                        . . . . . . . . . . seek
200   call mnseek
      go to 5000
!                                        . . . . . . . . . . simplex
300   call mnsimp
      if (isw(4).lt.1)  ierflg = 4
      go to 5000
!                                        . . . . . . migrad, minimize
400   continue
      nf = nfcn
      apsi = epsi
      call mnmigr
      call mnwerr
      if (isw(4).ge.1)         go to 5000
      ierflg = 4
      if (isw(1).eq.1)         go to 5000
      if (cword(1:3).eq.'MIG') go to 5000
      nfcnmx = nfcnmx + nf - nfcn
      nf = nfcn
      call mnsimp
      if (isw(1).eq.1)  go to 5000
      nfcnmx = nfcnmx + nf - nfcn
      call mnmigr
      if (isw(4).ge.1)  ierflg = 0
      call mnwerr
      go to 5000
!                                        . . . . . . . . . . minos
500   continue
      nsuper = nfcn + 2*(npar+1)*nfcnmx
!          possible loop over new minima
      epsi = 0.1 * up
510   continue
      call mncuve
      call mnmnos
      if (.not.lnewmn)  go to 5000
      call mnrset(0)
      call mnmigr
      call mnwerr
      if (nfcn.lt.nsuper)  go to 510
      write (minuit_fileout,'(/'' TOO MANY FUNCTION CALLS. MINOS GIVES UP''/)')
      ierflg = 4
      go to 5000
!                                        . . . . . . . . . .set, show
700   call mnset
      go to 5000
!                                        . . . . . . . . . . top of page
800   continue
      write (minuit_fileout,'(1H1)')
      go to 5000
!                                        . . . . . . . . . . fix
900   ltofix = .true.
!                                        . . (also release) ....
901   continue
      lfreed = .false.
      lfixed = .false.
      if (llist.eq.0)  then
         write (minuit_fileout,'(A,A)') cword,':  NO PARAMETERS REQUESTED '
         go to 5000
      endif
      do ilist= 1, llist
         iext = plist(ilist)
         chwhy = ' IS UNDEFINED.'
         if (iext.le.0)         go to 930
         if (iext.gt.nu)        go to 930
         if (nvarl(iext).lt.0)  go to 930
         chwhy = ' IS CONSTANT.  '
         if (nvarl(iext).eq.0)  go to 930
         liint = niofex(iext)
         if (ltofix) then
            chwhy = ' ALREADY FIXED.'
            if (liint.eq.0)      go to 930
            call mnfixp(liint,ierr)
            if (ierr.eq.0) then
               lfixed = .true.
            else
               ierflg = 4
            endif
         else
            chwhy = ' ALREADY VARIABLE.'
            if (liint.gt.0)      go to 930
            krl = -iabs(iext)
            call mnfree(krl)
            lfreed = .true.
         endif
         cycle
930      write (minuit_fileout,'(A,I4,A,A)') ' PARAMETER',iext,chwhy,' IGNORED.'
      end do
      if (lfreed.or.lfixed)  call mnrset(0)
      if (lfreed)  then
         isw(2) = 0
         dcovar = 1.
         edm = bigedm
         isw(4) = 0
      endif
      call mnwerr
      if (isw(5).gt.1)  call mnprin(5,amin)
      go to 5000
!                                        . . . . . . . . . . restore
1000  it = word7(1)
      if (it.gt.1.or.it.lt.0)  go to 1005
      lfreed = (npfix.gt.0)
      call mnfree(it)
      if (lfreed) then
         call mnrset(0)
         isw(2) = 0
         dcovar = 1.
         edm = bigedm
      endif
      go to 5000
1005  write (minuit_fileout,'(A,I4)') ' IGNORED.  UNKNOWN ARGUMENT:',it
      ierflg = 3
      go to 5000
!                                        . . . . . . . . . . release
1100  ltofix = .false.
      go to 901
!                                       . . . . . . . . . . scan . . .
1200  continue
      iext = word7(1)
      if (iext.le.0)  go to 1210
      it2 = 0
      if (iext.le.nu)  it2 = niofex(iext)
      if (it2.le.0)  go to 1250
1210  call mnscan
      go to 5000
1250  write (minuit_fileout,'(A,I4,A)') ' PARAMETER',iext,' NOT VARIABLE.'
      ierflg = 3
      go to 5000
!                                        . . . . . . . . . . contour
1300  continue
      ke1 = word7(1)
      ke2 = word7(2)
      if (ke1.eq.0)  then
         if (npar.eq.2)  then
            ke1 = nexofi(1)
            ke2 = nexofi(2)
         else
            write (minuit_fileout,'(A,A)') cword,':  NO PARAMETERS REQUESTED '
            ierflg = 3
            go to 5000
         endif
      endif
      nfcnmx = 1000
      call mncntr(ke1,ke2,ierrf)
      if (ierrf.gt.0)  ierflg = 3
      go to 5000
!                                        . . . . . . . . . . hesse
1400  continue
      call mnhess
      call mnwerr
      if (isw(5).ge.0)  call mnprin(2, amin)
      if (isw(5).ge.1)  call mnmatu(1)
      go to 5000
!                                        . . . . . . . . . . save
1500  continue
      call mnsave
      go to 5000
!                                        . . . . . . . . . . improve
1600  continue
      call mncuve
      call mnimpr
      if (lnewmn)  go to 400
      ierflg = 4
      go to 5000
!                                        . . . . . . . . . . call fcn
1700  iflag = word7(1)
      nparx = npar
      f = undefi
      call minuit_fcn(nparx,gin,f,u,iflag)
      nfcn = nfcn + 1
      nowprt = 0
      if (f.ne.undefi)  then
         if (amin.eq.undefi)  then
            amin = f
            nowprt = 1
         else if (f.lt.amin)  then
            amin = f
            nowprt = 1
         endif
         if (isw(5).ge.0.and.iflag.le.5.and.nowprt.eq.1)&
              &          call mnprin(5,amin)
         if (iflag.eq.3)  fval3=f
      endif
      if (iflag.gt.5)  call mnrset(1)
      go to 5000
!                                        . . . . . . . . . . standard
1800  call stand
      go to 5000
!                                       . . . return, stop, end, exit
1900  it = word7(1)
      if (fval3.ne.amin.and.it.eq.0)  then
         iflag = 3
         if (isw(5).ge.0)&
              &write (minuit_fileout,'(/A/)')&
              &' CALL TO USER FUNCTION WITH IFLAG = 3'
         nparx = npar
         call minuit_fcn(nparx,gin,f,u,iflag)
         nfcn = nfcn + 1
         fval3 = f
      endif
      ierflg = 11
      if (cword(1:3).eq.'END')  ierflg = 10
      if (cword(1:3).eq.'RET')  ierflg = 12
      go to 5000
!                                        . . . . . . . . . . clear
2200  continue
      call mncler
      if (isw(5).ge.1)  write (minuit_fileout,'(A)')&
           & ' MINUIT MEMORY CLEARED. NO PARAMETERS NOW DEFINED.'
      go to 5000
!                                        . . . . . . . . . . help
2300  continue
!ccc      if (index(cword,'SHO').gt.0)  go to 700
!ccc      if (index(cword,'SET').gt.0)  go to 700
      kcol = 0
      do icol= 5,lk
         if (cword(icol:icol).eq.' ') cycle
         kcol = icol
         exit
      end do
      !
      if (kcol.eq.0)  then
         comd = '*   '
      else
         comd = cword(kcol:lk)
      endif
      call mnhelp(comd,minuit_fileout)
      go to 5000
!                                       . . . . . . . . . . mncontour
2400  continue
      epsi = 0.05 * up
      ke1 = word7(1)
      ke2 = word7(2)
      if (ke1.eq.0.and.npar.eq.2) then
         ke1 = nexofi(1)
         ke2 = nexofi(2)
      endif
      nptu = word7(3)
      if (nptu.le.0)  nptu=20
      if (nptu.gt.mxpt)  nptu = mxpt
      nfcnmx =  100*(nptu+5)*(npar+1)
      call minuit_mncont(ke1,ke2,nptu,xptu,yptu,ierrf)
      if (ierrf.lt.nptu) ierflg = 4
      if (ierrf.eq.-1)   ierflg = 3
      go to 5000
!                                      . . . . . . . . . . jump
2600  continue
      step = word7(1)
      if (step.le.zero)  step = 2.
      rno = 0.
      izero = 0
      do i= 1, npar
         call mnrn15(rno,izero)
         rno = 2.0*rno - 1.0
         x(i) = x(i) + rno*step*werr(i)
      end do
      call mninex(x)
      call mnamin
      call mnrset(0)
      go to 5000
!                                      . . . . . . . . . . blank line
3300  continue
      write (minuit_fileout,'(10X,A)') ' BLANK COMMAND IGNORED.'
      ierflg = 1
      go to 5000
!  . . . . . . . . obsolete commands     . . . . . . . . . . . . . .
!                                      . . . . . . . . . . covariance
3400  continue
      write (minuit_fileout, '(A)') ' THE "COVARIANCE" COMMAND IS OSBSOLETE.',&
     & ' THE COVARIANCE MATRIX IS NOW SAVED IN A DIFFERENT FORMAT',&
     & ' WITH THE "SAVE" COMMAND AND READ IN WITH:"SET COVARIANCE"'
      ierflg = 3
      go to 5000
!                                        . . . . . . . . . . printout
3500  continue
      cneway = 'SET PRINT '
      go to 3100
!                                        . . . . . . . . . . gradient
3600  continue
      cneway = 'SET GRAD  '
      go to 3100
!                                        . . . . . . . . . . matout
3700  continue
      cneway = 'SHOW COVAR'
      go to 3100
!                                        . . . . . . . . . error def
3800  continue
      cneway = 'SET ERRDEF'
      go to 3100
!                                        . . . . . . . . . . limits
3900  continue
      cneway = 'SET LIMITS'
      go to 3100
!                                        . . . . . . . . . . punch
4000  continue
      cneway = 'SAVE      '
!                                ....... come from obsolete commands
3100  write (minuit_fileout, 3101) cword,cneway
3101  format (' OBSOLETE COMMAND:',1X,A10,5X,'PLEASE USE:',1X,A10)
      cword = cneway
      if (cword.eq.'SAVE      ') go to 1500
      go to 700
!                                 . . . . . . . . . . . . . . . . . .
5000  return
    end subroutine minuit_mnexcm

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnseek
!c   performs a rough (but global) minimization by monte carlo search.
!c        each time a new minimum is found, the search area is shifted
!c        to be centered at the best value.  random points are chosen
!c        uniformly over a hypercube determined by current step sizes.
!c   the metropolis algorithm accepts a worse point with probability
!c      exp(-d/up), where d is the degradation.  improved points
!c      are of course always accepted.  actual steps are random
!c      multiples of the nominal steps (dirin).
!c
    integer(FINT), parameter :: twopi=2.0*3.141593
    real(FREAL), dimension(mni) :: xbest, xmid !(mni)
    !
    integer(FINT) :: mxfail, mxstep, ifail, nparx, ipar, iext
    integer(FINT) :: istep, iseed, ib, j
    !
    real(FREAL) :: alpha, rnum, rnum1, rnum2, flast, dxdi, ftry, bar
    !
    mxfail = word7(1)
    if (mxfail.le.0)  mxfail=100+20*npar
    mxstep = 10*mxfail
    if (amin.eq.undefi)  call mnamin
    alpha = word7(2)
    if (alpha.le.zero)  alpha=3.
    if (isw(5).ge.1)  write (minuit_fileout, 3) mxfail,mxstep,alpha
3   format (' MNSEEK: MONTE CARLO MINIMIZATION USING METROPOLIS',&
         & ' ALGORITHM'/' TO STOP AFTER',I6,' SUCCESSIVE FAILURES, OR',&
         & I7,' STEPS'/' MAXIMUM STEP SIZE IS',F9.3,' ERROR BARS.')
    cstatu= 'INITIAL  '
    if (isw(5).ge.2)  call mnprin(2,amin)
    cstatu = 'UNCHANGED '
    ifail = 0
    rnum = zero
    rnum1 = zero
    rnum2 = zero
    nparx = npar
    flast = amin
!              set up step sizes, starting values
    do ipar =  1, npar
       iext = nexofi(ipar)
       dirin(ipar) = 2.0*alpha*werr(ipar)
       if (nvarl(iext).gt.1)  then
!              parameter with limits
          call mndxdi(x(ipar),ipar,dxdi)
          if (dxdi.eq.zero)  dxdi=1.
          dirin(ipar) = 2.0*alpha*werr(ipar)/dxdi
          if (abs(dirin(ipar)).gt.twopi)  dirin(ipar)=twopi
       endif
       xmid(ipar) = x(ipar)
       xbest(ipar) = x(ipar)
    end do
!                              search loop
    do istep= 1, mxstep
       if (ifail.ge.mxfail)  go to 600
       do ipar= 1, npar
          call mnrn15(rnum1,iseed)
          call mnrn15(rnum2,iseed)
          x(ipar) = xmid(ipar) + 0.5*(rnum1+rnum2-1.)*dirin(ipar)
       end do
       call mninex(x)
       call minuit_fcn(nparx,gin,ftry,u,4)
       nfcn = nfcn + 1
       if (ftry.lt.flast)  then
          if (ftry.lt.amin)  then
             cstatu = 'IMPROVEMNT'
             amin = ftry
             do ib= 1, npar
                xbest(ib) = x(ib)
             end do
             ifail = 0
             if (isw(5).ge.2) call mnprin(2,amin)
          endif
          go to 300
       else
          ifail = ifail + 1
!                   metropolis algorithm
          bar = (amin-ftry)/up
          call mnrn15(rnum,iseed)
          if (bar.lt.log(rnum))  cycle
       endif
!                    accept new point, move there
300    continue
       do j= 1, npar
          xmid(j) = x(j)
       end do
       flast = ftry
    end do
!                               end search loop
600 continue
    if (isw(5).gt.1) write (minuit_fileout,601) ifail
601 format(' MNSEEK:',I5,' SUCCESSIVE UNSUCCESSFUL TRIALS.')
    do ib= 1, npar
       x(ib) = xbest(ib)
    end do
    call mninex(x)
    if (isw(5).ge.1)  call mnprin(2,amin)
    if (isw(5).eq.0)  call mnprin(0,amin)
    return
  end subroutine mnseek

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnprin(inkode,fval)
!c        prints the values of the parameters at the time of the call.
!c        also prints other relevant information such as function value,
!c        estimated distance to minimum, parameter errors, step sizes.
!c
!         according to the value of ikode, the printout is:
!    ikode=inkode= 0    only info about function value
!                  1    parameter values, errors, limits
!                  2    values, errors, step sizes, internal values
!                  3    values, errors, step sizes, first derivs.
!                  4    values, parabolic errors, minos errors
!    when inkode=5, mnprin chooses ikode=1,2, or 3, according to isw(2)
!
!
    integer(FINT) :: inkode
    real(FREAL) :: fval
    !
    integer(FINT) :: i
    character*14 colhdu(6),colhdl(6), cx2,cx3
    character*11 cnambf, cblank
    character  chedm*10, cheval*15
    character(14), parameter :: cgetx='PLEASE GET X..'
    !
    integer(FINT) :: ikode, k, nc, m, ntrail, ic, lbl, nadd, ncol, kk, l
    real(FREAL) :: dcmax, dc, x1, x2, x3
    !
    data cblank/'          '/
!
    if (nu.eq.0)  then
       write (minuit_fileout,'(A)') ' THERE ARE CURRENTLY NO PARAMETERS DEFINED'
       go to 700
    endif
!                  get value of ikode based in inkode, isw(2)
    ikode = inkode
    if (inkode.eq.5) then
       ikode = isw(2)+1
       if (ikode.gt.3)  ikode=3
    endif
!                  set 'default' column headings
    do k= 1, 6
       colhdu(k) = 'UNDEFINED'
       colhdl(k) = 'COLUMN HEAD'
    end do
!              print title if minos errors, and title exists.
    if (ikode.eq.4.and.ctitl.ne.cundef)&
         &            write (minuit_fileout,'(/A,A)')  ' MINUIT TASK: ',ctitl
!              report function value and status
    if (fval.eq.undefi) then
       cheval = ' UNKNOWN       '
    else
       write (cheval,'(G15.7)') fval
    endif
    if (edm.eq.bigedm) then
       chedm = ' UNKNOWN  '
    else
       write (chedm, '(E10.2)') edm
    endif
    nc = nfcn-nfcnfr
    write (minuit_fileout,905)  cheval,cfrom,cstatu,nc,nfcn
905 format (/' FCN=',A,' FROM ',A8,'  STATUS=',A10,I6,' CALLS',&
         &         I9,' TOTAL')
    m = isw(2)
    if (m.eq.0.or.m.eq.2.or.dcovar.eq.zero) then
       write (minuit_fileout,907) chedm,istrat,covmes(m)
907    format (21X,'EDM=',A,'    STRATEGY=',I2,6X,A)
    else
       dcmax = 1.
       dc = min(dcovar,dcmax) * 100.0_FREAL
       write (minuit_fileout,908) chedm,istrat,dc
908    format (21X,'EDM=',A,'  STRATEGY=',I1,'  ERROR MATRIX',&
            &     ' UNCERTAINTY=',F5.1,'%')
    endif
!
    if (ikode.eq.0)  go to 700
!               find longest name (for rene!)
    ntrail = 10
    do i= 1, nu
       if (nvarl(i).lt.0) cycle
       do ic= 10,1,-1
          if (cpnam(i)(ic:ic).ne.' ') go to 16
       end do
       ic = 1
16     lbl = 10-ic
       if (lbl.lt.ntrail)  ntrail=lbl
    end do
    nadd = ntrail/2 + 1
    if (ikode.eq.1)  then
       colhdu(1) = '              '
       colhdl(1) = '      ERROR   '
       colhdu(2) = '      PHYSICAL'
       colhdu(3) = ' LIMITS       '
       colhdl(2) = '    NEGATIVE  '
       colhdl(3) = '    POSITIVE  '
    endif
    if (ikode.eq.2)  then
       colhdu(1) = '              '
       colhdl(1) = '      ERROR   '
       colhdu(2) = '    INTERNAL  '
       colhdl(2) = '    STEP SIZE '
       colhdu(3) = '    INTERNAL  '
       colhdl(3) = '      VALUE   '
    endif
    if (ikode.eq.3)  then
       colhdu(1) = '              '
       colhdl(1) = '      ERROR   '
       colhdu(2) = '       STEP   '
       colhdl(2) = '       SIZE   '
       colhdu(3) = '      FIRST   '
       colhdl(3) = '   DERIVATIVE '
    endif
    if (ikode.eq.4)  then
       colhdu(1) = '    PARABOLIC '
       colhdl(1) = '      ERROR   '
       colhdu(2) = '        MINOS '
       colhdu(3) = 'ERRORS        '
       colhdl(2) = '   NEGATIVE   '
       colhdl(3) = '   POSITIVE   '
    endif
!
    if (ikode.ne.4)  then
       if (isw(2).lt.3) colhdu(1)='  APPROXIMATE '
       if (isw(2).lt.1) colhdu(1)=' CURRENT GUESS'
    endif
    ncol = 3
    write (minuit_fileout, 910) (colhdu(kk),kk=1,ncol)
    write (minuit_fileout, 911) (colhdl(kk),kk=1,ncol)
910 format (/'  EXT PARAMETER ',     13X       ,6A14)
911 format ( '  NO.   NAME    ','    VALUE    ',6A14)
!
!                                        . . . loop over parameters . .
    do i= 1, nu
       if (nvarl(i).lt.0) cycle
       l = niofex(i)
       cnambf = cblank(1:nadd)//cpnam(i)
       if (l.eq.0)  go to 55
       !              variable parameter.
       x1 = werr(l)
       cx2 = cgetx
       cx3 = cgetx
       if (ikode.eq.1) then
          if (nvarl(i).le.1) then
             write (minuit_fileout, 952)  i,cnambf,u(i),x1
             cycle
          else
             x2 = alim(i)
             x3 = blim(i)
          endif
       endif
       if (ikode.eq.2) then
          x2 = dirin(l)
          x3 = x(l)
       endif
       if (ikode.eq.3) then
          x2 = dirin(l)
          x3 = grd(l)
          if (nvarl(i).gt.1.and.abs(cos(x(l))).lt.0.001)&
               &      cx3 = '** AT LIMIT **'
       endif
       if (ikode.eq.4) then
          x2 = ern(l)
          if (x2.eq.zero)   cx2=' '
          if (x2.eq.undefi) cx2='   AT LIMIT   '
          x3 = erp(l)
          if (x3.eq.zero)   cx3=' '
          if (x3.eq.undefi) cx3='   AT LIMIT   '
       endif
       if (cx2.eq.cgetx) write (cx2,'(G14.5)') x2
       if (cx3.eq.cgetx) write (cx3,'(G14.5)') x3
       write (minuit_fileout,952)   i,cnambf,u(i),x1,cx2,cx3
952    format (i4,1x,a11,2g14.5,2a)
       !               check if parameter is at limit
       if (nvarl(i).le.1.or.ikode.eq.3) cycle
       if (abs(cos(x(l))).lt.0.001)  write (minuit_fileout,1004)
1004   format (1h ,32x,42hwarning -   - above parameter is at limit.)
       cycle
!
!                                print constant or fixed parameter.
55     continue
       colhdu(1) = '   CONSTANT   '
       if (nvarl(i).gt.0)  colhdu(1) = '     FIXED    '
       if (nvarl(i).eq.4.and.ikode.eq.1) then
          write (minuit_fileout,'(I4,1X,A11,G14.5,A,2G14.5)')&
               &     i,cnambf,u(i),colhdu(1),alim(i),blim(i)
       else
          write (minuit_fileout,'(I4,1X,A11,G14.5,A)')  i,cnambf,u(i),colhdu(1)
       endif
    end do
!
    if (up.ne.updflt)  write (minuit_fileout,'(31X,A,G10.3)') 'ERR DEF=',up
700 continue
    return
  end subroutine mnprin

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnamin
!c        called  from many places.  initializes the value of amin by
!c        calling the user function. prints out the function value and
!c        parameter values if print flag value is high enough.
!c
    integer(FINT) :: nparx
    real(FREAL) :: fnew
    nparx = npar
    if (isw(5).ge.1) write (minuit_fileout,'(/A,A)') ' FIRST CALL TO ',&
         & 'USER FUNCTION AT NEW START POINT, WITH IFLAG=4.'
    call mnexin(x)
    call minuit_fcn(nparx,gin,fnew,u,4)
    nfcn = nfcn + 1
    amin = fnew
    edm = bigedm
    return
  end subroutine mnamin

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mndxdi(pint,ipar,dxdi)
!c        calculates the transformation factor between external and
!c        internal parameter values.     this factor is one for
!c        parameters which are not limited.     called from mnemat.
    !
    integer(FINT), intent(in) :: ipar
    real(FREAL), intent(in) :: pint
    real(FREAL), intent(out) :: dxdi
    !
    integer(FINT) :: i
    !
    i = nexofi(ipar)
    dxdi = 1.0
    if (nvarl(i).gt.1)&
         &      dxdi = 0.5 *abs((blim(i)-alim(i)) * cos(pint))
    return
  end subroutine mndxdi

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnrn15(val,inseed)
!         this is a super-portable random number generator.
!         it should not overflow on any 32-bit machine.
!         the cycle is only ~10**9, so use with care!
!         note especially that val must not be undefined on input.
!                    set default starting seed
    !
    real(FREAL), parameter :: three=3.0
    !
    real(FREAL) :: val
    integer(FINT) :: inseed, iseed, k
    data iseed/12345/
    if (val.eq.three)  go to 100
!
    inseed = iseed
    k = iseed/53668
    iseed = 40014*(iseed-k*53668) - k*12211
    if (iseed.lt.0) iseed = iseed + 2147483563
    val = real(iseed) * 4.656613e-10
    return
!               "entry" to set seed, flag is val=3.
100 iseed = inseed
    return
  end subroutine mnrn15

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mninex(pint)
!c        transforms from internal coordinates (pint) to external
!c        parameters (u).   the minimizing routines which work in
!c        internal coordinates call this routine before calling fcn.
    !
    real(FREAL), dimension(:) :: pint !(*)
    !
    integer(FINT) :: j,i
    !
    do j= 1, npar
       i = nexofi(j)
       if (nvarl(i).eq.1) then
          u(i) = pint(j)
       else
          u(i) = alim(i) + 0.5*(sin(pint(j)) +1.0) * (blim(i)-alim(i))
       endif
    end do
    return
  end subroutine mninex

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnsimp
!c        performs a minimization using the simplex method of nelder
!c        and mead (ref. -- comp. j. 7,308 (1965)).
!c
    real(FREAL), dimension(mni+1) :: y !(mni+1)
    real(FREAL) :: alpha, beta, gamma, rhomin, rhomax
    !
    integer(FINT) :: npfn, nparp1, nparx, jl, kg, ns, nf, k, jh, j, i
    real(FREAL) :: rho1, rho2, wg, dxdi, dmin, ynpp1, absmin, aming
    !real(FREAL) :: bestx, f, sig2, ncycl, pb, ystar, ystst, y1, y2
    real(FREAL) :: bestx, f, sig2, pb, ystar, ystst, y1, y2
    integer(FINT) :: ncycl
    real(FREAL) :: rho, yrho, jhold, ypbar
    !
    data alpha,beta,gamma,rhomin,rhomax / 1.0, 0.5, 2.0, 4.0, 8.0/
    !
    if (npar.le.0)  return
    if (amin.eq.undefi)  call mnamin
    cfrom = 'SIMPLEX '
    nfcnfr = nfcn
    cstatu= 'UNCHANGED '
    npfn=nfcn
    nparp1=npar+1
    nparx = npar
    rho1 = 1.0 + alpha
    rho2 = rho1 + alpha*gamma
    wg = 1.0/float(npar)
    if (isw(5).ge.0) write(minuit_fileout,100) epsi
100 format(' START SIMPLEX MINIMIZATION.    CONVERGENCE WHEN EDM .LT.'&
         &,E10.2 )
    do i= 1, npar
       dirin(i) = werr(i)
       call mndxdi(x(i),i,dxdi)
       if (dxdi.ne.zero) dirin(i)=werr(i)/dxdi
       dmin = epsma2*abs(x(i))
       if (dirin(i).lt.dmin)  dirin(i)=dmin
    end do
!**       choose the initial simplex using single-parameter searches
1   continue
    ynpp1 = amin
    jl = nparp1
    y(nparp1) = amin
    absmin = amin
    do i= 1, npar
       aming = amin
       pbar(i) = x(i)
       bestx = x(i)
       kg = 0
       ns = 0
       nf = 0
4      x(i) = bestx + dirin(i)
       call mninex(x)
       call minuit_fcn(nparx,gin, f, u, 4)
       nfcn = nfcn + 1
       if (f.lt.aming)  go to 6
!         failure
       if (kg.eq.1)  go to 8
       kg = -1
       nf = nf + 1
       dirin(i) = dirin(i) * (-0.4)
       if (nf.lt.3)  go to 4
!         stop after three failures
       bestx = x(i)
       dirin(i) = dirin(i) * 3.0
       aming = f
       go to 8
!
!         success
6      bestx = x(i)
       dirin(i) = dirin(i) * 3.0
       aming = f
       cstatu= 'PROGRESS  '
       kg = 1
       ns = ns + 1
       if (ns.lt.6)  go to 4
!
!         3 failures or 6 successes or
!         local minimum found in ith direction
8      y(i) = aming
       if (aming.lt.absmin)  jl = i
       if (aming.lt.absmin)  absmin = aming
       x(i) = bestx
       do k= 1, npar
          p(k,i) = x(k)
       end do
    end do
    jh = nparp1
    amin=y(jl)
    call mnrazz(ynpp1,pbar,y,jh,jl)
    do i= 1, npar
       x(i) = p(i,jl)
    end do
    call mninex(x)
    if (isw(5).ge.1)  call mnprin(5,amin)
    edm = bigedm
    sig2 = edm
    ncycl=0
!                                        . . . . .  start main loop
50  continue
    if (sig2.lt.epsi.and.edm.lt.epsi)     go to 76
    sig2 = edm
    if ((nfcn-npfn).gt.nfcnmx)  go to 78
!         calculate new point * by reflection
    do i= 1, npar
       pb = 0.
       do j= 1, nparp1
          pb = pb + wg * p(i,j)
       end do
       pbar(i) = pb - wg * p(i,jh)
       pstar(i)=(1.+alpha)*pbar(i)-alpha*p(i,jh)
    end do
    call mninex(pstar)
    call minuit_fcn(nparx,gin,ystar,u,4)
    nfcn=nfcn+1
    if(ystar.ge.amin) go to 70
!         point * better than jl, calculate new point **
    cstatu = 'PROGRESS  '
    do i=1,npar
       pstst(i)=gamma*pstar(i)+(1.-gamma)*pbar(i)
    end do
    call mninex(pstst)
    call minuit_fcn(nparx,gin,ystst,u,4)
    nfcn=nfcn+1
!         try a parabola through ph, pstar, pstst.  min = prho
    y1 = (ystar-y(jh)) * rho2
    y2 = (ystst-y(jh)) * rho1
    rho = 0.5 * (rho2*y1 -rho1*y2) / (y1 -y2)
    if (rho.lt.rhomin)  go to 66
    if (rho.gt.rhomax)  rho = rhomax
    do i= 1, npar
       prho(i) = rho*pbar(i) + (1.0-rho)*p(i,jh)
    end do
    call mninex(prho)
    call minuit_fcn(nparx,gin,yrho, u,4)
    nfcn = nfcn + 1
    if (yrho.lt.amin)     cstatu = 'PROGRESS  '
    if (yrho.lt.y(jl).and.yrho.lt.ystst)  go to 65
    if (ystst.lt.y(jl))  go to 67
    if (yrho.gt.y(jl))  go to 66
!         accept minimum point of parabola, prho
65  call mnrazz (yrho,prho,y,jh,jl)
    go to 68
66  if (ystst.lt.y(jl))  go to 67
    call mnrazz(ystar,pstar,y,jh,jl)
    go to 68
67  call mnrazz(ystst,pstst,y,jh,jl)
68  ncycl=ncycl+1
    if (isw(5).lt.2)  go to 50
    if (isw(5).ge.3.or.mod(ncycl, 10).eq.0) call mnprin(5,amin)
    go to 50
!         point * is not as good as jl
70  if (ystar.ge.y(jh))  go to 73
    jhold = jh
    call mnrazz(ystar,pstar,y,jh,jl)
    if (jhold.ne.jh)  go to 50
!         calculate new point **
73  do i=1,npar
       pstst(i)=beta*p(i,jh)+(1.-beta)*pbar(i)
    end do
    call mninex (pstst)
    call minuit_fcn(nparx,gin,ystst,u,4)
    nfcn=nfcn+1
    if(ystst.gt.y(jh)) go to 1
!     point ** is better than jh
    if (ystst.lt.amin)     cstatu = 'PROGRESS  '
    if (ystst.lt.amin)  go to 67
    call mnrazz(ystst,pstst,y,jh,jl)
    go to 50
!                                        . . . . . .  end main loop
76  if (isw(5).ge.0)  write(minuit_fileout,'(A)')&
         &                    ' SIMPLEX MINIMIZATION HAS CONVERGED.'
    isw(4) = 1
    go to 80
78  if (isw(5).ge.0)  write(minuit_fileout,'(A)')&
         &                    ' SIMPLEX TERMINATES WITHOUT CONVERGENCE.'
    cstatu= 'CALL LIMIT'
    isw(4) = -1
    isw(1) = 1
80  do i=1,npar
       pb = 0.
       do j=1,nparp1
          pb = pb + wg * p(i,j)
       end do
       pbar(i) = pb - wg * p(i,jh)
    end do
    call mninex(pbar)
    call minuit_fcn(nparx,gin,ypbar,u,4)
    nfcn=nfcn+1
    if (ypbar.lt.amin)  call mnrazz(ypbar,pbar,y,jh,jl)
    call mninex(x)
    if (nfcnmx+npfn-nfcn.lt.3*npar)  go to 90
    if (edm.gt.2.0*epsi)  go to 1
90  if (isw(5).ge.0)  call mnprin(5, amin)
    return
  end subroutine mnsimp

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnrazz(ynew,pnew,y,jh,jl)
!c        called only by mnsimp (and mnimpr) to add a new point
!c        and remove an old one from the current simplex, and get the
!c        estimated distance to minimum.
!c
    !
    real(FREAL), dimension(:) :: pnew, y !(*)
    real(FREAL) :: ynew
    integer(FINT) :: jh, jl
    !
    integer(FINT) :: i
    integer(FINT) :: nparp1, j
    real(FREAL) :: pbig, plit
    !
    do i=1,npar
       p(i,jh) = pnew(i)
    end do
    y(jh)=ynew
    if(ynew.lt.amin) then
       do i=1,npar
          x(i) = pnew(i)
       end do
       call mninex(x)
       amin = ynew
       cstatu = 'PROGRESS  '
       jl=jh
    endif
    jh = 1
    nparp1 = npar+1
20  do j=2,nparp1
       if (y(j).gt.y(jh))  jh = j
    end do
    edm = y(jh) - y(jl)
    if (edm.le.zero)  go to 45
    do i= 1, npar
       pbig = p(i,1)
       plit = pbig
       do j= 2, nparp1
          if (p(i,j).gt.pbig)  pbig = p(i,j)
          if (p(i,j).lt.plit)  plit = p(i,j)
       end do
       dirin(i) = pbig - plit
    end do
40  return
45  write (minuit_fileout, 1000)  npar
    go to 40
1000 format ('   FUNCTION VALUE DOES NOT SEEM TO DEPEND ON ANY OF THE',&
          &    I3,' VARIABLE PARAMETERS.' /10X,'VERIFY THAT STEP SIZES ARE',&
          &    ' BIG ENOUGH AND CHECK FCN LOGIC.'/1X,79(1H*)/1X,79(1H*)/)
  end subroutine mnrazz

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnmigr
!c        performs a local function minimization using basically the
!c        method of davidon-fletcher-powell as modified by fletcher
!c        ref. -- fletcher, comp.j. 13,317 (1970)   "switching method"
!c
    !
    !
    integer(FINT) :: i
    real(FREAL), dimension(mni) :: gs, step,  xxs, flnu, vg !(mni)
    logical ldebug
    real(FREAL), parameter :: toler=0.05
    !
    integer(FINT) :: nfcnmg, iswtr, npfn, nparx, lenv, nrstrt, npsdf
    integer(FINT) :: lined2, iext, ndex, j, kk, iter, m, n
    real(FREAL) :: rhotol, fzero, gdel, fs, gssq, ri
    real(FREAL) :: gvg, delgam, gdgssq, vgi, gami, dsum, vsum, d
    !
    if (npar.le.0)  return
    if (amin.eq.undefi)  call mnamin
    ldebug = (idbg(4).ge.1)
    cfrom = 'MIGRAD  '
    nfcnfr = nfcn
    nfcnmg = nfcn
    cstatu= 'INITIATE  '
    iswtr = isw(5) - 2*itaur
    npfn = nfcn
    nparx = npar
    lenv = npar*(npar+1)/2
    nrstrt = 0
    npsdf = 0
    lined2 = 0
    isw(4) = -1
    rhotol = 1.0e-3*apsi
    if (iswtr.ge.1)  write (minuit_fileout,470) istrat,rhotol
470 format (' START MIGRAD MINIMIZATION.  STRATEGY',I2,&
         &'.  CONVERGENCE WHEN EDM .LT.',E9.2)
!                                           initialization strategy
    if (istrat.lt.2.or.isw(2).ge.3)  go to 2
!                                come (back) here to restart completely
1   continue
    if (nrstrt.gt.istrat)  then
       cstatu= 'FAILED    '
       isw(4) = -1
       go to 230
    endif
!                                      . get full covariance and gradien
    call mnhess
    call mnwerr
    npsdf = 0
    if (isw(2).ge.1)  go to 10
!                                        . get gradient at start point
2   continue
    call mninex(x)
    if (isw(3).eq.1) then
       call minuit_fcn(nparx,gin,fzero,u,2)
       nfcn = nfcn + 1
    endif
    call mnderi
    if (isw(2).ge.1)  go to 10
!                                   sometimes start with diagonal matrix
    do i= 1, npar
       xxs(i) = x(i)
       step(i) = zero
    end do
!                           do line search if second derivative negative
    lined2 = lined2 + 1
    if (lined2.lt.(istrat+1)*npar) then
       do i= 1, npar
          if (g2(i).gt.zero) cycle
          step(i) = -sign(gstep(i),grd(i))
          gdel = step(i)*grd(i)
          fs = amin
          call mnline(xxs,fs,step,gdel,toler)
          call mnwarn('D','MNMIGR','NEGATIVE G2 LINE SEARCH')
          iext = nexofi(i)
          if (ldebug) write (minuit_fileout,'(A,I3,2G13.3)')&
               &    ' NEGATIVE G2 LINE SEARCH, PARAM ',iext,fs,amin
          go to 2
       end do
    endif
!                           make diagonal error matrix
    do i=1,npar
       ndex = i*(i-1)/2
       do j=1,i-1
          ndex = ndex + 1
          vhmat(ndex) = 0.
       end do
       ndex = ndex + 1
       if (g2(i).le.zero)  g2(i) = 1.
       vhmat(ndex) = 2./g2(i)
    end do
    dcovar = 1.
    if (ldebug) write (minuit_fileout,'(A,A/(1X,10G10.2))') ' DEBUG MNMIGR,',&
         &  ' STARTING MATRIX DIAGONAL,  VHMAT=', (vhmat(kk),kk=1,lenv)
!                                         ready to start first iteration
10  continue
    nrstrt = nrstrt + 1
    if (nrstrt.gt.istrat+1)  then
       cstatu= 'FAILED    '
       go to 230
    endif
    fs = amin
!                                        . . . get edm and set up loop
    edm = 0.
    do i= 1, npar
       gs(i) = grd(i)
       xxs(i) = x(i)
       ndex = i*(i-1)/2
       do j= 1, i-1
          ndex = ndex + 1
          edm = edm + gs(i)*vhmat(ndex)*gs(j)
       end do
       ndex = ndex + 1
       edm = edm + 0.5 * gs(i)**2 *vhmat(ndex)
    end do
    edm = edm * 0.5 * (1.0+3.0*dcovar)
    if (edm.lt.zero)  then
       call mnwarn('W','MIGRAD','STARTING MATRIX NOT POS-DEFINITE.')
       isw(2) = 0
       dcovar = 1.
       go to 2
    endif
    if (isw(2).eq.0)  edm=bigedm
    iter = 0
    call mninex(x)
    call mnwerr
    if (iswtr.ge.1)  call mnprin(3,amin)
    if (iswtr.ge.2)  call mnmatu(0)
!                                        . . . . .  start main loop
24  continue
    if (nfcn-npfn.ge.nfcnmx)  go to 190
    gdel = 0.
    gssq = 0.
    do i=1,npar
       ri = 0.
       gssq = gssq + gs(i)**2
       do j=1,npar
          m = max(i,j)
          n = min(i,j)
          ndex = m*(m-1)/2 + n
          ri = ri + vhmat(ndex) *gs(j)
       end do
       step(i) = -0.5*ri
       gdel = gdel + step(i)*gs(i)
    end do
    if (gssq.eq.zero)  then
       call mnwarn('D','MIGRAD',&
            &             ' FIRST DERIVATIVES OF FCN ARE ALL ZERO')
       go to 300
    endif
!                 if gdel positive, v not posdef
    if (gdel.ge.zero)  then
       call mnwarn('D','MIGRAD',' NEWTON STEP NOT DESCENT.')
       if (npsdf.eq.1)  go to 1
       call mnpsdf
       npsdf = 1
       go to 24
    endif
!                                        . . . . do line search
    call mnline(xxs,fs,step,gdel,toler)
    if (amin.eq.fs) go to 200
    cfrom  = 'MIGRAD  '
    nfcnfr = nfcnmg
    cstatu= 'PROGRESS  '
!                                        . get gradient at new point
    call mninex(x)
    if (isw(3).eq.1) then
       call minuit_fcn(nparx,gin,fzero,u,2)
       nfcn = nfcn + 1
    endif
    call mnderi
!                                         . calculate new edm
    npsdf = 0
81  edm = 0.
    gvg = 0.
    delgam = 0.
    gdgssq = 0.
    do i= 1, npar
       ri = 0.
       vgi = 0.
       do j= 1, npar
          m = max(i,j)
          n = min(i,j)
          ndex = m*(m-1)/2 + n
          vgi = vgi + vhmat(ndex)*(grd(j)-gs(j))
          ri  =  ri + vhmat(ndex)* grd(j)
       end do
       vg(i) = vgi*0.5
       gami = grd(i) - gs(i)
       gdgssq = gdgssq + gami**2
       gvg = gvg + gami*vg(i)
       delgam = delgam + dirin(i)*gami
       edm = edm + grd(i)*ri*0.5
    end do
    edm = edm * 0.5 * (1.0 + 3.0*dcovar)
!                          . if edm negative,  not positive-definite
    if (edm.lt.zero.or.gvg.le.zero)  then
       call mnwarn('D','MIGRAD','NOT POS-DEF. EDM OR GVG NEGATIVE.')
       cstatu = 'NOT POSDEF'
       if (npsdf.eq.1)  go to 230
       call mnpsdf
       npsdf = 1
       go to 81
    endif
!                            print information about this iteration
    iter = iter + 1
    if (iswtr.ge.3.or.(iswtr.eq.2.and.mod(iter,10).eq.1)) then
       call mnwerr
       call mnprin(3,amin)
    endif
    if (gdgssq.eq.zero)  call mnwarn('D','MIGRAD',&
         &           'NO CHANGE IN FIRST DERIVATIVES OVER LAST STEP')
    if (delgam.lt.zero) call mnwarn('D','MIGRAD',&
         &          'FIRST DERIVATIVES INCREASING ALONG SEARCH LINE')
!                                        .  update covariance matrix
    cstatu = 'IMPROVEMNT'
    if (ldebug) write (minuit_fileout,'(A,(1X,10G10.3))') ' VHMAT 1 =',&
         &             (vhmat(kk),kk=1,10)
    dsum = 0.
    vsum = 0.
    do i=1, npar
       do j=1, i
          d = dirin(i)*dirin(j)/delgam - vg(i)*vg(j)/gvg
          dsum = dsum + abs(d)
          ndex = i*(i-1)/2 + j
          vhmat(ndex) = vhmat(ndex) + 2.0*d
          vsum = vsum + abs(vhmat(ndex))
       end do
    end do
!                smooth local fluctuations by averaging dcovar
    dcovar = 0.5*(dcovar + dsum/vsum)
    if (iswtr.ge.3.or.ldebug) write (minuit_fileout,'(A,F5.1,A)')&
         &      ' RELATIVE CHANGE IN COV. MATRIX=',dcovar*100.,'%'
    if (ldebug) write (minuit_fileout,'(A,(1X,10G10.3))') ' VHMAT 2 =',&
         &             (vhmat(kk),kk=1,10)
    if (delgam.le.gvg)  go to 135
    do i= 1, npar
       flnu(i) = dirin(i)/delgam - vg(i)/gvg
    end do
    do i= 1, npar
       do j= 1, i
          ndex = i*(i-1)/2 + j
          vhmat(ndex) = vhmat(ndex) + 2.0*gvg*flnu(i)*flnu(j)
       end do
    end do
135 continue
!                                              and see if converged
    if (edm.lt.0.1*rhotol)  go to 300
!                                    if not, prepare next iteration
    do i= 1, npar
       xxs(i) = x(i)
       gs(i) = grd(i)
    end do
    fs = amin
    if (isw(2).eq.0.and.dcovar.lt.0.5 )  isw(2) = 1
    if (isw(2).eq.3.and.dcovar.gt.0.1 )  isw(2) = 1
    if (isw(2).eq.1.and.dcovar.lt.0.05)  isw(2) = 3
    go to 24
!                                        . . . . .  end main loop
!                                         . . call limit in mnmigr
190 isw(1) = 1
    if (isw(5).ge.0)&
         &     write (minuit_fileout,'(A)')  ' CALL LIMIT EXCEEDED IN MIGRAD.'
    cstatu = 'CALL LIMIT'
    go to 230
!                                         . . fails to improve . .
200 if (iswtr.ge.1)  write (minuit_fileout,'(A)')&
         &           ' MIGRAD FAILS TO FIND IMPROVEMENT'
    do i= 1, npar
       x(i) = xxs(i)
    end do
    if (edm.lt.rhotol)  go to 300
    if (edm.lt.abs(epsma2*amin))  then
       if (iswtr.ge.0)  write (minuit_fileout, '(A)')&
            &      ' MACHINE ACCURACY LIMITS FURTHER IMPROVEMENT.'
       go to 300
    endif
    if (istrat.lt.1)  then
       if (isw(5).ge.0) write (minuit_fileout, '(A)')&
            &    ' MIGRAD FAILS WITH STRATEGY=0.   WILL TRY WITH STRATEGY=1.'
       istrat = 1
    endif
    go to 1
!                                         . . fails to converge
230 if (iswtr.ge.0)  write (minuit_fileout,'(A)')&
         &    ' MIGRAD TERMINATED WITHOUT CONVERGENCE.'
    if (isw(2).eq.3)  isw(2) = 1
    isw(4) = -1
    go to 400
!                                         . . apparent convergence
300 if (iswtr.ge.0) write(minuit_fileout,'(/A)')&
         &   ' MIGRAD MINIMIZATION HAS CONVERGED.'
    if (itaur.eq.0) then
       if (istrat.ge.2.or.(istrat.eq.1.and.isw(2).lt.3)) then
          if (isw(5).ge.0)  write (minuit_fileout, '(/A)')&
               &      ' MIGRAD WILL VERIFY CONVERGENCE AND ERROR MATRIX.'
          call mnhess
          call mnwerr
          npsdf = 0
          if (edm.gt.rhotol) go to 10
       endif
    endif
    cstatu='CONVERGED '
    isw(4) = 1
!                                           come here in any case
400 continue
    cfrom = 'MIGRAD  '
    nfcnfr = nfcnmg
    call  mninex(x)
    call mnwerr
    if (iswtr.ge.0)  call mnprin (3,amin)
    if (iswtr.ge.1)  call mnmatu(1)
    return
  end subroutine mnmigr

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnhess
!c        calculates the full second-derivative matrix of fcn
!c        by taking finite differences. when calculating diagonal
!c        elements, it may iterate so that step size is nearly that
!c        which gives function change= up/10. the first derivatives
!c        of course come as a free side effect, but with a smaller
!c        step size in order to obtain a known accuracy.
!c
    !
    real(FREAL), dimension(mni) :: yy !(mni)
    logical ldebug
    character cbf1*22
    !
    integer(FINT) :: ncyc, npard, nparx, npar2, idrv, id, iext, i
    integer(FINT) :: icyc, multpy, ndex, j, ifail
    real(FREAL) :: tlrstp, tlrg2, fs1, df, aimsag, wint, dxdi, xtf, dmin, d
    real(FREAL) :: fs2, sag, g2bfor, dlast, stpinm, xti, xtj, elem, ztemp, g2i
    !
!
    ldebug = (idbg(3).ge.1)
    if (amin.eq.undefi)  call mnamin
    if (istrat.le.0) then
       ncyc = 3
       tlrstp = 0.5
       tlrg2  = 0.1
    else if (istrat.eq.1) then
       ncyc = 5
       tlrstp = 0.3
       tlrg2  = 0.05
    else
       ncyc = 7
       tlrstp = 0.1
       tlrg2  = 0.02
    endif
    if (isw(5).ge.2.or.ldebug)  write (minuit_fileout,'(A)')&
         &   '   START COVARIANCE MATRIX CALCULATION.'
    cfrom = 'HESSE   '
    nfcnfr = nfcn
    cstatu= 'OK        '
    npard = npar
!                 make sure starting at the right place
    call mninex(x)
    nparx = npar
    call minuit_fcn(nparx,gin,fs1,u,4)
    nfcn = nfcn + 1
    if (fs1.ne.amin) then
       df = amin - fs1
       write (cbf1(1:12),'(G12.3)') df
       call mnwarn('D','MNHESS',&
            &       'FUNCTION VALUE DIFFERS FROM AMIN BY '//cbf1(1:12) )
    endif
    amin = fs1
    if (ldebug) write (minuit_fileout,'(A,A)') ' PAR D   GSTEP          ',&
         &' D          G2         GRD         SAG    '
!                                        . . . . . . diagonal elements .
!         isw(2) = 1 if approx, 2 if not posdef, 3 if ok
!         aimsag is the sagitta we are aiming for in second deriv calc.
    aimsag = sqrt(epsma2)*(abs(amin)+up)
!         zero the second derivative matrix
    npar2 = npar*(npar+1)/2
    do i= 1,npar2
       vhmat(i) = 0.
    end do
!
!         loop over variable parameters for second derivatives
    idrv = 2
    do id= 1, npard
       i = id + npar - npard
       iext = nexofi(i)
       if (g2(i).eq.zero) then
          write (cbf1(1:4),'(I4)') iext
          call mnwarn('W','HESSE',&
               &      'SECOND DERIVATIVE ENTERS ZERO, PARAM '//cbf1(1:4) )
          wint = werr(i)
          if (nvarl(iext).gt.1) then
             call mndxdi(x(i),i,dxdi)
             if (abs(dxdi).lt..001) then
                wint = .01
             else
                wint = wint/abs(dxdi)
             endif
          endif
          g2(i) = up/wint**2
       endif
       xtf = x(i)
       dmin = 8.*epsma2*abs(xtf)
!
!                               find step which gives sagitta = aimsag
       d = abs(gstep(i))
       do icyc= 1, ncyc
!                               loop here only if sag=0.
          do multpy= 1, 5
!           take two steps
             x(i) = xtf + d
             call mninex(x)
             nparx = npar
             call minuit_fcn(nparx,gin,fs1,u,4)
             nfcn = nfcn + 1
             x(i) = xtf - d
             call mninex(x)
             call minuit_fcn(nparx,gin,fs2,u,4)
             nfcn = nfcn + 1
             x(i) = xtf
             sag = 0.5*(fs1+fs2-2.0*amin)
             if (sag.ne.zero) go to 30
             if (gstep(i).lt.zero) then
                if (d.ge..5) exit
                d = 10.*d
                if (d.gt.0.5)  d = 0.51
                cycle
             endif
             d = 10.*d
          end do
          write (cbf1(1:4),'(I4)') iext
          call mnwarn('W','HESSE',&
     &      'SECOND DERIVATIVE ZERO FOR PARAMETER'//cbf1(1:4) )
          go to 390
!                             sag is not zero
30        g2bfor = g2(i)
          g2(i) = 2.*sag/d**2
          grd(i) = (fs1-fs2)/(2.*d)
          if (ldebug) write (minuit_fileout,31) i,idrv,gstep(i),d,g2(i),grd(i),sag
31        format (i4,i2,6g12.5)
          gstep(i) = sign(d,gstep(i))
          dirin(i) = d
          yy(i) = fs1
          dlast = d
          d = sqrt(2.0*aimsag/abs(g2(i)))
!         if parameter has limits, max int step size = 0.5
          stpinm = 0.5
          if (gstep(i).lt.zero)  d = min(d,stpinm)
          if (d.lt.dmin)  d = dmin
!           see if converged
          if (abs((d-dlast)/d)         .lt.tlrstp)  go to 50
          if (abs((g2(i)-g2bfor)/g2(i)).lt.tlrg2 )  go to 50
          d = min(d, 10.*dlast)
          d = max(d, 0.1*dlast)
       end do
!                       end of step size loop
       write (cbf1,'(I2,2E10.2)') iext,sag,aimsag
       call mnwarn('D','MNHESS','SECOND DERIV. SAG,AIM= '//cbf1)
!
50     continue
       ndex = i*(i+1)/2
       vhmat(ndex) = g2(i)
    end do
!                              end of diagonal second derivative loop
    call mninex(x)
!                                     refine the first derivatives
    if (istrat.gt.0) call mnhes1
    isw(2) = 3
    dcovar = 0.
!                                        . . . .  off-diagonal elements
    if (npar.eq.1)  go to 214
    do i= 1, npar
       do j= 1, i-1
          xti = x(i)
          xtj = x(j)
          x(i) = xti + dirin(i)
          x(j) = xtj + dirin(j)
          call mninex(x)
          call minuit_fcn(nparx,gin,fs1,u,4)
          nfcn = nfcn + 1
          x(i) = xti
          x(j) = xtj
          elem = (fs1+amin-yy(i)-yy(j)) / (dirin(i)*dirin(j))
          ndex = i*(i-1)/2 + j
          vhmat(ndex) = elem
       end do
    end do
214 call mninex(x)
!                  verify matrix positive-definite
    call mnpsdf
    do i= 1, npar
       do j= 1, i
          ndex = i*(i-1)/2 + j
          p(i,j) = vhmat(ndex)
          p(j,i) = p(i,j)
       end do
    end do
    call mnvert(p,maxint,maxint,npar,ifail)
    if (ifail.gt.0)  then
       call mnwarn('W','HESSE', 'MATRIX INVERSION FAILS.')
       go to 390
    endif
!                                        . . . . . . .  calculate  e d m
    edm = 0.
    do i= 1, npar
!                              off-diagonal elements
       ndex = i*(i-1)/2
       do j= 1, i-1
          ndex = ndex + 1
          ztemp = 2.0 * p(i,j)
          edm = edm + grd(i)*ztemp*grd(j)
          vhmat(ndex) = ztemp
       end do
!                              diagonal elements
       ndex = ndex + 1
       vhmat(ndex) = 2.0 * p(i,i)
       edm = edm  + p(i,i) * grd(i)**2
    end do
    if (isw(5).ge.1.and.isw(2).eq.3.and.itaur.eq.0)&
         & write(minuit_fileout,'(A)')' COVARIANCE MATRIX CALCULATED SUCCESSFULLY'
    go to 900
!                              failure to invert 2nd deriv matrix
390 isw(2) = 1
    dcovar = 1.
    cstatu = 'FAILED    '
    if (isw(5).ge.0) write (minuit_fileout,'(A)')&
         &        '  MNHESS FAILS AND WILL RETURN DIAGONAL MATRIX. '
    do i= 1, npar
       ndex = i*(i-1)/2
       do j= 1, i-1
          ndex = ndex + 1
          vhmat(ndex) = 0.0
       end do
       ndex = ndex +1
       g2i = g2(i)
       if (g2i.le.zero)  g2i = 1.0
       vhmat(ndex) = 2.0/g2i
    end do
900 return
  end subroutine mnhess

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnhes1
!c      called from mnhess and mngrad
!c      calculate first derivatives (grd) and uncertainties (dgrd)
!c         and appropriate step sizes gstep
    !
    !
    logical ldebug
    character cbf1*22
    !
    integer(FINT) :: i
    integer(FINT) :: ncyc, idrv, nparx, icyc
    real(FREAL) :: dfmin, xtf, dmin, epspri, optstp, d, chgold
    real(FREAL) :: fs1, fs2, sag, grdold, grdnew, dgmin, change
    !
    ldebug = (idbg(5).ge.1)
    if (istrat.le.0) ncyc = 1
    if (istrat.eq.1) ncyc = 2
    if (istrat.gt.1) ncyc = 6
    idrv = 1
    nparx = npar
    dfmin = 4.*epsma2*(abs(amin)+up)
!                                     main loop over parameters
    do i= 1, npar
       xtf = x(i)
       dmin = 4.*epsma2*abs(xtf)
       epspri = epsma2 + abs(grd(i)*epsma2)
       optstp = sqrt(dfmin/(abs(g2(i))+epspri))
       d = 0.2 * abs(gstep(i))
       if (d.gt.optstp)  d = optstp
       if (d.lt.dmin)  d = dmin
       chgold = 10000.
!                                       iterate reducing step size
       do icyc= 1, ncyc
          x(i) = xtf + d
          call mninex(x)
          call minuit_fcn(nparx,gin,fs1,u,4)
          nfcn = nfcn + 1
          x(i) = xtf - d
          call mninex(x)
          call minuit_fcn(nparx,gin,fs2,u,4)
          nfcn = nfcn + 1
          x(i) = xtf
!                                       check if step sizes appropriate
          sag = 0.5*(fs1+fs2-2.0*amin)
          grdold = grd(i)
          grdnew = (fs1-fs2)/(2.0*d)
          dgmin = epsmac*(abs(fs1)+abs(fs2))/d
          if (ldebug) write (minuit_fileout,11) i,idrv,gstep(i),d,g2(i),grdnew,sag
11        format (i4,i2,6g12.5)
          if (grdnew.eq.zero)  go to 60
          change = abs((grdold-grdnew)/grdnew)
          if (change.gt.chgold.and.icyc.gt.1)  go to 60
          chgold = change
          grd(i) = grdnew
          gstep(i) = sign(d,gstep(i))
!                  decrease step until first derivative changes by <5%
          if (change.lt.0.05) go to 60
          if (abs(grdold-grdnew).lt.dgmin)  go to 60
          if (d.lt.dmin)  then
             call mnwarn('D','MNHES1','STEP SIZE TOO SMALL FOR 1ST DRV.')
             go to 60
          endif
          d = 0.2*d
       end do
!                                       loop satisfied = too many iter
       write (cbf1,'(2G11.3)') grdold,grdnew
       call mnwarn('D','MNHES1','TOO MANY ITERATIONS ON D1.'//cbf1)
60     continue
       dgrd(i) = max(dgmin,abs(grdold-grdnew))
    end do
!                                        end of first deriv. loop
    call mninex(x)
    return
  end subroutine mnhes1

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnwerr
!c          calculates the werr, external parameter errors,
!c      and the global correlation coefficients, to be called
!c      whenever a new covariance matrix is available.
!c
!                         calculate external error if v exists
    !
    integer(FINT) :: l, ndex, k1, j, k, ierr, iin, ndiag, i
    real(FREAL) :: dx, al, ba, du1, du2, denom
    !
    if (isw(2).ge.1) then
       do l= 1, npar
          ndex = l*(l+1)/2
          dx = sqrt(abs(vhmat(ndex)*up))
          i = nexofi(l)
          if (nvarl(i).gt.1)  then
             al = alim(i)
             ba = blim(i) - al
             du1 = al + 0.5 *(sin(x(l)+dx) +1.0) * ba - u(i)
             du2 = al + 0.5 *(sin(x(l)-dx) +1.0) * ba - u(i)
             if (dx.gt.1.0)  du1 = ba
             dx = 0.5 * (abs(du1) + abs(du2))
          endif
          werr(l) = dx
       end do
    endif
!                          global correlation coefficients
    if (isw(2).ge.1) then
       do i= 1, npar
          globcc(i) = 0.
          k1 = i*(i-1)/2
          do j= 1, i
             k = k1 + j
             p(i,j) = vhmat(k)
             p(j,i) = p(i,j)
          end do
       end do
       call mnvert(p,maxint,maxint,npar,ierr)
       if (ierr.eq.0)   then
          do iin= 1, npar
             ndiag = iin*(iin+1)/2
             denom = p(iin,iin)*vhmat(ndiag)
             if (denom.le.one.and.denom.ge.zero)  then
                globcc(iin) = 0.
             else
                globcc(iin) = sqrt(1.0-1.0/denom)
             endif
          end do
       endif
    endif
    return
  end subroutine mnwerr

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnvert(a,l,m,n,ifail)
!c        inverts a symmetric matrix.   matrix is first scaled to
!c        have all ones on the diagonal (equivalent to change of units)
!c        but no pivoting is done since matrix is positive-definite.
!c
    !
    integer(FINT) :: l, m, n, ifail
    !
    integer(FINT) :: i
    real(FREAL), dimension(l,m) :: a !(l,m)
    real(FREAL), dimension(mni) :: pp, q,  s !(mni)
    !
    real(FREAL) :: si
    !
    integer(FINT) :: j, k, kp1, km1
    !
    ifail=0
    if (n.lt.1)  go to 100
    if (n.gt.maxint)  go to 100
!                   scale matrix by sqrt of diag elements
    do i=1,n
       si = a(i,i)
       !if (si) 100,100,8
       if (si.gt.0.0_FREAL) then
          s(i) = 1.0/sqrt(si)
       else
          go to 100
       end if
    end do
    !
    do i= 1, n
       do j= 1, n
          a(i,j) = a(i,j) *s(i)*s(j)
       end do
    end do
!                                        . . . start main loop . . . .
    do i=1,n
       k = i
!                   preparation for elimination step1
       if (a(k,k).eq.zero)  go to 100
       q(k)=1./a(k,k)
       pp(k) = 1.0
       a(k,k)=0.0
       kp1=k+1
       km1=k-1
       if(km1)100,50,40
40     do j=1,km1
          pp(j)=a(j,k)
          q(j)=a(j,k)*q(k)
          a(j,k)=0.
       end do
50     if(k-n)51,60,100
51     do j=kp1,n
          pp(j)=a(k,j)
          q(j)=-a(k,j)*q(k)
          a(k,j)=0.0
       end do
!                   elimination proper
60     do j=1,n
          do k=j,n
             a(j,k)=a(j,k)+pp(j)*q(k)
          end do
       end do
    end do

!                   elements of left diagonal and unscaling
    do j= 1, n
       do k= 1, j
          a(k,j) = a(k,j) *s(k)*s(j)
          a(j,k) = a(k,j)
       end do
    end do
    return
!                   failure return
100 ifail=1
    return
  end subroutine mnvert

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnpsdf
!c        calculates the eigenvalues of v to see if positive-def.
!c        if not, adds constant along diagonal to make positive.
    !
    character chbuff*12
    real(FREAL), dimension(mni) :: s !(mni)
    !
    integer(FINT) :: ndex, ndexd, j, ifault, ip, i
    real(FREAL) :: epsmin, epspdf, dgmin, dg, pmin, pmax, padd
    !
    epsmin = 1.0e-6
    epspdf = max(epsmin, epsma2)
    dgmin = vhmat(1)
!                        check if negative or zero on diagonal
    do i= 1, npar
       ndex = i*(i+1)/2
       if (vhmat(ndex).le.zero) then
          write (chbuff(1:3),'(I3)') i
          call mnwarn('W',cfrom,&
               &'NEGATIVE DIAGONAL ELEMENT'//chbuff(1:3)//' IN ERROR MATRIX')
       endif
       if (vhmat(ndex).lt.dgmin)  dgmin = vhmat(ndex)
    end do
    if (dgmin.le.zero) then
       dg = (one+epspdf) - dgmin
       write (chbuff,'(E12.2)') dg
       call mnwarn('W',cfrom,&
            &     chbuff//' ADDED TO DIAGONAL OF ERROR MATRIX')
    else
       dg = zero
    endif
!                    store vhmat in p, make sure diagonal pos.
    do i= 1, npar
       ndex = i*(i-1)/2
       ndexd = ndex + i
       vhmat(ndexd) = vhmat(ndexd) + dg
       if (vhmat(ndexd).le.zero)   vhmat(ndexd) = 1.0
       s(i) = 1.0/sqrt(vhmat(ndexd))
       do j= 1, i
          ndex =  ndex + 1
          p(i,j) = vhmat(ndex) * s(i)*s(j)
       end do
    end do
!      call eigen (p,p,maxint,npar,pstar,-npar)
    call mneig(p,maxint,npar,maxint,pstar,epspdf,ifault)
    pmin = pstar(1)
    pmax = pstar(1)
    do ip= 2, npar
       if (pstar(ip).lt.pmin)  pmin = pstar(ip)
       if (pstar(ip).gt.pmax)  pmax = pstar(ip)
    end do
    pmax = max(abs(pmax), one)
    if ((pmin.le.zero.and.lwarn).or.isw(5).ge.2) then
       write (minuit_fileout,550)
       write (minuit_fileout,551) (pstar(ip),ip=1,npar)
    endif
    if (pmin.gt.epspdf*pmax)  go to 217
    if (isw(2).eq.3)  isw(2)=2
    padd = 1.0e-3*pmax - pmin
    do ip= 1, npar
       ndex = ip*(ip+1)/2
       vhmat(ndex) = vhmat(ndex) *(1.0 + padd)
    end do
    cstatu= 'NOT POSDEF'
    write (chbuff,'(G12.5)') padd
    call mnwarn('W',cfrom,&
         &   'MATRIX FORCED POS-DEF BY ADDING '//chbuff//' TO DIAGONAL.')
217 continue
    !
550 format (' EIGENVALUES OF SECOND-DERIVATIVE MATRIX:' )
551 format (7X,6E12.4)
    return
  end subroutine mnpsdf

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mneig(a,ndima,n,mits,work,precis,ifault)
    !
    real(FREAL), parameter :: zero=0.0, one=1.0, two=2.0, tol=1.0e-35
    integer(FINT) :: ndima, n, mits, ifault
    real(FREAL) :: precis
    real(FREAL), dimension(:,:) :: a !(ndima,*)
    real(FREAL), dimension(:) :: work !(*)
    !
    integer(FINT) :: i1, l, k, j, j1, n1, i0, m1, m, i
    real(FREAL) :: f, gl, hh, b, pt, r, pr, c, s, h
    !
    if(ndima.ne.size(a,1)) then
       stop "internal inconsistency. csm_minuit:3073"
    endif
!          precis is the machine precision epsmac
    ifault = 1
!
    i = n
    do i1 = 2,n
       l = i-2
       f = a(i,i-1)
       gl = zero
!
       if(l.lt.1) go to 25
!
       do k = 1,l
          gl = gl+a(i,k)**2
       end do
25     h = gl + f**2
!
       if(gl.gt.tol) go to 30
!
       work(i) = zero
       work(n+i) = f
       go to 65
30     l = l+1
!
       gl = sqrt(h)
!
       if(f.ge.zero) gl = -gl
!
       work(n+i) = gl
       h = h-f*gl
       a(i,i-1) = f-gl
       f = zero
       do j = 1,l
          a(j,i) = a(i,j)/h
          gl = zero
          do k = 1,j
             gl = gl+a(j,k)*a(i,k)
          end do
!
          if(j.ge.l) go to 47
!
          j1 = j+1
          do k = j1,l
             gl = gl+a(k,j)*a(i,k)
          end do
47        work(n+j) = gl/h
          f = f+gl*a(j,i)
       end do
       hh = f/(h+h)
       do j = 1,l
          f = a(i,j)
          gl = work(n+j)-hh*f
          work(n+j) = gl
          do k = 1,j
             a(j,k) = a(j,k)-f*work(n+k)-gl*a(i,k)
          end do
       end do
       work(i) = h
65     i = i-1
    end do
    work(1) = zero
    work(n+1) = zero
    do i = 1,n
       l = i-1
!
       if(work(i).eq.zero.or.l.eq.0) go to 100
!
       do j = 1,l
          gl = zero
          do k = 1,l
             gl = gl+a(i,k)*a(k,j)
          end do
          do k = 1,l
             a(k,j) = a(k,j)-gl*a(k,i)
          end do
       end do
100    work(i) = a(i,i)
       a(i,i) = one
!
       if(l.eq.0) cycle
!
       do j = 1,l
          a(i,j) = zero
          a(j,i) = zero
       end do
    end do
!
!
    n1 = n-1
    do i = 2,n
       i0 = n+i-1
       work(i0) = work(i0+1)
    end do
    work(n+n) = zero
    b = zero
    f = zero
    do l = 1,n
       j = 0
       h = precis*(abs(work(l))+abs(work(n+l)))
!
       if(b.lt.h) b = h
!
       do m1 = l,n
          m = m1
!
          if(abs(work(n+m)).le.b) exit
!
       end do
!
       if(m.eq.l) go to 205
!
160    if(j.eq.mits) return
!
       j = j+1
       pt = (work(l+1)-work(l))/(two*work(n+l))
       r = sqrt(pt*pt+one)
       pr = pt+r
!
       if(pt.lt.zero) pr=pt-r
!
       h = work(l)-work(n+l)/pr
       do i=l,n
          work(i) = work(i)-h
       end do
       f = f+h
       pt = work(m)
       c = one
       s = zero
       m1 = m-1
       i = m
       do i1 = l,m1
          j = i
          i = i-1
          gl = c*work(n+i)
          h = c*pt
!
          if(abs(pt).ge.abs(work(n+i))) go to 180
!
          c = pt/work(n+i)
          r = sqrt(c*c+one)
          work(n+j) = s*work(n+i)*r
          s = one/r
          c = c/r
          go to 190
180       c = work(n+i)/pt
          r = sqrt(c*c+one)
          work(n+j) = s*pt*r
          s = c/r
          c = one/r
190       pt = c*work(i)-s*gl
          work(j) = h+s*(c*gl+s*work(i))
          do k = 1,n
             h = a(k,j)
             a(k,j) = s*a(k,i)+c*h
             a(k,i) = c*a(k,i)-s*h
          end do
       end do
       work(n+l) = s*pt
       work(l) = c*pt
!
       if(abs(work(n+l)).gt.b) go to 160
!
205    work(l) = work(l)+f
    end do
    do i=1,n1
       k = i
       pt = work(i)
       i1 = i+1
       do j = i1,n
!
          if(work(j).ge.pt) cycle
!
          k = j
          pt = work(j)
       end do
!
       if(k.eq.i) cycle
!
       work(k) = work(i)
       work(i) = pt
       do j=1,n
          pt = a(j,i)
          a(j,i) = a(j,k)
          a(j,k) = pt
       end do
    end do
    ifault = 0
!
    return
  end subroutine mneig

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnderi
!c        calculates the first derivatives of fcn (grd),
!c        either by finite differences or by transforming the user-
!c        supplied derivatives to internal coordinates,
!c        according to whether isw(3) is zero or one.
!c
    !
    !
    logical ldebug
    character cbf1*22
    !
    integer(FINT) :: nparx, ncyc, icyc, liint, iext, i
    real(FREAL) :: fs1, df, dfmin, vrysml, tlrstp, tlrgrd, epspri
    real(FREAL) :: xtf, stepb4, optstp, step, stpmax, stpmin, fs2
    real(FREAL) :: grbfor, d1d2, dd
    !
    nparx = npar
    ldebug = (idbg(2).ge.1)
    if (amin.eq.undefi)  call mnamin
    if (isw(3).eq.1)  go to 100
    if (ldebug) then
!                       make sure starting at the right place
       call mninex(x)
       nparx = npar
       call minuit_fcn(nparx,gin,fs1,u,4)
       nfcn = nfcn + 1
       if (fs1.ne.amin) then
          df = amin - fs1
          write (cbf1(1:12),'(G12.3)') df
          call mnwarn('D','MNDERI',&
               &         'FUNCTION VALUE DIFFERS FROM AMIN BY '//cbf1(1:12) )
          amin = fs1
       endif
       write&
            &   (minuit_fileout,'(/''  FIRST DERIVATIVE DEBUG PRINTOUT.  MNDERI''/&
            &   '' PAR    DERIV     STEP      MINSTEP   OPTSTEP '',&
            &   '' D1-D2    2ND DRV'')')
    endif
    dfmin = 8. * epsma2*(abs(amin)+up)
    vrysml = 8.* epsmac**2
    if (istrat.le.0) then
       ncyc = 2
       tlrstp = 0.5
       tlrgrd = 0.1
    else if (istrat.eq.1) then
       ncyc = 3
       tlrstp = 0.3
       tlrgrd = 0.05
    else
       ncyc = 5
       tlrstp = 0.1
       tlrgrd = 0.02
    endif
!                                loop over variable parameters
    do i=1,npar
       epspri = epsma2 + abs(grd(i)*epsma2)
!         two-point derivatives always assumed necessary
!         maximum number of cycles over step size depends on strategy
       xtf = x(i)
       stepb4 = 0.
!                               loop as little as possible here!
       do icyc= 1, ncyc
!                 ........ theoretically best step
          optstp = sqrt(dfmin/(abs(g2(i))+epspri))
!                     step cannot decrease by more than a factor of ten
          step = max(optstp, abs(0.1*gstep(i)))
!                 but if parameter has limits, max step size = 0.5
          if (gstep(i).lt.zero.and.step.gt.0.5)  step=0.5
!                 and not more than ten times the previous step
          stpmax = 10.*abs(gstep(i))
          if (step.gt.stpmax)  step = stpmax
!                 minimum step size allowed by machine precision
          stpmin = max(vrysml, 8.*abs(epsma2*x(i)))
          if (step.lt.stpmin)  step = stpmin
!                 end of iterations if step change less than factor 2
          if (abs((step-stepb4)/step).lt.tlrstp)  go to 50
!         take step positive
          gstep(i) = sign(step, gstep(i))
          stepb4 = step
          x(i) = xtf + step
          call mninex(x)
          call minuit_fcn(nparx,gin,fs1,u,4)
          nfcn=nfcn+1
!         take step negative
          x(i) = xtf - step
          call mninex(x)
          call minuit_fcn(nparx,gin,fs2,u,4)
          nfcn=nfcn+1
          grbfor = grd(i)
          grd(i) = (fs1-fs2)/(2.0*step)
          g2(i) = (fs1+fs2-2.0*amin)/(step**2)
          x(i) = xtf
          if (ldebug) then
             d1d2 = (fs1+fs2-2.0*amin)/step
             write (minuit_fileout,41) i,grd(i),step,stpmin,optstp,d1d2,g2(i)
41           format (i4,2g11.3,5g10.2)
          endif
!         see if another iteration is necessary
          if (abs(grbfor-grd(i))/(abs(grd(i))+dfmin/step).lt.tlrgrd)&
               &        go to 50
       end do
!                           end of icyc loop. too many iterations
       if (ncyc.eq.1)  go to 50
       write (cbf1,'(2E11.3)')  grd(i),grbfor
       call mnwarn('D','MNDERI',&
            &         'FIRST DERIVATIVE NOT CONVERGED. '//cbf1)
50     continue
!
    end do
    call mninex(x)
    return
!                                        .  derivatives calc by fcn
100 do liint= 1, npar
       iext = nexofi(liint)
       if (nvarl(iext).gt.1)  go to 120
       grd(liint) = gin(iext)
       cycle
120    dd = (blim(iext)-alim(iext))*0.5 *cos(x(liint))
       grd(liint) = gin(iext)*dd
    end do
200 return
  end subroutine mnderi

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnline(start,fstart,step,slope,toler)
!c        perform a line search from position start
!c        along direction step, where the length of vector step
!c                   gives the expected position of minimum.
!c        fstart is value of function at start
!c        slope (if non-zero) is df/dx along step at start
!c        toler is initial tolerance of minimum in direction step
    !
    integer(FINT), parameter :: maxpt=12
    real(FREAL), parameter :: slambg=5., alpha=2.
    !
    integer(FINT) :: i
    real(FREAL), dimension(:) :: start, step !(*)
    real(FREAL), dimension(maxpt) :: xpq, ypq !(maxpt)
    character*1 chpq(maxpt)
    real(FREAL), dimension(3) :: xvals, fvals, coeff !(3)
    character*26 charal
    character*60 cmess
    !
! slambg and alpha control the maximum individual steps allowed.
! the first step is always =1. the max length of second step is slambg.
! the max size of subsequent steps is the maximum previous successful
!   step multiplied by alpha + the size of most recent successful step,
!   but cannot be smaller than slambg.
    !
    logical ldebug
    !
    integer(FINT) :: nparx, kk, nxypt, ipt, nvmax
    real(FREAL) :: fstart, slope, toler, overal, undral, f1
    real(FREAL) :: fvmin, xvmin, slamin, ratio, slam, toler8
    real(FREAL) :: slamax, flast, denom, f2, sdev, slopem
    real(FREAL) :: toler9, f3, fvmax
    !
    data charal / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
    !
    ldebug = (idbg(1).ge.1)
!                  starting values for overall limits on total step slam
    overal = 1000.
    undral = -100.
!                              debug check if start is ok
    if (ldebug)  then
       call mninex(start)
       call minuit_fcn(nparx,gin,f1,u,4)
       nfcn=nfcn+1
       if (f1.ne.fstart) then
          write (minuit_fileout,'(A/2E14.5/2X,10F10.5)')&
               & ' MNLINE START POINT NOT CONSISTENT, F VALUES, PARAMETERS=',&
               &  (x(kk),kk=1,npar)
       endif
    endif
!                                      . set up linear search along step
    fvmin = fstart
    xvmin = zero
    nxypt = 1
    chpq(1) = charal(1:1)
    xpq(1) = 0.
    ypq(1) = fstart
!               slamin = smallest possible value of abs(slam)
    slamin = zero
    do i= 1, npar
       if (step(i).eq.zero) cycle
       ratio = abs(start(i)/step(i))
       if (slamin.eq.zero)     slamin = ratio
       if (ratio.lt.slamin)  slamin = ratio
       x(i) = start(i) + step(i)
    end do
    if (slamin.eq.zero)  slamin = epsmac
    slamin = slamin*epsma2
    nparx = npar
!
    call mninex(x)
    call minuit_fcn(nparx,gin,f1,u,4)
    nfcn=nfcn+1
    nxypt = nxypt + 1
    chpq(nxypt) = charal(nxypt:nxypt)
    xpq(nxypt) = 1.
    ypq(nxypt) = f1
    if (f1.lt.fstart) then
       fvmin = f1
       xvmin = 1.0
    endif
!                         . quadr interp using slope gdel and two points
    slam = 1.
    toler8 = toler
    slamax = slambg
    flast = f1
!                         can iterate on two-points (cut) if no imprvmnt
25  continue
    denom = 2.0*(flast-fstart-slope*slam)/slam**2
!     if (denom.eq.zero)  denom = -0.1*slope
    slam  = 1.
    if (denom.ne.zero)  slam = -slope/denom
    if (slam.lt.zero)  slam = slamax
    if (slam.gt.slamax)  slam = slamax
    if (slam.lt.toler8)  slam = toler8
    if (slam.lt.slamin)  go to 80
    if (abs(slam-1.0).lt.toler8.and.f1.lt.fstart)  go to 70
    if (abs(slam-1.0).lt.toler8) slam = 1.0+toler8
    if (nxypt.ge.maxpt) go to 65
    do i= 1, npar
       x(i) = start(i) + slam*step(i)
    end do
    call mninex(x)
    call minuit_fcn(npar,gin,f2,u,4)
    nfcn = nfcn + 1
    nxypt = nxypt + 1
    chpq(nxypt) = charal(nxypt:nxypt)
    xpq(nxypt) = slam
    ypq(nxypt) = f2
    if (f2.lt.fvmin)  then
       fvmin = f2
       xvmin = slam
    endif
    if (fstart.eq.fvmin) then
       flast = f2
       toler8 = toler*slam
       overal = slam-toler8
       slamax = overal
       go to 25
    endif
!                                        . quadr interp using 3 points
    xvals(1) = xpq(1)
    fvals(1) = ypq(1)
    xvals(2) = xpq(nxypt-1)
    fvals(2) = ypq(nxypt-1)
    xvals(3) = xpq(nxypt)
    fvals(3) = ypq(nxypt)
!                             begin iteration, calculate desired step
50  continue
    slamax = max(slamax,alpha*abs(xvmin))
    call mnpfit(xvals,fvals,3,coeff,sdev)
    if (coeff(3).le.zero)  then
       slopem = 2.0*coeff(3)*xvmin + coeff(2)
       if (slopem.le.zero) then
          slam = xvmin + slamax
       else
          slam = xvmin - slamax
       endif
    else
       slam = -coeff(2)/(2.0*coeff(3))
       if (slam.gt.xvmin+slamax)  slam = xvmin+slamax
       if (slam.lt.xvmin-slamax)  slam = xvmin-slamax
    endif
    if (slam.gt.zero) then
       if (slam.gt.overal) slam = overal
    else
       if (slam.lt.undral) slam = undral
    endif
!               come here if step was cut below
52  continue
    toler9 = max(toler8,abs(toler8*slam))
    do ipt= 1, 3
       if (abs(slam-xvals(ipt)).lt.toler9)  go to 70
    end do
!                take the step
    if (nxypt.ge.maxpt) go to 65
    do i= 1, npar
       x(i) = start(i)+slam*step(i)
    end do
    call mninex(x)
    call minuit_fcn(nparx,gin,f3,u,4)
    nfcn = nfcn + 1
    nxypt = nxypt + 1
    chpq(nxypt) = charal(nxypt:nxypt)
    xpq(nxypt) = slam
    ypq(nxypt) = f3
!             find worst previous point out of three
    fvmax = fvals(1)
    nvmax = 1
    if (fvals(2).gt.fvmax) then
       fvmax = fvals(2)
       nvmax = 2
    endif
    if (fvals(3).gt.fvmax) then
       fvmax = fvals(3)
       nvmax = 3
    endif
!              if latest point worse than all three previous, cut step
    if (f3.ge.fvmax)  then
       if (nxypt.ge.maxpt) go to 65
       if (slam.gt.xvmin) overal = min(overal,slam-toler8)
       if (slam.lt.xvmin) undral = max(undral,slam+toler8)
       slam = 0.5*(slam+xvmin)
       go to 52
    endif
!              prepare another iteration, replace worst previous point
    xvals(nvmax) = slam
    fvals(nvmax) = f3
    if (f3.lt.fvmin)  then
       fvmin = f3
       xvmin = slam
    else
       if (slam.gt.xvmin) overal = min(overal,slam-toler8)
       if (slam.lt.xvmin) undral = max(undral,slam+toler8)
    endif
    if (nxypt.lt.maxpt)  go to 50
!                                            . . end of iteration . . .
!            stop because too many iterations
65  cmess = ' LINE SEARCH HAS EXHAUSTED THE LIMIT OF FUNCTION CALLS '
    if (ldebug) then
       write (minuit_fileout,'(A/(2X,6G12.4))') ' MNLINE DEBUG: STEPS=',&
            &    (step(kk),kk=1,npar)
    endif
    go to 100
!            stop because within tolerance
70  continue
    cmess = ' LINE SEARCH HAS ATTAINED TOLERANCE '
    go to 100
80  continue
    cmess = ' STEP SIZE AT ARITHMETICALLY ALLOWED MINIMUM'
100 continue
    amin = fvmin
    do i= 1, npar
       dirin(i) = step(i)*xvmin
       x(i) = start(i) + dirin(i)
    end do
    call mninex(x)
    if (xvmin.lt.0.)      call mnwarn('D','MNLINE',&
         &                   ' LINE MINIMUM IN BACKWARDS DIRECTION')
    if (fvmin.eq.fstart)  call mnwarn('D','MNLINE',&
         &                     ' LINE SEARCH FINDS NO IMPROVEMENT ')
    if (ldebug)  then
       write (minuit_fileout,'('' AFTER'',I3,'' POINTS,'',A)') nxypt,cmess
       call mnplot(xpq,ypq,chpq,nxypt,minuit_fileout,npagwd,npagln)
    endif
    return
  end subroutine mnline

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnplot(xpt,ypt,chpt,nxypt,nunit,npagwd,npagln)
!c        plots points in array xypt onto one page with labelled axes
!c        nxypt is the number of points to be plotted
!c        xpt(i) = x-coord. of ith point
!c        ypt(i) = y-coord. of ith point
!c        chpt(i) = character to be plotted at this position
!c        the input point arrays xpt, ypt, chpt are destroyed.
!c
    !
    integer(FINT), parameter :: maxwid=100
    !
    integer(FINT) :: nxypt, nunit, npagwd, npagln
    !
    real(FREAL), dimension(:) :: xpt, ypt !(*)
    character*1 chpt(*) ,  chsav,  chbest, cdot, cslash, cblank
    character cline*100, chmess*30
    real(FREAL), dimension(12) :: xvalus !(12)
    logical overpr
    !
    integer(FINT) :: i
    integer(FINT) :: maxnx, maxny, km1, iquit, ni, j, nx, ny, nxbest
    integer(FINT) :: nybest, isp1, linodd, ibk, k, ks, ix, iten
    real(FREAL) :: xbest, ybest, savx, savy, xmax, xmin, dxx, yprt
    real(FREAL) :: bwidx, ymax, ymin, dyy, bwidy, rlany, ax, ay, bx, by
    !
    data cdot,cslash,cblank/ '.' , '/' , ' '/
    !
    maxnx = min(npagwd-20,maxwid)
    if (maxnx.lt.10)  maxnx = 10
    maxny = npagln
    if (maxny.lt.10)  maxny = 10
    if (nxypt.le.1)  return
    xbest = xpt(1)
    ybest = ypt(1)
    chbest = chpt(1)
!         order the points by decreasing y
    km1 = nxypt - 1
    do i= 1, km1
       iquit = 0
       ni = nxypt - i
       do j= 1, ni
          if (ypt(j).gt.ypt(j+1)) cycle
          savx = xpt(j)
          xpt(j) = xpt(j+1)
          xpt(j+1) = savx
          savy = ypt(j)
          ypt(j) = ypt(j+1)
          ypt(j+1) = savy
          chsav = chpt(j)
          chpt(j) = chpt(j+1)
          chpt(j+1) = chsav
          iquit = 1
       end do
       if (iquit.eq.0) go to 160
    end do
160 continue
!         find extreme values
    xmax = xpt(1)
    xmin = xmax
    do i= 1, nxypt
       if (xpt(i).gt.xmax)  xmax = xpt(i)
       if (xpt(i).lt.xmin)  xmin = xpt(i)
    end do
    dxx = 0.001*(xmax-xmin)
    xmax = xmax + dxx
    xmin = xmin - dxx
    call mnbins(xmin,xmax,maxnx,xmin,xmax,nx,bwidx)
    ymax = ypt(1)
    ymin = ypt(nxypt)
    if (ymax.eq.ymin)  ymax=ymin+1.0
    dyy = 0.001*(ymax-ymin)
    ymax = ymax + dyy
    ymin = ymin - dyy
    call mnbins(ymin,ymax,maxny,ymin,ymax,ny,bwidy)
    rlany = ny
!         if first point is blank, it is an 'origin'
    if (chbest.eq.cblank)  go to 50
    xbest = 0.5 * (xmax+xmin)
    ybest = 0.5 * (ymax+ymin)
50  continue
!         find scale constants
    ax = 1.0/bwidx
    ay = 1.0/bwidy
    bx = -ax*xmin + 2.0
    by = -ay*ymin - 2.0
!         convert points to grid positions
    do i= 1, nxypt
       xpt(i) = ax*xpt(i) + bx
       ypt(i) = rlany-ay*ypt(i) - by
    end do
    nxbest = ax*xbest + bx
    nybest = rlany  - ay*ybest - by
!         print the points
    ny = ny + 2
    nx = nx + 2
    isp1 = 1
    linodd = 1
    overpr=.false.
    do i= 1, ny
       do ibk= 1, nx
          cline (ibk:ibk) = cblank
       end do
       cline(1:1) = cdot
       cline(nx:nx) = cdot
       cline(nxbest:nxbest) = cdot
       if (i.ne.1.and.i.ne.nybest.and.i.ne.ny)  go to 320
       do j= 1, nx
          cline(j:j) = cdot
       end do
320    continue
       yprt = ymax - float(i-1)*bwidy
       if (isp1.gt.nxypt)  go to 350
!         find the points to be plotted on this line
       do k= isp1,nxypt
          ks = ypt(k)
          if (ks.gt.i)  go to 345
          ix = xpt(k)
          if (cline(ix:ix).eq.cdot)  go to 340
          if (cline(ix:ix).eq.cblank)  go to 340
          if (cline(ix:ix).eq.chpt(k)) cycle
          overpr = .true.
!         overpr is true if one or more positions contains more than
!            one point
          cline(ix:ix) = '&'
          cycle
340       cline(ix:ix) = chpt(k)
       end do
       isp1 = nxypt + 1
       go to 350
345    isp1 = k
350    continue
       if (linodd.eq.1.or.i.eq.ny)  go to 380
       linodd = 1
       write (nunit, '(18X,A)')       cline(:nx)
       cycle
380    write (nunit,'(1X,G14.7,A,A)') yprt, ' ..', cline(:nx)
       linodd = 0
    end do
!         print labels on x-axis every ten columns
    do ibk= 1, nx
       cline(ibk:ibk) = cblank
       if (mod(ibk,10).eq.1)  cline(ibk:ibk) = cslash
    end do
    write (nunit, '(18X,A)')       cline(:nx)
!
    do ibk= 1, 12
       xvalus(ibk) = xmin + float(ibk-1)*10.*bwidx
    end do
    iten = (nx+9) / 10
    write (nunit,'(12X,12G10.4)')  (xvalus(ibk), ibk=1,iten)
    chmess = ' '
    if (overpr) chmess='   OVERPRINT CHARACTER IS&'
    write (nunit,'(25X,A,G13.7,A)') 'ONE COLUMN=',bwidx, chmess
500 return
  end subroutine mnplot

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnbins(a1,a2,naa,bl,bh,nb,bwid)
!         subroutine to determine reasonable histogram intervals
!         given absolute upper and lower bounds  a1 and a2
!         and desired maximum number of bins naa
!         program makes reasonable binning from bl to bh of width bwid
!         f. james,   august, 1974 , stolen for minuit, 1988
    !
    real(FREAL), parameter :: zero=0.0, one=1.0
    !
    real(FREAL) :: a1, a2, bl, bh, bwid
    integer(FINT) :: naa, nb
    !
    integer(FINT) :: na, llog, lwid, kwid
    real(FREAL) :: al, ah, awid, sigfig, sigrnd, alb
    !
    al = min(a1,a2)
    ah = max(a1,a2)
    if (al.eq.ah)  ah = al + 1.
!         if naa.eq.-1 , program uses bwid input from calling routine
    if (naa.eq.-1)  go to 150
10  na = naa - 1
    if (na.lt.1)  na = 1
!          get nominal bin width in expon form
20  awid = (ah-al)/float(na)
    llog = int(dlog10(dble(awid)))
    if (awid.le.one)  llog=llog-1
    sigfig = awid * (10.00 **(-llog))
!         round mantissa up to 2, 2.5, 5, or 10
    if(sigfig.gt.2.0)  go to 40
    sigrnd = 2.0
    go to 100
40  if (sigfig.gt.2.5)  go to 50
    sigrnd = 2.5
    go to 100
50  if(sigfig.gt.5.0)  go to 60
    sigrnd =5.0
    go to 100
60  sigrnd = 1.0
    llog = llog + 1
100 continue
    bwid = sigrnd*10.0**llog
    go to 200
!         get new bounds from new width bwid
150 if (bwid.le.zero)  go to 10
200 continue
    alb = al/bwid
    lwid=alb
    if (alb.lt.zero)  lwid=lwid-1
    bl = bwid*float(lwid)
    alb = ah/bwid + 1.0
    kwid = alb
    if (alb.lt.zero)  kwid=kwid-1
    bh = bwid*float(kwid)
    nb = kwid-lwid
    if (naa.gt.5)  go to 240
    if (naa.eq.-1)  return
!          request for one bin is difficult case
    if (naa.gt.1.or.nb.eq.1)  return
    bwid =  bwid*2.0
    nb  = 1
    return
240 if (2*nb.ne.naa)  return
    na = na + 1
    go to 20
  end subroutine mnbins

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnpfit(parx2p,pary2p,npar2p,coef2p,sdev2p)
!
!     to fit a parabola to npar2p points
!
!   npar2p   no. of points
!   parx2p(i)   x value of point i
!   pary2p(i)   y value of point i
!
!   coef2p(1...3)  coefficients of the fitted parabola
!   y=coef2p(1) + coef2p(2)*x + coef2p(3)*x**2
!   sdev2p= variance
!   method : chi**2 = min equation solved explicitly
    !
    integer(FINT) :: npar2p
    real(FREAL) :: sdev2p
    real(FREAL), dimension(npar2p) :: parx2p, pary2p, coef2p !(npar2p)
    real(FREAL), dimension(3) :: cz !(3)
    !
    integer(FINT) :: i
    real(FREAL) :: f, xm, x2, x3, x4, y, y2, xy, x2y, s, t, s2, a
    !
!
    do i=1,3
       cz(i)=0.
    end do
    sdev2p=0.
    if(npar2p.lt.3) go to 10
    f=npar2p
!--- center x values for reasons of machine precision
    xm=0.
    do i=1,npar2p
       xm=xm+parx2p(i)
    end do
    xm=xm/f
    x2=0.
    x3=0.
    x4=0.
    y=0.
    y2=0.
    xy=0.
    x2y=0.
    do i=1,npar2p
       s=parx2p(i)-xm
       t=pary2p(i)
       s2=s*s
       x2=x2+s2
       x3=x3+s*s2
       x4=x4+s2*s2
       y=y+t
       y2=y2+t*t
       xy=xy+s*t
       x2y=x2y+s2*t
    end do
    a=(f*x4-x2**2)*x2-f*x3**2
    if(a.eq.0.)  goto 10
    cz(3)=(x2*(f*x2y-x2*y)-f*x3*xy)/a
    cz(2)=(xy-x3*cz(3))/x2
    cz(1)=(y-x2*cz(3))/f
    if(npar2p.eq.3)  goto 6
    sdev2p=y2-(cz(1)*y+cz(2)*xy+cz(3)*x2y)
    if(sdev2p.lt.0.)  sdev2p=0.
    sdev2p=sdev2p/(f-3.)
6   cz(1)=cz(1)+xm*(xm*cz(3)-cz(2))
    cz(2)=cz(2)-2.*xm*cz(3)
10  continue
    do i=1,3
       coef2p(i)=cz(i)
    end do
    return
  end subroutine mnpfit

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnmatu(kode)
!c        prints the covariance matrix v when kode=1.
!c        always prints the global correlations, and
!c        calculates and prints the individual correlation coefficients
!c
    !
    integer(FINT) :: kode
    real(FREAL), dimension(mni) :: vline !(mni)
    !
    integer(FINT) :: i
    integer(FINT) :: isw2, isw5, ncoef, nparm, id, ix, ndi, j, m, n
    integer(FINT) :: ndex, ndj, it, iso, nsofar
    !
    isw2 = isw(2)
    if (isw2.lt.1)  then
       write (minuit_fileout,'(1X,A)')  covmes(isw2)
       go to 500
    endif
    if (npar.eq.0)  then
       write (minuit_fileout,'('' MNMATU: NPAR=0'')')
       go to 500
    endif
!                                       . . . . .external error matrix
    if (kode.eq.1)  then
       isw5 = isw(5)
       isw(5) = 2
       call minuit_mnemat(p,maxint)
       if (isw2.lt.3)  write (minuit_fileout,'(1X,A)')  covmes(isw2)
       isw(5) = isw5
    endif
!                                       . . . . . correlation coeffs. .
    if (npar.le.1)   go to 500
    call mnwerr
!     ncoef is number of coeff. that fit on one line, not to exceed 20
    ncoef = (npagwd-19)/6
    ncoef = min(ncoef,20)
    nparm = min(npar,ncoef)
    write (minuit_fileout, 150) (nexofi(id),id=1,nparm)
150 format (/36h parameter  correlation coefficients  /&
         &         18h       no.  global   ,20i6)
    do i= 1, npar
       ix = nexofi(i)
       ndi = i*(i+1)/2
       do j= 1, npar
          m = max(i,j)
          n = min(i,j)
          ndex = m*(m-1)/2 + n
          ndj = j*(j+1)/2
          vline(j) = vhmat(ndex)/sqrt(abs(vhmat(ndi)*vhmat(ndj)))
       end do
       nparm = min(npar,ncoef)
       write (minuit_fileout,171)   ix, globcc(i), (vline(it),it=1,nparm)
171    format (6x,i3,2x,f7.5,1x,20f6.3)
       if (i.le.nparm) cycle
       do iso= 1, 10
          nsofar = nparm
          nparm = min(npar,nsofar+ncoef)
          write (minuit_fileout,181)  (vline(it),it=nsofar+1,nparm)
181       format (19x,20f6.3)
          if (i.le.nparm) exit
       end do
192    continue
    end do
    if (isw2.lt.3)  write (minuit_fileout,'(1X,A)')  covmes(isw2)
500 return
  end subroutine mnmatu

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine minuit_mnemat(emat,ndim)
    !
    integer(FINT) :: ndim
    real(FREAL), dimension(ndim,ndim) :: emat !(ndim,ndim)
!c        calculates the external error matrix from the internal
!c        to be called by user, who must dimension emat at (ndim,ndim)
    !
    integer(FINT) :: i
    integer(FINT) :: npard, nperln, kga, j, kgb, iz, k, k2, kk
    real(FREAL) :: dxdi, dxdj
    !
    if (isw(2).lt.1)  return
    if (isw(5).ge.2)  write (minuit_fileout,'(/A,I4,A,I3,A,G10.3)')&
         &    ' EXTERNAL ERROR MATRIX.    NDIM=',ndim,'    NPAR=',npar,&
         &    '    ERR DEF=',up
!                    size of matrix to be printed
    npard = npar
    if (ndim.lt.npar)  then
       npard = ndim
       if (isw(5).ge.0) write (minuit_fileout,'(A,A)') ' USER-DIMENSIONED ',&
            &      ' ARRAY EMAT NOT BIG ENOUGH. REDUCED MATRIX CALCULATED.'
    endif
!                 nperln is the number of elements that fit on one line
    nperln = (npagwd-5)/10
    nperln = min(nperln,13)
    if (isw(5).ge.1.and.npard.gt.nperln)  write (minuit_fileout,'(A)')&
         &     ' ELEMENTS ABOVE DIAGONAL ARE NOT PRINTED.'
!                 i counts the rows of the matrix
    do i= 1, npard
       call mndxdi(x(i),i,dxdi)
       kga = i*(i-1)/2
       do j= 1, i
          call mndxdi(x(j),j,dxdj)
          kgb = kga + j
          emat(i,j) = dxdi * vhmat(kgb) * dxdj * up
          emat(j,i) = emat(i,j)
       end do
    end do
!                    iz is number of columns to be printed in row i
    if (isw(5).ge.2)  then
       do i= 1, npard
          iz = npard
          if (npard.ge.nperln)  iz = i
          do k= 1, iz, nperln
             k2 = k + nperln - 1
             if (k2.gt.iz)  k2=iz
             write (minuit_fileout,'(1X,13E10.3)')  (emat(i,kk),kk=k,k2)
          end do
       end do
    endif
    return
  end subroutine minuit_mnemat

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mncuve
    !
!c        makes sure that the current point is a local
!c        minimum and that the error matrix exists,
!c        or at least something good enough for minos and mncont
!c
    !
    !
    integer(FINT) :: i
    integer(FINT) :: ndex, j, iext
    real(FREAL) :: wint, dxdi
    !
    if (isw(4).lt.1) then
       write (minuit_fileout,'(/A,A)')&
            &    ' FUNCTION MUST BE MINIMIZED BEFORE CALLING ',cfrom
       apsi = epsi
       call mnmigr
    endif
    if (isw(2).lt.3)  then
       call mnhess
       if (isw(2).lt.1)  then
          call mnwarn('W',cfrom,'NO ERROR MATRIX.  WILL IMPROVISE.')
          do i=1,npar
             ndex = i*(i-1)/2
             do j=1,i-1
                ndex = ndex + 1
                vhmat(ndex) = 0.
             end do
             ndex = ndex + 1
             if (g2(i).le.zero)  then
                wint = werr(i)
                iext = nexofi(i)
                if (nvarl(iext).gt.1) then
                   call mndxdi(x(i),i,dxdi)
                   if (abs(dxdi).lt..001) then
                      wint = .01
                   else
                      wint = wint/abs(dxdi)
                   endif
                endif
                g2(i) = up/wint**2
             endif
             vhmat(ndex) = 2./g2(i)
          end do
          isw(2) = 1
          dcovar = 1.
       else
          call mnwerr
       endif
    endif
    return
  end subroutine mncuve

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnmnos
!c        performs a minos error analysis on those parameters for
!c        which it is requested on the minos command by calling
!c        mnmnot for each parameter requested.
    !
    !
    integer(FINT) :: ngood, nbad, nfcnmi, knt, ilax, ilax2, iin
    real(FREAL) :: val2pl, val2mi
    !
    if (npar.le.0)  go to 700
    ngood = 0
    nbad = 0
    nfcnmi = nfcn
!                                      . loop over parameters requested
    do knt= 1, npar
       if (int(word7(2)).eq.0) then
          ilax = nexofi(knt)
       else
          if (knt.ge.7)  go to 580
          ilax = int(word7(knt+1))
          if (ilax.eq.0)  go to 580
          if (ilax.gt.0.and.ilax.le.nu) then
             if (niofex(ilax).gt.0)  go to 565
          endif
          write (minuit_fileout,564) ilax
564       format (' PARAMETER NUMBER ',I5,' NOT VARIABLE. IGNORED.')
          cycle
       endif
565    continue
!                                         calculate one pair of m e's
       ilax2 = 0
       call mnmnot(ilax,ilax2,val2pl,val2mi)
       if (lnewmn)  go to 650
!                                          update ngood and nbad
       iin = niofex(ilax)
       if (erp(iin).gt.zero) then
          ngood=ngood+1
       else
          nbad=nbad+1
       endif
       if (ern(iin).lt.zero) then
          ngood=ngood+1
       else
          nbad=nbad+1
       endif
    end do
!                                           end of loop . . . . . . .
580 continue
!                                        . . . . printout final values .
    cfrom = 'MINOS   '
    nfcnfr = nfcnmi
    cstatu= 'UNCHANGED '
    if (ngood.eq.0.and.nbad.eq.0) go to 700
    if (ngood.gt.0.and.nbad.eq.0) cstatu='SUCCESSFUL'
    if (ngood.eq.0.and.nbad.gt.0) cstatu='FAILURE   '
    if (ngood.gt.0.and.nbad.gt.0) cstatu='PROBLEMS  '
    if (isw(5).ge.0) call mnprin(4,amin)
    if (isw(5).ge.2) call mnmatu(0)
    go to 900
!                                        . . . new minimum found . . . .
650 continue
    cfrom = 'MINOS   '
    nfcnfr = nfcnmi
    cstatu= 'NEW MINIMU'
    if (isw(5).ge.0) call mnprin(4,amin)
    write (minuit_fileout,675)
675 format(/50H NEW MINIMUM FOUND.  GO BACK TO MINIMIZATION STEP./1H ,&
         &60(1H=)/60X,1HV/60X,1HV/60X,1HV/57X,7HVVVVVVV/58X,5HVVVVV/59X,&
         &3HVVV/60X,1HV//)
    go to 900
700 write (minuit_fileout,'(A)') ' THERE ARE NO MINOS ERRORS TO CALCULATE.'
900 return
  end subroutine mnmnos

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnmnot(ilax,ilax2,val2pl,val2mi)
!c        performs a minos error analysis on one parameter.
!c        the parameter ilax is varied, and the minimum of the
!c        function with respect to the other parameters is followed
!c        until it crosses the value fmin+up.
!c
    !
    !
    character(4), parameter :: cpos='POSI',cneg='NEGA'
    integer(FINT) :: ilax, ilax2
    real(FREAL) :: val2pl, val2mi
    real(FREAL), dimension(mni) :: xdev, w, gcc !(mni)
    character*4 csig
    !
    integer(FINT) :: i
    integer(FINT) :: isw2, isw4, istrav, mpar, nfmxin, j, it
    integer(FINT) :: ndex, marc, imax, indx, ierr, isig, iercr
    real(FREAL) :: sigsav, dc, abest, ut, xunit, sig, nlimit
    real(FREAL) :: du1, delu, fac, aopt, eros, sav
    !
!                                        . . save and prepare start vals
    isw2 = isw(2)
    isw4 = isw(4)
    sigsav = edm
    istrav = istrat
    dc = dcovar
    lnewmn = .false.
    apsi  = epsi*0.5
    abest=amin
    mpar=npar
    nfmxin = nfcnmx
    do i= 1, mpar
       xt(i) = x(i)
    end do
    do j= 1, mpar*(mpar+1)/2
       vthmat(j) = vhmat(j)
    end do
    do i= 1, mpar
       gcc(i) = globcc(i)
       w(i) = werr(i)
    end do
    it = niofex(ilax)
    erp(it) = 0.
    ern(it) = 0.
    call mninex(xt)
    ut = u(ilax)
    if (nvarl(ilax).eq.1) then
       alim(ilax) = ut -100.*w(it)
       blim(ilax) = ut +100.*w(it)
    endif
    ndex = it*(it+1)/2
    xunit = sqrt(up/vthmat(ndex))
    marc = 0
    do i= 1, mpar
       if (i.eq.it) cycle
       marc = marc + 1
       imax = max(it,i)
       indx = imax*(imax-1)/2 + min(it,i)
       xdev(marc) = xunit*vthmat(indx)
    end do
!                           fix the parameter in question
    call mnfixp (it,ierr)
    if (ierr.gt.0)  then
       write (minuit_fileout,'(A,I5,A,I5)')&
            &    ' MINUIT ERROR. CANNOT FIX PARAMETER',ilax,'    INTERNAL',it
       go to 700
    endif
!                       . . . . . nota bene: from here on, npar=mpar-1
!      remember: mnfixp squeezes it out of x, xt, werr, and vhmat,
!                                                    not w, vthmat
    do isig= 1,2
       if (isig.eq.1) then
          sig = 1.0
          csig = cpos
       else
          sig = -1.0
          csig = cneg
       endif
!                                        . sig=sign of error being calcd
       if (isw(5).gt.1) write (minuit_fileout,806)  csig,ilax,cpnam(ilax)
806    format (/' DETERMINATION OF ',A4,'TIVE MINOS ERROR FOR PARAMETER',&
            &    I3, 2X ,A)
       if (isw(2).le.0) call mnwarn('D','MINOS','NO COVARIANCE MATRIX.')
       nlimit = nfcn + nfmxin
       istrat = max(istrav-1,0)
       du1 = w(it)
       u(ilax) = ut + sig*du1
       u(ilax) = min(u(ilax),blim(ilax))
       u(ilax) = max(u(ilax),alim(ilax))
       delu = u(ilax) - ut
!         stop if already at limit with negligible step size
       if (abs(delu)/(abs(ut)+abs(u(ilax))).lt.epsmac)  go to 440
       fac = delu/w(it)
       do i= 1, npar
          x(i) = xt(i) + fac*xdev(i)
       end do
       if (isw(5).gt.1) write (minuit_fileout,801)  ilax,ut,delu,u(ilax)
801    format (/' PARAMETER',I4,' SET TO',E11.3,' + ',E10.3,' = ',E12.3)
!                                        loop to hit amin+up
       ke1cr = ilax
       ke2cr = 0
       xmidcr = u(ilax)
       xdircr = delu
!
       amin = abest
       nfcnmx = nlimit - nfcn
       call mncros(aopt,iercr)
       if (abest-amin.gt.0.01*up)  go to 650
       if (iercr.eq.1)  go to 440
       if (iercr.eq.2)  go to 450
       if (iercr.eq.3)  go to 460
!                                        . error successfully calculated
       eros = xmidcr-ut + aopt*xdircr
       if (isw(5).gt.1) write (minuit_fileout,808)  csig,ilax,cpnam(ilax),eros
808    format (/9x,4hthe ,a4,  29htive minos error of parameter,i3,   2h&
            &, ,a10,      4h, is ,e12.4)
       go to 480
!                                        . . . . . . . . failure returns
440    if (isw(5).ge.1) write(minuit_fileout,807)  csig,ilax,cpnam(ilax)
807    format (5X,'THE ',A4,'TIVE MINOS ERROR OF PARAMETER',I3,', ',A,&
            &', EXCEEDS ITS LIMIT.'/)
       eros = undefi
       go to 480
450    if (isw(5).ge.1) write (minuit_fileout, 802)  csig,ilax,nfmxin
802    format (9X,'THE ',A,'TIVE MINOS ERROR',I4,' REQUIRES MORE THAN',&
            &   I5,' FUNCTION CALLS.'/)
       eros = 0.
       go to 480
460    if (isw(5).ge.1) write (minuit_fileout, 805) csig,ilax
805    format (25X,A,'TIVE MINOS ERROR NOT CALCULATED FOR PARAMETER',I4/)
       eros = 0.
!
480    if (isw(5).gt.1) write (minuit_fileout,'(5X, 74(1H*))')
       if (sig.lt.zero)  then
          ern(it) = eros
          if (ilax2.gt.0.and.ilax2.le.nu)  val2mi = u(ilax2)
       else
          erp(it) = eros
          if (ilax2.gt.0.and.ilax2.le.nu)  val2pl = u(ilax2)
       endif
    end do
!                                        . . parameter finished. reset v
!                       normal termination
    itaur = 1
    call mnfree(1)
    do j= 1, mpar*(mpar+1)/2
       vhmat(j) = vthmat(j)
    end do
    do i= 1, mpar
       werr(i) = w(i)
       globcc(i) = gcc(i)
       x(i) = xt(i)
    end do
    call mninex (x)
    edm = sigsav
    amin = abest
    isw(2) = isw2
    isw(4) = isw4
    dcovar = dc
    go to 700
!                       new minimum
650 lnewmn = .true.
    isw(2) = 0
    dcovar = 1.
    isw(4) = 0
    sav = u(ilax)
    itaur = 1
    call mnfree(1)
    u(ilax) = sav
    call mnexin(x)
    edm = bigedm
!                       in any case
700 continue
    itaur = 0
    nfcnmx = nfmxin
    istrat = istrav
    return
  end subroutine mnmnot

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mncros(aopt,iercr)
!c       find point where mneval=amin+up, along the line through
!c       xmidcr,ymidcr with direction xdircr,ydircr,   where x and y
!c       are parameters ke1cr and ke2cr.  if ke2cr=0 (from minos),
!c       only ke1cr is varied.  from mncont, both are varied.
!c       crossing point is at
!c        (u(ke1),u(ke2)) = (xmid,ymid) + aopt*(xdir,ydir)
!c
    !
    character(10), parameter :: chere='MNCROS    '
    integer(FINT), parameter :: mlsb=3, maxitr=15
    real(FREAL), parameter :: tlr=0.01
    character charal*28, chsign*4
    integer(FINT) :: i
    real(FREAL), dimension(mlsb) :: flsb, alsb !(mlsb)
    real(FREAL), dimension(3) :: coeff !(3)
    logical ldebug
    !
    integer(FINT) :: iercr, ipt, ik, kex, ierev, it, inew
    integer(FINT) :: ibest, iworst, ileft, iright, iout, itoohi

    real(FREAL) :: aopt, aminsv, aim, tlf, tla, aulim, zmid, zdir
    real(FREAL) :: zlim, anext, fnext, dfda, fdist, adist
    real(FREAL) :: bmin, bmax, ecarmn, ecarmx, noless, ecart, sdev
    real(FREAL) :: determ, rt, x1, x2, s1, s2, slope, smalla, aleft
    real(FREAL) :: aright
    integer(FINT) :: maxlk
    !
    data  charal/' .ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
    !
    ldebug = (idbg(6).ge.1)
    aminsv = amin
!        convergence when f is within tlf of aim and next prediction
!        of aopt is within tla of previous value of aopt
    aim = amin + up
    tlf = tlr*up
    tla = tlr
    xpt(1) = 0.0
    ypt(1) = aim
    chpt(1) = ' '
    ipt = 1
    if (ke2cr.eq.0) then
       xpt(2) = -1.0
       ypt(2) = amin
       chpt(2) = '.'
       ipt = 2
    endif
!                    find the largest allowed a
    aulim = 100.
    do ik= 1, 2
       if (ik.eq.1)  then
          kex = ke1cr
          zmid = xmidcr
          zdir = xdircr
       else
          if (ke2cr.eq.0) cycle
          kex = ke2cr
          zmid = ymidcr
          zdir = ydircr
       endif
       if (nvarl(kex).le.1) cycle
       if (zdir.eq.zero) cycle
       zlim = alim(kex)
       if (zdir.gt.zero) zlim = blim(kex)
       aulim = min(aulim,(zlim-zmid)/zdir)
    end do
!                  lsb = line search buffer
!          first point
    anext = 0.
    aopt = anext
    limset = .false.
    if (aulim.lt.aopt+tla)  limset = .true.
    call mneval(anext,fnext,ierev)
! debug printout:
    if (ldebug) write (minuit_fileout,'(A,I8,A,F10.5,A,2F10.5)')&
         & ' MNCROS: CALLS=',nfcn,'   AIM=',aim,'  F,A=',fnext,aopt
    if (ierev.gt.0)  go to 900
    if (limset.and.fnext.le.aim)  go to 930
    ipt = ipt + 1
    xpt(ipt) = anext
    ypt(ipt) = fnext
    chpt(ipt)= charal(ipt:ipt)
    alsb(1) = anext
    flsb(1) = fnext
    fnext = max(fnext,aminsv+0.1*up)
    aopt =  sqrt((up)/(fnext-aminsv)) - 1.0
    if (abs(fnext-aim).lt.tlf)  go to 800
!
    if (aopt.lt.-half)  aopt = -half
    if (aopt.gt.one)    aopt = one
    limset = .false.
    if (aopt.gt.aulim)  then
       aopt = aulim
       limset = .true.
    endif
    call mneval(aopt,fnext,ierev)
! debug printout:
    if (ldebug) write (minuit_fileout,'(A,I8,A,F10.5,A,2F10.5)')&
         & ' MNCROS: CALLS=',nfcn,'   AIM=',aim,'  F,A=',fnext,aopt
    if (ierev.gt.0)  go to 900
    if (limset.and.fnext.le.aim)  go to 930
    alsb(2) = aopt
    ipt = ipt + 1
    xpt(ipt) = alsb(2)
    ypt(ipt) = fnext
    chpt(ipt)= charal(ipt:ipt)
    flsb(2) = fnext
    dfda = (flsb(2)-flsb(1))/ (alsb(2)-alsb(1))
!                   dfda must be positive on the contour
    if (dfda.gt.zero)  go to 460
300 call mnwarn('D',chere,'LOOKING FOR SLOPE OF THE RIGHT SIGN')
    maxlk = maxitr - ipt
    do it= 1, maxlk
       alsb(1) = alsb(2)
       flsb(1) = flsb(2)
       aopt = alsb(1) + 0.2*real(it)
       limset = .false.
       if (aopt.gt.aulim)  then
          aopt = aulim
          limset = .true.
       endif
       call mneval(aopt,fnext,ierev)
! debug printout:
       if (ldebug) write (minuit_fileout,'(A,I8,A,F10.5,A,2F10.5)')&
            & ' MNCROS: CALLS=',nfcn,'   AIM=',aim,'  F,A=',fnext,aopt
       if (ierev.gt.0)  go to 900
       if (limset.and.fnext.le.aim)  go to 930
       alsb(2) = aopt
       ipt = ipt + 1
       xpt(ipt) = alsb(2)
       ypt(ipt) = fnext
       chpt(ipt)= charal(ipt:ipt)
       flsb(2) = fnext
       dfda = (flsb(2)-flsb(1))/ (alsb(2)-alsb(1))
       if (dfda.gt.zero)  go to 450
    end do
    call mnwarn('W',chere,'CANNOT FIND SLOPE OF THE RIGHT SIGN')
    go to 950
450 continue
!                    we have two points with the right slope
460 aopt = alsb(2) + (aim-flsb(2))/dfda
    fdist = min(abs(aim -flsb(1)),abs(aim -flsb(2)))
    adist = min(abs(aopt-alsb(1)),abs(aopt-alsb(2)))
    tla = tlr
    if (abs(aopt).gt.one)  tla = tlr*abs(aopt)
    if (adist.lt.tla.and.fdist.lt.tlf) go to 800
    if (ipt.ge.maxitr)  go to 950
    bmin = min(alsb(1),alsb(2)) - 1.0
    if (aopt.lt.bmin)  aopt = bmin
    bmax = max(alsb(1),alsb(2)) + 1.0
    if (aopt.gt.bmax)  aopt = bmax
!                    try a third point
    limset = .false.
    if (aopt.gt.aulim) then
       aopt = aulim
       limset = .true.
    endif
    call mneval(aopt,fnext,ierev)
! debug printout:
    if (ldebug) write (minuit_fileout,'(A,I8,A,F10.5,A,2F10.5)')&
         & ' MNCROS: CALLS=',nfcn,'   AIM=',aim,'  F,A=',fnext,aopt
    if (ierev.gt.0)  go to 900
    if (limset.and.fnext.le.aim)  go to 930
    alsb(3) = aopt
    ipt = ipt + 1
    xpt(ipt) = alsb(3)
    ypt(ipt) = fnext
    chpt(ipt)= charal(ipt:ipt)
    flsb(3) = fnext
    inew = 3
!                now we have three points, ask how many <aim
    ecarmn = abs(fnext-aim)
    ibest = 3
    ecarmx = 0.
    noless = 0
    do i= 1, 3
       ecart = abs(flsb(i) - aim)
       if (ecart.gt.ecarmx) then
          ecarmx = ecart
          iworst = i
       endif
       if (ecart.lt.ecarmn) then
          ecarmn = ecart
          ibest = i
       endif
       if (flsb(i).lt.aim) noless = noless + 1
    end do
    inew = ibest
!           if at least one on each side of aim, fit a parabola
    if (noless.eq.1.or.noless.eq.2) go to 500
!           if all three are above aim, third must be closest to aim
    if (noless.eq.0.and.ibest.ne.3)  go to 950
!           if all three below, and third is not best, then slope
!             has again gone negative, look for positive slope.
    if (noless.eq.3.and.ibest.ne.3) then
       alsb(2) = alsb(3)
       flsb(2) = flsb(3)
       go to 300
    endif
!           in other cases, new straight line thru last two points
    alsb(iworst) = alsb(3)
    flsb(iworst) = flsb(3)
    dfda = (flsb(2)-flsb(1))/ (alsb(2)-alsb(1))
    go to 460
!                parabola fit
500 call mnpfit(alsb,flsb,3,coeff,sdev)
    if (coeff(3).le.zero)  call mnwarn ('D',chere,&
         &             'CURVATURE IS NEGATIVE NEAR CONTOUR LINE.')
    determ =  coeff(2)**2 - 4.*coeff(3)*(coeff(1)-aim)
    if (determ.le.zero)   then
       call mnwarn('D',chere,'PROBLEM 2, IMPOSSIBLE DETERMINANT')
       go to 950
    endif
!                find which root is the right one
    rt = sqrt(determ)
    x1 = (-coeff(2) + rt)/(2.*coeff(3))
    x2 = (-coeff(2) - rt)/(2.*coeff(3))
    s1 = coeff(2) + 2.*x1*coeff(3)
    s2 = coeff(2) + 2.*x2*coeff(3)
    if (s1*s2.gt.zero) write (minuit_fileout,'(A)') ' MNCONTOUR PROBLEM 1'
    aopt = x1
    slope = s1
    if (s2.gt.zero)  then
       aopt = x2
       slope = s2
    endif
!         ask if converged
    tla = tlr
    if (abs(aopt).gt.one)  tla = tlr*abs(aopt)
    if (abs(aopt-alsb(ibest)).lt.tla.and.&
         &    abs(flsb(ibest)-aim).lt.tlf)  go to 800
    if (ipt.ge.maxitr)  go to 950
!         see if proposed point is in acceptable zone between l and r
!         first find ileft, iright, iout and ibest
    ileft = 0
    iright = 0
    ibest = 1
    ecarmx = 0.
    ecarmn = abs(aim-flsb(1))
    do i= 1, 3
       ecart = abs(flsb(i) - aim)
       if (ecart.lt.ecarmn) then
          ecarmn = ecart
          ibest = i
       endif
       if (ecart.gt.ecarmx) ecarmx = ecart
       if (flsb(i).gt.aim)  then
          if (iright.eq.0)  then
             iright = i
          else if (flsb(i).gt.flsb(iright)) then
             iout = i
          else
             iout = iright
             iright = i
          endif
       else if (ileft.eq.0)  then
          ileft = i
       else if (flsb(i).lt.flsb(ileft)) then
          iout = i
       else
          iout = ileft
          ileft = i
       endif
    end do
!       avoid keeping a very bad point next time around
    if (ecarmx.gt.10.*abs(flsb(iout)-aim))&
         &    aopt = half*aopt + half*half*(alsb(iright)+alsb(ileft))
!         knowing ileft and iright, get acceptable window
    smalla = 0.1*tla
    if (slope*smalla.gt.tlf)  smalla = tlf/slope
    aleft  = alsb(ileft)  + smalla
    aright = alsb(iright) - smalla
!         move proposed point aopt into window if necessary
    if (aopt.lt.aleft)  aopt = aleft
    if (aopt.gt.aright) aopt = aright
    if (aleft.gt.aright)aopt = half*(aleft + aright)
!         see if proposed point outside limits (should be impossible!)
    limset = .false.
    if (aopt.gt.aulim)  then
       aopt = aulim
       limset = .true.
    endif
!                  evaluate function at new point aopt
    call mneval(aopt,fnext,ierev)
! debug printout:
    if (ldebug) write (minuit_fileout,'(A,I8,A,F10.5,A,2F10.5)')&
         & ' MNCROS: CALLS=',nfcn,'   AIM=',aim,'  F,A=',fnext,aopt
    if (ierev.gt.0)  go to 900
    if (limset.and.fnext.le.aim)  go to 930
    ipt = ipt + 1
    xpt(ipt) = aopt
    ypt(ipt) = fnext
    chpt(ipt)= charal(ipt:ipt)
!                replace odd point by new one
    alsb(iout) = aopt
    flsb(iout) = fnext
!          the new point may not be the best, but it is the only one
!          which could be good enough to pass convergence criteria
    ibest = iout
    go to 500
!
!       contour has been located, return point to mncont or minos
800 continue
    iercr = 0
    go to 1000
!                error in the minimization
900 if (ierev.eq.1)  go to 940
    go to 950
!                parameter up against limit
930 iercr = 1
    go to 1000
!                too many calls to fcn
940 iercr = 2
    go to 1000
!                cannot find next point
950 iercr = 3
!                in any case
1000 continue
    if (ldebug) then
       itoohi = 0
       do i= 1, ipt
          if (ypt(i).gt.aim+up) then
             ypt(i) = aim+up
             chpt(i) = '+'
             itoohi = 1
          endif
       end do
       chsign = 'POSI'
       if (xdircr.lt.zero)  chsign = 'NEGA'
       if (ke2cr.eq.0)  write (minuit_fileout, '(2X,A,A,I3)')&
            &            chsign,'TIVE MINOS ERROR, PARAMETER ',ke1cr
       if (itoohi.eq.1)  write (minuit_fileout, '(10X,A)')&
            &            'POINTS LABELLED "+" WERE TOO HIGH TO PLOT.'
       if (iercr.eq.1) write (minuit_fileout,'(10X,A)')&
            &            'RIGHTMOST POINT IS UP AGAINST LIMIT.'
       call mnplot(xpt,ypt,chpt,ipt,minuit_fileout,npagwd,npagln)
    endif
    return
  end subroutine mncros

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mneval(anext,fnext,ierev)
!c      evaluates the function being analyzed by mncros, which is
!c      generally the minimum of fcn with respect to all remaining
!c      variable parameters.  common block /mn7xcr/ contains the
!c      data necessary to know the values of u(ke1cr) and u(ke2cr)
!c      to be used, namely     u(ke1cr) = xmidcr + anext*xdircr
!c      and (if ke2cr.ne.0)  u(ke2cr) = ymidcr + anext*ydircr
!c
    !
    !
    real(FREAL) :: anext, fnext
    integer(FINT) :: ierev
    !
    integer(FINT) :: nparx
    !
    u(ke1cr) = xmidcr + anext*xdircr
    if ( ke2cr.ne.0)  u(ke2cr) = ymidcr + anext*ydircr
    call mninex(x)
    nparx = npar
    call minuit_fcn(nparx,gin,fnext,u,4)
    nfcn = nfcn + 1
    ierev = 0
    if (npar.gt.0)  then
       itaur = 1
       amin = fnext
       isw(1) = 0
       call mnmigr
       itaur = 0
       fnext = amin
       if (isw(1).ge.1)  ierev = 1
       if (isw(4).lt.1)  ierev = 2
    endif
    return
  end subroutine mneval

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine MNSET
!C        Called from MNEXCM
!C        Interprets the commands that start with SET and SHOW
!C
!
!        file characteristics for SET INPUT
    integer(FINT), parameter :: NUMDBG = 6
    integer(FINT) :: i
    logical LNAME
    character CFNAME*64, CMODE*16
!       'SET ' or 'SHOW',  'ON ' or 'OFF', 'SUPPRESSED' or 'REPORTED  '
    character CKIND*4,    COPT*3,         CWARN*10
!        explanation of print level numbers -1:3  and strategies 0:2
    character CPRLEV(-1:3)*34 ,CSTRAT(0:2)*44
!        identification of debug options
    character*40 CDBOPT(0:NUMDBG)
!        things that can be set or shown
    character*10 CNAME(30)
    !
    integer(FINT) :: nname, nntot, kname, iprm, isw2, jseed, iunit
    integer(FINT) :: iset, idbopt, id, igrain, ikseed, iseed, iswsav
    integer(FINT) :: ii, kk
    real(FREAL) :: val
    !
    data CNAME( 1)/'FCN value '/
    data CNAME( 2)/'PARameters'/
    data CNAME( 3)/'LIMits    '/
    data CNAME( 4)/'COVariance'/
    data CNAME( 5)/'CORrelatio'/
    data CNAME( 6)/'PRInt levl'/
    data CNAME( 7)/'NOGradient'/
    data CNAME( 8)/'GRAdient  '/
    data CNAME( 9)/'ERRor def '/
    data CNAME(10)/'INPut file'/
    data CNAME(11)/'WIDth page'/
    data CNAME(12)/'LINes page'/
    data CNAME(13)/'NOWarnings'/
    data CNAME(14)/'WARnings  '/
    data CNAME(15)/'RANdom gen'/
    data CNAME(16)/'TITle     '/
    data CNAME(17)/'STRategy  '/
    data CNAME(18)/'EIGenvalue'/
    data CNAME(19)/'PAGe throw'/
    data CNAME(20)/'MINos errs'/
    data CNAME(21)/'EPSmachine'/
    data CNAME(22)/'OUTputfile'/
    data CNAME(23)/'BATch     '/
    data CNAME(24)/'INTeractiv'/
    data CNAME(25)/'VERsion   '/
    data NNAME/25/
!        options not intended for normal users
    data CNAME(26)/'reserve   '/
    data CNAME(27)/'NODebug   '/
    data CNAME(28)/'DEBug     '/
    data CNAME(29)/'SHOw      '/
    data CNAME(30)/'SET       '/
    data NNTOT/30/
!
    data CPRLEV(-1)/'-1: NO OUTPUT EXCEPT FROM "SHOW"  '/
    data CPRLEV( 0)/' 0: REDUCED OUTPUT                '/
    data CPRLEV( 1)/' 1: NORMAL OUTPUT                 '/
    data CPRLEV( 2)/' 2: EXTRA OUTPUT FOR PROBLEM CASES'/
    data CPRLEV( 3)/' 3: MAXIMUM OUTPUT                '/
!
    data CSTRAT( 0)/' 0: MINIMIZE THE NUMBER OF CALLS TO FUNCTION'/
    data CSTRAT( 1)/' 1: TRY TO BALANCE SPEED AGAINST RELIABILITY'/
    data CSTRAT( 2)/' 2: MAKE SURE MINIMUM TRUE, ERRORS CORRECT  '/
    !
    data CDBOPT(0)/'REPORT ALL EXCEPTIONAL CONDITIONS      '/
    data CDBOPT(1)/'MNLINE: LINE SEARCH MINIMIZATION       '/
    data CDBOPT(2)/'MNDERI: FIRST DERIVATIVE CALCULATIONS  '/
    data CDBOPT(3)/'MNHESS: SECOND DERIVATIVE CALCULATIONS '/
    data CDBOPT(4)/'MNMIGR: COVARIANCE MATRIX UPDATES      '/
    data CDBOPT(5)/'MNHES1: FIRST DERIVATIVE UNCERTAINTIES '/
    data CDBOPT(6)/'MNCONT: MNCONTOUR PLOT (MNCROS SEARCH) '/
!
!
    do I= 1, NNTOT
       if (index(CWORD(4:10),CNAME(I)(1:3)).gt.0)  goto 5
    end do
    I = 0
5   KNAME = I
!
!           Command could be SET xxx, SHOW xxx,  HELP SET or HELP SHOW
    if (index(CWORD(1:4),'HEL').gt.0)  goto 2000
    if (index(CWORD(1:4),'SHO').gt.0)  goto 1000
    if (index(CWORD(1:4),'SET').eq.0)  goto 1900
!                           ---
    CKIND = 'SET '
!                                        . . . . . . . . . . set unknown
    if (KNAME.le.0)  goto 1900
!                                        . . . . . . . . . . set known
    goto(3000,  20,  30,  40,3000,  60,  70,  80,  90, 100,&
         &110, 120, 130, 140, 150, 160, 170,3000, 190,3000,&
         &210, 220, 230, 240,3000,1900, 270, 280, 290, 300) , KNAME
!
!                                        . . . . . . . . . . set param
20  continue
    IPRM = WORD7(1)
    if (IPRM.gt.NU)  goto 25
    if (IPRM.le.0)   goto 25
    if (NVARL(IPRM).lt.0)  goto 25
    U(IPRM) = WORD7(2)
    call MNEXIN(X)
    ISW2 = ISW(2)
    call MNRSET(1)
!        Keep approximate covariance matrix, even if new param value
    ISW(2) = min(ISW2,1)
    CFROM = 'SET PARM'
    NFCNFR = NFCN
    CSTATU = 'NEW VALUES'
    goto 4000
25  write (minuit_fileout,'(A/)') ' UNDEFINED PARAMETER NUMBER.  IGNORED.'
    goto 4000
!                                        . . . . . . . . . . set limits
30  call MNLIMS
    goto 4000
!                                        . . . . . . . . . . set covar
40  continue
!   this command must be handled by MNREAD, and is not Fortran-callable
    goto 3000
!                                        . . . . . . . . . . set print
60  ISW(5) = WORD7(1)
    goto 4000
!                                        . . . . . . . . . . set nograd
70  ISW(3) = 0
    goto 4000
!                                        . . . . . . . . . . set grad
80  call MNGRAD
    goto 4000
!                                        . . . . . . . . . . set errdef
90  if (WORD7(1).eq.UP)  goto 4000
    if (WORD7(1).le.ZERO)  then
       if (UP.eq.UPDFLT)  goto 4000
       UP = UPDFLT
    else
       UP = WORD7(1)
    endif
    do I= 1, NPAR
       ERN(I) = 0.
       ERP(I) = 0.
    end do
    call MNWERR
    goto 4000
!                                        . . . . . . . . . . set input
! This command must be handled by MNREAD. If it gets this far,
!         it is illegal.
100 continue
    goto 3000
!                                        . . . . . . . . . . set width
110 NPAGWD = WORD7(1)
    NPAGWD = max(NPAGWD,50)
    goto 4000
!                                        . . . . . . . . . . set lines
120 NPAGLN = WORD7(1)
    goto 4000
!                                        . . . . . . . . . . set nowarn
130 LWARN = .false.
    goto 4000
!                                        . . . . . . . . . . set warn
140 LWARN = .true.
    call MNWARN('W','SHO','SHO')
    goto 4000
!                                        . . . . . . . . . . set random
150 JSEED = int(WORD7(1))
    VAL = 3.
    call MNRN15(VAL, JSEED)
    if (ISW(5).gt.0) write (minuit_fileout, 151) JSEED
151 format (' MINUIT RANDOM NUMBER SEED SET TO ',I10)
    goto 4000
!                                        . . . . . . . . . . set title
160 continue
!   this command must be handled by MNREAD, and is not Fortran-callable
    goto 3000
!                                        . . . . . . . . . set strategy
170 ISTRAT = WORD7(1)
    ISTRAT = max(ISTRAT,0)
    ISTRAT = min(ISTRAT,2)
    if (ISW(5).gt.0)  goto 1172
    goto 4000
!                                       . . . . . . . . . set page throw
190 NEWPAG = WORD7(1)
    goto 1190
!                                        . . . . . . . . . . set epsmac
210 if (WORD7(1).gt.ZERO.and.WORD7(1).lt.0.1) EPSMAC = WORD7(1)
    EPSMA2 = sqrt(EPSMAC)
    goto 1210
!                                        . . . . . . . . . . set outputf
220 continue
    IUNIT = WORD7(1)
    minuit_fileout = IUNIT
    ISTKWR(1) = IUNIT
    if (ISW(5).ge.0) goto 1220
    goto 4000
!                                        . . . . . . . . . . set batch
230 ISW(6) = 0
    if (ISW(5).ge.0)  goto 1100
    goto 4000
!                                        . . . . . . . . . . set interac
240 ISW(6) = 1
    if (ISW(5).ge.0)  goto 1100
    goto 4000
!                                        . . . . . . . . . . set nodebug
270 ISET = 0
    goto 281
!                                        . . . . . . . . . . set debug
280 ISET = 1
281 continue
    IDBOPT = WORD7(1)
    if (IDBOPT.gt.NUMDBG) goto 288
    if (IDBOPT.ge.0) then
       IDBG(IDBOPT) = ISET
       if (ISET.eq.1)  IDBG(0) = 1
    else
!             SET DEBUG -1  sets all debug options
       do ID= 0, NUMDBG
          IDBG(ID) = ISET
       end do
    endif
    LREPOR = (IDBG(0).ge.1)
    call MNWARN('D','SHO','SHO')
    goto 4000
288 write (minuit_fileout,289) IDBOPT
289 format (' UNKNOWN DEBUG OPTION',I6,' REQUESTED. IGNORED')
    goto 4000
!                                        . . . . . . . . . . set show
290 continue
!                                        . . . . . . . . . . set set
300 continue
    goto 3000
!                -----------------------------------------------------
1000 continue
!               at this point, CWORD must be 'SHOW'
    CKIND = 'SHOW'
    if (KNAME.le.0)  goto 1900
    goto (1010,1020,1030,1040,1050,1060,1070,1070,1090,1100,&
         &1110,1120,1130,1130,1150,1160,1170,1180,1190,1200,&
         &1210,1220,1100,1100,1250,1900,1270,1270,1290,1300),KNAME
!
!                                        . . . . . . . . . . show fcn
1010 continue
    if (AMIN.eq.UNDEFI)  call MNAMIN
    call MNPRIN (0,AMIN)
    goto 4000
!                                        . . . . . . . . . . show param
1020 continue
    if (AMIN.eq.UNDEFI)  call MNAMIN
    call MNPRIN (5,AMIN)
    goto 4000
!                                        . . . . . . . . . . show limits
1030 continue
    if (AMIN.eq.UNDEFI)  call MNAMIN
    call MNPRIN (1,AMIN)
    goto 4000
!                                        . . . . . . . . . . show covar
1040 call MNMATU(1)
    goto 4000
!                                        . . . . . . . . . . show corre
1050 call MNMATU(0)
    goto 4000
!                                        . . . . . . . . . . show print
1060 continue
    if (ISW(5).lt.-1)  ISW(5) = -1
    if (ISW(5).gt.3)  ISW(5) = 3
    write (minuit_fileout,'(A)') ' ALLOWED PRINT LEVELS ARE:'
    write (minuit_fileout,'(27X,A)') CPRLEV
    write (minuit_fileout,1061)  CPRLEV(ISW(5))
1061 format (/' CURRENT PRINTOUT LEVEL IS ',A)
    goto 4000
!                                        . . . . . . . show nograd, grad
1070 continue
    if (ISW(3).le.0) then
       write (minuit_fileout, 1081)
1081   format(' NOGRAD IS SET.  DERIVATIVES NOT COMPUTED IN FCN.')
    else
       write (minuit_fileout, 1082)
1082   format('   GRAD IS SET.  USER COMPUTES DERIVATIVES IN FCN.')
    endif
    goto 4000
!                                       . . . . . . . . . . show errdef
1090 write (minuit_fileout, 1091)  UP
1091 format (' ERRORS CORRESPOND TO FUNCTION CHANGE OF',G13.5)
    goto 4000
!                                       . . . . . . . . . . show input,
!                                                batch, or interactive
1100 continue
    inquire(UNIT=minuit_filein,NAMED=LNAME,NAME=CFNAME)
    CMODE = 'BATCH MODE      '
    if (ISW(6).eq.1)  CMODE = 'INTERACTIVE MODE'
    if (.not.LNAME)  CFNAME='unknown'
    write (minuit_fileout,1002) CMODE,minuit_filein,CFNAME
1002 format (' INPUT NOW BEING READ IN ',A,' FROM UNIT NO.',I3/&
          & ' FILENAME: ',A)
    goto 4000
!                                       . . . . . . . . . . show width
1110 write (minuit_fileout,1111) NPAGWD
1111 format (10X,'PAGE WIDTH IS SET TO',I4,' COLUMNS')
    goto 4000
!                                       . . . . . . . . . . show lines
1120 write (minuit_fileout,1121) NPAGLN
1121 format (10X,'PAGE LENGTH IS SET TO',I4,' LINES')
    goto 4000
!                                       . . . . . . .show nowarn, warn
1130 continue
    CWARN = 'SUPPRESSED'
    if (LWARN) CWARN = 'REPORTED  '
    write (minuit_fileout,1141) CWARN
1141 format (' MINUIT WARNING MESSAGES ARE ',A)
    if (.not.LWARN) call MNWARN('W','SHO','SHO')
    goto 4000
!                                      . . . . . . . . . . show random
1150 VAL = 0.
    call MNRN15(VAL,IGRAIN)
    IKSEED = IGRAIN
    write (minuit_fileout, 1151)  IKSEED
1151 format (' MINUIT RNDM SEED IS CURRENTLY=',I10/)
    VAL = 3.0
    ISEED = IKSEED
    call MNRN15(VAL,ISEED)
    goto 4000
!                                        . . . . . . . . . show title
1160 write (minuit_fileout,'(A,A)') ' TITLE OF CURRENT TASK IS:',CTITL
    goto 4000
!                                        . . . . . . . show strategy
1170 write (minuit_fileout, '(A)') ' ALLOWED STRATEGIES ARE:'
    write (minuit_fileout, '(20X,A)') CSTRAT
1172 write (minuit_fileout, 1175) CSTRAT(ISTRAT)
1175 format (/' NOW USING STRATEGY ',A/)
    goto 4000
!                                          . . . . . show eigenvalues
1180 continue
    ISWSAV = ISW(5)
    ISW(5) = 3
    if (ISW(2).lt.1)  then
       write (minuit_fileout,'(1X,A)') COVMES(0)
    else
       call MNPSDF
    endif
    ISW(5) = ISWSAV
    goto 4000
!                                            . . . . . show page throw
1190 write (minuit_fileout,'(A,I3)') ' PAGE THROW CARRIAGE CONTROL =',NEWPAG
    if (NEWPAG.eq.0)&
         &    write (minuit_fileout,'(A)') ' NO PAGE THROWS IN MINUIT OUTPUT'
    goto 4000
!                                        . . . . . . show minos errors
1200 continue
    do II= 1, NPAR
       if (ERP(II).gt.ZERO.or.ERN(II).lt.ZERO)  goto 1204
    end do
    write (minuit_fileout,'(A)')&
         &   '       THERE ARE NO MINOS ERRORS CURRENTLY VALID.'
    goto 4000
1204 continue
    call MNPRIN(4,AMIN)
    goto 4000
!                                        . . . . . . . . . show epsmac
1210 write (minuit_fileout,'(A,E12.3)')&
          &  ' FLOATING-POINT NUMBERS ASSUMED ACCURATE TO',EPSMAC
    goto 4000
!                                        . . . . . . show outputfiles
1220 continue
    write (minuit_fileout,'(A,I4)') '  MINUIT PRIMARY OUTPUT TO UNIT',minuit_fileout
    goto 4000
!                                        . . . . . . show version
1250 continue
    write (minuit_fileout,'(A,A)') ' THIS IS MINUIT VERSION:',CVRSN
    goto 4000
!                                        . . . . . . show nodebug, debug
1270 continue
    do ID= 0, NUMDBG
       COPT = 'OFF'
       if (IDBG(ID).ge.1)  COPT = 'ON '
       write (minuit_fileout,1286) ID, COPT, CDBOPT(ID)
    end do
1286 format (10X,'DEBUG OPTION',I3,' IS ',A3,' :',A)
    if (.not.LREPOR) call MNWARN('D','SHO','SHO')
    goto 4000
!                                        . . . . . . . . . . show show
1290 CKIND = 'SHOW'
    goto 2100
!                                        . . . . . . . . . . show set
1300 CKIND = 'SET '
    goto 2100
!                -----------------------------------------------------
!                              UNKNOWN COMMAND
1900 write (minuit_fileout, 1901) CWORD
1901 format (' THE COMMAND:',A10,' IS UNKNOWN.'/)
    goto 2100
!                -----------------------------------------------------
!                    HELP SHOW,  HELP SET,  SHOW SET, or SHOW SHOW
2000 CKIND = 'SET '
    if (index(CWORD(4:10),'SHO').gt.0)  CKIND = 'SHOW'
2100 write (minuit_fileout, 2101)  CKIND,CKIND, (CNAME(KK),KK=1,NNAME)
2101 format (' THE FORMAT OF THE ',A4,' COMMAND IS:'//&
          &   1X,A4,' xxx    [numerical arguments if any]'//&
          &   ' WHERE xxx MAY BE ONE OF THE FOLLOWING:'/&
          &   (7X,6A12))
    goto 4000
!                -----------------------------------------------------
!                               ILLEGAL COMMAND
3000 write (minuit_fileout,'('' ABOVE COMMAND IS ILLEGAL.   IGNORED'')')
4000 return
  end subroutine MNSET

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnlims
!c       called from mnset
!c       interprets the set lim command, to reset the parameter limits
!c
!
    !
    integer(FINT) :: i2, newcod, inu, ifx, kint
    real(FREAL) :: dxdi, snew
    !
    cfrom = 'SET LIM '
    nfcnfr = nfcn
    cstatu= 'NO CHANGE '
    i2 = word7(1)
    if (i2.gt.maxext.or.i2.lt.0)  go to 900
    if (i2.gt.0)  go to 30
!                                     set limits on all parameters
    newcod = 4
    if (word7(2).eq.word7(3))  newcod = 1
    do inu= 1, nu
       if (nvarl(inu).le.0) cycle
       if (nvarl(inu).eq.1.and.newcod.eq.1) cycle
       kint = niofex(inu)
!             see if parameter has been fixed
       if (kint.le.0)  then
          if (isw(5).ge.0)  write (minuit_fileout,'(11X,A,I3)')&
               &      ' LIMITS NOT CHANGED FOR FIXED PARAMETER:',inu
          cycle
       endif
       if (newcod.eq.1)  then
!            remove limits from parameter
          if (isw(5).gt.0)     write (minuit_fileout,134)  inu
          cstatu = 'NEW LIMITS'
          call mndxdi(x(kint),kint,dxdi)
          snew = gstep(kint)*dxdi
          gstep(kint) = abs(snew)
          nvarl(inu) = 1
       else
!             put limits on parameter
          alim(inu) = min(word7(2),word7(3))
          blim(inu) = max(word7(2),word7(3))
          if (isw(5).gt.0) write (minuit_fileout,237)  inu,alim(inu),blim(inu)
          nvarl(inu) = 4
          cstatu = 'NEW LIMITS'
          gstep(kint) = -0.1
       endif
    end do
    go to 900
!                                       set limits on one parameter
30  if (nvarl(i2).le.0)  then
       write (minuit_fileout,'(A,I3,A)') ' PARAMETER ',i2,' IS NOT VARIABLE.'
       go to 900
    endif
    kint = niofex(i2)
!                                       see if parameter was fixed
    if (kint.eq.0)  then
       write (minuit_fileout,'(A,I3)')&
            &     ' REQUEST TO CHANGE LIMITS ON FIXED PARAMETER:',i2
       do ifx= 1, npfix
          if (i2.eq.ipfix(ifx)) go to 92
       end do
       write (minuit_fileout,'(A)') ' MINUIT BUG IN MNLIMS. SEE F. JAMES'
92     continue
    endif
    if (word7(2).ne.word7(3))  go to 235
!                                       remove limits
    if (nvarl(i2).ne.1)  then
       if (isw(5).gt.0)  write (minuit_fileout,134)  i2
134    format (30h limits removed from parameter  ,i4)
       cstatu = 'NEW LIMITS'
       if (kint.le.0)  then
          gsteps(ifx) = abs(gsteps(ifx))
       else
          call mndxdi(x(kint),kint,dxdi)
          if (abs(dxdi).lt.0.01)  dxdi=0.01
          gstep(kint) = abs(gstep(kint)*dxdi)
          grd(kint) = grd(kint)*dxdi
       endif
       nvarl(i2) = 1
    else
       write (minuit_fileout,'(A,I3)') ' NO LIMITS SPECIFIED.  PARAMETER ',&
            &        i2,' IS ALREADY UNLIMITED.  NO CHANGE.'
    endif
    go to 900
!                                        put on limits
235 alim(i2) = min(word7(2),word7(3))
    blim(i2) = max(word7(2),word7(3))
    nvarl(i2) = 4
    if (isw(5).gt.0)   write (minuit_fileout,237)  i2,alim(i2),blim(i2)
237 format (10h parameter ,i3, 14h limits set to  ,2g15.5)
    cstatu = 'NEW LIMITS'
    if (kint.le.0)  then
       gsteps(ifx) = -0.1
    else
       gstep(kint) = -0.1
    endif
!
900 continue
    if (cstatu.ne.'NO CHANGE ')  then
       call mnexin(x)
       call mnrset(1)
    endif
    return
  end subroutine mnlims

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mngrad
!c       called from mnset
!c       interprets the set grad command, which informs minuit whether
!c       the first derivatives of fcn will be calculated by the user
!c       inside fcn.  it can check the user's derivative calculation
!c       by comparing it with a finite difference approximation.
!c
!
    !
    !
    character(4), parameter :: cgood='GOOD',cbad=' BAD',cnone='NONE'
    character*4 cwd
    logical lnone
    integer(FINT) :: i
    real(FREAL), dimension(mni) :: gf !(mni)
    !
    integer(FINT) :: nparx, istsav, lc
    real(FREAL) :: fzero, err
    !
!
    isw(3) = 1
    nparx = npar
    if (word7(1).gt.zero)  go to 2000
!                  get user-calculated first derivatives from fcn
    do i= 1, nu
       gin(i) = undefi
    end do
    call mninex(x)
    call minuit_fcn(nparx,gin,fzero,u,2)
    nfcn = nfcn + 1
    call mnderi
    do i= 1, npar
       gf(i) = grd(i)
    end do
!                    get minuit-calculated first derivatives
    isw(3) = 0
    istsav = istrat
    istrat = 2
    call mnhes1
    istrat = istsav
    write (minuit_fileout,51)
51  format(/' CHECK OF GRADIENT CALCULATION IN FCN'/12X,'PARAMETER',&
         & 6X,9HG(IN FCN) ,3X,9HG(MINUIT) ,2X,'DG(MINUIT)',3X,9HAGREEMENT)
    isw(3) = 1
    lnone = .false.
    do lc = 1, npar
       i = nexofi(lc)
       cwd = cgood
       err = dgrd(lc)
       if (abs(gf(lc)-grd(lc)).gt.err)  cwd = cbad
       if (gin(i).eq.undefi)  then
          cwd = cnone
          lnone = .true.
          gf(lc) = 0.
       endif
       if (cwd.ne.cgood)  isw(3) = 0
       write (minuit_fileout,99) i,cpnam(i),gf(lc),grd(lc),err,cwd
99     format (7x,i5,2x ,a10,3e12.4,4x ,a4)
    end do
    if (lnone) write (minuit_fileout,'(A)')&
         &  '  AGREEMENT=NONE  MEANS FCN DID NOT CALCULATE THE DERIVATIVE'
    if (isw(3).eq.0)  write (minuit_fileout,1003)
1003 format(/' MINUIT DOES NOT ACCEPT DERIVATIVE CALCULATIONS BY FCN'/&
          & ' TO FORCE ACCEPTANCE, ENTER "SET GRAD    1"'/)
!
2000 continue
    return
  end subroutine mngrad

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnscan
!c        scans the values of fcn as a function of one parameter
!c        and plots the resulting values as a curve using mnplot.
!c        it may be called to scan one parameter or all parameters.
!c        retains the best function and parameter values found.
    !
    !
    integer(FINT) :: ncall, nccall, iparwd, ipar, liint, nxypt
    integer(FINT) :: icall, nparx, nbins, nunit
    real(FREAL) :: xlreq, xhreq, ubest, unext, step, xl, xh, uhigh, fnext
    !
    xlreq = min(word7(3),word7(4))
    xhreq = max(word7(3),word7(4))
    ncall = word7(2) + 0.01
    if (ncall.le.1)  ncall = 41
    if (ncall.gt.maxcpt)  ncall = maxcpt
    nccall = ncall
    if (amin.eq.undefi)  call mnamin
    iparwd = word7(1) + 0.1
    ipar = max(iparwd, 0)
    liint = niofex(ipar)
    cstatu = 'NO CHANGE'
    if (iparwd.gt.0)  go to 200
!
!         equivalent to a loop over parameters requested
100 ipar = ipar + 1
    if (ipar.gt.nu)  go to 900
    liint = niofex(ipar)
    if (liint.le.0)  go to 100
!         set up range for parameter ipar
200 continue
    ubest = u(ipar)
    xpt(1) = ubest
    ypt(1) = amin
    chpt(1)= ' '
    xpt(2) = ubest
    ypt(2) = amin
    chpt(2)= 'X'
    nxypt = 2
    if (nvarl(ipar).gt.1)  go to 300
!         no limits on parameter
    if (xlreq.eq.xhreq)  go to 250
    unext = xlreq
    step = (xhreq-xlreq)/float(ncall-1)
    go to 500
250 continue
    xl = ubest - werr(liint)
    xh = ubest+  werr(liint)
    call mnbins(xl,xh,ncall, unext,uhigh,nbins,step)
    nccall = nbins + 1
    go to 500
!         limits on parameter
300 continue
    if (xlreq.eq.xhreq)  go to 350
    xl = max(xlreq,alim(ipar))
    xh = min(xhreq,blim(ipar))
    if (xl.ge.xh)  go to 700
    unext = xl
    step = (xh-xl)/float(ncall-1)
    go to 500
350 continue
    unext = alim(ipar)
    step = (blim(ipar)-alim(ipar))/float(ncall-1)
!         main scanning loop over parameter ipar
500 continue
    do icall = 1, nccall
       u(ipar) = unext
       nparx = npar
       call minuit_fcn(nparx,gin,fnext,u,4)
       nfcn = nfcn + 1
       nxypt = nxypt + 1
       xpt(nxypt) = unext
       ypt(nxypt) = fnext
       chpt(nxypt) = '*'
       if (fnext.lt.amin)  then
          amin = fnext
          ubest = unext
          cstatu= 'IMPROVED  '
       endif
530    continue
       unext = unext + step
    end do
!         finished with scan of parameter ipar
    u(ipar) = ubest
    call mnexin(x)
    if (isw(5).ge.1)  then
       write (minuit_fileout,1001)  newpag,ipar,cpnam(ipar)
       nunit = minuit_fileout
       call mnplot(xpt,ypt,chpt,nxypt,nunit,npagwd,npagln)
    endif
    go to 800
700 continue
    write (minuit_fileout,1000) ipar
800 continue
    if (iparwd.le.0)  go to 100
!         finished with all parameters
900 continue
    if (isw(5).ge.0) call mnprin(5,amin)
    return
1000 format (46H REQUESTED RANGE OUTSIDE LIMITS FOR PARAMETER  ,I3/)
1001 format (I1,'SCAN OF PARAMETER NO.',I3,3H,   ,A10)
  end subroutine mnscan

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mncntr(ke1,ke2,ierrf)
!c       to print function contours in two variables, on line printer
!c
    !
    !
    integer(FINT), parameter :: numbcs=20,nxmax=115
    !
    integer(FINT) :: ke1, ke2, ierrf
    !
    integer(FINT) :: i
    real(FREAL), dimension(numbcs) :: contur !(numbcs)
    real(FREAL), dimension(nxmax) :: fcna, fcnb !(nxmax)
    character clabel*(numbcs)
    character chln*(nxmax),chmid*(nxmax),chzero*(nxmax)
    data clabel/'0123456789ABCDEFGHIJ'/
    !
    integer(FINT) :: ki1, ki2, nparx, ngrid, nx, ny, ixmid, ixzero
    integer(FINT) :: ix, iy, ics, nl, nl2
    real(FREAL) :: xsav, ysav, devs, xlo, xup, ylo, yup, bwidx, bwidy
    real(FREAL) :: xb4, ff, unext, ylabel, fmx, fmn
    !
!                 input arguments: parx, pary, devs, ngrid
    if (ke1.le.0.or.ke2.le.0)  go to 1350
    if (ke1.gt.nu.or.ke2.gt.nu)  go to 1350
    ki1 = niofex(ke1)
    ki2 = niofex(ke2)
    if (ki1.le.0.or.ki2.le.0)  go to 1350
    if (ki1.eq.ki2)  go to 1350
!
    if (isw(2).lt.1)  then
       call mnhess
       call mnwerr
    endif
    nparx = npar
    xsav = u(ke1)
    ysav = u(ke2)
    devs = word7(3)
    if (devs.le.zero)  devs=2.
    xlo = u(ke1) - devs*werr(ki1)
    xup = u(ke1) + devs*werr(ki1)
    ylo = u(ke2) - devs*werr(ki2)
    yup = u(ke2) + devs*werr(ki2)
    ngrid = word7(4)
    if (ngrid.le.0)  then
       ngrid=25
       nx = min(npagwd-15,ngrid)
       ny = min(npagln-7, ngrid)
    else
       nx = ngrid
       ny = ngrid
    endif
    if (nx.lt.11) nx=11
    if (ny.lt.11) ny=11
    if (nx.ge.nxmax)  nx=nxmax-1
!         ask if parameter outside limits
    if (nvarl(ke1).gt.1)  then
       if (xlo.lt.alim(ke1))  xlo = alim(ke1)
       if (xup.gt.blim(ke1))  xup = blim(ke1)
    endif
    if (nvarl(ke2).gt.1)   then
       if (ylo.lt.alim(ke2))  ylo = alim(ke2)
       if (yup.gt.blim(ke2))  yup = blim(ke2)
    endif
    bwidx = (xup-xlo)/real(nx)
    bwidy = (yup-ylo)/real(ny)
    ixmid = int((xsav-xlo)*real(nx)/(xup-xlo)) + 1
    if (amin.eq.undefi)  call mnamin
    do i= 1, numbcs
       contur(i) = amin + up*float(i-1)**2
    end do
    contur(1) = contur(1) + 0.01*up
!                fill fcnb to prepare first row, and find column zero
    u(ke2) = yup
    ixzero = 0
    xb4 = one
    do ix= 1, nx+1
       u(ke1) = xlo + real(ix-1)*bwidx
       call minuit_fcn(nparx,gin,ff,u,4)
       fcnb(ix) = ff
       if (xb4.lt.zero.and.u(ke1).gt.zero)  ixzero = ix-1
       xb4 = u(ke1)
       chmid(ix:ix) = '*'
       chzero(ix:ix)= '-'
    end do
    write (minuit_fileout,'(A,I3,A,A)') ' Y-AXIS: PARAMETER ',&
         &      ke2,': ',cpnam(ke2)
    if (ixzero.gt.0)  then
       chzero(ixzero:ixzero) = '+'
       chln = ' '
       write (minuit_fileout,'(12X,A,A)') chln(1:ixzero),'X=0'
    endif
!                 loop over rows
    do iy= 1, ny
       unext = u(ke2) - bwidy
!                 prepare this line's background pattern for contour
       chln = ' '
       chln(ixmid:ixmid) = '*'
       if (ixzero.ne.0) chln(ixzero:ixzero) = ':'
       if (u(ke2).gt.ysav.and.unext.lt.ysav) chln=chmid
       if (u(ke2).gt.zero.and.unext.lt.zero) chln=chzero
       u(ke2) = unext
       ylabel = u(ke2) + 0.5*bwidy
!                 move fcnb to fcna and fill fcnb with next row
       do ix= 1, nx+1
          fcna(ix) = fcnb(ix)
          u(ke1) = xlo + real(ix-1)*bwidx
          call minuit_fcn(nparx,gin,ff,u,4)
          fcnb(ix) = ff
       end do
!                 look for contours crossing the fcnxy squares
       do ix= 1, nx
          fmx = max(fcna(ix),fcnb(ix),fcna(ix+1),fcnb(ix+1))
          fmn = min(fcna(ix),fcnb(ix),fcna(ix+1),fcnb(ix+1))
          do ics= 1, numbcs
             if (contur(ics).gt.fmn)  go to 240
          end do
          cycle
240       if (contur(ics).lt.fmx) chln(ix:ix)=clabel(ics:ics)
       end do
!                 print a row of the contour plot
       write (minuit_fileout,'(1X,G12.4,1X,A)') ylabel,chln(1:nx)
    end do
!                 contours printed, label x-axis
    chln = ' '
    chln( 1: 1) = 'I'
    chln(ixmid:ixmid) = 'I'
    chln(nx:nx) = 'I'
    write (minuit_fileout,'(14X,A)') chln(1:nx)
!                the hardest of all: print x-axis scale!
    chln = ' '
    if (nx.le.26) then
       nl = max(nx-12,2)
       nl2 = nl/2
       write (minuit_fileout,'(8X,G12.4,A,G12.4)') xlo,chln(1:nl),xup
       write (minuit_fileout,'(14X,A,G12.4)')   chln(1:nl2),xsav
    else
       nl = max(nx-24,2)/2
       nl2 = nl
       if (nl.gt.10) nl2=nl-6
       write (minuit_fileout,'(8X,G12.4,A,G12.4,A,G12.4)')  xlo,&
            &      chln(1:nl),xsav,chln(1:nl2),xup
    endif
    write (minuit_fileout,'(6X,A,I3,A,A,A,G12.4)') ' X-AXIS: PARAMETER',&
         &    ke1,': ',cpnam(ke1),'  ONE COLUMN=',bwidx
    write (minuit_fileout,'(A,G12.4,A,G12.4,A)') ' FUNCTION VALUES: F(I)=',&
         &    amin,' +',up,' *I**2'
!                 finished.  reset input values
    u(ke1) = xsav
    u(ke2) = ysav
    ierrf = 0
    return
1350 write (minuit_fileout,1351)
1351 format (' INVALID PARAMETER NUMBER(S) REQUESTED.  IGNORED.' /)
    ierrf = 1
    return
  end subroutine mncntr

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnsave
!c       writes current parameter values and step sizes onto file minuit_filesave
!c          in format which can be reread by minuit for restarting.
!c       the covariance matrix is also output if it exists.
!c
    !
    logical lopen,lname
    character cgname*64, cfname*64, canswr*1
    !
    integer(FINT) :: i
    integer(FINT) :: nlines, liint, npar2, ncovar
    !
!
    inquire(unit=minuit_filesave,opened=lopen,named=lname,name=cgname)
    if (lopen) then
       if (.not.lname) cgname='UNNAMED FILE'
       write (minuit_fileout,32) minuit_filesave,cgname
32     format (' CURRENT VALUES WILL BE SAVED ON UNIT',I3,': ',A/)
    else
!                new file, open it
       write (minuit_fileout,35) minuit_filesave
35     format (' UNIT',I3,' IS NOT OPENED.')
       if (isw(6).eq.1) then
          write (minuit_fileout,'(A)') ' PLEASE GIVE FILE NAME:'
          read (minuit_filein,'(A)') cfname
          open (unit=minuit_filesave,file=cfname,status='NEW',err=600)
          cgname = cfname
       else
          go to 650
       endif
    endif
!                               file is now correctly opened
    if (isw(6).eq.1)  then
       write (minuit_fileout,37)  minuit_filesave
37     format (' SHOULD UNIT',I3,' BE REWOUND BEFORE WRITING TO IT?' )
       read  (minuit_filein,'(A)')  canswr
       if (canswr.eq.'Y'.or.canswr.eq.'y') rewind(minuit_filesave)
    endif
!                               and rewound if requested
    write (minuit_filesave,'(10HSET TITLE )',err=700)
    write (minuit_filesave,'(A)')  ctitl
    write (minuit_filesave,'(10HPARAMETERS)')
    nlines = 3
!                                write out parameter values
    do i= 1, nu
       if (nvarl(i).lt.0) cycle
       nlines = nlines + 1
       liint = niofex(i)
       if (nvarl(i).gt.1)  go to 100
!         parameter without limits
       write (minuit_filesave,1001)  i,cpnam(i),u(i),werr(liint)
       cycle
!         parameter with limits
100    continue
       write (minuit_filesave,1001) i,cpnam(i),u(i),werr(liint),alim(i),blim(i)
1001   format (1X,I5,1H',A10,1H',4E13.5)
    end do
    write (minuit_filesave,'(A)')  ' '
    nlines = nlines + 1
!                                  write out covariance matrix, if any
    if (isw(2).lt.1)  go to 750
    write (minuit_filesave,1003,err=700)  npar
1003 format ('SET COVARIANCE',I6)
    npar2 = npar*(npar+1)/2
    write (minuit_filesave,1004) (vhmat(i),i=1,npar2)
1004 format (bn,7e11.4,3x)
    ncovar = npar2/7 + 1
    if (mod(npar2,7).gt.0)  ncovar = ncovar + 1
    nlines = nlines + ncovar
    write (minuit_fileout, 501) nlines,minuit_filesave,cgname(1:45)
501 format (1X,I5,' RECORDS WRITTEN TO UNIT',I4,':',A)
    if (ncovar.gt.0) write (minuit_fileout, 502) ncovar
502 format (' INCLUDING',I5,' RECORDS FOR THE COVARIANCE MATRIX.'/)
    go to 900
!                                           some error conditions
600 write (minuit_fileout,'(A,I4)') ' I/O ERROR: UNABLE TO OPEN UNIT',minuit_filesave
    go to 900
650 write (minuit_fileout,'(A,I4,A)') ' UNIT',minuit_filesave,' IS NOT OPENED.'
    go to 900
700 write (minuit_fileout,'(A,I4)') ' ERROR: UNABLE TO WRITE TO UNIT',minuit_filesave
    go to 900
750 write (minuit_fileout,'(A)') ' THERE IS NO COVARIANCE MATRIX TO SAVE.'
!
900 return
  end subroutine mnsave

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnimpr
!c        attempts to improve on a good local minimum by finding a
!c        better one.   the quadratic part of fcn is removed by mncalf
!c        and this transformed function is minimized using the simplex
!c        method from several random starting points.
!c        ref. -- goldstein and price, math.comp. 25, 569 (1971)
!c
    !
    !
    real(FREAL), parameter :: alpha=1., beta=0.5, gamma=2.0
    real(FREAL), dimension(mni) :: dsav !(mni)
    real(FREAL), dimension(mni+1) :: y !(mni+1)
    !
    real(FREAL) :: rnum
    !
    integer(FINT) :: i
    integer(FINT) :: npfn, nloop, nparx, nparp1, j, ndex, ifail
    integer(FINT) :: loop, iseed, jl, jh
    real(FREAL) :: wg, sigsav, reg, ycalf, amax, xi, sig2, ep, pb, ystar
    real(FREAL) :: ystst, jhold
    !
    data rnum/0./
    if (npar.le.0)  return
    if (amin.eq.undefi)  call mnamin
    cstatu = 'UNCHANGED '
    itaur = 1
    epsi = 0.1*up
    npfn=nfcn
    nloop = word7(2)
    if (nloop.le.0)  nloop = npar + 4
    nparx = npar
    nparp1=npar+1
    wg = 1.0/float(npar)
    sigsav = edm
    apsi = amin
    do i= 1, npar
       xt(i) = x(i)
       dsav(i) = werr(i)
       do j = 1, i
          ndex = i*(i-1)/2 + j
          p(i,j) = vhmat(ndex)
          p(j,i) = p(i,j)
       end do
    end do
    call mnvert(p,maxint,maxint,npar,ifail)
    if (ifail.ge.1)  go to 280
!               save inverted matrix in vt
    do i= 1, npar
       ndex = i*(i-1)/2
       do j= 1, i
          ndex = ndex + 1
          vthmat(ndex) = p(i,j)
       end do
    end do
    loop = 0
!
20  continue
    do i= 1, npar
       dirin(i) = 2.0*dsav(i)
       call mnrn15(rnum,iseed)
       x(i) = xt(i) + 2.0*dirin(i)*(rnum-0.5)
    end do
    loop = loop + 1
    reg = 2.0
    if (isw(5).ge.0)   write (minuit_fileout, 1040) loop
30  call  mncalf(x,ycalf)
    amin = ycalf
!                                        . . . . set up  random simplex
    jl = nparp1
    jh = nparp1
    y(nparp1) = amin
    amax = amin
    do i= 1, npar
       xi = x(i)
       call mnrn15(rnum,iseed)
       x(i) = xi - dirin(i) *(rnum-0.5)
       call mncalf(x,ycalf)
       y(i) = ycalf
       if (y(i).lt.amin)  then
          amin = y(i)
          jl = i
       else if (y(i).gt.amax)  then
          amax = y(i)
          jh = i
       endif
       do j= 1, npar
          p(j,i) = x(j)
       end do
       p(i,nparp1) = xi
       x(i) = xi
    end do
!
    edm = amin
    sig2 = edm
!                                        . . . . . . .  start main loop
50  continue
    if (amin.lt.zero)  go to 95
    if (isw(2).le.2)  go to 280
    ep = 0.1*amin
    if (sig2.lt.ep.and.edm.lt.ep  )     go to 100
    sig2 = edm
    if ((nfcn-npfn).gt.nfcnmx)  go to 300
!         calculate new point * by reflection
    do i= 1, npar
       pb = 0.
       do j= 1, nparp1
          pb = pb + wg * p(i,j)
       end do
       pbar(i) = pb - wg * p(i,jh)
       pstar(i)=(1.+alpha)*pbar(i)-alpha*p(i,jh)
    end do
    call mncalf(pstar,ycalf)
    ystar = ycalf
    if(ystar.ge.amin) go to 70
!         point * better than jl, calculate new point **
    do i=1,npar
       pstst(i)=gamma*pstar(i)+(1.-gamma)*pbar(i)
    end do
    call mncalf(pstst,ycalf)
    ystst = ycalf
66  if (ystst.lt.y(jl))  go to 67
    call mnrazz(ystar,pstar,y,jh,jl)
    go to 50
67  call mnrazz(ystst,pstst,y,jh,jl)
    go to 50
!         point * is not as good as jl
70  if (ystar.ge.y(jh))  go to 73
    jhold = jh
    call mnrazz(ystar,pstar,y,jh,jl)
    if (jhold.ne.jh)  go to 50
!         calculate new point **
73  do i=1,npar
       pstst(i)=beta*p(i,jh)+(1.-beta)*pbar(i)
    end do
    call mncalf(pstst,ycalf)
    ystst = ycalf
    if(ystst.gt.y(jh)) go to 30
!     point ** is better than jh
    if (ystst.lt.amin)  go to 67
    call mnrazz(ystst,pstst,y,jh,jl)
    go to 50
!                                        . . . . . .  end main loop
95  if (isw(5).ge.0)  write (minuit_fileout,1000)
    reg = 0.1
!                                        . . . . . ask if point is new
100 call mninex(x)
    call minuit_fcn(nparx,gin,amin,u,4)
    nfcn = nfcn + 1
    do i= 1, npar
       dirin(i) = reg*dsav(i)
       if (abs(x(i)-xt(i)).gt.dirin(i)) go to 150
    end do
    go to 230
150 nfcnmx = nfcnmx + npfn - nfcn
    npfn = nfcn
    call mnsimp
    if (amin.ge.apsi)  go to 325
    do i= 1, npar
       dirin(i) = 0.1 *dsav(i)
       if (abs(x(i)-xt(i)).gt.dirin(i)) go to 250
    end do
230 if (amin.lt.apsi)  go to 350
    go to 325
!                                        . . . . . . truly new minimum
250 lnewmn = .true.
    if (isw(2).ge.1) then
       isw(2) = 1
       dcovar = max(dcovar,half)
    else
       dcovar = 1.
    endif
    itaur = 0
    nfcnmx = nfcnmx + npfn - nfcn
    cstatu = 'NEW MINIMU'
    if (isw(5).ge.0)      write (minuit_fileout,1030)
    return
!                                        . . . return to previous region
280 if (isw(5).gt.0) write (minuit_fileout,1020)
    go to 325
300 isw(1) = 1
325 do i= 1, npar
       dirin(i) = 0.01*dsav(i)
       x(i) = xt(i)
    end do
    amin = apsi
    edm = sigsav
350 call mninex(x)
    if (isw(5).gt.0)    write (minuit_fileout,1010)
    cstatu= 'UNCHANGED '
    call mnrset(0)
    if (isw(2).lt.2)  go to 380
    if (loop.lt.nloop.and.isw(1).lt.1)  go to 20
380 call mnprin (5,amin)
    itaur = 0
    return
1000 format (54H AN IMPROVEMENT ON THE PREVIOUS MINIMUM HAS BEEN FOUND)
1010 format (51H IMPROVE HAS RETURNED TO REGION OF ORIGINAL MINIMUM)
1020 format (/44H COVARIANCE MATRIX WAS NOT POSITIVE-DEFINITE)
1030 format (/38H IMPROVE HAS FOUND A TRULY NEW MINIMUM/1H ,37(1H*)/)
1040 format (/18H START ATTEMPT NO.,I2,  20H TO FIND NEW MINIMUM)
  end subroutine mnimpr

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mncalf(pvec,ycalf)
!c        called only from mnimpr.  transforms the function fcn
!c        by dividing out the quadratic part in order to find further
!c        minima.    calculates  ycalf = (f-fmin)/(x-xmin)*v*(x-xmin)
!c
    !
    !
    real(FREAL), dimension(15) :: pvec !(15)
    real(FREAL) :: ycalf, f, denom
    integer(FINT) :: i
    integer(FINT) :: nparx, j, m, n, ndex
    !
    !
    nparx = npar
    call mninex(pvec)
    call minuit_fcn(nparx,gin,f,u,4)
    nfcn = nfcn + 1
    do i= 1, npar
       grd(i) = 0.
       do j= 1, npar
          m = max(i,j)
          n = min(i,j)
          ndex = m*(m-1)/2 + n
          grd(i) = grd(i) + vthmat(ndex) * (xt(j)-pvec(j))
       end do
    end do
    denom = 0.
    do i= 1, npar
       denom = denom + grd(i) * (xt(i)-pvec(i))
    end do
    if (denom.le.zero)  then
       dcovar = 1.
       isw(2) = 0
       denom = 1.0
    endif
    ycalf = (f-apsi) / denom
    return
  end subroutine mncalf

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine minuit_mncont(ke1,ke2,nptu,xptu,yptu,ierrf)
!C       Find NPTU points along a contour where the function
!C             FMIN (X(KE1),X(KE2)) =  AMIN+UP
!C       where FMIN is the minimum of FCN with respect to all
!C       the other NPAR-2 variable parameters (if any).
!C   IERRF on return will be equal to the number of points found:
!C     NPTU if normal termination with NPTU points found
!C     -1   if errors in the calling sequence (KE1, KE2 not variable)
!C      0   if less than four points can be found (using MNMNOT)
!C     n>3  if only n points can be found (n < NPTU)
!C
    !
    character(10), parameter :: CHERE='MNContour '
    !
    integer(FINT) :: i
    integer(FINT) :: ke1,ke2,nptu,ierrf
    !
    real(FREAL), dimension(nptu) :: XPTU, YPTU !(NPTU)
    real(FREAL), dimension(mni) :: W, GCC !(MNI)
    logical LDEBUG
    !
    integer(FINT) :: ki1, ki2, nfcnco, ki3, ke3, nowpts, next, nall
    integer(FINT) :: isw2, isw4, istrav, mpar, nfmxin, j, kints, ierr
    integer(FINT) :: inew, iold, i2, idist, i1, iercr, move, npcol
    integer(FINT) :: nfcol, line, lr
    real(FREAL) :: u1min, u2min, val2pl, val2mi, scalx, scaly
    real(FREAL) :: sigsav, dc, abest, bigdis, dist, a1, a2
    real(FREAL) :: xdir, ydir, sclfac, aopt
    !
!                 input arguments: parx, pary, devs, ngrid
    LDEBUG = (IDBG(6).ge.1)
    if (KE1.le.0.or.KE2.le.0)  GO TO 1350
    if (KE1.gt.NU.or.KE2.gt.NU)  GO TO 1350
    KI1 = NIOFEX(KE1)
    KI2 = NIOFEX(KE2)
    if (KI1.le.0.or.KI2.le.0)  GO TO 1350
    if (KI1.eq.KI2)  GO TO 1350
    if (NPTU.lt.4)  GO TO 1400
!
    NFCNCO = NFCN
    NFCNMX = 100*(NPTU+5)*(NPAR+1)
!           The minimum
    call MNCUVE
    U1MIN = U(KE1)
    U2MIN = U(KE2)
    IERRF = 0
    CFROM = CHERE
    NFCNFR = NFCNCO
    if (ISW(5).ge.0)  then
       write (minuit_fileout,'(1X,A,I4,A)')&
            &   'START MNCONTOUR CALCULATION OF',NPTU,' POINTS ON CONTOUR.'
       if (NPAR.gt.2) then
          if (NPAR.eq.3) then
             KI3 = 6 - KI1 - KI2
             KE3 = NEXOFI(KI3)
             write (minuit_fileout,'(1X,A,I3,2X,A)')&
                  &    'EACH POINT IS A MINIMUM WITH RESPECT TO PARAMETER ',&
                  &        KE3, CPNAM(KE3)
          else
             write (minuit_fileout,'(1X,A,I3,A)')&
                  &  'EACH POINT IS A MINIMUM WITH RESPECT TO THE OTHER',&
                  &        NPAR-2, ' VARIABLE PARAMETERS.'
          endif
       endif
    endif
!
!           Find the first four points using MNMNOT
!              ........................ first two points
    call MNMNOT(KE1,KE2,VAL2PL,VAL2MI)
    if (ERN(KI1).eq.UNDEFI)  then
       XPTU(1) = ALIM(KE1)
       call MNWARN('W',CHERE,'Contour squeezed by parameter limits.')
    else
       if (ERN(KI1).ge.ZERO)  GO TO 1500
       XPTU(1) = U1MIN+ERN(KI1)
    endif
    YPTU(1) = VAL2MI
!
    if (ERP(KI1).eq.UNDEFI)  then
       XPTU(3) = BLIM(KE1)
       call MNWARN('W',CHERE,'Contour squeezed by parameter limits.')
    else
       if (ERP(KI1).le.ZERO)  GO TO 1500
       XPTU(3) = U1MIN+ERP(KI1)
    endif
    YPTU(3) = VAL2PL
    SCALX = 1.0/(XPTU(3) - XPTU(1))
!              ........................... next two points
    call MNMNOT(KE2,KE1,VAL2PL,VAL2MI)
    if (ERN(KI2).eq.UNDEFI)  then
       YPTU(2) = ALIM(KE2)
       call MNWARN('W',CHERE,'Contour squeezed by parameter limits.')
    else
       if (ERN(KI2).ge.ZERO)  GO TO 1500
       YPTU(2) = U2MIN+ERN(KI2)
    endif
    XPTU(2) = VAL2MI
    if (ERP(KI2).eq.UNDEFI)  then
       YPTU(4) = BLIM(KE2)
       call MNWARN('W',CHERE,'Contour squeezed by parameter limits.')
    else
       if (ERP(KI2).le.ZERO)  GO TO 1500
       YPTU(4) = U2MIN+ERP(KI2)
    endif
    XPTU(4) = VAL2PL
    SCALY = 1.0/(YPTU(4) - YPTU(2))
    NOWPTS = 4
    NEXT = 5
    if (LDEBUG) then
       write (minuit_fileout,'(A)') ' Plot of four points found by MINOS'
       XPT(1) = U1MIN
       YPT(1) = U2MIN
       CHPT(1) = ' '
       NALL = min(NOWPTS+1,MAXCPT)
       do I= 2, NALL
          XPT(I) = XPTU(I-1)
          YPT(I) = YPTU(I-1)
       end do
       CHPT(2)= 'A'
       CHPT(3)= 'B'
       CHPT(4)= 'C'
       CHPT(5)= 'D'
       call MNPLOT(XPT,YPT,CHPT,NALL,minuit_fileout,NPAGWD,NPAGLN)
    endif
!
!               ..................... save some values before fixing
    ISW2 = ISW(2)
    ISW4 = ISW(4)
    SIGSAV = EDM
    ISTRAV = ISTRAT
    DC = DCOVAR
    APSI  = EPSI*0.5
    ABEST=AMIN
    MPAR=NPAR
    NFMXIN = NFCNMX
    do I= 1, MPAR
       XT(I) = X(I)
    end do
    do J= 1, MPAR*(MPAR+1)/2
       VTHMAT(J) = VHMAT(J)
    end do
    do I= 1, MPAR
       GCC(I) = GLOBCC(I)
       W(I) = WERR(I)
    end do
!                           fix the two parameters in question
    KINTS = NIOFEX(KE1)
    call MNFIXP (KINTS,IERR)
    KINTS = NIOFEX(KE2)
    call MNFIXP (KINTS,IERR)
!               ......................Fill in the rest of the points
    do INEW= NEXT, NPTU
!            find the two neighbouring points with largest separation
       BIGDIS = 0.
       do IOLD = 1, INEW-1
          I2 = IOLD + 1
          if (I2.eq.INEW) I2 = 1
          DIST = (SCALX*(XPTU(IOLD)-XPTU(I2)))**2 +                      &
               &          (SCALY*(YPTU(IOLD)-YPTU(I2)))**2
          if (DIST.gt.BIGDIS) then
             BIGDIS = DIST
             IDIST = IOLD
          endif
       end do
       I1 = IDIST
       I2 = I1 + 1
       if (I2.eq.INEW) I2 = 1
!                   next point goes between I1 and I2
       A1 = HALF
       A2 = HALF
300    XMIDCR = A1*XPTU(I1) + A2*XPTU(I2)
       YMIDCR = A1*YPTU(I1) + A2*YPTU(I2)
       XDIR = YPTU(I2) - YPTU(I1)
       YDIR = XPTU(I1) - XPTU(I2)
       SCLFAC = max(abs(XDIR*SCALX), abs(YDIR*SCALY))
       XDIRCR = XDIR/SCLFAC
       YDIRCR = YDIR/SCLFAC
       KE1CR = KE1
       KE2CR = KE2
!                Find the contour crossing point along DIR
       AMIN = ABEST
       call MNCROS(AOPT,IERCR)
       if (IERCR.gt.1)  then
!              If cannot find mid-point, try closer to point 1
          if (A1.gt.HALF) then
             if (ISW(5).ge.0)&
                  & write (minuit_fileout,'(A,A,I3,A)') ' MNCONT CANNOT FIND NEXT',&
                  &  ' POINT ON CONTOUR.  ONLY ',NOWPTS,' POINTS FOUND.'
             GO TO 950
          endif
          call MNWARN('W',CHERE,'Cannot find midpoint, try closer.')
          A1 = 0.75
          A2 = 0.25
          GO TO 300
       endif
!                Contour has been located, insert new point in list
       do MOVE= NOWPTS,I1+1,-1
          XPTU(MOVE+1) = XPTU(MOVE)
          YPTU(MOVE+1) = YPTU(MOVE)
       end do
       NOWPTS = NOWPTS + 1
       XPTU(I1+1) = XMIDCR + XDIRCR*AOPT
       YPTU(I1+1) = YMIDCR + YDIRCR*AOPT
    end do
950 continue
!
    IERRF = NOWPTS
    CSTATU = 'SUCCESSFUL'
    if (NOWPTS.lt.NPTU)  CSTATU = 'INCOMPLETE'
!                make a lineprinter plot of the contour
    if (ISW(5).ge.0) then
       XPT(1) = U1MIN
       YPT(1) = U2MIN
       CHPT(1) = ' '
       NALL = min(NOWPTS+1,MAXCPT)
       do I= 2, NALL
          XPT(I) = XPTU(I-1)
          YPT(I) = YPTU(I-1)
          CHPT(I)= 'X'
       end do
       write (minuit_fileout,'(A,I3,2X,A)') ' Y-AXIS: PARAMETER ',KE2,&
            &        CPNAM(KE2)
       call MNPLOT(XPT,YPT,CHPT,NALL,minuit_fileout,NPAGWD,NPAGLN)
       write (minuit_fileout,'(25X,A,I3,2X,A)') 'X-AXIS: PARAMETER ',&
            &         KE1,CPNAM(KE1)
    endif
!                 print out the coordinates around the contour
    if (ISW(5).ge.1)  then
       NPCOL = (NOWPTS+1)/2
       NFCOL = NOWPTS/2
       write (minuit_fileout,'(/I5,A,G13.5,A,G11.3)') NOWPTS,&
            &    ' POINTS ON CONTOUR.   FMIN=',ABEST,'   ERRDEF=',UP
       write (minuit_fileout,'(9X,A,3X,A,18X,A,3X,A)')&
            &         CPNAM(KE1),CPNAM(KE2),CPNAM(KE1),CPNAM(KE2)
       do LINE = 1, NFCOL
          LR = LINE + NPCOL
          write (minuit_fileout,'(1X,I5,2G13.5,10X,I5,2G13.5)')&
               &     LINE,XPTU(LINE),YPTU(LINE),LR,XPTU(LR),YPTU(LR)
       end do
       if (NFCOL.lt.NPCOL) write (minuit_fileout,'(1X,I5,2G13.5)')&
            &                         NPCOL,XPTU(NPCOL),YPTU(NPCOL)
    endif
!                                    . . contour finished. reset v
    ITAUR = 1
    call MNFREE(1)
    call MNFREE(1)
    do J= 1, MPAR*(MPAR+1)/2
       VHMAT(J) = VTHMAT(J)
    end do
    do I= 1, MPAR
       GLOBCC(I) = GCC(I)
       WERR(I) = W(I)
       X(I) = XT(I)
    end do
    call MNINEX (X)
    EDM = SIGSAV
    AMIN = ABEST
    ISW(2) = ISW2
    ISW(4) = ISW4
    DCOVAR = DC
    ITAUR = 0
    NFCNMX = NFMXIN
    ISTRAT = ISTRAV
    U(KE1) = U1MIN
    U(KE2) = U2MIN
    GO TO 2000
!                                     Error returns
1350 write (minuit_fileout,'(A)') ' INVALID PARAMETER NUMBERS.'
    GO TO 1450
1400 write (minuit_fileout,'(A)') ' LESS THAN FOUR POINTS REQUESTED.'
1450 IERRF = -1
    CSTATU = 'USER ERROR'
    GO TO 2000
1500 write (minuit_fileout,'(A)') ' MNCONT UNABLE TO FIND FOUR POINTS.'
    U(KE1) = U1MIN
    U(KE2) = U2MIN
    IERRF = 0
    CSTATU = 'FAILED'
2000 continue
    CFROM = CHERE
    NFCNFR = NFCNCO
    return
  end subroutine minuit_mncont

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine MNHELP(COMD,LOUT)
    character*(*) COMD
    character*3 CMD3
    integer(FINT) :: lout
    if(COMD(1:1).eq.'*')then
       write(LOUT,10000)
       write(LOUT,10001)
       GO TO 99
    endif
10000 format('   ==>List of MINUIT Interactive commands:',/,&
     &' CLEar     Reset all parameter names and values undefined',/,&
     &' CONtour   Make contour map of the user function',/,&
     &' EXIT      Exit from Interactive Minuit',/,&
     &' FIX       Cause parameter(s) to remain constant',/,&
     &' HESse     Calculate the Hessian or error matrix.',/,&
     &' IMPROVE   Search for a new minimum around current minimum',/,&
     &' MIGrad    Minimize by the method of Migrad',/,&
     &' MINImize  MIGRAD + SIMPLEX method if Migrad fails',/,&
     &' MINOs     Exact (non-linear) parameter error analysis')
10001 format(' MNContour Calculate one MINOS function contour',/,&
     &' PARameter Define or redefine new parameters and values',/,&
     &' RELease   Make previously FIXed parameters variable again',/,&
     &' REStore   Release last parameter fixed',/,&
     &' SAVe      Save current parameter values on a file',/,&
     &' SCAn      Scan the user function by varying parameters',/,&
     &' SEEk      Minimize by the method of Monte Carlo',/,&
     &' SET       Set various MINUIT constants or conditions',/,&
     &' SHOw      Show values of current constants or conditions',/,&
     &' SIMplex   Minimize by the method of Simplex')
      CMD3=COMD(1:3)
      if(CMD3.eq.'CLE')then
         write(LOUT,10100)
         GO TO 99
      endif
10100 format(' ***>CLEAR',/,&
     &' Resets all parameter names and values to undefined.',/,&
     &' Must normally be followed by a PARameters command or ',/,&
     &' equivalent, in order to define parameter values.')
      if(CMD3.eq.'CON')then
         write(LOUT,10200)
         GO TO 99
      endif
10200 format(' ***>CONTOUR <par1>  <par2>  [devs]  [ngrid]',/,&
     &' Instructs Minuit to trace contour lines of the user function',/,&
     &' with respect to the two parameters whose external numbers',/,&
     &' are <par1> and <par2>.',/,&
     &' Other variable parameters of the function, if any, will have',/,&
     &' their values fixed at the current values during the contour',/,&
     &' tracing. The optional parameter [devs] (default value 2.)',/,&
     &' gives the number of standard deviations in each parameter',/,&
     &' which should lie entirely within the plotting area.',/,&
     &' Optional parameter [ngrid] (default value 25 unless page',/,&
     &' size is too small) determines the resolution of the plot,',/,&
     &' i.e. the number of rows and columns of the grid at which the',/,&
     &' function will be evaluated. [See also MNContour.]')
      if(CMD3.eq.'END')then
         write(LOUT,10300)
         GO TO 99
      endif
10300 format(' ***>END',/,&
     &' Signals the end of a data block (i.e., the end of a fit),',/,&
     &' and implies that execution should continue, because another',/,&
     &' Data Block follows. A Data Block is a set of Minuit data',/,&
     &' consisting of',/,&
     &'     (1) A Title,',/,&
     &'     (2) One or more Parameter Definitions,',/,&
     &'     (3) A blank line, and',/,&
     &'     (4) A set of Minuit Commands.',/,&
     &' The END command is used when more than one Data Block is to',/,&
     &' be used with the same FCN function. It first causes Minuit',/,&
     &' to issue a CALL FCN with IFLAG=3, in order to allow FCN to',/,&
     &' perform any calculations associated with the final fitted',/,&
     &' parameter values, unless a CALL FCN 3 command has already',/,&
     &' been executed at the current FCN value.')
      if(CMD3.eq.'EXI')then
         write(LOUT,10400)
         GO TO 99
      endif
10400 format(' ***>EXIT',/,&
     &' Signals the end of execution.',/,&
     &' The EXIT command first causes Minuit to issue a CALL FCN',/,&
     &' with IFLAG=3, to allow FCN to perform any calculations',/,&
     &' associated with the final fitted parameter values, unless a',/,&
     &' CALL FCN 3 command has already been executed.')
      if(CMD3.eq.'FIX')then
         write(LOUT,10500)
         GO TO 99
      endif
10500 format(' ***>FIX} <parno> [parno] ... [parno]',/,&
     &' Causes parameter(s) <parno> to be removed from the list of',/,&
     &' variable parameters, and their value(s) will remain constant',/,&
     &' during subsequent minimizations, etc., until another command',/,&
     &' changes their value(s) or status.')
      if(CMD3.eq.'HES')then
         write(LOUT,10600)
         GO TO 99
      endif
10600 format(' ***>HESse  [maxcalls]',/,&
     &' Calculate, by finite differences, the Hessian or error matrix.',&
     &/,'  That is, it calculates the full matrix of second derivatives'&
     &,/,' of the function with respect to the currently variable',/,&
     &' parameters, and inverts it, printing out the resulting error',/,&
     &' matrix. The optional argument [maxcalls] specifies the',/,&
     &' (approximate) maximum number of function calls after which',/,&
     &' the calculation will be stopped.')
      if(CMD3.eq.'IMP')then
         write(LOUT,10700)
         GO TO 99
      endif
10700 format(' ***>IMPROVE  [maxcalls]',/,&
     &' If a previous minimization has converged, and the current',/,&
     &' values of the parameters therefore correspond to a local',/,&
     &' minimum of the function, this command requests a search for',/,&
     &' additional distinct local minima.',/,&
     &' The optional argument [maxcalls] specifies the (approximate)',/,&
     &' maximum number of function calls after which the calculation',/,&
     &' will be stopped.')
      if(CMD3.eq.'MIG')then
         write(LOUT,10800)
         GO TO 99
      endif
10800 format(' ***>MIGrad  [maxcalls]  [tolerance]',/,&
     &' Causes minimization of the function by the method of Migrad,',/,&
     &' the most efficient and complete single method, recommended',/,&
     &' for general functions (see also MINImize).',/,&
     &' The minimization produces as a by-product the error matrix',/,&
     &' of the parameters, which is usually reliable unless warning',/,&
     &' messages are produced.',/,&
     &' The optional argument [maxcalls] specifies the (approximate)',/,&
     &' maximum number of function calls after which the calculation',/,&
     &' will be stopped even if it has not yet converged.',/,&
     &' The optional argument [tolerance] specifies required tolerance',&
     &/,' on the function value at the minimum.',/,&
     &' The default tolerance is 0.1, and the minimization will stop',/,&
     &' when the estimated vertical distance to the minimum (EDM) is',/,&
     &' less than 0.001*[tolerance]*UP (see [SET ERRordef]).')
      if(COMD(1:4).eq.'MINI')then
         write(LOUT,10900)
         GO TO 99
      endif
10900 format(' ***>MINImize  [maxcalls] [tolerance]',/,&
     &' Causes minimization of the function by the method of Migrad,',/,&
     &' as does the MIGrad command, but switches to the SIMplex method',&
     &/,' if Migrad fails to converge. Arguments are as for MIGrad.',/,&
     &' Note that command requires four characters to be unambiguous.')
      if(COMD(1:4).eq.'MINO')then
         write(LOUT,11000)
         GO TO 99
      endif
11000 format(' ***>MINOs  [maxcalls]  [parno] [parno] ...',/,&
     &' Causes a Minos error analysis to be performed on the parameters'&
     &,/,' whose numbers [parno] are specified. If none are specified,',&
     &/,' Minos errors are calculated for all variable parameters.',/,&
     &' Minos errors may be expensive to calculate, but are very',/,&
     &' reliable since they take account of non-linearities in the',/,&
     &' problem as well as parameter correlations, and are in general',/&
     &' asymmetric.',/,&
     &' The optional argument [maxcalls] specifies the (approximate)',/,&
     &' maximum number of function calls per parameter requested,',/,&
     &' after which the calculation will stop for that parameter.')
      if(CMD3.eq.'MNC')then
         write(LOUT,11100)
         GO TO 99
      endif
11100 format(' ***>MNContour  <par1> <par2> [npts]',/,&
     &' Calculates one function contour of FCN with respect to',/,&
     &' parameters par1 and par2, with FCN minimized always with',/,&
     &' respect to all other NPAR-2 variable parameters (if any).',/,&
     &' Minuit will try to find npts points on the contour (default 20)'&
     &,/,' If only two parameters are variable at the time, it is not',&
     &/,' necessary to specify their numbers. To calculate more than',/,&
     &' one contour, it is necessary to SET ERRordef to the appropriate'&
     &,/,' value and issue the MNContour command for each contour.')
      if(CMD3.eq.'PAR')then
         write(LOUT,11150)
         GO TO 99
      endif
11150 format(' ***>PARameters',/,&
     &' followed by one or more parameter definitions.',/,&
     &' Parameter definitions are of the form:',/,&
     &'   <number>  ''name''  <value>  <step>  [lolim] [uplim] ',/,&
     &' for example:',/,&
     &'  3  ''K width''  1.2   0.1' ,/,&
     &' the last definition is followed by a blank line or a zero.')
      if(CMD3.eq.'REL')then
         write(LOUT,11200)
         GO TO 99
      endif
11200 format(' ***>RELease  <parno> [parno] ... [parno]',/,&
     &' If <parno> is the number of a previously variable parameter',/,&
     &' which has been fixed by a command: FIX <parno>, then that',/,&
     &' parameter will return to variable status.  Otherwise a warning' &
     &,/,' message is printed and the command is ignored.',/,&
     &' Note that this command operates only on parameters which were',/&
     &' at one time variable and have been FIXed. It cannot make',/,&
     &' constant parameters variable; that must be done by redefining',/&
     &' the parameter with a PARameters command.')
      if(CMD3.eq.'RES')then
         write(LOUT,11300)
         GO TO 99
      endif
11300 format(' ***>REStore  [code]',/,&
     &' If no [code] is specified, this command restores all previously'&
     &,/,' FIXed parameters to variable status. If [code]=1, then only',&
     &/,' the last parameter FIXed is restored to variable status.',/,&
     &' If code is neither zero nor one, the command is ignored.')
      if(CMD3.eq.'RET')then
         write(LOUT,11400)
         GO TO 99
      endif
11400 format(' ***>RETURN',/,&
     &' Signals the end of a data block, and instructs Minuit to return'&
     &,/,' to the program which called it. The RETurn command first',/,&
     &' causes Minuit to CALL FCN with IFLAG=3, in order to allow FCN',/&
     &,' to perform any calculations associated with the final fitted',/&
     &,' parameter values, unless a CALL FCN 3 command has already been'&
     &,/,' executed at the current FCN value.')
      if(CMD3.eq.'SAV')then
         write(LOUT,11500)
         GO TO 99
      endif
11500 format(' ***>SAVe',/,&
     &' Causes the current parameter values to be saved on a file in',/,&
     &' such a format that they can be read in again as Minuit',/,&
     &' parameter definitions. If the covariance matrix exists, it is',/&
     &,' also output in such a format. The unit number is by default 7,'&
     &,/,' or that specified by the user in his call to MINTIO or',/,&
     &' MNINIT. The user is responsible for opening the file previous'  &
     &,/,' to issuing the [SAVe] command (except where this can be done'&
     &,/,' interactively).')
      if(CMD3.eq.'SCA')then
         write(LOUT,11600)
         GO TO 99
      endif
11600 format(' ***>SCAn  [parno]  [numpts] [from]  [to]',/,&
     &' Scans the value of the user function by varying parameter',/,&
     &' number [parno], leaving all other parameters fixed at the',/,&
     &' current value. If [parno] is not specified, all variable',/,&
     &' parameters are scanned in sequence.',/,&
     &' The number of points [numpts] in the scan is 40 by default,',/,&
     &' and cannot exceed 100. The range of the scan is by default',/,&
     &' 2 standard deviations on each side of the current best value,',&
     &/,' but can be specified as from [from] to [to].',/,&
     &' After each scan, if a new minimum is found, the best parameter' &
     &,/,' values are retained as start values for future scans or',/,&
     &' minimizations. The curve resulting from each scan is plotted',/&
     &,' on the output unit in order to show the approximate behaviour' &
     &,/,' of the function.',/,&
     &' This command is not intended for minimization, but is sometimes'&
     &,/,' useful for debugging the user function or finding a',/,&
     &' reasonable starting point.')
      if(CMD3.eq.'SEE')then
         write(LOUT,11700)
         GO TO 99
      endif
11700 format(' ***>SEEk  [maxcalls]  [devs]',/,&
     &' Causes a Monte Carlo minimization of the function, by choosing',&
     &/,' random values of the variable parameters, chosen uniformly',/,&
     &' over a hypercube centered at the current best value.',/,&
     &' The region size is by default 3 standard deviations on each',/,&
     &' side, but can be changed by specifying the value of [devs].')
      if(CMD3.eq.'SET')then
         write(LOUT,11800)
         write(LOUT,11801)
         write(LOUT,11802)
         write(LOUT,11803)
         write(LOUT,11804)
         write(LOUT,11805)
         write(LOUT,11806)
         write(LOUT,11807)
         write(LOUT,11808)
         write(LOUT,11809)
         write(LOUT,11810)
         write(LOUT,11811)
         write(LOUT,11812)
         write(LOUT,11813)
         write(LOUT,11814)
         write(LOUT,11815)
         write(LOUT,11816)
         write(LOUT,11817)
         GO TO 99
      endif
11800 format(' ***>SET <option_name>',/,/,&
     &'  SET BATch',/,&
     &'    Informs Minuit that it is running in batch mode.',//,&
     &'  SET EPSmachine  <accuracy>',/,&
     &'    Informs Minuit that the relative floating point arithmetic',/&
     &'    precision is <accuracy>. Minuit determines the nominal',/,&
     &'    precision itself, but the SET EPSmachine command can be',/,&
     &'    used to override Minuit own determination, when the user',/,&
     &'    knows that the FCN function value is not calculated to',/,&
     &'    the nominal machine accuracy. Typical values of <accuracy>',/&
     &'    are between 10**-5 and 10**-14.')


11801 format(/,'  SET ERRordef  <up>',/,&
     &'    Sets the value of UP (default value= 1.), defining parameter'&
     &,/,'    errors. Minuit defines parameter errors as the change',/,&
     &'    in parameter value required to change the function value',/,&
     &'    by UP. Normally, for chisquared fits UP=1, and for negative' &
     &,/,'    log likelihood, UP=0.5.')

11802 format(/,'   SET GRAdient  [force]',/,&
     &'    Informs Minuit that the user function is prepared to',/,&
     &'    calculate its own first derivatives and return their values' &
     &,/,'    in the array GRAD when IFLAG=2 (see specs of FCN).',/,&
     &'    If [force] is not specified, Minuit will calculate',/,&
     &'    the FCN derivatives by finite differences at the current',/,&
     &'    point and compare with the user calculation at that point,'  &
     &,/,'    accepting the user values only if they agree.',/,&
     &'    If [force]=1, Minuit does not do its own derivative',/,&
     &'    calculation, and uses the derivatives calculated in FCN.')

11803 format(/,'   SET INPut  [unitno]  [filename]',/,&
     &'    Causes Minuit, in data-driven mode only, to read subsequent',&
     &/,'    commands (or parameter definitions) from a different input'&
     &,/,'    file. If no [unitno] is specified, reading reverts to the'&
     &,/,'    previous input file, assuming that there was one.',/,&
     &'    If [unitno] is specified, and that unit has not been opened,'&
     &,/,'    then Minuit attempts to open the file [filename]} if a',/,&
     &'    name is specified. If running in interactive mode and',/,&
     &'    [filename] is not specified and [unitno] is not opened,',/,&
     &'    Minuit prompts the user to enter a file name.',/,&
     &'    If the word REWIND is added to the command (note:no blanks',/&
     &'    between INPUT and REWIND), the file is rewound before',/,&
     &'    reading. Note that this command is implemented in standard',/&
     &'    Fortran 77 and the results may depend on the  system;',/,&
     &'    for example, if a filename is given under VM/CMS, it must',/,&
     &'    be preceeded by a slash.')

11804 format(/,'   SET INTeractive',/,&
     &'    Informs Minuit that it is running interactively.')

11805 format(/,'   SET LIMits  [parno]  [lolim]  [uplim]',/,&
     &'    Allows the user to change the limits on one or all',/,&
     &'    parameters. If no arguments are specified, all limits are',/,&
     &'    removed from all parameters. If [parno] alone is specified,',&
     &/,'    limits are removed from parameter [parno].',/,&
     &'    If all arguments are specified, then parameter [parno] will',&
     &/,'    be bounded between [lolim] and [uplim].',/,&
     &'    Limits can be specified in either order, Minuit will take',/,&
     &'    the smaller as [lolim] and the larger as [uplim].',/,&
     &'    However, if [lolim] is equal to [uplim], an error condition',&
     &/,'    results.')

11806 format(/,'   SET LINesperpage',/,&
     &'     Sets the number of lines for one page of output.',/,&
     &'     Default value is 24 for interactive mode')

11807 format(/,'   SET NOGradient',/,&
     &'    The inverse of SET GRAdient, instructs Minuit not to',&
     &/,'    use the first derivatives calculated by the user in FCN.')

11808 format(/,'   SET NOWarnings',/,&
     &'    Supresses Minuit warning messages.')

11809 format(/,'   SET OUTputfile  <unitno>',/,&
     &'    Instructs Minuit to write further output to unit <unitno>.')

11810 format(/,'   SET PAGethrow  <integer>',/,&
     &'    Sets the carriage control character for ``new page'' to',/,&
     &'    <integer>. Thus the value 1 produces a new page, and 0',/,&
     &'    produces a blank line, on some devices (see TOPofpage)')


11811 format(/,'   SET PARameter  <parno>  <value>',/,&
     &'    Sets the value of parameter <parno> to <value>.',/,&
     &'    The parameter in question may be variable, fixed, or',/,&
     &'    constant, but must be defined.')

11812 format(/,'   SET PRIntout  <level>',/,&
     &'    Sets the print level, determining how much output will be',/,&
     &'    produced. Allowed values and their meanings are displayed',/,&
     &'    after a SHOw PRInt command, and are currently <level>=:',/,&
     &'      [-1]  no output except from SHOW commands',/,&
     &'       [0]  minimum output',/,&
     &'       [1]  default value, normal output',/,&
     &'       [2]  additional output giving intermediate results.',/,&
     &'       [3]  maximum output, showing progress of minimizations.',/&
     &'    Note: See also the SET WARnings command.')

11813 format(/,'   SET RANdomgenerator  <seed>',/,&
     &'    Sets the seed of the random number generator used in SEEk.',/&
     &'    This can be any integer between 10000 and 900000000, for',/,&
     &'    example one which was output from a SHOw RANdom command of',/&
     &'    a previous run.')

11814 format(/,'   SET STRategy  <level>',/,&
     &'    Sets the strategy to be used in calculating first and second'&
     &,/,'    derivatives and in certain minimization methods.',/,&
     &'    In general, low values of <level> mean fewer function calls',&
     &/,'    and high values mean more reliable minimization.',/,&
     &'    Currently allowed values are 0, 1 (default), and 2.')

11815 format(/,'   SET TITle',/,&
     &'    Informs Minuit that the next input line is to be considered',&
     &/,'    the (new) title for this task or sub-task.  This is for',/,&
     &'    the convenience of the user in reading his output.')

11816 format(/,'   SET WARnings',/,&
     &'    Instructs Minuit to output warning messages when suspicious',&
     &/,'    conditions arise which may indicate unreliable results.',/&
     &'    This is the default.')

11817 format(/,'    SET WIDthpage',/,&
     &'    Informs Minuit of the output page width.',/,&
     &'    Default values are 80 for interactive jobs')
      if(CMD3.eq.'SHO')then
         write(LOUT,11900)
         write(LOUT,11901)
         write(LOUT,11902)
         write(LOUT,11903)
         write(LOUT,11904)
         GO TO 99
      endif
11900 format(' ***>SHOw  <option_name>',/,&
     &'  All SET XXXX commands have a corresponding SHOw XXXX command.',&
     &/,'  In addition, the SHOw commands listed starting here have no',&
     &/,'  corresponding SET command for obvious reasons.')

11901 format(/,'   SHOw CORrelations',/,&
     &'    Calculates and prints the parameter correlations from the',/,&
     &'    error matrix.')

11902 format(/,'   SHOw COVariance',/,&
     &'    Prints the (external) covariance (error) matrix.')

11903 format(/,'   SHOw EIGenvalues',/,&
     &'    Calculates and prints the eigenvalues of the covariance',/,&
     &'    matrix.')

11904 format(/,'   SHOw FCNvalue',/,&
     &'    Prints the current value of FCN.')
      if(CMD3.eq.'SIM')then
         write(LOUT,12000)
         GO TO 99
      endif
12000 format(' ***>SIMplex  [maxcalls]  [tolerance]',/,&
     &' Performs a function minimization using the simplex method of',/&
     &' Nelder and Mead. Minimization terminates either when the',/,&
     &' function has been called (approximately) [maxcalls] times,',/,&
     &' or when the estimated vertical distance to minimum (EDM) is',/,&
     &' less than [tolerance].',/,&
     &' The default value of [tolerance] is 0.1*UP(see SET ERRordef).')
      if(CMD3.eq.'STA')then
         write(LOUT,12100)
         GO TO 99
      endif
12100 format(' ***>STAndard',/,&
     &' Causes Minuit to execute the Fortran instruction CALL STAND',/,&
     &' where STAND is a subroutine supplied by the user.')
      if(CMD3.eq.'STO')then
         write(LOUT,12200)
         GO TO 99
      endif
12200 format(' ***>STOP',/,&
     &' Same as EXIT.')
      if(CMD3.eq.'TOP')then
         write(LOUT,12300)
         GO TO 99
      endif
12300 format(' ***>TOPofpage',/,&
     &' Causes Minuit to write the character specified in a',/,&
     &' SET PAGethrow command (default = 1) to column 1 of the output'  &
     &,/,' file, which may or may not position your output medium to',&
     &/,' the top of a page depending on the device and system.')
      write(LOUT,13000)
13000 format(' Unknown MINUIT command. Type HELP for list of commands.')
   99 return
  end subroutine MNHELP

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine stand
!c        optional user-supplied subroutine is called whenever the
!c        command "standard" appears.
!c
    write (*,*) 'STANDARD NOT IMPLEMENTED.'
    return
  end subroutine stand

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine minuit_minuit
!
!  CPNAM   Parameter name (10 characters)
!  U       External (visible to user in minuit_fcn) value of parameter
!  ALIM, BLIM Lower and upper parameter limits. If both zero, no limits.
!  ERP,ERN Positive and negative MINOS errors, if calculated.
!  WERR    External parameter error (standard deviation, defined by UP)
!  GLOBCC  Global Correlation Coefficient
!  NVARL   =-1 if parameter undefined,      =0 if constant,
!          = 1 if variable without limits,  =4 if variable with limits
!   (Note that if parameter has been fixed, NVARL=1 or =4, and NIOFEX=0)
!  NIOFEX  Internal parameter number, or zero if not currently variable
!  NEXOFI  External parameter number for currently variable parameters
!  X, XT   Internal parameter values (X are sometimes saved in XT)
!  DIRIN   (Internal) step sizes for current step
!  variables with names ending in ..S are saved values for fixed params
!  VHMAT   (Internal) error matrix stored as Half MATrix, since
!                it is symmetric
!  VTHMAT  VHMAT is sometimes saved in VTHMAT, especially in MNMNOT
!
!  ISW definitions:
!      ISW(1) =0 normally, =1 means CALL LIMIT EXCEEDED
!      ISW(2) =0 means no error matrix
!             =1 means only approximate error matrix
!             =2 means full error matrix, but forced pos-def.
!             =3 means good normal full error matrix exists
!      ISW(3) =0 if Minuit is calculating the first derivatives
!             =1 if first derivatives calculated inside minuit_fcn
!      ISW(4) =-1 if most recent minimization did not converge.
!             = 0 if problem redefined since most recent minimization.
!             =+1 if most recent minimization did converge.
!      ISW(5) is the PRInt level.  See SHO PRIntlevel
!      ISW(6) = 0 for batch mode, =1 for interactive mode
!                      =-1 for originally interactive temporarily batch
!
!  LWARN is true if warning messges are to be put out (default=true)
!            SET WARN turns it on, set NOWarn turns it off
!  LREPOR is true if exceptional conditions are put out (default=false)
!            SET DEBUG turns it on, SET NODebug turns it off
!  LIMSET is true if a parameter is up against limits (for MINOS)
!  LNOLIM is true if there are no limits on any parameters (not yet used
!  LNEWMN is true if the previous process has unexpectedly improved minuit_fcn
!  LPHEAD is true if a heading should be put out for the next parameter
!        definition, false if a parameter has just been defined
!
    !
    integer(FINT) :: irc, iflgut, nparx
    real(FREAL) :: fzero, first, fnew
    !
    character(40) :: CWHYXT
    data CWHYXT/'FOR UNKNOWN REASONS                     '/
    integer(FINT) :: JSYSRD,JSYSWR,JSYSSA
    data JSYSRD,JSYSWR,JSYSSA/5,6,7/
!                                 . . . . . . . . . . initialize minuit
    if(.not.minuit_isinit()) then
       ! Set defaults
       irc=minuit_init(JSYSRD,JSYSWR,JSYSSA)
    end if
    write (minuit_fileout,'(1X,75(1H*))')
!                                      . . . . initialize new data block
100 continue
    write (minuit_fileout,'(1X,75(1H*))')
    NBLOCK = NBLOCK + 1
    write (minuit_fileout,'(26X,A,I4)')  'MINUIT DATA BLOCK NO.',NBLOCK
    write (minuit_fileout,'(1X,75(1H*))')
!               . . . . . . . . . . .   set parameter lists to undefined
    call MNCLER
!                                             . . . . . . . . read title
    call MNREAD(1,IFLGUT)
    if (IFLGUT .eq. 2)  GO TO 500
    if (IFLGUT .eq. 3)  GO TO 600
!                                        . . . . . . . . read parameters
    call MNREAD(2,IFLGUT)
    if (IFLGUT .eq. 2)  GO TO 500
    if (IFLGUT .eq. 3)  GO TO 600
    if (IFLGUT .eq. 4)  GO TO 700
!                              . . . . . . verify FCN not time-dependent
    write (minuit_fileout,'(/A,A)') ' MINUIT: FIRST CALL TO USER FUNCTION,',  &
         &    ' WITH IFLAG=1'
    NPARX = NPAR
    call MNINEX(X)
    FZERO = UNDEFI
    call minuit_fcn(NPARX,GIN,FZERO,U,1)
    FIRST = UNDEFI
    call minuit_fcn(NPARX,GIN,FIRST,U,4)
    NFCN = 2
    if (FZERO.eq.UNDEFI .and. FIRST.eq.UNDEFI)  then
       CWHYXT = 'BY ERROR IN USER FUNCTION.  '
       write (minuit_fileout,'(/A,A/)') ' USER HAS NOT CALCULATED FUNCTION', &
            &    ' VALUE WHEN IFLAG=1 OR 4'
       GO TO 800
    endif
    AMIN = FIRST
    if (FIRST .eq. UNDEFI) AMIN=FZERO
    call MNPRIN(1,AMIN)
    NFCN = 2
    if (FIRST .eq. FZERO)  GO TO 300
    FNEW = 0.0
    call minuit_fcn(NPARX,GIN,FNEW,U,4)
    if  (FNEW .ne. AMIN) write (minuit_fileout,280) AMIN, FNEW
280 format (/' MINUIT WARNING: PROBABLE ERROR IN USER FUNCTION.'/     &
         &         ' FOR FIXED VALUES OF PARAMETERS, FCN IS TIME-DEPENDENT'/&
         &         ' F =',E22.14,' FOR FIRST CALL'/                         &
         &         ' F =',E22.14,' FOR SECOND CALL.'/)
    NFCN = 3
300 FVAL3 = 2.0*AMIN+1.0
    !                                   . . . . . . . . . . . read commands
    call MNREAD(3,IFLGUT)
    if (IFLGUT .eq. 2)  GO TO 500
    if (IFLGUT .eq. 3)  GO TO 600
    if (IFLGUT .eq. 4)  GO TO 700
    !Warning! some implicit dimension assumed!
    CWHYXT = 'BY MINUIT COMMAND: '//CWORD(1:21)
    if (index(CWORD,'STOP').gt. 0)  GO TO 800
    if (index(CWORD,'EXI') .gt. 0)  GO TO 800
    if (index(CWORD,'RET') .eq. 0)  GO TO 100
    CWHYXT = 'AND RETURNS TO USER PROGRAM.    '
    write (minuit_fileout,'(A,A)')  ' ..........MINUIT TERMINATED ',CWHYXT
    return
    !                                           . . . . . . stop conditions
500 continue
    CWHYXT = 'BY END-OF-DATA ON PRIMARY INPUT FILE.   '
    GO TO 800
600 continue
    CWHYXT = 'BY UNRECOVERABLE READ ERROR ON INPUT.   '
    GO TO 800
700 continue
    CWHYXT = ': FATAL ERROR IN PARAMETER DEFINITIONS. '
800 write (minuit_fileout,'(A,A)')  ' ..........MINUIT TERMINATED ',CWHYXT
    stop
!
!  ......................entry to set unit numbers  - - - - - - - - - -
  end subroutine minuit_minuit

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine mnit_mintio(ia,ib,ic)
    !
    integer(FINT) :: ia, ib ,ic
    !
    minuit_filein = ia
    minuit_fileout = ib
    minuit_filesave = ic
    !
  end subroutine mnit_mintio
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine MNREAD(IFLGIN,IFLGUT)
!C        Called from MINUIT.  Reads all user input to MINUIT.
!C     This routine is highly unstructured and defies normal logic.
!C
!C     IFLGIN indicates the function originally requested:
!C           = 1: read one-line title
!C             2: read parameter definitions
!C             3: read MINUIT commands
!C
!C     IFLGUT= 1: reading terminated normally
!C             2: end-of-data on input
!C             3: unrecoverable read error
!C             4: unable to process parameter requests
!C             5: more than 100 incomprehensible commands
!C internally,
!C     IFLGDO indicates the subfunction to be performed on the next
!C         input record: 1: read a one-line title
!C                       2: read a parameter definition
!C                       3: read a command
!C                       4: read in covariance matrix
!C     for example, when IFLGIN=3, but IFLGDO=1, then it should read
!C       a title, but this was requested by a command, not by MINUIT.
!C
    !
    !
    integer(FINT) :: i
    integer(FINT) :: iflgin, iflgut, iflgdo, incomp
    integer(FINT) :: ic, ierr, icondp, icondn, npar2
    !
    !
    character(80) :: CRDBUF !*80
    character(10) :: CUPBUF !*10
    character(26) :: CLOWER, CUPPER !*26
    character(40), dimension(3) :: CPROMT !(3)*40
    logical LEOF
    data CPROMT/' ENTER MINUIT TITLE, or "SET INPUT n" : ',           &
         &      ' ENTER MINUIT PARAMETER DEFINITION:     ',           &
         &      ' ENTER MINUIT COMMAND:                  '/
!
    data CLOWER/'abcdefghijklmnopqrstuvwxyz'/
    data CUPPER/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
!
    IFLGUT = 1
    IFLGDO = IFLGIN
    LEOF = .false.
    INCOMP = 0
!                                           . . . . read next record
10  continue
    if (ISW(6) .eq. 1) then
       write (minuit_fileout,'(A)') CPROMT(IFLGDO)
       if (IFLGDO .eq. 2)  LPHEAD = .false.
    endif
    CRDBUF = '   '
    read (minuit_filein,'(A)',ERR=500,end=45)  CRDBUF
!
!                 CUPBUF is the first few characters in upper case
    CUPBUF(1:10) = CRDBUF(1:10)
    do I= 1, 10
       if (CRDBUF(I:I) .eq. '''') GO TO 13
       do IC= 1, 26
          if (CRDBUF(I:I) .eq. CLOWER(IC:IC)) CUPBUF(I:I)=CUPPER(IC:IC)
       end do
    end do
13  continue
!                                           . .   preemptive commands
    LEOF = .false.
    if (index(CUPBUF,'*EOF') .eq. 1)    then
       write (minuit_fileout,'(A,I3)') ' *EOF ENCOUNTERED ON UNIT NO.',minuit_filein
       LPHEAD = .true.
       GO TO 50
    endif
    if (index(CUPBUF,'SET INP') .eq. 1)    then
       ICOMND = ICOMND + 1
       write (minuit_fileout, 21) ICOMND,CRDBUF(1:50)
21     format (' **********'/' **',I5,' **',A/' **********')
       LPHEAD = .true.
       GO TO 50
    endif
    GO TO 80
!                                    . . hardware EOF on current minuit_filein
45  CRDBUF = '*EOF '
    write (minuit_fileout,'(A,I3)') ' END OF DATA ON UNIT NO.',minuit_filein
!                                     or SET INPUT command
50  continue
    call MNSTIN(CRDBUF,IERR)
    if (IERR .eq. 0)  GO TO 10
    if (IERR .eq. 2)  then
       if (.not. LEOF) then
          write (minuit_fileout,'(A,A/)') ' TWO CONSECUTIVE EOFs ON ',     &
               &              'PRIMARY INPUT FILE WILL TERMINATE EXECUTION.'
          LEOF = .true.
          GO TO 10
       endif
    endif
    IFLGUT = IERR
    GO TO 900
80  if (IFLGDO .gt. 1) GO TO 100
!                            read title        . . . . .   IFLGDO = 1
!              if title is 'SET TITLE', skip and read again
    if (index(CUPBUF,'SET TIT') .eq. 1)  GO TO 10
    call minuit_mnseti(CRDBUF(1:50))
    write (minuit_fileout,'(1X,A50)')  CTITL
    write (minuit_fileout,'(1X,78(1H*))')
    LPHEAD = .true.
    if (IFLGIN .eq. IFLGDO)  GO TO 900
    IFLGDO = IFLGIN
    GO TO 10
!                            data record is not a title.
100 continue
    if (IFLGDO .gt. 2)  GO TO 300
!                          expect parameter definitions.   IFLGDO = 2
!              if parameter def is 'PARAMETER', skip and read again
    if (index(CUPBUF,'PAR') .eq. 1)  GO TO 10
!              if line starts with SET TITLE, read a title first
    if (index(CUPBUF,'SET TIT') .eq. 1)  then
       IFLGDO = 1
       GO TO 10
    endif
!                      we really have parameter definitions now
    call MNPARS(CRDBUF,ICONDP)
    if (ICONDP .eq. 0)  GO TO 10
!          format error
    if (ICONDP .eq. 1)  then
       if (ISW(6) .eq. 1)  then
          write (minuit_fileout,'(A)') ' FORMAT ERROR.  IGNORED.  ENTER AGAIN.'
          GO TO 10
       else
          write (minuit_fileout,'(A)') ' ERROR IN PARAMETER DEFINITION'
          IFLGUT = 4
          GO TO 900
       endif
    endif
!                     ICONDP = 2            . . . end parameter requests
    if (ISW(5).ge.0 .and. ISW(6).lt.1) write (minuit_fileout,'(4X,75(1H*))')
    LPHEAD = .true.
    if (IFLGIN .eq. IFLGDO)  GO TO 900
    IFLGDO = IFLGIN
    GO TO 10
!                                              . . . . .   IFLGDO = 3
!                                           read commands
300 continue
    call MNCOMD(CRDBUF,ICONDN)
!C     ICONDN = 0: command executed normally
!C              1: command is blank, ignored
!C              2: command line unreadable, ignored
!C              3: unknown command, ignored
!C              4: abnormal termination (e.g., MIGRAD not converged)
!C              5: command is a request to read PARAMETER definitions
!C              6: 'SET INPUT' command
!C              7: 'SET TITLE' command
!C              8: 'SET COVAR' command
!C              9: reserved
!C             10: END command
!C             11: EXIT or STOP command
!C             12: RETURN command
    if (ICONDN .eq. 2 .or. ICONDN .eq. 3) then
       INCOMP = INCOMP + 1
       if (INCOMP .gt. 100) then
          IFLGUT = 5
          GO TO 900
       endif
    endif
!                         parameter
    if (ICONDN .eq. 5)  IFLGDO = 2
!                         SET INPUT
    if (ICONDN .eq. 6)  GO TO 50
!                         SET TITLE
    if (ICONDN .eq. 7)  IFLGDO = 1
!                                        . . . . . . . . . . set covar
    if (ICONDN .eq. 8) then
       ICOMND = ICOMND + 1
       write (minuit_fileout,405) ICOMND,CRDBUF(1:50)
405    format (1H ,10(1H*)/' **',I5,' **',A)
       write (minuit_fileout, '(1H ,10(1H*))' )
       NPAR2 = NPAR*(NPAR+1)/2
       read (minuit_filein,420,ERR=500,end=45)  (VHMAT(I),I=1,NPAR2)
420    format (BN,7E11.4,3X)
       ISW(2) = 3
       DCOVAR = 0.0
       if (ISW(5) .ge. 0)  call MNMATU(1)
       if (ISW(5) .ge. 1)  call MNPRIN(2,AMIN)
       GO TO 10
    endif
    if (ICONDN .lt. 10) GO TO 10
    GO TO 900
!                                              . . . . error conditions
500 IFLGUT = 3
900 return
  end subroutine MNREAD
!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine MNSTIN(CRDBUF,IERR)
!C Called from MNREAD.
!C Implements the SET INPUT command to change input units.
!C If command is: 'SET INPUT'   'SET INPUT 0'   or  '*EOF',
!C                 or 'SET INPUT , ,  ',
!C                reverts to previous input unit number,if any.
!C
!C      If it is: 'SET INPUT n'  or  'SET INPUT n filename',
!C                changes to new input file, added to stack
!C
!C      IERR = 0: reading terminated normally
!C             2: end-of-data on primary input file
!C             3: unrecoverable read error
!C             4: unable to process request
!C
    !
    !
    integer(FINT) :: ierr, lend, ic, icol, ic1, ic2, iunit
    real(FREAL) :: funit
    !
    !
    character(*) :: CRDBUF !*(*)
    character(10) :: CUNIT !*10
    character(64) :: CFNAME, CGNAME !*64
    character(1) :: CANSWR !*1
    character(16) :: CMODE !*16
    logical :: LOPEN, LREWIN, NONAME, LNAME
    NONAME = .true.
    IERR = 0
    if (index(CRDBUF,'*EOF') .eq. 1) goto 190
    if (index(CRDBUF,'*eof') .eq. 1) goto 190
    LEND = len(CRDBUF)
!                               look for end of SET INPUT command
    do IC= 8,LEND
       if (CRDBUF(IC:IC) .eq. ' ') goto 25
       if (CRDBUF(IC:IC) .eq. ',') goto 53
    end do
    goto 200
25  continue
!         look for end of separator between command and first argument
    ICOL = IC+1
    do IC= ICOL,LEND
       if (CRDBUF(IC:IC) .eq. ' ') cycle
       if (CRDBUF(IC:IC) .eq. ',') goto 53
       goto 55
    end do
    goto 200
53  IC = IC + 1
55  IC1 = IC
!                      see if "REWIND" was requested in command
    LREWIN = .false.
    if (index(CRDBUF(1:IC1),'REW') .gt. 5)  LREWIN=.true.
    if (index(CRDBUF(1:IC1),'rew') .gt. 5)  LREWIN=.true.
!                      first argument begins in or after col IC1
    do IC= IC1,LEND
       if (CRDBUF(IC:IC) .eq. ' ') cycle
       if (CRDBUF(IC:IC) .eq. ',') goto 200
       goto 80
    end do
    goto 200
80  IC1 = IC
!                        first argument really begins in col IC1
    do IC= IC1+1,LEND
       if (CRDBUF(IC:IC) .eq. ' ') goto 108
       if (CRDBUF(IC:IC) .eq. ',') goto 108
    end do
    IC = LEND + 1
108 IC2 = IC-1
!                            end of first argument is in col IC2
110 continue
    CUNIT = CRDBUF(IC1:IC2)
    write (minuit_fileout,'(A,A)') ' UNIT NO. :',CUNIT
    read (CUNIT,'(BN,F10.0)',ERR=500) funit
    IUNIT = int(funit,FINT)
    if (IUNIT .eq. 0)  goto 200
!                             skip blanks and commas, find file name
    do IC= IC2+1,LEND
       if (CRDBUF(IC:IC) .eq. ' ') cycle
       if (CRDBUF(IC:IC) .eq. ',') cycle
       goto 130
    end do
    goto 131
130 continue
    CFNAME = CRDBUF(IC:LEND)
    NONAME = .false.
    write (minuit_fileout, '(A,A)') ' FILE NAME IS:',CFNAME
!              ask if file exists, if not ask for name and open it
131 continue
    inquire(UNIT=IUNIT,OPENED=LOPEN,NAMED=LNAME,NAME=CGNAME)
    if (LOPEN) then
       if (NONAME) then
          goto 136
       else
          if (.not.LNAME) CGNAME='unknown'
          write (minuit_fileout,132) IUNIT,CGNAME,CFNAME
132       format (' UNIT',I3,' ALREADY OPENED WITH NAME:',A/         &
               &                  '                 NEW NAME IGNORED:',A)
       endif
    else
!                new file, open it
       write (minuit_fileout,135) IUNIT
135    format (' UNIT',I3,' IS NOT OPENED.')
       if (NONAME) then
          write (minuit_fileout,'(A)') ' NO FILE NAME GIVEN IN COMMAND.'
          if (ISW(6) .lt. 1)  goto 800
          write (minuit_fileout,'(A)') ' PLEASE GIVE FILE NAME:'
          read (minuit_filein,'(A)') CFNAME
       endif
       open (UNIT=IUNIT,FILE=CFNAME,STATUS='OLD',ERR=600)
       write (minuit_fileout,'(A)') ' FILE OPENED SUCCESSFULLY.'
    endif
!                                     . .   file is correctly opened
136 if (LREWIN) goto 150
    if (ISW(6) .lt. 1)  goto 300
    write (minuit_fileout,137)  IUNIT
137 format (' SHOULD UNIT',I3,' BE REWOUND?' )
    read  (minuit_filein,'(A)')  CANSWR
    if (CANSWR.ne.'Y' .and. CANSWR.ne.'y') goto 300
150 rewind IUNIT
    goto 300
!                      *EOF
190 continue
    if (NSTKRD .eq. 0)  then
       IERR = 2
       goto 900
    endif
!                      revert to previous input file
200 continue
    if (NSTKRD .eq. 0)  then
       write (minuit_fileout, '(A,A)') ' COMMAND IGNORED:',CRDBUF
       write (minuit_fileout, '(A)') ' ALREADY READING FROM PRIMARY INPUT'
    else
       minuit_filein = ISTKRD(NSTKRD)
       NSTKRD = NSTKRD - 1
       if (NSTKRD .eq. 0)  ISW(6) = IABS(ISW(6))
       if (ISW(5) .ge. 0)  then
          inquire(UNIT=minuit_filein,NAMED=LNAME,NAME=CFNAME)
          CMODE = 'BATCH MODE      '
          if (ISW(6) .eq. 1)  CMODE = 'INTERACTIVE MODE'
          if (.not.LNAME) CFNAME='unknown'
          if (MNUNPT(CFNAME))  CFNAME='unprintable'
          write (minuit_fileout,290) CMODE,minuit_filein,CFNAME
290       format (' INPUT WILL NOW BE READ IN ',A,' FROM UNIT NO.',I3/  &
               &    ' FILENAME: ',A)
       endif
    endif
    goto 900
!                      switch to new input file, add to stack
300 continue
    if (NSTKRD .ge. MAXSTK)  then
       write (minuit_fileout, '(A)') ' INPUT FILE STACK SIZE EXCEEDED.'
       goto 800
    endif
    NSTKRD = NSTKRD + 1
    ISTKRD(NSTKRD) = minuit_filein
    minuit_filein = IUNIT
!                   ISW(6) = 0 for batch, =1 for interactive, and
!                      =-1 for originally interactive temporarily batch
    if (ISW(6) .eq. 1)  ISW(6) = -1
    goto 900
!                      format error
500 continue
    write (minuit_fileout,'(A,A)') ' CANNOT READ FOLLOWING AS INTEGER:',CUNIT
    goto 800
600 continue
    write (minuit_fileout, 601) CFNAME
601 format (' SYSTEM IS UNABLE TO OPEN FILE:',A)
!                      serious error
800 continue
    IERR = 3
900 continue
    return
  end subroutine MNSTIN

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine MNPARS(CRDBUF,ICONDN)
!C        Called from MNREAD and user-callable
!C    Implements one parameter definition, that is:
!C       parses the string CRDBUF and calls MNPARM
!
! output conditions:
!        ICONDN = 0    all OK
!        ICONDN = 1    error, attempt to define parameter is ignored
!        ICONDN = 2    end of parameter definitions
!
!
    !
    !
    integer(FINT) :: icondn, lenbuf, kapo1, kapo2, istart, k, icy, ibegin, lnc, llist, ierr
    real(FREAL) ::fk, uk, wk, a, b, xk
    !
    !
    real(FREAL), dimension(MAXP) :: PLIST !(MAXP)
    character(10) :: CNAMK !*10
    character(*) :: CRDBUF !*(*)
    character(20) :: CELMNT !*20
    character(MAXCWD) :: COMAND !*(MAXCWD)
!
    LENBUF = len(CRDBUF)
!                     find out whether fixed or free-field format
    KAPO1 = index(CRDBUF,'''')
    if (KAPO1 .eq. 0)  goto 150
    KAPO2 = index(CRDBUF(KAPO1+1:),'''')
    if (KAPO2 .eq. 0)  goto 150
!          new (free-field) format
    KAPO2 = KAPO2 + KAPO1
!                             skip leading blanks if any
    do ISTART=1, KAPO1-1
       if (CRDBUF(ISTART:ISTART) .ne. ' ')  goto 120
    end do
    goto 210
120 continue
!                               parameter number integer
    CELMNT = CRDBUF(ISTART:KAPO1-1)
    read (CELMNT,'(BN,F20.0)',ERR=180) FK
    K = FK
    if (K .le. 0)  goto 210
    !Warning: some assumed dimension here!
    !CNAMK = 'PARAM '//CELMNT
    CNAMK = 'PARAM '//CELMNT(1:4)
    if (KAPO2-KAPO1 .gt. 1) CNAMK = CRDBUF(KAPO1+1:KAPO2-1)
!  special handling if comma or blanks and a comma follow 'name'
    do ICY= KAPO2+1,LENBUF
       if (CRDBUF(ICY:ICY) .eq. ',') goto 139
       if (CRDBUF(ICY:ICY) .ne. ' ') goto 140
    end do
    UK = 0.
    WK = 0.
    A  = 0.
    B = 0.
    goto 170
139 continue
    ICY = ICY+1
140 continue
    IBEGIN = ICY
    call MNCRCK(CRDBUF(IBEGIN:),MAXCWD,COMAND,LNC,MAXP,PLIST,LLIST, IERR)
    if (IERR .gt. 0)  goto 180
    UK = PLIST(1)
    WK = 0.
    if (LLIST .ge. 2)  WK = PLIST(2)
    A = 0.
    if (LLIST .ge. 3)  A = PLIST(3)
    B = 0.
    if (LLIST .ge. 4)  B = PLIST(4)
    goto 170
!          old (fixed-field) format
150 continue
    read (CRDBUF, 158,ERR=180)  XK,CNAMK,UK,WK,A,B
158 format (BN,F10.0, A10, 4F10.0)
    K = XK
    if (K .eq. 0)  goto 210
!          parameter format cracked, implement parameter definition
170 call minuit_mnparm(K,CNAMK,UK,WK,A,B,IERR)
    ICONDN = IERR
    return
!          format or other error
180 continue
    ICONDN = 1
    return
!        end of data
210 continue
    ICONDN = 2
    return
  end subroutine MNPARS

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine MNCRCK(CRDBUF,MAXCWD,COMAND,LNC,MXP,PLIST, LLIST,IERR)
!C
!C       Called from MNREAD.
!C       Cracks the free-format input, expecting zero or more
!C         alphanumeric fields (which it joins into COMAND(1:LNC))
!C         followed by one or more numeric fields separated by
!C         blanks and/or one comma.  The numeric fields are put into
!C         the LLIST (but at most MXP) elements of PLIST.
!C      IERR = 0 if no errors,
!C           = 1 if error(s).
!C      Diagnostic messages are written to minuit_fileout
!C
    integer(FINT), parameter :: MAXELM=25
    integer(FINT), parameter :: MXLNEL=19
    !
    !
    integer(FINT) :: lnc, llist, ierr, maxcwd, ielmnt, lend, nextb, ipos, ibegin
    integer(fint) :: iend, nelmnt, kcmnd, ic, left, ltoadd, ifld, nreq
    !
    !
    integer(fint), intent(in) :: mxp
    character(*) :: comand, crdbuf
    character(13) ::  cnumer !*13
    character(mxlnel), dimension(maxelm) :: celmnt !(maxelm)*(mxlnel)
    character(15) :: cnull !*15
    integer(fint), dimension(maxelm) :: lelmnt !(maxelm)
    real(freal), dimension(mxp) :: plist !(mxp)
    !dimension lelmnt(maxelm),plist(mxp)
    !
    data cnull /')null string   '/
    data cnumer/'123456789-.0+'/
    !
    ielmnt = 0
    lend = len(crdbuf)
    nextb = 1
    ierr = 0
!                                   . . . .  loop over words celmnt
10  continue
    do ipos= nextb,lend
       ibegin = ipos
       if (crdbuf(ipos:ipos).eq.' ') cycle
       if (crdbuf(ipos:ipos).eq.',')  goto 250
       goto 150
    end do
    goto 300
150 continue
!               found beginning of word, look for end
    do ipos = ibegin+1,lend
       if (crdbuf(ipos:ipos).eq.' ')  goto 250
       if (crdbuf(ipos:ipos).eq.',')  goto 250
    end do
    ipos = lend+1
250 iend = ipos-1
    ielmnt = ielmnt + 1
    if (iend .ge. ibegin) then
       celmnt(ielmnt) = crdbuf(ibegin:iend)
    else
       celmnt(ielmnt) = cnull
    endif
    lelmnt(ielmnt) = iend-ibegin+1
    if (lelmnt(ielmnt) .gt. mxlnel)  then
       write (minuit_fileout, 253) crdbuf(ibegin:iend),celmnt(ielmnt)
253    format (' minuit warning: input data word too long.'           &
            &   /'     original:',a                                            &
            &   /' truncated to:',a)
       lelmnt(ielmnt) = mxlnel
    endif
    if (ipos .ge. lend) goto 300
    if (ielmnt .ge. maxelm)  goto 300
!                     look for comma or beginning of next word
    do ipos= iend+1,lend
       if (crdbuf(ipos:ipos) .eq. ' ') cycle
       nextb = ipos
       if (crdbuf(ipos:ipos) .eq. ',') nextb = ipos+1
       goto 10
    end do
!                 all elements found, join the alphabetic ones to
!                                form a command
300 continue
    nelmnt = ielmnt
    comand = ' '
    lnc = 1
    plist(1) = 0.
    llist = 0
    if (ielmnt .eq. 0)  goto 900
    kcmnd = 0
    do ielmnt = 1, nelmnt
       if (celmnt(ielmnt) .eq. cnull)  goto 450
       do ic= 1, 13
          if (celmnt(ielmnt)(1:1) .eq. cnumer(ic:ic)) goto 450
       end do
       if (kcmnd .ge. maxcwd) cycle
       left = maxcwd-kcmnd
       ltoadd = lelmnt(ielmnt)
       if (ltoadd .gt. left) ltoadd=left
       comand(kcmnd+1:kcmnd+ltoadd) = celmnt(ielmnt)(1:ltoadd)
       kcmnd = kcmnd + ltoadd
       if (kcmnd .eq. maxcwd) cycle
       kcmnd = kcmnd + 1
       comand(kcmnd:kcmnd) = ' '
    end do
    lnc = kcmnd
    goto 900
450 continue
    lnc = kcmnd
!                      . . . .  we have come to a numeric field
    llist = 0
    do ifld= ielmnt,nelmnt
       llist = llist + 1
       if (llist .gt. mxp) then
          nreq = nelmnt-ielmnt+1
          write (minuit_fileout,511) nreq,mxp
511       format (/' minuit warning in mncrck: '/ ' command has input',i5,  &
               & ' numeric fields, but minuit can accept only',i3)
          goto 900
       endif
       if (celmnt(ifld) .eq. cnull)  then
          plist(llist) = 0.
       else
          read (celmnt(ifld), '(bn,f19.0)',err=575) plist(llist)
       endif
       cycle
575    write (minuit_fileout,'(A,A,A)') ' FORMAT ERROR IN NUMERIC FIELD: "',     &
            & CELMNT(IFLD)(1:LELMNT(IFLD)),'"'
       IERR = 1
       PLIST(LLIST) = 0.
    end do
!                                  end loop over numeric fields
900 continue
    if (LNC .le. 0)  LNC=1
    return
  end subroutine MNCRCK

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine MNCOMD(CRDBIN,ICONDN)
!C        Called by user.  'Reads' a command string and executes.
!C     Equivalent to MNEXCM except that the command is given as a
!C          character string.
!C
!C     ICONDN = 0: command executed normally
!C              1: command is blank, ignored
!C              2: command line unreadable, ignored
!C              3: unknown command, ignored
!C              4: abnormal termination (e.g., MIGRAD not converged)
!C              5: command is a request to read PARAMETER definitions
!C              6: 'SET INPUT' command
!C              7: 'SET TITLE' command
!C              8: 'SET COVAR' command
!C              9: reserved
!C             10: END command
!C             11: EXIT or STOP command
!C             12: RETURN command
!C
    integer(FINT) :: i
    real(FREAL), dimension(MAXP) :: PLIST !(MAXP)
    character(MAXCWD) :: COMAND !*(MAXCWD)
    character(26) ::  CLOWER, CUPPER !*26
    logical LEADER
    character(*) ::CRDBIN !*(*)
    character(100) :: CRDBUF !*100
    !
    !
    !integer(FINT) :: icondn, lenbuf, ipos, ic, isyswr, lnc, llist, ierr
    integer(FINT) :: icondn, lenbuf, ipos, ic, lnc, llist, ierr
    !
    !
    data CLOWER/'abcdefghijklmnopqrstuvwxyz'/
    data CUPPER/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
!
    LENBUF = len(CRDBIN)
    CRDBUF = CRDBIN
    ICONDN = 0
!     record not case-sensitive, get upper case, strip leading blanks
    LEADER = .true.
    IPOS = 1
    do I= 1, min(MAXCWD,LENBUF)
       if (CRDBUF(I:I) .eq. '''') exit
       if (CRDBUF(I:I) .eq. ' ')  then
          if (LEADER) IPOS = IPOS + 1
          cycle
       endif
       LEADER = .false.
       do IC= 1, 26
          if (CRDBUF(I:I) .eq. CLOWER(IC:IC)) CRDBUF(I:I)=CUPPER(IC:IC)
       end do
    end do
!                     blank or null command
    if (IPOS .gt. LENBUF)  then
       write (minuit_fileout,'(A)') ' BLANK COMMAND IGNORED.'
       ICONDN = 1
       GO TO 900
    endif
!                                           . .   preemptive commands
!               if command is 'PARAMETER'
    if (CRDBUF(IPOS:IPOS+2) .eq. 'PAR')    then
       ICONDN = 5
       LPHEAD = .true.
       GO TO 900
    endif
!               if command is 'SET INPUT'
    if (CRDBUF(IPOS:IPOS+6) .eq. 'SET INP')  then
       ICONDN = 6
       LPHEAD = .true.
       GO TO 900
    endif
!              if command is 'SET TITLE'
    if (CRDBUF(IPOS:IPOS+6) .eq. 'SET TIT')  then
       ICONDN = 7
       LPHEAD = .true.
       GO TO 900
    endif
!               if command is 'SET COVARIANCE'
    if (CRDBUF(IPOS:IPOS+6) .eq. 'SET COV')   then
       ICONDN = 8
       LPHEAD = .true.
       GO TO 900
    endif
!               crack the command . . . . . . . . . . . . . . . .
    call MNCRCK(CRDBUF(IPOS:LENBUF),MAXCWD,COMAND,LNC,MAXP,PLIST,LLIST,IERR)
    if (IERR .gt. 0) then
       write (minuit_fileout,'(A)') ' COMMAND CANNOT BE INTERPRETED'
       ICONDN = 2
       GO TO 900
    endif
!
    call minuit_mnexcm(COMAND(1:LNC),PLIST,LLIST,IERR)
    ICONDN = IERR
900 return
  end subroutine MNCOMD

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  logical function mnunpt(cfname) 
!           is .true. if cfname contains unprintable characters.        
    character(*),intent(in) :: cfname !*(*) 
    character(40), parameter :: cp1=' ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklm'
    character(40), parameter :: cp2='nopqrstuvwxyz1234567890./;:[]$%*_!@#&+()'
    character(80) :: cpt
    integer(FINT) :: i, l, ic
    !character cpt*80, cp1*40,cp2*40 
    !parameter (cp1=' ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklm') 
    !parameter (cp2='nopqrstuvwxyz1234567890./;:[]$%*_!@#&+()') 
    cpt=cp1//cp2 
    mnunpt = .false. 
    l = len(cfname) 
    do 100 i= 1, l 
       do ic= 1, 80 
          if (cfname(i:i) .eq. cpt(ic:ic))  go to 100 
       end do
       mnunpt = .true. 
       exit
100 end do
    return 
  end function mnunpt

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
  subroutine minuit_mnpout(iuext,chnam,val,err,xlolim,xuplim,iuint)
!c     user-called
!c   provides the user with information concerning the current status
!c          of parameter number iuext. namely, it returns:
!c        chnam: the name of the parameter
!c        val: the current (external) value of the parameter
!c        err: the current estimate of the parameter uncertainty
!c        xlolim: the lower bound (or zero if no limits)
!c        xuplim: the upper bound (or zero if no limits)
!c        iuint: the internal parameter number (or zero if not variable,
!c           or negative if undefined).
!c  note also:  if iuext is negative, then it is -internal parameter
!c           number, and iuint is returned as the external number.
!c     except for iuint, this is exactly the inverse of mnparm
!c
    !
    character(*) :: chnam
    integer(FINT), intent(in) :: iuext
    integer(FINT), intent(out) :: iuint
    real(FREAL), intent(out) :: val, err, xlolim, xuplim
    !
    integer(FINT) :: liint, iext, nvl
    !
    xlolim = 0.
    xuplim = 0.
    err = 0.
    if (iuext .eq. 0)  go to 100
    if (iuext .lt. 0)  then
!                   internal parameter number specified
       liint = -iuext
       if (liint .gt. npar) go to 100
       iext = nexofi(liint)
       iuint = iext
    else
!                    external parameter number specified
       iext = iuext
       if (iext .eq. 0)   go to 100
       if (iext .gt. nu)  go to 100
       liint = niofex(iext)
       iuint = liint
    endif
!                     in both cases
    nvl = nvarl(iext)
    if (nvl .lt. 0) go to 100
    chnam = cpnam(iext)
    val = u(iext)
    if (liint .gt. 0)  err = werr(liint)
    if (nvl .eq. 4) then
       xlolim = alim(iext)
       xuplim = blim(iext)
    endif
    return
!                parameter is undefined
100 iuint = -1
    chnam = 'undefined'
    val = 0.
    return
  end subroutine minuit_mnpout

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!

!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!


!
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!



  !
  !
end module minuit
