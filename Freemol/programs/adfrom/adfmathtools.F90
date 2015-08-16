!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H MODULE MATHTOOLS (Freemol1999 F.Mariotti Nov 1999)
!H To be inserted in inlcudes
!H-----------------------------------------------------------------------------
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H 
!
module adfmathtools
  !
  use vartypes
  !
  implicit none
  private
  !
  logical, save :: mathtools_debug = .true.
  integer(FINT), save :: mathtools_error = 0
  real(FREAL), private, parameter :: zero = 0.D0, one = 1.D0
  real(FREAL), private, parameter :: eps = epsilon(1.0_FREAL)
  !
  public math_a4pi50 !Partial port not tested
  public math_ang050 !Partial port not tested
  public math_ang110 !Partial port not tested
  public math_ang194 !Partial port not tested
  public math_euler  !convert port Not tested
  public math_ylm4pi !convert port Not tested
  public math_ylmw2  !convert port Not tested
  public math_invidx !convert port Not tested
  public math_permut
  public math_getorigin
  public math_lm
  public math_getdgt
  public math_factorials
  public math_getfact
  public math_bincoe
  public math_masscenter
  public quicksort
  public quicksortidx
  !
  !
  interface math_masscenter
     module procedure math_masscenter_r
     module procedure math_masscenter_i
  end interface
  !
  !
  interface quicksort
     module procedure quicksort_r
     module procedure quicksort_i
  end interface
  !
  !
  interface quicksortidx
     module procedure quicksortidx_r
     module procedure quicksortidx_i
     module procedure quicksortidx_nr
     module procedure quicksortidx_ni
  end interface
  !
contains
  !
!
!H-----------------------------------------------------------------------------
!H  subroutine math_masscenter_r(xyz,we,rvec,doit)
!H-----------------------------------------------------------------------------
!
  subroutine math_masscenter_r(xyz,we,rvec,doit)
!
!H-----------------------------------------------------------------------------
!H  
!H  We are right now limited to 3 dimensional vectors!!!!
!H  xyz       xyz cartesian
!H  we        Weight of each center !! you should give it!!!
!H  rvec      Position of the center of mass
!H  doit      Optinally do the translation!
!H  
!H-----------------------------------------------------------------------------
!
    !
    real(FREAL), dimension(:,:) :: xyz
    real(FREAL), dimension(:) :: we
    real(FREAL), dimension(3) :: rvec

    logical, optional :: doit

    rvec = 0.0_FREAL

    rvec(1) = sum(we*xyz(1,:))
    rvec(2) = sum(we*xyz(2,:))
    rvec(3) = sum(we*xyz(3,:))

    rvec = rvec / sum(we)

    if(rvec(1).lt.eps) rvec(1) = 0.0_FREAL
    if(rvec(2).lt.eps) rvec(2) = 0.0_FREAL
    if(rvec(3).lt.eps) rvec(3) = 0.0_FREAL

    if(present(doit)) then
       if (doit) then
          xyz(1,:) = xyz(1,:) - rvec(1)
          xyz(2,:) = xyz(2,:) - rvec(2)
          xyz(3,:) = xyz(3,:) - rvec(3)
       endif
    endif

  end subroutine math_masscenter_r
!
!H-----------------------------------------------------------------------------
!H  subroutine math_masscenter_i(xyz,we,rvec,doit)
!H-----------------------------------------------------------------------------
!H  
!H  Optionally use integer weights!!
!H  
!H  
!H-----------------------------------------------------------------------------
!H  
!
  subroutine math_masscenter_i(xyz,we,rvec,doit)

    integer(FINT), dimension(:) :: we
    real(FREAL), dimension(:,:) :: xyz
    real(FREAL), dimension(3) :: rvec

    logical, optional :: doit

    rvec = 0.0_FREAL

    rvec(1) = sum(real(we,KIND=FREAL)*xyz(1,:))
    rvec(2) = sum(real(we,KIND=FREAL)*xyz(2,:))
    rvec(3) = sum(real(we,KIND=FREAL)*xyz(3,:))

    rvec = rvec / real(sum(we),KIND=FREAL)

    if(rvec(1).lt.eps) rvec(1) = 0.0_FREAL
    if(rvec(2).lt.eps) rvec(2) = 0.0_FREAL
    if(rvec(3).lt.eps) rvec(3) = 0.0_FREAL

    if(present(doit)) then
       if (doit) then
          xyz(1,:) = xyz(1,:) - rvec(1)
          xyz(2,:) = xyz(2,:) - rvec(2)
          xyz(3,:) = xyz(3,:) - rvec(3)
       endif
    endif

  end subroutine math_masscenter_i
!
!H      subroutine math_a4pi50(x,y,z,w)
!H
!W    TODO!!
!
  subroutine math_a4pi50(x, y, z, w)
    !
    real(FREAL), dimension(:), intent(out) :: x,y,z
    real(FREAL), dimension(:), intent(out) :: w
    !
    real(FREAL) :: pi4
    !
    real(FREAL) :: wsum
    !
    !
    !-----subroutine arguments are   x(n),y(n),z(n),w(n)
    !-----where x,y,z are coordinates of mesh points (on 4pi sphere)
    !-----and w is the associated integration weight (sum{w(j),j=1,50}=4pi).
    !
    call math_ang050(x,y,z,w)
    !
    pi4 = 16.D0*atan(1.D0)
    w = pi4 * w
    wsum = sum(w)
    !
    if(mathtools_debug) then
       if((wsum-pi4).gt.(1000.0D0*epsilon(1.0D0))) then
          write (0,*) ' Error in math_a4pi50 ',wsum-pi4
       endif
    endif
    !
    return
  end subroutine math_a4pi50
!
!
!
!
!***********************************************************************
!
!  gauss integration meshes on surface of unit sphere,
!  formulae of v.i. lebedev.
!
!  programmed by a. d. becke, 1986.
!
!  some constants extended to ca. 16 decimal digits by ross m. dickson,
!  july-aug 1991.
!
!  literature citations are:
!
!  v. i. lebedev, "values of the nodes and weights of ninth to
!  seventeenth order gauss-markov quadrature formulae invariant under
!  the octahedron groups with inversion", zh. vychisl. mat. mat. fiz.
!  15(1), 48--54 (1975).
!
!  v. i. lebedev, "quadratures on a sphere", zh. vychisl. mat. mat. fiz.
!  16(2), 293--306 (1976).
!
!  v. i. lebedev, "spherical quadrature formulas exact to orders 25--29"
!  sibirsk. mat. zh. 18(1), 132--142 (1977).
!
!  a. d. becke, "a fully numerical integration scheme for polyatomic
!  molecules", j. chem. phys. 88(4), 2547--2553 (1988).
!
!
!     subroutine math_ang050 : 50 points, order l=11.
!     subroutine math_ang110 : 110 points, order l=17.
!     subroutine math_ang194 : 194 points, order l=23.
!     subroutine math_ang302 : 302 points, order l=29.
!
!     subroutine arguments are   x(n),y(n),z(n),w(n)
!     where x,y,z are coordinates of mesh points (on unit sphere)
!     and w is the associated integration weight.
!
  subroutine math_ang050(x,y,z,w)

    real(FREAL), dimension(:), intent(out) :: x, y, z
    real(FREAL), dimension(:), intent(out) :: w

    real(FREAL), parameter :: a1 = 1.26984126984126984d-02
    real(FREAL), parameter :: a2 = 2.25749559082892416d-02
    real(FREAL), parameter :: a3 = 2.10937500000000000d-02
    real(FREAL), parameter :: b1 = 2.01733355379188713d-02

    integer(FINT), dimension(50) :: ix, iy, iz
    real(FINT), dimension(50) :: wi

    real(FREAL) :: l
    real(FREAL) :: m
    real(FREAL) :: rt2
    real(FREAL) :: rt3
    integer(FINT) :: i

    data l /0.301511344577763574d0/
    data m /0.904534033733290888d0/

    data(wi(i),i= 1,50) /6 * a1, 12 * a2, 8 * a3, 24 * b1/

    data(ix(i),i= 1, 6) /1,-1, 0, 0, 0, 0/
    data(iy(i),i= 1, 6) /0, 0, 1,-1, 0, 0/
    data(iz(i),i= 1, 6) /0, 0, 0, 0, 1,-1/

    data(ix(i),i= 7,18) /0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1/
    data(iy(i),i= 7,18) /1,-1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1/
    data(iz(i),i= 7,18) /1, 1,-1,-1, 1, 1,-1,-1, 0, 0, 0, 0/
    data(ix(i),i=19,26) /1,-1, 1,-1, 1,-1, 1,-1/
    data(iy(i),i=19,26) /1, 1,-1,-1, 1, 1,-1,-1/
    data(iz(i),i=19,26) /1, 1, 1, 1,-1,-1,-1,-1/

    data(ix(i),i=27,50)/1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1/
    data(iy(i),i=27,50)/1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1/
    data(iz(i),i=27,50)/1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1/

    rt2 = 1.d0/sqrt(2.d0)
    rt3 = 1.d0/sqrt(3.d0)

!do i=1,6
    x(1:6) = ix(1:6)
    y(1:6) = iy(1:6)
    z(1:6) = iz(1:6)

!do i=7,18
    x(7:18) = rt2 * ix(7:18)
    y(7:18) = rt2 * iy(7:18)
    z(7:18) = rt2 * iz(7:18)

!do i=19,26
    x(19:26) = rt3 * ix(19:26)
    y(19:26) = rt3 * iy(19:26)
    z(19:26) = rt3 * iz(19:26)
!
    x(27:42) = l * ix(27:42)
    y(27:34) = l * iy(27:34)
    z(27:34) = m * iz(27:34)

    y(35:42) = m * iy(35:42)
    z(35:50) = l * iz(35:50)

    x(43:50) = m * ix(43:50)
    y(43:50) = l * iy(43:50)

!weights
    w(1:50) = wi

    return
  end subroutine math_ang050
!

  subroutine math_ang110(x, y, z, w)

    real(FREAL), dimension(:) :: x, y, z
    real(FREAL), dimension(:), intent(inout) :: w

    integer(FINT), dimension(110) :: ix, iy, iz
    real(FINT), dimension(110) :: wi

    integer(FINT) :: i

    real(FREAL) :: rt3, p1, q1
    real(FREAL), dimension(3) :: l, m

    real(FREAL), parameter :: a1 = 3.82827049493716160d-03
    real(FREAL), parameter :: a3 = 9.79373751248751249d-03
    real(FREAL), parameter :: b1 = 8.21173728319111098d-03
    real(FREAL), parameter :: b2 = 9.59547133607096285d-03
    real(FREAL), parameter :: b3 = 9.94281489117810328d-03
    real(FREAL), parameter :: c1 = 9.69499636166302833d-03

    data(wi(i),i=1,110) /6 * a1, 8 * a3, 24 * b1, 24 * b2, 24 *b3, 24 * c1/

    data l(1) /0.185115635344736169d+00/
    data l(2) /0.395689473055941909d+00/
    data l(3) /0.690421048382292178d+00/
    data m(1) /0.965124035086594106d+00/
    data m(2) /0.828769981252592211d+00/
    data m(3) /0.215957291845848832d+00/
    data p1 /0.878158910604066133d0/
    data q1 /0.478369028812150197d0/

    data(ix(i),i= 1, 6) /1,-1, 0, 0, 0, 0/
    data(iy(i),i= 1, 6) /0, 0, 1,-1, 0, 0/
    data(iz(i),i= 1, 6) /0, 0, 0, 0, 1,-1/

    data(ix(i), i = 7, 14) / 1,-1, 1,-1, 1,-1, 1,-1 /
    data(iy(i), i = 7, 14) / 1, 1,-1,-1, 1, 1,-1,-1 /
    data(iz(i), i = 7, 14) / 1, 1, 1, 1,-1,-1,-1,-1 /

    data(ix(i),i=15,86)/&
         1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,&
         1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,&
         1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1/
    data(iy(i),i=15,86)/&
         1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,&
         1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,&
         1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1/
    data(iz(i),i=15,86)/&
         1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,&
         1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,&
         1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1/

    data(ix(i),i=87,110)/1,-1,1,-1,1,-1,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0/
    data(iy(i),i=87,110)/1,1,-1,-1,0,0,0,0,1,-1,1,-1,1,1,-1,-1,0,0,0,0,1,-1,1,-1/
    data(iz(i),i=87,110)/0,0,0,0,1,1,-1,-1,1,1,-1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1/

    rt3 = 1.d0/sqrt(3.d0)

!do i=1,6
    x(1:6) = ix(1:6)
    y(1:6) = iy(1:6)
    z(1:6) = iz(1:6)

!do i=7,14
    x(7:14) = rt3 * ix(7:14)
    y(7:14) = rt3 * iy(7:14)
    z(7:14) = rt3 * iz(7:14)

!15:22,39:46,63:70
    x(15:22) = l(1) * ix(15:22)
    y(15:22) = l(1) * iy(15:22)
    z(15:22) = m(1) * iz(15:22)

    x(39:46) = l(1) * ix(39:46)
    y(39:46) = l(1) * iy(39:46)
    z(39:46) = m(1) * iz(39:46)

    x(63:70) = l(1) * ix(63:70)
    y(63:70) = l(1) * iy(63:70)
    z(63:70) = m(1) * iz(63:70)

!23:30,47:54,71:78
    x(23:30) = l(2) * ix(23:30)
    y(23:30) = m(2) * iy(23:30)
    z(23:30) = l(2) * iz(23:30)

    x(47:54) = l(2) * ix(47:54)
    y(47:54) = m(2) * iy(47:54)
    z(47:54) = l(2) * iz(47:54)

    x(71:78) = l(2) * ix(71:78)
    y(71:78) = m(2) * iy(71:78)
    z(71:78) = l(2) * iz(71:78)

!31:38,55:62,79:86
    x(31:38) = m(3) * ix(31:38)
    y(31:38) = l(3) * iy(31:38)
    z(31:38) = l(3) * iz(31:38)

    x(55:62) = m(3) * ix(55:62)
    y(55:62) = l(3) * iy(55:62)
    z(55:62) = l(3) * iz(55:62)

    x(79:86) = m(3) * ix(79:86)
    y(79:86) = l(3) * iy(79:86)
    z(79:86) = l(3) * iz(79:86)

!87:90
    x(87:90) = p1 * ix(87:90)
    y(87:90) = q1 * iy(87:90)
    z(87:90) = 0.d0

!91:94
    x(91:94) = p1 * ix(91:94)
    y(91:94) = 0.d0
    z(91:94) = q1 * iz(91:94)

!95:98
    x(95:98) = 0.d0
    y(95:98) = p1 * iy(95:98)
    z(95:98) = q1 * iz(95:98)

!99:102
    x(99:102) = q1 * ix(99:102)
    y(99:102) = p1 * iy(99:102)
    z(99:102) = 0.d0

!103:106
    x(103:106) = q1 * ix(103:106)
    y(103:106) = 0.d0
    z(103:106) = p1 * iz(103:106)

!107:110
    x(107:110) = 0.d0
    y(107:110) = q1 * iy(107:110)
    z(107:110) = p1 * iz(107:110)

    w(1:110) = wi

    return
  end subroutine math_ang110
!
  subroutine math_ang194(x, y, z, w)

    real(FREAL), dimension(:) :: x, y, z
    real(FREAL), dimension(:), intent(inout) :: w

    integer(FINT), dimension(194) :: ix, iy, iz
    real(FINT), dimension(194) :: wi

    integer(FINT) :: i

    real(FREAL) :: rt2, rt3, p1, q1, r1, s1, t1
    real(FREAL), dimension(4) :: l, m

    real(FREAL), parameter :: a1 = 1.78234044724461116d-03
    real(FREAL), parameter :: a2 = 5.71690594997710189d-03
    real(FREAL), parameter :: a3 = 5.57338317884873797d-03
    real(FREAL), parameter :: b1 = 5.51877146727d-03
    real(FREAL), parameter :: b2 = 5.15823771181d-03
    real(FREAL), parameter :: b3 = 5.60870408259d-03
    real(FREAL), parameter :: b4 = 4.10677702817d-03
    real(FREAL), parameter :: c1 = 5.05184606461480848d-03
!
!-----note the typographical error in lebedev's second paper:  the
!-----rational expression on page 18 for d1 should contain 2 to the
!-----ninth power rather than 2 to the sixth power in the denominator.
!
    real(FREAL), parameter :: d1 = 5.530248916233093701d-03

    data(wi(i),i=1,194)/6*a1,12*a2,8*a3,24*b1,24*b2,24*b3,24*b4,24*c1,48*d1/

    data(l(i),i=1,4)/&
         0.444693317871743737d0,0.289246562757543879d0,&
         0.671297344269522658d0,0.129933544765006692d0/
    data(m(i),i=1,4)/&
         0.777493219314767121d0,0.912509096867473678d0,&
         0.314196994182585949d0,0.982972302707253295d0/

    data p1/0.938319218137591521d0/
    data q1/0.345770219761128271d0/
    data r1/0.836036015482458872d0/
    data s1/0.159041710538352976d0/
    data t1/0.525118572443642018d0/

    data(ix(i),i=1,6)/1,-1,0,0,0,0/
    data(iy(i),i=1,6)/0,0,1,-1,0,0/
    data(iz(i),i=1,6)/0,0,0,0,1,-1/

    data(ix(i),i=7,18)/0,0,0,0,1,-1,1,-1,1,-1,1,-1/
    data(iy(i),i=7,18)/1,-1,1,-1,0,0,0,0,1,1,-1,-1/
    data(iz(i),i=7,18)/1,1,-1,-1,1,1,-1,-1,0,0,0,0/

    data(ix(i),i=19,26)/1,-1,1,-1,1,-1,1,-1/
    data(iy(i),i=19,26)/1,1,-1,-1,1,1,-1,-1/
    data(iz(i),i=19,26)/1,1,1,1,-1,-1,-1,-1/

    data(ix(i),i=27,122)/&
         1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,&
         1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,&
         1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,&
         1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1/
    data(iy(i),i=27,122)/&
         1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,&
         1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,&
         1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,&
         1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1/
    data(iz(i),i=27,122)/&
         1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,&
         1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,&
         1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,&
         1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1/
    data(ix(i),i=123,146)/1,-1,1,-1,1,-1,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0/
    data(iy(i),i=123,146)/1,1,-1,-1,0,0,0,0,1,-1,1,-1,1,1,-1,-1,0,0,0,0,1,-1,1,-1/
    data(iz(i),i=123,146)/0,0,0,0,1,1,-1,-1,1,1,-1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1/
    data(ix(i),i=147,194)/&
         1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,&
         1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,1,-1,1,-1/
    data(iy(i),i=147,194)/&
         1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,&
         1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1,  1,1,-1,-1,1,1,-1,-1/
    data(iz(i),i=147,194)/&
         1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,&
         1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1,  1,1,1,1,-1,-1,-1,-1/

    rt2 = 1.d0/sqrt(2.d0)
    rt3 = 1.d0/sqrt(3.d0)

!do i=1,6
    x(1:6) = ix(1:6)
    y(1:6) = iy(1:6)
    z(1:6) = iz(1:6)

!do i=7,18
    x(7:18) = rt2 * ix(7:18)
    y(7:18) = rt2 * iy(7:18)
    z(7:18) = rt2 * iz(7:18)

!do i=19,26
    x(19:26) = rt3 * ix(19:26)
    y(19:26) = rt3 * iy(19:26)
    z(19:26) = rt3 * iz(19:26)
!
    x(27:42) = l(1) * ix(27:42)
    y(27:34) = l(1) * iy(27:34)
    z(27:34) = m(1) * iz(27:34)

    y(35:42) = m(1) * iy(35:42)
    z(35:50) = l(1) * iz(35:50)

    x(43:50) = m(1) * ix(43:50)
    y(43:50) = l(1) * iy(43:50)
!
    x(51:66) = l(2) * ix(51:66)
    y(51:58) = l(2) * iy(51:58)
    z(51:58) = m(2) * iz(51:58)

    y(59:66) = m(2) * iy(59:66)
    z(59:74) = l(2) * iz(59:74)

    x(67:74) = m(2) * ix(67:74)
    y(67:74) = l(2) * iy(67:74)
!
    x(75:90) = l(3) * ix(75:90)
    y(75:82) = l(3) * iy(75:82)
    z(75:82) = m(3) * iz(75:82)

    y(83:90) = m(3) * iy(83:90)
    z(83:98) = l(3) * iz(83:98)

    x(91:98) = m(3) * ix(91:98)
    y(91:98) = l(3) * iy(91:98)
!
    x(99:114) = l(4) * ix(99:114)
    y(99:106) = l(4) * iy(99:106)
    z(99:106) = m(4) * iz(99:106)

    y(107:114) = m(4) * iy(107:114)
    z(107:122) = l(4) * iz(107:122)

    x(115:122) = m(4) * ix(115:122)
    y(115:122) = l(4) * iy(115:122)

!123:126
    x(123:126) = p1 * ix(123:126)
    y(123:126) = q1 * iy(123:126)
    z(123:126) = 0.d0

!127:130
    x(127:130) = p1 * ix(127:130)
    y(127:130) = 0.d0
    z(127:130) = q1 * iz(127:130)

!131:134
    x(131:134) = 0.d0
    y(131:134) = p1 * iy(131:134)
    z(131:134) = q1 * iz(131:134)

!135:138
    x(135:138) = q1 * ix(135:138)
    y(135:138) = p1 * iy(135:138)
    z(135:138) = 0.d0

!139:142
    x(139:142) = q1 * ix(139:142)
    y(139:142) = 0.d0
    z(139:142) = p1 * iz(139:142)

!143:146
    x(143:146) = 0.d0
    y(143:146) = q1 * iy(143:146)
    z(143:146) = p1 * iz(143:146)

!147:154
    x(147:154) = r1 * ix(147:154)
    y(147:154) = s1 * iy(147:154)
    z(147:154) = t1 * iz(147:154)

!155:162
    x(155:162) = r1 * ix(155:162)
    y(155:162) = t1 * iy(155:162)
    z(155:162) = s1 * iz(155:162)

!163:170
    x(163:170) = s1 * ix(163:170)
    y(163:170) = r1 * iy(163:170)
    z(163:170) = t1 * iz(163:170)

!171:178
    x(171:178) = s1 * ix(171:178)
    y(171:178) = t1 * iy(171:178)
    z(171:178) = r1 * iz(171:178)

!179:186
    x(179:186) = t1 * ix(179:186)
    y(179:186) = r1 * iy(179:186)
    z(179:186) = s1 * iz(179:186)

!187:194
    x(187:194) = t1 * ix(187:194)
    y(187:194) = s1 * iy(187:194)
    z(187:194) = r1 * iz(187:194)

    w(1:194) = wi

    return
  end subroutine math_ang194
!
!
!H
!H----------------------------------------------------------------------
!H    SUBROUTINE MATH_EULER(ANGLE,ROTAX,INV,ALPHA,BETA,GAMMA,RMAT)
!H----------------------------------------------------------------------
!H
!
  subroutine math_euler(angle, rotax, inv, alpha, beta, gamma, rmat,leps)
    !
    !H----------------------------------------------------------------------
    !H
    !H    calculates the euler angles 'alpha,beta,gamma' and the 3-d
    !H    rotation matrix 'rmat' corresponding to a rotation by 'angle'
    !H    radians about the rotation axis through point 'x,y,z' and
    !H    the origin. if 'inv'=.true. , an inversion is performed
    !H    following the pure rotation.
    !H
    !H    alpha = 1st rotation: about the z-axis,    -pi.lt.alpha.le.pi
    !H    beta  = 2nd rotation: about the y'-axis,    0 .le.beta .le.pi
    !H    gamma = 3rd rotation: about the z"-axis,   -pi.lt.gamma.le.pi
    !H
    !H    alpha = 1st rotation: about the z-axis,    -pi.< .alpha.<=.pi
    !H    beta  = 2nd rotation: about the y'-axis,    0 .<=.beta .<=.pi
    !H    gamma = 3rd rotation: about the z"-axis,   -pi.< .gamma.<=.pi
    !H
    !H    all rotations are passive (rotations of the coord system
    !H    rather than the object vectors). the coord system is assumed to
    !H    be right-handed, and rotations are positive if counterclockwise.
    !H    the axis vector 'rotax' is returned normalized and oriented so
    !H    that z is positive. (if z=0, so that y is pos; if z & y =0, so
    !H    that x is pos.) angle is returned with the appropriate sign,
    !H    within the range -pi.lt.angle.le.pi
    !H
    !H    the euler angles are extracted from the rotation matrix which
    !H    transforms z-axis into the rot'n axis, rotates by angle, and
    !H    then inverts the first transformation.
    !H
    !H    cf. rose,  p. 50ff, p.65
    !H
    !H    called by: symops,closur,molrot
    !H
    !H----------------------------------------------------------------------
    !
    !W    F.Mariotti LEPS parameter in this case can be the machine precision
    !W    and not the symmetry precision!!
    !
    !
    !
    !
    logical, intent(in) :: inv
    real(FREAL) :: angle
    real(FREAL) :: leps
    real(FREAL) :: alpha
    real(FREAL) :: beta
    real(FREAL) :: gamma
    real(FREAL), dimension(3) :: rotax
    real(FREAL), dimension(3,3) :: rmat
    !
    !     FROM INLINE CODE
    !
    integer(FINT) :: i
    integer(FINT) :: j
    integer(FINT) :: nmod
    real(FREAL) :: a
    real(FREAL) :: b
    real(FREAL) :: ca
    real(FREAL) :: cb
    real(FREAL) :: cosa
    real(FREAL) :: cosb
    real(FREAL) :: cosg
    real(FREAL) :: ct
    real(FREAL) :: rax
    real(FREAL) :: sa
    real(FREAL) :: sb
    real(FREAL) :: sina
    real(FREAL) :: sinb
    real(FREAL) :: sing
    real(FREAL) :: st
    real(FREAL) :: t
    real(FREAL) :: xd1
    real(FREAL) :: xd2
    real(FREAL) :: xn1
    real(FREAL) :: xn2
    real(FREAL) :: xyproj
    real(FREAL) :: x
    real(FREAL) :: y
    real(FREAL) :: z
    real(FREAL) :: pi
    real(FREAL) :: twopi
    !
    !
    !
    pi = 4.0D0 * atan(1.0D0)
    twopi = 2.D0 * pi
    !
    a = zero
    b = zero
    t = zero
    x = rotax(1)
    y = rotax(2)
    z = rotax(3)
    !
    !-----normalize the rotation axis vector ------------------------------
    !
    rax = sqrt(x * x + y * y + z * z)
    !
    !-----if no axis of rotation was specified, error stop ----------------
    !
    if (rax.lt.leps) then
       print *,'DEBUG Euler: undefined rotational axis.'
    endif
    !
    x = x / rax
    y = y / rax
    z = z / rax
    if (abs(x).lt.leps) x = zero
    if(abs(y).lt.leps) y = zero
    if (abs(z).lt.leps) z = zero
    !
    !-----bring angle within the range -pi to pi --------------------------
    !
    nmod = int(angle/twopi)
    angle = angle-real(nmod,KIND=FREAL) * twopi
    if (angle.gt.pi) angle = -(twopi - angle)
    if (angle.lt. - pi) angle = twopi + angle
    !
    !-----direct the rotation axis vector into the upper
    !     hemisphere if it points toward negative z ---------------------
    !
    !
    !OLDCODE
    !      if (z.gt.zero)                             go to 10
    !      if (z.eq.zero.and.y.gt.zero)               go to 10
    !      if (z.eq.zero.and.y.eq.zero.and.x.gt.zero) go to 10
    !
    !      x=-x
    !      y=-y
    !      z=-z
    !      angle=-angle
    ! 10   if (abs(abs(angle)-pi).lt.leps) angle=pi
    !NEWLINES
    if (z.le.zero) then
       if (.not.(z.eq.zero.and.y.gt.zero) ) then
          if (.not.(z.eq.zero.and.y.eq.zero.and.x.gt.zero) ) then
             x = - x
             y = - y
             z = - z
             angle = - angle
          endif
       endif
    endif
    !
    if (abs(abs(angle)-pi).lt.leps) angle = pi
    !NEWLINES
    !
    !
    if (abs(angle).lt.leps) angle = zero
    xyproj = sqrt(x * x + y * y)
    rotax(1) = x
    rotax(2) = y
    rotax(3) = z
    !
    !-----find rot angles to transform z-axis into rotax,
    !     rotate by angle, & transform back: (ab)**-1 * t * (ab) --------
    !
    if (xyproj.ne.zero) then
       if (y.eq.zero.and.x.lt.zero) a = pi
       if (y.ne.zero) a = - atan2( - y, x)
       b = atan2(xyproj, z)
    endif
    !
    t = angle
    !
    !-----extract euler angles for the net rotation -----------------------
    !
    sa = sin(a)
    ca = cos(a)
    sb = sin(b)
    cb = cos(b)
    st = sin(t)
    ct = cos(t)
    beta = acos(one-(one-ct) * sb * sb)
    !
    !-----special cases: beta=zero or beta=pi ( sin(beta)=zero ) ----------
    !
    if (beta.lt.leps) then
       alpha = zero
       beta = zero
       gamma = angle
    elseif ((pi - beta) .lt.leps) then
       alpha = - atan2(x, y)
       beta = pi
       gamma = - alpha
    else
       xn1 = sa * cb *(one-ct)
       xn2 = ca * st
       xd1 = ca * cb *(one-ct)
       xd2 = sa * st
       alpha = atan2((xn1 - xn2),(xd1 + xd2))
       gamma = atan2((xn1 + xn2),(xd2 - xd1))
    endif
    !
    !-----------------------------------------------------------------------
    !     construct the 3-d rotation matrix for this set of euler angles
    !-----------------------------------------------------------------------
    !
    sina = sin(alpha)
    cosa = cos(alpha)
    sinb = sin(beta)
    cosb = cos(beta)
    sing = sin(gamma)
    cosg = cos(gamma)
    rmat(1, 1) = cosa * cosb * cosg - sina * sing
    rmat(1, 2) = sina * cosb * cosg + cosa * sing
    rmat(1, 3) = - sinb * cosg
    rmat(2, 1) = - cosa * cosb * sing - sina * cosg
    rmat(2, 2) = - sina * cosb * sing + cosa * cosg
    rmat(2, 3) = sinb * sing
    rmat(3, 1) = cosa * sinb
    rmat(3, 2) = sina * sinb
    rmat(3, 3) = cosb
    !
    !-----if this is an improper rotation, invert through the origin ------
    !
    if (inv) then
       rmat = rmat * (-1.0D0)
    endif
    !
    !
    !
    return
  end subroutine math_euler
!
!
!
!
!
!H    ******************************************************************
!H    *       DOUBLE PRECISION FUNCTION GETDGT(A)                      *
!H    ******************************************************************
!H
  real(FREAL) function math_getdgt(a)
    !
    real(FREAL), intent(in):: a
    !
    math_getdgt = a - real(int(a),KIND=FREAL)
    return
  end function math_getdgt
!
!
!
!
!
!
!
!H
!H    SUBROUTINE MATH_YLM4PI(LMAX,Z,PHI,YL)
!H
!H
!H
!H
!H
!
  subroutine math_ylm4pi(lmax, z, phi, yl)
!
!
!H
!H    ******************************************************************
!H
!H       generates the real spherical harmonics from l=0 to lmax,
!H     m=-l to l, with cos(theta)=z. m varies most rapidly.
!H       cos(mtheta) form has m>0, sin(mtheta) form has m<0.
!H======
!Hnote_1: sqrt{4pi/(2l+1)}*ylm are calculated! i.e. spherical harmonics
!H======  normalized to 4pi instead of (2l+1)/4pi or 1 !
!H        see comments in program.
!H======
!Hnote_2: phase convention: (-1)**m*ylm(abramowitz & stegun)
!H======  i.e.:  y(1,-1)=p_y=cte*y and y(1,1)=p_x=cte*x ; etc. ...
!H        see:     abramowitz & stegun eq. 8.5.3
!H        see comments in program.
!H
!H    ******************************************************************
!H
!

    real(FREAL), dimension(:), intent(out) :: yl
    real(FREAL) :: pi2in
    real(FREAL) :: y0, zero, one, two
    real(FREAL) :: z, phi
!
    integer(FINT) :: lmax
    real(FREAL) :: x
    real(FREAL) :: sinphi
    real(FREAL) :: cosphi
    real(FREAL) :: sinmp
    real(FREAL) :: cosmp
    real(FREAL) :: pmm
    real(FREAL) :: facm
!Normalization to one
!      real(FREAL) :: pi4
    real(FREAL) :: faclm
    integer(FINT) :: mfac2
    integer(FINT) :: mfac
    integer(FINT) :: m
    integer(FINT) :: istor1
    integer(FINT) :: istor2
    integer(FINT) :: mp1
    integer(FINT) :: lfac
    real(FREAL) :: plm
    integer(FINT) :: l
    integer(FINT) :: lmm
    integer(FINT) :: lpm
    real(FREAL) :: plm1
    real(FREAL) :: plm0
    !integer(FINT) :: yl1
    real(FREAL) :: yl1
    real(FREAL) :: costp
    integer(FINT) :: iy
!
    data pi2in / 1.591549430918954d-1 /
    data y0, zero, one, two / 0.282094791773878d0, 0.d0, 1.d0, 2.d0 /
!
!-----Choose normalization
!
!-----nomalisation to 1/sqrt(4pi)
    yl(1) = y0
!-----nomalisation to one
!H    yl(1)=one
!
!
    if (lmax.eq.0) return
    x = sqrt(one-z * z)
    sinphi = sin(phi)
    cosphi = cos(phi)
    sinmp = zero
    cosmp = one
    pmm = one
    facm = pi2in
!Normalization to one
!      pi4 = two / pi2in
    faclm = facm / two
    mfac2 = 0
    mfac = 1
    m = 0
!
!-----p(m,m) calculated by recursion ----------------------------------
!
    istor1 = 1
    istor2 = 1
    do mp1 = 1, lmax
       lfac = mfac
       plm = pmm
       mfac2 = mfac2 + 2
!
!-----recursion for p(l,m), l=m+1 to lmax -----------------------------
!
       do l = mp1, lmax
          istor1 = istor1 + l + l
          istor2 = istor2 + l + l
          lmm = l - m
          lpm = l + m
          plm1 = plm * z * real(lfac,KIND=FREAL)
          lfac = lfac + 2
          if (l.ne.mp1) then
             plm1 = (plm1 -(real(lpm-1,KIND=FREAL))*plm0) / real(lmm,KIND=FREAL)
          endif
          if (m.ne.0) then
             faclm = real(lmm,KIND=FREAL) * faclm /real(lpm,KIND=FREAL)
          endif
!     nomalisation to 1/sqrt(4pi)
          yl1 = sqrt(real(lfac,KIND=FREAL) * faclm) * plm1
!     nomalisation to one
!     yl1=sqrt(real(lfac,KIND=FREAL)*pi4*faclm)*plm1
          yl(istor2) = yl1 * cosmp
          if (m.ne.0) yl(istor1) = yl1 * sinmp
          yl(istor2) = yl1 * cosmp
          if (m.ne.0) yl(istor1) = yl1 * sinmp
          plm0 = plm
          plm = plm1
       end do
       m = mp1
       istor1 = m * m + 1
       istor2 = istor1 + m + m
!
!-----calculate l=m element and sin, cos mphi for next m iteration ----
!
       facm = facm / (real(mfac * mfac2,KIND=FREAL) )
       faclm = facm
       pmm = - pmm * real(mfac,KIND=FREAL) * x
       costp = cosmp * cosphi - sinmp * sinphi
       sinmp = sinmp * cosphi + cosmp * sinphi
       cosmp = costp
       mfac = mfac + 2
!     nomalisation to 1/sqrt(4pi)
       yl1 = pmm * sqrt(faclm * real(mfac,KIND=FREAL) )
!     nomalisation to one
!     yl1=pmm*sqrt(faclm*pi4*real(mfac,KIND=FREAL))
       yl(istor2) = yl1 * cosmp
       yl(istor1) = yl1 * sinmp
    end do
!
!---  phase is according to abramowitz and stegun or frank harris ----
!---  adapt phase to modern convention
!
    iy = 0
    do l = 0, lmax
       do m = - l, l
          iy = iy + 1
          if (mod(abs(m),2).eq.1) then
             yl(iy) = - yl(iy)
          endif
       enddo
    enddo
    return
!
  end subroutine math_ylm4pi
!
!
!H    SUBROUTINE MATH_YLMW2(Z,PHI,YL,LMAX)
!H
!
  subroutine math_ylmw2(z,phi,yl,lmax)
!
!H********************************************************************
!H                                                                   *
!H     ylmw calculates real spherical harmonics directly rather      *
!H     than recursively as in ylm.  it is faster but is currently    *
!H     limited to lmax < 9 ; it is suspected that numerical          *
!H     difficulties may arise for greater l values.                  *
!H                                                                   *
!H     the real spherical harmonics are stored from l=0 to lmax,     *
!H     m= 0 to l, with cos(theta)=z. m varies most rapidly.          *
!H     cos(mtheta) are stored in yl(...,1), sin(mtheta) in yl(...,2) *
!H                                                                   *
!H     note: ylm are normalized to (2l+1)/4pi !                      *
!H                                                                   *
!H********************************************************************
!H
!
    integer(FINT), intent(in) :: lmax
    real(FREAL), dimension(:), intent(out) :: yl

    integer(FINT) :: i
    real(FREAL), dimension(45) :: fac

!
! From Inline Code
!
    real(FREAL) :: n,z,phi,a,x,zz,xx,sinphi,cosphi,y,sn2phi,cs2phi
    real(FREAL) :: sn3phi,cs3phi,sn4phi,cs4phi,sn5phi,cs5phi,z3,z4
    real(FREAL) :: x3,x4,sn6phi,cs6phi,sn7phi,cs7phi,sn8phi,cs8phi


    !
    !     fac contains the normalization factors of the associated
    !     legendre polynomials.
    !
    data fac / .282094792d0, .488602512d0, .488602512d0, .315391565d0,&
              1.092548431d0, .546274215d0, .373176333d0, .457045800d0,&
              1.445305721d0, .590043590d0, .105785547d0, .669046544d0,&
               .473087348d0,1.770130770d0, .625835736d0, .116950321d0,&
               .452946651d0,2.396768392d0, .489238300d0,2.075662315d0,&
               .656382057d0, .063569202d0, .582621362d0, .460602630d0,&
               .921205260d0, .504564901d0,2.366619162d0, .683184105d0,&
               .068284277d0, .090331608d0, .221266346d0, .156458934d0,&
              1.037831158d0, .518915579d0,2.645960663d0, .707162733d0,&
               .009086770d0, .012115694d0, .065164608d0, .176466594d0,&
               .068345218d0, .037911105d0, .005849811d0,2.915706834d0,&
               .728926709d0 /
    do i = 1, lmax*lmax
       yl(i) = 0.d0
       yl(i) = 0.d0
    enddo
    yl(1) = 0.282094792d0
!
      if (lmax.eq.0) return
      a = 1.d0 - z * z
      if (abs(a) .lt.1.d-5) a = 0.d0
!
!     the following statement accounts for the phase convention
!     adopted in ylm1, e.g. px ~ -x, py ~ -y, pz ~ +z.
!
      x = - sqrt(a)
      zz = z * z
      xx = x * x
      sinphi = sin(phi)
      cosphi = cos(phi)
      yl(3) = fac(2) * z
      y = fac(3) * x
      yl(4) = y * cosphi
      yl(2) = y * sinphi
!
      if (lmax.eq.1) return
      sn2phi = sin(2.d0 * phi)
      cs2phi = cos(2.d0 * phi)
      yl(7) = fac(4) * (3.d0 * zz - 1.d0)
      y = fac(5) * z * x
      yl(8) = y * cosphi
      yl(6) = y * sinphi
      y = fac(6) * xx
      yl(9) = y * cs2phi
      yl(5) = y * sn2phi
!
      if (lmax.eq.2) return
      sn3phi = sin(3.d0 * phi)
      cs3phi = cos(3.d0 * phi)
      yl(13) = fac(7) * (5.d0 * zz * z - 3.d0 * z)
      y = fac(8) * x * (5.d0 * zz - 1.d0)
      yl(14) = y * cosphi
      yl(12) = y * sinphi
      y = fac(9) * z * xx
      yl(15) = y * cs2phi
      yl(11) = y * sn2phi
      y = fac(10) * x * xx
      yl(16) = y * cs3phi
      yl(10) = y * sn3phi
!
      if (lmax.eq.3) return
      sn4phi = sin(4.d0 * phi)
      cs4phi = cos(4.d0 * phi)
      yl(21) = fac(11) * (35.d0 * zz * zz - 30.d0 * zz + 3.d0)
      y = fac(12) * x * (7.d0 * zz * z - 3.d0 * z)
      yl(22) = y * cosphi
      yl(20) = y * sinphi
      y = fac(13) * xx * (7.d0 * zz - 1.d0)
      yl(23) = y * cs2phi
      yl(19) = y * sn2phi
      y = fac(14) * xx * x * z
      yl(24) = y * cs3phi
      yl(18) = y * sn3phi
      y = fac(15) * xx * xx
      yl(25) = y * cs4phi
      yl(17) = y * sn4phi
!
      if (lmax.eq.4) return
      sn5phi = sin(5.d0 * phi)
      cs5phi = cos(5.d0 * phi)
      z3 = zz * z
      z4 = zz * zz
      x3 = xx * x
      x4 = xx * xx
      yl(31) = fac(16) * (63.d0 * z4 * z - 70.d0 * z3 + 15.d0 * z)
      y = fac(17) * x * (21.d0 * z4 - 14.d0 * zz + 1.d0)
      yl(32) = y * cosphi
      yl(30) = y * sinphi
      y = fac(18) * xx * (3.d0 * z3 - z)
      yl(33) = y * cs2phi
      yl(29) = y * sn2phi
      y = fac(19) * x3 * (9.d0 * zz - 1.d0)
      yl(34) = y * cs3phi
      yl(28) = y * sn3phi
      y = fac(20) * x4 * z
      yl(35) = y * cs4phi
      yl(27) = y * sn4phi
      y = fac(21) * x4 * x
      yl(36) = y * cs5phi
      yl(26) = y * sn5phi
!
      if (lmax.eq.5) return
      sn6phi = sin(6.d0 * phi)
      cs6phi = cos(6.d0 * phi)
      yl(43) = fac(22) * (231.d0 * z4 * zz - 315.d0 * z4 + 105.d0 *zz - 5.d0)
      y = fac(23) * x * (33.d0 * z4 * z - 30.d0 * z3 + 5.d0 * z)
      yl(44) = y * cosphi
      yl(42) = y * sinphi
      y = fac(24) * xx * (33.d0 * z4 - 18.d0 * zz + 1.d0)
      yl(45) = y * cs2phi
      yl(41) = y * sn2phi
      y = fac(25) * x3 * (11.d0 * z3 - 3.d0 * z)
      yl(46) = y * cs3phi
      yl(40) = y * sn3phi
      y = fac(26) * x4 * (11.d0 * zz - 1.d0)
      yl(47) = y * cs4phi
      yl(39) = y * sn4phi
      y = fac(27) * x4 * x * z
      yl(48) = y * cs5phi
      yl(38) = y * sn5phi
      y = fac(28) * x4 * xx
      yl(49) = y * cs6phi
      yl(37) = y * sn6phi
!
      if (lmax.eq.6) return
      sn7phi = sin(7.d0 * phi)
      cs7phi = cos(7.d0 * phi)
      yl(57) = fac(29) * (429.d0 * z4 * z3 - 693.d0 * z4 * z + 315.d0* z3 - 35.d0 * z)
      y = fac(30) * x * (429.d0 * z4 * zz - 495.d0 * z4 + 135.d0 * zz -5.d0)
      yl(58) = y * cosphi
      yl(56) = y * sinphi
      y = fac(31) * xx * (143.d0 * z3 * zz - 110.d0 * z3 + 15.d0 * z)
      yl(59) = y * cs2phi
      yl(55) = y * sn2phi
      y = fac(32) * x3 * (143.d0 * z4 - 66.d0 * zz + 3.d0)
      yl(60) = y * cs3phi
      yl(54) = y * sn3phi
      y = fac(33) * x4 * (13.d0 * z3 - 3.d0 * z)
      yl(61) = y * cs4phi
      yl(53) = y * sn4phi
      y = fac(34) * x4 * x * (13.d0 * zz - 1.d0)
      yl(62) = y * cs5phi
      yl(52) = y * sn5phi
      y = fac(35) * x3 * x3 * z
      yl(63) = y * cs6phi
      yl(51) = y * sn6phi
      y = fac(36) * x3 * x4
      yl(64) = y * cs7phi
      yl(50) = y * sn7phi
!
      if (lmax.eq.7) return
      sn8phi = sin(8.d0 * phi)
      cs8phi = cos(8.d0 * phi)
      yl(73) = fac(37) * (6435.d0 * z4 * z4 - 12012.d0 * z4 * zz +6930.d0 * z4 - 1260.d0 * zz + 35.d0)
      y = fac(38) * x * (6435.d0 * z4 * z3 - 9009.d0 * z3 * zz +3465.d0 * z3 - 315.d0 * z)
      yl(74) = y * cosphi
      yl(72) = y * sinphi
      y = fac(39) * xx * (1001.d0 * z4 * zz - 1001.d0 * z4 + 231.d0 *zz - 7.d0)
      yl(75) = y * cs2phi
      yl(71) = y * sn2phi
      y = fac(40) * x3 * (273.d0 * z3 * zz - 182.d0 * z3 + 21.d0 * z)
      yl(76) = y * cs3phi
      yl(70) = y * sn3phi
      y = fac(41) * x4 * (455.d0 * z4 - 182.d0 * zz + 7.d0)
      yl(77) = y * cs4phi
      yl(69) = y * sn4phi
      y = fac(42) * x3 * xx * (455.d0 * z3 - 91.d0 * z)
      yl(78) = y * cs5phi
      yl(68) = y * sn5phi
      y = fac(43) * x3 * x3 * (1365.d0 * zz - 91.d0)
      yl(79) = y * cs6phi
      yl(67) = y * sn6phi
      y = fac(44) * x3 * x4 * z
      yl(80) = y * cs7phi
      yl(66) = y * sn7phi
      y = fac(45) * x4 * x4
      yl(81) = y * cs8phi
      yl(65) = y * sn8phi
!
      if (lmax.eq.8) return
      print 1000
      stop
 1000 format(' ********* program stop -- lmax > 8 in ylmw **********')
      end subroutine math_ylmw2


!
!-----
!H
!H      function math_lm(l,m)
!H
!###
  integer(FINT) function math_lm(l, m)
!
!
    integer(FINT),intent(in) :: l
    integer(FINT),intent(in) :: m
!
    math_lm = l *(l + 1) + m + 1
!
    return
  end function math_lm
!
!
!
!H    ******************************************************************
!H    *              SUBROUTINE MATH_INVIDX(IDMN,INDX,IRANK)                *
!H    ******************************************************************
!H
!
  subroutine math_invidx(indx, irank)
!
    integer(FINT),intent(in), dimension(:) :: indx
    integer(FINT),intent(out), dimension(:) :: irank
!
    integer(FINT) :: j
!
    do j = 1, size(indx)
       irank(indx(j) ) = j
    enddo
!
    return
  end subroutine math_invidx


  recursive subroutine quicksort_r(list, order)

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.
    ! Included in Utilities by F.Mariotti.

    implicit none

    real(FREAL), dimension (:), intent(INOUT)  :: list
    integer(FINT), dimension (:), intent(OUT) :: order

    ! Local variable
    integer(FINT) :: i

    do i = 1, size(list)
       order(i) = i
    end do

    call quick_sort_1(1, size(list))

  contains

    recursive subroutine quick_sort_1(left_end, right_end)

      integer(FINT), intent(IN) :: left_end, right_end

      !     Local variables
      integer(FINT)             :: i, j, itemp
      real(FREAL)               :: reference, temp
      integer(FINT), parameter  :: max_simple_sort_size = 6

      if (right_end < left_end + max_simple_sort_size) then
         ! Use interchange sort for small lists
         call interchange_sort(left_end, right_end)

      else
         ! Use partition ("quick") sort
         itemp = (left_end + right_end)/2 !Correct a compiler bug
         reference = list(itemp)
         i = left_end - 1; j = right_end + 1

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


            if (i < j) then
               ! Swap two out-of-order elements
               temp = list(i); list(i) = list(j); list(j) = temp
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
      real(FREAL)               :: temp

      do i = left_end, right_end - 1
         do j = i+1, right_end
            if (list(i) > list(j)) then
               temp = list(i); list(i) = list(j); list(j) = temp
               itemp = order(i); order(i) = order(j); order(j) = itemp
            end if
         end do
      end do

    end subroutine interchange_sort

  end subroutine quicksort_r

  recursive subroutine quicksort_i(list, order)

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.
    ! Included in Utilities by F.Mariotti.

    implicit none

    integer(FINT), dimension (:), intent(INOUT)  :: list
    integer(FINT), dimension (:), intent(OUT) :: order

    ! Local variable
    integer(FINT) :: i

    do i = 1, size(list)
       order(i) = i
    end do

    call quick_sort_1(1, size(list))

  contains

    recursive subroutine quick_sort_1(left_end, right_end)

      integer(FINT), intent(IN) :: left_end, right_end

      !     Local variables
      integer(FINT)             :: reference, temp
      integer(FINT)             :: i, j, itemp
      integer(FINT), parameter  :: max_simple_sort_size = 6

      if (right_end < left_end + max_simple_sort_size) then
         ! Use interchange sort for small lists
         call interchange_sort(left_end, right_end)

      else
         ! Use partition ("quick") sort
         reference = (left_end + right_end)/2
         reference = list(reference)
         i = left_end - 1; j = right_end + 1

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


            if (i < j) then
               ! Swap two out-of-order elements
               temp = list(i); list(i) = list(j); list(j) = temp
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
      integer(FINT)             :: i, j, itemp, temp

      do i = left_end, right_end - 1
         do j = i+1, right_end
            if (list(i) > list(j)) then
               temp = list(i); list(i) = list(j); list(j) = temp
               itemp = order(i); order(i) = order(j); order(j) = itemp
            end if
         end do
      end do

    end subroutine interchange_sort

  end subroutine quicksort_i

  recursive subroutine quicksortidx_i(list, order)

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.
    !
    ! Included in Utilities and modified in order to do not sort list
    ! but just to return the index (by F.Mariotti).

    implicit none

    integer(FINT), dimension (:), intent(INOUT)  :: list
    integer(FINT), dimension (:), intent(OUT) :: order

    ! Local variable
    integer(FINT) :: i

    do i = 1, size(list)
       order(i) = i
    end do

    call quick_sort_1(1, size(list))

  contains

    recursive subroutine quick_sort_1(left_end, right_end)

      integer(FINT), intent(IN) :: left_end, right_end

      !     Local variables
      integer(FINT)             :: reference, temp
      integer(FINT)             :: i, j, itemp
      integer(FINT), parameter  :: max_simple_sort_size = 6

      if (right_end < left_end + max_simple_sort_size) then
         ! Use interchange sort for small lists
         call interchange_sort(left_end, right_end)

      else
         ! Use partition ("quick") sort
         itemp = (left_end + right_end)/2
         itemp = order(itemp)
         reference = list(itemp)
         i = left_end - 1; j = right_end + 1

         do
            ! Scan list from left end until element >= reference is found
            do
               i = i + 1
               if (list(order(i)) >= reference) exit
            end do
            ! Scan list from right end until element <= reference is found
            do
               j = j - 1
               if (list(order(j)) <= reference) exit
            end do


            if (i < j) then
               ! Swap two out-of-order elements
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
      integer(FINT)             :: i, j, itemp, temp

      do i = left_end, right_end - 1
         do j = i+1, right_end
            if (list(order(i)) > list(order(j))) then
               itemp = order(i); order(i) = order(j); order(j) = itemp
            end if
         end do
      end do

    end subroutine interchange_sort

  end subroutine quicksortidx_i
!
! Quick Sort Index Intengers with dimension
!
  recursive subroutine quicksortidx_ni(idmn, list, order)

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.
    !
    ! Included in Utilities and modified in order to do not sort list
    ! but just to return the index (by F.Mariotti).

    implicit none

    integer(FINT), intent(in) :: idmn
    integer(FINT), dimension (:), intent(INOUT)  :: list
    integer(FINT), dimension (:), intent(OUT) :: order

    ! Local variable
    integer(FINT) :: i
!
    if(idmn.gt.size(list)) then
       stop ' Internal error quicksortidx_ni'
    endif
    do i = 1, idmn
       order(i) = i
    end do

    call quick_sort_1(1, idmn)

  contains

    recursive subroutine quick_sort_1(left_end, right_end)

      integer(FINT), intent(IN) :: left_end, right_end

      !     Local variables
      integer(FINT)             :: reference, temp
      integer(FINT)             :: i, j, itemp
      integer(FINT), parameter  :: max_simple_sort_size = 6

      if (right_end < left_end + max_simple_sort_size) then
         ! Use interchange sort for small lists
         call interchange_sort(left_end, right_end)

      else
         ! Use partition ("quick") sort
         itemp = (left_end + right_end)/2
         itemp = order(itemp)
         reference = list(itemp)
         i = left_end - 1; j = right_end + 1

         do
            ! Scan list from left end until element >= reference is found
            do
               i = i + 1
               if (list(order(i)) >= reference) exit
            end do
            ! Scan list from right end until element <= reference is found
            do
               j = j - 1
               if (list(order(j)) <= reference) exit
            end do


            if (i < j) then
               ! Swap two out-of-order elements
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
      integer(FINT)             :: i, j, itemp, temp

      do i = left_end, right_end - 1
         do j = i+1, right_end
            if (list(order(i)) > list(order(j))) then
               itemp = order(i); order(i) = order(j); order(j) = itemp
            end if
         end do
      end do

    end subroutine interchange_sort

  end subroutine quicksortidx_ni
!
!
!
  recursive subroutine quicksortidx_r(list, order)

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.
    !
    ! Included in Utilities and modified in order to do not sort list
    ! but just to return the index (by F.Mariotti).

    implicit none

    real(FREAL), dimension (:), intent(INOUT)  :: list
    integer(FINT), dimension (:), intent(OUT) :: order

    ! Local variable
    integer(FINT) :: i

    do i = 1, size(list)
       order(i) = i
    end do

    call quick_sort_1(1, size(list))

  contains

    recursive subroutine quick_sort_1(left_end, right_end)

      integer(FINT), intent(IN) :: left_end, right_end

      !     Local variables
      real(FREAL)               :: reference
      integer(FINT)             :: i, j, itemp
      integer(FINT), parameter  :: max_simple_sort_size = 6

      if (right_end < left_end + max_simple_sort_size) then
         ! Use interchange sort for small lists
         call interchange_sort(left_end, right_end)

      else
         ! Use partition ("quick") sort
         itemp = (left_end + right_end)/2
         itemp = order(itemp) !F90 Sun compiler2.0 bug .. that's why
         reference = list(itemp)
         i = left_end - 1; j = right_end + 1

         do
            ! Scan list from left end until element >= reference is found
            do
               i = i + 1
               itemp = order(i)
               if (list(itemp) >= reference) exit
            end do
            ! Scan list from right end until element <= reference is found
            do
               j = j - 1
               itemp = order(j)
               if (list(itemp) <= reference) exit
            end do


            if (i < j) then
               ! Swap two out-of-order elements
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

      do i = left_end, right_end - 1
         do j = i+1, right_end
            if (list(order(i)) > list(order(j))) then
               itemp = order(i); order(i) = order(j); order(j) = itemp
            end if
         end do
      end do

    end subroutine interchange_sort

  end subroutine quicksortidx_r
!
! Quick Sort Index IntengersReals with dimension
!
  recursive subroutine quicksortidx_nr(idmn, list, order)

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.
    !
    ! Included in Utilities and modified in order to do not sort list
    ! but just to return the index (by F.Mariotti).

    implicit none

    integer(FINT), intent(in) :: idmn
    real(FREAL), dimension (:), intent(INOUT)  :: list
    integer(FINT), dimension (:), intent(OUT) :: order

    ! Local variable
    integer(FINT) :: i
!
    if(idmn.gt.size(list)) then
       stop ' Internal error quicksortidx_ni'
    endif
!
    do i = 1, size(list)
       order(i) = i
    end do

    call quick_sort_1(1, size(list))

  contains

    recursive subroutine quick_sort_1(left_end, right_end)

      integer(FINT), intent(IN) :: left_end, right_end

      !     Local variables
      real(FREAL)               :: reference
      integer(FINT)             :: i, j, itemp
      integer(FINT), parameter  :: max_simple_sort_size = 6

      if (right_end < left_end + max_simple_sort_size) then
         ! Use interchange sort for small lists
         call interchange_sort(left_end, right_end)

      else
         ! Use partition ("quick") sort
         itemp = (left_end + right_end)/2
         itemp = order(itemp) !F90 Sun compiler2.0 bug .. that's why
         reference = list(itemp)
         i = left_end - 1; j = right_end + 1

         do
            ! Scan list from left end until element >= reference is found
            do
               i = i + 1
               itemp = order(i)
               if (list(itemp) >= reference) exit
            end do
            ! Scan list from right end until element <= reference is found
            do
               j = j - 1
               itemp = order(j)
               if (list(itemp) <= reference) exit
            end do


            if (i < j) then
               ! Swap two out-of-order elements
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

      do i = left_end, right_end - 1
         do j = i+1, right_end
            if (list(order(i)) > list(order(j))) then
               itemp = order(i); order(i) = order(j); order(j) = itemp
            end if
         end do
      end do

    end subroutine interchange_sort

  end subroutine quicksortidx_nr
!
!
!
  subroutine math_permut(neqs,s,rmat,naperm,invrnt,leps,match)
    !
    !***********************************************************************
    !
    !     permutes the 'neqs' members of a symmetry-equivalent set of
    !     atoms by applying rotation matrix 'rmat' to the atomic
    !     coordinate vectors 's', producing rotated vectors 't'.
    !     'naperm(na)' is returned as the atom into which na is permuted
    !     by the symmetry operation.( passive rotation, possibly
    !     including inversion ).
    !
    !     if the set of neqs vectors is not closed under the symmetry
    !     operation, return invrnt = .false.
    !
    !     called by: addop
    !
    !***********************************************************************
    !
    use vartypes
    implicit none
    !integer, parameter :: KINT  = kind(1)
    !integer, parameter :: FREAL = kind(1.0d0)
    !integer, parameter :: LCHARS = 160


    !
    logical :: invrnt
    logical, dimension(:) :: match
    real(FREAL) :: leps
    !
    integer(FINT), dimension(:) :: naperm
    real(FREAL), dimension(3) :: t
    real(FREAL), dimension(3,3) :: rmat
    real(FREAL), dimension(:,:) :: s
    !
    !     FROM INLINE CODE
    !
    !integer(FINT) :: i
    !integer(FINT) :: j
    !integer(FINT) :: na
    logical :: notdone
    integer(FINT) :: na1
    integer(FINT) :: na2
    integer(FINT) :: neqs
    real(FREAL) :: rab
    !real(FREAL) :: xd
    !real(FREAL) :: yd
    !real(FREAL) :: zd
    !
    invrnt = .false.
    notdone = .true.
    !
    match = .false.  !! Vector
    !


!XXX diff version....

    do na1=1,neqs
       t = matmul(rmat,s(1:3,na1))
       invrnt = .false.
       do na2=1,neqs
          rab = sqrt(sum((t-s(1:3,na2))*(t-s(1:3,na2))))
          if(rab.lt.leps) then
             invrnt = .true.
             exit
          endif
       enddo
       if(.not.invrnt) then
          exit
       endif
    enddo

    return

!    do na1 = 1,neqs
    na1 = 1
    do while(notdone)
       !
       !--------rotate coord vector of atom na1 with rmat ---------------------
       t = matmul(rmat,s(1:3,na1))
       !
       !--------test rotated vector for congruency with one of the s vectors --
       do na2 = 1,neqs
          if(.not.match(na2)) then
             rab = sqrt(sum((t-s(1:3,na2))*(t-s(1:3,na2))))
             if(rab.lt.leps) then
                invrnt = .true.
                match(na2) = .true.
                naperm(na1) = na2
                exit
             endif
          endif
       enddo
       !
       na1 = na1 + 1
       if(na1.gt.neqs)then
          exit
       endif
       !if(.not.invrnt) exit
       invrnt = .false.
    enddo
    !
    !
    !
    return
  end subroutine math_permut
!
!H------------------------------------------------------------------------------
!H Function math_getorigin
!H------------------------------------------------------------------------------
!H
!
  integer(FINT) function math_getorigin(xyz, nat, leps)
!
!H
!H    Find the orgin in the xyz vector
!H     n   return the index of the origin if any
!H     0   if no origin
!H    -1   on error
!H    -2   if more then one origin
!H
!H    Version with all loop on atoms (i.e. check more then one origin)
!H
!
    integer(FINT), intent(in) :: nat
    real(FREAL), intent(in) :: leps
    real(FREAL), intent(in), dimension(:,:) :: xyz
    !
    integer(FINT) :: i
    real(FREAL) :: r
    !
    math_getorigin = 0
    i = 0

    do i = 1, nat
       r=sqrt(xyz(1,i)*xyz(1,i)+xyz(2,i)*xyz(2,i)+xyz(3,i)*xyz(3,i))
       if (r.lt.leps) then
          if (math_getorigin.ne.0) then
             math_getorigin = - 2
          endif
          math_getorigin = i
       endif
    enddo
    !
    return
  end function math_getorigin
  !
  !
  !
  !
  !H--------------------------------------------------------------------------
  !H Function math_factorials(fact,idmn)
  !H--------------------------------------------------------------------------
  !H
  !
  subroutine math_factorials(fact,idmn)
    !
    implicit none
    !
    !H  
    !H  Compute factorials
    !H  
    !H  fact    pointer to a vector containing factorials.
    !H          Factorials are stored as N! = fact(N+1) i.e. 0!=fact(1)
    !H  
    !H  idmn    The biggest factorial required. If idmn is less
    !H          then the previous computed vector the factorials
    !H           are not recoputed.
    !H
    !H          idmn is the vector dimension: to store N! idmn=N+1
    !H          idmn.lt.0 require deallocation
    !H  
    !H  
    !
    real(FREAL), dimension(:), pointer :: fact
    integer(FINT), intent(in) :: idmn
    !
    logical, save :: iscomputed = .false.
    integer(FINT), save :: max_fact
    real(FREAL), dimension(:), save, allocatable, target :: factorials
    !
    integer(FINT) :: irc
    integer(FINT) :: i
    !
    ! Check for deallocate request!
    !
    if(idmn.lt.0) then !idmn < 0 means deallocate
       if(allocated(factorials)) then
          deallocate(factorials,STAT=irc)
          iscomputed = .false. !set it not computed
       endif
       nullify(fact) !We nullify it anyway
       return
    endif
    !
    ! Catch Errors
    !
    ! our factorials is N! = fact(N+1) ie N=1 is 0!
    if(idmn.eq.0) then
       mathtools_error = -1
       return
    endif
    !
    ! Compute factorials
    !
    if(.not.iscomputed) then
       if(allocated(factorials)) then
          deallocate(factorials,STAT=irc)
          nullify(fact)
          if(irc.ne.0) then
             mathtools_error = -1
             return
          endif
       endif
    else
       ! it has been computed
       if(idmn.le.size(factorials)) then
          return !we had stored it already
       endif
       !
       deallocate(factorials,STAT=irc)
       if(irc.ne.0) then
          mathtools_error = -1
          return
       endif
       nullify(fact)
    endif
    !
    ! it has nor been computed or computed at lower dim= do it again
    allocate(factorials(idmn),STAT=irc)
    if(irc.ne.0) then
       mathtools_error = -1
       return
    endif
    factorials(1) = 1.0_FREAL
    !
    ! Compute factorials
    !
    do i=2,idmn
       factorials(i) = real(i-1,KIND=FREAL) * factorials(i-1)
    enddo
    !
    ! Associate the pointer
    !
    fact => factorials(1:idmn)
    !
  end subroutine math_factorials
  !
  !
  !
  recursive function math_getfact(fact) result(rfact)
    !
    real(FREAL), intent(in) :: fact
    real(FREAL) :: rfact
    !
    if(fact.gt.eps)then
       rfact = fact * math_getfact(fact-1.0_FREAL)
    else
       rfact = 1.0_FREAL
    endif
    !
  end function math_getfact
  !
  !
  !
  !
  subroutine math_bincoe(i,j,b)
    !
    !     this subroutine computes binomial coefficients
    !
    integer(FINT), intent(in) :: i
    integer(FINT), intent(in) :: j
    real(FREAL), intent(out) :: b
    !
    integer(FINT) :: k
    !real(FREAL), dimension(:), pointer :: fact
    real(FREAL), dimension(60) :: fact
    !
    fact(1)=1.d0
    do k=2,60
       fact(k)=fact(k-1)*real(k-1,KIND=FREAL)
    enddo
    !
    k = max(i,j) + 1
    !call math_factorials(fact,k)
    !
    k=j-i
    !write (0,*) 'B',i,j,k
    !call flush(0)
    b=fact(j+1)/(fact(i+1)/fact(k+1))
    !
  end subroutine math_bincoe
  !
  !
  !
end module adfmathtools
!
!
!
