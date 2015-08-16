!
!H
!H-----------------------------------------------------------------------------
!H This File is part of the Freemol package
!H (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
!H Please read the full copyright statment in the COPYING file
!H-----------------------------------------------------------------------------
!H MODULE mathtools  Freemol by F.Mariotti: (c) F.Mariotti
!H-----------------------------------------------------------------------------
!H $Id: mathtools.F90,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
!H-----------------------------------------------------------------------------
!H 
!H DESCRIPTION: We store here maths routines with a genral meaning.
!H              This module will split up in different other modules
!H              as needed but there will be one point when it will get
!H              a better defined picture.
!H
!H mathtools_XXX is the public name
!H mt_XXX is the private one.
!H 
!H 
!H TODO: 
!H 
!H 
!H WARN: 
!H 
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H 
!
module mathtools
  !
  use vartypes
  use messages
  !use modlapack
  !
  implicit none
  private
  !
  ! PUBLIC DECLARATIONS:
  !--------------------
  
  ! Public Procedures
  !------------------
  public mathtools_init
  public mathtools_isinit
  public mathtools_geterror
  !
  public mathtools_masscenter
  public mathtools_inertia
  public mathtools_euler
  !public mathtools_
  !public mathtools_
  !
  ! Interfaces
  !-----------
  interface mathtools_masscenter
     module procedure mathtools_masscenter_r
     module procedure mathtools_masscenter_i
  end interface

  !
  ! Error Interface
  !----------------
  ! We define two error interfaces. Each one return the latest error
  ! code or string value. The String error function can be used to inqure
  ! about a given error code.
  !----------------------------------------------------------------------
  interface mathtools_geterror
     module procedure mathtools_geterror_i
     module procedure mathtools_geterror_c
  end interface !mathtools_geterror

  ! LOCAL DECLARATIONS:
  !--------------------
  !
  !
  ! error variables
  !----------------
  integer(FINT), save, private :: mt_error
  character(FLCHARS), save, private :: mt_error_message
  logical :: mt_isinit = .false.
  !
contains
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H Function mathtools_init()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function mathtools_init()
!
!H
!H-----------------------------------------------------------------------------
!H
!H This routine initialize the mathtools module.
!H
!H-----------------------------------------------------------------------------
!H
!
    !
    ! WE do not need yet initialization so
    !-------------------------------------
    mt_error = 0
    mt_error_message = ""
    mt_isinit = .true.
    mathtools_init = 0
    !
    return 
  end function mathtools_init
!
!H
!H-----------------------------------------------------------------------------
!H logical function mathtools_isinit()
!H-----------------------------------------------------------------------------
!H
!
  logical function mathtools_isinit()
!
!H
!H-----------------------------------------------------------------------------
!H We merelyreturn the initialization flag.
!H-----------------------------------------------------------------------------
!H
!
    !
    mathtools_isinit = mt_isinit
    !
  end function mathtools_isinit
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H integer(FINT) function mathtools_geterror_i()
!H-----------------------------------------------------------------------------
!H
!
  integer(FINT) function mathtools_geterror_i()
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
    mathtools_geterror_i = mt_error
    !
  end function mathtools_geterror_i
!
!H-----------------------------------------------------------------------------
!H
!
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H character(FLCHARS) function mathtools_geterror_c(core)
!H-----------------------------------------------------------------------------
!H
!
  character(FLCHARS) function mathtools_geterror_c(code)
    !
    integer(FINT), intent(in) :: code
    !integer(FINT), optional, intent(in) :: code
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
    !integer(FINT) :: lcode !the local used error code
    !
    !lcode = mt_error
    !if(present(code)) then
    !   lcode = code
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
       mathtools_geterror_c = "No Errors."
    case(1_FINT)
       mathtools_geterror_c = "Not documented error or generic error.."
    case(2_FINT)
       mathtools_geterror_c = "Double Error: Why did you get this message?."
    case default
       mathtools_geterror_c = "Not documented error or generic error.."
    end select
    !
    return
    !
  end function mathtools_geterror_c
!
!H
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H START
!H-----------------------------------------------------------------------------
!H-----------------------------------------------------------------------------
!H
!
!H
!H-----------------------------------------------------------------------------
!H subroutine mathtools_masscenter_r(xyz,we,rvec,doit)
!H-----------------------------------------------------------------------------
!H
!
  subroutine mathtools_masscenter_r(xyz,we,rvec,doit)
!
!H-----------------------------------------------------------------------------
!H  
!H  xyz       xyz cartesian
!H  we        Weight of each center !! you should give it!!!
!H  rvec      Position of the center of mass
!H  doit      Optinally do the translation!
!H  
!H-----------------------------------------------------------------------------
!
    !
    real(FREAL), dimension(:,:), intent(inout) :: xyz
    real(FREAL), dimension(:), intent(in) :: we
    real(FREAL), dimension(:), intent(out) :: rvec
    !
    integer(FINT) :: i
    !
    logical, optional :: doit
    !
    real(FREAL) :: wesum
    !
    ! Check compatibility
    !--------------------
    ! We assume first index in xyz as the dimension of the space
    ! and the second index the coordinates: in 3D space xyz
    ! is xyz(1:3,particles). Which means to have rvec as much as rvec(1:3).
    !----------------------------------------------------------------------
    if(size(rvec,1).ne.size(xyz,1)) then
       mt_error = 1
       mt_error_message = "Conformation error in: mass_center"
       call message(MESERRO,"[mathtools_masscenter]: Conformation error:1.")
       call message(MESDEBG,"[mathtools_masscenter]: Conformation error:1.")
       stop 1
    end if
    ! We check now the other side ...
    !--------------------------------
    if(size(we,1).ne.size(xyz,2)) then
       mt_error = 1
       mt_error_message = "Conformation error in: mass_center"
       write (*,*) 's1', size(we,1)
       write (*,*) 's2', size(xyz,2)
       call message(MESERRO,"[mathtools_masscenter]: Conformation error:2.")
       call message(MESDEBG,"[mathtools_masscenter]: Conformation error:2.")
       stop 1
    end if
    rvec(:) = 0.0_FREAL
    !
    ! We get sum weithed stuff! there should be abetter way in f90
    ! but at present I do not care!
    !-------------------------------------------------------------
    do i=1,size(rvec,1)
       rvec(i) = sum(we(:)*xyz(i,:))
    end do
    !
    ! Vectorial normalization
    !------------------------
    wesum=sum(we(:))
    rvec(:) = rvec(:) / wesum
    !
    !We remove zeros errors
    !----------------------
    do i=1,size(rvec,1)
       if(abs(rvec(i)).lt.F_EPS) rvec(i) = 0.0_FREAL
    end do
    !
    if(present(doit)) then
       if (doit) then
          do i=1,size(rvec,1)
             xyz(i,:) = xyz(i,:) - rvec(i)
          end do
          !call message (MESOUT,"[mathtools_masscenter]: performed center of mass translation.")
       endif
    endif


  end subroutine mathtools_masscenter_r
!
!H
!H-----------------------------------------------------------------------------
!H subroutine mathtools_masscenter_r(xyz,we,rvec,doit)
!H-----------------------------------------------------------------------------
!H
!
  subroutine mathtools_masscenter_i(xyz,we,rvec,doit)
!
!H-----------------------------------------------------------------------------
!H  
!H  xyz       xyz cartesian
!H  we        Weight of each center !! you should give it!!!
!H  rvec      Position of the center of mass
!H  doit      Optinally do the translation!
!H  
!H-----------------------------------------------------------------------------
!
    !
    real(FREAL), dimension(:,:), intent(inout) :: xyz
    integer(FINT), dimension(:), intent(in) :: we
    real(FREAL), dimension(:), intent(out) :: rvec
    !
    integer(FINT) :: i
    !
    logical, optional :: doit
    !
    ! Check compatibility
    !--------------------
    ! We assume first index in xyz as the dimension of the space
    ! and the second index the coordinates: in 3D space xyz
    ! is xyz(1:3,particles). Which means to have rvec as much as rvec(1:3).
    !----------------------------------------------------------------------
    if(size(rvec,1).lt.size(xyz,1)) then
       mt_error = 1
       mt_error_message = "Conformation error in: mass_center"
       call message(MESDEBG,"[, intent(in)]: Conformation error.")
       stop 1
    end if
    ! We check now the other side ...
    !--------------------------------
    if(size(we,1).lt.size(xyz,2)) then
       mt_error = 1
       mt_error_message = "Conformation error in: mass_center"
       call message(MESDEBG,"[mathtools_masscenter]: Conformation error:2.")
       stop 1
    end if
    rvec(:) = 0.0_FREAL
    !
    ! We get sum weithed stuff! there should be abetter way in f90
    ! but at present I do not care!
    !-------------------------------------------------------------
    do i=1,size(rvec,1)
       rvec(i) = sum(we*xyz(i,:))
    end do
    !
    ! Vectorial normalization
    !------------------------
    rvec(:) = rvec(:) / sum(we(:))
    !
    !We remove zeros errors
    !----------------------
    do i=1,size(rvec,1)
       if(rvec(i).lt.F_EPS) rvec(i) = 0.0_FREAL
    end do
    !
    if(present(doit)) then
       if (doit) then
          do i=1,size(rvec,1)
             xyz(i,:) = xyz(i,:) - rvec(i)
          end do
          !call message (MESOUT,"[mathtools_masscenter]: performed center of mass translation.")
       endif
    endif
    !
    !
  end subroutine mathtools_masscenter_i
  !
  !
!
!H
!H-----------------------------------------------------------------------------
!H subroutine mathtools_inertia(xyz,we,rvec,doit)
!H-----------------------------------------------------------------------------
!H
!
  subroutine mathtools_inertia(xyz,we,rvec,doit)
!
!H-----------------------------------------------------------------------------
!H  
!H  xyz       xyz cartesian
!H  we        Weight of each center !! you should give it!!!
!H  rvec      Rotational Matrix
!H  doit      Optinally do the translation!
!H  
!H-----------------------------------------------------------------------------
!
    !
    real(FREAL), dimension(:,:), intent(inout) :: xyz
    real(FREAL), dimension(:), intent(in) :: we
    real(FREAL), dimension(:,:), intent(out) :: rvec
    !
    integer(FINT) :: i,num
    !
    logical, optional :: doit
    !
    real(FREAL), dimension(6) :: rimat
    real(FREAL), dimension(3,3) :: arimat, trvec
    real(FREAL), dimension(3) :: eigval
    real(FREAL), dimension(9) :: rwork
    real(FREAL) :: dsum
    !integer(FINT) :: info,ilow,iup !old for diagla
    integer(FINT) :: info
    logical, save :: dojakobi
    logical :: debug
    debug = .true.
    dojakobi = .false.
    !dojakobi = .true.
    !
    ! Check compatibility
    !--------------------
    ! We assume first index in xyz as the dimension of the space
    ! and the second index the coordinates: in 3D space xyz
    ! is xyz(1:3,particles). Which means to have rvec as much as rvec(1:3).
    !----------------------------------------------------------------------
    if(size(rvec,1).lt.size(xyz,1)) then
       mt_error = 1
       mt_error_message = "Conformation error in: mass_center"
       call message(MESDEBG,"[mathtools_masscenter]: Conformation error:1.")
       stop 1
    end if
    ! We check now the other side ...
    !--------------------------------
    if(size(we,1).lt.size(xyz,2)) then
       mt_error = 1
       mt_error_message = "Conformation error in: mass_center"
       call message(MESDEBG,"[mathtools_masscenter]: Conformation error:2.")
       stop 1
    end if
    ! Now we can do it only for 3D objects
    !-------------------------------------
    if(size(xyz,1).gt.3) then
       call message(MESDEBG,"[mathtools_masscenter]: We handle only 3D now")
       stop 3
    endif
    num = size(we)
    rvec(:,:) = 0.0_FREAL
    trvec(:,:) = 0.0_FREAL
    !
    !
    !dsum =0.0_FREAL
    !do i=1,num
    !   dsum = dsum + we(i) * (xyz(1,i) * xyz(3,i))
    !   write (*,*) ' Step:',i,' V ',dsum,' Val: ',we(i) * (xyz(1,i) * xyz(3,i)),'terms: ',we(i),(xyz(1,i),xyz(3,i))
    !end do
    !write (*,*) ' END: V',dsum
    !
    rimat(1) = - sum(we(1:num) * (xyz(1,1:num) * xyz(1,1:num)))
    rimat(2) = - sum(we(1:num) * (xyz(1,1:num) * xyz(2,1:num)))
    rimat(3) = - sum(we(1:num) * (xyz(2,1:num) * xyz(2,1:num)))
    rimat(4) = - sum(we(1:num) * (xyz(1,1:num) * xyz(3,1:num)))
    rimat(5) = - sum(we(1:num) * (xyz(2,1:num) * xyz(3,1:num)))
    rimat(6) = - sum(we(1:num) * (xyz(3,1:num) * xyz(3,1:num)))
    !
    dsum = rimat(1) + rimat(3) + rimat(6)
    rimat(1) = rimat(1) - dsum
    rimat(3) = rimat(3) - dsum
    rimat(6) = rimat(6) - dsum
    !
    arimat(1,1)=rimat(1)
    arimat(1,2)=rimat(2)
    arimat(1,3)=rimat(4)
    arimat(2,1)=rimat(2)
    arimat(2,2)=rimat(3)
    arimat(2,3)=rimat(5)
    arimat(3,1)=rimat(4)
    arimat(3,2)=rimat(5)
    arimat(3,3)=rimat(6)
    !
    if(debug) then
       write(6,'("# Inertia Tensor: (Before Diagonalization)")') 
       write(6,'("# Up Diagonal Columnwise")') 
       write(6,'("#",1X,6(1X,D16.8))') rimat
       write(6,'("# Inertia Tensor: (Full Matrix)")') 
       write(6,'("# 1,1 - 1,3      ",2X,3(1X,D16.8))') arimat(1,1:3)
       write(6,'("# 2,1 - 2,3      ",2X,3(1X,D16.8))') arimat(2,1:3)
       write(6,'("# 3,1 - 3,3      ",2X,3(1X,D16.8))') arimat(3,1:3)
    endif
    !
    !
    if (dojakobi) then
       arimat(1,1)=rimat(1)
       arimat(1,2)=rimat(2)
       arimat(1,3)=rimat(4)
       arimat(2,1)=rimat(2)
       arimat(2,2)=rimat(3)
       arimat(2,3)=rimat(5)
       arimat(3,1)=rimat(4)
       arimat(3,2)=rimat(5)
       arimat(3,3)=rimat(6)
       call jakobi(3,3,arimat(1:3,1:3),eigval(3),rvec(1:3,1:3),info)
       if (info.gt.0) info=0
       if(debug) then
          write(6,'("# (Jakobi) Inertia Tensor")') 
          write(6,'("# Inertia EigenValues:",2X,3(1X,D16.8))') eigval(1:3)
          write(6,'("# Eigen Vectors:")') 
          write(6,'("# 1,1 - 1,3      ",2X,3(1X,D16.8))') rvec(1,1:3)
          write(6,'("# 2,1 - 2,3      ",2X,3(1X,D16.8))') rvec(2,1:3)
          write(6,'("# 3,1 - 3,3      ",2X,3(1X,D16.8))') rvec(3,1:3)
          write(6,'("# USING Jakobi")') 
       endif
    else
       !
       rwork(:)=0.0_FREAL
       if (.false.) then
          call dspev('V','U',3,rimat(1:6),eigval(1:3),rvec(1:3,1:3),3,rwork(1:9),info)
       else
          call dsyev('V','U',3,arimat(1:3,1:3),3,eigval(1:3),rwork(1:9),9,info)
          rvec(:,:)=arimat(:,:)
       end if
       if(debug) then
          write(6,'("# Inertia Tensor")') 
          write(6,'("# Inertia EigenValues:",2X,3(1X,D16.8))') eigval(1:3)
          write(6,'("# Eigen Vectors:")') 
          write(6,'("# 1,1 - 1,3      ",2X,3(1X,D16.8))') rvec(1,1:3)
          write(6,'("# 2,1 - 2,3      ",2X,3(1X,D16.8))') rvec(2,1:3)
          write(6,'("# 3,1 - 3,3      ",2X,3(1X,D16.8))') rvec(3,1:3)
       endif
       !
    end if
    !
    if (info.ne.0) then
       call message(MESERRO,'[mathtools_inertia] Diag error in inertia.')
       return
    endif
    !
    !
    if(present(doit)) then
       if(doit) then
       trvec(1,1:3) = rvec(1:3,1)
       trvec(2,1:3) = rvec(1:3,2)
       trvec(3,1:3) = rvec(1:3,3)

          do i=1,num
             xyz(1:3,i) = matmul(trvec(1:3,1:3),xyz(1:3,i))
          enddo
       endif
    endif
    !
    !
    !
  end subroutine mathtools_inertia
!
!H
!H----------------------------------------------------------------------
!H subroutine mathtools_euler(angle,rotax,inv,alpha,beta,gamma,rmat,leps)
!H----------------------------------------------------------------------
!H
!
  subroutine mathtools_euler(angle, rotax, inv, alpha, beta, gamma, rmat,leps)
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
!    integer(FINT) :: i
!    integer(FINT) :: j
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
    real(FREAL) :: zero
    real(FREAL) :: one
    !
    !
    !
    pi = 4.0D0 * atan(1.0D0)
    twopi = 2.D0 * pi
    zero = 0.0_FREAL
    one = 1.0_FREAL
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
  end subroutine mathtools_euler
  !
  !
!
!H
!H----------------------------------------------------------------------
!H subroutine jakobi.. to test
!H----------------------------------------------------------------------
!H
!
      subroutine jakobi(N,M,A,U,V,NROT) 
!                                                                       
! EIGENVALUES AND EIGENVECTORS OF A real SYMMETRIC MATRIX A OF logical  
! SIZE A(N,N) STORED IN AN ARRAY OF PHYSICAL SIZE A(M,M), where         
! M .ge. N.                                                             
!                                                                       
! ON OUTPUT:                                                            
! (1) ELEMENTS OF A ABOVE THE DIAGONAL ARE DESTROYED.                   
! (2) THE ARRAY U(M) RETURNS THE EIGENVALUES OF A(N,N) IN ITS FIRST N   
!     ELEMENTS.                                                         
! (3) THE COLUMNS OF THE MATRIX V(N,N), STORED IN THE ARRAY V(M,M),     
!     CONTAIN THE EIGENVECTORS OF A.                                    
! (4) THE VARIABLE NROT RETURNS THE NUMBER OF ITERATIONS OF JACOBI      
!     ROTATION THAT WERE REQUIRED TO ANNIHILATE THE OFF-DIAGONAL        
!     ELEMENTS OF A(N,N) TO MACHINE precision.                          
!                                                                       
! FORTRAN CODE ADAPTED FROM WILLIAM H. PRESS, BRIAN P. FLANNERY, SAUL A.
! TEUKOLSKY, AND WILLIAM T. VETTERLING (1986).  NUMERICAL RECIPIES:  THE
!                                               --------- --------   ---
! ART OF SCIENTIFIC COMPUTING, PP. 335-349.  CAMBRIDGE, ENGLAND:        
! --- -- ---------- ---------                                           
! CAMBRIDGE UNIVERSITY PRESS.                                           
!                                                                       
        integer(FINT), parameter:: NMAX=999
        integer(FINT) :: n,m,nrot
        real(FREAL), dimension(1:M,1:M) :: A !(M,M)
        real(FREAL), dimension(1:M) :: U !(M)
        real(FREAL), dimension(1:M,1:M) :: V !(M,M)
        real(FREAL), dimension(1:NMAX) :: B !(NMAX)
        real(FREAL), dimension(1:NMAX) :: Z !(NMAX) 
        !double precision A,U,V,B,Z,T,THRESH,E,AII,AJJ,THETA,C,S,TAU,P,Q 
        real(FREAL) :: T,THRESH,E,AII,AJJ,THETA,C,S,TAU,P,Q 
        integer(FINT) :: i,j,k,ncycle
!                                                                       
! INITIALIZE EIGENVECTORS MATRIX TO AN IDENTITY MATRIX.                 
!                                                                       
      do I=1,N 
      do J=1,N 
        V(I,J)=0 
      end do 
        V(I,I)=1 
      end do 
!                                                                       
! INITIALIZE EIGENVALUES U AND WORK VECTOR B TO A(I,I) AND ZERO WORK    
! VECTOR Z.                                                             
!                                                                       
      do I=1,N 
        T=A(I,I) 
        U(I)=T 
        B(I)=T 
        Z(I)=0 
      end do 
!                                                                       
! PERFORM UP TO 50 ITERATIONS OF UP TO N*(N - 1)/2 JACOBI ROTATIONS.    
!                                                                       
      NROT=0 
      NCYCLE=1 
      do while (NCYCLE.le.50) 
!                                                                       
! TEST FOR NORMAL return WHEN MAXIMUM MAGNITUDE OF AN OFF-DIAGONAL      
! ELEMENT EQUALS ZERO TO MACHINE precision.  THE TEST PRESUMES THAT     
! ARITHMETIC UNDERFLOW VALUES ARE SET TO ZERO.                          
!                                                                       
        T=0 
        do I=1,N-1 
        do J=I+1,N 
          T=T+abs(A(I,J)) 
        end do 
        end do 
        if (T.eq.0) return 
!                                                                       
! SET OFF-DIAGONAL THRESHOLD.                                           
!                                                                       
        if (NCYCLE.lt.4) then 
          THRESH=T/(5*N**2) 
        else 
          THRESH=0 
        end if 
!                                                                       
! ROTATE TO ANNIHILATE OFF-DIAGONAL ELEMENTS.                           
!                                                                       
        do I=1,N-1 
        do J=I+1,N 
          T=abs(A(I,J)) 
          E=100*T 
          AII=abs(U(I)) 
          AJJ=abs(U(J)) 
          if (NCYCLE.gt.4.and.AII+E.eq.AII.and.AJJ+E.eq.AJJ) then 
!                                                                       
! AFTER FOUR CYCLES, SKIP THE ROTATION if abs(A(I,J)) IS SMALL COMPARED 
! TO BOTH abs(A(I,I)) AND abs(A(J,J)).                                  
!                                                                       
            A(I,J)=0 
          else if (T.gt.THRESH) then 
            T=abs(AJJ-AII) 
            if (T+E.eq.T) then 
!                                                                       
! EFFECTIVELY, SET T = 1/(2*THETA)                                      
!                                                                       
              T=A(I,J)/(AJJ-AII) 
            else 
              THETA=(AJJ-AII)/(2*A(I,J)) 
              T=1/(abs(THETA)+sqrt(1+THETA**2)) 
              if (THETA.lt.0) T=-T 
            end if 
            C=1/sqrt(1+T**2) 
            S=T*C 
            TAU=S/(1+C) 
!                                                                       
! ADJUST EIGENVALUES U AND WORK VECTOR Z.                               
!                                                                       
            E=T*A(I,J) 
            Z(I)=Z(I)-E 
            Z(J)=Z(J)+E 
            U(I)=U(I)-E 
            U(J)=U(J)+E 
            A(I,J)=0 
!                                                                       
! ROTATIONS  1 .le. K .lt. I.                                           
!                                                                       
            do K=1,I-1 
              P=A(K,I) 
              Q=A(K,J) 
              A(K,I)=P-S*(Q+P*TAU) 
              A(K,J)=Q+S*(P-Q*TAU) 
            end do 
!                                                                       
! ROTATIONS  I .lt. K .lt. J.                                           
!                                                                       
            do K=I+1,J-1 
              P=A(I,K) 
              Q=A(K,J) 
              A(I,K)=P-S*(Q+P*TAU) 
              A(K,J)=Q+S*(P-Q*TAU) 
            end do 
!                                                                       
! ROTATIONS  J .lt. K .le. N.                                           
!                                                                       
            do K=J+1,N 
              P=A(I,K) 
              Q=A(J,K) 
              A(I,K)=P-S*(Q+P*TAU) 
              A(J,K)=Q+S*(P-Q*TAU) 
            end do 
!                                                                       
! COMPUTE AND STORE EIGENVECTORS.                                       
!                                                                       
            do K=1,N 
              P=V(K,I) 
              Q=V(K,J) 
              V(K,I)=P-S*(Q+P*TAU) 
              V(K,J)=Q+S*(P-Q*TAU) 
            end do 
            NROT=NROT+1 
          end if 
        end do 
        end do 
!                                                                       
! ADJUST EIGENVALUES AND WORK VECTORS.                                  
!                                                                       
        do I=1,N 
          B(I)=B(I)+Z(I) 
          U(I)=B(I) 
          Z(I)=0 
        end do 
        NCYCLE=NCYCLE+1 
      end do 
      stop '50 CYCLES OF JACOBI ROTATIONS SHOULD NEVER BE NECESSARY.' 
    end subroutine jakobi

  

  !
  !
end module mathtools






