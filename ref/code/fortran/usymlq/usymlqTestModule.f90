!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File usymlqTestModule.f90
!
!    Hprod   Aprod1   Aprod2    lstp     usymlqtest
!
! These routines define a class of least-squares test problems
! for testing algorithm LSQR (Paige and Saunders, ACM TOMS, 1982).
! Aprod1 and Aprod2 define matrix-vector products required by
! subroutine USYMLQ for a test matrix of the form  A = Y*D*Z,
! where Y and Z are Householder transformations and D is diagonal.
!
! This file illustrates how USYMLQ can call Aprod1 and Aprod2 with a
! short fixed parameter list, even if they need arbitrary other data.
! Here, they need arrays d, hx, hy created by subroutine usymlqtest
! in this module before USYMLQ is called (see below).
!
! 25 May 2014: usymlqTestModule.f90 derived from lsqrTestModule.f90.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module usymlqTestModule

  use  usymlqDataModule,    only : ip, dp
  use  usymlqblasInterface, only : dnrm2, dscal
  use  usymlqModule,        only : USYMLQ
  use  usymlqCheckModule,   only : Acheck, xcheck

  implicit none
  private
  public      :: usymlqtest
  private     :: Hprod, Aprod1, Aprod2, lstp
  intrinsic   :: dot_product

  ! DYNAMIC WORKSPACE DEFINED HERE.
  ! They are allocated in usymlqtest and used by Aprod1, Aprod2.

  real(dp), allocatable :: d(:), hy(:), hz(:) ! These define A = Y D Z.
  real(dp), allocatable :: wm(:), wn(:)       ! Work vectors.

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Hprod (n,z,x)

    integer(ip), intent(in)    :: n
    real(dp),    intent(in)    :: z(n)
    real(dp),    intent(inout) :: x(n)

    !-------------------------------------------------------------------
    ! Hprod  applies a Householder transformation stored in z
    ! to return x = (I - 2*z*z')*x.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------

    real(dp) :: s

    s = dot_product( z,x )
    s = s + s
    x = x - s*z
  end subroutine Hprod

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Aprod1(m,n,x,y)

    integer(ip), intent(in)    :: m,n
    real(dp),    intent(in)    :: x(n)
    real(dp),    intent(inout) :: y(m)

    !-------------------------------------------------------------------
    ! Aprod1 computes y = y + A*x without altering x,
    ! where A is a test matrix of the form  A = Y*D*Z,
    ! and the matrices D, Y, Z are represented by
    ! the allocatable vectors d, hy, hz in this module.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------

    intrinsic           :: min
    integer(ip)         :: minmn
    real(dp), parameter :: zero = 0.0_dp

    minmn = min(m,n)
    wn    = x
    call Hprod (n,hz,wn)
    wm(1:minmn) = d(1:minmn)*wn(1:minmn)
    wm(n+1:m)   = zero
    call Hprod (m,hy,wm)
    y = y + wm

  end subroutine Aprod1

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine Aprod2(m,n,x,y)

    integer(ip), intent(in)    :: m,n
    real(dp),    intent(inout) :: x(n)
    real(dp),    intent(in)    :: y(m)

    !-------------------------------------------------------------------
    ! Aprod2 computes x = x + A'*y without altering y,
    ! where A is a test matrix of the form  A = Y*D*Z,
    ! and the matrices D, Y, Z are represented by
    ! the allocatable vectors d, hy, hz in this module.
    !
    ! 23 Sep 2007: Fortran 90 version.
    !-------------------------------------------------------------------

    intrinsic           :: min
    integer(ip)         :: minmn
    real(dp), parameter :: zero = 0.0_dp

    minmn = min(m,n)
    wm    = y
    call Hprod (m,hy,wm)
    wn(1:minmn) = d(1:minmn)*wm(1:minmn)
    wn(m+1:n)   = zero
    call Hprod (n,hz,wn)
    x = x + wn

  end subroutine Aprod2

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine lstp  (m,n,nduplc,npower,damp,x,b,Acond,rnorm)

    integer(ip), intent(in)    :: m, n, nduplc, npower
    real(dp),    intent(in)    :: damp
    real(dp),    intent(inout) :: x(n)
    real(dp),    intent(out)   :: b(m)
    real(dp),    intent(out)   :: Acond, rnorm

    !-------------------------------------------------------------------
    ! lstp  generates a sparse least-squares test problem of the form
    !           (   A    )*x = ( b ) 
    !           ( damp*I )     ( 0 )
    ! for solution by USYMLQ, or a sparse underdetermined system
    !            Ax + damp*s = b
    ! for solution by USYMLQ or CRAIG.  The matrix A is m by n and is
    ! constructed in the form  A = Y*D*Z,  where D is an m by n
    ! diagonal matrix, and Y and Z are Householder transformations.
    !
    ! m and n may contain any positive values.
    ! If m >= n  or  damp = 0, the true solution is x as given.
    ! Otherwise, x is modified to contain the true solution.
    !
    ! 1982---1991: Various versions implemented.
    ! 06 Feb 1992: lstp generalized to allow any m and n.
    ! 07 Sep 2007: Line by line translation for Fortran 90 compilers
    !              by Eric Badel <badel@nancy.inra.fr>.
    ! 23 Sep 2007: Fortran 90 version with modules.
    !-------------------------------------------------------------------

    intrinsic           :: min, cos, sin, sqrt
    integer(ip)         :: i, j, minmn
    real(dp)            :: alfa, beta, dampsq, fourpi, t
    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp

    !-------------------------------------------------------------------
    ! Make two vectors of norm 1.0 for the Householder transformations.
    ! fourpi  need not be exact.
    !-------------------------------------------------------------------
    minmn  = min(m,n)
    dampsq = damp**2
    fourpi = 4.0_dp * 3.141592_dp
    alfa   = fourpi / m
    beta   = fourpi / n

    do i = 1,m
       hy(i) = sin( alfa*i )
    end do

    do i = 1,n
       hz(i) = cos( beta*i )
    end do

    alfa   = dnrm2 ( m, hy, 1 )
    beta   = dnrm2 ( n, hz, 1 )
    call dscal ( m, (one/alfa), hy, 1 )
    call dscal ( n, (one/beta), hz, 1 )

    !-------------------------------------------------------------------
    ! Set the diagonal matrix D.  These are the singular values of A.
    !-------------------------------------------------------------------
    do i = 1,minmn
       j    = (i-1+nduplc) / nduplc
       t    =  j*nduplc
       t    =  t / minmn
       d(i) =  t**npower
    end do

    Acond  = (d(minmn)**2 + dampsq) / (d(1)**2 + dampsq)
    Acond  = sqrt( Acond )

    !-------------------------------------------------------------------
    ! If m >=n, the input x will be the true solution.
    ! If m < n, reset x = different true solution.
    ! It must be of the form x = Z(w1) for some w1.
    !                             (0 )
    ! We get w1 from the top of w = Zx.
    !-------------------------------------------------------------------
    wn = x
    call Hprod (n,hz,wn)
    if (m < n) then
       wn(m+1:n) = zero
       call Hprod (n,hz,wn)
       x = wn
    end if

    ! Let r = Y rbar and xbar = Z x.
    ! Solve D r1bar = damp^2 x1bar, where r1bar is in wm(1:minmn).

    wm(1:minmn)   = dampsq * wn(1:minmn)/d(1:minmn)

    wm(minmn+1:m) = one      ! Set r2bar to be anything (empty if m <= n).
    call Hprod (m,hy,wm)     ! Form r = Y rbar  in wm.

    rnorm  = dnrm2 (m,wm,1)  ! Compute norm(r)
    b      = wm              ! and  b = r + Ax.
    call Aprod1(m,n,x,b)

  end subroutine lstp

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine usymlqtest(m,n,nduplc,npower,damp,nout)

    integer(ip), intent(in)    :: m, n, nduplc, npower, nout
    real(dp),    intent(in)    :: damp

    !-------------------------------------------------------------------
    ! This is an example driver routine for running USYMLQ.
    ! It generates a test problem, solves it, and examines the results.
    !
    ! 1982---1991: Various versions implemented.
    ! 23 Sep 2007: f90 version.
    ! 27 Sep 2007: Reformatted a little.
    !------------------------------------------------------------------------

    real(dp)                :: b(m), c(n), x(n), xtrue(n)
    integer(ip)             :: inform, istop, itn, itnlim, j, minmn, nprint
    real(dp)                :: atol, btol, conlim, Anorm, Acond,     &
                               rnorm, enorm, etol, wnorm, xnorm

    real(dp),     parameter :: one = 1.0_dp
    character(*), parameter :: line = '----------------------------------'

    !------------------------------------------------------------------------
    ! Generate the specified test problem
    ! and check that Aprod1, Aprod2 form y + Ax and x + A'y consistently.
    !------------------------------------------------------------------------
    do j = 1,n                              ! Set the desired solution xtrue.
    !  xtrue(j) = real(n-j,dp)              ! Same as f77 usymlq.for 
       xtrue(j) = 0.1_dp*j                  ! For least-squares problems, this is it.
    end do                                  ! If m<n, lstp will alter it.

    minmn  = min(m,n)                       ! Allocate arrays.
    allocate( d(minmn), hy(m), hz(n) )      ! Vectors defining A = Y D Z.
    allocate( wm(m), wn(n) )                ! Work vectors for Aprod1, Aprod2.

    call lstp  ( m,n,nduplc,npower,damp, &  ! Generate test problem.
                 xtrue,b,Acond,rnorm )      ! If m<n, xtrue is altered.

    write(nout,1000) line,line,m,n,nduplc,npower,damp,Acond,rnorm, &
                     line,line

    call Acheck( m,n,Aprod1,Aprod2,nout,inform )   ! Check Aprod1, Aprod2.
          
    if (inform > 0) then
       write(nout,'(a)') 'Check tol in subroutine Acheck in usymlqCheckModule'
       stop
    end if

    !------------------------------------------------------------------------
    ! Set input parameters for USYMLQ
    ! and solve the problem defined by Aprod1, Aprod2, b, damp.
    !------------------------------------------------------------------------
    atol   = 1e-12_dp
    btol   = atol
    conlim = 1000.0_dp * Acond
    itnlim = 1*(m + n + 50)

    c      = one  ! An arbitrary vector.  c = x is best but usually impossible!

    call USYMLQ( m, n, Aprod1, Aprod2, b, c, x,         &
                 atol, btol, conlim, itnlim, nout,      &
                 istop, itn, Anorm, Acond, rnorm, xnorm )

    call xcheck( m, n, Aprod1, Aprod2, b, damp, x, &
                 Anorm, atol, nout, inform )

    nprint = min(m,n,8)
    write(nout,2500)    (j, x(j), j=1,nprint)    ! Print some of the solution

    wn     = x - xtrue                           ! Print a clue about whether
    wnorm  = dnrm2 (n,wn,1)                      ! the solution looks OK.
    xnorm  = dnrm2 (n,xtrue,1)
    enorm  = wnorm/(one + xnorm)
    etol   = 1.0e-3_dp
    if (inform <= 3) then
       write(nout,3000) itn, rnorm, enorm
    else
       write(nout,3100) itn, rnorm, enorm
    end if

    deallocate( wn, wm, hz, hy, d )              ! Free arrays (in reverse)
    return

 1000 format(1p &
      // 1x, 2a &
      /  ' Least-Squares Test Problem      P(', 4i5, e12.2, ' )'     &
      /  ' Condition no. =', e12.4, '     Residual function =', e17.9&
      /  1x, 2a)
 2500 format(/' Solution  x:'          / 4(i6, g14.6))
 3000 format(1p, " USYMLQ  appears to be successful.  Itns =", i7,  &
         "  norm(r) =", e9.2, "  x rel err =", e9.2)
 3100 format(1p, " USYMLQ  appears to have failed.    Itns =", i7,  &
         "  norm(r) =", e9.2, "  x rel err =", e9.2)

  end subroutine usymlqtest

end module usymlqTestModule
