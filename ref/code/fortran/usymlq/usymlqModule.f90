!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File usymlqModule.f90
!
!     USYMLQ   normlz
!
! USYMLQ solves Ax = b 
! using the orthogonal tridiagonalization algorithm described in
!
!    M. A. Saunders, H. D. Simon, and E. L. Yip,
!    Two conjugate-gradient-type methods for unsymmetric linear equations,
!    SIAM J. Numer. Anal. 25(4) 927-940 (1988)
!
! generalized from n x n to m x n systems as described in
!
!    L. Reichel and Q. Ye,
!    A generalized LSQR algorithm,
!    Numer. Linear Algebra Appl. 15 643-660 (2008).
!
! Michael Saunders <saunders@stanford.edu>
! Systems Optimization Laboratory
! Dept of Management Science and Engineering
! Stanford University
!
! Dominique Orban  <dominique.orban@gerad.ca>
! GERAD and Mathematics and Industrial Engineering Dept
! Ecole Polytechnique
! Montreal, Canada
!
! 19 Jul 1984: Original f66 version.
! 17 Oct 1988: Converted to f77.
! 24 May 2014: f90 version derived from f77 version.
!              Generalized to m x n systems.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module usymlqModule

  use  usymlqDataModule,  only : ip, dp
  implicit none
  private
  public  :: USYMLQ
  private :: normlz

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine USYMLQ( m, n, Aprod1, Aprod2, b, c, x,         &
                     atol, btol, conlim, itnlim, nout,      &
                     istop, itn, Anorm, Acond, rnorm, xnorm )

    integer(ip), intent(in)  :: m, n, itnlim, nout
    integer(ip), intent(out) :: istop, itn
    real(dp),    intent(in)  :: b(m), c(n)
    real(dp),    intent(out) :: x(n)
    real(dp),    intent(in)  :: atol, btol, conlim
    real(dp),    intent(out) :: Anorm, Acond, rnorm, xnorm

    interface
       subroutine Aprod1(m,n,x,y)                   ! y := y + A*x
         use usymlqDataModule, only : ip, dp
         integer(ip), intent(in)    :: m,n
         real(dp),    intent(in)    :: x(n)
         real(dp),    intent(inout) :: y(m)
       end subroutine Aprod1

       subroutine Aprod2(m,n,x,y)                   ! x := x + A'*y
         use usymlqDataModule, only : ip, dp
         integer(ip), intent(in)    :: m,n
         real(dp),    intent(inout) :: x(n)
         real(dp),    intent(in)    :: y(m)
       end subroutine Aprod2
    end interface

    !------------------------------------------------------------------
    ! USYMLQ finds a solution x to the system of linear equations Ax = b
    ! where A is a real nonsingular matrix with m rows and n columns,
    ! and b is a real m-vector.
    ! The matrix A is treated as a linear operator.  It is accessed
    ! by means of subroutine calls with the following purpose:
    !
    ! call Aprod1(m,n,x,y)  must compute y = y + A*x  without altering x.
    ! call Aprod2(m,n,x,y)  must compute x = x + A'*y without altering y.
    !
    ! USYMLQ uses an iterative method to approximate the solution.
    ! The number of iterations required to reach a certain accuracy
    ! depends strongly on the scaling of the problem.  Poor scaling of
    ! the rows or columns of A should therefore be avoided where possible.
    ! In the absence of better information, the nonzero columns of A
    ! should be scaled so that they all have the same Euclidean norm.
    !
    ! istop   Output     An integer giving the reason for termination...
    !
    !            0       x = 0 is the exact solution.
    !                    No iterations were performed.
    !
    !            1       norm(Ax - b) is sufficiently small,
    !                    given the values of atol and btol.
    !
    !            3       An estimate of cond(A) has exceeded conlim.
    !                    The system Ax = b appears to be ill-conditioned.
    !                    Otherwise, there could be an error in the Aprods.
    !
    !            4       norm(Ax - b) is as small as seems reasonable.
    !
    !            6       cond(A) seems to be so large that there is
    !                    not much point in doing further iterations,
    !                    given the precision of this machine.
    !                    There could be an error in the Aprods.
    !
    !            7       The iteration limit itnlim was reached.
    !-------------------------------------------------------------------

    intrinsic              :: abs, mod, sqrt
    integer(ip)            :: i
    real(dp)               :: alfa, beta, beta1, gama, gama1, gamma, delta, epsln, &
                              AAnorm, bnorm, cgnorm, cs, ctol, dbar, eta, gbar,    &
                              oldg, qrnorm, rhs1, rhs2, rtol,                      &
                              s, sn, t, test, t1, xxnorm, x1cg, z, zbar
    real(dp)               :: u1(m), u2(m), v1(n), v2(n), w(n)
    real(dp), parameter    :: zero = 0.0, one = 1.0

    ! Initialize.

    if (nout > 0) then
       write(nout, 1000) m, n, atol, conlim, btol, itnlim
    end if
    ctol   = zero
    if (conlim > zero) ctol = one / conlim
    Acond  = zero
    xxnorm = zero
    itn    = 0
    istop  = 0

    u2     = zero
    v2     = zero
    x      = zero

    ! Set up the first vectors for the tridiagonalization.

    u1 = b
    v1 = c
    call normlz( m, u1, beta1 )
    call normlz( n, v1, gama1 )
    w  = v1
    call Aprod1( m, n, v1, u2 )
    call Aprod2( m, n, v2, u1 )
    alfa = dot_product( u1, u2 )
    u2   = u2 - alfa*u1
    v2   = v2 - alfa*v1
    call normlz( m, u2, beta )
    call normlz( n, v2, gama )

    AAnorm = alfa**2 + gama**2
    gbar   = alfa
    dbar   = beta
    bnorm  = beta1
    qrnorm = beta1
    rhs1   = beta1
    rhs2   = zero
    if (nout > 0) then
       write(nout, 1200)
       test   = one
       write(nout, 1500) itn, x(1), qrnorm, qrnorm, test
    end if

    !===================================================================
    ! Main iteration loop.
    !===================================================================
    do
       itn = itn + 1

       ! Perform the next step of the tridiagonalization.

       oldg   = gama
       u1     = (-gama)*u1
       call Aprod1( m, n, v2, u1 )
       v1     = (-beta)*v1
       call Aprod2( m, n, v1, u2 )
       alfa   = dot_product( u1, u2 )

       do i = 1, m
          t     = u2(i)
          u2(i) = u1(i) - alfa*t
          u1(i) = t
       end do

       call normlz( m, u2, beta )
       alfa   = dot_product( v1, v2 )

       do i = 1, n
          t     = v2(i)
          v2(i) = v1(i) - alfa*t
          v1(i) = t
       end do

       call normlz( n, v2, gama )
                                                      
       ! Use a plane rotation on the right to eliminate the
       ! super-diagonal element (oldg) of the tridiagonal matrix.

       gamma  =   sqrt( gbar**2 + oldg**2 )
       cs     =   gbar / gamma
       sn     =   oldg / gamma
       delta  =   cs*dbar + sn*alfa
       gbar   =   sn*dbar - cs*alfa
       epsln  =   sn*beta
       dbar   = - cs*beta

       ! Update x and w.

       z      =   rhs1 / gamma
       s      =   z * cs
       t      =   z * sn

       do i = 1, n
          x(i)  = (w(i)*s + v1(i)*t) + x(i)
          w(i)  = w(i)*sn - v1(i)*cs
       end do

       ! Test for convergence.

       rhs1   = rhs2 - delta*z
       rhs2   =      - epsln*z
       zbar   = rhs1 / gbar
       eta    = sn*z - cs*zbar
       AAnorm = AAnorm + alfa**2 + beta**2 + gama**2
       Anorm  = sqrt( AAnorm )
       cgnorm = beta*abs( eta )
       qrnorm = qrnorm*sn
       xnorm  = sqrt( xxnorm + zbar**2 )
       xxnorm = xxnorm + z**2

       ! Now use these norms to estimate certain other quantities,
       ! some of which will be small near a solution.

       test   = cgnorm / bnorm
       t1     = test   / (one + anorm*xnorm/bnorm)
       rtol   = btol   + atol * anorm*xnorm/bnorm
       t1     = one + t1
       if (itn >= itnlim) istop = 7
       if (t1  <= one   ) istop = 4

       ! Allow for tolerances set by the user.

       if (test  .le. rtol) istop = 1
       !==================================================================

       ! See if it is time to print something.

       if (nout <= 0) go to 600
       if (n   <= 40         ) go to 400
       if (itn <= 10         ) go to 400
       if (itn >= itnlim - 10) go to 400
       if (mod(itn,10)  ==  0) go to 400
       if (test  <= 10.0*rtol) go to 400
       go to 600

       ! Print a line for this iteration.

400    if (mod(itn,10) == 1) write(nout, 1400)
       x1cg   = x(1) + zbar*w(1)
       write(nout, 1500) itn, x1cg, cgnorm, qrnorm, test, anorm

       ! Stop if appropriate.

600    if (istop > 0) exit
    end do
    !===================================================================
    ! End of iteration loop.
    !===================================================================

    ! Move x to xCG, the solution of T*y = beta1*e1, x = V*y.

    x = x + zbar*w

    ! Print the stopping condition.

    if (nout <= 0) go to 900
    write(nout, 1900) istop, itn, Anorm, bnorm, rnorm, xnorm
    if (istop == 0) write(nout, 2000)
    if (istop == 1) write(nout, 2100)
    if (istop == 4) write(nout, 2400)
    if (istop == 7) write(nout, 2700)
900 rnorm  = cgnorm

1000 format(// 9x, ' USYMLQ     --     solution of unsymmetric Ax = b'    &
             / 9x, ' The matrix A has', i8, ' rows   and', i8,' cols' &
             / 9x, ' atol   =', 1p, e10.2, 12x, 'conlim =', e10.2     &
             / 9x, ' btol   =',     e10.2, 12x, 'itnlim =', i10)
1200 format(// 3x, 'itn', 9x, 'x(1)', 11x, 'cgnorm', 7x, &
                   'qrnorm       test        norm(A)' /)
1400 format(1x)
1500 format(i6, 1p, e20.10, e13.3, 2e13.3, 2e11.2)
1900 format(/ ' istop  =', i2,    15x, 'itn    =', i8     &
            / ' Anorm  =', es12.5, 5x, 'bnorm  =', es12.5 &
            / ' rnorm  =', es12.5, 5x, 'xnorm  =', es12.5)
2000 format(/ ' The exact solution is x = 0.'           )
2100 format(/ ' Ax - b is small enough given atol, btol')
2400 format(/ ' Ax - b is small enough for this machine')
2700 format(/ ' The iteration limit has been reached'   )

  end subroutine USYMLQ

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine normlz( n, x, beta )

    integer(ip), intent(in)    :: n
    real(dp),    intent(inout) :: x(n)
    real(dp),    intent(out)   :: beta

    !------------------------------------------------------------------
    ! normlz computes the Euclidean norm of x and returns it in beta.
    ! If x is nonzero, it is scaled so that norm(X) = 1.
    ! ------------------------------------------------------------------

    real(dp),    parameter :: zero = 0.0, one = 1.0

    beta   = dot_product( x, x )
    beta   = sqrt(beta)
    if (beta > zero) then
       x = (one/beta) * x
    end if

  end subroutine normlz

end module usymlqModule
