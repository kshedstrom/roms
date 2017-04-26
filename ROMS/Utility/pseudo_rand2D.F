#include "cppdefs.h"
#ifdef STOPERTURB
      MODULE pseudo_rand2D
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================

      USE mod_fftw3
      USE mod_parallel

      CONTAINS

      SUBROUTINE initfftdim(nx,ny,fnx,fny)
!
!=======================================================================
! Get FFT dimensions (next power of two)
!=======================================================================
!
      implicit none
      integer, intent(in)  :: nx,ny
      integer, intent(out) :: fnx,fny
      fnx = ceiling(log(float(nx))/log(2.))
      fnx = 2**fnx
      fny = ceiling(log(float(ny))/log(2.))
      fny = 2**fny

      END SUBROUTINE initfftdim




      SUBROUTINE pseudo2D(Randmat,nx,ny,rh,dx,dy,plan)
!
!=======================================================================
! This routine calculates the pseudo random filds using
! the procedure outlined in Evensen (1994).
!=======================================================================
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in)    :: nx,ny           ! horizontal dimensions
      real*8, intent(in)     :: rh              ! Horizontal decorrelation length
      real*8, intent(in)     :: dx,dy           ! horizontal size of input grid
      real*8, intent(out)    :: Randmat(nx,ny)  ! generated random fields
      TYPE(C_PTR),intent(in) :: plan            ! plan for FFTW
!
!  Local variable declarations.
!
      integer :: n1,n2                 ! horizontal dimensions in fft grid
      integer :: l,p,n,m,i,j1,j2       ! Loop dummies
      real*8  :: kappa2,lambda2,pi2    ! Squared parameters
      real*8  :: kappa,lambda,deltak   ! Equation parameters
      real*8  :: asum                  ! Summing dummy
      real*8  :: adjust                ! Adjusting coefficient
      real*8  :: a1,b1,tol,fval        ! Parameters of equation solver
      real*8, allocatable :: phi(:,:)  ! Wave phases
      real*8, parameter   :: pi=3.141592653589     
!
!  Saved variable declarations.
!
      real*8, save  :: rh_save=0.0  ! saving rh used in preprocessing.  Allow for new call to
                                    ! zeroin if rh changes.
      real*8, save  :: sigma,sigma2
      real*8, save  :: c
      integer, save :: n1_save=0
      integer, save :: n2_save=0
!
!  Variables for FFTW 
!
      real         :: real_part, imag_part
      type(C_PTR)  :: q,s             ! pointers for FFTW
      complex(C_DOUBLE_COMPLEX), pointer :: Qinv(:,:),x(:,:)

!
!-----------------------------------------------------------------------
!  Get the dimensions of FFT
!-----------------------------------------------------------------------
! 
      CALL initfftdim(nx,ny,n1,n2)
!
!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------
!
      ALLOCATE(phi(-n1/2:n1/2,-n2/2:n2/2))
      ! For FFTW, we allocate arrays with an external C function :
      ! fftw_alloc_complex
      ! (http://www.fftw.org/fftw3_doc/Allocating-aligned-memory-in-Fortran.html)
      q = fftw_alloc_complex(int(n1*n2, C_SIZE_T))
      s = fftw_alloc_complex(int(n1*n2, C_SIZE_T))
      CALL c_f_pointer(q, x, [n1,n2])
      CALL c_f_pointer(s, Qinv, [n1,n2])
!
!-----------------------------------------------------------------------
!  Compute equation parameters
!-----------------------------------------------------------------------
! 
      pi2=2.0*pi
      deltak=pi2**2/(float(n1*n2)*dx*dy)
      kappa=pi2/(float(n1)*dx)
      kappa2=kappa**2
      lambda=pi2/(float(n2)*dy)
      lambda2=lambda**2
!
!-----------------------------------------------------------------------
!  Solve the sigma equation (only the first time or if rh/n1/n2 change)
!-----------------------------------------------------------------------
!
      IF (rh /= rh_save .OR. n1 /= n1_save .OR. n2 /= n2_save) THEN
         rh_save=rh
         n1_save=n1
         n2_save=n2
         a1=0.1e-07
         b1=0.1e-06
         tol=0.1e-10
         CALL zeroin(func2D,sigma,a1,b1,tol,rh,dx,dy,fval,n1,n2)
   
         sigma2=sigma**2
         asum=0.0
         DO p=-n2/2+1,n2/2
            DO l=-n1/2+1,n1/2
               asum=asum+exp(-2.0*(kappa2*float(l*l)+&
                    lambda2*float(p*p))/sigma2)
            ENDDO
         ENDDO
         c=sqrt(1.0/(deltak*asum))
      ENDIF
!
!-----------------------------------------------------------------------
!  Calculate the random wave phases
!-----------------------------------------------------------------------
!
      call RANDOM_NUMBER(phi)
      phi=pi2*phi
!
!-----------------------------------------------------------------------
!  Calculating the wave amplitudes (A3)
!-----------------------------------------------------------------------
!
      DO p=-n2/2,n2/2-1
       DO l=-n1/2,n1/2-1
         real_part =                                                    &
     &      EXP(-(kappa2*FLOAT(l*l)+lambda2*FLOAT(p*p))/sigma2)*        &
     &      COS(phi(l,p))*SQRT(deltak)*c
         imag_part =                                                    &
     &      EXP(-(kappa2*FLOAT(l*l)+lambda2*FLOAT(p*p))/sigma2)*        &
     &      SIN(phi(l,p))*SQRT(deltak)*c
         Qinv(MOD(l+n1,n1)+1,MOD(p+n2,n2)+1) =                    &
     &      CMPLX(real_part,imag_part)
       END DO
      END DO
!
!-----------------------------------------------------------------------
!  Compute the 2-dimensional inverse discrete fourier transfrom (FFTW)
!-----------------------------------------------------------------------
!
      CALL fftw_execute_dft(plan,Qinv, x)
!
!-----------------------------------------------------------------------
!  Save data in array and compute variance
!-----------------------------------------------------------------------
!
      Randmat(1:nx,1:ny) = x(1:nx,1:ny)
      asum = SUM(Randmat*Randmat) / (nx*ny)
!
!-----------------------------------------------------------------------
!  Check if variance of the pseudo random field is one
!  If it is not one, amplify the magnitude of Randmat
!-----------------------------------------------------------------------
!
      adjust = SQRT( 1.00 / asum )
      Randmat(:,:) = Randmat(:,:) * adjust
!
!-----------------------------------------------------------------------     
!  Deallocate arrays
!-----------------------------------------------------------------------
!
      DEALLOCATE(phi)
      call fftw_free(q)
      call fftw_free(s)

      END SUBROUTINE pseudo2D


      SUBROUTINE random2D(Randmat,nx,ny)
!
!=======================================================================
! This routine calculates the uniform random fields
!=======================================================================
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: nx,ny           ! horizontal dimensions
      real*8, intent(out) :: Randmat(nx,ny)  ! generated random fields
!
!  Local variable declarations.
!
      integer :: l,p                   ! Loop dummies
      real*8  :: asum                  ! Summing dummy
      real*8  :: adjust                ! Adjusting coefficient

      CALL RANDOM_NUMBER(Randmat)
      asum=SUM(Randmat*Randmat)
      asum = asum / (nx*ny)
      adjust = SQRT( 1.00 / asum )
      Randmat(:,:) = Randmat(:,:) * adjust

      END SUBROUTINE random2D





      SUBROUTINE zeroin(func,zeropkt,ax,bx,tol,length,dx,dy,fval,n1,n2)
!
!=======================================================================
! Finds zero of function f.
! A zero of the function  $func(x,length,dx,n1,n2)$
! is computed in the interval $[ax,bx]$.
! Zeroin| returns a zero $x$ in the given interval
! to within a tolerance  $4*macheps*abs(x) + tol$, where macheps
! is the relative machine precision.

! This function subprogram is a slightly  modified  translation  of
! the algol 60 procedure  zero  given in  richard brent, algorithms for
! minimization without derivatives, prentice - hall, inc. (1973).
!=======================================================================
!
      real*8, external :: func
      integer     :: n1,n2
      real*8      :: zeropkt
      real*8      :: length,dx,dy
      real*8      :: ax   ! left endpoint of initial interval
      real*8      :: bx   ! right endpoint of initial interval
      real*8      :: tol  !  desired length of the interval of uncertainty of the
      real*8      :: a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      real*8      :: fval
   
     ! compute eps, the relative machine precision
      icorr=0
      eps = 1.0
   10 eps = eps/2.0
      tol1 = 1.0 + eps
      IF (tol1 .GT. 1.0) GO TO 10
      ! initialization
   77 a = ax
      b = bx
      fa = func(a,length,dx,dy,n1,n2)
      fb = func(b,length,dx,dy,n1,n2)
      IF (fa*fb .GT. 0.0) THEN
         ax=0.1*ax
         bx=10.0*bx
         icorr=icorr+1
         IF (icorr < 20) THEN
           GOTO 77
         ELSE
           CALL xcstop('(zeroin)')
           STOP  '(zeroin)'
         ENDIF
      ENDIF
     ! begin step
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 IF (ABS(fc) .GE. ABS(fb)) GO TO 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
     ! convergence test
   40 tol1 = 2.0*eps*ABS(b) + 0.5*tol
      xm = .5*(c - b)
      IF (ABS(xm) .LE. tol1) GO TO 90
      IF (fb .EQ. 0.0) GO TO 90
     ! is bisection necessary
      IF (ABS(e) .LT. tol1) GO TO 70
      IF (ABS(fa) .LE. ABS(fb)) GO TO 70
     ! is quadratic interpolation possible
      IF (a .NE. c) GO TO 50
     ! linear interpolation
      s = fb/fa
      p = 2.0*xm*s
      q = 1.0 - s
      GO TO 60
     ! inverse quadratic interpolation
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0))
      q = (q - 1.0)*(r - 1.0)*(s - 1.0)
     ! adjust signs
   60 IF (p .GT. 0.0) q = -q
      p = ABS(p)
     ! is interpolation acceptable
      IF ((2.0*p) .GE. (3.0*xm*q - ABS(tol1*q))) GO TO 70
      IF (p .GE. ABS(0.5*e*q)) GO TO 70
      e = d
      d = p/q
      GO TO 80
     ! bisection
   70 d = xm
      e = d
     ! complete step
   80 a = b
      fa = fb
      IF (ABS(d) .GT. tol1) b = b + d
      IF (ABS(d) .LE. tol1) b = b + SIGN(tol1, xm)
      fb = func(b,length,dx,dy,n1,n2)
      IF ((fb*(fc/ABS(fc))) .GT. 0.0) GO TO 20
      GO TO 30
     ! done
   90 zeropkt = b
      fval=func(b,length,dx,dy,n1,n2)
      END SUBROUTINE zeroin



      real*8 FUNCTION func2D(sigma,length,dx,dy,n1,n2)
!=======================================================================
! Function used to calculate $sigma$ and $c$.
!=======================================================================
      implicit none
      real*8  :: sum1,sum2,sigma,length
      real*8  :: sigma2,pi2,kappa,kappa2,lambda,lambda2
      real*8  :: dx,dy
      integer :: l,p,n1,n2
      
      real*8, parameter :: pi=3.141592653589
   
      sigma2=sigma**2
      pi2=2.0*pi
      kappa=pi2/(float(n1)*dx)
      kappa2=kappa**2
      lambda=pi2/(float(n2)*dy)
      lambda2=lambda**2

      ! Calculate sum1
      sum1=0.0
      DO p=-n2/2+1,n2/2
         DO l=-n1/2+1,n1/2
            sum1=sum1+exp(-2.0*(kappa2*float(l*l)                       &
     &          +lambda2*float(p*p))/sigma2)*cos(kappa*float(l)*length)
         ENDDO
      ENDDO
   
      ! Calculate sum2
      sum2=0.0
      DO p=-n2/2+1,n2/2
         DO l=-n1/2+1,n1/2
            sum2=sum2+exp(-2.0*(kappa2*float(l*l)                       &
     &          +lambda2*float(p*p))/sigma2)
         ENDDO
      ENDDO
      func2D = sum1/sum2 - exp(-1.0)
      END FUNCTION func2D



      SUBROUTINE random_seed_fixed(seed_user)
!=======================================================================
! Routine to initialize seed
!=======================================================================
      implicit none
      ! ----- variables for portable seed setting -----
      integer, intent (in)               :: seed_user
      integer                            :: i_seed
      integer, dimension(:), allocatable :: a_seed
      !integer, dimension(1:8) :: dt_seed
      ! ----- end of variables for seed setting -----
      ! ----- Set up random seed portably -----
      CALL RANDOM_SEED(size=i_seed)
      ALLOCATE(a_seed(1:i_seed))
      CALL RANDOM_SEED(get=a_seed)
      !CALL DATE_AND_TIME(values=dt_seed)
      a_seed(i_seed)=seed_user; a_seed(1)=seed_user*60+100
      CALL RANDOM_SEED(put=a_seed)
      DEALLOCATE(a_seed)

      END SUBROUTINE random_seed_fixed


      SUBROUTINE xcstop(cerror)
!=======================================================================
! Routine to handle error
!=======================================================================
      implicit none
      character*(*), intent(in) :: cerror
      IF     (cerror .NE. ' ') THEN
         WRITE(6,*) '*********************************************'
         WRITE(6,*) cerror
         WRITE(6,*) '*********************************************'
      ENDIF
      STOP '(xcstop)'
      END SUBROUTINE xcstop

      END MODULE pseudo_rand2D
#endif
