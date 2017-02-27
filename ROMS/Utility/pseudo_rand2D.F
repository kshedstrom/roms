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
      USE mod_iounits

      integer,parameter :: mnproc = 1

      CONTAINS

      SUBROUTINE initfftdim(nx,ny,fnx,fny)
!
!=======================================================================
! Get FFT dimensions (next power of two)
!=======================================================================
!
      !Initialize fft dimensions and get fnx and fny
      implicit none
      integer, intent(in)  :: nx,ny
      integer, intent(out) :: fnx,fny
      fnx = ceiling(log(float(nx))/log(2.))
      fnx = 2**fnx
      fny = ceiling(log(float(ny))/log(2.))
      fny = 2**fny
      IF (mnproc==1) THEN
         WRITE(stdout,&
         '(" initfftdim: Fourier transform dimensions ",2i6)') fnx,fny
      END IF

      END SUBROUTINE


      SUBROUTINE pseudo2D(Amat,nx,ny,rh,dx,dy)
!
!=======================================================================
! This routine calculates the pseudo random filds using
! the procedure outlined in Evensen (1994).
!=======================================================================
!
      implicit none
      integer, intent(in) :: nx,ny           ! horizontal dimensions
      real*8, intent(out) :: Amat(nx,ny)     ! generated random fields
      real, intent(in)    :: rh              ! Horizontal decorrelation length
      real, intent(in)    :: dx,dy           ! horizontal size of input grid

      integer             :: n1,n2           ! horizontal dimensions in fft grid
      integer             :: nmax, j1, j2, n1n2max
      real*8, allocatable :: w(:)
      real, save          :: rh_save=0.0  ! saving rh used in preprocessing.  Allow for new call to
                                          ! zeroin if rh changes.
      real, save    :: sigma,sigma2
      real, save    :: c
      integer, save :: n1_save=0
      integer, save :: n2_save=0
      integer       :: l,p,n,m,i
      real          :: kappa2,lambda2,kappa,lambda
      real          :: pi2,deltak,asum,adjust
      real          :: a1,b1,tol,fval
      real, allocatable    :: phi(:,:)
      real, allocatable    :: y(:,:)   ! Physical field
      real, parameter :: pi=3.141592653589
      
      ! Variables for FFTW 
      TYPE(C_PTR)                   :: plan
      real                          :: real_part, imag_part
      complex(C_DOUBLE_COMPLEX),&
            & allocatable           :: Qinv(:,:),x(:,:)
  

      CALL initfftdim(nx,ny,n1,n2)
 
      ALLOCATE(phi(-n1/2:n1/2,-n2/2:n2/2))
      ALLOCATE(y(0:2*n1-1,0:n2-1))
      ALLOCATE(x(n1,n2))
      ALLOCATE(Qinv(-n1/2:n1/2,-n2/2:n2/2))
   
      y=0.0
   
      n1n2max = 2 + ceiling(sqrt(max(n1*1.0,n2*1.0)))
      n1n2max = max(n1/2,n2/2)
      ALLOCATE(w(0:n1n2max))
   
      nmax=nx
      pi2=2.0*pi
      deltak=pi2**2/(float(n1*n2)*dx*dy)
      kappa=pi2/(float(n1)*dx)
      kappa2=kappa**2
      lambda=pi2/(float(n2)*dy)
      lambda2=lambda**2
   
      IF (rh /= rh_save .OR. n1 /= n1_save .OR. n2 /= n2_save) THEN
         rh_save=rh
         n1_save=n1
         n2_save=n2
   
         IF (mnproc==1) PRINT *,&
            'pseudo2D: Solving for sigma. rh, dx and dy =',rh,dx,dy
         a1=0.1e-07
         b1=0.1e-06
         tol=0.1e-10
         call zeroin(func2D,sigma,a1,b1,tol,rh,dx,dy,fval,n1,n2)
   
         sigma2=sigma**2
         asum=0.0
         DO p=-n2/2+1,n2/2
            DO l=-n1/2+1,n1/2
               asum=asum+exp(-2.0*(kappa2*float(l*l)+&
                    lambda2*float(p*p))/sigma2)
            ENDDO
         ENDDO
         c=sqrt(1.0/(deltak*asum))
   
         IF (mnproc==1) PRINT *,'Ranfield: sigma  ',sigma
         IF (mnproc==1) PRINT *,'Ranfield: c=     ',c
      ENDIF
   
      !call random_seed() <-- Do not call random_seed here ! bjchoi. 2014.12.16.

      ! Calculating the random wave phases
      call RANDOM_NUMBER(phi)
      phi=pi2*phi

      ! Calculating the wave amplitudes (A3)
      DO p=-n2/2,n2/2-1
       DO l=-n1/2,n1/2-1

         real_part = &
            EXP(-(kappa2*FLOAT(l*l)+lambda2*FLOAT(p*p))/sigma2)*&
            COS(phi(l,p))*SQRT(deltak)*c

         imag_part = &
            EXP(-(kappa2*FLOAT(l*l)+lambda2*FLOAT(p*p))/sigma2)*&
            SIN(phi(l,p))*SQRT(deltak)*c


         Qinv(MOD(l+n1,n1)-n1/2,MOD(p+n2,n2)-n2/2) = &
              CMPLX(real_part,imag_part)
       ENDDO
      ENDDO
     !write(*,*) ' max, min of x amplitude = ', maxval(x),minval(x)

     ! 2-dimensional complex discrete fourier transfrom
     ! call cdft2d(nmax, 2*n1, n2, 1, x, ip, w)
     plan = fftw_plan_dft_2d(n1,n2, Qinv, x,FFTW_BACKWARD,&
            FFTW_ESTIMATE)
     CALL fftw_execute_dft(plan,Qinv, x)
     CALL dfftw_destroy_plan(plan)

      ! store the result on the physical grid y(j1,j2)
      DO j2 = 0, n2 - 1
        DO j1 = 0, n1 - 1
              y(j1, j2) = x(j1, j2)
        END DO
      END DO

      ! save the data
      asum = 0.0
      DO m=1,ny
         DO i=1,nx
            Amat(i,m)=y(i-1,m-1)
            asum = asum + Amat(i,m)**2
         ENDDO
      ENDDO

      ! check if variance of the pseudo ran field is one
      ! if it is not one, amplify the magnitude of Amat
      ! added by bjchoi 2014. 12. 19.
       asum = asum / (nx*ny)
       adjust = SQRT( 1.00/ asum )
       Amat(:,:)=Amat(:,:)*adjust

      !print *,'min,max y   :',minval(y),maxval(y)

      DEALLOCATE(phi, y, x, Qinv)

      END SUBROUTINE pseudo2D



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
      real, external :: func
      integer n1,n2
      real zeropkt,length,dx,dy
      real ax   ! left endpoint of initial interval
      real bx   ! right endpoint of initial interval
      real tol  !  desired length of the interval of uncertainty of the
      real  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      real  abs,sign,fval
   
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
         !if (mnproc==1) write(*,'(2(a,g13.5))')'zeroin: No convergence, ax=',ax,'
         !bx=',bx
           CALL xcstop('(zeroin)')
           STOP  '(zeroin)'
         ENDIF
      ENDIF
     ! begin step
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 IF (abs(fc) .GE. abs(fb)) GO TO 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
     ! convergence test
   40 tol1 = 2.0*eps*abs(b) + 0.5*tol
      xm = .5*(c - b)
      IF (abs(xm) .LE. tol1) GO TO 90
      IF (fb .EQ. 0.0) GO TO 90
     ! is bisection necessary
      IF (abs(e) .LT. tol1) GO TO 70
      IF (abs(fa) .LE. abs(fb)) GO TO 70
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
      p = abs(p)
     ! is interpolation acceptable
      IF ((2.0*p) .GE. (3.0*xm*q - abs(tol1*q))) GO TO 70
      IF (p .GE. abs(0.5*e*q)) GO TO 70
      e = d
      d = p/q
      GO TO 80
     ! bisection
   70 d = xm
      e = d
     ! complete step
   80 a = b
      fa = fb
      IF (abs(d) .GT. tol1) b = b + d
      IF (abs(d) .LE. tol1) b = b + sign(tol1, xm)
      fb = func(b,length,dx,dy,n1,n2)
      IF ((fb*(fc/abs(fc))) .GT. 0.0) GO TO 20
      GO TO 30
     ! done
   90 zeropkt = b
      fval=func(b,length,dx,dy,n1,n2)
      END SUBROUTINE zeroin



      real FUNCTION func2D(sigma,length,dx,dy,n1,n2)
!=======================================================================
! Function used to calculate $sigma$ and $c$.
!=======================================================================
      implicit none
      real sum1,sum2,sigma,length
      real sigma2,pi2,kappa,kappa2,lambda,lambda2,dx,dy
      integer l,p,n1,n2
      real, parameter :: pi=3.141592653589
   
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
            sum1=sum1+exp(-2.0*(kappa2*float(l*l)&
            +lambda2*float(p*p))/sigma2)*cos(kappa*float(l)*length)
         ENDDO
      ENDDO
   
      ! Calculate sum2
      sum2=0.0
      DO p=-n2/2+1,n2/2
         DO l=-n1/2+1,n1/2
            sum2=sum2+exp(-2.0*(kappa2*float(l*l)&
                 +lambda2*float(p*p))/sigma2)
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
      WRITE(stdout,*) 'same pseudo random field will be generated for',&
           ' the given  seed'
      WRITE(stdout,*) 'your random seeds = ', a_seed(8), a_seed(1)
      DEALLOCATE(a_seed)

      END SUBROUTINE random_seed_fixed


      SUBROUTINE xcstop(cerror)
!=======================================================================
! Routine to handle error
!=======================================================================
      implicit none
      character*(*), intent(in) :: cerror
      IF     (cerror .NE. ' ') THEN
         WRITE(stdout,*) '*********************************************'
         WRITE(stdout,*) cerror
         WRITE(stdout,*) '*********************************************'
      ENDIF
      STOP '(xcstop)'
      END SUBROUTINE xcstop

      END MODULE pseudo_rand2D
#endif
