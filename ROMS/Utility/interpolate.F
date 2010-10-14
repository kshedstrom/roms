#include "cppdefs.h"

#undef CUBIC_MASKED    /* Bicubic interpolation with mask, needs work */

      MODULE interpolate_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This module contains several all purpuse generic routines:          !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    cinterp2d     Bicubic  interpolation for any 2D field.            !
!    linterp2d     Bilinear interpolation for any 2D field.            !
!    hindices      Finds model grid cell for any datum.                !
!    try_range     Binary search of model grid cell for any datum.     !
!    inside        Closed polygon datum search.                        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds

      implicit none

      CONTAINS

      SUBROUTINE linterp2d (ng, LBx, UBx, LBy, UBy,                     &
     &                      Xinp, Yinp, Finp,                           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Istr, Iend, Jstr, Jend,                     &
     &                      Iout, Jout,                                 &
     &                      Xout, Yout,                                 &
     &                      Fout, Fmin, Fmax)
!
!=======================================================================
!                                                                      !
!  Given any gridded 2D field, Finp, this routine linearly interpolate !
!  to locations (Xout,Yout).  To facilitate the  interpolation  within !
!  any irregularly gridded 2D field,  the fractional grid cell indices !
!  (Iout,Jout) with respect Finp are needed at input.  Notice that the !
!  routine "hindices" can be used to compute these indices.            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     LBx        I-dimension lower bound of gridded field, Finp.       !
!     UBx        I-dimension upper bound of gridded field, Finp.       !
!     LBy        J-dimension lower bound of gridded field, Finp.       !
!     UBy        J-dimension upper bound of gridded field, Finp.       !
!     Xinp       X-locations of gridded field, Finp.                   !
!     Yinp       Y-locations of gridded field, Finp.                   !
!     Finp       2D field to interpolate from.                         !
!     LBi        I-dimension Lower bound of data to interpolate, Fout. !
!     UBi        I-dimension Upper bound of data to interpolate, Fout. !
!     LBj        J-dimension Lower bound of data to interpolate, Fout. !
!     UBj        J-dimension Upper bound of data to interpolate, Fout. !
!     Istr       Starting data I-index to interpolate, Fout.           !
!     Iend       Ending   data I-index to interpolate, Fout.           !
!     Jstr       Starting data J-index to interpolate, Fout.           !
!     Jend       Ending   data J-index to interpolate, Fout.           !
!     Iout       I-fractional Xinp grid cell containing Xout.          !
!     Jout       J-fractional Yinp grid cell containing Yout.          !
!     Xout       X-locations to interpolate, Fout.                     !
!     Yout       Y-locations to interpolate, Fout.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Interpolated 2D field.                                !
!     Fmin       Interpolated field minimum value.                     !
!     Fmax       Interpolated field maximum value.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBx, UBx, LBy, UBy
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
      real(r8), intent(in) :: Xinp(LBx:UBx,LBy:UBy)
      real(r8), intent(in) :: Yinp(LBx:UBx,LBy:UBy)
      real(r8), intent(in) :: Finp(LBx:UBx,LBy:UBy)

      real(r8), intent(in) :: Iout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Jout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Xout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Yout(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: Fout(LBi:UBi,LBj:UBj)
!
      real(r8), intent(out) :: Fmin, Fmax
!
!  Local variable declarations.
!
      integer :: i, i1, i2, j, j1, j2

      real(r8) :: p1, p2, q1, q2
!
!-----------------------------------------------------------------------
!  Linearly interpolate requested field
!-----------------------------------------------------------------------
!
      Fmin=1.0E+35_r8
      Fmax=-1.0E+35_r8
      DO j=Jstr,Jend
        DO i=Istr,Iend
          i1=INT(Iout(i,j))
          i2=i1+1
          j1=INT(Jout(i,j))
          j2=j1+1
          IF (((LBx.le.i1).and.(i1.le.UBx)).and.                        &
     &        ((LBy.le.j1).and.(j1.le.UBy))) THEN
            p2=REAL(i2-i1,r8)*(Iout(i,j)-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(Jout(i,j)-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
            Fout(i,j)=p1*q1*Finp(i1,j1)+                                &
     &                p2*q1*Finp(i2,j1)+                                &
     &                p2*q2*Finp(i2,j2)+                                &
     &                p1*q2*Finp(i1,j2)
            Fmin=MIN(Fmin,Fout(i,j))
            Fmax=MAX(Fmax,Fout(i,j))
          END IF
        END DO
      END DO

      END SUBROUTINE linterp2d

      SUBROUTINE cinterp2d (ng, LBx, UBx, LBy, UBy,                     &
     &                      Xinp, Yinp, Finp,                           &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      Istr, Iend, Jstr, Jend,                     &
     &                      Iout, Jout,                                 &
     &                      Xout, Yout,                                 &
     &                      Fout, Fmin, Fmax)
!
!=======================================================================
!                                                                      !
!  Given any gridded 2D field,  Finp, at locations (Xinp,Yinp) this    !
!  routine performs bicubic interpolation at locations (Xout,Yout).    !
!  To facilitate the interpolation within any  irregularly  gridded    !
!  field, the fractional grid cell indices (Iout,Jout) with respect    !
!  Finp are needed at input. Notice that the routine "hindices" can    !
!  be used to compute these indices.                                   !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     LBx        I-dimension lower bound of gridded field, Finp.       !
!     UBx        I-dimension upper bound of gridded field, Finp.       !
!     LBy        J-dimension lower bound of gridded field, Finp.       !
!     UBy        J-dimension upper bound of gridded field, Finp.       !
!     Xinp       X-locations of gridded field, Finp.                   !
!     Yinp       Y-locations of gridded field, Finp.                   !
!     Finp       2D field to interpolate from.                         !
!     LBi        I-dimension Lower bound of data to interpolate, Fout. !
!     UBi        I-dimension Upper bound of data to interpolate, Fout. !
!     LBj        J-dimension Lower bound of data to interpolate, Fout. !
!     UBj        J-dimension Upper bound of data to interpolate, Fout. !
!     Istr       Starting data I-index to interpolate, Fout.           !
!     Iend       Ending   data I-index to interpolate, Fout.           !
!     Jstr       Starting data J-index to interpolate, Fout.           !
!     Jend       Ending   data J-index to interpolate, Fout.           !
!     Iout       I-fractional Xinp grid cell containing Xout.          !
!     Jout       J-fractional Yinp grid cell containing Yout.          !
!     Xout       X-locations to interpolate, Fout.                     !
!     Yout       Y-locations to interpolate, Fout.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Fout       Interpolated 2D field.                                !
!     Fmin       Interpolated field minimum value.                     !
!     Fmax       Interpolated field maximum value.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBx, UBx, LBy, UBy
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend
!
      real(r8), intent(in) :: Xinp(LBx:UBx,LBy:UBy)
      real(r8), intent(in) :: Yinp(LBx:UBx,LBy:UBy)
      real(r8), intent(in) :: Finp(LBx:UBx,LBy:UBy)

      real(r8), intent(in) :: Iout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Jout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Xout(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Yout(LBi:UBi,LBj:UBj)

      real(r8), intent(out) :: Fout(LBi:UBi,LBj:UBj)
!
      real(r8), intent(out) :: Fmin, Fmax
!
!  Local variable declarations.
!
      integer i, ic, iter, i1, i2, j, jc, j1, j2

      real(r8) :: a11, a12, a21, a22
      real(r8) :: e11, e12, e21, e22
      real(r8) :: cff, d1, d2, dfc, dx, dy, eta, xi, xy, yx
      real(r8) :: f0, fx, fxx, fxxx, fxxy, fxy, fxyy, fy, fyy, fyyy

      real(r8), parameter :: C01 = 1.0_r8/48.0_r8
      real(r8), parameter :: C02 = 1.0_r8/32.0_r8
      real(r8), parameter :: C03 = 0.0625_r8                  ! 1/16
      real(r8), parameter :: C04 = 1.0_r8/6.0_r8
      real(r8), parameter :: C05 = 0.25_r8
      real(r8), parameter :: C06 = 0.5_r8
      real(r8), parameter :: C07 = 0.3125_r8                  ! 5/16
      real(r8), parameter :: C08 = 0.625_r8                   ! 5/8
      real(r8), parameter :: C09 = 1.5_r8
      real(r8), parameter :: C10 = 13.0_r8/24.0_r8

      real(r8), parameter :: LIMTR = 3.0_r8
      real(r8), parameter :: spv = 0.0_r8            ! HGA need work

      real(r8), dimension(-1:2,-1:2) :: dfx, dfy, f
!
!-----------------------------------------------------------------------
!  Interpolates requested field locations (Xout,Yout).
!-----------------------------------------------------------------------
!
      Fmin=1.0E+35_r8
      Fmax=-1.0E+35_r8
      DO j=Jstr,Jend
        DO i=Istr,Iend
          i1=INT(Iout(i,j))
          i2=i1+1
          j1=INT(Jout(i,j))
          j2=j1+1
          IF (((LBx.le.i1).and.(i1.le.UBx)).and.                        &
     &        ((LBy.le.j1).and.(j1.le.UBy))) THEN
!
!  Determine local fractional coordinates (xi,eta) corresponding to
!  the target point (Xout,Yout) on the grid (Xinp,Yinp). Here, "xi"
!  and "eta" are defined, in such a way, that xi=eta=0 corresponds
!  to the middle of the cell (i1:i1+1,j1:j1+1), while xi=+/-1/2 and
!  eta=+/-1/2 (any combination +/- signs) corresponds to the four
!  corner points of the cell. Inside the cell it is assumed that
!  (Xout,Yout) are expressed via bi-linear functions of (xi,eta),
!  where term proportional to xi*eta does not vanish because
!  coordinate transformation may be at least weakly non-orthogonal
!  due to discretization errors. The associated non-linear system
!  is solved by iterative method of Newton.
!
            xy=Xinp(i2,j2)-Xinp(i1,j2)-Xinp(i2,j1)+Xinp(i1,j1)
            yx=Yinp(i2,j2)-Yinp(i1,j2)-Yinp(i2,j1)+Yinp(i1,j1)
            dx=Xout(i,j)-0.25_r8*(Xinp(i2,j2)+Xinp(i1,j2)+              &
     &                            Xinp(i2,j1)+Xinp(i1,j1))
            dy=Yout(i,j)-0.25_r8*(Yinp(i2,j2)+Yinp(i1,j2)+              &
     &                            Yinp(i2,j1)+Yinp(i1,j1))
!
!  The coordinate transformation matrix:
!
!           e11 e12
!           e21 e22
!
!  contains derivatives of (Xinp,Yinp) with respect to (xi,eta). Because
!  the coordinates may be non-orthogonal (at least due to discretization
!  errors), the nonlinear system
!
!           e11*xi+e12*eta+xy*xi*eta=dx
!           e21*xi+e22*eta+yx*xi*eta=dy
!
!  needs to be solved in order to retain symmetry.
!
            e11=0.5_r8*(Xinp(i2,j2)-Xinp(i1,j2)+Xinp(i2,j1)-Xinp(i1,j1))
            e12=0.5_r8*(Xinp(i2,j2)+Xinp(i1,j2)-Xinp(i2,j1)-Xinp(i1,j1))
            e21=0.5_r8*(Yinp(i2,j2)-Yinp(i1,j2)+Yinp(i2,j1)-Yinp(i1,j1))
            e22=0.5_r8*(Yinp(i2,j2)+Yinp(i1,j2)-Yinp(i2,j1)-Yinp(i1,j1))
!
            cff=1.0_r8/(e11*e22-e12*e21)
            xi=cff*(e22*dx-e12*dy)
            eta=cff*(e11*dy-e21*dx)
!
            DO iter=1,4
              d1=dx-e11*xi-e12*eta-xy*xi*eta
              d2=dy-e21*xi-e22*eta-yx*xi*eta
              a11=e11+xy*eta
              a12=e12+xy*xi
              a21=e21+yx*eta
              a22=e22+yx*xi
              cff=1.0_r8/(a11*a22-a12*a21)
              xi =xi +cff*(a22*d1-a12*d2)
              eta=eta+cff*(a11*d2-a21*d1)
            END DO

#ifndef CUBIC_MASKED
!
!  Genuinely two-dimensional, isotropic cubic interpolation scheme
!  using 12-point stencil.  In the code below the interpolated field,
!  Fout, is expanded into two-dimensional Taylor series of local
!  fractional coordinates "xi" and "eta", retaining all terms of
!  combined power up to third order (that is, xi, eta, xi^2, eta^2,
!  xi*eta, xi^3, eta^3, xi^2*eta, and xi*eta^2), with all
!  coefficients (i.e, derivatives) computed via           x  x
!  two-dimensional finite difference expressions          |  |
!  of "natural" order of accuracy: 4th-order for       x--x--x--x
!  the field itself and its first derivatives in          |  |
!  both directions; and 2nd-order for all higher-      x--x--x--x
!  order derivatives. The permissible range of            |  |
!  of coordinates is -1/2 < xi,eta < +1/2, which          x--x
!  covers the central cell on the stencil, while
!  xi=eta=0 corresponds to its center. This interpolation scheme has
!  the property that if xi,eta=+/-1/2 (any combination of +/- signs)
!  it reproduces exactly value of the function at the corresponding
!  corner of the central "working" cell. However, it does not pass
!  exactly through the  extreme points of the stencil, where either
!  xi=+/-3/2 or eta+/-3/2. And, unlike a split-directional scheme,
!  when interpolating along the line eta=+/-1/2 (similarly xi=+/-1/2),
!  it has non-zero contribution from points on the side from the line,
!  except if xi=-1/2; 0; +1/2 (similarly eta=-1/2; 0; +1/2).
!
            DO jc=-1,2
              DO ic=-1,2
                f(ic,jc)=Finp(MAX(1,MIN(UBx,i1+ic)),                    &
     &                        MAX(1,MIN(UBy,j1+jc)))
              END DO
            END DO

            f0=C07*(f(1,1)+f(1,0)+f(0,1)+f(0,0))-                       &
     &         C02*(f(2,0)+f(2,1)+f(1,2)+f(0,2)+                        &
     &              f(-1,1)+f(-1,0)+f(0,-1)+f(1,-1))

            fx=C08*(f(1,1)+f(1,0)-f(0,1)-f(0,0))-                       &
     &         C01*(f(2,1)+f(2,0)-f(-1,1)-f(-1,0))-                     &
     &         C03*(f(1,2)-f(0,2)+f(1,-1)-f(0,-1))

            fy=C08*(f(1,1)-f(1,0)+f(0,1)-f(0,0))-                       &
     &         C01*(f(1,2)+f(0,2)-f(1,-1)-f(0,-1))-                     &
     &         C03*(f(2,1)-f(2,0)+f(-1,1)-f(-1,0))

            fxy=f(1,1)-f(1,0)-f(0,1)+f(0,0)

            fxx=C05*(f(2,1)-f(1,1)-f(0,1)+f(-1,1)+                      &
     &               f(2,0)-f(1,0)-f(0,0)+f(-1,0))

            fyy=C05*(f(1,2)-f(1,1)-f(1,0)+f(1,-1)+                      &
     &               f(0,2)-f(0,1)-f(0,0)+f(0,-1))

            fxxx=C06*(f(2,1)+f(2,0)-f(-1,1)-f(-1,0))-                   &
     &           C09*(f(1,1)+f(1,0)-f(0,1)-f(0,0))

            fyyy=C06*(f(1,2)+f(0,2)-f(1,-1)-f(0,-1))-                   &
     &           C09*(f(1,1)-f(1,0)+f(0,1)-f(0,0))

            fxxy=C06*(f(2,1)-f(1,1)-f(0,1)+f(-1,1)-                     &
     &                f(2,0)+f(1,0)+f(0,0)-f(-1,0))

            fxyy=C06*(f(1,2)-f(1,1)-f(1,0)+f(1,-1)-                     &
     &                f(0,2)+f(0,1)+f(0,0)-f(0,-1))
#else
!
!  Algorithm below is equivalent to the one above, except that special
!  care is taken to avoid interpolation accross land. This is achieved
!  by shortening the stencil and reducing order of polynomial, if
!  extreme points of the stencil touch land. This is achieved by
!  expressing all f0,fx,fy,...,fxyy in terms of values of interpolated
!  field at the four corners of central cell (which already checked to
!  stay away from land), and eight one-sided differences dfx,dfy (see
!  below) in such a way that field values at the extreme points of the
!  12-point stencil do not participate directly into f0,fx,...,fxyy.
!  Should an extreme point of the stencil touch land, thus making it
!  impossible to compute the corresponding one-sided difference, this
!  difference is retracted toward the center of the stencil.
!
!  Optionally, a slope-limiting algorithm may be employed to prevent
!  spurious oscillations of the interpolant. This is a valuable property,
!  if dealing with rough data, however, as a side effect, it turns off
!  high-order interpolation in the vicinity of extrema.
!
!  The slope-limiting algorithm employed here checks that two consecutive
!  elementary differences, "dfx" and "dfc" have the same sign and differ
!  in magnitude by no more than factor of 3.
!
            f(0,0)=Finp(i1,j1)
            f(1,0)=Finp(i2,j1)
            f(0,1)=Finp(i1,j2)
            f(1,1)=Finp(i2,j2)
!
            dfc=f(1,1)-f(0,1)
            IF (i1+2.gt.UBx) THEN
              dfx(1,1)=dfc
            ELSE IF (Finp(i1+2,j2).eq.spv) THEN
              dfx(1,1)=dfc
            ELSE
              dfx(1,1)=Finp(i1+2,j2)-f(1,1)
# ifdef LIMTR
              IF ((dfx(1,1)*dfc).lt.0.0_r8) THEN
                dfx(1,1)=0.0_r8
              ELSE IF (ABS(dfx(1,1)).gt.(LIMTR*ABS(dfc))) THEN
                dfx(1,1)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(1,0)-f(0,0)
            IF ((i1+2).gt.UBx) THEN
              dfx(1,0)=dfc
            ELSE IF (Finp(i1+2,j1).eq.spv) THEN
              dfx(1,0)=dfc
            ELSE
              dfx(1,0)=Finp(i1+2,j1)-f(1,0)
# ifdef LIMTR
              IF ((dfx(1,0)*dfc).lt.0.0_r8) THEN
                dfx(1,0)=0.0_r8
              ELSE IF (ABS(dfx(1,0)).gt.(LIMTR*ABS(dfc))) THEN
                dfx(1,0)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(1,1)-f(0,1)
            IF (i1-1.lt.1) THEN
              dfx(0,1)=dfc
            ELSE IF (Finp(i1-1,j2).eq.spv) THEN
              dfx(0,1)=dfc
            ELSE
              dfx(0,1)=f(0,1)-Finp(i1-1,j2)
# ifdef LIMTR
              IF ((dfx(0,1)*dfc).lt.0.0_r8) THEN
                dfx(0,1)=0.0_r8
              ELSE IF (ABS(dfx(0,1)).gt.(LIMTR*ABS(dfc))) THEN
                dfx(0,1)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(1,0)-f(0,0)
            IF (i1-1.lt.1) THEN
              dfx(0,0)=dfc
            ELSE IF (Finp(i1-1,j1).eq.spv) THEN
              dfx(0,0)=dfc
            ELSE
              dfx(0,0)=f(0,0)-Finp(i1-1,j1)
# ifdef LIMTR
              IF ((dfx(0,0)*dfc).lt.0.0_r8) THEN
                dfx(0,0)=0.0_r8
              ELSE IF (ABS(dfx(0,0)).gt.(LIMTR*ABS(dfc))) THEN
                dfx(0,0)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(1,1)-f(1,0)
            IF (j1+2.gt.UBy) THEN
              dfy(1,1)=dfc
            ELSE IF (Finp(i2,j1+2).eq.spv) THEN
              dfy(1,1)=dfc
            ELSE
              dfy(1,1)=Finp(i2,j1+2)-f(1,1)
# ifdef LIMTR
              IF ((dfy(1,1)*dfc).lt.0.0_r8) THEN
                dfy(1,1)=0.0_r8
              ELSEIF (ABS(dfy(1,1)).gt.(LIMTR*ABS(dfc))) THEN
                dfy(1,1)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(0,1)-f(0,0)
            IF (j1+2.gt.UBy) THEN
              dfy(0,1)=dfc
            ELSE IF (Finp(i1,j1+2).eq.spv) THEN
              dfy(0,1)=dfc
            ELSE
              dfy(0,1)=Finp(i1,j1+2)-f(0,1)
# ifdef LIMTR
              IF ((dfy(0,1)*dfc).lt.0.0_r8) THEN
                dfy(0,1)=0.0_r8
              ELSE IF (ABS(dfy(0,1)).gt.(LIMTR*ABS(dfc))) THEN
                dfy(0,1)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(1,1)-f(1,0)
            IF (j1-1.lt.1) THEN
              dfy(1,0)=dfc
            ELSE IF (Finp(i2,j1-1).eq.spv) THEN
              dfy(1,0)=dfc
            ELSE
              dfy(1,0)=f(1,0)-Finp(i2,j1-1)
# ifdef LIMTR
              IF ((dfy(1,0)*dfc).lt.0.0_r8) THEN
                dfy(1,0)=0.0_r8
              ELSE IF (ABS(dfy(1,0)).gt.(LIMTR*ABS(dfc))) THEN
                dfy(1,0)=LIMTR*dfc
              END IF
# endif
            END IF
!
            dfc=f(0,1)-f(0,0)
            IF (j1-1.lt.1) THEN
              dfy(0,0)=dfc
            ELSE IF (Finp(i1,j1-1).eq.spv) THEN
              dfy(0,0)=dfc
            ELSE
              dfy(0,0)=f(0,0)-Finp(i1,j1-1)
# ifdef LIMTR
              IF ((dfy(0,0)*dfc).lt.0.0_r8) THEN
                dfy(0,0)=0.0_r8
              ELSEIF (ABS(dfy(0,0)).gt.(LIMTR*ABS(dfc))) THEN
                dfy(0,0)=LIMTR*dfc
              END IF
# endif
            END IF
!
            f0=C05*(f(1,1)+f(1,0)+f(0,1)+f(0,0))-                       &
     &         C02*(dfx(1,1)+dfx(1,0)-dfx(0,1)-dfx(0,0)+                &
     &              dfy(1,1)-dfy(1,0)+dfy(0,1)-dfy(0,0))

            fx=C10*(f(1,1)-f(0,1)+f(1,0)-f(0,0))-                       &
     &         C01*(dfx(1,1)+dfx(1,0)+dfx(0,1)+dfx(0,0))-               &
     &         C03*(dfy(1,1)-dfy(0,1)-dfy(1,0)+dfy(0,0))

            fy=C10*(f(1,1)-f(1,0)+f(0,1)-f(0,0))-                       &
     &         C01*(dfy(1,1)+dfy(0,1)+dfy(1,0)+dfy(0,0))-               &
     &         C03*(dfx(1,1)-dfx(1,0)-dfx(0,1)+dfx(0,0))

            fxy=f(1,1)-f(1,0)-f(0,1)+f(0,0)

            fxx=C05*(dfx(1,1)-dfx(0,1)+dfx(1,0)-dfx(0,0))

            fyy=C05*(dfy(1,1)-dfy(1,0)+dfy(0,1)-dfy(0,0))

            fxxx=C06*(dfx(1,1)+dfx(1,0)+dfx(0,1)+dfx(0,0))-             &
     &           f(1,1)+f(0,1)-f(1,0)+f(0,0)

            fyyy=C06*(dfy(1,1)+dfy(0,1)+dfy(1,0)+dfy(0,0))-             &
     &           f(1,1)+f(1,0)-f(0,1)+f(0,0)

            fxxy=C06*(dfx(1,1)-dfx(0,1)-dfx(1,0)+dfx(0,0))

            fxyy=C06*(dfy(1,1)-dfy(1,0)-dfy(0,1)+dfy(0,0))
#endif
            Fout(i,j)=f0+                                               &
     &                fx*xi+                                            &
     &                fy*eta+                                           &
     &                C06*fxx*xi*xi+                                    &
     &                fxy*xi*eta+                                       &
     &                C06*fyy*eta*eta+                                  &
     &                C04*fxxx*xi*xi*xi+                                &
     &                C06*fxxy*xi*xi*eta+                               &
     &                C04*fyyy*eta*eta*eta+                             &
     &                C06*fxyy*xi*eta*eta
            Fmin=MIN(Fmin,Fout(i,j))
            Fmax=MAX(Fmax,Fout(i,j))
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE cinterp2d

      SUBROUTINE hindices (ng, LBi, UBi, LBj, UBj,                      &
     &                     Is, Ie, Js, Je,                              &
     &                     angler, Xgrd, Ygrd,                          &
     &                     LBm, UBm, LBn, UBn,                          &
     &                     Ms, Me, Ns, Ne,                              &
     &                     Xpos, Ypos, Ipos, Jpos,                      &
     &                     IJspv, rectangular)
!
!=======================================================================
!                                                                      !
!  Given any geographical locations Xpos and Ypos, this routine finds  !
!  the corresponding array cell indices (Ipos, Jpos) of gridded  data  !
!  Xgrd and Ygrd containing each requested location. This indices are  !
!  usually used for interpolation.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng          Nested grid number.                                  !
!     LBi         I-dimension Lower bound of gridded data.             !
!     UBi         I-dimension Upper bound of gridded data.             !
!     LBj         J-dimension Lower bound of gridded data.             !
!     UBj         J-dimension Upper bound of gridded data.             !
!     Is          Starting gridded data I-index to search.             !
!     Ie          Ending   gridded data I-index to search.             !
!     Js          Starting gridded data J-index to search.             !
!     Je          Ending   gridded data J-index to search.             !
!     angler      Gridded data angle between X-axis and true EAST      !
!                   (radians).                                         !
!     Xgrd        Gridded data X-locations (usually, longitude).       !
!     Ygrd        Gridded data Y-locations (usually, latitude).        !
!     LBm         I-dimension Lower bound of requested locations.      !
!     UBm         I-dimension Upper bound of requested locations.      !
!     LBn         J-dimension Lower bound of requested locations.      !
!     UBn         J-dimension Upper bound of requested locations.      !
!     Ms          Starting requested locations I-index to search.      !
!     Me          Ending   requested locations I-index to search.      !
!     Ns          Starting requested locations J-index to search.      !
!     Ne          Ending   requested locations J-index to search.      !
!     Xpos        Requested X-locations to process (usually longitude).!
!     Ypos        Requested Y-locations to process (usually latitude). !
!     IJspv       Unbounded special value to assign.                   !
!     rectangular Logical switch indicating that gridded data has a    !
!                   plaid distribution.                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Ipos       Fractional I-cell index containing locations in data. !
!     Jpos       Fractional J-cell index containing locations in data. !
!                                                                      !
!  Calls:    Try_Range                                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
#ifdef OFFLINE_FLOATS
      USE mod_iounits
#endif
!
!  Imported variable declarations.
!
      logical, intent(in) :: rectangular
      integer, intent(in) :: ng
      integer, intent(in) :: LBi, UBi, LBj, UBj, Is, Ie, Js, Je
      integer, intent(in) :: LBm, UBm, LBn, UBn, Ms, Me, Ns, Ne

      real(r8), intent(in) :: IJspv
!
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Xgrd(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Ygrd(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: Xpos(LBm:UBm,LBn:UBn)
      real(r8), intent(in) :: Ypos(LBm:UBm,LBn:UBn)

      real(r8), intent(out) :: Ipos(LBm:UBm,LBn:UBn)
      real(r8), intent(out) :: Jpos(LBm:UBm,LBn:UBn)
!
!  Local variable declarations.
!
      logical :: found, foundi, foundj

      integer :: Imax, Imin, Jmax, Jmin, i, i0, j, j0, mp, np

      real(r8) :: aa2, ang, bb2, diag2, dx, dy, phi
      real(r8) :: xfac, xpp, yfac, ypp
#ifdef OFFLINE_FLOATS
      real(r8) :: xh(5), yh(5), xh0, yh0
#endif
!
!-----------------------------------------------------------------------
!  Determine grid cell indices containing requested position points.
!  Then, interpolate to fractional cell position.
!-----------------------------------------------------------------------
!
      DO np=Ns,Ne
        DO mp=Ms,Me
          Ipos(mp,np)=IJspv
          Jpos(mp,np)=IJspv
!
!  The gridded data has a plaid distribution so the search is trivial.
!
          IF (rectangular) THEN
            foundi=.FALSE.
            I_LOOP : DO i=LBi,UBi-1
              IF ((Xgrd(i  ,1).le.Xpos(mp,np)).and.                     &
     &            (Xgrd(i+1,1).gt.Xpos(mp,np))) THEN
                Imin=i
                foundi=.TRUE.
                EXIT I_LOOP
              END IF
            END DO I_LOOP
            foundj=.FALSE.
            J_LOOP : DO j=LBj,UBj-1
              IF ((Ygrd(1,j  ).le.Ypos(mp,np)).and.                     &
     &            (Ygrd(1,j+1).gt.Ypos(mp,np))) THEN
                Jmin=j
                foundj=.TRUE.
                EXIT J_LOOP
              END IF
            END DO J_LOOP
            found=foundi.and.foundj
!
!  Check each position to find if it falls inside the whole domain.
!  Once it is stablished that it inside, find the exact cell to which
!  it belongs by successively dividing the domain by a half (binary
!  search).
!
          ELSE
            found=try_range(ng, LBi, UBi, LBj, UBj,                     &
     &                      Xgrd, Ygrd,                                 &
     &                      Is, Ie, Js, Je,                             &
     &                      Xpos(mp,np), Ypos(mp,np))
            IF (found) THEN
              Imin=Is
              Imax=Ie
              Jmin=Js
              Jmax=Je
              DO while (((Imax-Imin).gt.1).or.((Jmax-Jmin).gt.1))
                IF ((Imax-Imin).gt.1) THEN
                  i0=(Imin+Imax)/2
                  found=try_range(ng, LBi, UBi, LBj, UBj,               &
     &                            Xgrd, Ygrd,                           &
     &                            Imin, i0, Jmin, Jmax,                 &
     &                            Xpos(mp,np), Ypos(mp,np))
                  IF (found) THEN
                    Imax=i0
                  ELSE
                    Imin=i0
                  END IF
                END IF
                IF ((Jmax-Jmin).gt.1) THEN
                  j0=(Jmin+Jmax)/2
                  found=try_range(ng, LBi, UBi, LBj, UBj,               &
     &                            Xgrd, Ygrd,                           &
     &                            Imin, Imax, Jmin, j0,                 &
     &                            Xpos(mp,np), Ypos(mp,np))
                  IF (found) THEN
                    Jmax=j0
                  ELSE
                    Jmin=j0
                  END IF
                END IF
              END DO
              found=(Is.le.Imin).and.(Imin.le.Ie).and.                  &
     &              (Is.le.Imax).and.(Imax.le.Ie).and.                  &
     &              (Js.le.Jmin).and.(Jmin.le.Je).and.                  &
     &              (Js.le.Jmax).and.(Jmax.le.Je)
            END IF
          END IF

#ifdef OFFLINE_FLOATS

          found=try_range(ng, LBi, UBi, LBj, UBj,                       &
     &                    Xgrd, Ygrd,                                   &
     &                    Imin, Imax, Jmin, Jmax,                       &
     &                    Xpos(mp,np), Ypos(mp,np))
           
          IF (.not. found) THEN
            WRITE(stdout,*) 'Initial float ', mp,                       &
     &                      ' location not found in grid.  ... '
            WRITE(stdout,*) '(x,y)=',Xpos(mp,np),Ypos(mp,np)
!            WRITE(stdout,*) '(I,J)min=',Imin,Jmin
!            WRITE(stdout,*) '(I,J)max=',Imax,Jmax
!            WRITE(stdout,*) 'grid (x(0,0),y(0,0))=',                   &
!     &                       Xgrd(Imin,Jmin),Ygrd(Imin,Jmin)
!            WRITE(stdout,*) 'grid (x(0,1),y(0,1))=',                   &
!     &                       Xgrd(Imin,Jmax),Ygrd(Imin,Jmax)
!            WRITE(stdout,*) 'grid (x(1,1),y(1,1))=',                   &
!     &                       Xgrd(Imax,Jmax),Ygrd(Imax,Jmax)
!            WRITE(stdout,*) 'grid (x(1,0),y(1,0))=',                   &
!     &                       Xgrd(Imax,Jmin),Ygrd(Imax,Jmin)
            WRITE(stdout,*) 'Check for float location error or ',       &
     &                      'missing float.'
            WRITE(stdout,*) 'Float grid coords will be default (0.5,0.5).'
          END IF

! New code by HPB, modified by Patrick Drake, added Feb. 2007
!
! knowing the exact cell, calculate the exact indices
!       CAE adding this if statement for mpi case where proc does not have float
        IF (found) THEN
          yh(1) = Ygrd(Imin,Jmin+1)
          yh(2) = Ygrd(Imin+1,Jmin+1)
          yh(3) = Ygrd(Imin+1,Jmin)
          yh(4) = Ygrd(Imin,Jmin)
          xh(1) = Xgrd(Imin,Jmin+1)
          xh(2) = Xgrd(Imin+1,Jmin+1)
          xh(3) = Xgrd(Imin+1,Jmin)
          xh(4) = Xgrd(Imin,Jmin)

          yh0 = Ypos(mp,np)
          xh0 = Xpos(mp,np)

          call quadbisect(mp,xh0,yh0,xh,yh,dx,dy)
          Ipos(mp,np)=REAL(Imin,r8)+dx
          Jpos(mp,np)=REAL(Jmin,r8)+dy
        END IF

! End New code by HPB
#else /* OFFLINE_FLOATS */
!
!  Knowing the correct cell, calculate the exact indices, accounting
!  for a possibly rotated grid.  If spherical, convert all positions
!  to meters first.
!
          IF (found) THEN
            IF (spherical) THEN
              yfac=Eradius*deg2rad
              xfac=yfac*COS(Ypos(mp,np)*deg2rad)
              xpp=(Xpos(mp,np)-Xgrd(Imin,Jmin))*xfac
              ypp=(Ypos(mp,np)-Ygrd(Imin,Jmin))*yfac
            ELSE
              xfac=1.0_r8
              yfac=1.0_r8
              xpp=Xpos(mp,np)-Xgrd(Imin,Jmin)
              ypp=Ypos(mp,np)-Ygrd(Imin,Jmin)
            END IF
!
!  Use Law of Cosines to get cell parallelogram "shear" angle.
!
            diag2=((Xgrd(Imin+1,Jmin)-Xgrd(Imin,Jmin+1))*xfac)**2+      &
     &            ((Ygrd(Imin+1,Jmin)-Ygrd(Imin,Jmin+1))*yfac)**2
            aa2=((Xgrd(Imin,Jmin)-Xgrd(Imin+1,Jmin))*xfac)**2+          &
     &          ((Ygrd(Imin,Jmin)-Ygrd(Imin+1,Jmin))*yfac)**2
            bb2=((Xgrd(Imin,Jmin)-Xgrd(Imin,Jmin+1))*xfac)**2+          &
     &          ((Ygrd(Imin,Jmin)-Ygrd(Imin,Jmin+1))*yfac)**2
            phi=ASIN((diag2-aa2-bb2)/(2.0_r8*SQRT(aa2*bb2)))
!
!  Transform float position into curvilinear coordinates. Assume the
!  cell is rectanglar, for now.
!
            ang=angler(Imin,Jmin)
            dx=xpp*COS(ang)+ypp*SIN(ang)
            dy=ypp*COS(ang)-xpp*SIN(ang)
!
!  Correct for parallelogram.
!
            dx=dx+dy*TAN(phi)
            dy=dy/COS(phi)
!
!  Scale with cell side lengths to translate into cell indices.
!
            dx=MIN(MAX(0.0_r8,dx/SQRT(aa2)),1.0_r8)
            dy=MIN(MAX(0.0_r8,dy/SQRT(bb2)),1.0_r8)
            Ipos(mp,np)=REAL(Imin,r8)+dx
            Jpos(mp,np)=REAL(Jmin,r8)+dy
          END IF
#endif /* OFFLINE_FLOATS */
        END DO
      END DO

      RETURN
      END SUBROUTINE hindices

      LOGICAL FUNCTION try_range (ng, LBi, UBi, LBj, UBj, Xgrd, Ygrd,   &
     &                            Imin, Imax, Jmin, Jmax, Xo, Yo)
!
!=======================================================================
!                                                                      !
!  Given a grided domain with matrix coordinates Xgrd and Ygrd, this   !
!  function finds if the point (Xo,Yo)  is inside the box defined by   !
!  the requested corners (Imin,Jmin) and (Imax,Jmax). It will return   !
!  logical switch  try_range=.TRUE.  if (Xo,Yo) is inside, otherwise   !
!  it will return false.                                               !
!                                                                      !
!  Calls:   inside                                                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
      integer, intent(in) :: Imin, Imax, Jmin, Jmax

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Xgrd(LBi:,LBj:)
      real(r8), intent(in) :: Ygrd(LBi:,LBj:)
#else
      real(r8), intent(in) :: Xgrd(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Ygrd(LBi:UBi,LBj:UBj)
#endif

      real(r8), intent(in) :: Xo, Yo
!
!  Local variable declarations.
!
      integer ::  Nb, i, j, shft, ic

      real(r8), dimension(2*(Jmax-Jmin+Imax-Imin)+1) :: Xb, Yb
!
!-----------------------------------------------------------------------
!  Define closed polygon.
!-----------------------------------------------------------------------
!
!  Note that the last point (Xb(Nb),Yb(Nb)) does not repeat first
!  point (Xb(1),Yb(1)).  Instead, in function inside, it is implied
!  that the closing segment is (Xb(Nb),Yb(Nb))-->(Xb(1),Yb(1)). In
!  fact, function inside sets Xb(Nb+1)=Xb(1) and Yb(Nb+1)=Yb(1).
!
      Nb=2*(Jmax-Jmin+Imax-Imin)
      shft=1-Imin
      DO i=Imin,Imax-1
        Xb(i+shft)=Xgrd(i,Jmin)
        Yb(i+shft)=Ygrd(i,Jmin)
      END DO
      shft=1-Jmin+Imax-Imin
      DO j=Jmin,Jmax-1
        Xb(j+shft)=Xgrd(Imax,j)
        Yb(j+shft)=Ygrd(Imax,j)
      END DO
      shft=1+Jmax-Jmin+2*Imax-Imin
      DO i=Imax,Imin+1,-1
        Xb(shft-i)=Xgrd(i,Jmax)
        Yb(shft-i)=Ygrd(i,Jmax)
      END DO
      shft=1+2*Jmax-Jmin+2*(Imax-Imin)
      DO j=Jmax,Jmin+1,-1
        Xb(shft-j)=Xgrd(Imin,j)
        Yb(shft-j)=Ygrd(Imin,j)
      END DO
!
!-----------------------------------------------------------------------
!  Check if point (Xo,Yo) is inside of the defined polygon.
!-----------------------------------------------------------------------
!
      try_range=inside(Nb, Xb, Yb, Xo, Yo)
      RETURN
      END FUNCTION try_range

#ifdef OFFLINE_FLOATS

      LOGICAL FUNCTION inside (Nb, Xb, Yb, Xo, Yo)
!
!=======================================================================
!                                                                      !
!  Given the vectors Xb and Yb of size Nb, defining the coordinates    !
!  of a closed polygon,  this function find if the point (Xo,Yo) is    !
!  inside the polygon.  If the point  (Xo,Yo)  falls exactly on the    !
!  boundary of the polygon, it is still considered inside.             !
!                                                                      !
!  This algorithm does not rely on the setting of  Xb(Nb)=Xb(1) and    !
!  Yb(Nb)=Yb(1).  Instead, it assumes that the last closing segment    !
!  is (Xb(Nb),Yb(Nb)) --> (Xb(1),Yb(1)).                               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Reid, C., 1969: A long way from Euclid. Oceanography EMR,         !
!      page 174.                                                       !
!                                                                      !
!  Algorithm:                                                          !
!                                                                      !
!  The decision whether the point is  inside or outside the polygon    !
!  is done by counting the number of crossings from the ray (Xo,Yo)    !
!  to (Xo,-infinity), hereafter called meridian, by the boundary of    !
!  the polygon.  In this counting procedure,  a crossing is counted    !
!  as +2 if the crossing happens from "left to right" or -2 if from    !
!  "right to left". If the counting adds up to zero, then the point    !
!  is outside.  Otherwise,  it is either inside or on the boundary.    !
!                                                                      !
!  This routine is a modified version of the Reid (1969) algorithm,    !
!  where all crossings were counted as positive and the decision is    !
!  made  based on  whether the  number of crossings is even or odd.    !
!  This new algorithm may produce different results  in cases where    !
!  Xo accidentally coinsides with one of the (Xb(k),k=1:Nb) points.    !
!  In this case, the crossing is counted here as +1 or -1 depending    !
!  of the sign of (Xb(k+1)-Xb(k)).  Crossings  are  not  counted if    !
!  Xo=Xb(k)=Xb(k+1).  Therefore, if Xo=Xb(k0) and Yo>Yb(k0), and if    !
!  Xb(k0-1) < Xb(k0) < Xb(k0+1),  the crossing is counted twice but    !
!  with weight +1 (for segments with k=k0-1 and k=k0). Similarly if    !
!  Xb(k0-1) > Xb(k0) > Xb(k0+1), the crossing is counted twice with    !
!  weight -1 each time.  If,  on the other hand,  the meridian only    !
!  touches the boundary, that is, for example, Xb(k0-1) < Xb(k0)=Xo    !
!  and Xb(k0+1) < Xb(k0)=Xo, then the crossing is counted as +1 for    !
!  segment k=k0-1 and -1 for segment k=k0, resulting in no crossing.   !
!                                                                      !
!  Note 1: (Explanation of the logical condition)                      !
!                                                                      !
!  Suppose  that there exist two points  (x1,y1)=(Xb(k),Yb(k))  and    !
!  (x2,y2)=(Xb(k+1),Yb(k+1)),  such that,  either (x1 < Xo < x2) or    !
!  (x1 > Xo > x2).  Therefore, meridian x=Xo intersects the segment    !
!  (x1,y1) -> (x2,x2) and the ordinate of the point of intersection    !
!  is:                                                                 !
!                                                                      !
!                 y1*(x2-Xo) + y2*(Xo-x1)                              !
!             y = -----------------------                              !
!                          x2-x1                                       !
!                                                                      !
!  The mathematical statement that point  (Xo,Yo)  either coinsides    !
!  with the point of intersection or lies to the north (Yo>=y) from    !
!  it is, therefore, equivalent to the statement:                      !
!                                                                      !
!         Yo*(x2-x1) >= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 > 0      !
!  or                                                                  !
!         Yo*(x2-x1) <= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 < 0      !
!                                                                      !
!  which, after noting that  Yo*(x2-x1) = Yo*(x2-Xo + Xo-x1) may be    !
!  rewritten as:                                                       !
!                                                                      !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) >= 0,   if   x2-x1 > 0      !
!  or                                                                  !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) <= 0,   if   x2-x1 < 0      !
!                                                                      !
!  and both versions can be merged into  essentially  the condition    !
!  that (Yo-y1)*(x2-Xo)+(Yo-y2)*(Xo-x1) has the same sign as x2-x1.    !
!  That is, the product of these two must be positive or zero.         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: Nb

      real(r8), intent(in) :: Xo, Yo

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: Xb(:), Yb(:)
#else
      real(r8), intent(inout) :: Xb(Nb+1), Yb(Nb+1)
#endif
!
!  Local variable declarations.
!
      integer, parameter :: Nstep =128

      integer :: crossings, i, inc, k, kk, nc

      integer, dimension(Nstep) :: Sindex

      real(r8) :: dx1, dx2, dxy
!
!-----------------------------------------------------------------------
!  Find intersections.
!-----------------------------------------------------------------------
!
!  Set crossings counter and close the contour of the polygon.
!
      crossings=0
      Xb(Nb+1)=Xb(1)
      Yb(Nb+1)=Yb(1)
!
!  The search is optimized.  First select the indices of segments
!  where Xb(k) is different from Xb(k+1) and Xo falls between them.
!  Then, further investigate these segments in a separate loop.
!  Doing it in two stages takes less time because the first loop is
!  pipelined.
!
      DO kk=0,Nb-1,Nstep
        nc=0
        DO k=kk+1,MIN(kk+Nstep,Nb)
          !Fix missing search step, added by Patrick Drake, Feb.  2007
          !Does meridian coincide with a vertical side of the polygon?
          !If so, count it as inside with 1 crossing
          IF ( ((Xb(k).eq.Xo).and.(Xb(k+1).eq.Xo)).and.                 &
     &        (((Yb(k+1)-Yo)*(Yo-Yb(k))) .gt.0.0_r8) ) THEN
             crossings = 1
             goto 10
          END IF
          !End PTD Fix
          !Only investigate points for crossings if they fall
          !between Xb(k+1) and Xb(k) 
          IF (((Xb(k+1)-Xo)*(Xo-Xb(k)).ge.0.0_r8).and.                  &
     &        (Xb(k).ne.Xb(k+1))) THEN
            nc=nc+1
            Sindex(nc)=k
          END IF
        END DO
        !Now investigate crossings
        DO i=1,nc
          k=Sindex(i)
          IF (Xb(k).ne.Xb(k+1)) THEN
            dx1=Xo-Xb(k)
            dx2=Xb(k+1)-Xo
            dxy=dx2*(Yo-Yb(k))-dx1*(Yb(k+1)-Yo)
            inc=0

            !This is the case where (Xo,Yo) is on a vertex
            IF ((Xb(k).eq.Xo).and.(Yb(k).eq.Yo)) THEN
              crossings=1
              goto 10
            !This is the case where the meridian intersects a vertex
            ELSE IF (((dx1.eq.0.0_r8).and.(Yo.gt.Yb(k))) .or.           &
     &              ((dx2.eq.0.0_r8).and.(Yo.gt.Yb(k+1))) ) THEN
              inc=1
            !This is the test for whether the meridian is "north" of segment
            ELSE IF ((dx1*dx2.gt.0.0_r8).and.                           &
     &              ((Xb(k+1)-Xb(k))*dxy.ge.0.0_r8)) THEN  ! see note 1
              inc=2
            END IF
            !Account for the Meridian crossing polygon twice, so point is
            !outside and net crossings equal 0
            IF (Xb(k+1).gt.Xb(k)) THEN
              crossings=crossings+inc
            ELSE
              crossings=crossings-inc
            END IF
          ELSE
            print *, '  Not sure if this should happen.  Check logic.'
            print *, 'CAE in interpolate.F . Exiting.'
            stop
          END IF
        END DO
      END DO
!
!  Determine if point (Xo,Yo) is inside of closed polygon.
!
  10  IF (crossings.eq.0) THEN
        inside=.FALSE.
      ELSE
        inside=.TRUE.
      END IF
      !END IF !End HPB test for vertex
      RETURN
      END FUNCTION inside

       subroutine quadbisect(mp,x0,y0,x,y,xfrac,yfrac)
!=======================================================================
! Copyright (c) ?
!=======================================================================
! QUADBISECT determines the fractional loaction (0-1)
! of a point within a quadrilateral by iterative bisection
! in both axis
!
! x0 and y0 are coordinates of the point
! x and y are the coordinates of the vertices of the quadrilateral
! xfrac and yfrac are the fractional coordinates (0..1)
!
! routine calls INSIDE, to determine whether a point is within a
! quadrilateral
!
! Harold P. Batchelder, 11 July 2006,
! Modified by Patrick Drake, Feb. 2007
!========================================================================

      USE mod_iounits

      real(r8), intent(in) :: x0, y0
      real(r8), intent(inout) :: x(5),y(5) !add an extra element that inside will use
      real(r8), intent(out) :: xfrac, yfrac
      real(r8) :: tol, xmin, xmax, ymin, ymax
      real(r8) :: mx(9), my(9)
      real(r8) :: xur(5), yur(5), xlr(5), ylr(5)
      real(r8) :: xul(5), yul(5), xll(5), yll(5)
      real(r8) :: xcn, ycn, xd, yd, eps
      logical :: done
      logical :: result
      integer :: nbv, count, quad, try, mp

! some setup

       tol = 0.000001_r8
       xmin = 0.0_r8
       xmax = 1.0_r8
       ymin = 0.0_r8
       ymax = 1.0_r8
       done = .FALSE.
       count = 0
       nbv = 4
       eps = epsilon(x0)

      try = 0
      DO WHILE (.not. done)
        count = count + 1
        mx(1) = x(1)
        mx(3) = x(2)
        mx(5) = x(3)
        mx(7) = x(4)
        my(1) = y(1)
        my(3) = y(2)
        my(5) = y(3)
        my(7) = y(4)
!
! find midpoints of edges
        mx(2) = 0.5_r8 * (mx(1) + mx(3))
        mx(4) = 0.5_r8 * (mx(3) + mx(5))
        mx(6) = 0.5_r8 * (mx(5) + mx(7))
        mx(8) = 0.5_r8 * (mx(7) + mx(1))
        my(2) = 0.5_r8 * (my(1) + my(3))
        my(4) = 0.5_r8 * (my(3) + my(5))
        my(6) = 0.5_r8 * (my(5) + my(7))
        my(8) = 0.5_r8 * (my(7) + my(1))

! find centroid 
        mx(9) = 0.25_r8 * (mx(1)+mx(3)+mx(5)+mx(7));
        my(9) = 0.25_r8 * (my(1)+my(3)+my(5)+my(7));
!
! Define upper left quadrant
        xul(1:4) = (/mx(1), mx(2), mx(9), mx(8)/)
        yul(1:4) = (/my(1), my(2), my(9), my(8)/)
!
! Define upper right quadrant
        xur(1:4) = (/mx(2), mx(3), mx(4), mx(9)/)
        yur(1:4) = (/my(2), my(3), my(4), my(9)/)
!
! Define lower right quadrant
        xlr(1:4) = (/mx(9), mx(4), mx(5), mx(6)/)
        ylr(1:4) = (/my(9), my(4), my(5), my(6)/)
!
! Define lower left quadrant
        xll(1:4) = (/mx(8), mx(9), mx(6), mx(7)/)
        yll(1:4) = (/my(8), my(9), my(6), my(7)/)

!cae 11/07 moving this outside loop        try = 0
        IF (inside(nbv, xul, yul, x0, y0)) THEN
!
! Try upper left quadrant  
          quad = 1
          xmax = 0.5_r8 * (xmin + xmax)
          ymin = 0.5_r8 * (ymin + ymax)
          x = xul
          y = yul
          try = 0
        ELSE IF (inside(nbv, xur, yur, x0, y0)) THEN
!
! Try upper right quadrant
          quad = 2
          xmin = 0.5_r8 * (xmin + xmax)
          ymin = 0.5_r8 * (ymin + ymax)
          x = xur
          y = yur
          try = 0
        ELSE IF (inside(nbv, xlr, ylr, x0, y0)) THEN
!
! Try lower right quadrant
          quad = 3
          xmin = 0.5_r8 * (xmin + xmax)
          ymax = 0.5_r8 * (ymin + ymax)
          x = xlr
          y = ylr
          try = 0
        ELSE IF (inside(nbv, xll, yll, x0, y0)) THEN
!
! Try lower left quadrant
          quad = 4 
          xmax = 0.5_r8 * (xmin + xmax)
          ymax = 0.5_r8 * (ymin + ymax)
          x = xll
          y = yll
          try = 0
        ELSE
!
! Point lost inside box, probably due to machine error,
! so bloat box a little, one side at a time
          try = try + 1
          IF (try.eq.1) THEN
            WRITE(stdout,*) 'Bloated right side of bounding box ',      &
     &                      'during grid coords. Search for float ', mp
            x(2) = x(2) + 2*eps
            x(3) = x(3) + 2*eps
          ELSE IF (try.eq.2) THEN 
!
! Unbloat right side
            x(2) = x(2) - 2*eps
            x(3) = x(3) - 2*eps
            WRITE(stdout,*) 'Bloated bottom of bounding box ',          &
     &                      'during grid coords. Search for float ', mp
            y(3) = y(3) - 2*eps
            y(4) = y(4) - 2*eps
          ELSE IF (try.eq.3) THEN 
!
! Unbloat bottom
            y(3) = y(3) + 2*eps
            y(4) = y(4) + 2*eps
            WRITE(stdout,*) 'Bloated left side of bounding box ',       &
     &                      'during grid coords. Search for float ', mp
            x(4) = x(4) - 2*eps
            x(1) = x(1) - 2*eps
          ELSE IF (try.eq.4) THEN 
!
! Unbloat left side
            x(4) = x(4) + 2*eps
            x(1) = x(1) + 2*eps
            WRITE(stdout,*) 'Bloated top of bounding box during ',      &
     &                      'grid coords. Search for float ', mp
            y(1) = y(1) + 2*eps
            y(2) = y(2) + 2*eps
          ELSE
!
! Float lost, break program
           WRITE(stdout,*) 'Float ', mp, ' lost during grid coords.  ', &
     &                      'Search, giving up and stopping program.'
           WRITE(stdout,*) 'Last bounding box:'
           WRITE(stdout,*) 'x1, x2, x3, x4: ', x(1:4)
           WRITE(stdout,*) 'y1, y2, y3, y4: ', y(1:4)
           WRITE(stdout,*) 'Float(x0,y0): ', x0, y0
           STOP
          END IF

        END IF

!       evaluate convergence using relative change of centroids
        xcn = SUM(x(1:4))/4.0_r8 ! centroid of new subregion
        ycn = SUM(y(1:4))/4.0_r8
        xd = ABS((mx(9)-xcn))
        yd = ABS((my(9)-ycn))
        IF ((try.eq.0) .AND. ((xd .LT. tol)) .AND. (yd .LT. tol))       &
     &      done = .TRUE.

      END DO

      xfrac = 0.5_r8 * (xmin+xmax)
      yfrac = 0.5_r8 * (ymin+ymax)

     END SUBROUTINE quadbisect

#else /* OFFLINE_FLOATS */

      LOGICAL FUNCTION inside (Nb, Xb, Yb, Xo, Yo)
!
!=======================================================================
!                                                                      !
!  Given the vectors Xb and Yb of size Nb, defining the coordinates    !
!  of a closed polygon,  this function find if the point (Xo,Yo) is    !
!  inside the polygon.  If the point  (Xo,Yo)  falls exactly on the    !
!  boundary of the polygon, it still considered inside.                !
!                                                                      !
!  This algorithm does not rely on the setting of  Xb(Nb)=Xb(1) and    !
!  Yb(Nb)=Yb(1).  Instead, it assumes that the last closing segment    !
!  is (Xb(Nb),Yb(Nb)) --> (Xb(1),Yb(1)).                               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Reid, C., 1969: A long way from Euclid. Oceanography EMR,         !
!      page 174.                                                       !
!                                                                      !
!  Algorithm:                                                          !
!                                                                      !
!  The decision whether the point is  inside or outside the polygon    !
!  is done by counting the number of crossings from the ray (Xo,Yo)    !
!  to (Xo,-infinity), hereafter called meridian, by the boundary of    !
!  the polygon.  In this counting procedure,  a crossing is counted    !
!  as +2 if the crossing happens from "left to right" or -2 if from    !
!  "right to left". If the counting adds up to zero, then the point    !
!  is outside.  Otherwise,  it is either inside or on the boundary.    !
!                                                                      !
!  This routine is a modified version of the Reid (1969) algorithm,    !
!  where all crossings were counted as positive and the decision is    !
!  made  based on  whether the  number of crossings is even or odd.    !
!  This new algorithm may produce different results  in cases where    !
!  Xo accidentally coinsides with one of the (Xb(k),k=1:Nb) points.    !
!  In this case, the crossing is counted here as +1 or -1 depending    !
!  of the sign of (Xb(k+1)-Xb(k)).  Crossings  are  not  counted if    !
!  Xo=Xb(k)=Xb(k+1).  Therefore, if Xo=Xb(k0) and Yo>Yb(k0), and if    !
!  Xb(k0-1) < Xb(k0) < Xb(k0+1),  the crossing is counted twice but    !
!  with weight +1 (for segments with k=k0-1 and k=k0). Similarly if    !
!  Xb(k0-1) > Xb(k0) > Xb(k0+1), the crossing is counted twice with    !
!  weight -1 each time.  If,  on the other hand,  the meridian only    !
!  touches the boundary, that is, for example, Xb(k0-1) < Xb(k0)=Xo    !
!  and Xb(k0+1) < Xb(k0)=Xo, then the crossing is counted as +1 for    !
!  segment k=k0-1 and -1 for segment k=k0, resulting in no crossing.   !
!                                                                      !
!  Note 1: (Explanation of the logical condition)                      !
!                                                                      !
!  Suppose  that there exist two points  (x1,y1)=(Xb(k),Yb(k))  and    !
!  (x2,y2)=(Xb(k+1),Yb(k+1)),  such that,  either (x1 < Xo < x2) or    !
!  (x1 > Xo > x2).  Therefore, meridian x=Xo intersects the segment    !
!  (x1,y1) -> (x2,x2) and the ordinate of the point of intersection    !
!  is:                                                                 !
!                                                                      !
!                 y1*(x2-Xo) + y2*(Xo-x1)                              !
!             y = -----------------------                              !
!                          x2-x1                                       !
!                                                                      !
!  The mathematical statement that point  (Xo,Yo)  either coinsides    !
!  with the point of intersection or lies to the north (Yo>=y) from    !
!  it is, therefore, equivalent to the statement:                      !
!                                                                      !
!         Yo*(x2-x1) >= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 > 0      !
!  or                                                                  !
!         Yo*(x2-x1) <= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 < 0      !
!                                                                      !
!  which, after noting that  Yo*(x2-x1) = Yo*(x2-Xo + Xo-x1) may be    !
!  rewritten as:                                                       !
!                                                                      !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) >= 0,   if   x2-x1 > 0      !
!  or                                                                  !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) <= 0,   if   x2-x1 < 0      !
!                                                                      !
!  and both versions can be merged into  essentially  the condition    !
!  that (Yo-y1)*(x2-Xo)+(Yo-y2)*(Xo-x1) has the same sign as x2-x1.    !
!  That is, the product of these two must be positive or zero.         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: Nb

      real(r8), intent(in) :: Xo, Yo

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: Xb(:), Yb(:)
#else
      real(r8), intent(inout) :: Xb(Nb+1), Yb(Nb+1)
#endif
!
!  Local variable declarations.
!
      integer, parameter :: Nstep =128

      integer :: crossings, i, inc, k, kk, nc

      integer, dimension(Nstep) :: Sindex

      real(r8) :: dx1, dx2, dxy
!
!-----------------------------------------------------------------------
!  Find intersections.
!-----------------------------------------------------------------------
!
!  Set crossings counter and close the contour of the polygon.
!
      crossings=0
      Xb(Nb+1)=Xb(1)
      Yb(Nb+1)=Yb(1)
!
!  The search is optimized.  First select the indices of segments
!  where Xb(k) is different from Xb(k+1) and Xo falls between them.
!  Then, further investigate these segments in a separate loop.
!  Doing it in two stages takes less time because the first loop is
!  pipelined.
!
      DO kk=0,Nb-1,Nstep
        nc=0
        DO k=kk+1,MIN(kk+Nstep,Nb)
          IF (((Xb(k+1)-Xo)*(Xo-Xb(k)).ge.0.0_r8).and.                  &
     &        (Xb(k).ne.Xb(k+1))) THEN
            nc=nc+1
            Sindex(nc)=k
          END IF
        END DO
        DO i=1,nc
          k=Sindex(i)
          IF (Xb(k).ne.Xb(k+1)) THEN
            dx1=Xo-Xb(k)
            dx2=Xb(k+1)-Xo
            dxy=dx2*(Yo-Yb(k))-dx1*(Yb(k+1)-Yo)
            inc=0
            IF ((Xb(k).eq.Xo).and.(Yb(k).eq.Yo)) THEN
              crossings=1
              goto 10
            ELSE IF (((dx1.eq.0.0_r8).and.(Yo.ge.Yb(k  ))).or.          &
     &              ((dx2.eq.0.0_r8).and.(Yo.ge.Yb(k+1)))) THEN
              inc=1
            ELSE IF ((dx1*dx2.gt.0.0_r8).and.                           &
     &              ((Xb(k+1)-Xb(k))*dxy.ge.0.0_r8)) THEN  ! see note 1
              inc=2
            END IF
            IF (Xb(k+1).gt.Xb(k)) THEN
              crossings=crossings+inc
            ELSE
              crossings=crossings-inc
            END IF
          END IF
        END DO
      END DO
!
!  Determine if point (Xo,Yo) is inside of closed polygon.
!
  10  IF (crossings.eq.0) THEN
        inside=.FALSE.
      ELSE
        inside=.TRUE.
      END IF
      RETURN
      END FUNCTION inside
#endif /* OFFLINE_FLOATS */

      END MODULE interpolate_mod
