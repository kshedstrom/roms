#include "cppdefs.h"
      MODULE interp_fish_mod
#if defined NONLINEAR && defined NEMURO_SAN
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group         Mark Hadfield   !
!    Licensed under a MIT/X style license             John M. Klinck   !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine interpolates requested field at the float trajectory   !
!  locations.                                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     LBk        K-dimension Lower bound.                              !
!     UBk        K-dimension Upper bound.                              !
!     Lstr       Starting float index to process.                      !
!     Lend       Ending   float index to process.                      !
!     itime      Fish time level to process.                           !
!     ifield     ID of field to compute.                               !
!     gtype      Grid type. If negative, interpolate fish slopes.      !
!     maskit     Should the field be masked? Ignored if Land/Sea       !
!                 masking is not active.                               !
!     pm         Inverse grid spacing (1/m) in the XI-direction.       !
!     pn         Inverse grid spacing (1/m) in the ETA-direction.      !
!     Amask      Field Land/Sea mask.                                  !
!     A          Field to interpolate from.                            !
!     fishthread   Float parallel thread bounded switch.                 !
!     bounded    Float grid bounded status switch.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     track      Interpolated field: track(ifield,itime,:).            !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: interp_fish

      CONTAINS
!
!***********************************************************************
      SUBROUTINE interp_fish (ng, LBi, UBi, LBj, UBj, LBk, UBk,         &
     &                          Lstr, Lend, itime, ifield, numflt,      &
     &                          gtype, maskit,                          &
     &                          pm, pn,                                 &
# ifdef MASKING
     &                          Amask,                                  &
# endif
     &                          A, fishthread, bounded, track, bioenergy)
!***********************************************************************
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Lstr, Lend, itime, ifield, gtype, numflt

      logical, intent(in) :: maskit
      logical, intent(in) :: fishthread(Lstr:Lend)
      logical, intent(in) :: bounded(numflt)

      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk)

      real(r8), intent(in) :: track(NFV(ng),0:NFT,numflt)
      real(r8), intent(inout) :: bioenergy(NFishV(ng),numflt)
!
!  Local variable declarations.
!
      logical :: Irvar, Iuvar, Jrvar, Jvvar, Krvar, Kwvar, Lmask
      logical :: halo

      integer :: Ir, Iu, Jr, Jv, Kr, Kw, l
      integer :: i1, i2, j1, j2, k1, k2, khm, khp, vtype

      real(r8) :: p1, p2, q1, q2, r1, r2
      real(r8) :: s111, s211, s121, s221, s112, s212, s122, s222

# ifdef MASKING
      integer :: Irn, Irnm1, Irnp1, Jrn, Jrnm1, Jrnp1

      real(r8) :: cff1, cff2
# endif
!
!-----------------------------------------------------------------------
!  Initialize various internal variables.
!-----------------------------------------------------------------------
!
!  Determine variable type switches.
!
      vtype=ABS(gtype)
      Irvar=(vtype.eq.r2dvar).or.(vtype.eq.r3dvar).or.                  &
     &      (vtype.eq.v2dvar).or.(vtype.eq.v3dvar).or.                  &
     &      (vtype.eq.w3dvar)
      Jrvar=(vtype.eq.r2dvar).or.(vtype.eq.r3dvar).or.                  &
     &      (vtype.eq.u2dvar).or.(vtype.eq.u3dvar).or.                  &
     &      (vtype.eq.w3dvar)
      Iuvar=(vtype.eq.u2dvar).or.(vtype.eq.u3dvar)
      Jvvar=(vtype.eq.v2dvar).or.(vtype.eq.v3dvar)
      Krvar=(vtype.eq.r3dvar).or.(vtype.eq.u3dvar).or.                  &
     &      (vtype.eq.v3dvar)
      Kwvar=(vtype.eq.w3dvar)
!
!  Determine whether to allow for masking in horizontal interpolation.
!
# ifdef MASKING
      Lmask=maskit
# else
      Lmask=.FALSE.
# endif
!
!  If not interpolating slope, set multipliers to 1.
!
      IF (gtype.ge.0) THEN
        s111=1.0_r8
        s121=1.0_r8
        s211=1.0_r8
        s221=1.0_r8
        s112=1.0_r8
        s122=1.0_r8
        s212=1.0_r8
        s222=1.0_r8
      END IF
!
!-----------------------------------------------------------------------
!  Loop through fish.
!-----------------------------------------------------------------------
!
      DO l=Lstr,Lend
        IF (fishthread(l)) THEN
          IF (.not.bounded(l)) THEN
            bioenergy(ifield,l)=spval
          ELSE
!
!  Calculate indices and weights for vertical interpolation, if any.
!
            IF (Krvar) THEN
              Kr=INT(track(izgrd,itime,l)+0.5_r8)
              k1=MIN(MAX(Kr  ,1),N(ng))
              k2=MIN(MAX(Kr+1,1),N(ng))
              r2=REAL(k2-k1,r8)*(track(izgrd,itime,l)+                  &
     &                           0.5_r8-REAL(k1,r8))
            ELSE IF (Kwvar) THEN
              Kw=INT(track(izgrd,itime,l))
              k1=MIN(MAX(Kw  ,0),N(ng))
              k2=MIN(MAX(Kw+1,0),N(ng))
              r2=REAL(k2-k1,r8)*(track(izgrd,itime,l)-REAL(k1,r8))
            ELSE
              k1=1
              k2=1
              r2=0.0_r8
            END IF
            r1=1.0_r8-r2
!
!-----------------------------------------------------------------------
!  Interpolation on RHO-columns.
!-----------------------------------------------------------------------
!
            IF (Irvar.and.Jrvar) THEN
!
              Ir=INT(track(ixgrd,itime,l))
              Jr=INT(track(iygrd,itime,l))
!
              i1=MIN(MAX(Ir  ,0),Lm(ng)+1)
              i2=MIN(MAX(Ir+1,1),Lm(ng)+1)
              j1=MIN(MAX(Jr  ,0),Mm(ng)+1)
              j2=MIN(MAX(Jr+1,1),Mm(ng)+1)
!
              p2=REAL(i2-i1,r8)*(track(ixgrd,itime,l)-REAL(i1,r8))
              q2=REAL(j2-j1,r8)*(track(iygrd,itime,l)-REAL(j1,r8))
              p1=1.0_r8-p2
              q1=1.0_r8-q2
!
              IF (Lmask) THEN
# ifdef MASKING
                cff1=p1*q1*r1*Amask(i1,j1)*s111*A(i1,j1,k1)+            &
     &               p2*q1*r1*Amask(i2,j1)*s211*A(i2,j1,k1)+            &
     &               p1*q2*r1*Amask(i1,j2)*s121*A(i1,j2,k1)+            &
     &               p2*q2*r1*Amask(i2,j2)*s221*A(i2,j2,k1)+            &
     &               p1*q1*r2*Amask(i1,j1)*s112*A(i1,j1,k2)+            &
     &               p2*q1*r2*Amask(i2,j1)*s212*A(i2,j1,k2)+            &
     &               p1*q2*r2*Amask(i1,j2)*s122*A(i1,j2,k2)+            &
     &               p2*q2*r2*Amask(i2,j2)*s222*A(i2,j2,k2)
                cff2=p1*q1*r1*Amask(i1,j1)+                             &
     &               p2*q1*r1*Amask(i2,j1)+                             &
     &               p1*q2*r1*Amask(i1,j2)+                             &
     &               p2*q2*r1*Amask(i2,j2)+                             &
     &               p1*q1*r2*Amask(i1,j1)+                             &
     &               p2*q1*r2*Amask(i2,j1)+                             &
     &               p1*q2*r2*Amask(i1,j2)+                             &
     &               p2*q2*r2*Amask(i2,j2)
                IF (cff2.gt.0.0_r8) THEN
                  bioenergy(ifield,l)=cff1/cff2
                ELSE
                  bioenergy(ifield,l)=0.0_r8
                END IF
# endif
              ELSE
                bioenergy(ifield,l)=p1*q1*r1*s111*A(i1,j1,k1)+          &
     &                              p2*q1*r1*s211*A(i2,j1,k1)+          &
     &                              p1*q2*r1*s121*A(i1,j2,k1)+          &
     &                              p2*q2*r1*s221*A(i2,j2,k1)+          &
     &                              p1*q1*r2*s112*A(i1,j1,k2)+          &
     &                              p2*q1*r2*s212*A(i2,j1,k2)+          &
     &                              p1*q2*r2*s122*A(i1,j2,k2)+          &
     &                              p2*q2*r2*s222*A(i2,j2,k2)
              END IF
!
!-----------------------------------------------------------------------
!  Interpolation at horizontal velocity points.
!-----------------------------------------------------------------------
!
            ELSE
              Ir=INT(track(ixgrd,itime,l))
              Jr=INT(track(iygrd,itime,l))
              Iu=INT(track(ixgrd,itime,l)+0.5_r8)
              Jv=INT(track(iygrd,itime,l)+0.5_r8)
!
              halo=.FALSE.
# ifdef MASKING
!
!  Is the float inside the halo (i.e. inside a masked grid cell or
!  within 0.5*delta of one)?  Depending on which portion of the cell
!  the float is in,  we look for adjacent and diagonally adjacent
!  masked cells.
!
!  Note that the halo-checking code is evaluated every time
!  interp_fish is called for a velocity variable with masking
!  activated. This is at least twice as often as it needs to be
!  called, as the float does not move between interpolating U and V.
!
!  The special-case code for periodic boundary conditions may be
!  unnecessary
!
              IF (Lmask) THEN
                Irn=NINT(track(ixgrd,itime,l))
                Jrn=NINT(track(iygrd,itime,l))
#  ifdef EW_PERIODIC
                IF (Irn.ge.Lm(ng)) THEN
                  Irnp1=Irn+1-Lm(ng)
                ELSE
                  Irnp1=Irn+1
                END IF
                IF (Irn.le.1) THEN
                  Irnm1=Irn-1+Lm(ng)
                ELSE
                  Irnm1=Irn-1
                END IF
#  else
                Irnm1=Irn-1
                Irnp1=Irn+1
#  endif
#  ifdef NS_PERIODIC
                IF (Jrn.ge.Mm(ng)) THEN
                  Jrnp1=Jrn+1-Mm(ng)
                ELSE
                  Jrnp1=Jrn+1
                END IF
                IF (Jrn.le.1) THEN
                  Jrnm1=Jrn-1+Mm(ng)
                ELSE
                  Jrnm1=Jrn-1
                END IF
#  else
                Jrnm1=Jrn-1
                Jrnp1=Jrn+1
#  endif
                IF (Amask(Irn,Jrn).lt.0.5_r8) THEN
                  halo=.TRUE.
                ELSE IF ((Ir.lt.Irn).and.                               &
     &                   (Amask(Irn-1,Jrn).lt.0.5_r8)) THEN
                  halo=.TRUE.
                ELSE IF ((Ir.eq.Irn).and.                               &
     &                   (Amask(Irn+1,Jrn).lt.0.5_r8)) THEN
                  halo=.TRUE.
                ELSE IF ((Jr.lt.Jrn).and.                               &
     &                   (Amask(Irn,Jrn-1).lt.0.5_r8)) THEN
                  halo=.TRUE.
                ELSE IF ((Jr.eq.Jrn).and.                               &
     &                   (Amask(Irn,Jrn+1).lt.0.5_r8)) THEN
                  halo=.TRUE.
                ELSE IF ((Ir.lt.Irn).and.(Jr.lt.Jrn).and.               &
     &                   (Amask(Irn-1,Jrn-1).lt.0.5_r8)) THEN
                  halo=.TRUE.
                ELSE IF ((Ir.eq.Irn).and.(Jr.lt.Jrn).and.               &
     &                   (Amask(Irn+1,Jrn-1).lt.0.5_r8)) THEN
                  halo=.TRUE.
                ELSE IF ((Ir.lt.Irn).and.(Jr.eq.Jrn).and.               &
     &                   (Amask(Irn-1,Jrn+1).lt.0.5_r8)) THEN
                  halo=.TRUE.
                ELSE IF ((Ir.eq.Irn).and.(Jr.eq.Jrn).and.               &
     &                   (Amask(Irn+1,Jrn+1).lt.0.5_r8)) THEN
                  halo=.TRUE.
                END IF
              END IF
# endif
!
!-----------------------------------------------------------------------
!  Interpolation at U-points.
!-----------------------------------------------------------------------
!
              IF (Iuvar) THEN
                IF (halo) THEN
# ifdef MASKING
!
!  Velocity interpolation inside the halo is linear in the parallel
!  direction and nearest-neighbour in the perpendicular direction.
!  The latter ensures that the perpendicular velocity is zero
!  everywhere on the perimeter of a masked cell.
!
                  i1=MIN(MAX(Iu  ,1),Lm(ng)+1)
                  i2=MIN(MAX(Iu+1,1),Lm(ng)+1)
                  j1=Jrn
!
                  p2=REAL(i2-i1,r8)*                                    &
     &               (track(ixgrd,itime,l)-REAL(i1,r8)+0.5_r8)
                  p1=1.0_r8-p2
                  q1=1.0_r8
!
                  IF (gtype.lt.0) THEN
                    s111=0.5_r8*(pm(i1-1,j1)+pm(i1,j1))
                    s211=0.5_r8*(pm(i2-1,j1)+pm(i2,j1))
                    s112=s111
                    s212=s112
                  END IF
!
                  bioenergy(ifield,l)=p1*q1*r1*s111*A(i1,j1,k1)+        &
     &                                p2*q1*r1*s211*A(i2,j1,k1)+        &
     &                                p1*q1*r2*s112*A(i1,j1,k2)+        &
     &                                p2*q1*r2*s212*A(i2,j1,k2)
# endif
                ELSE
!
!  Bilinear interpolation outside halo.
!
                  i1=MIN(MAX(Iu  ,1),Lm(ng)+1)
                  i2=MIN(MAX(Iu+1,1),Lm(ng)+1)
                  j1=MIN(MAX(Jr  ,0),Mm(ng)+1)
                  j2=MIN(MAX(Jr+1,0),Mm(ng)+1)
!
                  p2=REAL(i2-i1,r8)*                                    &
     &               (track(ixgrd,itime,l)-REAL(i1,r8)+0.5_r8)
                  q2=REAL(j2-j1,r8)*                                    &
     &               (track(iygrd,itime,l)-REAL(j1,r8))
                  p1=1.0_r8-p2
                  q1=1.0_r8-q2
!
                  IF (gtype.lt.0) THEN
                    s111=0.5_r8*(pm(i1-1,j1)+pm(i1,j1))
                    s211=0.5_r8*(pm(i2-1,j1)+pm(i2,j1))
                    s121=0.5_r8*(pm(i1-1,j2)+pm(i1,j2))
                    s221=0.5_r8*(pm(i2-1,j2)+pm(i2,j2))
                    s112=s111
                    s212=s112
                    s122=s121
                    s222=s221
                  END IF
!
                  bioenergy(ifield,l)=p1*q1*r1*s111*A(i1,j1,k1)+        &
     &                                p2*q1*r1*s211*A(i2,j1,k1)+        &
     &                                p1*q2*r1*s121*A(i1,j2,k1)+        &
     &                                p2*q2*r1*s221*A(i2,j2,k1)+        &
     &                                p1*q1*r2*s112*A(i1,j1,k2)+        &
     &                                p2*q1*r2*s212*A(i2,j1,k2)+        &
     &                                p1*q2*r2*s122*A(i1,j2,k2)+        &
     &                                p2*q2*r2*s222*A(i2,j2,k2)
                END IF
!
!-----------------------------------------------------------------------
!  Interpolation at V-points.
!-----------------------------------------------------------------------
!
              ELSE IF (Jvvar) THEN
                IF (halo) THEN
# ifdef MASKING
!
!  Velocity interpolation inside the halo is linear in the parallel
!  direction and nearest-neighbour in the perpendicular direction.
!  The latter ensures that the perpendicular velocity is zero
!  everywhere on the perimeter of a masked cell.
!
                  i1=Irn
                  j1=MIN(MAX(Jv  ,1),Mm(ng)+1)
                  j2=MIN(MAX(Jv+1,1),Mm(ng)+1)
!
                  q2=REAL(j2-j1,r8)*                                    &
     &               (track(iygrd,itime,l)-REAL(j1,r8)+0.5_r8)
                  p1=1.0_r8
                  q1=1.0_r8-q2
!
                  IF (gtype.lt.0) THEN
                    s111=0.5_r8*(pn(i1,j1-1)+pn(i1,j1))
                    s121=0.5_r8*(pn(i1,j2-1)+pn(i1,j2))
                    s112=s111
                    s122=s121
                  END IF
!
                  bioenergy(ifield,l)=p1*q1*r1*s111*A(i1,j1,k1)+        &
     &                                p1*q2*r1*s121*A(i1,j2,k1)+        &
     &                                p1*q1*r2*s112*A(i1,j1,k2)+        &
     &                                p1*q2*r2*s122*A(i1,j2,k2)
# endif
                ELSE
!
!  Bilinear interpolation outside halo.
!
                  i1=MIN(MAX(Ir  ,0),Lm(ng)+1)
                  i2=MIN(MAX(Ir+1,1),Lm(ng)+1)
                  j1=MIN(MAX(Jv  ,1),Mm(ng)+1)
                  j2=MIN(MAX(Jv+1,1),Mm(ng)+1)
!
                  p2=REAL(i2-i1,r8)*                                    &
     &               (track(ixgrd,itime,l)-REAL(i1,r8))
                  q2=REAL(j2-j1,r8)*                                    &
     &               (track(iygrd,itime,l)-REAL(j1,r8)+0.5_r8)
                  p1=1.0_r8-p2
                  q1=1.0_r8-q2
!
                  IF (gtype.lt.0) THEN
                    s111=0.5_r8*(pn(i1,j1-1)+pn(i1,j1))
                    s211=0.5_r8*(pn(i2,j1-1)+pn(i2,j1))
                    s121=0.5_r8*(pn(i1,j2-1)+pn(i1,j2))
                    s221=0.5_r8*(pn(i2,j2-1)+pn(i2,j2))
                    s112=s111
                    s212=s112
                    s122=s121
                    s222=s221
                  END IF
!
                  bioenergy(ifield,l)=p1*q1*r1*s111*A(i1,j1,k1)+        &
     &                                p2*q1*r1*s211*A(i2,j1,k1)+        &
     &                                p1*q2*r1*s121*A(i1,j2,k1)+        &
     &                                p2*q2*r1*s221*A(i2,j2,k1)+        &
     &                                p1*q1*r2*s112*A(i1,j1,k2)+        &
     &                                p2*q1*r2*s212*A(i2,j1,k2)+        &
     &                                p1*q2*r2*s122*A(i1,j2,k2)+        &
     &                                p2*q2*r2*s222*A(i2,j2,k2)
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE interp_fish
#endif
      END MODULE interp_fish_mod
