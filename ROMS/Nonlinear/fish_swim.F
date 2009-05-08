#include "cppdefs.h"
#if defined NONLINEAR && defined FLOATS && defined NEMURO_SAN
      SUBROUTINE fish_swim (ng, LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          lflt, itime, nnew,                      &
     &                          diZ, djZ, kZmax,                        &
     &                          pm, pn,                                 &
# ifdef SOLVE3D
     &                          Hz,                                     &
# endif
# ifdef MASKING
     &                          Amask,                                  &
# endif
     &                          A, track, bioenergy)
!
!svn $Id: interp_floats.F 294 2009-01-09 21:37:26Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group         Mark Hadfield   !
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
!     itime      Floats time level to process.                         !
!     ifield     ID of field to compute.                               !
!     gtype      Grid type. If negative, interpolate floats slopes.    !
!     maskit     Should the field be masked? Ignored if Land/Sea       !
!                 masking is not active.                               !
!     nudg       Vertical random walk term to be added to the field.   !
!     pm         Inverse grid spacing (1/m) in the XI-direction.       !
!     pn         Inverse grid spacing (1/m) in the ETA-direction.      !
!     Hz         Vertical thicknesses (m).                             !
!     Amask      Field Land/Sea mask.                                  !
!     A          Field to interpolate from.                            !
!     MyThread   Float parallel thread bounded switch.                 !
!     bounded    Float grid bounded status switch.                     !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     track      Interpolated field: track(ifield,itime,:).            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_biology
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: lflt, itime, nnew

      real(r8), intent(inout) :: diZ, djZ
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
# endif
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk,3,1:NT(ng))
      real(r8), intent(in) :: bioenergy(NFishV(ng),Nfloats(ng))

      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfloats(ng))
!
!  Local variable declarations.
!
      integer :: Ir, Jr, Kr, IrMax, JrMax, KrMax, kZmax
      integer :: i1, i2, j1, j2, i, j, k

      real(r8) :: SZooC, LZooC, PZooC, TZooC, ZooMax
      real(r8) :: a_A, b_A, Uswim, dxyG, fac
!
!
!-----------------------------------------------------------------------
!  Compute swimming velocity based on food concentration
!-----------------------------------------------------------------------
!
! Swimming velocity based on fish weight
! Uswim relationship is in cm/s, so multiply by 0.01 for m/s)
      a_A=3.9_r8
      b_A=0.13_r8
      Uswim=a_A*bioenergy(ifwwt,lflt)**b_A
      Uswim=0.01_r8*Uswim

      Ir=INT(track(ixgrd,itime,lflt)+0.5_r8)
      Jr=INT(track(iygrd,itime,lflt)+0.5_r8)
!      Kr=INT(track(izgrd,itime,lflt)+0.5_r8)
!
! Find k for cell with max food
!   Surface food concentration
# ifdef MASKING
      SZooC=ZSpref(ng)*Amask(Ir,Jr)*A(Ir,Jr,N(ng),nnew,iSzoo)
      LZooC=ZLpref(ng)*Amask(Ir,Jr)*A(Ir,Jr,N(ng),nnew,iLzoo)
      PZooC=ZPpref(ng)*Amask(Ir,Jr)*A(Ir,Jr,N(ng),nnew,iPzoo)
# else
      SZooC=ZSpref(ng)*A(Ir,Jr,N(ng),nnew,iSzoo)
      LZooC=ZLpref(ng)*A(Ir,Jr,N(ng),nnew,iLzoo)
      PZooC=ZPpref(ng)*A(Ir,Jr,N(ng),nnew,iPzoo)
# endif
      TZooC=SZooC+LZooC+PZooC
      ZooMax=TZooC
      KrMax=N(ng)
!   Subsurface food concentration
      DO k=N(ng)-1,1,-1
# ifdef MASKING
        SZooC=ZSpref(ng)*Amask(Ir,Jr)*A(Ir,Jr,k,nnew,iSzoo)
        LZooC=ZLpref(ng)*Amask(Ir,Jr)*A(Ir,Jr,k,nnew,iLzoo)
        PZooC=ZPpref(ng)*Amask(Ir,Jr)*A(Ir,Jr,k,nnew,iPzoo)
# else
        SZooC=ZSpref(ng)*A(Ir,Jr,k,nnew,iSzoo)
        LZooC=ZLpref(ng)*A(Ir,Jr,k,nnew,iLzoo)
        PZooC=ZPpref(ng)*A(Ir,Jr,k,nnew,iPzoo)
# endif
        TZooC=SZooC+LZooC+PZooC
        IF (TZooC.gt.ZooMax) THEN
          ZooMax=TZooC
          KrMax=k
        END IF
      END DO
!
! Find i,j for cell with max food
!      ZooMax=TZooC
      IrMax=Ir
      JrMax=Jr
!      KrMax=Kr
!
      DO j=Jr-1,Jr+1
        DO i=Ir-1,Ir+1
# ifdef MASKING
          SZooC=ZSpref(ng)*Amask(i,j)*A(i,j,KrMax,nnew,iSzoo)
          LZooC=ZLpref(ng)*Amask(i,j)*A(i,j,KrMax,nnew,iLzoo)
          PZooC=ZPpref(ng)*Amask(i,j)*A(i,j,KrMax,nnew,iPzoo)
# else
          SZooC=ZSpref(ng)*A(i,j,KrMax,nnew,iSzoo)
          LZooC=ZLpref(ng)*A(i,j,KrMax,nnew,iLzoo)
          PZooC=ZPpref(ng)*A(i,j,KrMax,nnew,iPzoo)
# endif
          TZooC=SZooC+LZooC+PZooC
          IF (TZooC.gt.ZooMax) THEN
            ZooMax=TZooC
            IrMax=i
            JrMax=j
          END IF
        END DO
      END DO
!
! If present cell is not best location, swim to best location
!      dxyG=(pm(Ir,Jr)**2+pn(Ir,Jr)**2)**.5
      diZ=Uswim*REAL(IrMax-Ir)*pm(Ir,Jr)
      djZ=Uswim*REAL(JrMax-Jr)*pn(Ir,Jr)
!
! If present cell is best location, swim against currents to stay
      IF ((diZ.eq.0.0_r8).and.(djZ.eq.0.0_r8)) THEN
        IF (track(ixrhs,itime,lflt).ge.0.0_r8) THEN
          diZ=-MIN(Uswim*pm(Ir,Jr),track(ixrhs,itime,lflt))
        ELSE
          diZ=MIN(Uswim*pm(Ir,Jr),-track(ixrhs,itime,lflt))
        END IF
        IF (track(iyrhs,itime,lflt).ge.0.0_r8) THEN
          djZ=-MIN(Uswim*pn(Ir,Jr),track(iyrhs,itime,lflt))
        ELSE
          djZ=MIN(Uswim*pn(Ir,Jr),-track(iyrhs,itime,lflt))
        END IF
      END IF
!
! k index of cell with max food to move fish in step_float
      kZmax=KrMax
!
#else
      SUBROUTINE fish_swim
#endif
      RETURN
      END SUBROUTINE fish_swim
