#include "cppdefs.h"
      MODULE fish_swim_mod
#if defined NONLINEAR && defined NEMURO_SAN
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

      implicit none

      PRIVATE
      PUBLIC  :: fish_swim

      CONTAINS
!
!***********************************************************************
      SUBROUTINE fish_swim (ng, LBi, UBi, LBj, UBj, LBk, UBk,           &
     &                          Lstr, Lend, itime, nnew,                &
     &                          pm, pn,                                 &
# ifdef SOLVE3D
     &                          Hz,                                     &
# endif
     &                          MyThread, bounded, track,               &
     &                          bioenergy, alive)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_biology
      USE mod_grid
      USE utility_mod, ONLY : nrng
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: Lstr, Lend, itime, nnew

      logical, intent(in) :: MyThread(Nfish(ng))
      logical, intent(in) :: bounded(Nfish(ng))
      logical, intent(in) :: alive(Nfish(ng))

      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
# endif
      real(r8), intent(in) :: bioenergy(NFishV(ng),Nfish(ng))
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfish(ng))
!
!  Local variable declarations.
!
      logical :: lflag=.TRUE.

      integer :: ierr
      integer :: iseed = 149876

      integer :: Ir, Jr, Kr, IrMax, JrMax, KrMax, kZmax
      integer :: i1, i2, j1, j2, i, j, k, itm1, l

      real(r8) :: SZooC, LZooC, PZooC, TZooC, ZooMax
      real(r8) :: Fweight, Fgrowth, Fmort, Ftemp, Fdens, Fpval
      real(r8) :: score, scoreMax, a_A, b_A, Uswim
      real(r8) :: d_i, d_j, d_ij, d_k, theta_ij, fac, snudg
      real(r8) :: vnudg(1)
      real(r8) :: Fvar, Fvar_opt, Fvar_sigma, var_f, var_g
      real(r8) :: dist, epsd, epsd_sigma, ee, epsx, epsy, h1, h2
      real(r8) :: f_xdis, g_epsx, f_ydis, g_epsy, xdis, ydis
      real(r8) :: f_zdis, g_epsz, epsz, zdis, zpos, cff1, cff2
      real(r8) :: tcheck
!
!-----------------------------------------------------------------------
!  Compute swimming velocity based on Railsback or Humston behavior
!-----------------------------------------------------------------------
!
      DO l=Lstr,Lend
        IF (MyThread(l).and.bounded(l).and.alive(l)) THEN
!
          Ir=INT(track(ixgrd,itime,l))
          Jr=INT(track(iygrd,itime,l))
          Kr=INT(track(izgrd,itime,l))
!
          Fweight=bioenergy(ifwwt,l)
!
          IF (Fweight.ge.0.0_r8) THEN
!
! Swimming velocity based on fish weight
! Uswim relationship is in cm/s
            a_A=15.0_r8
            b_A=0.13_r8
            Uswim=a_A*Fweight**b_A

# ifdef RAILSBACK
!
! MISSING HERE: CELL SEARCH RANDOMIZATION
!
            scoreMax=-99.0_r8
            IrMax=Ir
            JrMax=Jr
            KrMax=Kr
            DO j=Jr-1,Jr+1
              DO i=Ir-1,Ir+1
                DO k=MAX(Kr-1,1),MIN(Kr+1,N(ng))
                  CALL fish_growth_ijk (ng, i, j, k, nnew,                  &
     &                                      Fweight, Uswim, Fgrowth)
! MISSING HERE: CALL FOR MORTALITY SUBROUTINE
                  Fmort=0.0_r8
                  score=exp(Fgrowth)*exp(Fmort)
                  IF (score.gt.scoreMax) THEN
                    scoreMax=score
                    IrMax=i
                    JrMax=j
                    KrMax=k
                  END IF
                END DO
              END DO
            END DO
            CALL ran1_s (snudg)
            fac=2.0_r8*(snudg-0.5_r8)
            d_i=REAL(IrMax,r8)-track(ixgrd,itime,l)
            d_j=REAL(JrMax,r8)-track(iygrd,itime,l)
            d_k=REAL(KrMax,r8)-track(izgrd,itime,l)
! Angle with +-0.5 radians randomness
            theta_ij=atan2(d_j,d_i)+0.5_r8*fac
! Distance with 30% randomness
            d_ij=((d_i/pm(Ir,Jr))**2+(d_j/pn(Ir,Jr))**2)**0.5_r8
! Uswim is in cm/s, so multiply by 0.01 for m/s)
            d_ij=MIN(0.01_r8*Uswim*dt(ng),d_ij)
            d_ij=d_ij+0.3_r8*fac*d_ij
! Swim towards best location
            track(ixrhs,itime,l)=d_ij*cos(theta_ij)*pm(Ir,Jr)/dt(ng)
            track(iyrhs,itime,l)=d_ij*sin(theta_ij)*pn(Ir,Jr)/dt(ng)
            track(izrhs,itime,l)=d_k/dt(ng)
! Temporary overwrite until vertical motion is working
            track(izrhs,itime,l)=0.0_r8
!
# endif
!
# ifdef HUMSTON
!
            Ftemp=track(itemp+NFV(ng)-NT(ng),itime,l)
!            Fpval=bioenergy(ifpvl,l)
            Fdens=track(ifden,itime,l)
            tcheck=MOD(time(ng)/86400.0_r8,24.0_r8)
!
! Behavior based on temperature
            Fvar=Ftemp
            Fvar_opt=15.0_r8
            Fvar_sigma=3.0_r8
! Behavior based on p-value
!            Fvar=Fpval
!            Fvar_opt=1.0_r8
!            Fvar_sigma=0.2_r8
! Use commented lines below to check behavior based on longitude gradient
!            Fvar=track(ixgrd,itime,l)
!            Fvar_opt=15.0_r8
!            Fvar_sigma=5.0_r8
!
            itm1=MOD(itime-1+NFT+1,NFT+1)
            xdis=track(ixrhs,itm1,l)*dt(ng)/pm(Ir,Jr)
            ydis=track(iyrhs,itm1,l)*dt(ng)/pn(Ir,Jr)
! Use velx, vely below for purely random motion
!              xdis=0.0_r8
!              ydis=0.0_r8
!
! Horizontal behavior (updated once daily or when hits land)
            IF ((tcheck.eq.REAL(INT(tcheck)))) THEN
! Uswim is in cm/s, so multiply by 0.01 for m/s)
              dist=0.01_r8*Uswim*dt(ng)
! Compute random displacements
              epsd=(0.5_r8*dist**2)**0.5_r8
              epsd_sigma=0.5_r8*dist
              lflag=.TRUE.
              DO WHILE (lflag)
                CALL nrng (iseed, vnudg(:), 1, ierr)
                IF (ABS(vnudg(1)).le.3.0_r8) lflag=.FALSE.
              END DO
              ee=vnudg(1)*epsd_sigma+epsd
              CALL ran1_s (snudg)
              IF (snudg.lt.0.5_r8) epsx=-ee
              IF (snudg.ge.0.5_r8) epsx=ee
              lflag=.TRUE.
              DO WHILE (lflag)
                CALL nrng (iseed, vnudg(:), 1, ierr)
                IF (ABS(vnudg(1)).le.3.0_r8) lflag=.FALSE.
              END DO
              ee=vnudg(1)*epsd_sigma+epsd
              CALL ran1_s (snudg)
              IF (snudg.lt.0.5_r8) epsy=-ee
              IF (snudg.ge.0.5_r8) epsy=ee
! Compute happiness-based  displacements
              fac=EXP(-0.5_r8*((Fvar-Fvar_opt)/Fvar_sigma)**2)
              h1=0.75_r8
              h2=0.9_r8
              var_f=h1*fac
              var_g=1.0_r8-h2*fac
              f_xdis=xdis*var_f
              g_epsx=epsx*var_g
              f_ydis=ydis*var_f
              g_epsy=epsy*var_g
              xdis=f_xdis+g_epsx
              ydis=f_ydis+g_epsy
            END IF
! Reflect if hits land using density anomaly
            IF (Fdens.eq.0.0_r8) THEN
              xdis=-1.0_r8*xdis
              ydis=-1.0_r8*ydis
            END IF
! Swim towards best location
            track(ixrhs,itime,l)=xdis*pm(Ir,Jr)/dt(ng)
            track(iyrhs,itime,l)=ydis*pn(Ir,Jr)/dt(ng)
!
! Vertical behavior (update at every time step)
! DOES NOT WORK YET
            zdis=track(izgrd,itm1,l)*dt(ng)
! Uswim is in cm/s, so multiply by 0.01 for m/s)
!            dist=0.01_r8*Uswim*dt(ng)
! Max distance is 100m, or total water depth (whichever is less)
!            dist=min(1.0_r8,-GRID(ng)%z_w(Ir,Jr,1))
            dist=1.0_r8
! Compute random displacements
            epsd=(0.5_r8*dist**2)**0.5_r8
            epsd_sigma=0.5_r8*dist
            lflag=.TRUE.
            DO WHILE (lflag)
              CALL nrng (iseed, vnudg(:), 1, ierr)
              IF (ABS(vnudg(1)).le.3.0_r8) lflag=.FALSE.
            END DO
            ee=vnudg(1)*epsd_sigma+epsd
            CALL ran1_s (snudg)
            IF (snudg.lt.0.5_r8) epsz=-ee
            IF (snudg.ge.0.5_r8) epsz=ee
! Compute happiness-based  displacements
            fac=EXP(-0.5_r8*((Fvar-Fvar_opt)/Fvar_sigma)**2)
            h1=0.75_r8
            h2=0.9_r8
            var_f=h1*fac
            var_g=1.0_r8-h2*fac
            f_zdis=zdis*var_f
            g_epsz=epsz*var_g
            zdis=f_zdis+g_epsz
! Swim towards best location
!            track(izrhs,itime,l)=zdis/dt(ng)
! Temporary overwrite until vertical motion is working
            track(izrhs,itime,l)=0.0_r8
!
# endif
!
          END IF
!
        END IF
      END DO
!
      RETURN
      END SUBROUTINE fish_swim
#endif
      END MODULE fish_swim_mod
