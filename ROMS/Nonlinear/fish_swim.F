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
     &                          lflt, itime, nnew,                      &
     &                          diS, djS, kZmax, Fswim,                 &
     &                          pm, pn,                                 &
# ifdef SOLVE3D
     &                          Hz,                                     &
# endif
# ifdef MASKING
     &                          Amask,                                  &
# endif
     &                          A, track, bioenergy)
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
      logical, intent(inout) :: Fswim

      integer, intent(in) :: ng, LBi, UBi, LBj, UBj, LBk, UBk
      integer, intent(in) :: lflt, itime, nnew

      real(r8), intent(inout) :: diS, djS
      real(r8), intent(in) :: pm(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: pn(LBi:UBi,LBj:UBj)
# ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
# endif
# ifdef MASKING
      real(r8), intent(in) :: Amask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: A(LBi:UBi,LBj:UBj,LBk:UBk,3,1:NT(ng))
      real(r8), intent(in) :: bioenergy(NFishV(ng),Nfish(ng))

      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfish(ng))
!
!  Local variable declarations.
!
      logical :: lflag=.TRUE.

      integer :: ierr
      integer :: iseed = 149876

      integer :: Ir, Jr, Kr, IrMax, JrMax, KrMax, kZmax
      integer :: i1, i2, j1, j2, i, j, k

      real(r8) :: SZooC, LZooC, PZooC, TZooC, ZooMax
      real(r8) :: Fweight, Fgrowth, Fmort, Ftemp, Fdens, Fpval
      real(r8) :: score, scoreMax, a_A, b_A, Uswim
      real(r8) :: d_i, d_j, d_ij, theta_ij, fac, snudg
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
      Ir=INT(track(ixgrd,itime,lflt)-0.5_r8)
      Jr=INT(track(iygrd,itime,lflt)-0.5_r8)
      Kr=INT(track(izgrd,itime,lflt)-0.5_r8)
!
      Fweight=bioenergy(ifwwt,lflt)
!
      IF (Fweight.ge.0.0_r8) Fswim=.TRUE.
!
      IF (Fswim) THEN
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
        IrMax=0
        JrMax=0
        KrMax=0
        DO j=Jr,Jr+1
          DO i=Ir,Ir+1
            DO k=1,N(ng)
              CALL fish_growth_ijk (ng, Fweight, Uswim,                 &
#  ifdef MASKING
     &                              Amask(i,j),                         &
#  endif
     &                              A(i,j,k,nnew,:), Fgrowth)
!
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
!
        CALL ran1_s (snudg)
        fac=2.0_r8*(snudg-0.5_r8)
!
        d_i=REAL(IrMax,r8)-track(ixgrd,itime,lflt)
        d_j=REAL(JrMax,r8)-track(iygrd,itime,lflt)
! Angle with +-0.5 radians randomness
        theta_ij=atan2(d_j,d_i)+0.5_r8*fac
! Distance with 30% randomness
        d_ij=1.0_r8+0.3_r8*fac
!
! Swim towards best location
! Uswim is in cm/s, so multiply by 0.01 for m/s)
        diS=0.01_r8*Uswim*d_ij*cos(theta_ij)*pm(Ir,Jr)
        djS=0.01_r8*Uswim*d_ij*sin(theta_ij)*pn(Ir,Jr)
! k index of cell with max food to move fish in step_float
        kZmax=KrMax
!
# endif
!
# ifdef HUMSTON
!
        Ftemp=track(itemp+NFV(ng)-NT(ng),itime,lflt)
        Fpval=bioenergy(ifpvl,lflt)
        Fdens=track(ifden,itime,lflt)
        tcheck=MOD(time(ng)/86400.0_r8,24.0_r8)
!
! Behavior based on temperature
        Fvar=Ftemp
        Fvar_opt=15.0_r8
        Fvar_sigma=3.0_r8
! Behavior based on p-value
!        Fvar=Fpval
!        Fvar_opt=1.0_r8
!        Fvar_sigma=0.2_r8
! Use commented lines below to check behavior based on longitude gradient
!        Fvar=track(ixgrd,itime,lflt)
!        Fvar_opt=15.0_r8
!        Fvar_sigma=5.0_r8
!
! Horizontal behavior (updated once daily or when hits land)
        IF ((tcheck.eq.REAL(INT(tcheck))).or.(Fdens.eq.0.0_r8)) THEN
          xdis=bioenergy(ixgrd,lflt)
          ydis=bioenergy(iygrd,lflt)
! Use velx, vely below for purely random motion
!          xdis=0.0_r8
!          ydis=0.0_r8
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
! Reflect if hits land using density anomaly
          IF (Fdens.eq.0.0_r8) THEN
            xdis=-1.0_r8*xdis
            ydis=-1.0_r8*ydis
          END IF
! Swim towards best location
          diS=xdis*pm(Ir,Jr)/dt(ng)
          djS=ydis*pn(Ir,Jr)/dt(ng)
          bioenergy(ixgrd,lflt)=xdis
          bioenergy(iygrd,lflt)=ydis
        ELSE
! Swim towards best location
          xdis=bioenergy(ixgrd,lflt)
          ydis=bioenergy(iygrd,lflt)
          diS=xdis*pm(Ir,Jr)/dt(ng)
          djS=ydis*pn(Ir,Jr)/dt(ng)
        END IF
!
! Vertical behavior (update at every time step)
! DOES NOT WORK YET
        zdis=bioenergy(izgrd,lflt)
! Uswim is in cm/s, so multiply by 0.01 for m/s)
!        dist=0.01_r8*Uswim*dt(ng)
! Max distance is 100m, or total water depth (whichever is less)
!        dist=min(1.0_r8,-GRID(ng)%z_w(Ir,Jr,1))
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
        zpos=GRID(ng)%z_w(Ir,Jr,Kr)+zdis
        cff1=GRID(ng)%z_w(Ir,Jr,N(ng))
        kZmax=Kr
        DO k=N(ng)-1,1,-1
          cff2=GRID(ng)%z_w(Ir,Jr,k)
          IF ((zpos-cff2)*(cff1-zpos).ge.0.0_r8) THEN
            kZmax=k
          END IF
          cff1=cff2
        END DO
        bioenergy(izgrd,lflt)=zdis
!        write(6,*) 'HUMSTON   ',zdis,zpos
! Temporary overwrite until vertical motion is working
        kZmax=Kr
!
# endif
!
      END IF
!
      RETURN
      END SUBROUTINE fish_swim
#endif
      END MODULE fish_swim_mod
