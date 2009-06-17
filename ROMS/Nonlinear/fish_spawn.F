#include "cppdefs.h"
      MODULE fish_spawn_mod
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
      PUBLIC  :: fish_spawn

      CONTAINS
!
!***********************************************************************
      SUBROUTINE fish_spawn (ng, Lstr, Lend, itime,                     &
     &                           MyThread, bounded, track,              &
     &                           bioenergy, alive)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_biology
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend, itime

      logical, intent(in) :: MyThread(Nfish(ng))
      logical, intent(in) :: bounded(Nfish(ng))
      logical, intent(in) :: alive(Nfish(ng))

      real(r8), intent(in) :: track(NFV(ng),0:NFT,Nfish(ng))
      real(r8), intent(inout) :: bioenergy(NFishV(ng),Nfish(ng))
!
!  Local variable declarations.
!
      integer :: l
!
      real(r8) :: Ftemp, temp_opt, temp_sigma
      real(r8) :: sp_start, sp_end, sp_time, tcheck
!
!
!-----------------------------------------------------------------------
!  Compute feeding terms on zooplankton to update growth
!-----------------------------------------------------------------------
!
      temp_opt=15.0_r8
      temp_sigma=3.0_r8
!
      sp_start=120.0_r8
      sp_end=210.0_r8
! sp_time: time in days, modulo 360 days
      sp_time=REAL(INT(time(ng)/86400.0_r8/360.0_r8))
      sp_time=time(ng)/86400.0_r8-360.0_r8*sp_time
! tcheck: spawn at midnight
      tcheck=MOD(time(ng)/86400.0_r8,24.0_r8)
!
      IF ((sp_time.ge.sp_start).and.(sp_time.le.sp_end).and.            &
     &    (tcheck.eq.REAL(INT(tcheck)))) THEN
!
        DO l=Lstr,Lend
          IF (MyThread(l).and.bounded(l)) THEN
!
! Worth should be replaced by alive in check
!
            Ftemp=track(itemp+NFV(ng)-NT(ng),itime,l)
!
            IF ((alive(l)).and.                                         &
     &          (bioenergy(ifage,l).ge.Fage0(ng)).and.                  &
     &          (INT(bioenergy(ifbatch,l)).lt.Nbatch(ng)).and.          &
     &          (ABS(Ftemp-temp_opt).le.temp_sigma).and.                &
     &          ((sp_time-bioenergy(iftspwn,l)).ge.7.0_r8)) THEN
!
              bioenergy(ifbatch,l)=bioenergy(ifbatch,l)+1.0_r8
              bioenergy(iftspwn,l)=sp_time
!
            END IF
!
          END IF
        END DO
!
      END IF

      RETURN
      END SUBROUTINE fish_spawn
#endif
      END MODULE fish_spawn_mod
