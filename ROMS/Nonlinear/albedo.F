#include "cppdefs.h"
      MODULE albedo_mod
#if defined ALBEDO && !defined ALBEDO_FILE && !defined ANA_ALBEDO
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the albedo                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: albedo_eval

contains
!
!***********************************************************************
      SUBROUTINE albedo_eval (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
# ifdef ICE_MODEL
      USE mod_ice
# endif
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 17, __LINE__, __FILE__)
# endif
      CALL albedo_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
# ifdef ICE_MODEL
     &                  liold(ng),                                      &
     &                  linew(ng),                                      &
# endif
# if defined SHORTWAVE && defined ALBEDO_CURVE
     &                  GRID(ng) % latr,                                &
# endif
# ifdef ICE_MODEL
     &                  ICE(ng) % ai,                                   &
     &                  ICE(ng) % hi,                                   &
     &                  FORCES(ng) % albedo_ice,                        &
#  ifdef ICE_THERMO
     &                  ICE(ng) % hsn,                                  &
     &                  ICE(ng) % tis,                                  &
#  endif
#  ifdef ALBEDO_CSIM
     &                  FORCES(ng) % frswpen,                           &
#  endif
# endif
# ifdef ALBEDO_HACK
     &                  GRID(ng) % mask_albedo,                         &
# endif
     &                  FORCES(ng) % albedo                             &
     &                  )
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 17, __LINE__, __FILE__)
# endif
      RETURN
      END SUBROUTINE albedo_eval
!
!***********************************************************************
      SUBROUTINE albedo_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
# ifdef ICE_MODEL
     &                        liold, linew,                             &
# endif
# if defined SHORTWAVE && defined ALBEDO_CURVE
     &                        latr,                                     &
# endif
# ifdef ICE_MODEL
     &                        ai, hi, albedo_ice,                       &
#  ifdef ICE_THERMO
     &                        hsn, tis,                                 &
#  endif
#  ifdef ALBEDO_CSIM
     &                        frswpen,                                  &
#  endif
# endif
# ifdef ALBEDO_HACK
     &                        mask_albedo,                              &
# endif
     &                        albedo                                    &
     &                        )
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
# if defined BEST_NPZ
      USE mod_biology
# endif
!
# if defined BEST_NPZ && defined CLIM_ICE_1D
      USE mod_clima
# endif

      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
# ifdef ICE_MODEL
      integer, intent(in) :: liold
      integer, intent(in) :: linew
# endif
!
# ifdef ASSUMED_SHAPE
#  if defined SHORTWAVE && defined ALBEDO_CURVE
      real(r8), intent(in) :: latr(LBi:,LBj:)
#  endif
#  ifdef ICE_MODEL
      real(r8), intent(in) :: ai(LBi:,LBj:,:)
      real(r8), intent(in) :: hi(LBi:,LBj:,:)
      real(r8), intent(out) :: albedo_ice(LBi:,LBj:)
#   ifdef ICE_THERMO
      real(r8), intent(in) :: hsn(LBi:,LBj:,:)
      real(r8), intent(in) :: tis(LBi:,LBj:)
#   endif
#   ifdef ALBEDO_CSIM
      real(r8), intent(out) :: frswpen(LBi:,LBj:)
#   endif
#  endif
#  ifdef ALBEDO_HACK
      real(r8), intent(in) :: mask_albedo(LBi:,LBj:)
#  endif
      real(r8), intent(out) :: albedo(LBi:,LBj:)

# else
#  if defined SHORTWAVE && defined ALBEDO_CURVE
      real(r8), intent(in) :: latr(LBi:UBi,LBj:UBj)
#  endif
#  ifdef ICE_MODEL
      real(r8), intent(in) :: ai(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: hi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(out) :: albedo_ice(LBi:UBi,LBj:UBj)
#   ifdef ICE_THERMO
      real(r8), intent(in) :: hsn(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in) :: tis(LBi:UBi,LBj:UBj)
#   endif
#   ifdef ALBEDO_CSIM
      real(r8), intent(out) :: frswpen(LBi:UBi,LBj:UBj)
#   endif
#  endif
#ifdef ALBEDO_HACK
      real(r8), intent(in) :: mask_albedo(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(out) :: albedo(LBi:UBi,LBj:UBj)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, listp

      real(r8) :: cff1, cff2
      real(r8) :: cff
# ifdef ICE_BULK_FLUXES
#  ifdef ALBEDO_CSIM
      real(r8), parameter :: albicev=0.78_r8 !visible ice albedo for h > ahmax; CESM2/CICE5.1 default, DuVivier et al., 2018 (cf. 0.73 in Briegleb et al., 2004)
      real(r8), parameter :: albicei=0.36_r8 !near-ir ice albedo for h > ahmax; CESM2/CICE5.1 default, DuVivier et al., 2018 (cf. 0.33 in Briegleb et al., 2004)
      real(r8), parameter :: albsnowv=0.98_r8 !cold snow albedo, visible; CESM2/CICE5.1 default, DuVivier et al., 2018 (cf. 0.96 in Briegleb et al., 2004)
      real(r8), parameter :: albsnowi=0.70_r8 !cold snow albedo, near IR; CESM2/CICE5.1 default, DuVivier et al., 2018 (cf. 0.68 in Briegleb et al., 2004)
      real(r8), parameter :: ahmax=0.3_r8 !thickness above which ice albedo is constant (m); CESM2/CICE5.1 default, DuVivier et al., 2018 (cf. 0.5 m in Briegleb et al., 2004)
      real(r8), parameter :: snowpatch=0.02_r8 !parameter for fractional snow area (m); CESM2/CICE5.1 default, DuVivier et al., 2018 (cf. 0.02 m in Briegleb et al., 2004)
      real(r8), parameter :: dT_melt=1.0_r8 !change in temp to give dalb_mlt albedo change (C); CESM2/CICE5.1 default, DuVivier et al., 2018 (cf. 1.0 C in Briegleb et al., 2004)
      real(r8), parameter :: dalb_mlt=-0.075_r8 !albedo change per dT_melt change in temp for ice; CESM2/CICE5.1 default, DuVivier et al., 2018 (cf. -0.075 in Briegleb et al., 2004)
      real(r8), parameter :: dalb_mltv=-0.1_r8 !albedo vis change per dT_melt change in temp for snow; CESM2/CICE5.1 default, DuVivier et al., 2018 (cf. -0.1 in Briegleb et al., 2004)
      real(r8), parameter :: dalb_mlti=-0.15_r8 !albedo ir change per dT_melt change in temp for snow; CESM2/CICE5.1 default, DuVivier et al., 2018 (cf. -0.15 in Briegleb et al., 2004)
      real(r8), parameter :: Timelt=0.0_r8 !melting temperature, ice top surface (C)
      real(r8), parameter :: Tsmelt=0.0_r8 !melting temperature, snow top surface (C)
      real(r8), parameter :: puny=1.0e-11_r8
      real(r8), parameter :: fv=0.43_r8 !visual fraction of total shortwave (Aumont et al., 2015)
      real(r8), parameter :: i0vis=0.70_r8 !fraction of penetrating solar rad (visible); CESM2/CICE5.1 default, DuVivier et al., 2018 (cf. 0.7 in Briegleb et al., 2004)

      real(r8) :: fhtan, fh, albo, alvni, alini, dTs, fT
      real(r8) :: alvns, alins, asnow, alvn, alin
! Was these values
!     real(r8), parameter :: alb_i_thick=0.54_r8
!     real(r8), parameter :: alb_s_dry=0.83_r8
!     real(r8), parameter :: alb_s_wet=0.70_r8
! Values from Ungermann et al, 2017
!     real(r8), parameter :: alb_i_thick=0.71_r8
!     real(r8), parameter :: alb_s_dry=0.86_r8
!     real(r8), parameter :: alb_s_wet=0.79_r8
! Try something intermediate
!     real(r8), parameter :: alb_i_thick=0.65_r8
!     real(r8), parameter :: alb_s_dry=0.86_r8
!     real(r8), parameter :: alb_s_wet=0.72_r8
#  else
#ifdef ICE_BOX
      real(r8), parameter :: alb_i_dry=0.75_r8
      real(r8), parameter :: alb_i_wet=0.64_r8
      real(r8), parameter :: alb_s_dry=0.85_r8
      real(r8), parameter :: alb_s_wet=0.82_r8
#else
      real(r8), parameter :: alb_i_dry=0.65_r8
      real(r8), parameter :: alb_i_wet=0.60_r8
      real(r8), parameter :: alb_s_dry=0.85_r8
      real(r8), parameter :: alb_s_wet=0.72_r8
#endif
#  endif
      real(r8) :: albs, qlwi, qlh_i, qsh_i
      real(r8) :: le_i, dq_i,fqlat1, slp, Qsati
      real(r8) :: vap_p_i
#  ifdef ICE_ALB_EC92
      real(r8) :: albi, albsn, thki_n, thksn_n
#  endif
# endif
      real(r8), parameter :: alb_w=0.06_r8
# if defined ICE_MODEL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ice_thick
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: snow_thick
# endif

# include "set_bounds.h"
!-------------------------------------------------------------------------------
! PURPOSE:
!   computes albedo over snow/ice/water
!-------------------------------------------------------------------------------

# if defined ICE_MODEL
      IF (PerfectRST(ng) .and. iic(ng).eq.ntstart(ng)) THEN
        listp = liold
      ELSE
        listp = linew
      END IF
#  ifdef ALBEDO_CSIM
      fhtan = atan(ahmax*4.0_r8)
#  endif
# endif

! Note that this loop needs to be cleaned of all global arrays for
! OpenMP.
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
# if defined ICE_MODEL
! Calculate the ice/snow albedo
          ice_thick(i,j) = hi(i,j,listp)/(ai(i,j,listp)+0.001)
          snow_thick(i,j) = hsn(i,j,listp)/(ai(i,j,listp)+0.001)
!
#  ifdef ICE_ALB_EC92
! Ice and snow albedo is calculated from Ebert and Curry,1992b
          albi=0.06_r8
          albsn=0.06_r8
          IF (ai(i,j,listp) .ge. min_a(ng)) THEN
            thki_n = ice_thick(i,j)
            thki_n = MAX(thki_n,0.00001_r8)
            thksn_n = snow_thick(i,j)
            albi=0.082409_r8*LOG(thki_n)+0.485472_r8
            IF (thki_n.GE.1._r8) albi=0.07616_r8*thki_n+0.414492_r8
            IF (thki_n.GE.2._r8) albi=0.561632_r8
!approximated values for albsn(depends on COSZ, but small variation)
            albsn=0.83_r8
            albedo_ice(i,j)=albi
            IF (hsn(i,j,listp).GT.0._r8) albedo_ice(i,j)=albsn
!            IF (sfwat(i,j,listp).GT.0._r8) albs=0.10737_r8              &
!     &         +0.518_r8*EXP(-8.1_r8 *sfwat(i,j,listp)-0.47_r8)         &
!     &         +0.341_r8*EXP(-31.8_r8*sfwat(i,j,listp)-0.94_r8)         &
!     &         +0.131_r8*EXP(-2.6_r8 *sfwat(i,j,listp)-3.82_r8)
!            ENDIF
          ELSE
            albedo_ice(i,j)=alb_w
          ENDIF
#  elif defined ALBEDO_CSIM
! Here we follow ice_shortwave.F90/compute_albedos subroutine (CESM2/CICE5.1)
          IF (ai(i,j,listp) .gt. min_a(ng)) THEN
            ! bare ice, thickness dependence
            fh = min(atan(ice_thick(i,j)*4.0_r8)/fhtan, 1.0_r8)
            albo = alb_w*(1.0_r8-fh)
            alvni = albicev*fh + albo
            alini = albicei*fh + albo
            ! bare ice, temperature dependence
            dTs = Timelt - tis(i,j)
            fT = min(dTs/dT_melt-1.0_r8, 0.0_r8) !=0 for dTs>=dT_melt, =-1 for dTs=0
            alvni = max(alb_w, alvni - dalb_mlt*fT) !Ice albedo decreases with warming, capped by ocean albedo
            alini = max(alb_w, alini - dalb_mlt*fT)

            IF (snow_thick(i,j) > puny) THEN
              ! snow on ice, temperature dependence
              alvns = max(alb_w, albsnowv - dalb_mltv*fT)
              alins = max(alb_w, albsnowi - dalb_mlti*fT)
              !Note: the CESM2/CICE5.1 code (ice_shortwave.F90/compute_albedos) does not impose lower cap here,
              !      but it seems safer to do so for the case of tis>>Timelt (=> fT<<0)
              asnow = snow_thick(i,j) / (snow_thick(i,j) + snowpatch)
              alvn = alvni*(1.0_r8-asnow) + alvns*asnow !Combined ice/snow visual albedo
              alin = alini*(1.0_r8-asnow) + alins*asnow !Combined ice/snow near IR albedo
            ELSE
              asnow = 0.0_r8
              alvn = alvni !Combined ice/snow visual albedo
              alin = alini !Combined ice/snow near IR albedo
            ENDIF

            albedo_ice(i,j) = fv*alvn + (1.0_r8-fv)*alin !Combined ice/snow albedo wrt total incident SW (absorbed SW = (1-albedo_ice)*SW)
            frswpen(i,j) = (1.0_r8-alvni)*(1.0_r8-asnow)*i0vis*fv !Penetrating fraction of total incident SW (penetrating SW qi_o_n = frswpen*SW)
            !Note: We assume that only the visual part incident on snow-uncovered ice can penetrate through the surface layer (if not reflected)
          ELSE
            albedo_ice(i,j) = alb_w
            frswpen(i,j) = (1.0_r8-alb_w)*i0vis*fv !This should be arbitrary
          ENDIF
#  elif defined ICE_BOX
          IF (ai(i,j,listp) .gt. min_a(ng)) THEN
            IF (hsn(i,j,listp).gt.0._r8) THEN
              IF (tis(i,j) .gt. -1.0_r8) THEN
                albedo_ice(i,j) = alb_s_wet
              ELSE
                albedo_ice(i,j) = alb_s_dry
              ENDIF
            ELSE
              IF (tis(i,j) .gt. -1.0_r8) THEN
                albedo_ice(i,j) = alb_i_wet
              ELSE
                albedo_ice(i,j) = alb_i_dry
              ENDIF
            ENDIF
          ELSE
            albedo_ice(i,j)=alb_w
          ENDIF
#  else
          cff1 = alb_s_wet - alb_s_dry
          cff2 = alb_i_wet - alb_i_dry
          IF (ai(i,j,listp) .gt. min_a(ng)) THEN
            IF (hsn(i,j,listp).gt.0._r8) THEN
              IF (tis(i,j) .gt. -1.0_r8) THEN
                albedo_ice(i,j) = cff1*(tis(i,j)+1.0_r8)+alb_s_dry
              ELSE
                albedo_ice(i,j) = alb_s_dry
              ENDIF
            ELSE
              IF (tis(i,j) .gt. -1.0_r8) THEN
                albedo_ice(i,j) = cff2*(tis(i,j)+1.0_r8)+alb_i_dry
              ELSE
                albedo_ice(i,j) = alb_i_dry
              ENDIF
            ENDIF
          ELSE
            albedo_ice(i,j)=alb_w
          ENDIF
#  endif
#  ifdef ALBEDO_HACK
          albedo_ice(i,j) = max(albedo_ice(i,j), 0.95*mask_albedo(i,j))
#  endif
# endif
! Compute ocean albedo
# ifdef ALBEDO_CURVE
#  ifdef BIO_1D
!using lat for M2 for whole domain
          albedo(i,j) = (0.069_r8 - 0.011_r8*                           &
     &                        cos(2*deg2rad*56.877))
#  else
          albedo(i,j) = (0.069_r8 - 0.011_r8*                           &
     &                        cos(2*deg2rad*latr(i,j)))
#  endif
# else
          albedo(i,j)=alb_w
# endif
# ifdef ALBEDO_HACK
          albedo(i,j) = max(albedo(i,j), 0.95*mask_albedo(i,j))
# endif
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          albedo)
# ifdef ICE_MODEL
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          albedo_ice)
#  ifdef ALBEDO_CSIM
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          frswpen)
#  endif
# endif
      END IF
# ifdef DISTRIBUTE
#  ifdef ICE_MODEL
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    albedo, albedo_ice)
#   ifdef ALBEDO_CSIM
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    frswpen)
#   endif
#  else
      CALL mp_exchange2d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    albedo)
#  endif
# endif

      RETURN
      END SUBROUTINE albedo_tile
#endif

      END module albedo_mod
