      SUBROUTINE biology (ng,tile)
!
!svn $Id$
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !

!TO DO : clean CPP keys OXYGEN/CARBON
!TO DO : clean/customize IO, we will probably need some other stuff

      USE mod_param
      USE mod_forces
      USE mod_ncparam
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
      USE mod_diags
      USE mod_mixing
      USE shapiro_mod
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15)
#endif
      CALL biology_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   GRID(ng) % h,                                  &
     &                   GRID(ng) % omn,                                &
     &                   FORCES(ng) % srflx,                            &
#ifdef OPTIC_MANIZZA
     &                   OCEAN(ng) % decayW,                            &
#endif
#ifdef BULK_FLUXES
     &                   OCEAN(ng) % u10_neutral,                       &
     &                   FORCES(ng) % Uwind,                            &
     &                   FORCES(ng) % Vwind,                            &
     &                   FORCES(ng) % Pair,                             &
#else
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
#endif
#ifdef COBALT_CARBON
     &                   FORCES(ng) % atmCO2,                           &
#endif
#ifdef COBALT_IRON
     &                   FORCES(ng) % soluble_fe,                       &
     &                   FORCES(ng) % mineral_fe,                       &
     &                   FORCES(ng) % ironsed,                          &
#endif
     &                   MIXING(ng) % hsbl,                             &
     &                   MIXING(ng) % ksbl,                             &
#ifdef DIAGNOSTICS_BIO
     &                   DIAGS(ng) % DiaBio2d,                          &
     &                   DIAGS(ng) % DiaBio3d,                          &
#endif
     &                   OCEAN(ng) % rho,                               &
     &                   OCEAN(ng) % zeta,                              &
     &                   OCEAN(ng) % obgc,                              &
#ifdef BENTHIC
     &                   OCEAN(ng) % bt,                                &
#endif
#ifdef TIMESERIES
     &                   OCEAN(ng) % tms,                               &
#endif
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif
      RETURN
      END SUBROUTINE biology

!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
#ifdef MASKING
     &                         rmask,                                   &
#endif
     &                         Hz, z_r, z_w, h, omn, srflx,             &
#ifdef OPTIC_MANIZZA
     &                         decayW,                                  &
#endif
#ifdef BULK_FLUXES
     &                         u10_neutral,                             &
     &                         Uwind, Vwind, Pair,                      &
#else
     &                         sustr, svstr,                            &
#endif
#ifdef COBALT_CARBON
     &                         atmCO2,                                  &
#endif
#ifdef COBALT_IRON
     &                         soluble_fe, mineral_fe, ironsed,         &
#endif
     &                         hsbl, ksbl,                              &
#ifdef DIAGNOSTICS_BIO
     &                         DiaBio2d, DiaBio3d,                      &
#endif
     &                         rho, zeta,                               &
     &                         obgc,                                    &
#ifdef BENTHIC
     &                         bt,                                      &
#endif
#ifdef TIMESERIES
     &                         tms,                                     &
#endif
     &                         t)
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mod_parallel
      USE shapiro_mod
      USE FMS_ocmip2_co2calc_mod, only : FMS_ocmip2_co2calc, CO2_dope_vector
      USE mod_iounits, only : stdout
# ifdef SOLVE3D
      USE bc_3d_mod
# endif

#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_reduce
      USE distribute_mod, ONLY : mp_reduce2
#endif
!  Imported variable declarations.
!

      implicit none

      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

! RD: cleaning of CPP  keys needed ASSUME_SHAPE, OXYGEN, CARBON !!!
#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: omn(LBi:,LBj:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
# ifdef OPTIC_MANIZZA
      real(r8), intent(in) :: decayW(LBi:,LBj:,0:,:)
#endif
# ifdef BULK_FLUXES
      real(r8), intent(in) :: u10_neutral(LBi:,LBj:)
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
      real(r8), intent(in) :: Pair(LBi:,LBj:)
# else
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
# endif
# ifdef COBALT_CARBON
      real(r8), intent(in) :: atmCO2(LBi:,LBj:)
#endif
# ifdef COBALT_IRON
      real(r8), intent(in) :: soluble_fe(LBi:,LBj:)
      real(r8), intent(in) :: mineral_fe(LBi:,LBj:)
      real(r8), intent(in) :: ironsed(LBi:,LBj:)
#endif
      real(r8), intent(in)    :: hsbl(LBi:,LBj:)
      integer,  intent(in)    :: ksbl(LBi:,LBj:)
#ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:,LBj:,:)
      real(r8), intent(inout) :: DiaBio3d(LBi:,LBj:,:,:)
#endif
      real(r8), intent(in)    :: rho(LBi:,LBj:,:)
      real(r8), intent(in)    :: zeta(LBi:,LBj:,:)
      real(r8), intent(inout) :: obgc(LBi:,LBj:,:,:,:)
#ifdef BENTHIC
      real(r8), intent(inout) :: bt(LBi:,LBj:,:,:,:)
#endif
#ifdef TIMESERIES
      real(r8), intent(inout) :: tms(:,:,:)
#endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: omn(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
# ifdef OPTIC_MANIZZA
      real(r8), intent(in) :: decayW(LBi:UBi,LBj:UBj,0:UBk,4)
#endif
# ifdef BULK_FLUXES
      real(r8), intent(in) :: u10_neutral(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Pair(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
# endif
# ifdef COBALT_CARBON
      real(r8), intent(in) :: atmCO2(LBi:UBi,LBj:UBj)
# endif
# ifdef COBALT_IRON
      real(r8), intent(in) :: soluble_fe(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: mineral_fe(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ironsed(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in)    :: hsbl(LBi:UBi,LBj:UBj)
      integer,  intent(in)    :: ksbl(LBi:UBi,LBj:UBj)
#ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:UBi,LBj:UBj,NDbio2d)
      real(r8), intent(inout) :: DiaBio3d(LBi:UBi,LBj:UBj,UBk,NDbio3d)
#endif
      real(r8), intent(in)    :: rho(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in)    :: zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: obgc(LBi:UBi,LBj:UBj,N(ng),3,NOBGC) ! RD: not sure it'd work
#ifdef BENTHIC
      real(r8), intent(inout) :: bt(LBi:UBi,LBj:UBj,NBL(ng),3,UBt) ! RD: not sure it'd work
#endif
#ifdef TIMESERIES
      real(r8), intent(inout) :: tms(1,1,NT) ! RD: not sure it'd work
#endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif

!--------------------------------------------------------------------------------
!  Local variable declarations.
!
    integer :: i, j, k, nphyto, nzoo, nprey ! loop indices
    integer :: ntau, kblt, m, k_100, k_200 
    integer :: tau

    !real, dimension(:,:,:) ,pointer :: grid_tmask
    !integer, dimension(:,:),pointer :: mask_coast,grid_kmt
    ! convert to ROMS arrays
    ! loop index for sinking
    integer :: isink, ibio

    !
    !------------------------------------------------------------------------
    ! Local Variables
    !------------------------------------------------------------------------
    !
!    real, dimension(:),   allocatable     :: tmp_irr_band
    real, dimension(:,:), allocatable     :: rho_dzt_100, rho_dzt_200
    real, dimension(1:NUM_ZOO,1:NUM_PREY) :: ipa_matrix,pa_matrix,ingest_matrix
    ! RD debugging
    real, dimension(1:NUM_ZOO,1:NUM_PREY) :: ingest_matrix_glo
    real, dimension(1:NUM_PREY)           :: hp_ipa_vec,hp_pa_vec,hp_ingest_vec
    real, dimension(1:NUM_PREY)           :: prey_vec,prey_p2n_vec,prey_fe2n_vec,prey_si2n_vec
    real, dimension(1:NUM_ZOO)            :: tot_prey
    real                                  :: r_dt
    real                                  :: feprime
    real                                  :: juptake_di_tot2nterm
    real                                  :: log_btm_flx
    real                                  :: P_C_m
    real                                  :: p_lim_nhet
    real                                  :: TK, PRESS, PKSPA, PKSPC
    real                                  :: tmp_hblt, tmp_irrad, tmp_irrad_ML,tmp_opacity
    real                                  :: drho_dzt
    real                                  :: tot_prey_hp, sw_fac_denom, ingest_p2n, refuge_conc 
    real                                  :: bact_ldon_lim, bact_uptake_ratio, vmax_bact
    real                                  :: fpoc_btm, log_fpoc_btm
    logical                               :: used
    integer                               :: nb

    real(8), dimension(IminS:ImaxS,JminS:JmaxS) :: twodim_int_diag

    real(8), dimension(IminS:ImaxS,JminS:JmaxS,1:N(ng)) :: rmask3d
    real(8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: FM,swdk3, chl_conc
    real(8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng),4) :: DK
    real(8), dimension(IminS:ImaxS,JminS:JmaxS,1:N(ng)) :: tmp_decay
    real(8), dimension(IminS:ImaxS,JminS:JmaxS) :: mxl_depth, mxl_subgrid
    real(8), dimension(IminS:ImaxS,JminS:JmaxS) :: rho_ref
    integer(4), dimension(IminS:ImaxS,JminS:JmaxS) :: mxl_blev

    integer(4) :: overflow

    real(8) :: rho_crit = 0.03d0
    real(8) :: rho_ref_depth = 10.0d0
    real(8) :: drho_crit 
    real(8) :: cff_ts, cff_btf

    type(CO2_dope_vector) :: CO2_dope_vec

    ! Sinking routine
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: ndet_sinking
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: sidet_sinking
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: cadet_calc_sinking
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: cadet_arag_sinking
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: pdet_sinking
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: fedet_sinking
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: lithdet_sinking
    real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC

    real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
    real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
    real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
    real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
    real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
    real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
    real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
    real(r8), dimension(IminS:ImaxS,N(ng)) :: qc

    real(r8) :: cff, cff1, cff2, cff3, cff4, cff5, cff6, cff7, cff8
    real(r8) :: fac, fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8
    real(r8) :: cffL, cffR, cu, dltL, dltR

    integer, dimension(IminS:ImaxS,N(ng)) :: ksource
    integer :: ks
    real(r8) :: dtdays

    real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: sed_losses

    ! diagnostics 3d reduced
    integer                         ::  ndiag_int3
    integer                         :: it
    real(r8), allocatable           :: total_tracer(:)
    character (len=3), allocatable  :: op_handle(:)
    real(r8)                        :: sum_t
    real(r8)                        :: zvolume

    ! gas transfer air/ocean
    real(r8)                                     :: sal, ST
    real(r8)                                     :: co2_trans_vel
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: airsea_co2_flx
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: delta_co2star
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: co2_sc_no, co2_alpha, co2_csurf
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: u10squ, co2_piston_vel, pCO2atm, pCO2surf

    real(r8)                                     :: o2_saturation, tt, ts, ts2, ts3, ts4, ts5
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: airsea_o2_flx, o2_sc_no, o2_piston_vel, o2_csurf

    ! dust deposition
    real(r8)                                         :: frac_iron=0.035d0
    real(r8)                                         :: extra_frac_iron=0.03d0
    real(r8)                                         :: iron_dust_solub=0.02d0
    real(r8)                                         :: molar_feO2=87.8438d0    ! g/mol
    real(r8)                                         :: molar_fe=55.8450d0      ! g/mol
    real(r8)                                         :: dust_pen_depth=600.0d0  ! m
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,UBk) :: iron_dust_src
    real(r8), dimension(IminS:ImaxS,JminS:JmaxS,UBk) :: lith_dust_src
    

!-----------------------------------------------------------------------

#include "set_bounds.h"
!
#ifdef DISTRIBUTE
!  DISTRIBUTE is TRUE and model is needed
#endif

#ifdef DEBUG_COBALT
IF ( Master ) WRITE(stdout,*) '>>>   --------------- Cobalt debugging prints ------------'
#endif


!----------------------------------------------------------------------------
! reset loop index to this value to force overflow if any i,j,k out of loops

  overflow=1e9 

!----------------------------------------------------------------------------

  refuge_conc = 1.0e-9

#ifndef COBALT_DO_NOTHING
!----------------------------------------------------------------------------
! for volume integral of tracers

  ndiag_int3 = NT(ng) ! for now, we may want to increase that number then some
                      ! modifications will be needed in mod_ocean
  IF (.NOT. ALLOCATED(total_tracer))    ALLOCATE (total_tracer(ndiag_int3))
  IF (.NOT. ALLOCATED(op_handle))       ALLOCATE (op_handle(ndiag_int3))

  DO it=1,ndiag_int3
     op_handle(it) = 'SUM'
  ENDDO

!----------------------------------------------------------------------------
! indices for dope vector (used in FMS)
! RD: double-check me

!  CO2_dope_vec%isd = LBi ; CO2_dope_vec%isc = LBi
!  CO2_dope_vec%ied = UBi ; CO2_dope_vec%iec = UBi
!  CO2_dope_vec%jsd = LBj ; CO2_dope_vec%jsc = LBj
!  CO2_dope_vec%jed = UBj ; CO2_dope_vec%jec = UBj
  CO2_dope_vec%isd = IminS ; CO2_dope_vec%isc = IminS
  CO2_dope_vec%ied = ImaxS ; CO2_dope_vec%iec = ImaxS
  CO2_dope_vec%jsd = JminS ; CO2_dope_vec%jsc = JminS
  CO2_dope_vec%jed = JmaxS ; CO2_dope_vec%jec = JmaxS

!----------------------------------------------------------------------------
!
!  Create a 3D-mask to avoid loops in tracer update 
!
!----------------------------------------------------------------------------

  rmask3d = 0.

  DO k=1,UBk
   DO j=Jstr,Jend
    DO i=Istr,Iend
     rmask3d(i,j,k) = rmask(i,j)
    ENDDO
   ENDDO
  ENDDO
  k=overflow

!!!! RD DO NOTHING
#endif

!----------------------------------------------------------------------------
!
!  Time-stepping 
!
!----------------------------------------------------------------------------

  r_dt = 1.0d0 / dt(ng)
  tau = 1

!----------------------------------------------------------------------------
!
! At first step, allocate the arrays and copy the parameters to cobalt types
!
!----------------------------------------------------------------------------

      IF (iic(ng).eq.ntstart(ng)) THEN

        CALL cobalt_alloc_arrays(ng, tile,                          &
     &                               IminS, ImaxS, JminS, JmaxS, UBk )
        CALL cobalt_add_params(ng)

      ENDIF

#ifndef COBALT_DO_NOTHING

    ! RD dev notes:
    ! read temperature and salinity to cobalt array and units
    DO k=1,UBk
      DO j=Jstr,Jend
        DO i=Istr,Iend
           cobalt%f_temp(i,j,k) = max( t(i,j,k,nstp,itemp) ,  0.0d0 )
           cobalt%f_salt(i,j,k) = max( t(i,j,k,nstp,isalt) ,  0.0d0 )
        ENDDO
      ENDDO
    ENDDO
    i=overflow ; j=overflow ; k=overflow

#ifdef DEBUG_COBALT
IF ( Master ) WRITE(stdout,*) '>>>    Before CALL FMS surface min/max(temp)  =', MINVAL(cobalt%f_temp(:,:,:)),  MAXVAL(cobalt%f_temp(:,:,:))
IF ( Master ) WRITE(stdout,*) '>>>    Before CALL FMS surface min/max(salt) =', MINVAL(cobalt%f_salt(:,:,:)), MAXVAL(cobalt%f_salt(:,:,:))
#endif

#ifdef DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
! If appropriate, initialize time-averaged diagnostic arrays.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsDIA(ng)).and.                                 &
     &     (MOD(iic(ng),nDIA(ng)).eq.1)).or.                            &
     &    ((iic(ng).ge.ntsDIA(ng)).and.(nDIA(ng).eq.1)).or.             &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
        DO it=1,NDbio2d
          DO j=Jstr,Jend
            DO i=Istr,Iend
              DiaBio2d(i,j,it)=0.0_r8
            END DO
          END DO
        END DO
        DO it=1,NDbio3d
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                DiaBio3d(i,j,k,it)=0.0_r8
              END DO
            END DO
          END DO
        END DO
      END IF
#endif

#ifdef COBALT_CARBON
!---------------------------------------------------------------------
!Calculate co3_ion
!Also calculate co2 fluxes csurf and alpha for the next round of exchange
!---------------------------------------------------------------------
!   
! RD dev notes :
! htotal and co3_ion are read from restart even though it is not prognostic so that we
! don't have a glitch at the begining of each job
       DO k=1,UBk
         DO j=Jstr,Jend
           DO i=Istr,Iend
              cobalt%f_htotal(i,j,k)  = max( obgc(i,j,k,nstp,iohtotal) , 0.0d0 )
              cobalt%f_co3_ion(i,j,k) = max( obgc(i,j,k,nstp,ioco3_ion) , 0.0d0 )
           ENDDO
         ENDDO
       ENDDO
       i=overflow ; j=overflow ; k=overflow
!
! RD dev notes :
! Copy some tracer variables (po4, sio4, dic and alk) to cobalt arrays
! these will be used for the carbon chemistry routines
! variables will be copied again after carbon routines
    DO k=1,UBk
      DO j=Jstr,Jend
        DO i=Istr,Iend
           cobalt%f_po4(i,j,k)    = max( t(i,j,k,nstp,ipo4) ,  0.0d0 )
           cobalt%f_sio4(i,j,k)   = max( t(i,j,k,nstp,isio4),  0.0d0 )
           cobalt%f_dic(i,j,k)    = max( t(i,j,k,nstp,idic) ,  0.0d0 )
           cobalt%f_alk(i,j,k)    = max( t(i,j,k,nstp,ialk) ,  0.0d0 )
           cobalt%f_o2(i,j,k)     = max( t(i,j,k,nstp,io2)  ,  0.0d0 )
        ENDDO
      ENDDO
    ENDDO
    i=overflow ; j=overflow ; k=overflow
!
! RD dev notes :
! carbon chemistry routines are first applied to surface ( k = UBk )
    k=UBk
    cobalt%htotallo(:,:) = cobalt%htotal_scale_lo * cobalt%f_htotal(:,:,k)
    cobalt%htotalhi(:,:) = cobalt%htotal_scale_hi * cobalt%f_htotal(:,:,k)

#ifdef DEBUG_COBALT
IF ( Master ) WRITE(stdout,*) '>>>    Before CALL FMS surface min/max(htotal)  =', MINVAL(cobalt%f_htotal(:,:,k)),  MAXVAL(cobalt%f_htotal(:,:,k))
IF ( Master ) WRITE(stdout,*) '>>>    Before CALL FMS surface min/max(co3_ion) =', MINVAL(cobalt%f_co3_ion(:,:,k)), MAXVAL(cobalt%f_co3_ion(:,:,k))
#endif


  CALL FMS_ocmip2_co2calc(CO2_dope_vec,rmask3d(:,:,k), &
 &       cobalt%f_temp(:,:,k),                         &
 &       cobalt%f_salt(:,:,k),                         &
 &       cobalt%f_dic(:,:,k),                          &
 &       cobalt%f_po4(:,:,k),                          &
 &       cobalt%f_sio4(:,:,k),                         &
 &       cobalt%f_alk(:,:,k),                          &
 &       cobalt%htotallo, cobalt%htotalhi,             &
                                !InOut
 &       cobalt%f_htotal(:,:,k),                       &
                                !OUT
 &       co2star=cobalt%co2_csurf(:,:),                &
 &       alpha=cobalt%co2_alpha(:,:),                  &
 &       pCO2surf=cobalt%pco2_csurf(:,:),              &
 &       co3_ion=cobalt%f_co3_ion(:,:,k))
!
#ifdef DEBUG_COBALT
IF ( Master ) WRITE(stdout,*) '>>>    After CALL FMS surface min/max(htotal)  =', MINVAL(cobalt%f_htotal(:,:,k)),  MAXVAL(cobalt%f_htotal(:,:,k))
IF ( Master ) WRITE(stdout,*) '>>>    After CALL FMS surface min/max(co3_ion) =', MINVAL(cobalt%f_co3_ion(:,:,k)), MAXVAL(cobalt%f_co3_ion(:,:,k))
#endif
!
#ifdef DIAGNOSTICS_BIO
! RD dev notes:
! save diagnostics on carbon chemistry
   DO j=Jstr,Jend
     DO i=Istr,Iend
        DiaBio2d(i,j,ialpha)    = DiaBio2d(i,j,ialpha)    + cobalt%co2_alpha(i,j)
        DiaBio2d(i,j,ico2star)  = DiaBio2d(i,j,ico2star)  + cobalt%co2_csurf(i,j)
        DiaBio2d(i,j,ipco2surf) = DiaBio2d(i,j,ipco2surf) + cobalt%pco2_csurf(i,j)
     ENDDO
   ENDDO
#endif
!
! RD dev notes : gas transfer for CO2
! code adapted from cobalt set_boundary_values
!
   DO j=Jstr,Jend
     DO i=Istr,Iend

       !This calculation needs an input of SST and SSS
       sal = cobalt%f_salt(i,j,UBk)
       ST  = cobalt%f_temp(i,j,UBk)

       !---------------------------------------------------------------------
       !     CO2
       !---------------------------------------------------------------------

       !---------------------------------------------------------------------
       !  Compute the Schmidt number of CO2 in seawater using the
       !  formulation presented by Wanninkhof (1992, J. Geophys. Res., 97,
       !  7373-7382).
       !---------------------------------------------------------------------

       ! RD dev notes : replaced namelist values by values in TABLE A1 from Wanninkhof (1992)
       co2_sc_no(i,j) = ( 2073.1d0 - 125.62d0 * ST + 3.6276d0 * ST**2 - 0.043219d0 * ST**3 ) * rmask(i,j)

       ! RD dev notes : use 10 meters neutral wind and not wind from forcing
       ! files, this field is computed by ccsm_flux.F (Large/Yeager) or bulk_flux.F (Fairall)
       u10squ(i,j) = u10_neutral(i,j)*u10_neutral(i,j) 

       ! RD dev notes : piston/transfer velocity also from Wanninkhof (1992) eq. 3 (valid for short-term
       ! or steady winds), this piston velocity is given in cm/hour
       co2_piston_vel(i,j) = 0.31d0 * u10squ(i,j) * sqrt(660.0 / (co2_sc_no(i,j) + epsln))

       ! convert piston velocity in m.s-1
       co2_piston_vel(i,j) = co2_piston_vel(i,j) * 0.01d0 / 3600.0d0

       ! CO2 solubility computed in FMS module is in mol/kg/atm, convert it in
       ! mol/m3/atm
       co2_alpha(i,j) = cobalt%co2_alpha(i,j) * cobalt%Rho_0 

       ! convert ppmv to partial pressure in atm
       pCO2atm(i,j)  = atmCO2(i,j) * 1.0e-6             
       pCO2surf(i,j) = cobalt%pco2_csurf(i,j) * 1.0e-6 

       ! Wanninkhov 1992 equation A2 
       ! here the flux is taken positive when CO2 goes from atmosphere -> ocean
       ! co2_piston_vel = k [m.s-1]
       ! co2_alpha = L [mol.m-3.atm-1]
       ! pCO2atm and pCO2surf [atm]
       airsea_co2_flx(i,j) = co2_piston_vel(i,j) * co2_alpha(i,j) * (pCO2atm(i,j) - pCO2surf(i,j))

       ! convert airsea_co2_flx from mol.m-2.s-1 to mol.kg-1.s-1 
       airsea_co2_flx(i,j) = airsea_co2_flx(i,j) / cobalt%Rho_0 / Hz(i,j,UBk)

#ifdef COBALT_CONSERVATION_TEST
       airsea_co2_flx(i,j) = 0.0d0
#endif

#ifdef DIAGNOSTICS_BIO
       DiaBio2d(i,j,ico2_flx) = DiaBio2d(i,j,ico2_flx) + airsea_co2_flx(i,j) * dt(ng)
#endif

       !---------------------------------------------------------------------
       !     O2
       !---------------------------------------------------------------------
       !  Compute the oxygen saturation concentration at 1 atm total
       !  pressure in mol/kg given the temperature (t, in deg C) and
       !  the salinity (s, in permil)
       !
       !  From Garcia and Gosrdon (1992), Limnology and Oceonography.
       !  The formula used is from page 1310, eq (8).
       !
       !  *** Note: the "a3*ts^2" term (in the paper) is incorrect. ***
       !  *** It shouldn't be there.                                ***
       !
       !  o2_saturation is defined between T(freezing) <= T <= 40 deg C and
       !                                   0 permil <= S <= 42 permil
       !
       ! check value: T = 10 deg C, S = 35 permil,
       !              o2_saturation = 0.282015 mol m-3
       !---------------------------------------------------------------------
       !

       tt = 298.15 - ST
       tk = 273.15 + ST
       ts = log(tt / tk)
       ts2 = ts  * ts
       ts3 = ts2 * ts
       ts4 = ts3 * ts
       ts5 = ts4 * ts

       o2_saturation = (1000.0/22391.6) * rmask(i,j) *                   & !convert from ml/l to mol m-3
 &           exp( cobalt%a_0 + cobalt%a_1*ts + cobalt%a_2*ts2 +          &
 &           cobalt%a_3*ts3 + cobalt%a_4*ts4 + cobalt%a_5*ts5 +          &
 &           (cobalt%b_0 + cobalt%b_1*ts + cobalt%b_2*ts2 +              &
 &           cobalt%b_3*ts3 + cobalt%c_0*sal)*sal)

       !---------------------------------------------------------------------
       !  Compute the Schmidt number of O2 in seawater using the
       !  formulation proposed by Keeling et al. (1998, Global Biogeochem.
       !  Cycles, 12, 141-163).
       !---------------------------------------------------------------------
       !

       ! RD dev notes : took values from Sarmiento/Gruber table 3.3.1 (verif reference above)
       o2_sc_no(i,j)  = (1638.0d0 - 81.83d0 * ST + 1.483d0 * ST**2 - 0.008004 * ST**3) * rmask(i,j)

       ! RD dev notes : piston/transfer velocity also from Wanninkhof (1992) eq. 3 (valid for short-term
       ! or steady winds), this piston velocity is given in cm/hour
       o2_piston_vel(i,j) = 0.31d0 * u10squ(i,j) * sqrt(660.0 / (o2_sc_no(i,j) + epsln))

       ! convert piston velocity in m.s-1
       o2_piston_vel(i,j) = o2_piston_vel(i,j) * 0.01d0 / 3600.0d0

       ! convert o2 from mol/kg to mol/m3
       o2_csurf(i,j) = cobalt%f_o2(i,j,UBk) * cobalt%Rho_0 

       ! flux is k.(O2sat - [O2]) and is in mol.m-2.s-1
       airsea_o2_flx(i,j) = o2_piston_vel(i,j) * (o2_saturation - o2_csurf(i,j))

       ! convert flux in mol.kg-1.s-1
       airsea_o2_flx(i,j) = airsea_o2_flx(i,j) / cobalt%Rho_0 / Hz(i,j,UBk)

       ! don't think I need this now
       !o2_alpha(i,j) = (o2_saturation / 0.21)

#ifdef COBALT_CONSERVATION_TEST
        airsea_o2_flx(i,j) = 0.0d0
#endif
#ifdef DIAGNOSTICS_BIO
        DiaBio2d(i,j,io2_flx) = DiaBio2d(i,j,io2_flx) + airsea_o2_flx(i,j) * dt(ng)
#endif
     ENDDO
   ENDDO

!   IF( Master ) WRITE(stdout,*) 'Atmospheric CO2 is ', MINVAL(atmCO2) , MAXVAL(atmCO2)

   ! RD dev notes :
   ! carbon chemistry routines are then applied from susurface to bottom layer
   DO k=UBk-1,1,-1

      cobalt%htotallo(:,:) = cobalt%htotal_scale_lo * cobalt%f_htotal(:,:,k)
      cobalt%htotalhi(:,:) = cobalt%htotal_scale_hi * cobalt%f_htotal(:,:,k)

  ! RD needs cleaning : I don't like having LBi:UBi,LBj:UBj in t array

      CALL FMS_ocmip2_co2calc(CO2_dope_vec,rmask3d(:,:,k), &
     &       cobalt%f_temp(:,:,k),                         &
     &       cobalt%f_salt(:,:,k),                         &
     &       cobalt%f_dic(:,:,k),                          &
     &       cobalt%f_po4(:,:,k),                          &
     &       cobalt%f_sio4(:,:,k),                         &
     &       cobalt%f_alk(:,:,k),                          &
     &       cobalt%htotallo, cobalt%htotalhi,             &
                                    !InOut
     &       cobalt%f_htotal(:,:,k),                       &
                                    !OUT
     &       co3_ion=cobalt%f_co3_ion(:,:,k))

   ENDDO
   k=overflow

#endif

!!!! DO_NOTHING  key
#endif 
     
!----------------------------------------------------------------------------
! Copy the prognostics variables from ROMS tracers array to a first set of
! cobalt variables used to compute (something else)
!----------------------------------------------------------------------------

!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend
      ! Nitrogen dynamics
      phyto(SMALL)%f_n(i,j,k)    = max( t(i,j,k,nstp,insm) ,   0.0d0 )
      phyto(LARGE)%f_n(i,j,k)    = max( t(i,j,k,nstp,inlg) ,   0.0d0 )
      phyto(DIAZO)%f_n(i,j,k)    = max( t(i,j,k,nstp,indi) ,   0.0d0 )
      zoo(1)%f_n(i,j,k)          = max( t(i,j,k,nstp,insmz) ,  0.0d0 )
      zoo(2)%f_n(i,j,k)          = max( t(i,j,k,nstp,inmdz) ,  0.0d0 )
      zoo(3)%f_n(i,j,k)          = max( t(i,j,k,nstp,inlgz) ,  0.0d0 )
      cobalt%f_ldon(i,j,k)       = max( t(i,j,k,nstp,ildon) ,  0.0d0 )
      cobalt%f_sldon(i,j,k)      = max( t(i,j,k,nstp,isldon) , 0.0d0 )
      cobalt%f_srdon(i,j,k)      = max( t(i,j,k,nstp,isrdon) , 0.0d0 )
      bact(1)%f_n(i,j,k)         = max( t(i,j,k,nstp,inbact) , 0.0d0 )
      cobalt%f_nh4(i,j,k)        = max( t(i,j,k,nstp,inh4) ,   0.0d0 )
      cobalt%f_no3(i,j,k)        = max( t(i,j,k,nstp,ino3) ,   0.0d0 )
      cobalt%f_ndet(i,j,k)       = max( t(i,j,k,nstp,indet) ,  0.0d0 )
#ifdef COBALT_MINERALS
      ! Biogenic Minerals and Lithogenic Materials
      cobalt%f_sio4(i,j,k)       = max( t(i,j,k,nstp,isio4) ,       0.0d0 )
      cobalt%f_silg(i,j,k)       = max( t(i,j,k,nstp,isilg) ,       0.0d0 )
      cobalt%f_sidet(i,j,k)      = max( t(i,j,k,nstp,isidet) ,      0.0d0 )
      cobalt%f_cadet_calc(i,j,k) = max( t(i,j,k,nstp,icadet_calc) , 0.0d0 )
      cobalt%f_cadet_arag(i,j,k) = max( t(i,j,k,nstp,icadet_arag) , 0.0d0 )
      cobalt%f_lith(i,j,k)       = max( t(i,j,k,nstp,ilith) ,       0.0d0 )
      cobalt%f_lithdet(i,j,k)    = max( t(i,j,k,nstp,ilithdet) ,    0.0d0 )
#endif
#ifdef COBALT_PHOSPHORUS
      ! Phosporus Dynamics
      cobalt%f_ldop(i,j,k)       = max( t(i,j,k,nstp,ildop) ,  0.0d0 )
      cobalt%f_sldop(i,j,k)      = max( t(i,j,k,nstp,isldop) , 0.0d0 )
      cobalt%f_srdop(i,j,k)      = max( t(i,j,k,nstp,isrdop) , 0.0d0 )
      cobalt%f_po4(i,j,k)        = max( t(i,j,k,nstp,ipo4) ,   0.0d0 )
      cobalt%f_pdet(i,j,k)       = max( t(i,j,k,nstp,ipdet) ,  0.0d0 )
#endif
#ifdef COBALT_IRON
      ! Iron Dynamics
      phyto(SMALL)%f_fe(i,j,k)   = max( t(i,j,k,nstp,ifesm) ,  0.0d0 )
      phyto(DIAZO)%f_fe(i,j,k)   = max( t(i,j,k,nstp,ifedi) ,  0.0d0 )
      phyto(LARGE)%f_fe(i,j,k)   = max( t(i,j,k,nstp,ifelg) ,  0.0d0 )
      cobalt%f_fed(i,j,k)        = max( t(i,j,k,nstp,ifed) ,   0.0d0 )
      cobalt%f_fedet(i,j,k)      = max( t(i,j,k,nstp,ifedet) , 0.0d0 )
#endif
#ifdef COBALT_CARBON
      ! Oxygen, Carbon and Alkalinity
      cobalt%f_o2(i,j,k)         = max( t(i,j,k,nstp,io2) ,  0.0d0 )
      cobalt%f_dic(i,j,k)        = max( t(i,j,k,nstp,idic) , 0.0d0 )
      cobalt%f_alk(i,j,k)        = max( t(i,j,k,nstp,ialk) , 0.0d0 )
#endif
  ENDDO ; ENDDO ; ENDDO

! diagnostic tracers that are passed between time steps (except chlorophyll)

#ifndef COBALT_DO_NOTHING

#ifdef BENTHIC
! RD dev notes :
! RE DO THIS !!!
! bottom fluxes are passed into ROMS in the bt array, it is well suited for 2d 
! variables, is not advected/diffused and has a restart implementation
! irr_mem, hotal and co3_ion are in DiaBio3d (3d diagnostics variables) and
! don't have this restart implementation

! DO WE EVER NEED THOSE ????

  DO j=Jstr,Jend
    DO i=Istr,Iend
    cobalt%f_cased(i,j,1)          = max( bt(i,j,1,nstp,icased) , 0.0d0 )
    ! fluxes don't need to be positive definite
!    cobalt%f_cadet_arag_btf(i,j,1) = bt(i,j,1,nstp,icadet_arag_btf)
!    cobalt%f_cadet_calc_btf(i,j,1) = bt(i,j,1,nstp,icadet_calc_btf)
!    cobalt%f_ndet_btf(i,j,1)       = bt(i,j,1,nstp,indet_btf)
!    cobalt%f_pdet_btf(i,j,1)       = bt(i,j,1,nstp,ipdet_btf)
!    cobalt%f_sidet_btf(i,j,1)      = bt(i,j,1,nstp,isidet_btf)

    ENDDO
  ENDDO
  i=overflow ; j=overflow 
#endif

! RD dev notes:
! irradiance memory read from restart to avoid glitch at the begining of a new
! job
  DO k=1,UBk
    DO j=Jstr,Jend
      DO i=Istr,Iend
      cobalt%f_irr_mem(i,j,k)      = max( obgc(i,j,k,nstp,ioirr_mem) , 0.0d0 )
      ENDDO
    ENDDO
  ENDDO
  i=overflow ; j=overflow ; k=overflow
!
! RD dev notes : dust deposition to iron and lith surface source
  DO j=Jstr,Jend
    DO i=Istr,Iend

       iron_dust_src(i,j,UBk) = soluble_fe(i,j) / (rho0 * Hz(i,j,UBk) )
       ! units mol/kg/s       =  [mol.m-2.s-1]  / [kg.m-3] * [m]

       lith_dust_src(i,j,UBk) = mineral_fe(i,j) / (rho0 * Hz(i,j,UBk) )
       ! units g/kg/s         = [g.m-2.s-1]     / [kg.m-3] * [m]

#ifdef COBALT_CONSERVATION_TEST
        iron_dust_src(i,j,UBk) = 0.0d0
        lith_dust_src(i,j,UBk) = 0.0d0
#endif
#ifdef DIAGNOSTICS_BIO
        DiaBio3d(i,j,UBk,ife_bulk_flx) = DiaBio3d(i,j,UBk,ife_bulk_flx) + iron_dust_src(i,j,UBk) * dt(ng)
#endif

    ENDDO
  ENDDO
  i=overflow ; j=overflow
!
! RD dev notes : dust deposition to iron and lith subsurface source
! set to zero for now, some other models add iron and lith in the upper 600
! meters below the surface. Here for future developments
   DO k=1,UBk-1
     DO j=Jstr,Jend
       DO i=Istr,Iend

          iron_dust_src(i,j,k) = 0.0d0
          lith_dust_src(i,j,k) = 0.0d0
#ifdef COBALT_CONSERVATION_TEST
          iron_dust_src(i,j,k) = 0.0d0
          lith_dust_src(i,j,k) = 0.0d0
#endif
#ifdef DIAGNOSTICS_BIO
          DiaBio3d(i,j,k,ife_bulk_flx) = DiaBio3d(i,j,k,ife_bulk_flx) + iron_dust_src(i,j,k) * dt(ng)
#endif

       ENDDO
     ENDDO
   ENDDO
   i=overflow ; j=overflow ; k=overflow

!
!-----------------------------------------------------------------------------------
! 1: Phytoplankton growth and nutrient uptake calculations
!-----------------------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------------------
! 1.1: Nutrient Limitation Calculations
!-----------------------------------------------------------------------------------
!
! Calculate iron cell quota
!

  DO nphyto=1,NUM_PHYTO
     phyto(nphyto)%q_fe_2_n(:,:,:)=max(0.0d0,phyto(nphyto)%f_fe(:,:,:) / &
   & max(epsln,phyto(nphyto)%f_n(:,:,:)))
 
     phyto(nphyto)%q_p_2_n(:,:,:) = phyto(nphyto)%p_2_n_static
  ENDDO
  nphyto=overflow

!
! N limitation with NH4 inhibition after Frost and Franzen (1992)
!

  DO nphyto=2,NUM_PHYTO
     phyto(nphyto)%no3lim(:,:,:) = cobalt%f_no3(:,:,:) /   &
   & ( (phyto(nphyto)%k_no3 + cobalt%f_no3(:,:,:))     *   &
   & (1.0d0 + cobalt%f_nh4(:,:,:) / phyto(nphyto)%k_nh4) )

     phyto(nphyto)%nh4lim(:,:,:) = cobalt%f_nh4(:,:,:) /   &
   & (phyto(nphyto)%k_nh4 + cobalt%f_nh4(:,:,:))
  ENDDO
  nphyto=overflow

!
! O2 inhibition term for diazotrophs
!

  nphyto=DIAZO
  phyto(nphyto)%o2lim(:,:,:) = (1.0d0 -               &
  & cobalt%f_o2(:,:,:)**cobalt%o2_inhib_Di_pow /      &
  & (cobalt%f_o2(:,:,:)**cobalt%o2_inhib_Di_pow +     &
  & cobalt%o2_inhib_Di_sat**cobalt%o2_inhib_Di_pow))

  nphyto=overflow
!
! SiO4, PO4 and Fe uptake limitation with Michaelis-Mentin 
!

  phyto(LARGE)%silim(:,:,:) = cobalt%f_sio4(:,:,:) /  &
  & (phyto(LARGE)%k_sio4 + cobalt%f_sio4(:,:,:))

  DO nphyto=1,NUM_PHYTO

     phyto(nphyto)%po4lim(:,:,:) = cobalt%f_po4(:,:,:) / &
   & (phyto(nphyto)%k_po4 + cobalt%f_po4(:,:,:))

     phyto(nphyto)%felim(:,:,:) = cobalt%f_fed(:,:,:) / &
   & (phyto(nphyto)%k_fed + cobalt%f_fed(:,:,:))

     phyto(nphyto)%def_fe(:,:,:)=phyto(nphyto)%q_fe_2_n(:,:,:)**2.0d0 / &
   & (phyto(nphyto)%k_fe_2_n**2.0d0+phyto(nphyto)%q_fe_2_n(:,:,:)**2.0d0)

  ENDDO
  nphyto=overflow

!
! Calculate nutrient limitation based on the most limiting nutrient (liebig_lim)
!

  nphyto=DIAZO
  phyto(nphyto)%liebig_lim(:,:,:) = phyto(nphyto)%o2lim(:,:,:) * &
  & min( phyto(nphyto)%po4lim(:,:,:), phyto(nphyto)%def_fe(:,:,:) )

  DO nphyto=2,NUM_PHYTO
     phyto(nphyto)%liebig_lim(:,:,:) = min( phyto(nphyto)%no3lim(:,:,:)+ &
   &                                        phyto(nphyto)%nh4lim(:,:,:), &
   &                                        phyto(nphyto)%po4lim,        &
   &                                        phyto(nphyto)%def_fe )
  ENDDO
  nphyto=overflow

!
!-----------------------------------------------------------------------
! 1.2: Light Limitation/Growth Calculations
!-----------------------------------------------------------------------
!
! RD dev notes : copy fractional decay for red and blue bands previously
! computed in optic_manizza

#ifdef OPTIC_MANIZZA
   DO k=0,UBk
     DO j=Jstr,Jend
       DO i=Istr,Iend
          swdk3(i,j,k) = decayW(i,j,k,3) + decayW(i,j,k,4) ! use only visible bands
       END DO
     END DO
   ENDDO
   i=overflow ; j=overflow ; k=overflow
#else
???
# maybe put something with lmd_swfrac but does it make sense really ?
#endif

#ifdef DIAGNOSTICS_BIO
   DO k=1,UBk
     DO j=Jstr,Jend
       DO i=Istr,Iend
          DiaBio3d(i,j,k,iswdk) = DiaBio3d(i,j,k,iswdk) + swdk3(i,j,k)
       END DO
     END DO
   ENDDO
   i=overflow ; j=overflow ; k=overflow
#endif


   !------------ Compute mixed layer depth -----------------------------

   ! RD dev notes :
   ! array initialization :
   ! mxl bottom level set to model's bottom level -1
   ! mxl_depth is set to the depth of the water column
   ! rho_ref is set to the surface value
   DO j=Jstr,Jend
     DO i=Istr,Iend
        mxl_blev(i,j)  = 1
        mxl_depth(i,j) = -h(i,j)
        rho_ref(i,j)   = rho(i,j,UBk)
     ENDDO
   ENDDO
   i=overflow ; j=overflow

   DO k=0,UBk
     DO j=Jstr,Jend
       DO i=Istr,Iend
          FM(i,j,k) = z_w(i,j,UBk) - z_w(i,j,k)
       ENDDO
     ENDDO
   ENDDO

   ! update the reference density to the density at z = rho_ref_depth
   DO k=1,UBk
     DO j=Jstr,Jend
       DO i=Istr,Iend
          IF (FM(i,j,k-1)  > rho_ref_depth ) rho_ref(i,j) = rho(i,j,k)
       ENDDO
     ENDDO
   ENDDO
   i=overflow ; j=overflow ; k=overflow

   ! find which is the level just below the mixed layer and update the array
   ! bottom to top, if rho > rho_ref + 0.03 increase k until not true anymore
   DO k=1,UBk
     DO j=Jstr,Jend
       DO i=Istr,Iend
          IF (rho(i,j,k) > rho_ref(i,j) + rho_crit)   mxl_blev(i,j) = k
       ENDDO
     ENDDO
   ENDDO
   i=overflow ; j=overflow ; k=overflow

   ! find the real mxl with linear interpolation 
   DO j=Jstr,Jend
     DO i=Istr,Iend

      IF ( mxl_blev(i,j) > 1) THEN
          ! this would be the depth if we consider all the level
          mxl_depth(i,j) =  FM(i,j,mxl_blev(i,j))

          ! we approach the real mixed later depth with linear interp between
          ! vertical levels
          !drho_crit = rho_ref(i,j) + rho_crit - rho(i,j,mxl_blev(i,j)+1)

!          mxl_subgrid(i,j) =  drho_crit *                                &
!  &       ( ( z_w(i,j,mxl_blev(i,j)+1) - z_w(i,j,mxl_blev(i,j)) ) /      &
!  &         max( rho(i,j,mxl_blev(i,j)) - rho(i,j,mxl_blev(i,j)+1), 1.0e-9 ) )

          ! mixed layer depth is then shallower than depth of the bottom level
!          mxl_depth(i,j) = FM(i,j,mxl_blev(i,j)) - mxl_subgrid(i,j)

      ELSE
          mxl_depth(i,j) = h(i,j)
      ENDIF

     ENDDO
   ENDDO
   i=overflow ; j=overflow

! I should use something like in pre_step3d to treat heating in layers with fluxes
! at interfaces.

!   DO j=Jstr,Jend
!     DO i=Istr,Iend
!      ! init the irradiance in the mixed layer and last level of mixed layer
!      tmp_irrad_ML = 0.0d0 ; kblt = mxl_blev(i,j) + 1
!
!      DO k=1,N(ng)
!         ! light penetration decay function at rho-point
!         tmp_decay(i,j,k) = 0.5 * ( swdk3(i,j,k-1) + swdk3(i,j,k) )
!         ! compute the instant irradiance 
!         cobalt%irr_inst(i,j,k) = rho0 * Cp * srflx(i,j) * tmp_decay(i,j,k)
!
!         ! RD test :
!         !cobalt%irr_inst(i,j,k) = cobalt%irr_inst(i,j,k) * 0.5
!
!         ! integrate along depth irradiances for all levels in mixed layer
!         IF ( k .ge. kblt ) THEN
!            tmp_irrad_ML = tmp_irrad_ML + ( cobalt%irr_inst(i,j,k) * Hz(i,j,k) )
!         ENDIF
!         ! set base value of irr_mix to irr_inst (Charles Stock told me to do so)
!         cobalt%irr_mix(i,j,k) = cobalt%irr_inst(i,j,k)
!      ENDDO
!      ! update irr_mix in the mix layer to the integrated value
!      cobalt%irr_mix(i,j,kblt:UBk) = tmp_irrad_ML * rmask(i,j) / mxl_depth(i,j)
!
!     ENDDO
!   ENDDO
!   i=overflow ; j=overflow

   DO j=Jstr,Jend
     DO i=Istr,Iend
      DO k=1,N(ng)
         ! compute the instant irradiance in layer k as a difference of fluxes
         ! between interfaces k and k-1 (like in pre_step3d)
         cobalt%irr_inst(i,j,k) = rho0 * Cp * srflx(i,j) * &
       & ( swdk3(i,j,k) - swdk3(i,j,k-1) )

         ! set base value of irr_mix to irr_inst (Charles Stock told me to do so)
         cobalt%irr_mix(i,j,k) = cobalt%irr_inst(i,j,k)

      ENDDO
     ENDDO
   ENDDO
   i=overflow ; j=overflow ; k=overflow

   DO j=Jstr,Jend
     DO i=Istr,Iend

      ! init the irradiance in the mixed layer and last level of mixed layer
      ! we are working on a layer-based integration hence +1 for last level of
      ! mixed layer
      tmp_irrad_ML = 0.0d0 ; kblt = mxl_blev(i,j) + 1

      ! integrate irradiance in last level of mixed layer with correction for
      ! depth
!      k=kblt
!      tmp_irrad_ML = tmp_irrad_ML + cobalt%irr_inst(i,j,k) * &
!    &                max( Hz(i,j,k) - mxl_subgrid(i,j), 0.0d0 )
      
      DO k=kblt,N(ng)
!      DO k=kblt+1,N(ng)
         tmp_irrad_ML = tmp_irrad_ML + cobalt%irr_inst(i,j,k) * Hz(i,j,k)
      ENDDO

!      DO k=kblt+1,N(ng)
      DO k=kblt,N(ng)
         cobalt%irr_mix(i,j,k) = tmp_irrad_ML * rmask(i,j) / mxl_depth(i,j)
      ENDDO

     ENDDO
   ENDDO
   i=overflow ; j=overflow



#ifdef DEBUG_COBALT
IF( Master ) WRITE(stdout,*) '>>>   max srflx is = ', MAXVAL(srflx)
IF( Master ) WRITE(stdout,*) '>>>   max swdk is = ', MAXVAL(swdk3)
IF( Master ) WRITE(stdout,*) '>>>   max irr_inst is = ', MAXVAL(cobalt%irr_inst)
IF( Master ) WRITE(stdout,*) '>>>   max irr_mem is = ', MAXVAL(cobalt%f_irr_mem)
IF( Master ) WRITE(stdout,*) '>>>   max irr_mix is = ', MAXVAL(cobalt%irr_mix)
#endif 

  !
  ! Calculate the temperature limitation (expkT) and the time integrated
  ! irradiance (f_irr_mem) to which the Chl:C ratio responds (~24 hours)
  !

  DO k=1,UBk
    DO j=Jstr,Jend
      DO i=Istr,Iend
         cobalt%expkT(i,j,k) = exp(cobalt%kappa_eppley * cobalt%f_temp(i,j,k))
         cobalt%f_irr_mem(i,j,k) = (cobalt%f_irr_mem(i,j,k) + &
       & (cobalt%irr_mix(i,j,k)  - cobalt%f_irr_mem(i,j,k)) * &
       & min(1.0d0,cobalt%gamma_irr_mem * dt(ng))) * rmask(i,j)
      ENDDO
    ENDDO
  ENDDO
  i=overflow ; j=overflow ; k=overflow

!
! Phytoplankton growth rate calculation based on Geider et al. (1997)
!

  DO k=1,UBk
    DO j=Jstr,Jend
      DO i=Istr,Iend

! RD dev notes : charlie notes suggested to get rid of gross_prim_prod
!         cobalt%gross_prim_prod(i,j,k) = 0.0
         cobalt%f_chl(i,j,k) = 0.0d0

         DO nphyto=1,NUM_PHYTO
            P_C_m = phyto(nphyto)%liebig_lim(i,j,k) *                    &
          & phyto(nphyto)%P_C_max * cobalt%expkT(i,j,k) + epsln

            phyto(nphyto)%theta(i,j,k) = (phyto(nphyto)%thetamax -       &
          & cobalt%thetamin) / (1.0 + phyto(nphyto)%thetamax *           &
          & phyto(nphyto)%alpha * cobalt%f_irr_mem(i,j,k)*0.5 / P_C_m) + &
          & cobalt%thetamin

            cobalt%f_chl(i,j,k) = cobalt%f_chl(i,j,k) + cobalt%c_2_n *   &
          & 12.0e6 * phyto(nphyto)%theta(i,j,k) * phyto(nphyto)%f_n(i,j,k)

            phyto(nphyto)%irrlim(i,j,k) = (1.0-exp(-phyto(nphyto)%alpha* &
          & cobalt%irr_inst(i,j,k) * phyto(nphyto)%theta(i,j,k) / P_C_m))

          ! calculate the growth rate
            phyto(nphyto)%mu(i,j,k) = P_C_m / (1.0 + cobalt%zeta) *      &
          & phyto(nphyto)%irrlim(i,j,k) - cobalt%expkT(i,j,k)     *      &
          & phyto(nphyto)%bresp * phyto(nphyto)%f_n(i,j,k)        /      &
          & (refuge_conc + phyto(nphyto)%f_n(i,j,k))

!            cobalt%gross_prim_prod(i,j,k)=cobalt%gross_prim_prod(i,j,k)+ &
!          & P_C_m*phyto(nphyto)%irrlim(i,j,k) * phyto(nphyto)%f_n(i,j,k)

          ! Negative growth assumed to go to cell death rather than respiration
          ! (see manual)
!            cobalt%net_prim_prod(i,j,k)=max(phyto(nphyto)%mu(i,j,k),0.)* &
!          & phyto(nphyto)%f_n(i,j,k)
         ENDDO
         nphyto=overflow

!      cobalt%gross_prim_prod(i,j,k) = cobalt%gross_prim_prod(i,j,k) *    &
!    & cobalt%c_2_n*spery

!      cobalt%net_prim_prod(i,j,k) = cobalt%net_prim_prod(i,j,k) *        &
!    & cobalt%c_2_n*spery

      ENDDO
    ENDDO
  ENDDO
  i=overflow ; j=overflow ; k=overflow

!-----------------------------------------------------------------------
! 1.3: Nutrient uptake calculations 
!-----------------------------------------------------------------------
!
! Uptake of nitrate and ammonia
!

  nphyto=DIAZO
  
   phyto(nphyto)%juptake_n2(:,:,:) = max( 0.0,                           &
 & (1.0 - phyto(LARGE)%no3lim(:,:,:) - phyto(LARGE)%nh4lim(:,:,:)) *     &
 & phyto(nphyto)%mu(:,:,:) * phyto(nphyto)%f_n(:,:,:) )

   phyto(nphyto)%juptake_nh4(:,:,:)=max(0.0,phyto(LARGE)%nh4lim(:,:,:) * &
 & phyto(nphyto)%mu(:,:,:)*phyto(nphyto)%f_n(:,:,:))

   phyto(nphyto)%juptake_no3(:,:,:)=max(0.0,phyto(LARGE)%no3lim(:,:,:) * &
 & phyto(nphyto)%mu(:,:,:)*phyto(nphyto)%f_n(:,:,:))

   ! If growth is negative, net remineralization of organic material
   phyto(nphyto)%juptake_nh4(:,:,:) = phyto(nphyto)%juptake_nh4(:,:,:) + &
 & min(0.0,phyto(nphyto)%mu(:,:,:)  * phyto(nphyto)%f_n(:,:,:))

   phyto(nphyto)%jprod_n(:,:,:) = phyto(nphyto)%juptake_nh4(:,:,:)     + &
 & phyto(nphyto)%juptake_no3(:,:,:) + phyto(nphyto)%juptake_n2(:,:,:)

   DO nphyto=2,NUM_PHYTO

      phyto(nphyto)%juptake_no3(:,:,:)=max(0.0,phyto(nphyto)%mu(:,:,:) * &
    & phyto(nphyto)%f_n(:,:,:) * phyto(nphyto)%no3lim(:,:,:)           / &
    & (phyto(nphyto)%no3lim(:,:,:) + phyto(nphyto)%nh4lim(:,:,:)+epsln) )
      
      phyto(nphyto)%juptake_nh4(:,:,:)=max(0.0,phyto(nphyto)%mu(:,:,:) * &
    & phyto(nphyto)%f_n(:,:,:) * phyto(nphyto)%nh4lim(:,:,:)           / &
    & (phyto(nphyto)%no3lim(:,:,:) + phyto(nphyto)%nh4lim(:,:,:)+epsln) )

      ! If growth is negative, net remineralization of organic material
      phyto(nphyto)%juptake_nh4(:,:,:)=phyto(nphyto)%juptake_nh4(:,:,:)+ &
    & min(0.0,phyto(nphyto)%mu(:,:,:) *phyto(nphyto)%f_n(:,:,:))

      phyto(nphyto)%jprod_n(:,:,:) = phyto(nphyto)%juptake_nh4(:,:,:) +  &
    & phyto(nphyto)%juptake_no3(:,:,:)

   ENDDO
   nphyto=overflow

!
! Phosphorous uptake
! 

  nphyto=DIAZO
  phyto(nphyto)%juptake_po4(:,:,:) = (phyto(nphyto)%juptake_n2(:,:,:)  + &
  & phyto(nphyto)%juptake_nh4(:,:,:)+phyto(nphyto)%juptake_no3(:,:,:)) * &
  & phyto(nphyto)%p_2_n_static

  DO nphyto=2,NUM_PHYTO
     phyto(nphyto)%juptake_po4(:,:,:)=(phyto(nphyto)%juptake_no3(:,:,:)+ & 
     phyto(nphyto)%juptake_nh4(:,:,:)) * phyto(nphyto)%p_2_n_static
  ENDDO
  nphyto=overflow

!
! Iron uptake
! 

  DO k=1,UBk
    DO j=Jstr,Jend
      DO i=Istr,Iend

        DO nphyto=1,NUM_PHYTO
          IF (phyto(nphyto)%q_fe_2_n(i,j,k).lt.phyto(nphyto)%fe_2_n_max) THEN
             phyto(nphyto)%juptake_fe(i,j,k) = phyto(nphyto)%P_C_max * &
           & cobalt%expkT(i,j,k) * phyto(nphyto)%f_n(i,j,k)          * &
           & phyto(nphyto)%felim(i,j,k) * cobalt%fe_2_n_upt_fac
          ELSE
             phyto(nphyto)%juptake_fe(i,j,k) = 0.0
          ENDIF
        ENDDO
        nphyto=overflow

      ENDDO
    ENDDO
  ENDDO 
  i=overflow ; j=overflow ; k=overflow

!
! Silicate uptake
!

   cobalt%nlg_diatoms(:,:,:) = phyto(LARGE)%f_n(:,:,:) *                 &
 & phyto(LARGE)%silim(:,:,:)

   cobalt%q_si_2_n_lg_diatoms(:,:,:) = cobalt%f_silg(:,:,:) /            &
 & (cobalt%nlg_diatoms(:,:,:) + epsln)

   phyto(LARGE)%juptake_sio4(:,:,:)=max(phyto(LARGE)%juptake_no3(:,:,:)+ &
 & phyto(LARGE)%juptake_nh4(:,:,:),0.0) * phyto(LARGE)%silim(:,:,:)    * &
 & phyto(LARGE)%silim(:,:,:) * phyto(LARGE)%si_2_n_max

! CAS: set q_si_2_n values for each of the phyto groups for consumption calculations
! Note that this is si_2_n in large phytoplankton pool, not in diatoms themselves 

! RD : forgot this line (oups...)
   phyto(LARGE)%q_si_2_n(:,:,:) = cobalt%f_silg(:,:,:)/(phyto(LARGE)%f_n(:,:,:)+epsln)

!
!-----------------------------------------------------------------------
! 2: Bacterial Growth and Uptake Calculations 
!-----------------------------------------------------------------------
!
!
! calculate an effective maximum ldon uptake rate (at 0 deg. C) for bacteria
! from specified values of bact(1)%gge_max, bact(1)%mu_max and bact(1)%bresp
!

  vmax_bact = (1.0d0/bact(1)%gge_max)*(bact(1)%mu_max + bact(1)%bresp)

  DO k=1,Ubk
    DO j=Jstr,Jend
      DO i=Istr,Iend

         bact(1)%temp_lim(i,j,k) = exp(bact(1)%ktemp*cobalt%f_temp(i,j,k))

         bact_ldon_lim = cobalt%f_ldon(i,j,k)/(bact(1)%k_ldon + cobalt%f_ldon(i,j,k))

         bact(1)%juptake_ldon(i,j,k)=vmax_bact*bact(1)%temp_lim(i,j,k) * &
       & bact_ldon_lim * bact(1)%f_n(i,j,k)

         bact_uptake_ratio = ( cobalt%f_ldop(i,j,k)/max(cobalt%f_ldon(i,j,k),epsln) )

         bact(1)%juptake_ldop(i,j,k) = bact(1)%juptake_ldon(i,j,k)*bact_uptake_ratio

         IF (bact_uptake_ratio.lt.bact(1)%q_p_2_n) THEN

            bact(1)%jprod_n(i,j,k) = bact(1)%gge_max                   * &
          & bact(1)%juptake_ldop(i,j,k)*16.0 - bact(1)%f_n(i,j,k)      / &
          & (refuge_conc + bact(1)%f_n(i,j,k))*bact(1)%temp_lim(i,j,k) * &
          & bact(1)%bresp * bact(1)%f_n(i,j,k)

         ELSE

            bact(1)%jprod_n(i,j,k) = bact(1)%gge_max                   * &
          & bact(1)%juptake_ldon(i,j,k) - bact(1)%f_n(i,j,k)           / &
          & (refuge_conc + bact(1)%f_n(i,j,k))*bact(1)%temp_lim(i,j,k) * &
          & bact(1)%bresp * bact(1)%f_n(i,j,k)
         ENDIF

      ENDDO
    ENDDO
  ENDDO
  i=overflow ; j=overflow ; k=overflow

!
!-----------------------------------------------------------------------
! 3: Plankton foodweb dynamics
!-----------------------------------------------------------------------
!
!
! 3.1 Plankton foodweb dynamics: consumption by zooplankton and higher predators
!
!
! Set-up local matrices for calculating zooplankton ingestion of
! multiple prey types.  The rows are consumers (i.e., NUM_ZOO zooplankton
! groups), the columns are food sources (i.e., NUM_PREY potential food sources)
!
! ipa_matrix = the innate prey availability matrix
! pa_matrix = prey availability matrix after accounting for switching 
! ingest_matrix = ingestion matrix
! tot_prey = total prey available to predator m 
!
! The definition of predator-prey matrices is intended to allow for
! efficient experimentation with predator-prey interconnections.
! However, we are still working to reduce the runtime required to
! include this feature.  The matrix structures are thus included,
! but the standard COBALT interactions have been hard-coded such
! that changing linkages requires changing the prey availability
! values and adding additional code to handle the new linkages.
!
! With regard to stoichiometry, the primary ingestion calculations
! (i.e., those within the i, j, k loops) are coded to allow for 
! variable stoichiometry.  Several sections of the code corresponding
! to predator-prey and other linkages not in included in the
! default COBALT parameterizations have been commented out to
! avoid unnecessary calculations.
!

    DO m=1,NUM_ZOO 
       ipa_matrix(m,1) = zoo(m)%ipa_diaz
       ipa_matrix(m,2) = zoo(m)%ipa_lgp
       ipa_matrix(m,3) = zoo(m)%ipa_smp
       ipa_matrix(m,4) = zoo(m)%ipa_bact
       ipa_matrix(m,5) = zoo(m)%ipa_smz
       ipa_matrix(m,6) = zoo(m)%ipa_mdz
       ipa_matrix(m,7) = zoo(m)%ipa_lgz
       ipa_matrix(m,8) = zoo(m)%ipa_det
       tot_prey(m)     = 0.0
       DO nprey=1,NUM_PREY 
           pa_matrix(m,nprey)     = 0.0
           ingest_matrix(m,nprey) = 0.0
       ENDDO 
    ENDDO 
    m=overflow ; nprey=overflow

!
! Set-up local matrices for calculating higher predator ingestion
! of multiple prey types
!

    hp_ipa_vec(1) = cobalt%hp_ipa_diaz
    hp_ipa_vec(2) = cobalt%hp_ipa_lgp
    hp_ipa_vec(3) = cobalt%hp_ipa_smp
    hp_ipa_vec(4) = cobalt%hp_ipa_bact
    hp_ipa_vec(5) = cobalt%hp_ipa_smz
    hp_ipa_vec(6) = cobalt%hp_ipa_mdz
    hp_ipa_vec(7) = cobalt%hp_ipa_lgz
    hp_ipa_vec(8) = cobalt%hp_ipa_det
    tot_prey_hp   = 0.0
    DO nprey=1,NUM_PREY  
       hp_pa_vec(nprey)     = 0.0                  
       hp_ingest_vec(nprey) = 0.0              
    ENDDO 
    nprey=overflow

! 
! Set all static stoichiometric ratios outside k,j,i loop
!

    prey_p2n_vec(1) = phyto(DIAZO)%p_2_n_static
    prey_p2n_vec(2) = phyto(LARGE)%p_2_n_static
    prey_p2n_vec(3) = phyto(SMALL)%p_2_n_static
    prey_p2n_vec(4) = bact(1)%q_p_2_n
    prey_p2n_vec(5) = zoo(1)%q_p_2_n
    prey_p2n_vec(6) = zoo(2)%q_p_2_n
    prey_p2n_vec(7) = zoo(3)%q_p_2_n

    prey_fe2n_vec(4) = 0.0
    prey_fe2n_vec(5) = 0.0
    prey_fe2n_vec(6) = 0.0
    prey_fe2n_vec(7) = 0.0

    prey_si2n_vec(1) = 0.0
    prey_si2n_vec(3) = 0.0
    prey_si2n_vec(4) = 0.0
    prey_si2n_vec(5) = 0.0
    prey_si2n_vec(6) = 0.0
    prey_si2n_vec(7) = 0.0


!RD debug : does not breaks anything, just a diag
DO nzoo=1,NUM_ZOO
  DO m=1,NUM_PREY
     ingest_matrix_glo(nzoo,m) = 0.
  ENDDO
ENDDO


  DO k=1,UBk
    DO j=Jstr,Jend
      DO i=Istr,Iend

       !
       ! 3.1.1: Calculate zooplankton ingestion fluxes
       !

       ! Prey vectors for ingestion and loss calculations 
       ! (note: ordering of phytoplankton must be consistent with
       !  DIAZO, LARGE, SMALL ordering inherited from TOPAZ)
       !
       prey_vec(1) = max(phyto(DIAZO)%f_n(i,j,k) - refuge_conc,0.0d0)
       prey_vec(2) = max(phyto(LARGE)%f_n(i,j,k) - refuge_conc,0.0d0)
       prey_vec(3) = max(phyto(SMALL)%f_n(i,j,k) - refuge_conc,0.0d0)
       prey_vec(4) = max(bact(1)%f_n(i,j,k)      - refuge_conc,0.0d0)
       prey_vec(5) = max(zoo(1)%f_n(i,j,k)       - refuge_conc,0.0d0)
       prey_vec(6) = max(zoo(2)%f_n(i,j,k)       - refuge_conc,0.0d0)
       prey_vec(7) = max(zoo(3)%f_n(i,j,k)       - refuge_conc,0.0d0)
       prey_vec(8) = max(cobalt%f_ndet(i,j,k)    - refuge_conc,0.0d0)
       ! 
       ! Set dynamic stoichiometric rations inside k,j,i loop
       prey_p2n_vec(8)  = cobalt%f_pdet(i,j,k)/(cobalt%f_ndet(i,j,k)+epsln)
       prey_fe2n_vec(1) = phyto(DIAZO)%q_fe_2_n(i,j,k)
       prey_fe2n_vec(2) = phyto(LARGE)%q_fe_2_n(i,j,k)
       prey_fe2n_vec(3) = phyto(SMALL)%q_fe_2_n(i,j,k)
       prey_fe2n_vec(8) = cobalt%f_fedet(i,j,k)/(cobalt%f_ndet(i,j,k)+epsln)
       prey_si2n_vec(2) = phyto(LARGE)%q_si_2_n(i,j,k)
       prey_si2n_vec(8) = cobalt%f_sidet(i,j,k)/(cobalt%f_ndet(i,j,k)+epsln)

       !
       ! Calculate zooplankton ingestion
       !
       ! Small zooplankton (m = 1) consuming small phytoplankton (3) and
       ! bacteria (4).  sw_fac_denom is the denominator of the abundance-
       ! based switching factor, tot_prey is the total available prey 
       ! after accounting for switching.
       !

       m = 1 
       zoo(m)%temp_lim(i,j,k) = exp(zoo(m)%ktemp*cobalt%f_temp(i,j,k)) 

       sw_fac_denom = (ipa_matrix(m,3)*prey_vec(3))**zoo(m)%nswitch +    &
     &                (ipa_matrix(m,4)*prey_vec(4))**zoo(m)%nswitch

       pa_matrix(m,3) = ipa_matrix(m,3)*                                 &
     &                  ((ipa_matrix(m,3)*prey_vec(3))**zoo(m)%nswitch / &
     &                    (sw_fac_denom+epsln))**(1.0/zoo(m)%mswitch)

       pa_matrix(m,4) = ipa_matrix(m,4)*                                 &
     &                  ((ipa_matrix(m,4)*prey_vec(4))**zoo(m)%nswitch / &
     &                    (sw_fac_denom+epsln))**(1.0/zoo(m)%mswitch)

       tot_prey(m) = pa_matrix(m,3) * prey_vec(3) + pa_matrix(m,4) *     &
     &                prey_vec(4)

       ingest_matrix(m,3) = zoo(m)%temp_lim(i,j,k) * zoo(m)%imax *       &
     &                pa_matrix(m,3)* prey_vec(3)*zoo(m)%f_n(i,j,k) /    &
     &                (zoo(m)%ki+tot_prey(m))


       ingest_matrix(m,4) = zoo(m)%temp_lim(i,j,k) * zoo(m)%imax *       &
     &                pa_matrix(m,4) * prey_vec(4) * zoo(m)%f_n(i,j,k) / &
     &                (zoo(m)%ki+tot_prey(m))


       zoo(m)%jingest_n(i,j,k) = ingest_matrix(m,3) + ingest_matrix(m,4)
       zoo(m)%jingest_p(i,j,k) = ingest_matrix(m,3)*prey_p2n_vec(3) +    &
     &                           ingest_matrix(m,4)*prey_p2n_vec(4)
       zoo(m)%jingest_fe(i,j,k)= ingest_matrix(m,3)*prey_fe2n_vec(3)

       !
       ! Medium zooplankton (m = 2) consuming diazotrophs (1), large
       ! phytoplankton (2), and small zooplankton (5) 
       !

       m = 2 
       zoo(m)%temp_lim(i,j,k) = exp(zoo(m)%ktemp*cobalt%f_temp(i,j,k))
       sw_fac_denom = (ipa_matrix(m,1)*prey_vec(1))**zoo(m)%nswitch +    &
     &                (ipa_matrix(m,2)*prey_vec(2))**zoo(m)%nswitch +    &
     &                (ipa_matrix(m,5)*prey_vec(5))**zoo(m)%nswitch
       pa_matrix(m,1) = ipa_matrix(m,1)*                                 &
     &                 ( (ipa_matrix(m,1)*prey_vec(1))**zoo(m)%nswitch / &
     &                   (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,2) = ipa_matrix(m,2)*                                 & 
     &                 ( (ipa_matrix(m,2)*prey_vec(2))**zoo(m)%nswitch / &
     &                   (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,5) = ipa_matrix(m,5)*                                 & 
     &                 ( (ipa_matrix(m,5)*prey_vec(5))**zoo(m)%nswitch / &
     &                   (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       tot_prey(m) =  pa_matrix(m,1)*prey_vec(1) +                       &
     &                pa_matrix(m,2)*prey_vec(2) +                       &
     &                pa_matrix(m,5)*prey_vec(5)
       ingest_matrix(m,1) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*          &
     &                      pa_matrix(m,1)*                              &
     &                      prey_vec(1)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,2) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*          &
     &                      pa_matrix(m,2)*                              &
     &                      prey_vec(2)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,5) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*          &
     &                      pa_matrix(m,5)*                              &
     &                      prey_vec(5)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))

       zoo(m)%jingest_n(i,j,k) = ingest_matrix(m,1)+ingest_matrix(m,2) + &
     &                           ingest_matrix(m,5)
       zoo(m)%jingest_p(i,j,k) = ingest_matrix(m,1)*prey_p2n_vec(1) +    &
     &                           ingest_matrix(m,2)*prey_p2n_vec(2) +    &
     &                           ingest_matrix(m,5)*prey_p2n_vec(5)
       zoo(m)%jingest_fe(i,j,k) = ingest_matrix(m,1)*prey_fe2n_vec(1) +  &
     &                            ingest_matrix(m,2)*prey_fe2n_vec(2)
       zoo(m)%jingest_sio2(i,j,k) = ingest_matrix(m,2)*prey_si2n_vec(2)

       !
       ! Large zooplankton (m = 3) consuming diazotrophs (1), large phytoplankton (2)
       ! and medium zooplankton (6)
       !

       m = 3
       zoo(m)%temp_lim(i,j,k) = exp(zoo(m)%ktemp*cobalt%f_temp(i,j,k))
       sw_fac_denom = (ipa_matrix(m,1)*prey_vec(1))**zoo(m)%nswitch +    &
     &                (ipa_matrix(m,2)*prey_vec(2))**zoo(m)%nswitch +    &
     &                (ipa_matrix(m,6)*prey_vec(6))**zoo(m)%nswitch
       pa_matrix(m,1) = ipa_matrix(m,1)*                                 &
     &                 ( (ipa_matrix(m,1)*prey_vec(1))**zoo(m)%nswitch / &
     &                   (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,2) = ipa_matrix(m,2)*                                 &
     &                 ( (ipa_matrix(m,2)*prey_vec(2))**zoo(m)%nswitch / &
     &                   (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       pa_matrix(m,6) = ipa_matrix(m,6)*                                 &
     &                 ( (ipa_matrix(m,6)*prey_vec(6))**zoo(m)%nswitch / &
     &                   (sw_fac_denom+epsln) )**(1.0/zoo(m)%mswitch)
       tot_prey(m) =  pa_matrix(m,1)*prey_vec(1) + pa_matrix(m,2)*       &
     &                prey_vec(2) +                                      &
     &                pa_matrix(m,6)*prey_vec(6)
       ingest_matrix(m,1) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*          &
     &                      pa_matrix(m,1)*                              &
     &                      prey_vec(1)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,2) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*          &
     &                      pa_matrix(m,2)*                              &
     &                      prey_vec(2)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       ingest_matrix(m,6) = zoo(m)%temp_lim(i,j,k)*zoo(m)%imax*          &
     &                      pa_matrix(m,6)*                              &
     &                      prey_vec(6)*zoo(m)%f_n(i,j,k)/(zoo(m)%ki+tot_prey(m))
       zoo(m)%jingest_n(i,j,k) = ingest_matrix(m,1)+ingest_matrix(m,2) + &
     &                           ingest_matrix(m,6)
       zoo(m)%jingest_p(i,j,k) = ingest_matrix(m,1)*prey_p2n_vec(1) +    &
     &                           ingest_matrix(m,2)*prey_p2n_vec(2) +    &
     &                           ingest_matrix(m,6)*prey_p2n_vec(6)
       zoo(m)%jingest_fe(i,j,k) = ingest_matrix(m,1)*prey_fe2n_vec(1) +  &
     &                            ingest_matrix(m,2)*prey_fe2n_vec(2)
       zoo(m)%jingest_sio2(i,j,k) = ingest_matrix(m,2)*prey_si2n_vec(2)

       cobalt%total_filter_feeding(i,j,k) = ingest_matrix(2,1) +         &
     &                                      ingest_matrix(2,2) +         &
     &    ingest_matrix(2,3) + ingest_matrix(3,1) + ingest_matrix(3,2) + & 
     &    ingest_matrix(3,3) !+ hp_ingest_vec(1) + hp_ingest_vec(2) + hp_ingest_vec(3)  ! charlie's notes

       !
       ! Calculate losses to zooplankton
       !

       DO nphyto=1,NUM_PHYTO
          phyto(nphyto)%jzloss_n(i,j,k) = 0.0
       ENDDO
       nphyto=overflow

       DO m=1,NUM_ZOO !{
          phyto(DIAZO)%jzloss_n(i,j,k) = phyto(DIAZO)%jzloss_n(i,j,k) + ingest_matrix(m,DIAZO)
          phyto(LARGE)%jzloss_n(i,j,k) = phyto(LARGE)%jzloss_n(i,j,k) + ingest_matrix(m,LARGE)
          phyto(SMALL)%jzloss_n(i,j,k) = phyto(SMALL)%jzloss_n(i,j,k) + ingest_matrix(m,SMALL)
       ENDDO !} m
       m=overflow

       DO nphyto=1,NUM_PHYTO !{
          phyto(nphyto)%jzloss_p(i,j,k) = phyto(nphyto)%jzloss_n(i,j,k)*prey_p2n_vec(nphyto)
          phyto(nphyto)%jzloss_fe(i,j,k) = phyto(nphyto)%jzloss_n(i,j,k)*prey_fe2n_vec(nphyto)
          phyto(nphyto)%jzloss_sio2(i,j,k) = phyto(nphyto)%jzloss_n(i,j,k)*prey_si2n_vec(nphyto)  
       ENDDO !} n
       nphyto=overflow

       !
       ! losses of bacteria to zooplankton 
       !

       bact(1)%jzloss_n(i,j,k) = 0.0
       DO m=1,NUM_ZOO !{
          bact(1)%jzloss_n(i,j,k) = bact(1)%jzloss_n(i,j,k) + ingest_matrix(m,4)
       ENDDO !} m
       m=overflow
       bact(1)%jzloss_p(i,j,k) = bact(1)%jzloss_n(i,j,k)*prey_p2n_vec(4)

       !
       ! losses of zooplankton to zooplankton
       !

       DO nzoo=1,NUM_ZOO !{
          zoo(nzoo)%jzloss_n(i,j,k) = 0.0

          DO m=1,NUM_ZOO !{
             zoo(nzoo)%jzloss_n(i,j,k) = zoo(nzoo)%jzloss_n(i,j,k) + ingest_matrix(m,NUM_PHYTO+1+nzoo)
          ENDDO !} m
          m=overflow

          zoo(nzoo)%jzloss_p(i,j,k) = zoo(nzoo)%jzloss_n(i,j,k)*prey_p2n_vec(NUM_PHYTO+1+nzoo)
       ENDDO !} n
       nzoo=overflow 

       !
       ! 3.1.2 Calculate ingestion by higher predators
       !

       ! The higher-predator ingestion calculations mirror those used for zooplankton
       !
       cobalt%hp_temp_lim(i,j,k) = exp(cobalt%ktemp_hp*cobalt%f_temp(i,j,k))
       sw_fac_denom = (hp_ipa_vec(6)*prey_vec(6))**cobalt%nswitch_hp +   &
     &                (hp_ipa_vec(7)*prey_vec(7))**cobalt%nswitch_hp
       hp_pa_vec(6) = hp_ipa_vec(6)*                                     &
     &                ( (hp_ipa_vec(6)*prey_vec(6))**cobalt%nswitch_hp / &
     &                  (sw_fac_denom+epsln) )**(1.0/cobalt%mswitch_hp)
       hp_pa_vec(7) = hp_ipa_vec(7)*                                     &
     &                ( (hp_ipa_vec(7)*prey_vec(7))**cobalt%nswitch_hp / &
     &                  (sw_fac_denom+epsln) )**(1.0/cobalt%mswitch_hp)
       tot_prey_hp = hp_pa_vec(6)*prey_vec(6) + hp_pa_vec(7)*prey_vec(7)
       hp_ingest_vec(6) = cobalt%hp_temp_lim(i,j,k)*cobalt%imax_hp*      &
     &                    hp_pa_vec(6)*                                  &
     &                    prey_vec(6)*tot_prey_hp**(cobalt%coef_hp-1)/   &
     &                    (cobalt%ki_hp+tot_prey_hp)
       hp_ingest_vec(7) = cobalt%hp_temp_lim(i,j,k)*cobalt%imax_hp*      &
     &                    hp_pa_vec(7)*                                  &
     &                    prey_vec(7)*tot_prey_hp**(cobalt%coef_hp-1)/   &
     &                    (cobalt%ki_hp+tot_prey_hp)
       cobalt%hp_jingest_n(i,j,k) = hp_ingest_vec(6) + hp_ingest_vec(7)
       cobalt%hp_jingest_p(i,j,k) = hp_ingest_vec(6)*prey_p2n_vec(6) +   &
     &                              hp_ingest_vec(7)*prey_p2n_vec(7)

       !
       ! losses of zooplankton to higher predators
       !
       DO nzoo=1,NUM_ZOO !{
         zoo(nzoo)%jhploss_n(i,j,k) = hp_ingest_vec(NUM_PHYTO+1+nzoo)
         zoo(nzoo)%jhploss_p(i,j,k) = zoo(nzoo)%jhploss_n(i,j,k)*prey_p2n_vec(NUM_PHYTO+1+nzoo)
       ENDDO !} n
       nzoo=overflow

      ! RD debug, just a diag 
      DO nzoo=1,NUM_ZOO
        DO m=1,NUM_PREY
           ingest_matrix_glo(nzoo,m) = MAX(ingest_matrix_glo(nzoo,m),ingest_matrix(nzoo,m))
        ENDDO
      ENDDO
      ! end

      ENDDO
    ENDDO
  ENDDO
  i=overflow ; j=overflow ; k=overflow

!  WRITE(stdout,*) 'Max of ingest matrix is ', MAXVAL(ingest_matrix_glo), ' at ', MAXLOC(ingest_matrix_glo)
!  WRITE(stdout,*) 'Min of ingest matrix is ', MINVAL(ingest_matrix_glo), ' at ', MINLOC(ingest_matrix_glo)

!
! 3.2: Plankton foodweb dynamics: Other mortality and loss terms
!
!  
! 3.2.1 Calculate losses of phytoplankton to aggregation 
!

  DO nphyto=1,NUM_PHYTO 
       phyto(nphyto)%jaggloss_n(:,:,:) = phyto(nphyto)%agg*phyto(nphyto)%f_n(:,:,:)**2.0 
       phyto(nphyto)%jaggloss_p(:,:,:) = phyto(nphyto)%jaggloss_n(:,:,:)*phyto(nphyto)%q_p_2_n(:,:,:)
       phyto(nphyto)%jaggloss_fe(:,:,:) = phyto(nphyto)%jaggloss_n(:,:,:)*phyto(nphyto)%q_fe_2_n(:,:,:)
       phyto(nphyto)%jaggloss_sio2(:,:,:) = phyto(nphyto)%jaggloss_n(:,:,:)*phyto(nphyto)%q_si_2_n(:,:,:)
  ENDDO 
  nphyto=overflow

!
! 3.2.2 Calculate phytoplankton and bacterial losses to viruses
!

  DO nphyto=1,NUM_PHYTO 
     phyto(nphyto)%jvirloss_n(:,:,:) = bact(1)%temp_lim(:,:,:)*phyto(nphyto)%vir*phyto(nphyto)%f_n(:,:,:)**2.0 
     phyto(nphyto)%jvirloss_p(:,:,:) = phyto(nphyto)%jvirloss_n(:,:,:)*phyto(nphyto)%q_p_2_n(:,:,:)
     phyto(nphyto)%jvirloss_fe(:,:,:) = phyto(nphyto)%jvirloss_n(:,:,:)*phyto(nphyto)%q_fe_2_n(:,:,:)
     phyto(nphyto)%jvirloss_sio2(:,:,:) = phyto(nphyto)%jvirloss_n(:,:,:)*phyto(nphyto)%q_si_2_n(:,:,:)
  ENDDO 
  nphyto=overflow

  bact(1)%jvirloss_n(:,:,:) = bact(1)%temp_lim(:,:,:)*bact(1)%vir*bact(1)%f_n(:,:,:)**2.0
  bact(1)%jvirloss_p(:,:,:) = bact(1)%jvirloss_n(:,:,:)*bact(1)%q_p_2_n

!
! 3.2.3 Calculate losses to exudation
!

  nphyto = DIAZO
  phyto(nphyto)%jexuloss_n(:,:,:) = phyto(nphyto)%exu*                   &
  &                            max(phyto(nphyto)%juptake_no3(:,:,:)+     &
  &                            phyto(nphyto)%juptake_nh4(:,:,:)+phyto(nphyto)%juptake_n2(:,:,:),0.0)
  phyto(nphyto)%jexuloss_p(:,:,:) = phyto(nphyto)%exu*max(phyto(nphyto)%juptake_po4(:,:,:),0.0)
  phyto(nphyto)%jexuloss_fe(:,:,:) = phyto(nphyto)%exu*max(phyto(nphyto)%juptake_fe(:,:,:),0.0)
  DO nphyto=2,NUM_PHYTO 
     phyto(nphyto)%jexuloss_n(:,:,:) = phyto(nphyto)%exu*                &
  &  max(phyto(nphyto)%juptake_no3(:,:,:)+phyto(nphyto)%juptake_nh4(:,:,:),0.0)
     phyto(nphyto)%jexuloss_p(:,:,:) = phyto(nphyto)%exu*max(phyto(nphyto)%juptake_po4(:,:,:),0.0)
     phyto(nphyto)%jexuloss_fe(:,:,:) = phyto(nphyto)%exu*max(phyto(nphyto)%juptake_fe(:,:,:),0.0)
  ENDDO
  nphyto=overflow

!
! 3.3: Plankton foodweb dynamics: Production calculations
!

  DO k=1,Ubk
    DO j=Jstr,Jend
      DO i=Istr,Iend

!
! 3.3.1: Calculate the production of detritus and dissolved organic material
!

       ! initialize some cumulative COBALT-wide production diagnostics
       cobalt%jprod_ndet(i,j,k)  = 0.0
       cobalt%jprod_pdet(i,j,k)  = 0.0
       cobalt%jprod_sldon(i,j,k) = 0.0
       cobalt%jprod_ldon(i,j,k)  = 0.0
       cobalt%jprod_srdon(i,j,k) = 0.0
       cobalt%jprod_sldop(i,j,k) = 0.0
       cobalt%jprod_ldop(i,j,k)  = 0.0
       cobalt%jprod_srdop(i,j,k) = 0.0
       cobalt%jprod_fedet(i,j,k) = 0.0
       cobalt%jprod_fed(i,j,k)   = 0.0
       cobalt%jprod_sidet(i,j,k) = 0.0
       cobalt%jprod_sio4(i,j,k)  = 0.0
       cobalt%jprod_po4(i,j,k)   = 0.0
       cobalt%jprod_nh4(i,j,k)   = 0.0

!
! Production of detritus and dissolved organic material from zooplankton egestion 
!   

       DO m=1,NUM_ZOO
           zoo(m)%jprod_ndet(i,j,k)  = zoo(m)%phi_det*zoo(m)%jingest_n(i,j,k)
           zoo(m)%jprod_pdet(i,j,k)  = zoo(m)%phi_det*zoo(m)%jingest_p(i,j,k)
           zoo(m)%jprod_sldon(i,j,k) = zoo(m)%phi_sldon*zoo(m)%jingest_n(i,j,k)
           zoo(m)%jprod_ldon(i,j,k)  = zoo(m)%phi_ldon*zoo(m)%jingest_n(i,j,k)
           zoo(m)%jprod_srdon(i,j,k) = zoo(m)%phi_srdon*zoo(m)%jingest_n(i,j,k)
           zoo(m)%jprod_sldop(i,j,k) = zoo(m)%phi_sldop*zoo(m)%jingest_p(i,j,k)
           zoo(m)%jprod_ldop(i,j,k)  = zoo(m)%phi_ldop*zoo(m)%jingest_p(i,j,k)
           zoo(m)%jprod_srdop(i,j,k) = zoo(m)%phi_srdop*zoo(m)%jingest_p(i,j,k)
           zoo(m)%jprod_fedet(i,j,k) = zoo(m)%phi_det*zoo(m)%jingest_fe(i,j,k)
           zoo(m)%jprod_sidet(i,j,k) = zoo(m)%phi_det*zoo(m)%jingest_sio2(i,j,k)


           ! augment cumulative production with zooplankton terms
           cobalt%jprod_ndet(i,j,k)  = cobalt%jprod_ndet(i,j,k)  + zoo(m)%jprod_ndet(i,j,k)
           cobalt%jprod_pdet(i,j,k)  = cobalt%jprod_pdet(i,j,k)  + zoo(m)%jprod_pdet(i,j,k)
           cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) + zoo(m)%jprod_sldon(i,j,k)
           cobalt%jprod_ldon(i,j,k)  = cobalt%jprod_ldon(i,j,k)  + zoo(m)%jprod_ldon(i,j,k)
           cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) + zoo(m)%jprod_srdon(i,j,k)
           cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) + zoo(m)%jprod_sldop(i,j,k)
           cobalt%jprod_ldop(i,j,k)  = cobalt%jprod_ldop(i,j,k)  + zoo(m)%jprod_ldop(i,j,k)
           cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) + zoo(m)%jprod_srdop(i,j,k)
           cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) + zoo(m)%jprod_fedet(i,j,k)
           cobalt%jprod_sidet(i,j,k) = cobalt%jprod_sidet(i,j,k) + zoo(m)%jprod_sidet(i,j,k)
       ENDDO !} m
       m=overflow

       !
       ! Production of detritus and dissolved organic material from higher predator egestion 
       ! (did not track individual terms, just add to cumulative total)
       !

       cobalt%jprod_ndet(i,j,k)  = cobalt%jprod_ndet(i,j,k)  + cobalt%hp_phi_det*cobalt%hp_jingest_n(i,j,k)
       cobalt%jprod_pdet(i,j,k)  = cobalt%jprod_pdet(i,j,k)  + cobalt%hp_phi_det*cobalt%hp_jingest_p(i,j,k)
       cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) + cobalt%hp_phi_sldon*cobalt%hp_jingest_n(i,j,k)
       cobalt%jprod_ldon(i,j,k)  = cobalt%jprod_ldon(i,j,k)  + cobalt%hp_phi_ldon*cobalt%hp_jingest_n(i,j,k)
       cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) + cobalt%hp_phi_srdon*cobalt%hp_jingest_n(i,j,k)
       cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) + cobalt%hp_phi_sldop*cobalt%hp_jingest_p(i,j,k)
       cobalt%jprod_ldop(i,j,k)  = cobalt%jprod_ldop(i,j,k)  + cobalt%hp_phi_ldop*cobalt%hp_jingest_p(i,j,k)
       cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) + cobalt%hp_phi_srdop*cobalt%hp_jingest_p(i,j,k)
       cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) + cobalt%hp_phi_det*cobalt%hp_jingest_fe(i,j,k)
       cobalt%jprod_sidet(i,j,k) = cobalt%jprod_sidet(i,j,k) + cobalt%hp_phi_det*cobalt%hp_jingest_sio2(i,j,k)
       
       !
       ! Sources from phytoplankton aggregation
       !

       do m=1,NUM_PHYTO
           cobalt%jprod_ndet(i,j,k)  = cobalt%jprod_ndet(i,j,k)  + phyto(m)%jaggloss_n(i,j,k)
           cobalt%jprod_pdet(i,j,k)  = cobalt%jprod_pdet(i,j,k)  + phyto(m)%jaggloss_p(i,j,k)
           cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) + phyto(m)%jaggloss_fe(i,j,k)
           cobalt%jprod_sidet(i,j,k) = cobalt%jprod_sidet(i,j,k) + phyto(m)%jaggloss_sio2(i,j,k)
       enddo !} m
       m=overflow

       !
       ! Sources from viral lysis of phytoplankton (0 in default formulation) and exudation
       !

       DO m=1,NUM_PHYTO
           cobalt%jprod_ldon(i,j,k) = cobalt%jprod_ldon(i,j,k) +         &
        &             cobalt%lysis_phi_ldon*phyto(m)%jvirloss_n(i,j,k) + &
        &             phyto(m)%jexuloss_n(i,j,k) 
           cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) +       &
        &                              cobalt%lysis_phi_sldon*phyto(m)%jvirloss_n(i,j,k)
           cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) +       &
        &                              cobalt%lysis_phi_srdon*phyto(m)%jvirloss_n(i,j,k)
           cobalt%jprod_ldop(i,j,k) = cobalt%jprod_ldop(i,j,k) +         &
        &             cobalt%lysis_phi_ldop*phyto(m)%jvirloss_p(i,j,k) + &
        &             phyto(m)%jexuloss_p(i,j,k)
           cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) +       &
        &                              cobalt%lysis_phi_sldop*phyto(m)%jvirloss_p(i,j,k)
           cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) +       &
        &                              cobalt%lysis_phi_srdop*phyto(m)%jvirloss_p(i,j,k)
           cobalt%jprod_fed(i,j,k)   = cobalt%jprod_fed(i,j,k)   +       &
        &                              phyto(m)%jvirloss_fe(i,j,k) +     &
        &                              phyto(m)%jexuloss_fe(i,j,k) 
           cobalt%jprod_sio4(i,j,k) = cobalt%jprod_sio4(i,j,k) + phyto(m)%jvirloss_sio2(i,j,k)
       ENDDO !} m
       m=overflow

!
! Sources of dissolved organic material from viral lysis due to bacteria 
!

       cobalt%jprod_ldon(i,j,k)  = cobalt%jprod_ldon(i,j,k)  + cobalt%lysis_phi_ldon*bact(1)%jvirloss_n(i,j,k)
       cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) + cobalt%lysis_phi_sldon*bact(1)%jvirloss_n(i,j,k)
       cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) + cobalt%lysis_phi_srdon*bact(1)%jvirloss_n(i,j,k)
       cobalt%jprod_ldop(i,j,k)  = cobalt%jprod_ldop(i,j,k)  + cobalt%lysis_phi_ldop*bact(1)%jvirloss_p(i,j,k)
       cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) + cobalt%lysis_phi_sldop*bact(1)%jvirloss_p(i,j,k)
       cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) + cobalt%lysis_phi_srdop*bact(1)%jvirloss_p(i,j,k)

!
! Sources of dissolved organic material from bacterial mortality (metabolic costs higher than food uptake).
! These conditions are assumed to lead to a lysis-like redistribution of bacteria organic matter.
!

       cobalt%jprod_ldon(i,j,k) = cobalt%jprod_ldon(i,j,k) -             &
     &                            cobalt%lysis_phi_ldon*                 &
     &                            min(bact(1)%jprod_n(i,j,k),0.0)
       cobalt%jprod_sldon(i,j,k) = cobalt%jprod_sldon(i,j,k) -           &
     &                             cobalt%lysis_phi_sldon*               &
     &                             min(bact(1)%jprod_n(i,j,k),0.0)
       cobalt%jprod_srdon(i,j,k) = cobalt%jprod_srdon(i,j,k) -           &
     &                             cobalt%lysis_phi_srdon*               &
     &                             min(bact(1)%jprod_n(i,j,k),0.0)
       cobalt%jprod_ldop(i,j,k) = cobalt%jprod_ldop(i,j,k) -             &
     &                            cobalt%lysis_phi_ldop*                 &
     &                            min(bact(1)%jprod_n(i,j,k)*bact(1)%q_p_2_n,0.0)
       cobalt%jprod_sldop(i,j,k) = cobalt%jprod_sldop(i,j,k) -           &
     &                             cobalt%lysis_phi_sldop*               &
     &                             min(bact(1)%jprod_n(i,j,k)*bact(1)%q_p_2_n,0.0)
       cobalt%jprod_srdop(i,j,k) = cobalt%jprod_srdop(i,j,k) -           &
     &                             cobalt%lysis_phi_srdop*               &
     &                             min(bact(1)%jprod_n(i,j,k)*bact(1)%q_p_2_n,0.0)
!
! 3.3.2: Calculate the remineralization of organic material by free-living bacteria
!

       bact(1)%jprod_nh4(i,j,k) = bact(1)%juptake_ldon(i,j,k) - max(bact(1)%jprod_n(i,j,k),0.0)
       bact(1)%jprod_po4(i,j,k) = bact(1)%juptake_ldop(i,j,k) -          &
     &                            max(bact(1)%jprod_n(i,j,k)*bact(1)%q_p_2_n,0.0)
       cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) + bact(1)%jprod_nh4(i,j,k)
       cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) + bact(1)%jprod_po4(i,j,k)

!
! 3.3.3: Zooplankton production and excretion calculations
!

       DO m=1,NUM_ZOO

          ingest_p2n = zoo(m)%jingest_p(i,j,k)/(zoo(m)%jingest_n(i,j,k)+epsln)

          IF (ingest_p2n .lt. zoo(m)%q_p_2_n) THEN
             zoo(m)%jprod_n(i,j,k) = zoo(m)%gge_max*zoo(m)%jingest_p(i,j,k)*(1.0/zoo(m)%q_p_2_n)
          ELSE
             zoo(m)%jprod_n(i,j,k) = zoo(m)%gge_max*zoo(m)%jingest_n(i,j,k)
          ENDIF

          ! adjust production terms for basal respiration costs
            zoo(m)%jprod_n(i,j,k) = zoo(m)%jprod_n(i,j,k) -              &
          &                         zoo(m)%f_n(i,j,k)/(refuge_conc +     &
          &                         zoo(m)%f_n(i,j,k))*                  &
          &                         zoo(m)%temp_lim(i,j,k)*zoo(m)%bresp*zoo(m)%f_n(i,j,k)
          !
          ! Ingested material that does not go to zooplankton production, detrital production
          ! or production of dissolved organic material is excreted as nh4 or po4.  If production
          ! is negative, zooplankton are lost to large detritus 
          !
          IF (zoo(m)%jprod_n(i,j,k) .gt. 0.0) THEN 
             zoo(m)%jprod_nh4(i,j,k) =  zoo(m)%jingest_n(i,j,k)   -      &
           &                            zoo(m)%jprod_ndet(i,j,k)  -      &
           &                            zoo(m)%jprod_n(i,j,k)     -      &
           &                            zoo(m)%jprod_ldon(i,j,k)  -      &
           &                            zoo(m)%jprod_sldon(i,j,k) -      &
           &                            zoo(m)%jprod_srdon(i,j,k)
             zoo(m)%jprod_po4(i,j,k) =  zoo(m)%jingest_p(i,j,k)   -      &
           &                            zoo(m)%jprod_pdet(i,j,k)  -      & 
           &                            zoo(m)%jprod_n(i,j,k)     *      &
           &                            zoo(m)%q_p_2_n            -      &
           &                            zoo(m)%jprod_ldop(i,j,k)  -      &
           &                            zoo(m)%jprod_sldop(i,j,k) -      &
           &                            zoo(m)%jprod_srdop(i,j,k)
          ELSE
             ! None of the ingestion material goes to zooplankton production
             zoo(m)%jprod_nh4(i,j,k) =  zoo(m)%jingest_n(i,j,k)   -      &
           &                            zoo(m)%jprod_ndet(i,j,k)  -      &
           &                            zoo(m)%jprod_ldon(i,j,k)  -      &
           &                            zoo(m)%jprod_sldon(i,j,k) -      &
           &                            zoo(m)%jprod_srdon(i,j,k)
             zoo(m)%jprod_po4(i,j,k) =  zoo(m)%jingest_p(i,j,k)   -      &
           &                            zoo(m)%jprod_pdet(i,j,k)  -      &
           &                            zoo(m)%jprod_ldop(i,j,k)  -      &
           &                            zoo(m)%jprod_sldop(i,j,k) -      &
           &                            zoo(m)%jprod_srdop(i,j,k)

             ! The negative production (i.e., mortality) is lost to large detritus. Update values
             ! for zooplankton and for total.

             zoo(m)%jprod_ndet(i,j,k) = zoo(m)%jprod_ndet(i,j,k) - zoo(m)%jprod_n(i,j,k)
             zoo(m)%jprod_pdet(i,j,k) = zoo(m)%jprod_pdet(i,j,k) - zoo(m)%jprod_n(i,j,k)*zoo(m)%q_p_2_n
             cobalt%jprod_ndet(i,j,k) = cobalt%jprod_ndet(i,j,k) - zoo(m)%jprod_n(i,j,k)
             cobalt%jprod_pdet(i,j,k) = cobalt%jprod_pdet(i,j,k) - zoo(m)%jprod_n(i,j,k)*zoo(m)%q_p_2_n 
          ENDIF

          ! cumulative production of inorganic nutrients 
          cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) + zoo(m)%jprod_nh4(i,j,k)
          cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) + zoo(m)%jprod_po4(i,j,k)

          !
          ! Any ingested iron that is not allocated to detritus is routed back to the
          ! dissolved pool.       
          !
          zoo(m)%jprod_fed(i,j,k) = (1.0 - zoo(m)%phi_det)*zoo(m)%jingest_fe(i,j,k)
          cobalt%jprod_fed(i,j,k) = cobalt%jprod_fed(i,j,k) + zoo(m)%jprod_fed(i,j,k)
          !
          ! Any ingested opal that is not allocated to detritus is assumed to undergo
          ! rapid dissolution to dissolved silica
          !
          zoo(m)%jprod_sio4(i,j,k) = (1.0 - zoo(m)%phi_det)*zoo(m)%jingest_sio2(i,j,k)
          cobalt%jprod_sio4(i,j,k) = cobalt%jprod_sio4(i,j,k) + zoo(m)%jprod_sio4(i,j,k)
 
       ENDDO !} m
       m=overflow

       !
       ! Excretion by higher predators
       !
       cobalt%jprod_fed(i,j,k)  = cobalt%jprod_fed(i,j,k) +              &
     &                            (1.0-cobalt%hp_phi_det)*cobalt%hp_jingest_fe(i,j,k)
       cobalt%jprod_sio4(i,j,k) = cobalt%jprod_sio4(i,j,k) +             &
     &                            (1.0-cobalt%hp_phi_det)*cobalt%hp_jingest_sio2(i,j,k)
       cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) +               &
     &                           cobalt%hp_phi_nh4*cobalt%hp_jingest_n(i,j,k)
       cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) +               &
     &                           cobalt%hp_phi_po4*cobalt%hp_jingest_p(i,j,k)


      ENDDO
    ENDDO
  ENDDO
  i=overflow ; j=overflow ; k=overflow


   ! RD dev notes :
   ! cobalt%zt is the depth at rho-point, here computed by the difference
   ! between the position of the free surface and the rho-point. This quantity
   ! must be positive because it is used then to compute a pressure.
   DO k=1,UBk
     DO j=Jstr,Jend
       DO i=Istr,Iend
       
          cobalt%zt(i,j,k)=( z_w(i,j,N(ng)) - z_r(i,j,k) )

       ENDDO
     ENDDO
   ENDDO
   i=overflow ; j=overflow ; k=overflow

!
!------------------------------------------------------------------------------------
! 4: Production of calcium carbonate (Calcite and Aragonite) and lithogenic
! material
!------------------------------------------------------------------------------------
!

  DO k=1,UBk
    DO j=Jstr,Jend
      DO i=Istr,Iend

!
! 4.1: Calculate aragonite and calcite saturation states
!

       TK = cobalt%f_temp(i,j,k) + 273.15
       PRESS = 0.1016 * cobalt%zt(i,j,k) + 1.013
       PKSPA = 171.945 + 0.077993 * TK - 2903.293 / TK -                 &
     &         71.595 * log10(TK) - (-0.068393 + 1.7276e-3 *             &
     &         TK + 88.135 / TK)*sqrt(max(epsln,cobalt%f_salt(i,j,k))) + &
     &         0.10018 * max(epsln, cobalt%f_salt(i,j,k)) -              &
     &         5.9415e-3 * max(epsln, cobalt%f_salt(i,j,k))**(1.5) -     &
     &         0.02 - (48.76 - 2.8 - 0.5304 * cobalt%f_temp(i,j,k)) *    &
     &         (PRESS - 1.013) / (191.46 * TK) + (1e-3 * (11.76-0.3692 * &
     &         cobalt%f_temp(i,j,k))) * (PRESS - 1.013) *                &
     &         (PRESS - 1.013) / (382.92 * TK)

       cobalt%co3_sol_arag(i,j,k) = 10**(-PKSPA) / (2.937d-4 * max(5.0, cobalt%f_salt(i,j,k)))
       cobalt%omega_arag(i,j,k) = cobalt%f_co3_ion(i,j,k) / cobalt%co3_sol_arag(i,j,k)

       PKSPC = 171.9065 + 0.077993 * TK - 2839.319 / TK - 71.595 *       &
     &         log10(TK) - (-0.77712 + 2.8426e-3 *                       &
     &         TK + 178.34 / TK)*sqrt(max(epsln,cobalt%f_salt(i,j,k))) + &
     &         0.07711 * max(epsln, cobalt%f_salt(i,j,k)) -              &
     &         4.1249e-3 * max(epsln, cobalt%f_salt(i,j,k))**(1.5) -     &
     &         0.02 - (48.76 - 0.5304 * cobalt%f_temp(i,j,k)) *          &
     &         (PRESS - 1.013) / (191.46 * TK) + (1e-3 * (11.76 -        &
     &         0.3692 * cobalt%f_temp(i,j,k))) * (PRESS - 1.013) *       &
     &         (PRESS - 1.013) / (382.92 * TK)

       cobalt%co3_sol_calc(i,j,k) = 10**(-PKSPC) / (2.937d-4 * max(5.0, cobalt%f_salt(i,j,k)))
       cobalt%omega_calc(i,j,k) = cobalt%f_co3_ion(i,j,k) / cobalt%co3_sol_calc(i,j,k)


#ifdef DIAGNOSTICS_BIO
      DiaBio3d(i,j,k,ico3_sol_arag)     = DiaBio3d(i,j,k,ico3_sol_arag)     + cobalt%co3_sol_arag(i,j,k)
      DiaBio3d(i,j,k,ico3_sol_calc)     = DiaBio3d(i,j,k,ico3_sol_calc)     + cobalt%co3_sol_calc(i,j,k)
      DiaBio3d(i,j,k,iomega_cadet_arag) = DiaBio3d(i,j,k,iomega_cadet_arag) + cobalt%omega_arag(i,j,k)
      DiaBio3d(i,j,k,iomega_cadet_calc) = DiaBio3d(i,j,k,iomega_cadet_calc) + cobalt%omega_calc(i,j,k)
#endif

      ENDDO
    ENDDO
  ENDDO
  i=overflow ; j=overflow ; k=overflow

!
! 4.2: Calculate the production rate of aragonite and calcite detritus 
!

  DO k=1,UBk
    DO j=Jstr,Jend
      DO i=Istr,Iend
        cobalt%jprod_cadet_arag(i,j,k) = ( zoo(2)%jzloss_n(i,j,k) +      &
      &                                    zoo(3)%jzloss_n(i,j,k) +      &
      &                                    zoo(2)%jhploss_n(i,j,k) +     &
      &                                    zoo(3)%jhploss_n(i,j,k) ) *   &
      & cobalt%ca_2_n_arag * min( cobalt%caco3_sat_max,                  &
      & max(0.0d0,cobalt%omega_arag(i,j,k) - 1.0d0) ) +                  &
      & epsln

        cobalt%jprod_cadet_calc(i,j,k) = (zoo(1)%jzloss_n(i,j,k) +       &
      &              phyto(SMALL)%jaggloss_n(i,j,k))*cobalt%ca_2_n_calc* &
      &              min(cobalt%caco3_sat_max, max(0.0d0, cobalt%omega_calc(i,j,k) - 1.0d0)) + epsln
      ENDDO
    ENDDO
  ENDDO
  i=overflow ; j=overflow ; k=overflow

!
! 4.3: Lithogenic detritus production (repackaged from f_lith during filter feeding)
!

  DO k=1,UBk
    DO j=Jstr,Jend
      DO i=Istr,Iend

       cobalt%jprod_lithdet(i,j,k)=( cobalt%total_filter_feeding(i,j,k)/ &
     &                             ( phyto(LARGE)%f_n(i,j,k) +           &
     &                               phyto(DIAZO)%f_n(i,j,k) + epsln ) * &  
     &                               cobalt%phi_lith + cobalt%k_lith ) * cobalt%f_lith(i,j,k)
      ENDDO
    ENDDO
  ENDDO
  i=overflow ; j=overflow ; k=overflow

!
!---------------------------------------------------------------------------------------------------------
! 5: Detrital dissolution and remineralization calculation
!---------------------------------------------------------------------------------------------------------
!
!
! 5.1: Dissolution of aragonite, calcite and opal detrital particles
!

   cobalt%jdiss_cadet_arag(:,:,:) = cobalt%gamma_cadet_arag * &
 &       max(0.0, 1.0 - cobalt%omega_arag(:,:,:)) * cobalt%f_cadet_arag(:,:,:)
   cobalt%jdiss_cadet_calc(:,:,:) = cobalt%gamma_cadet_calc * &
 &       max(0.0, 1.0 - cobalt%omega_calc(:,:,:)) * cobalt%f_cadet_calc(:,:,:)
   cobalt%jdiss_sidet(:,:,:) = cobalt%gamma_sidet * cobalt%f_sidet(:,:,:)
   cobalt%jprod_sio4(:,:,:) = cobalt%jprod_sio4(:,:,:) + cobalt%jdiss_sidet(:,:,:)

!
! 5.2: Remineralization of nitrogen, phosphorous and iron detritus accounting for oxygen 
!      and mineral protection 
!

  DO k=1,UBk
    DO j=Jstr,Jend
      DO i=Istr,Iend

       cobalt%jno3denit_wc(i,j,k) = 0.0
       !
       !   Under oxic conditions
       !
       IF (cobalt%f_o2(i,j,k) .gt. cobalt%o2_min) THEN  !{
          cobalt%jremin_ndet(i,j,k) = cobalt%gamma_ndet *                &
        &                             cobalt%f_o2(i,j,k) /               & 
        &                           ( cobalt%k_o2+cobalt%f_o2(i,j,k) )*  &
        &                             max( 0.0, cobalt%f_ndet(i,j,k) -   &
        &                             cobalt%rpcaco3*                    &
        &                            (cobalt%f_cadet_arag(i,j,k) +       &
        &                             cobalt%f_cadet_calc(i,j,k)) -      & 
        & cobalt%rplith*cobalt%f_lithdet(i,j,k) - cobalt%rpsio2*cobalt%f_sidet(i,j,k) )
       !
       ! Under sub-oxic conditions
       !
       ELSE !}{
          cobalt%jremin_ndet(i,j,k) = cobalt%gamma_ndet*cobalt%o2_min /  &
        &      (cobalt%k_o2 + cobalt%o2_min)*                            &
        &      cobalt%f_no3(i,j,k) / (phyto(SMALL)%k_no3 +               &
        &      cobalt%f_no3(i,j,k))*                                     &
        &      max(0.0, cobalt%f_ndet(i,j,k) -                           &
        &      cobalt%rpcaco3*(cobalt%f_cadet_arag(i,j,k) +              &
        &      cobalt%f_cadet_calc(i,j,k)) -                             &
        &      cobalt%rplith*cobalt%f_lithdet(i,j,k) - cobalt%rpsio2*cobalt%f_sidet(i,j,k) )
!RD dev notes : was + cobalt%rpsio2*cobalt%f_sidet(i,j,k) in original code but -
!in charlie's suppl material (eq (58))


          cobalt%jno3denit_wc(i,j,k) = cobalt%jremin_ndet(i,j,k) * cobalt%n_2_n_denit
#ifdef COBALT_CONSERVATION_TEST
          cobalt%jno3denit_wc(i,j,k) = 0.
#endif
       ENDIF !}

       ! RD testing a simpler version of remineralization
       !cobalt%jremin_ndet(i,j,k) = cobalt%gamma_ndet*cobalt%f_ndet(i,j,k)

       !
       ! P and Fe assumed to be protected similarly to N
       !
       cobalt%jremin_pdet(i,j,k) = cobalt%jremin_ndet(i,j,k)/            &
     &                            (cobalt%f_ndet(i,j,k) + epsln)*        &
     &                             cobalt%f_pdet(i,j,k)
       cobalt%jremin_fedet(i,j,k) = cobalt%jremin_ndet(i,j,k) /          &
     &                             (cobalt%f_ndet(i,j,k) + epsln) *      &
     &                              cobalt%remin_eff_fedet*cobalt%f_fedet(i,j,k)

      ENDDO
    ENDDO
  ENDDO
  i=overflow ; j=overflow ; k=overflow

!
!--------------------------------------------------------------------------------------------
! 6: Miscellaneous sources and sinks: Nitrification, Iron Scavenging, Coastal Iron inputs
!--------------------------------------------------------------------------------------------
!

  DO k=1,UBk
    DO j=Jstr,Jend
      DO i=Istr,Iend

       !
       !  Nitrification
       !

       cobalt%jnitrif(i,j,k) = cobalt%gamma_nitrif*cobalt%expkT(i,j,k) * &
     &                         cobalt%f_nh4(i,j,k) *                     &
     &                         phyto(SMALL)%nh4lim(i,j,k) *              &
     &                         (1.0 - cobalt%f_irr_mem(i,j,k) /          &
     &                         (cobalt%irr_inhibit + cobalt%f_irr_mem(i,j,k)))

       !
       ! Solve for free iron
       !

       cobalt%kfe_eq_lig(i,j,k) = min(cobalt%kfe_eq_lig_ll,              &
     &                                10**( log10(cobalt%kfe_eq_lig_hl)+ &
     &      max(0.0,cobalt%gamma_fescav*log10(cobalt%io_fescav/cobalt%irr_inst(i,j,k))) ) ) 

       feprime = 1.0 + cobalt%kfe_eq_lig(i,j,k) * (cobalt%felig_bkg +    &
     &                                             cobalt%felig_2_don *  &
     &           (cobalt%f_sldon(i,j,k) + cobalt%f_srdon(i,j,k)) - cobalt%f_fed(i,j,k))
       feprime = (-feprime + (feprime * feprime + 4.0 *                  &
     &             cobalt%kfe_eq_lig(i,j,k) *                            &
     &             cobalt%f_fed(i,j,k))**(0.5)) / (2.0 * cobalt%kfe_eq_lig(i,j,k))

       !
       ! Iron adsorption to detrital particles
       !

       cobalt%jfe_ads(i,j,k) = min(r_dt,cobalt%alpha_fescav*feprime)
       IF (cobalt%f_fed(i,j,k).gt.1.0e-9) THEN !{
          cobalt%jfe_ads(i,j,k) = min(r_dt,5.0*cobalt%alpha_fescav*cobalt%f_fed(i,j,k))
       ENDIF !}
       !
       ! Coastal iron inputs (proxy for sediment inputs for areas with poorly resolved shelves)
       !
!       cobalt%jfe_coast(i,j,k) = cobalt%fe_coast * mask_coast(i,j) *     &
!     &                           grid_tmask(i,j,k) / sqrt(grid_dat(i,j))
      ! RD dev notes : we are supposed to resolve the shelves in ROMS now
      cobalt%jfe_coast(i,j,k) = 0.

      ENDDO
    ENDDO
  ENDDO
  i=overflow ; j=overflow ; k=overflow


!
!-------------------------------------------------------------------------------------------------
! Sinking of detritus (adapted from nemuro)
!-------------------------------------------------------------------------------------------------
!

  DO j=Jstr,Jend

    ! compute inverse thickness,...
    DO k=1,N(ng)
      DO i=Istr,Iend
        Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
      END DO
    END DO
    DO k=1,N(ng)-1
      DO i=Istr,Iend
        Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
      END DO
    END DO
    DO k=2,N(ng)-1
      DO i=Istr,Iend
        Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
      END DO
    END DO

    DO isink=1,Nsink
       ibio=idsink(isink)

       ! reset losses to sediment
       sed_losses(:,j) = 0.

       DO k=1,N(ng)
         DO i=Istr,Iend
           !qc(i,k)=max(t(i,j,k,nstp,ibio),0.0d0) ! we don't want negative values
           qc(i,k)=t(i,j,k,nstp,ibio) ! we want negative values
           !qc(i,k)=max(t(i,j,k,nnew,ibio) * Hz_inv(i,k) ,0.0d0) ! we don't want negative values
! RD #TEST
!            qc(i,k)=0.0d0 ! set zeros everywhere and hope for the best
         END DO
       END DO
!
       DO k=N(ng)-1,1,-1
         DO i=Istr,Iend
           FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
         END DO
       END DO
       DO k=2,N(ng)-1
         DO i=Istr,Iend
           dltR=Hz(i,j,k)*FC(i,k)
           dltL=Hz(i,j,k)*FC(i,k-1)
           cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
           cffR=cff*FC(i,k)
           cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
           IF ((dltR*dltL).le.0.0_r8) THEN
             dltR=0.0_r8
             dltL=0.0_r8
           ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
             dltR=cffL
           ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
             dltL=cffR
           END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
           cff=(dltR-dltL)*Hz_inv3(i,k)
           dltR=dltR-cff*Hz(i,j,k+1)
           dltL=dltL+cff*Hz(i,j,k-1)
           bR(i,k)=qc(i,k)+dltR
           bL(i,k)=qc(i,k)-dltL
           WR(i,k)=(2.0_r8*dltR-dltL)**2
           WL(i,k)=(dltR-2.0_r8*dltL)**2
         END DO
       END DO
       cff=1.0E-14_r8
       DO k=2,N(ng)-2
         DO i=Istr,Iend
           dltL=MAX(cff,WL(i,k  ))
           dltR=MAX(cff,WR(i,k+1))
           bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
           bL(i,k+1)=bR(i,k)
         END DO
       END DO
       DO i=Istr,Iend
         FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
         bL(i,N(ng))=bR(i,N(ng)-1)
         bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
         bL(i,N(ng))=bR(i,N(ng)-1)
         bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
         bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
         bL(i,N(ng))=qc(i,N(ng))       ! conditions
         bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
         bR(i,1)=bL(i,2)
         bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
         bR(i,1)=bL(i,2)
         bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
         bL(i,2)=qc(i,1)               ! bottom grid boxes are
         bR(i,1)=qc(i,1)               ! re-assumed to be
         bL(i,1)=qc(i,1)               ! piecewise constant.
#endif
       END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
       DO k=1,N(ng)
         DO i=Istr,Iend
           dltR=bR(i,k)-qc(i,k)
           dltL=qc(i,k)-bL(i,k)
           cffR=2.0_r8*dltR
           cffL=2.0_r8*dltL
           IF ((dltR*dltL).lt.0.0_r8) THEN
             dltR=0.0_r8
             dltL=0.0_r8
           ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
             dltR=cffL
           ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
             dltL=cffR
           END IF
           bR(i,k)=qc(i,k)+dltR
           bL(i,k)=qc(i,k)-dltL
         END DO
       END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
       !cff=dtdays*ABS(wsink(ng))
       cff=dt(ng)*ABS(wsink(ng)) ! RD : in cobalt wsink is in m.s-1
       !cff=86400*dt(ng)*ABS(wsink(ng)) ! RD test
       DO k=1,N(ng)
         DO i=Istr,Iend
           FC(i,k-1)=0.0_r8
           WL(i,k)=z_w(i,j,k-1)+cff
           WR(i,k)=Hz(i,j,k)*qc(i,k)
           ksource(i,k)=k
         END DO
       END DO
       DO k=1,N(ng)
         DO ks=k,N(ng)-1
           DO i=Istr,Iend
             IF (WL(i,k).gt.z_w(i,j,ks)) THEN
               ksource(i,k)=ks+1
               FC(i,k-1)=FC(i,k-1)+WR(i,ks)
             END IF
           END DO
         END DO
       END DO
!
!  Finalize computation of flux: add fractional part.
!
       DO k=1,N(ng)
         DO i=Istr,Iend
           ks=ksource(i,k)
           cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
           FC(i,k-1)=FC(i,k-1)+                                    &
  &                  Hz(i,j,ks)*cu*                                &
  &                  (bL(i,ks)+                                    &
  &                   cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
  &                       (1.5_r8-cu)*                             &
  &                       (bR(i,ks)+bL(i,ks)-                      &
  &                        2.0_r8*qc(i,ks))))
         END DO
       END DO

       ! update with sinking flux
       ! RD TODO this is probably not the right arrays to change !!!
       !DO k=1,N(ng)
       DO k=1,N(ng)
         DO i=Istr,Iend

         !FC(i,0) = 0.0d0

         IF ( ibio == indet ) THEN
            ndet_sinking(i,j,k) = r_dt * (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
#ifdef DIAGNOSTICS_BIO
            !DiaBio3d(i,j,k,indet_b4sink) = Hz(i,j,k) * t(i,j,k,nstp,indet) 
            !DiaBio3d(i,j,k,indet_b4sink) = Hz(i,j,k) * qc(i,k)
            !DiaBio3d(i,j,k,indet_b4sink) = 0.0d0
            !DiaBio3d(i,j,k,indet_flx)    = FC(i,k-1) !FC(N) is zero anyway
            !DiaBio3d(i,j,k,indet_afsink) = Hz(i,j,k) * ( t(i,j,k,nstp,indet) + ndet_sinking(i,j,k) * dt(ng) )
            !DiaBio3d(i,j,k,indet_afsink) = Hz(i,j,k) * t(i,j,k,nstp,indet) + (FC(i,k)-FC(i,k-1))
            !DiaBio3d(i,j,k,indet_afsink) = Hz(i,j,k) * qc(i,k) + (FC(i,k)-FC(i,k-1))
            !DiaBio3d(i,j,k,indet_afsink) = (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
#endif
         ELSEIF ( ibio == isidet ) THEN
            sidet_sinking(i,j,k) = r_dt * (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
         ELSEIF ( ibio == icadet_calc ) THEN
            cadet_calc_sinking(i,j,k) = r_dt * (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
         ELSEIF ( ibio == icadet_arag ) THEN
            cadet_arag_sinking(i,j,k) = r_dt * (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
         ELSEIF ( ibio == ipdet ) THEN
            pdet_sinking(i,j,k) = r_dt * (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
            !DiaBio3d(i,j,k,ipdet_b4sink) = t(i,j,k,nstp,ipdet) 
            !DiaBio3d(i,j,k,ipdet_afsink) = t(i,j,k,nstp,ipdet) + pdet_sinking(i,j,k) * dt(ng)
         ELSEIF ( ibio == ifedet ) THEN
            fedet_sinking(i,j,k) = r_dt * (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
            !DiaBio3d(i,j,k,ifedet_b4sink) = t(i,j,k,nstp,ifedet) 
            !DiaBio3d(i,j,k,ifedet_afsink) = t(i,j,k,nstp,ifedet) + fedet_sinking(i,j,k) * dt(ng)
         ELSEIF ( ibio == ilithdet ) THEN
            lithdet_sinking(i,j,k) = r_dt * (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
         ENDIF

         END DO
       END DO

       ! Set bottom fluxes in Cobalt units :
       ! FC in [mol.kg-1.m] cos' (see (1)), rho0 in [kg.m-3], dt in [s]
       ! flux is then in [mol.m-2.s-1]

       DO i=Istr,Iend
         IF ( ibio == indet ) THEN
             cobalt%f_ndet_btf(i,j,1) = FC(i,0) * rho0 * r_dt
         ELSEIF ( ibio == isidet ) THEN
             cobalt%f_sidet_btf(i,j,1) = FC(i,0) * rho0 * r_dt
         ELSEIF ( ibio == ipdet ) THEN
             cobalt%f_pdet_btf(i,j,1) = FC(i,0) * rho0 * r_dt
         ELSEIF ( ibio == icadet_arag ) THEN
             cobalt%f_cadet_arag_btf(i,j,1) = FC(i,0) * rho0 * r_dt
         ELSEIF ( ibio == icadet_calc ) THEN
             cobalt%f_cadet_calc_btf(i,j,1) = FC(i,0) * rho0 * r_dt
         ELSEIF ( ibio == ifedet ) THEN
             cobalt%f_fedet_btf(i,j,1) = FC(i,0) * rho0 * r_dt
         ELSEIF ( ibio == ilithdet ) THEN
             cobalt%f_lithdet_btf(i,j,1) = FC(i,0) * rho0 * r_dt
         ELSE
             PRINT *, 'no such sinking flux' ; STOP
         ENDIF
       ENDDO

    ENDDO ! sink loop

 ENDDO ! j loop


 ! RD : check that all of those bottom fluxes get back in the system some way or
 ! another !!!

!
!-------------------------------------------------------------------------------------------------
! 7: Sedimentary fluxes/transformations
!-------------------------------------------------------------------------------------------------
!

!
! Nitrogen flux from the sediments
! 

  ! RD dev notes : bottom will always be k=1 so I got rid of grid_kmt and replaced by k=1
  k=1
  DO j=Jstr,Jend
    DO i=Istr,Iend

      IF (cobalt%f_ndet_btf(i,j,1) .gt. 0.0) THEN

      ! fpoc_bottom in mmoles C m-2 day-1 for burial relationship
      fpoc_btm = (cobalt%f_ndet_btf(i,j,1)*cobalt%c_2_n*sperd*1000.0)
      cobalt%frac_burial(i,j) = (0.013 + 0.53*fpoc_btm**2.0)/((7.0+fpoc_btm)**2.0)


#ifdef COBALT_CONSERVATION_TEST
      ! RD dev notes ; this is to test conservation
      cobalt%frac_burial(i,j) = 0.0d0
#endif


      cobalt%fndet_burial(i,j) = cobalt%frac_burial(i,j)*cobalt%f_ndet_btf(i,j,1)
      cobalt%fpdet_burial(i,j) = cobalt%frac_burial(i,j)*cobalt%f_pdet_btf(i,j,1)

      ! fpoc_bottom in micromoles C cm-2 day-1 for denitrification relationship, cap at 43
      ! to prevent anomalous extrapolation of the relationship
      log_fpoc_btm = log(min(43.0,0.1*fpoc_btm))
      cobalt%fno3denit_sed(i,j) = min(cobalt%f_no3(i,j,k)*cobalt%Rho_0*  &
    & r_dt,                                                              &
    & min((cobalt%f_ndet_btf(i,j,1)-cobalt%fndet_burial(i,j))*           &
    &     cobalt%n_2_n_denit,                                            &
    &     10.0**(-0.9543+0.7662*log_fpoc_btm - 0.235*log_fpoc_btm**2.0)/ &
    &     (cobalt%c_2_n*sperd*100.0)*                                    &
    &     cobalt%n_2_n_denit*cobalt%f_no3(i,j,k)/                        &
    &     (cobalt%k_no3_denit + cobalt%f_no3(i,j,k))))

#ifdef COBALT_CONSERVATION_TEST
      cobalt%fno3denit_sed(i,j) = 0.0d0
#endif

      IF (cobalt%f_o2(i,j,k) .gt. cobalt%o2_min) THEN  !{
          cobalt%fnoxic_sed(i,j) = max(0.0, min(cobalt%f_o2(i,j,k)*      &
    &              cobalt%Rho_0*r_dt*(1.0/cobalt%o2_2_nh4),              &
    &              cobalt%f_ndet_btf(i,j,1) - cobalt%fndet_burial(i,j) - &
    &              cobalt%fno3denit_sed(i,j)/cobalt%n_2_n_denit))
      ELSE
          cobalt%fnoxic_sed(i,j) = 0.0
      ENDIF !}

      cobalt%fno3denit_sed(i,j) = cobalt%fno3denit_sed(i,j) +            &
    &             min(cobalt%f_no3(i,j,k)*cobalt%Rho_0*r_dt-             &
    &                 cobalt%fno3denit_sed(i,j),                         &
    &                (cobalt%f_ndet_btf(i,j,1)-cobalt%fnoxic_sed(i,j)-   &
    &                 cobalt%fndet_burial(i,j) -                         &
    &                 cobalt%fno3denit_sed(i,j)/cobalt%n_2_n_denit)*     &
    &                 cobalt%n_2_n_denit)

#ifdef COBALT_CONSERVATION_TEST
      cobalt%fno3denit_sed(i,j) = 0.0d0
#endif

      cobalt%fnfeso4red_sed(i,j) = max(0.0, cobalt%f_ndet_btf(i,j,1)-    &
    &                              cobalt%fnoxic_sed(i,j)-               &
    &                              cobalt%fndet_burial(i,j)-cobalt%fno3denit_sed(i,j)/cobalt%n_2_n_denit)

     ELSE

      cobalt%fnfeso4red_sed(i,j) = 0.0d0
      cobalt%fno3denit_sed(i,j)  = 0.0d0
      cobalt%fnoxic_sed(i,j)     = 0.0d0

     ENDIF

    ! iron from sediment 
    ! read from file :
    cobalt%ffe_sed(i,j)  = ironsed(i,j) / rho0 / Hz(i,j,1)
    ! units : mol.kg-1.s-1 = mol.m-2.s-1 / ( kg.m-3 * m )

    ! THE original code from cobalt uses ndet bottom flux as a proxy to diagnose
    ! iron sources from sediments
    ! cobalt%ffe_sed(i,j) = cobalt%fe_2_n_sed * cobalt%f_ndet_btf(i,j,1)

!#ASKCHARLIE if we need this
#ifdef COBALT_CONSERVATION_TEST
    cobalt%ffe_sed(i,j) = 0.0d0
#endif
#ifdef DIAGNOSTICS_BIO
    ! output in diagnostic variable
    DiaBio2d(i,j,iironsed_flx) = DiaBio2d(i,j,iironsed_flx) + cobalt%ffe_sed(i,j) *  dt(ng) 
#endif


    !
    ! Calcium carbonate flux and burial
    !
!    cobalt%fcased_redis(i,j)=max(0.0,min(0.5*cobalt%f_cased(i,j,1)*r_dt, &
!  &                            min(0.5 * cobalt%f_cadet_calc_btf(i,j,1), &
!  &                   0.165 * cobalt%f_ndet_btf(i,j,1) * cobalt%c_2_n) + &
!  &           0.1244 / spery * max(0.0, 1.0 - cobalt%omega_calc(i,j,k) + &
!  &    4.38 * cobalt%f_ndet_btf(i,j,1) * cobalt%c_2_n * spery)**(2.91) * &
!  &                       max(1.0, cobalt%f_lithdet_btf(i,j,1) * spery + &
!  &                             cobalt%f_cadet_calc_btf(i,j,1) * 100.0 * &
!  &                             spery)**(-2.55) * cobalt%f_cased(i,j,1)))

    ! simpler bottom flux
    cobalt%fcased_redis(i,j)=cobalt%f_cadet_calc_btf(i,j,1)

    cobalt%fcased_burial(i,j) = max(0.0, cobalt%f_cadet_calc_btf(i,j,1)* &
  &                             cobalt%f_cased(i,j,1) / 8.1e3)




! RD : fcased_burial needs to be set at zero for conservation ?
#ifdef COBALT_CONSERVATION_TEST
    cobalt%fcased_burial(i,j) = 0.0d0
#endif

    cobalt%f_cased(i,j,1) = cobalt%f_cased(i,j,1) +                      &
  &                        (cobalt%f_cadet_calc_btf(i,j,1) -             &
  &                         cobalt%fcased_redis(i,j) -                   &
  &                         cobalt%fcased_burial(i,j)) / cobalt%z_sed *  &
  &                         dt(ng) * rmask(i,j)              


   !
   ! Bottom flux boundaries passed to the vertical mixing routine 
   !
    cobalt%b_alk(i,j) = - 2.0*(cobalt%fcased_redis(i,j)+                 &
   &                           cobalt%f_cadet_arag_btf(i,j,1)) -         &
   &                          (cobalt%f_ndet_btf(i,j,1) -                &
   &                           cobalt%fndet_burial(i,j)) +               &
   &                           cobalt%alk_2_n_denit * cobalt%fno3denit_sed(i,j)


    cobalt%b_dic(i,j) =  - cobalt%fcased_redis(i,j) -                    &
   &                       cobalt%f_cadet_arag_btf(i,j,1) -              &
   &                      (cobalt%f_ndet_btf(i,j,1) - cobalt%fndet_burial(i,j)) * cobalt%c_2_n


! RD :this simehow has to be passed to bt

#ifdef COBALT_CONSERVATION_TEST
    cobalt%b_fed(i,j) = - cobalt%f_fedet_btf(i,j,1)
#else
    cobalt%b_fed(i,j) = - cobalt%ffe_sed(i,j)
#endif
    cobalt%b_nh4(i,j) = - cobalt%f_ndet_btf(i,j,1) + cobalt%fndet_burial(i,j)
    cobalt%b_no3(i,j) = cobalt%fno3denit_sed(i,j)
    cobalt%b_o2(i,j)  = cobalt%o2_2_nh4 * (cobalt%fnoxic_sed(i,j) + cobalt%fnfeso4red_sed(i,j))
    cobalt%b_po4(i,j) = - cobalt%f_pdet_btf(i,j,1) + cobalt%fpdet_burial(i,j)
    cobalt%b_sio4(i,j)= - cobalt%f_sidet_btf(i,j,1)

   ! these have to be added to the bottom level but with opposite sign


! ndet_btf = ndet * wsink * dt
! put source/sinks at each level for Xsink = X * wsink *dt -------------- x1 * w * dt positive


!                                                        ---------------- x2 * w * dt negative
! X2 = X2 + (x1 * wsink * dt) - (x2 * wsink * dt)
! X_btf = X_bottom * wsinl * dt

    ENDDO
  ENDDO

  DO k=2,N(ng)
    DO j=Jstr,Jend
      DO i=Istr,Iend

      cobalt%f_cased(i,j,k) = 0.0

      ENDDO
    ENDDO
  ENDDO
  

#ifdef DIAGNOSTICS_BIO
  DO j=Jstr,Jend
    DO i=Istr,Iend
! array overflow !!!
!       CALL cobalt_vertical_integral( cobalt%jprod_ndet(i,j,:), &
! &          Hz(i,j,:), z_w(i,j,:), 1000.0d0, 100.0d0 , UBk, twodim_int_diag(i,j) )

!    DiaBio2d(i,j,ijprod_ndet_100)  = DiaBio2d(i,j,iironsed_flx) + twodim_int_diag(i,j)  * 86400 * 6.625 * 1000 * 12

!       CALL cobalt_vertical_integral( cobalt%jremin_ndet(i,j,:), &
! &          Hz(i,j,:), z_w(i,j,:), 1000.0d0, 100.0d0 , UBk, twodim_int_diag(i,j) )
! &          Hz(i,j,:), z_w(i,j,:), omn(i,j), 100.0d0 , UBk, cobalt%jremin_ndet_100(i,j) )

!    DiaBio2d(i,j,ijremin_ndet_100) = DiaBio2d(i,j,ijremin_ndet_100) + twodim_int_diag(i,j) * 86400 * 6.625 * 1000 * 12

    ENDDO
  ENDDO
#endif


! DO_NOTHING 
#endif

!
!-----------------------------------------------------------------------
! 8: Source/sink calculations 
!-----------------------------------------------------------------------
!
! * Phytoplankton Nitrogen and Phosphorus
!
!
!   *** Diazotrophic Phytoplankton Nitrogen
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jndi(i,j,k) = phyto(DIAZO)%mu(i,j,k)*phyto(DIAZO)%f_n(i,j,k) -  &
&                      phyto(DIAZO)%jzloss_n(i,j,k)   -                  &
&                      phyto(DIAZO)%jhploss_n(i,j,k)  -                  &
&                      phyto(DIAZO)%jaggloss_n(i,j,k) -                  &
&                      phyto(DIAZO)%jvirloss_n(i,j,k) - phyto(DIAZO)%jexuloss_n(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Large Phytoplankton Nitrogen
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jnlg(i,j,k) = phyto(LARGE)%mu(i,j,k)*phyto(LARGE)%f_n(i,j,k) -  &
&                      phyto(LARGE)%jzloss_n(i,j,k)   -                  &
&                      phyto(LARGE)%jhploss_n(i,j,k)  -                  &
&                      phyto(LARGE)%jaggloss_n(i,j,k) -                  &
&                      phyto(LARGE)%jvirloss_n(i,j,k) -                  &
&                      phyto(LARGE)%jexuloss_n(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Small Phytoplankton Nitrogen
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jnsm(i,j,k) = phyto(SMALL)%mu(i,j,k)*phyto(SMALL)%f_n(i,j,k) -  &
&                      phyto(SMALL)%jzloss_n(i,j,k) -                    &
&                      phyto(SMALL)%jhploss_n(i,j,k) -                   &
&                      phyto(SMALL)%jaggloss_n(i,j,k) -                  &
&                      phyto(SMALL)%jvirloss_n(i,j,k) -                  &
&                      phyto(SMALL)%jexuloss_n(i,j,k)                                         

  ENDDO ; ENDDO ; ENDDO
!
! * Phytoplankton Silicon and Iron
!
!
!   *** Large Phytoplankton Silicon
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jsilg(i,j,k) = phyto(LARGE)%juptake_sio4(i,j,k)  -              &
&                       phyto(LARGE)%jzloss_sio2(i,j,k)   -              &
&                       phyto(LARGE)%jhploss_sio2(i,j,k)  -              &
&                       phyto(LARGE)%jaggloss_sio2(i,j,k) -              &
&                       phyto(LARGE)%jvirloss_sio2(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Diazotrophic Phytoplankton Iron
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jfedi(i,j,k) = phyto(DIAZO)%juptake_fe(i,j,k)  -                &
&                       phyto(DIAZO)%jzloss_fe(i,j,k)   -                &
&                       phyto(DIAZO)%jhploss_fe(i,j,k)  -                &
&                       phyto(DIAZO)%jaggloss_fe(i,j,k) -                &
&                       phyto(DIAZO)%jvirloss_fe(i,j,k) -                &
&                       phyto(DIAZO)%jexuloss_fe(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Large Phytoplankton Iron
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jfelg(i,j,k) = phyto(LARGE)%juptake_fe(i,j,k)  -                & 
&                       phyto(LARGE)%jzloss_fe(i,j,k)   -                &
&                       phyto(LARGE)%jhploss_fe(i,j,k)  -                &
&                       phyto(LARGE)%jaggloss_fe(i,j,k) -                &
&                       phyto(LARGE)%jvirloss_fe(i,j,k) -                &
&                       phyto(LARGE)%jexuloss_fe(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Small Phytoplankton Iron
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jfesm(i,j,k) = phyto(SMALL)%juptake_fe(i,j,k)  -                &
&                       phyto(SMALL)%jzloss_fe(i,j,k)   -                &
&                       phyto(SMALL)%jhploss_fe(i,j,k)  -                &
&                       phyto(SMALL)%jaggloss_fe(i,j,k) -                &
&                       phyto(SMALL)%jvirloss_fe(i,j,k) -                &
&                       phyto(SMALL)%jexuloss_fe(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Bacteria
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jnbact(i,j,k) = bact(1)%jprod_n(i,j,k)    -                     &
&                        bact(1)%jzloss_n(i,j,k)   -                     &
&                        bact(1)%jvirloss_n(i,j,k) -                     &
&                        bact(1)%jhploss_n(i,j,k)  

  ENDDO ; ENDDO ; ENDDO
!
!
! * Zooplankton 
!
!
!   *** Small zooplankton
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jnsmz(i,j,k) = zoo(1)%jprod_n(i,j,k) - zoo(1)%jzloss_n(i,j,k) - &
&                       zoo(1)%jhploss_n(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Medium zooplankton
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jnmdz(i,j,k) = zoo(2)%jprod_n(i,j,k) - zoo(2)%jzloss_n(i,j,k) - &
&                       zoo(2)%jhploss_n(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Large zooplankton
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jnlgz(i,j,k) = zoo(3)%jprod_n(i,j,k) - zoo(3)%jzloss_n(i,j,k) - &
&                       zoo(3)%jhploss_n(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!
! * NO3
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jno3(i,j,k) =  cobalt%jnitrif(i,j,k) -                          &
&                       phyto(DIAZO)%juptake_no3(i,j,k) -                &
&                       phyto(LARGE)%juptake_no3(i,j,k) -                &
&                       phyto(SMALL)%juptake_no3(i,j,k) -                &
&                       cobalt%jno3denit_wc(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
! * Other nutrients
!
!
!   *** NH4
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jprod_nh4(i,j,k) = cobalt%jprod_nh4(i,j,k) + cobalt%jremin_ndet(i,j,k)

  cobalt%jnh4(i,j,k) = cobalt%jprod_nh4(i,j,k)         -                 &
&                      phyto(DIAZO)%juptake_nh4(i,j,k) -                 &
&                      phyto(LARGE)%juptake_nh4(i,j,k) -                 &
&                      phyto(SMALL)%juptake_nh4(i,j,k) -                 &
&                      cobalt%jnitrif(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** PO4
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jprod_po4(i,j,k) = cobalt%jprod_po4(i,j,k) + cobalt%jremin_pdet(i,j,k) 

  cobalt%jpo4(i,j,k) = cobalt%jprod_po4(i,j,k)         -                 &
&                      phyto(DIAZO)%juptake_po4(i,j,k) -                 &
&                      phyto(LARGE)%juptake_po4(i,j,k) -                 &
&                      phyto(SMALL)%juptake_po4(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** SiO4
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jsio4(i,j,k) = cobalt%jprod_sio4(i,j,k) - phyto(LARGE)%juptake_sio4(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Fed
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jprod_fed(i,j,k) = cobalt%jprod_fed(i,j,k)    +                 &
&                           cobalt%jremin_fedet(i,j,k) +                 &
&                           cobalt%jfe_coast(i,j,k)  
  cobalt%jfed(i,j,k) = cobalt%jprod_fed(i,j,k)        -                  &
&                      phyto(DIAZO)%juptake_fe(i,j,k) -                  &
&                      phyto(LARGE)%juptake_fe(i,j,k) -                  &
&                      phyto(SMALL)%juptake_fe(i,j,k) -                  &
&                      cobalt%jfe_ads(i,j,k)

  ! RD dev notes : add dust source from atmosphere
  cobalt%jfed(i,j,k) = cobalt%jfed(i,j,k) + iron_dust_src(i,j,k)

  ENDDO ; ENDDO ; ENDDO

!
!-----------------------------------------------------------------------
!     Detrital Components
!-----------------------------------------------------------------------
!
!
!   *** Cadet_arag
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jcadet_arag(i,j,k) = cobalt%jprod_cadet_arag(i,j,k) -           &
&                             cobalt%jdiss_cadet_arag(i,j,k) +           &
&                             cadet_arag_sinking(i,j,k)

#ifdef DIAGNOSTICS_BIO
  ! diag on production term
  DiaBio3d(i,j,k,ijprod_cadet_arag) = DiaBio3d(i,j,k,ijprod_cadet_arag) + cobalt%jprod_cadet_arag(i,j,k)
  ! diag on dissolution term
  DiaBio3d(i,j,k,ijdiss_cadet_arag) = DiaBio3d(i,j,k,ijdiss_cadet_arag) + cobalt%jdiss_cadet_arag(i,j,k)
#endif

  ENDDO ; ENDDO ; ENDDO
!
!   *** Cadet_calc
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jcadet_calc(i,j,k) = cobalt%jprod_cadet_calc(i,j,k) -           &
&                             cobalt%jdiss_cadet_calc(i,j,k) +           &
&                             cadet_calc_sinking(i,j,k)

#ifdef DIAGNOSTICS_BIO
  ! diag on production term
  DiaBio3d(i,j,k,ijprod_cadet_calc) = DiaBio3d(i,j,k,ijprod_cadet_calc) + cobalt%jprod_cadet_calc(i,j,k)
  ! diag on dissolution term
  DiaBio3d(i,j,k,ijdiss_cadet_calc) = DiaBio3d(i,j,k,ijdiss_cadet_calc) + cobalt%jdiss_cadet_calc(i,j,k)
#endif

  ENDDO ; ENDDO ; ENDDO
!
!   *** Fedet
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend


  cobalt%jprod_fedet(i,j,k) = cobalt%jprod_fedet(i,j,k) + cobalt%jfe_ads(i,j,k)

  cobalt%jfedet(i,j,k) = cobalt%jprod_fedet(i,j,k)   -                   &
&                        cobalt%jremin_fedet(i,j,k)  -                   &
&                        cobalt%det_jzloss_fe(i,j,k) -                   & 
&                        cobalt%det_jhploss_fe(i,j,k) +                  &
&                        fedet_sinking(i,j,k)

#ifdef DIAGNOSTICS_BIO
  ! diag on production term
  DiaBio3d(i,j,k,ijprod_fedet) = DiaBio3d(i,j,k,ijprod_fedet) + cobalt%jprod_fedet(i,j,k)
  ! diag on remineralization
  DiaBio3d(i,j,k,ijremin_fedet) = DiaBio3d(i,j,k,ijremin_fedet) + cobalt%jremin_fedet(i,j,k)
  ! diag
  DiaBio3d(i,j,k,idet_jzloss_fe) = DiaBio3d(i,j,k,idet_jzloss_fe) + cobalt%det_jzloss_fe(i,j,k)
  ! diag 
  DiaBio3d(i,j,k,idet_jhploss_fe) = DiaBio3d(i,j,k,idet_jhploss_fe) +cobalt%det_jhploss_fe(i,j,k)
#endif
  
  ENDDO ; ENDDO ; ENDDO
!
!   *** Lithdet and Lith
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jlithdet(i,j,k) = cobalt%jprod_lithdet(i,j,k) + &
&                          lithdet_sinking(i,j,k)

  cobalt%jlith(i,j,k) = - cobalt%jlithdet(i,j,k)

  ! RD dev notes : add dust source from atmosphere
  cobalt%jlith(i,j,k) = cobalt%jlith(i,j,k) + lith_dust_src(i,j,k)

#ifdef DIAGNOSTICS_BIO
  ! diag on production term
  DiaBio3d(i,j,k,ijprod_lithdet) = DiaBio3d(i,j,k,ijprod_lithdet) + cobalt%jprod_lithdet(i,j,k)
#endif

  ENDDO ; ENDDO ; ENDDO
!
!   *** Ndet
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jndet(i,j,k) = cobalt%jprod_ndet(i,j,k)   -                     &
&                       cobalt%jremin_ndet(i,j,k)  -                     &
&                       cobalt%det_jzloss_n(i,j,k) -                     &
&                       cobalt%det_jhploss_n(i,j,k) +                    &
&                       ndet_sinking(i,j,k)

#ifdef DIAGNOSTICS_BIO
  ! diag on production term
  DiaBio3d(i,j,k,ijprod_ndet) = DiaBio3d(i,j,k,ijprod_ndet) + cobalt%jprod_ndet(i,j,k)
  ! diag on remineralization
  DiaBio3d(i,j,k,ijremin_ndet) = DiaBio3d(i,j,k,ijremin_ndet) + cobalt%jremin_ndet(i,j,k)
  ! diag
  DiaBio3d(i,j,k,idet_jzloss_n) = DiaBio3d(i,j,k,idet_jzloss_n) + cobalt%det_jzloss_n(i,j,k)
  ! diag 
  DiaBio3d(i,j,k,idet_jhploss_n) = DiaBio3d(i,j,k,idet_jhploss_n) + cobalt%det_jhploss_n(i,j,k)
#endif

  ENDDO ; ENDDO ; ENDDO
!
!   *** Pdet
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jpdet(i,j,k) = cobalt%jprod_pdet(i,j,k)   -                     &
&                       cobalt%jremin_pdet(i,j,k)  -                     &
&                       cobalt%det_jzloss_p(i,j,k) -                     &
&                       cobalt%det_jhploss_p(i,j,k) +                    &
&                       pdet_sinking(i,j,k)

#ifdef DIAGNOSTICS_BIO
  ! diag on production term
  DiaBio3d(i,j,k,ijprod_pdet) = DiaBio3d(i,j,k,ijprod_pdet) + cobalt%jprod_pdet(i,j,k)
  ! diag on remineralization
  DiaBio3d(i,j,k,ijremin_pdet) = DiaBio3d(i,j,k,ijremin_pdet) + cobalt%jremin_pdet(i,j,k)
  ! diag
  DiaBio3d(i,j,k,idet_jzloss_p) = DiaBio3d(i,j,k,idet_jzloss_p) + cobalt%det_jzloss_p(i,j,k)
  ! diag 
  DiaBio3d(i,j,k,idet_jhploss_p) = DiaBio3d(i,j,k,idet_jhploss_p) + cobalt%det_jhploss_p(i,j,k)
#endif

  ENDDO ; ENDDO ; ENDDO
!
!   *** Sidet
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jsidet(i,j,k) = cobalt%jprod_sidet(i,j,k)   -                   & 
&                        cobalt%jdiss_sidet(i,j,k)   -                   &
&                        cobalt%det_jzloss_si(i,j,k) -                   &
&                        cobalt%det_jhploss_si(i,j,k) +                  &
&                        sidet_sinking(i,j,k)

#ifdef DIAGNOSTICS_BIO
  ! diag on production term
  DiaBio3d(i,j,k,ijprod_sidet) = DiaBio3d(i,j,k,ijprod_sidet) + cobalt%jprod_sidet(i,j,k)
#endif

  ENDDO ; ENDDO ; ENDDO
!
! * Dissolved Organic Matter
!
!
!   *** Labile Dissolved Organic Nitrogen
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jldon(i,j,k) = cobalt%jprod_ldon(i,j,k)                 +       &
&                       cobalt%gamma_sldon*cobalt%f_sldon(i,j,k) +       &
&                       cobalt%gamma_srdon*cobalt%f_srdon(i,j,k) -       &
&                       bact(1)%juptake_ldon(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Labile Dissolved Organic Phosphorous
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jldop(i,j,k) = cobalt%jprod_ldop(i,j,k)                 +       &
&                       cobalt%gamma_sldop*cobalt%f_sldop(i,j,k) +       &
&                       cobalt%gamma_srdop*cobalt%f_srdop(i,j,k) -       &
&                       bact(1)%juptake_ldop(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Semilabile Dissolved Organic Nitrogen
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jsldon(i,j,k) = cobalt%jprod_sldon(i,j,k) -                     &
&                        cobalt%gamma_sldon*cobalt%f_sldon(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Semilabile dissolved organic phosphorous  
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jsldop(i,j,k) = cobalt%jprod_sldop(i,j,k) -                     &
&                        cobalt%gamma_sldop*cobalt%f_sldop(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Refractory Dissolved Organic Nitrogen
! 
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jsrdon(i,j,k) = cobalt%jprod_srdon(i,j,k) -  cobalt%gamma_srdon * cobalt%f_srdon(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
!   *** Refractory dissolved organic phosphorous
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jsrdop(i,j,k) = cobalt%jprod_srdop(i,j,k) - cobalt%gamma_srdop * cobalt%f_srdop(i,j,k)

  ENDDO ; ENDDO ; ENDDO
!
! * O2
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jo2(i,j,k) =(cobalt%o2_2_no3*(phyto(DIAZO)%juptake_no3(i,j,k) + &
&                    phyto(LARGE)%juptake_no3(i,j,k)                   + &
&                    phyto(SMALL)%juptake_no3(i,j,k))                  + & 
&                    cobalt%o2_2_nh4*(phyto(DIAZO)%juptake_nh4(i,j,k)  + &
&                    phyto(LARGE)%juptake_nh4(i,j,k)                   + &
&                    phyto(SMALL)%juptake_nh4(i,j,k)                   + &  
&                    phyto(DIAZO)%juptake_n2(i,j,k))) * rmask3d(i,j,k)

       IF (cobalt%f_o2(i,j,k) .gt. cobalt%o2_min) THEN 
          cobalt%jo2(i,j,k) = cobalt%jo2(i,j,k)                       -  &
        &                     cobalt%o2_2_nh4*cobalt%jprod_nh4(i,j,k) -  &
        &                     cobalt%o2_2_nitrif*cobalt%jnitrif(i,j,k) 
       ENDIF  

  ENDDO ; ENDDO ; ENDDO

  ! RD dev notes : add contribution from O2 air/sea fluxes
  DO j=Jstr,Jend ; DO i=Istr,Iend

  cobalt%jo2(i,j,UBk) = cobalt%jo2(i,j,UBk) + airsea_o2_flx(i,j)

  ENDDO ; ENDDO 

!
! * The Carbon system
!
!
!   *** Alkalinity
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jalk(i,j,k) = (2.0 * (cobalt%jdiss_cadet_arag(i,j,k)          + &
&    cobalt%jdiss_cadet_calc(i,j,k) - cobalt%jprod_cadet_arag(i,j,k)   - &
&    cobalt%jprod_cadet_calc(i,j,k)) + phyto(DIAZO)%juptake_no3(i,j,k) + &
&    phyto(LARGE)%juptake_no3(i,j,k) + phyto(SMALL)%juptake_no3(i,j,k) + &
&    cobalt%jprod_nh4(i,j,k) - phyto(DIAZO)%juptake_nh4(i,j,k)         - & 
&    phyto(LARGE)%juptake_nh4(i,j,k) - phyto(SMALL)%juptake_nh4(i,j,k) - &
&    2.0 * cobalt%jnitrif(i,j,k) + cobalt%alk_2_n_denit * cobalt%jno3denit_wc(i,j,k))

  ENDDO ; ENDDO ; ENDDO
!
!   *** Dissolved Inorganic Carbon
!
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

  cobalt%jdic(i,j,k) =(cobalt%c_2_n * (cobalt%jno3(i,j,k)              + &
&    cobalt%jnh4(i,j,k) + cobalt%jno3denit_wc(i,j,k)                   - &
&    phyto(DIAZO)%juptake_n2(i,j,k))                                   + &
&    cobalt%jdiss_cadet_arag(i,j,k) + cobalt%jdiss_cadet_calc(i,j,k)   - &
&    cobalt%jprod_cadet_arag(i,j,k) - cobalt%jprod_cadet_calc(i,j,k)) 

  ENDDO ; ENDDO ; ENDDO

  ! RD dev notes : add contribution from CO2 air/sea fluxes
  DO j=Jstr,Jend ; DO i=Istr,Iend

  cobalt%jdic(i,j,UBk) = cobalt%jdic(i,j,UBk) + airsea_co2_flx(i,j)

  ENDDO ; ENDDO 


  ! RD dev notes :
  ! adding the bottom flux b_var to the bottom level concentration
  ! cobalt%p_var in [mol.kg-1] , cobalt%b_var in [mol.m-2.s-1]
  DO j=Jstr,Jend
    DO i=Istr,Iend
      cff_btf = 1.0d0 / ( rho0 * Hz(i,j,1) )
      !
      cobalt%jalk(i,j,1)  = cobalt%jalk(i,j,1)  - cobalt%b_alk(i,j)  * cff_btf
      cobalt%jdic(i,j,1)  = cobalt%jdic(i,j,1)  - cobalt%b_dic(i,j)  * cff_btf
      cobalt%jfed(i,j,1)  = cobalt%jfed(i,j,1)  - cobalt%b_fed(i,j)  * cff_btf
      cobalt%jnh4(i,j,1)  = cobalt%jnh4(i,j,1)  - cobalt%b_nh4(i,j)  * cff_btf
      cobalt%jno3(i,j,1)  = cobalt%jno3(i,j,1)  - cobalt%b_no3(i,j)  * cff_btf
      cobalt%jo2(i,j,1)   = cobalt%jo2(i,j,1)   - cobalt%b_o2(i,j)   * cff_btf
      cobalt%jpo4(i,j,1)  = cobalt%jpo4(i,j,1)  - cobalt%b_po4(i,j)  * cff_btf 
      cobalt%jsio4(i,j,1) = cobalt%jsio4(i,j,1) - cobalt%b_sio4(i,j) * cff_btf
    ENDDO
  ENDDO

  ! RD dev notes : Copy other BGC variables to new step (this array is not passed to dynamics)
  ! RD : Sanity check with values = iic done and successful, array with index nstp is previous step
!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend
         obgc(i,j,k,nnew,iochl)          = cobalt%f_chl(i,j,k)
         obgc(i,j,k,nnew,ioco3_ion)      = cobalt%f_co3_ion(i,j,k)
         obgc(i,j,k,nnew,iohtotal)       = cobalt%f_htotal(i,j,k)
         obgc(i,j,k,nnew,ioirr_mem)      = cobalt%f_irr_mem(i,j,k)
  ENDDO ; ENDDO ; ENDDO

  ! Set value for boundaries
  DO it=1,NOBGC
     CALL bc_r3d_tile (ng, tile,                                     &
  &                      LBi, UBi, LBj, UBj, 1, UBk,                 &
  &                      obgc(:,:,:,nnew,it))
  END DO

#ifdef BENTHIC
!
!-----------------------------------------------------------------------
!     Saving the bottom fluxes to the bt array
!-----------------------------------------------------------------------
!

! RD test
!    bt(:,:,1,nnew,icased) = bt(:,:,1,nstp,icased) + 1
#ifdef DEBUG_COBALT
    IF (iic(ng) == ntstart(ng)) THEN 
       IF ( Master ) WRITE(stdout,*) '>>>    after restart max(bt(icased)) =', MAXVAL(bt(:,:,1,:,icased))
    ENDIF
#endif

    ! RD TODO maybe this can be changed to bt(...nstp ) = cobalt%...

    bt(:,:,1,nnew,icased)          = cobalt%f_cased(:,:,1)
    bt(:,:,1,nnew,icadet_arag_btf) = cobalt%f_cadet_arag_btf(:,:,1)
    bt(:,:,1,nnew,icadet_calc_btf) = cobalt%f_cadet_calc_btf(:,:,1)
    bt(:,:,1,nnew,indet_btf)       = cobalt%f_ndet_btf(:,:,1)
    bt(:,:,1,nnew,ipdet_btf)       = cobalt%f_pdet_btf(:,:,1)
    bt(:,:,1,nnew,isidet_btf)      = cobalt%f_sidet_btf(:,:,1)

    bt(:,:,1,nstp,icased)          = bt(:,:,1,nnew,icased) 
    bt(:,:,1,nstp,icadet_arag_btf) = bt(:,:,1,nnew,icadet_arag_btf)
    bt(:,:,1,nstp,icadet_calc_btf) = bt(:,:,1,nnew,icadet_calc_btf)
    bt(:,:,1,nstp,indet_btf)       = bt(:,:,1,nnew,indet_btf) 
    bt(:,:,1,nstp,ipdet_btf)       = bt(:,:,1,nnew,ipdet_btf)
    bt(:,:,1,nstp,isidet_btf)      = bt(:,:,1,nnew,isidet_btf)

#endif

#ifdef DIAGNOSTICS_BIO
    !
    !-----------------------------------------------------------------------
    !       Save variables for diagnostics
    !-----------------------------------------------------------------------
    !

    ! RD dev notes :
    ! update 2d diagnostics variables in ROMS
    DO j=Jstr,Jend
      DO i=Istr,Iend
         DiaBio2d(i,j,icased)          = DiaBio2d(i,j,icased) + cobalt%f_cased(i,j,1)
         DiaBio2d(i,j,icadet_arag_btf) = DiaBio2d(i,j,icadet_arag_btf) + cobalt%f_cadet_arag_btf(i,j,1)
         DiaBio2d(i,j,icadet_calc_btf) = DiaBio2d(i,j,icadet_calc_btf) + cobalt%f_cadet_calc_btf(i,j,1)
         DiaBio2d(i,j,indet_btf)       = DiaBio2d(i,j,indet_btf) + cobalt%f_ndet_btf(i,j,1)
         DiaBio2d(i,j,ipdet_btf)       = DiaBio2d(i,j,ipdet_btf) + cobalt%f_pdet_btf(i,j,1)
         DiaBio2d(i,j,isidet_btf)      = DiaBio2d(i,j,isidet_btf) + cobalt%f_sidet_btf(i,j,1)
         DiaBio2d(i,j,imxl_depth)      = DiaBio2d(i,j,imxl_depth) + mxl_depth(i,j)
         DiaBio2d(i,j,imxl_level)      = DiaBio2d(i,j,imxl_level) + mxl_blev(i,j)

         DiaBio2d(i,j,ialk_btf)  = DiaBio2d(i,j,ialk_btf) + cobalt%b_alk(i,j)
         DiaBio2d(i,j,idic_btf)  = DiaBio2d(i,j,idic_btf) + cobalt%b_dic(i,j)
         DiaBio2d(i,j,ifed_btf)  = DiaBio2d(i,j,ifed_btf) + cobalt%b_fed(i,j)
         DiaBio2d(i,j,inh4_btf)  = DiaBio2d(i,j,inh4_btf) + cobalt%b_nh4(i,j)
         DiaBio2d(i,j,ino3_btf)  = DiaBio2d(i,j,ino3_btf) + cobalt%b_no3(i,j)
         DiaBio2d(i,j,io2_btf)   = DiaBio2d(i,j,io2_btf)  + cobalt%b_o2(i,j)
         DiaBio2d(i,j,ipo4_btf)  = DiaBio2d(i,j,ipo4_btf) + cobalt%b_po4(i,j)
         DiaBio2d(i,j,isio4_btf) = DiaBio2d(i,j,isio4_btf)+ cobalt%b_sio4(i,j)
      ENDDO
    ENDDO
    i=overflow ; j=overflow

    ! RD dev notes :
    ! update 3d diagnostics variables in ROMS
    DO k=1,UBk
      DO j=Jstr,Jend
        DO i=Istr,Iend
           DiaBio3d(i,j,k,ichl)          = DiaBio3d(i,j,k,ichl)     + cobalt%f_chl(i,j,k)
           DiaBio3d(i,j,k,ico3_ion)      = DiaBio3d(i,j,k,ico3_ion) + cobalt%f_co3_ion(i,j,k)
           DiaBio3d(i,j,k,ihtotal)       = DiaBio3d(i,j,k,ihtotal)  + cobalt%f_htotal(i,j,k)
           DiaBio3d(i,j,k,iirr_mem)      = DiaBio3d(i,j,k,iirr_mem) + cobalt%f_irr_mem(i,j,k)
           DiaBio3d(i,j,k,iirr_mix)      = DiaBio3d(i,j,k,iirr_mix) + cobalt%irr_mix(i,j,k)
           DiaBio3d(i,j,k,iirr_inst)     = DiaBio3d(i,j,k,iirr_inst) + cobalt%irr_inst(i,j,k)
           ! RD : debuuging variables (remove or replace later)
           !DiaBio3d(i,j,k,itheta_small)  = DK(i,j,k,2)
           !DiaBio3d(i,j,k,itheta_large)  = DK(i,j,k,3)
           !DiaBio3d(i,j,k,itheta_diazo)  = DK(i,j,k,4)
           !DiaBio3d(i,j,k,ipcm_small)    = FM(i,j,k-1)
           !DiaBio3d(i,j,k,ipcm_large)    = tmp_decay(i,j,k)
        ENDDO
      ENDDO
    ENDDO
    i=overflow ; j=overflow ; k=overflow

#endif
!
!-----------------------------------------------------------------------------------
! End: Give back the tracers arrays to the dynamics
!-----------------------------------------------------------------------------------
!
! RD dev notes :
! the max(value, 0.) makes really a difference on the results
! Remark : Xa = Xa + (Xb -Xa) * Hz

!  DO k=1,UBk ; DO j=Jstr,Jend ; DO i=Istr,Iend
  DO j=Jstr,Jend ; DO k=1,UBk ; DO i=Istr,Iend

   cff_ts = dt(ng) * rmask(i,j) * Hz(i,j,k)
#ifdef COBALT_NOSOURCE
   cff_ts = 0.0d0
#endif
   ! Nitrogen dynamics
   t(i,j,k,nnew,insm)   = t(i,j,k,nnew,insm)   + cobalt%jnsm(i,j,k)   * cff_ts
   t(i,j,k,nnew,inlg)   = t(i,j,k,nnew,inlg)   + cobalt%jnlg(i,j,k)   * cff_ts
   t(i,j,k,nnew,indi)   = t(i,j,k,nnew,indi)   + cobalt%jndi(i,j,k)   * cff_ts
   t(i,j,k,nnew,insmz)  = t(i,j,k,nnew,insmz)  + cobalt%jnsmz(i,j,k)  * cff_ts
   t(i,j,k,nnew,inmdz)  = t(i,j,k,nnew,inmdz)  + cobalt%jnmdz(i,j,k)  * cff_ts
   t(i,j,k,nnew,inlgz)  = t(i,j,k,nnew,inlgz)  + cobalt%jnlgz(i,j,k)  * cff_ts
   t(i,j,k,nnew,ildon)  = t(i,j,k,nnew,ildon)  + cobalt%jldon(i,j,k)  * cff_ts
   t(i,j,k,nnew,isldon) = t(i,j,k,nnew,isldon) + cobalt%jsldon(i,j,k) * cff_ts
   t(i,j,k,nnew,isrdon) = t(i,j,k,nnew,isrdon) + cobalt%jsrdon(i,j,k) * cff_ts
   t(i,j,k,nnew,inbact) = t(i,j,k,nnew,inbact) + cobalt%jnbact(i,j,k) * cff_ts
   t(i,j,k,nnew,inh4)   = t(i,j,k,nnew,inh4)   + cobalt%jnh4(i,j,k)   * cff_ts
   t(i,j,k,nnew,ino3)   = t(i,j,k,nnew,ino3)   + cobalt%jno3(i,j,k)   * cff_ts
   t(i,j,k,nnew,indet)  = t(i,j,k,nnew,indet)  + cobalt%jndet(i,j,k)  * cff_ts
#ifdef COBALT_MINERALS
   ! Biogenic Minerals and Lithogenic Materials
   t(i,j,k,nnew,isio4)       = t(i,j,k,nnew,isio4)       + cobalt%jsio4(i,j,k)       * cff_ts
   t(i,j,k,nnew,isilg)       = t(i,j,k,nnew,isilg)       + cobalt%jsilg(i,j,k)       * cff_ts
   t(i,j,k,nnew,isidet)      = t(i,j,k,nnew,isidet)      + cobalt%jsidet(i,j,k)      * cff_ts
   t(i,j,k,nnew,icadet_calc) = t(i,j,k,nnew,icadet_calc) + cobalt%jcadet_calc(i,j,k) * cff_ts
   t(i,j,k,nnew,icadet_arag) = t(i,j,k,nnew,icadet_arag) + cobalt%jcadet_arag(i,j,k) * cff_ts
   t(i,j,k,nnew,ilith)       = t(i,j,k,nnew,ilith)       + cobalt%jlith(i,j,k)       * cff_ts
   t(i,j,k,nnew,ilithdet)    = t(i,j,k,nnew,ilithdet)    + cobalt%jlithdet(i,j,k)    * cff_ts
#endif
#ifdef COBALT_PHOSPHORUS
   ! Phosporus Dynamics
   t(i,j,k,nnew,ildop)  = t(i,j,k,nnew,ildop)  + cobalt%jldop(i,j,k)  * cff_ts
   t(i,j,k,nnew,isldop) = t(i,j,k,nnew,isldop) + cobalt%jsldop(i,j,k) * cff_ts
   t(i,j,k,nnew,isrdop) = t(i,j,k,nnew,isrdop) + cobalt%jsrdop(i,j,k) * cff_ts
   t(i,j,k,nnew,ipo4)   = t(i,j,k,nnew,ipo4)   + cobalt%jpo4(i,j,k)   * cff_ts
   t(i,j,k,nnew,ipdet)  = t(i,j,k,nnew,ipdet)  + cobalt%jpdet(i,j,k)  * cff_ts
#endif
#ifdef COBALT_IRON
   ! Iron Dynamics
   t(i,j,k,nnew,ifesm)  = t(i,j,k,nnew,ifesm)  + cobalt%jfesm(i,j,k)  * cff_ts
   t(i,j,k,nnew,ifedi)  = t(i,j,k,nnew,ifedi)  + cobalt%jfedi(i,j,k)  * cff_ts
   t(i,j,k,nnew,ifelg)  = t(i,j,k,nnew,ifelg)  + cobalt%jfelg(i,j,k)  * cff_ts
   t(i,j,k,nnew,ifed)   = t(i,j,k,nnew,ifed)   + cobalt%jfed(i,j,k)   * cff_ts
   t(i,j,k,nnew,ifedet) = t(i,j,k,nnew,ifedet) + cobalt%jfedet(i,j,k) * cff_ts
#endif
#ifdef COBALT_CARBON
   ! Oxygen, Carbon and Alkalinity
   t(i,j,k,nnew,io2)  = t(i,j,k,nnew,io2)  + cobalt%jo2(i,j,k)  * cff_ts
   t(i,j,k,nnew,idic) = t(i,j,k,nnew,idic) + cobalt%jdic(i,j,k) * cff_ts
   t(i,j,k,nnew,ialk) = t(i,j,k,nnew,ialk) + cobalt%jalk(i,j,k) * cff_ts
#endif

  ENDDO ; ENDDO ; ENDDO
  i=overflow ; j=overflow ; k=overflow

#ifdef TIMESERIES
!----------------------------------------------------------------------------------------
!
!      DIAGNOSTICS 3D
!
!----------------------------------------------------------------------------------------


   ! init the tracers to zero
   total_tracer(:) = 0.0d0

   DO it=1,NT(ng)
     DO k=1,UBk
       DO j=Jstr,Jend
         DO i=Istr,Iend
          !total_tracer(it) = total_tracer(it) + ( t(i,j,k,nstp,it) * omn(i,j) * rmask(i,j) * Hz(i,j,k) )
          total_tracer(it) = total_tracer(it) + ( t(i,j,k,nnew,it) * omn(i,j) * rmask(i,j) )
         ENDDO
       ENDDO
     ENDDO
   ENDDO

   ! sum over all tiles
   CALL mp_reduce (ng, iNLM, ndiag_int3, total_tracer, op_handle)


   ! same for volume
   zvolume = 0.0d0
   DO k=1,UBk
     DO j=Jstr,Jend
       DO i=Istr,Iend
          zvolume = zvolume + ( omn(i,j) * rmask(i,j) * Hz(i,j,k) ) 
       ENDDO
     ENDDO
   ENDDO

   CALL mp_reduce (ng, iNLM, 1, zvolume, 'SUM')

   DO it=1,NT(ng)
      tms(1,1,it) = total_tracer(it) / zvolume
   ENDDO

#endif

! DO_NOTHING
!!!!#endif


      RETURN
      END SUBROUTINE biology_tile
!-------------------------------------------------------------------------------

   SUBROUTINE check_overdrive(rhs,lhs,ratio,bioeqvar)

   USE mod_parallel

   IMPLICIT NONE

   REAL(8),INTENT(in) :: rhs, lhs ! right/left hand side
   REAL(8),INTENT(in)                  :: ratio    ! ratio allowed
   CHARACTER(len=*),INTENT(in)         :: bioeqvar ! equation for var

   IF( ABS(rhs) > ratio * ABS(lhs) ) THEN
     PRINT *, 'Message from core', MyRank, ': Boy you are in trouble here'
     PRINT *, 'Message from core', MyRank, ' : ', TRIM(bioeqvar), ' right hand side is way too high'
     STOP
   ENDIF

   END SUBROUTINE

!---------------------------------------------------------------------------------------------

   SUBROUTINE cobalt_vertical_integral(lar_tracer, lar_layerthick, lar_depth_w_point, \
              lsr_cell_area, lsr_depth_integration, lsi_number_levels, lsr_tracer_vert_integral)

   IMPLICIT NONE

   REAL(8),DIMENSION(:),INTENT(in)  :: lar_tracer, lar_layerthick, lar_depth_w_point
   REAL(8),INTENT(in)               :: lsr_cell_area, lsr_depth_integration
   INTEGER(4),INTENT(in)            :: lsi_number_levels
   REAL(8),INTENT(out)              :: lsr_tracer_vert_integral

   INTEGER :: k, kdi, overflow

   overflow=1e6

   ! while below depth of integration, keep incrementing kdi
   DO k=1,lsi_number_levels
      IF (lar_depth_w_point(k) < -1.0 * lsr_depth_integration )  kdi = k
   ENDDO
   k=overflow

   lsr_tracer_vert_integral = 0.
   DO k=lsi_number_levels,kdi,-1
      lsr_tracer_vert_integral = lsr_tracer_vert_integral + lar_tracer(k) * lar_layerthick(k) * lsr_cell_area
   ENDDO

   END SUBROUTINE

