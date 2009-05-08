      MODULE mod_forces
!
!svn $Id: mod_forces.F 975 2009-05-05 22:51:13Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Surface momentum stresses.                                          !
!                                                                      !
!  sustr        Kinematic surface momentum flux (wind stress) in       !
!                 the XI-direction (m2/s2) at horizontal U-points.     !
!  sustrG       Latest two-time snapshots of input "sustr" grided      !
!                 data used for interpolation.                         !
!  svstr        Kinematic surface momentum flux (wind stress) in       !
!                 the ETA-direction (m2/s2) at horizontal V-points.    !
!  svstrG       Latest two-time snapshots of input "svstr" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Bottom momentum stresses.                                           !
!                                                                      !
!  bustr        Kinematic bottom momentum flux (bottom stress) in      !
!                 the XI-direction (m2/s2) at horizontal U-points.     !
!  bvstr        Kinematic bottom momentum flux (bottom stress) in      !
!                 ETA-direction (m2/s2) at horizontal V-points.        !
!                                                                      !
!  Surface wind induced waves.                                         !
!                                                                      !
!  Hwave        Surface wind induced wave height (m).                  !
!  HwaveG       Latest two-time snapshots of input "Hwave" grided      !
!                 data used for interpolation.                         !
!  Dwave        Surface wind induced wave direction (radians).         !
!  DwaveG       Latest two-time snapshots of input "Dwave" grided      !
!                 data used for interpolation.                         !
!  Lwave        Mean surface wavelength read in from swan output       !
!  LwaveG       Latest two-time snapshots of input "Lwave" grided      !
!                 data used for interpolation.                         !
!  Pwave_top    Wind induced surface wave period (s).                  !
!  Pwave_topG   Latest two-time snapshots of input "Pwave_top" grided  !
!                 data used for interpolation.                         !
!  Pwave_bot    Wind induced bottom wave period (s).                   !
!  Pwave_botG   Latest two-time snapshots of input "Pwave_bot" grided  !
!                 data used for interpolation.                         !
!  Ub_swan      Bottom orbital velocity read in from swan output       !
!  Ub_swanG     Latest two-time snapshots of input "Ub_swan" grided    !
!                 data used for interpolation.                         !
!  wave_dissip  Wave dissipation                                       !
!  wave_dissipG Latest two-time snapshots of input "wave_dissip"       !
!                 gridded data used for interpolation.                 !
!  Wave_break   Percent of wave breaking for use with roller model.    !
!  Wave_breakG  Latest two-time snapshots of input "wave_break"        !
!                 gridded data used for interpolation.                 !
!                                                                      !
!  Solar shortwave radiation flux.                                     !
!                                                                      !
!  srflx        Kinematic surface shortwave solar radiation flux       !
!                 (Celsius m/s) at horizontal RHO-points               !
!  srflxG       Latest two-time snapshots of input "srflx" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Cloud fraction.                                                     !
!                                                                      !
!  cloud        Cloud fraction (percentage/100).                       !
!  cloudG       Latest two-time snapshots of input "cloud" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface heat fluxes, Atmosphere-Ocean bulk parameterization.        !
!                                                                      !
!  lhflx        Kinematic net sensible heat flux (degC m/s).           !
!  lrflx        Kinematic net longwave radiation (degC m/s).           !
!  shflx        Kinematic net sensible heat flux (degC m/s).           !
!                                                                      !
!  Surface air humidity.                                               !
!                                                                      !
!  Hair         Surface air specific (g/kg) or relative humidity       !
!                 (percentage).                                        !
!  HairG        Latest two-time snapshots of input "Hair" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface air pressure.                                               !
!                                                                      !
!  Pair         Surface air pressure (mb).                             !
!  PairG        Latest two-time snapshots of input "Pair" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface air temperature.                                            !
!                                                                      !
!  Tair         Surface air temperature (Celsius)                      !
!  TairG        Latest two-time snapshots of input "Tair" grided       !
!                 data used for interpolation.                         !
!  Surface Winds.                                                      !
!                                                                      !
!  Uwind        Surface wind in the XI-direction (m/s) at              !
!                 horizontal RHO-points.                               !
!  UwindG       Latest two-time snapshots of input "Uwind" grided      !
!                 data used for interpolation.                         !
!  Vwind        Surface wind in the ETA-direction (m/s) at             !
!                 horizontal RHO-points.                               !
!  VwindG       Latest two-time snapshots of input "Vwind" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Rain fall rate.                                                     !
!                                                                      !
!  evap         Evaporation rate (kg/m2/s).                            !
!  rain         Rain fall rate (kg/m2/s).                              !
!  rainG        Latest two-time snapshots of input "rain" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface tracer fluxes.                                              !
!                                                                      !
!  stflx        Kinematic surface flux of tracer type variables        !
!                 (temperature: degC m/s; salinity: PSU m/s) at        !
!                 horizontal RHO-points.                               !
!  stflxG       Latest two-time snapshots of input "stflx" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Bottom tracer fluxes.                                               !
!                                                                      !
!  btflx        Kinematic bottom flux of tracer type variables         !
!                 (temperature: degC m/s; salinity: PSU m/s) at        !
!                horizontal RHO-points.                                !
!  btflxG       Latest two-time snapshots of input "btflx" grided      !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface heat flux correction.                                       !
!                                                                      !
!  dqdt         Kinematic surface net heat flux sensitivity to SST,    !
!                 d(Q)/d(SST), (m/s).                                  !
!  dqdtG        Latest two-time snapshots of input "dqdt" grided       !
!                 data used for interpolation.                         !
!  sst          Sea surface temperature (Celsius).                     !
!  sstG         Latest two-time snapshots of input "sst" grided        !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface freshwater flux correction.                                 !
!                                                                      !
!  sss          Sea surface salinity (PSU).                            !
!  sssG         Latest two-time snapshots of input "sss" grided        !
!                 data used for interpolation.                         !
!  sssflx       Sea surface salinity flux correction.                  !
!  sssflxG      Latest two-time snapshots of input "sssflx" grided     !
!                 data used for interpolation.                         !
!                                                                      !
!  Surface spectral downwelling irradiance.                            !
!                                                                      !
!  SpecIr       Spectral irradiance (NBands) from 400-700 nm at        !
!                 5 nm bandwidth.                                      !
!  avcos        Cosine of average zenith angle of downwelling          !
!                 spectral photons.                                    !
!                                                                      !
!=======================================================================
!
        USE mod_kinds
        implicit none
        TYPE T_FORCES
!
!  Nonlinear model state.
!
          real(r8), pointer :: sustr(:,:)
          real(r8), pointer :: svstr(:,:)
          real(r8), pointer :: bustr(:,:)
          real(r8), pointer :: bvstr(:,:)
          real(r8), pointer :: srflx(:,:)
          real(r8), pointer :: Hair(:,:)
          real(r8), pointer :: HairG(:,:,:)
          real(r8), pointer :: Tair(:,:)
          real(r8), pointer :: TairG(:,:,:)
          real(r8), pointer :: stflx(:,:,:)
          real(r8), pointer :: btflx(:,:,:)
        END TYPE T_FORCES
        TYPE (T_FORCES), allocatable :: FORCES(:)
      CONTAINS
      SUBROUTINE allocate_forces (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( FORCES(Ngrids) )
!
!  Nonlinear model state
!
      allocate ( FORCES(ng) % sustr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % svstr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % bustr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % bvstr(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % srflx(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % Hair(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % HairG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % Tair(LBi:UBi,LBj:UBj) )
      allocate ( FORCES(ng) % TairG(LBi:UBi,LBj:UBj,2) )
      allocate ( FORCES(ng) % stflx(LBi:UBi,LBj:UBj,NT(ng)) )
      allocate ( FORCES(ng) % btflx(LBi:UBi,LBj:UBj,NT(ng)) )
      RETURN
      END SUBROUTINE allocate_forces
      SUBROUTINE initialize_forces (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, j, k
      integer :: itrc
      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrR, IstrT, IstrU, Iend, IendR, IendT
      integer :: Jstr, JstrR, JstrT, JstrV, Jend, JendR, JendT
!
      Istr =BOUNDS(ng)%Istr (tile)
      IstrR=BOUNDS(ng)%IstrR(tile)
      IstrT=BOUNDS(ng)%IstrT(tile)
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      IendR=BOUNDS(ng)%IendR(tile)
      IendT=BOUNDS(ng)%IendT(tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JstrT=BOUNDS(ng)%JstrT(tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
      JendR=BOUNDS(ng)%JendR(tile)
      JendT=BOUNDS(ng)%JendT(tile)
!
!  Set array initialization range.
!
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            FORCES(ng) % sustr(i,j) = IniVal
            FORCES(ng) % svstr(i,j) = IniVal
            FORCES(ng) % bustr(i,j) = IniVal
            FORCES(ng) % bvstr(i,j) = IniVal
            FORCES(ng) % srflx(i,j) = IniVal
            FORCES(ng) % Hair(i,j) = IniVal
            FORCES(ng) % Tair(i,j) = IniVal
            DO itrc=1,NT(ng)
              FORCES(ng) % stflx(i,j,itrc) = IniVal
              FORCES(ng) % btflx(i,j,itrc) = IniVal
            END DO
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE initialize_forces
      END MODULE mod_forces
