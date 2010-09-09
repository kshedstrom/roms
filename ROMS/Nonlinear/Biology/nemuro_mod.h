!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Nemuro ecosystem model:                              !
!                                                                      !
!  AlphaPS     Small Phytoplankton photochemical reaction coefficient: !
!                initial slope (low light) of the P-I curve,           !
!                [1/(W/m2) 1/day].                                     !
!  AlphaPL     Large Phytoplankton photochemical reaction coefficient: !
!                initial slope (low light) of the P-I curve,           !
!                [1/(W/m2) 1/day].                                     !
!  AlphaZL     Large Zooplankton assimilation efficiency,              !
!                [nondimemsional].                                     !
!  AlphaZP     Predator Zooplankton assimilation efficiency,           !
!                [nondimemsional].                                     !
!  AlphaZS     Small Zooplankton assimilation efficiency,              !
!                [nondimemsional].                                     !
!  AttPL       Light attenuation due to Large Phytoplankton, self-     !
!                shading coefficient, [m2/millimole_N].                !
!  AttPS       Light attenuation due to Small Phytoplankton, self-     !
!                shading coefficient, [m2/millimole_N].                !
!  AttSW       Light attenuation due to sea water, [1/m].              !
!  BetaPL      Large Phytoplankton photoinhibition coefficient,        !
!                [1/(W/m2) 1/day].                                     !
!  BetaPS      Small Phytoplankton photoinhibition coefficient,        !
!                [1/(W/m2) 1/day].                                     !
!  BetaZL      Large Zooplankton growth efficiency [nondimensional].   !
!  BetaZP      Predator Zooplankton growth efficiency [nondimensional].!
!  BetaZS      Small Zooplankton growth efficiency [nondimensional].   !
!  BioIter     Maximum number of iterations to achieve convergence of  !
!                the nonlinear solution.                               !
!  GammaL      Large Phytoplankton ratio of extracellular excretion to !
!                photosynthesis [nondimensional].                      !
!  GammaS      Small Phytoplankton ratio of extracellular excretion to !
!                photosynthesis [nondimensional].                      !
!  GRmaxLpl    Large Zooplankton maximum grazing rate on Large         !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxLps    Large Zooplankton maximum grazing rate on Small         !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxLzs    Small Zooplankton maximum grazing rate on Small         !
!                Zooplankton at 0 Celsius, [1/day].                    !
!  GRmaxPpl    Predator Zooplankton maximum grazing rate on Large      !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxPzl    Predator Zooplankton maximum grazing rate on Large      !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxPzs    Predator Zooplankton maximum grazing rate on Small      !
!                Zooplankton at 0 Celsius, [1/day].                    !
!  GRmaxSps    Small Zooplankton maximum grazing rate on Small         !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  GRmaxSpl    Small Zooplankton maximum grazing rate on Large         !
!                Phytoplankton at 0 Celsius, [1/day].                  !
!  KD2N        Temperature coefficient for DON to NH4 decomposition,   !
!                [1/Celsius].                                          !
!  KGppL       Large Phytoplankton temperature coefficient for         !
!                photosynthetic rate, [1/Celsius].                     !
!  KGppS       Small Phytoplankton temperature coefficient for         !
!                photosynthetic rate, [1/Celsius].                     !
!  KGraL       Large Zooplankton temperature coefficient for grazing,  !
!                [1/Celsius].                                          !
!  KGraP       Predator Zooplankton temperature coefficient for        !
!                grazing,[1/Celsius].                                  !
!  KGraS       Small Zooplankton temperature coefficient for grazing,  !
!                [1/Celsius].                                          !
!  KMorPL      Large Phytoplankton temperature coefficient for         !
!                mortality, [1/Celsius].                               !
!  KMorPS      Small Phytoplankton temperature coefficient for         !
!                mortality, [1/Celsius].                               !
!  KMorZL      Large Zooplankton temperature coefficient for           !
!                mortality, [1/Celsius].                               !
!  KMorZP      Predator Zooplankton temperature coefficient for        !
!                mortality, [1/Celsius].                               !
!  KMorZS      Small Zooplankton temperature coefficient for           !
!                mortality, [1/Celsius].                               !
!  KNit        Temperature coefficient for nitrification (NH4 to NO3)  !
!                decomposition, [1/Celsius].                           !
!  KNH4L       Large Phytoplankton half satuation constant for NH4,    !
!                [millimole_N/m3].                                     !
!  KNH4S       Small Phytoplankton half satuation constant for NH4,    !
!                [millimole_N/m3].                                     !
!  KNO3L       Large Phytoplankton half satuation constant for NO3,    !
!                [millimole_N/m3].                                     !
!  KNO3S       Small Phytoplankton half satuation constant for NO3,    !
!                [millimole_N/m3].                                     !
!  KO2S        Temperature coefficient for Opal to SiOH4 decomposition,!
!                [1/Celsius].                                          !
!  KP2D        Temperature coefficient for PON to DON decomposition,   !
!                [1/Celsius].                                          !
!  KP2N        Temperature coefficient for PON to NH4 decomposition,   !
!                [1/Celsius].                                          !
!  KPL2ZS      Small Zooplankton half-saturation coefficient for       !
!                ingestion on Large Phytoplankton [millimole_N/m3]^2.  !
!  KPL2ZL      Large Zooplankton half-saturation coefficient for       !
!                ingestion on Large Phytoplankton [millimole_N/m3]^2.  !
!  KPL2ZP      Predator Zooplankton half-saturation coefficient for    !
!                ingestion on Large Phytoplankton [millimole_N/m3]^2.  !
!  KPS2ZL      Larg Zooplankton half-saturation coefficient for        !
!                ingestion on Small Phytoplankton [millimole_N/m3]^2.  !
!  KPS2ZS      Small Zooplankton half-saturation coefficient for       !
!                ingestion on Small Phytoplankton [millimole_N/m3]^2.  !
!  KResPL      Large Phytoplankton temperature coefficient for         !
!                respiration, [1/Celsius].                             !
!  KResPS      Small Phytoplankton temperature coefficient for         !
!                respiration, [1/Celsius].                             !
!  KSiL        Large Phytoplankton half satuation constant for SiOH4,  !
!                [millimole_Si/m3].                                    !
!  KZL2ZP      Predator Zooplankton half-saturation coefficient for    !
!                ingestion on Large Zooplankton [millimole_N/m3]^2.    !
!  KZS2ZL      Large Zooplankton half-saturation coefficient for       !
!                ingestion on Small Phytoplankton [millimole_N/m3]^2.  !
!  KZS2ZP      Predator Zooplankton half-saturation coefficient for    !
!                ingestion on Small Zooplankton [millimole_N/m3]^2.    !
!  LamL        Large Zooplankton Ivlev constant, [m3/millimole_N].     !
!  LamP        Predator Zooplankton Ivlev constant, [m3/millimole_N].  !
!  LamS        Small Zooplankton Ivlev constant, [m3/millimole_N].     !
!  MorPL0      Large Phytoplankton mortality rate at 0 Celsius,        !
!                [m3/millimole_N 1/day].                               !
!  MorPS0      Small Phytoplankton mortality rate at 0 Celsius,        !
!                [m3/millimole_N 1/day].                               !
!  MorZL0      Large Zooplankton mortality rate at 0 Celsius,          !
!                [m3/millimole_N 1/day].                               !
!  MorZP0      Predator Zooplankton mortality rate at 0 Celsius,       !
!                [m3/millimole_N 1/day].                               !
!  MorZS0      Small Zooplankton mortality rate at 0 Celsius,          !
!                [m3/millimole_N 1/day].                               !
!  Nit0        Nitrification (NH4 to NO3) rate at 0 Celsius, [1/day].  !
!  PARfrac     Fraction of shortwave radiation that is available for   !
!                photosyntesis [nondimensional].                       !
!  PL2ZSstar   Small Zooplankton threshold value for grazing on        !
!                Large Phytoplankton, [millimole_N/m3].                !
!  PL2ZLstar   Large Zooplankton threshold value for grazing on        !
!                Large Phytoplankton, [millimole_N/m3].                !
!  PL2ZPstar   Predator Zooplankton threshold value for grazing on     !
!                Large Phytoplankton, [millimole_N/m3].                !
!  PS2ZLstar   Large Zooplankton threshold value for grazing on        !
!                Small Phytoplankton, [millimole_N/m3].                !
!  PS2ZSstar   Small Zooplankton threshold value for grazing on        !
!                Small Phytoplankton, [millimole_N/m3].                !
!  PusaiL      Large Phytoplankton Ammonium inhibition coefficient,    !
!                [m3/millimole_N].                                     !
!  PusaiPL     Predator Zooplankton grazing on Large Phytoplankton     !
!                inhibition coefficient, [m3/millimole_N].             !
!  PusaiS      Small Phytoplankton Ammonium inhibition coefficient,    !
!                [m3/millimole_N].                                     !
!  PusaiZS     Predator Zooplankton grazing on Small Zooplankton       !
!                inhibition coefficient, [m3/millimole_N].             !
!  ResPL0      Large Phytoplankton respiration rate at 0 Celsius,      !
!                [1/day].                                              !
!  ResPS0      Small Phytoplankton respiration rate at 0 Celsius,      !
!                [1/day].                                              !
!  RSiN        Si:N ratio [millimole_Si/millimole_N].                  !
!  setVOpal    Opal Settling (sinking) velocity [m/day].               !
!  setVPON     PON Settling (sinking) velocity [m/day].                !
!  VD2N0       DON to NH4 decomposition rate at 0 Celsius, [1/day].    !
!  VmaxL       Maximum Large Phytoplankton photosynthetic rate [1/day] !
!                in the absence of photoinhibition under optimal light.!
!  VmaxS       Maximum Small Phytoplankton photosynthetic rate [1/day] !
!                in the absence of photoinhibition under optimal light.!
!  VO2S0       Opal to Silicate decomposition rate at 0 Celsius,       !
!                [1/day].                                              !
!  VP2D0       PON to DON decomposition rate at 0 Celsius, [1/day].    !
!  VP2N0       PON to NH4 decomposition rate at 0 Celsius, [1/day].    !
!  ZL2ZPstar   Small Zooplankton threshold value for grazing on        !
!                Small Phytoplankton, [millimole_N/m3].                !
!  ZS2ZLstar   Large Zooplankton threshold value for grazing on        !
!                Small Zooplankton, [millimole_N/m3].                  !
!  ZS2ZPstar   Predator Zooplankton threshold value for grazing on     !
!                Small Zooplankton, [millimole_N/m3].                  !
!                                                                      !
! Parameters for iron limitation                                       !
!  T_Fe      Iron uptake timescale, [day].                             !
!  A_Fe      Empirical Fe:C power, [-].                                !
!  B_Fe      Empirical Fe:C coefficient, [1/M-C].                      !
!  SK_FeC    Small phytoplankton Fe:C at F=0.5, [muM-Fe/M-C].          !
!  LK_FeC    Large phytoplankton Fe:C at F=0.5, [muM-Fe/M-C].          !
!  FeRR      Fe remineralization rate, [1/day].                        !
!                                                                      !
#ifdef NEMURO_SAN
! parameters for the fish particle model
!  Nfishperyear Number of fish per year class per species
!  Nspecies     Number of fish species
!  Nyearclass   Number of year classes
!  Nfish        Total number of fish
# ifdef PREDATOR
!  Npredspecies      Number of predator species
!  Npredperspecies   Number of fish per year class per species
!  Npred             Total number of predators
# endif
# ifdef FLEET
!  Nboats        Total number of boats
!  Nports        Total number of ports
# endif
#endif
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Set biological tracer identification indices.
!
      integer, allocatable :: idbio(:)  ! Biological tracers
      integer :: iLphy                  ! Large Phytoplankton biomass
      integer :: iSphy                  ! Small Phytoplankton biomass
      integer :: iLzoo                  ! Large Zooplankton biomass
      integer :: iSzoo                  ! Small Zooplankton biomass
      integer :: iPzoo                  ! Predator Zooplankton biomass
      integer :: iNO3_                  ! Nitrate concentration
      integer :: iNH4_                  ! Ammonium concentration
      integer :: iPON_                  ! Particulate Organic Nitrogen
      integer :: iDON_                  ! Dissolved Organic Nitrogen
      integer :: iSiOH                  ! Silicate concentration
      integer :: iopal                  ! Particulate organic silica
# ifdef IRON_LIMIT
      integer :: iFeSp                  ! Small phytoplankton iron
      integer :: iFeLp                  ! Large phytoplankton iron
      integer :: iFeD_                  ! Available disolved iron
# endif
      integer, parameter :: max_species = 5
#  ifdef PREDATOR
      integer, parameter :: max_predspecies = 1
#  endif
#  ifdef FLEET
      integer, parameter :: max_boats = 100
      integer, parameter :: max_ports = 5
#  endif


!
!  Biological parameters.
!
      integer, dimension(Ngrids) :: BioIter

      real(r8), dimension(Ngrids) :: AlphaPL         ! 1/(W/m2) 1/day
      real(r8), dimension(Ngrids) :: AlphaPS         ! 1/(W/m2) 1/day
      real(r8), dimension(Ngrids) :: AlphaZL         ! nondimensional
      real(r8), dimension(Ngrids) :: AlphaZP         ! nondimensional
      real(r8), dimension(Ngrids) :: AlphaZS         ! nondimensional
      real(r8), dimension(Ngrids) :: AttPL           ! m2/mmole_N
      real(r8), dimension(Ngrids) :: AttPS           ! m2/mmole_N
      real(r8), dimension(Ngrids) :: AttSW           ! 1/m
      real(r8), dimension(Ngrids) :: BetaPL          ! 1/(W/m2) 1/day
      real(r8), dimension(Ngrids) :: BetaPS          ! 1/(W/m2) 1/day
      real(r8), dimension(Ngrids) :: BetaZS          ! nondimensional
      real(r8), dimension(Ngrids) :: BetaZL          ! nondimensional
      real(r8), dimension(Ngrids) :: BetaZP          ! nondimensional
      real(r8), dimension(Ngrids) :: GammaL          ! nondimensional
      real(r8), dimension(Ngrids) :: GammaS          ! nondimensional
      real(r8), dimension(Ngrids) :: GRmaxLpl        ! 1/day
      real(r8), dimension(Ngrids) :: GRmaxLps        ! 1/day
      real(r8), dimension(Ngrids) :: GRmaxLzs        ! 1/day
      real(r8), dimension(Ngrids) :: GRmaxPpl        ! 1/day
      real(r8), dimension(Ngrids) :: GRmaxPzl        ! 1/day
      real(r8), dimension(Ngrids) :: GRmaxPzs        ! 1/day
      real(r8), dimension(Ngrids) :: GRmaxSpl        ! 1/day
      real(r8), dimension(Ngrids) :: GRmaxSps        ! 1/day
      real(r8), dimension(Ngrids) :: KD2N            ! 1/Celsius
      real(r8), dimension(Ngrids) :: KGppL           ! 1/Celsius
      real(r8), dimension(Ngrids) :: KGppS           ! 1/Celsius
      real(r8), dimension(Ngrids) :: KGraL           ! 1/Celsius
      real(r8), dimension(Ngrids) :: KGraP           ! 1/Celsius
      real(r8), dimension(Ngrids) :: KGraS           ! 1/Celsius
      real(r8), dimension(Ngrids) :: KMorPL          ! 1/Celsius
      real(r8), dimension(Ngrids) :: KMorPS          ! 1/Celsius
      real(r8), dimension(Ngrids) :: KMorZL          ! 1/Celsius
      real(r8), dimension(Ngrids) :: KMorZP          ! 1/Celsius
      real(r8), dimension(Ngrids) :: KMorZS          ! 1/Celsius
      real(r8), dimension(Ngrids) :: KNH4L           ! mmole_N/m3
      real(r8), dimension(Ngrids) :: KNH4S           ! mmole_N/m3
      real(r8), dimension(Ngrids) :: KNit            ! 1/Celsius
      real(r8), dimension(Ngrids) :: KNO3L           ! mmole_N/m3
      real(r8), dimension(Ngrids) :: KNO3S           ! mmole_N/m3
      real(r8), dimension(Ngrids) :: KO2S            ! 1/Celsius
      real(r8), dimension(Ngrids) :: KP2D            ! 1/Celsius
      real(r8), dimension(Ngrids) :: KP2N            ! 1/Celsius
      real(r8), dimension(Ngrids) :: KPL2ZL          ! mmole_N/m3
      real(r8), dimension(Ngrids) :: KPL2ZS          ! mmole_N/m3
      real(r8), dimension(Ngrids) :: KPS2ZL          ! mmole_N/m3
      real(r8), dimension(Ngrids) :: KPS2ZS          ! mmole_N/m3
      real(r8), dimension(Ngrids) :: KPL2ZP          ! mmole_N/m3
      real(r8), dimension(Ngrids) :: KResPL          ! 1/Celsius
      real(r8), dimension(Ngrids) :: KResPS          ! 1/Celsius
      real(r8), dimension(Ngrids) :: KSiL            ! mmole_Si/m3
      real(r8), dimension(Ngrids) :: KZL2ZP          ! mmole_N/m3
      real(r8), dimension(Ngrids) :: KZS2ZL          ! mmole_N/m3
      real(r8), dimension(Ngrids) :: KZS2ZP          ! mmole_N/m3
      real(r8), dimension(Ngrids) :: LamL            ! m3/mmole_N
      real(r8), dimension(Ngrids) :: LamP            ! m3/mmole_N
      real(r8), dimension(Ngrids) :: LamS            ! m3/mmole_N
      real(r8), dimension(Ngrids) :: MorPL0          ! m3/mmole_N/day
      real(r8), dimension(Ngrids) :: MorPS0          ! m3/mmole_N/day
      real(r8), dimension(Ngrids) :: MorZL0          ! m3/mmole_N 1/day
      real(r8), dimension(Ngrids) :: MorZP0          ! m3/mmole_N 1/day
      real(r8), dimension(Ngrids) :: MorZS0          ! m3/mmole_N 1/day
      real(r8), dimension(Ngrids) :: Nit0            ! 1/day
      real(r8), dimension(Ngrids) :: PARfrac         ! nondimensional
      real(r8), dimension(Ngrids) :: PusaiL          ! m3/mmole_N
      real(r8), dimension(Ngrids) :: PusaiPL         ! m3/mmole_N
      real(r8), dimension(Ngrids) :: PusaiS          ! m3/mmole_N
      real(r8), dimension(Ngrids) :: PusaiZS         ! m3/mmole_N
      real(r8), dimension(Ngrids) :: PL2ZLstar       ! mmole_N/m3
      real(r8), dimension(Ngrids) :: PL2ZPstar       ! mmole_N/m3
      real(r8), dimension(Ngrids) :: PL2ZSstar       ! mmole_N/m3
      real(r8), dimension(Ngrids) :: PS2ZLstar       ! mmole_N/m3
      real(r8), dimension(Ngrids) :: PS2ZSstar       ! mmole_N/m3
      real(r8), dimension(Ngrids) :: ResPL0          ! 1/day
      real(r8), dimension(Ngrids) :: ResPS0          ! 1/day
      real(r8), dimension(Ngrids) :: RSiN            ! mmole_Si/mmole_N
      real(r8), dimension(Ngrids) :: setVOpal        ! m/day
      real(r8), dimension(Ngrids) :: setVPON         ! m/day
      real(r8), dimension(Ngrids) :: VD2N0           ! 1/day
      real(r8), dimension(Ngrids) :: VmaxL           ! 1/day
      real(r8), dimension(Ngrids) :: VmaxS           ! 1/day
      real(r8), dimension(Ngrids) :: VO2S0           ! 1/day
      real(r8), dimension(Ngrids) :: VP2D0           ! 1/day
      real(r8), dimension(Ngrids) :: VP2N0           ! 1/day
      real(r8), dimension(Ngrids) :: ZL2ZPstar       ! mmole_N/m3
      real(r8), dimension(Ngrids) :: ZS2ZLstar       ! mmole_N/m3
      real(r8), dimension(Ngrids) :: ZS2ZPstar       ! mmole_N/m3
#ifdef IRON_LIMIT
      real(r8), dimension(Ngrids) :: T_Fe            ! day
      real(r8), dimension(Ngrids) :: A_Fe            ! nondimensional
      real(r8), dimension(Ngrids) :: B_Fe            ! 1/M-C
      real(r8), dimension(Ngrids) :: SK_FeC          ! muM-Fe/M-C
      real(r8), dimension(Ngrids) :: LK_FeC          ! muM-Fe/M-C
      real(r8), dimension(Ngrids) :: FeRR            ! 1/day
#endif
#ifdef NEMURO_SAN
      integer, dimension(max_species, Ngrids)  :: Nbatch        !  nondimensional
      real(r8), dimension(max_species, Ngrids) :: Fwwt0         !  grams
      real(r8), dimension(max_species, Ngrids) :: Fwth0         !  millions of fish 
      real(r8), dimension(max_species, Ngrids) :: Fage0         !  years
      real(r8), dimension(max_species, Ngrids) :: Fspstr        !  yearday
      real(r8), dimension(max_species, Ngrids) :: Fspend        !  yearday
      integer, dimension(max_species, Ngrids)  :: Fbehave       !  1=Fitness, 2=Kinesis
      real(r8), dimension(max_species, Ngrids) :: ZSpref        !  nondimensional
      real(r8), dimension(max_species, Ngrids) :: ZLpref        !  nondimensional
      real(r8), dimension(max_species, Ngrids) :: ZPpref        !  nondimensional
      real(r8), dimension(max_species, Ngrids) :: CAL_ZF        !
      real(r8), dimension(max_species, Ngrids) :: K_ZS          !   
      real(r8), dimension(max_species, Ngrids) :: K_ZL          !
      real(r8), dimension(max_species, Ngrids) :: K_ZP          !
      real(r8), dimension(max_species, Ngrids) :: a_C           !
      real(r8), dimension(max_species, Ngrids) :: b_C           !
      real(r8), dimension(max_species, Ngrids) :: a_F           !
      real(r8), dimension(max_species, Ngrids) :: a_E           !
      real(r8), dimension(max_species, Ngrids) :: a_R           !
      real(r8), dimension(max_species, Ngrids) :: b_R           !
      real(r8), dimension(max_species, Ngrids) :: d_R           !
      real(r8), dimension(max_species, Ngrids) :: a_A           !
      real(r8), dimension(max_species, Ngrids) :: b_A           !
      real(r8), dimension(max_species, Ngrids) :: a_S           !
      real(r8), dimension(max_species, Ngrids) :: Nymort        !
      real(r8), dimension(max_species, Ngrids) :: Fymort        !
# ifdef PREDATOR
      real(r8), dimension(max_predspecies, Ngrids) :: Pwwt0     !  grams
      real(r8), dimension(max_predspecies, Ngrids) :: Pwth0     !  millions of fish
      real(r8), dimension(max_predspecies, Ngrids) :: K_Fish    !
      real(r8), dimension(max_predspecies, Ngrids) :: Fpref     !
      real(r8), dimension(max_predspecies, Ngrids) :: a_Cmax    !
      real(r8), dimension(max_predspecies, Ngrids) :: b_Cmax    !
      real(r8), dimension(max_predspecies, Ngrids) :: a_Swim    !
      real(r8), dimension(max_predspecies, Ngrids) :: b_Swim    !
# endif
# ifdef FLEET
      integer, dimension(max_ports, Ngrids) :: iPort            ! grid cell
      integer, dimension(max_ports, Ngrids) :: jPort            ! grid cell
      integer, dimension(Ngrids) :: EncMax                      ! nondimensional
      real(r8), dimension(Ngrids) :: CatchMax                   ! kg
      real(r8), dimension(Ngrids) :: TravCost                   ! $
      real(r8), dimension(Ngrids) :: BoatVel                    ! km/h
      real(r8), dimension(Ngrids) :: Qcatch                     ! nondimensional
      real(r8), dimension(Ngrids) :: FishTime                   ! hour
      real(r8), dimension(Ngrids) :: EncRate                    ! nondimensional
      real(r8), dimension(max_ports, Ngrids) :: CatchPrice      ! $/kg
# endif
#endif

#ifdef NEMURO_SAN
      integer, allocatable :: idfish(:) ! Fish species map
      integer, allocatable :: idfish_inv(:) ! Fish species inverse map
# ifdef PREDATOR
      integer, allocatable :: idpred(:) ! Pred species map
      integer, allocatable :: idpred_inv(:) ! Pred species inverse map
# endif
#endif

#ifdef NEMURO_SAN
!
! NOTE: When adding variables, must update NFishV in mod_param.F
!
      integer, dimension(Ngrids) :: Nfishperyear
      integer, dimension(Ngrids) :: Nspecies
      integer, dimension(Ngrids) :: Nyearclass
!   Nfish=Nfishperyear*Nspecies*Nyearclass
      integer, dimension(Ngrids) :: Nfish

# ifdef PREDATOR
      integer, dimension(Ngrids) :: Npredperspecies
      integer, dimension(Ngrids) :: Npredspecies
      integer, dimension(Ngrids) :: Npred
# endif
# ifdef FLEET
      integer, dimension(Ngrids) :: Nboats
      integer, dimension(Ngrids) :: Nports
# endif
#endif

      CONTAINS

      SUBROUTINE initialize_biology
!
!=======================================================================
!                                                                      !
!  This routine sets several variables needed by the biology model.    !
!  It allocates and assigns biological tracers indices.                !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: i, ic
!
!-----------------------------------------------------------------------
!  Set number of biological tracers.
!-----------------------------------------------------------------------
!
#  ifdef IRON_LIMIT
      NBT = 14
#  else
      NBT = 11
#  endif
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
      END IF
#  ifdef NEMURO_SAN
      IF (.not. allocated(idfish))       allocate ( idfish(max_species) )
      IF (.not. allocated(idfish_inv))   allocate ( idfish_inv(max_species) )
#   ifdef PREDATOR
      IF (.not. allocated(idpred))       allocate ( idpred(max_predspecies) )
      IF (.not. allocated(idpred_inv))   allocate ( idpred_inv(max_predspecies) )
#   endif
#  endif
!
!  Set identification indices.
!
      ic=NAT+NPT+NCS+NNS
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      iSphy=ic+1
      iLphy=ic+2
      iSzoo=ic+3
      iLzoo=ic+4
      iPzoo=ic+5
      iNO3_=ic+6
      iNH4_=ic+7
      iPON_=ic+8
      iDON_=ic+9
      iSiOH=ic+10
      iopal=ic+11
# ifdef IRON_LIMIT
      iFeSp=ic+12
      iFeLp=ic+13
      iFeD_=ic+14
# endif

      RETURN
      END SUBROUTINE initialize_biology
