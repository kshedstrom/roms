!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
# ifdef FISHING_FLEET
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
      integer :: iFeD_                  ! Available dissolved iron
# endif
      integer, parameter :: max_species = 5
      integer, parameter :: max_lstages = 6
#  ifdef PREDATOR
      integer, parameter :: max_predspecies = 1
#  endif
#  ifdef FISHING_FLEET
      integer, parameter :: max_boats = 100
      integer, parameter :: max_ports = 5
#  endif


!
!  Biological parameters.
!
      integer, allocatable :: BioIter(:)

      real(r8), allocatable :: AlphaPL(:)            ! 1/(W/m2) 1/day
      real(r8), allocatable :: AlphaPS(:)            ! 1/(W/m2) 1/day
      real(r8), allocatable :: AlphaZL(:)            ! nondimensional
      real(r8), allocatable :: AlphaZP(:)            ! nondimensional
      real(r8), allocatable :: AlphaZS(:)            ! nondimensional
      real(r8), allocatable :: AttPL(:)              ! m2/mmole_N
      real(r8), allocatable :: AttPS(:)              ! m2/mmole_N
      real(r8), allocatable :: AttSW(:)              ! 1/m
      real(r8), allocatable :: BetaPL(:)             ! 1/(W/m2) 1/day
      real(r8), allocatable :: BetaPS(:)             ! 1/(W/m2) 1/day
      real(r8), allocatable :: BetaZS(:)             ! nondimensional
      real(r8), allocatable :: BetaZL(:)             ! nondimensional
      real(r8), allocatable :: BetaZP(:)             ! nondimensional
      real(r8), allocatable :: GammaL(:)             ! nondimensional
      real(r8), allocatable :: GammaS(:)             ! nondimensional
      real(r8), allocatable :: GRmaxLpl(:)           ! 1/day
      real(r8), allocatable :: GRmaxLps(:)           ! 1/day
      real(r8), allocatable :: GRmaxLzs(:)           ! 1/day
      real(r8), allocatable :: GRmaxPpl(:)           ! 1/day
      real(r8), allocatable :: GRmaxPzl(:)           ! 1/day
      real(r8), allocatable :: GRmaxPzs(:)           ! 1/day
      real(r8), allocatable :: GRmaxSpl(:)           ! 1/day
      real(r8), allocatable :: GRmaxSps(:)           ! 1/day
      real(r8), allocatable :: KD2N(:)               ! 1/Celsius
      real(r8), allocatable :: KGppL(:)              ! 1/Celsius
      real(r8), allocatable :: KGppS(:)              ! 1/Celsius
      real(r8), allocatable :: KGraL(:)              ! 1/Celsius
      real(r8), allocatable :: KGraP(:)              ! 1/Celsius
      real(r8), allocatable :: KGraS(:)              ! 1/Celsius
      real(r8), allocatable :: KMorPL(:)             ! 1/Celsius
      real(r8), allocatable :: KMorPS(:)             ! 1/Celsius
      real(r8), allocatable :: KMorZL(:)             ! 1/Celsius
      real(r8), allocatable :: KMorZP(:)             ! 1/Celsius
      real(r8), allocatable :: KMorZS(:)             ! 1/Celsius
      real(r8), allocatable :: KNH4L(:)              ! mmole_N/m3
      real(r8), allocatable :: KNH4S(:)              ! mmole_N/m3
      real(r8), allocatable :: KNit(:)               ! 1/Celsius
      real(r8), allocatable :: KNO3L(:)              ! mmole_N/m3
      real(r8), allocatable :: KNO3S(:)              ! mmole_N/m3
      real(r8), allocatable :: KO2S(:)               ! 1/Celsius
      real(r8), allocatable :: KP2D(:)               ! 1/Celsius
      real(r8), allocatable :: KP2N(:)               ! 1/Celsius
      real(r8), allocatable :: KPL2ZL(:)             ! mmole_N/m3
      real(r8), allocatable :: KPL2ZS(:)             ! mmole_N/m3
      real(r8), allocatable :: KPS2ZL(:)             ! mmole_N/m3
      real(r8), allocatable :: KPS2ZS(:)             ! mmole_N/m3
      real(r8), allocatable :: KPL2ZP(:)             ! mmole_N/m3
      real(r8), allocatable :: KResPL(:)             ! 1/Celsius
      real(r8), allocatable :: KResPS(:)             ! 1/Celsius
      real(r8), allocatable :: KSiL(:)               ! mmole_Si/m3
      real(r8), allocatable :: KZL2ZP(:)             ! mmole_N/m3
      real(r8), allocatable :: KZS2ZL(:)             ! mmole_N/m3
      real(r8), allocatable :: KZS2ZP(:)             ! mmole_N/m3
      real(r8), allocatable :: LamL(:)               ! m3/mmole_N
      real(r8), allocatable :: LamP(:)               ! m3/mmole_N
      real(r8), allocatable :: LamS(:)               ! m3/mmole_N
      real(r8), allocatable :: MorPL0(:)             ! m3/mmole_N/day
      real(r8), allocatable :: MorPS0(:)             ! m3/mmole_N/day
      real(r8), allocatable :: MorZL0(:)             ! m3/mmole_N 1/day
      real(r8), allocatable :: MorZP0(:)             ! m3/mmole_N 1/day
      real(r8), allocatable :: MorZS0(:)             ! m3/mmole_N 1/day
      real(r8), allocatable :: Nit0(:)               ! 1/day
      real(r8), allocatable :: PARfrac(:)            ! nondimensional
      real(r8), allocatable :: PusaiL(:)             ! m3/mmole_N
      real(r8), allocatable :: PusaiPL(:)            ! m3/mmole_N
      real(r8), allocatable :: PusaiS(:)             ! m3/mmole_N
      real(r8), allocatable :: PusaiZS(:)            ! m3/mmole_N
      real(r8), allocatable :: PL2ZLstar(:)          ! mmole_N/m3
      real(r8), allocatable :: PL2ZPstar(:)          ! mmole_N/m3
      real(r8), allocatable :: PL2ZSstar(:)          ! mmole_N/m3
      real(r8), allocatable :: PS2ZLstar(:)          ! mmole_N/m3
      real(r8), allocatable :: PS2ZSstar(:)          ! mmole_N/m3
      real(r8), allocatable :: ResPL0(:)             ! 1/day
      real(r8), allocatable :: ResPS0(:)             ! 1/day
      real(r8), allocatable :: RSiN(:)               ! mmole_Si/mmole_N
      real(r8), allocatable :: setVOpal(:)           ! m/day
      real(r8), allocatable :: setVPON(:)            ! m/day
      real(r8), allocatable :: VD2N0(:)              ! 1/day
      real(r8), allocatable :: VmaxL(:)              ! 1/day
      real(r8), allocatable :: VmaxS(:)              ! 1/day
      real(r8), allocatable :: VO2S0(:)              ! 1/day
      real(r8), allocatable :: VP2D0(:)              ! 1/day
      real(r8), allocatable :: VP2N0(:)              ! 1/day
      real(r8), allocatable :: ZL2ZPstar(:)          ! mmole_N/m3
      real(r8), allocatable :: ZS2ZLstar(:)          ! mmole_N/m3
      real(r8), allocatable :: ZS2ZPstar(:)          ! mmole_N/m3
#ifdef IRON_LIMIT
      real(r8), allocatable :: T_Fe(:)               ! day
      real(r8), allocatable :: A_Fe(:)               ! nondimensional
      real(r8), allocatable :: B_Fe(:)               ! 1/M-C
      real(r8), allocatable :: SK_FeC(:)             ! muM-Fe/M-C
      real(r8), allocatable :: LK_FeC(:)             ! muM-Fe/M-C
      real(r8), allocatable :: FeRR(:)               ! 1/day
#endif
#ifdef NEMURO_SAN
      real(r8), allocatable :: Fwwt0(:,:)         !  grams
      real(r8), allocatable :: Fwth0(:,:)         !  millions of fish
      real(r8), allocatable :: Fspstr(:,:)        !  yearday
      real(r8), allocatable :: Fspend(:,:)        !  yearday
      real(r8), allocatable :: FspTmin(:,:)       !  deg.C
      real(r8), allocatable :: FspTmax(:,:)       !  deg.C
      integer, allocatable  :: Hbehave(:,:)       !  Hor.behave
      integer, allocatable  :: Vbehave(:,:)       !  Vert.behave
      real(r8), allocatable :: ZSpref(:,:,:)
      real(r8), allocatable :: ZSpref_L(:,:)
      real(r8), allocatable :: ZSpref_J(:,:)
      real(r8), allocatable :: ZSpref_A(:,:)
      real(r8), allocatable :: ZLpref(:,:,:)
      real(r8), allocatable :: ZLpref_L(:,:)
      real(r8), allocatable :: ZLpref_J(:,:)
      real(r8), allocatable :: ZLpref_A(:,:)
      real(r8), allocatable :: ZPpref(:,:,:)
      real(r8), allocatable :: ZPpref_L(:,:)
      real(r8), allocatable :: ZPpref_J(:,:)
      real(r8), allocatable :: ZPpref_A(:,:)
      real(r8), allocatable :: K_ZS(:,:,:)
      real(r8), allocatable :: K_ZS_L(:,:)
      real(r8), allocatable :: K_ZS_J(:,:)
      real(r8), allocatable :: K_ZS_A(:,:)
      real(r8), allocatable :: K_ZL(:,:,:)
      real(r8), allocatable :: K_ZL_L(:,:)
      real(r8), allocatable :: K_ZL_J(:,:)
      real(r8), allocatable :: K_ZL_A(:,:)
      real(r8), allocatable :: K_ZP(:,:,:)
      real(r8), allocatable :: K_ZP_L(:,:)
      real(r8), allocatable :: K_ZP_J(:,:)
      real(r8), allocatable :: K_ZP_A(:,:)
      real(r8), allocatable :: Cal_Z(:,:)
      real(r8), allocatable :: Cal_F(:,:)
      real(r8), allocatable :: a_C(:,:,:)
      real(r8), allocatable :: a_C_L(:,:)
      real(r8), allocatable :: a_C_J(:,:)
      real(r8), allocatable :: a_C_A(:,:)
      real(r8), allocatable :: b_C(:,:,:)
      real(r8), allocatable :: b_C_L(:,:)
      real(r8), allocatable :: b_C_J(:,:)
      real(r8), allocatable :: b_C_A(:,:)
      real(r8), allocatable :: pvalmax(:,:,:)
      real(r8), allocatable :: pvalmax_L(:,:)
      real(r8), allocatable :: pvalmax_J(:,:)
      real(r8), allocatable :: pvalmax_A(:,:)
      real(r8), allocatable :: a_R(:,:,:)
      real(r8), allocatable :: a_R_L(:,:)
      real(r8), allocatable :: a_R_J(:,:)
      real(r8), allocatable :: a_R_A(:,:)
      real(r8), allocatable :: b_R(:,:,:)
      real(r8), allocatable :: b_R_L(:,:)
      real(r8), allocatable :: b_R_J(:,:)
      real(r8), allocatable :: b_R_A(:,:)
      real(r8), allocatable :: activity(:,:,:)
      real(r8), allocatable :: activity_L(:,:)
      real(r8), allocatable :: activity_J(:,:)
      real(r8), allocatable :: activity_A(:,:)
      real(r8), allocatable :: d_R(:,:,:)
      real(r8), allocatable :: d_R_L(:,:)
      real(r8), allocatable :: d_R_J(:,:)
      real(r8), allocatable :: d_R_A(:,:)
      real(r8), allocatable :: Fswim(:,:,:)
      real(r8), allocatable :: Fswim_L(:,:)
      real(r8), allocatable :: Fswim_J(:,:)
      real(r8), allocatable :: Fswim_A(:,:)
      real(r8), allocatable :: a_AE(:,:,:)
      real(r8), allocatable :: a_AE_L(:,:)
      real(r8), allocatable :: a_AE_J(:,:)
      real(r8), allocatable :: a_AE_A(:,:)
      real(r8), allocatable :: b_AE(:,:,:)
      real(r8), allocatable :: b_AE_L(:,:)
      real(r8), allocatable :: b_AE_J(:,:)
      real(r8), allocatable :: b_AE_A(:,:)
      real(r8), allocatable :: AEmax(:,:,:)
      real(r8), allocatable :: AEmax_L(:,:)
      real(r8), allocatable :: AEmax_J(:,:)
      real(r8), allocatable :: AEmax_A(:,:)
      real(r8), allocatable :: te1(:,:,:)
      real(r8), allocatable :: te1_L(:,:)
      real(r8), allocatable :: te1_J(:,:)
      real(r8), allocatable :: te1_A(:,:)
      real(r8), allocatable :: te2(:,:,:)
      real(r8), allocatable :: te2_L(:,:)
      real(r8), allocatable :: te2_J(:,:)
      real(r8), allocatable :: te2_A(:,:)
      real(r8), allocatable :: te3(:,:,:)
      real(r8), allocatable :: te3_L(:,:)
      real(r8), allocatable :: te3_J(:,:)
      real(r8), allocatable :: te3_A(:,:)
      real(r8), allocatable :: te4(:,:,:)
      real(r8), allocatable :: te4_L(:,:)
      real(r8), allocatable :: te4_J(:,:)
      real(r8), allocatable :: te4_A(:,:)
      real(r8), allocatable :: xk1(:,:,:)
      real(r8), allocatable :: xk1_L(:,:)
      real(r8), allocatable :: xk1_J(:,:)
      real(r8), allocatable :: xk1_A(:,:)
      real(r8), allocatable :: xk2(:,:,:)
      real(r8), allocatable :: xk2_L(:,:)
      real(r8), allocatable :: xk2_J(:,:)
      real(r8), allocatable :: xk2_A(:,:)
      real(r8), allocatable :: xk3(:,:,:)
      real(r8), allocatable :: xk3_L(:,:)
      real(r8), allocatable :: xk3_J(:,:)
      real(r8), allocatable :: xk3_A(:,:)
      real(r8), allocatable :: xk4(:,:,:)
      real(r8), allocatable :: xk4_L(:,:)
      real(r8), allocatable :: xk4_J(:,:)
      real(r8), allocatable :: xk4_A(:,:)
      real(r8), allocatable :: cr(:,:,:)
      real(r8), allocatable :: cr_L(:,:)
      real(r8), allocatable :: cr_J(:,:)
      real(r8), allocatable :: cr_A(:,:)
      real(r8), allocatable :: tr(:,:,:)
      real(r8), allocatable :: tr_L(:,:)
      real(r8), allocatable :: tr_J(:,:)
      real(r8), allocatable :: tr_A(:,:)
      real(r8), allocatable :: Wffeed(:,:)
      real(r8), allocatable :: Lffeed(:,:)
      real(r8), allocatable :: WeightLJ(:,:)
      real(r8), allocatable :: LengthLJ(:,:)
      real(r8), allocatable :: WeightJA(:,:)
      real(r8), allocatable :: LengthJA(:,:)
      real(r8), allocatable :: aw2l(:,:,:)
      real(r8), allocatable :: aw2l_L(:,:)
      real(r8), allocatable :: aw2l_J(:,:)
      real(r8), allocatable :: aw2l_A(:,:)
      real(r8), allocatable :: bw2l(:,:,:)
      real(r8), allocatable :: bw2l_L(:,:)
      real(r8), allocatable :: bw2l_J(:,:)
      real(r8), allocatable :: bw2l_A(:,:)
      real(r8), allocatable :: dSLk(:,:,:)
      real(r8), allocatable :: dSLk_L(:,:)
      real(r8), allocatable :: dSLk_J(:,:)
      real(r8), allocatable :: dSLk_A(:,:)
      real(r8), allocatable :: dSLinf(:,:,:)
      real(r8), allocatable :: dSLinf_L(:,:)
      real(r8), allocatable :: dSLinf_J(:,:)
      real(r8), allocatable :: dSLinf_A(:,:)
      real(r8), allocatable :: al2w(:,:,:)
      real(r8), allocatable :: al2w_L(:,:)
      real(r8), allocatable :: al2w_J(:,:)
      real(r8), allocatable :: al2w_A(:,:)
      real(r8), allocatable :: bl2w(:,:,:)
      real(r8), allocatable :: bl2w_L(:,:)
      real(r8), allocatable :: bl2w_J(:,:)
      real(r8), allocatable :: bl2w_A(:,:)
      real(r8), allocatable :: abatch(:,:)
      real(r8), allocatable :: bbatch(:,:)
      real(r8), allocatable :: T0batch(:,:)
      real(r8), allocatable :: apof(:,:)
      real(r8), allocatable :: bpof(:,:)
      real(r8), allocatable :: T0pof(:,:)
      real(r8), allocatable :: epg(:,:)
      real(r8), allocatable :: eegg(:,:)
      real(r8), allocatable :: megg(:,:)
      real(r8), allocatable :: breed(:,:)
      real(r8), allocatable :: amature(:,:)
      real(r8), allocatable :: bmature(:,:)
      real(r8), allocatable :: pctxwt(:,:)
      real(r8), allocatable :: pctgain(:,:)
      real(r8), allocatable :: Nmort(:,:,:)
      real(r8), allocatable :: Nmort_E(:,:)
      real(r8), allocatable :: Nmort_Y(:,:)
      real(r8), allocatable :: Nmort_L(:,:)
      real(r8), allocatable :: Nmort_J(:,:)
      real(r8), allocatable :: Nmort_A(:,:)
!      real(r8), allocatable :: Nymort(:,:)
      real(r8), allocatable :: Fymort(:,:)
      real(r8), allocatable :: DDmort1_J(:,:)
      real(r8), allocatable :: DDmort2_J(:,:)
      real(r8), allocatable :: DDmort3_J(:,:)
      real(r8), allocatable :: DDmort4_J(:,:)
      real(r8), allocatable :: DDscale_J(:,:)
# ifdef PREDATOR
      real(r8), allocatable :: Pwwt0(:,:)     ! grams
      real(r8), allocatable :: Pwth0(:,:)     ! worth
      real(r8), allocatable :: Pswim(:,:)     ! cm s-1
      real(r8), allocatable :: HIF(:,:)       ! HIF
      real(r8), allocatable :: FEUE1(:,:)     ! FEUE1
      real(r8), allocatable :: FEUE2(:,:)     ! FEUE2
      real(r8), allocatable :: FEUE3(:,:)     ! FEUE3
      real(r8), allocatable :: FEUE4(:,:)     ! FEUE4
      real(r8), allocatable :: PED1(:,:)      ! PED1
      real(r8), allocatable :: PED2(:,:)      ! PED2
      real(r8), allocatable :: PED3(:,:)      ! PED3
      real(r8), allocatable :: PED4(:,:)      ! PED4
      real(r8), allocatable :: Eprot(:,:)     ! Eprot
      real(r8), allocatable :: Efat(:,:)      ! Efat
      real(r8), allocatable :: RefpFat(:,:)   ! RefpFat
      real(r8), allocatable :: RefpPr(:,:)    ! RefpPr
      real(r8), allocatable :: Linf(:,:)      ! Linf
      real(r8), allocatable :: vbk(:,:)       ! vbk
      real(r8), allocatable :: vbt0(:,:)      ! vbt0
      real(r8), allocatable :: al2m(:,:)      ! al2m
      real(r8), allocatable :: bl2m(:,:)      ! bl2m
      real(r8), allocatable :: Pmgstr(:,:)    ! day
      real(r8), allocatable :: Pmgend(:,:)    ! day
      real(r8), allocatable :: Fpref(:,:)
      real(r8), allocatable :: K_Fish(:,:)
      real(r8), allocatable :: Pcmax(:,:)     ! g g-1 day-1
# endif
# ifdef FISHING_FLEET
      integer, allocatable :: iPort(:,:)           ! grid cell
      integer, allocatable :: jPort(:,:)           ! grid cell
      integer, allocatable  :: EncMax(:)           ! nondimensional
      real(r8), allocatable :: CatchMax(:)         ! kg
      real(r8), allocatable :: TravCost(:)         ! $
      real(r8), allocatable :: BoatVel(:)          ! km/h
      real(r8), allocatable :: Qcatch(:)           ! nondimensional
      real(r8), allocatable :: FishTime(:)         ! hour
      real(r8), allocatable :: EncRate(:)          ! nondimensional
      real(r8), allocatable :: CatchPrice(:,:)     ! $/kg
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
      integer, allocatable :: Nfishperyear(:)
      integer, allocatable :: Nspecies(:)
      integer, allocatable :: Nyearclass(:)
!   Nfish=Nfishperyear*Nspecies*Nyearclass
      integer, allocatable :: Nfish(:)

# ifdef PREDATOR
      integer, allocatable :: Npredperspecies(:)
      integer, allocatable :: Npredspecies(:)
      integer, allocatable :: Npred(:)
# endif
# ifdef FISHING_FLEET
      integer, allocatable :: Nboats(:)
      integer, allocatable :: Nports(:)
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
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(AlphaPL)) THEN
        allocate ( AlphaPL(Ngrids) )
      END IF
      IF (.not.allocated(AlphaPS)) THEN
        allocate ( AlphaPS(Ngrids) )
      END IF
      IF (.not.allocated(AlphaZL)) THEN
        allocate ( AlphaZL(Ngrids) )
      END IF
      IF (.not.allocated(AlphaZP)) THEN
        allocate ( AlphaZP(Ngrids) )
      END IF
      IF (.not.allocated(AlphaZS)) THEN
        allocate ( AlphaZS(Ngrids) )
      END IF
      IF (.not.allocated(AttPL)) THEN
        allocate ( AttPL(Ngrids) )
      END IF
      IF (.not.allocated(AttPS)) THEN
        allocate ( AttPS(Ngrids) )
      END IF
      IF (.not.allocated(AttSW)) THEN
        allocate ( AttSW(Ngrids) )
      END IF
      IF (.not.allocated(BetaPL)) THEN
        allocate ( BetaPL(Ngrids) )
      END IF
      IF (.not.allocated(BetaPS)) THEN
        allocate ( BetaPS(Ngrids) )
      END IF
      IF (.not.allocated(BetaZS)) THEN
        allocate ( BetaZS(Ngrids) )
      END IF
      IF (.not.allocated(BetaZL)) THEN
        allocate ( BetaZL(Ngrids) )
      END IF
      IF (.not.allocated(BetaZP)) THEN
        allocate ( BetaZP(Ngrids) )
      END IF
      IF (.not.allocated(GammaL)) THEN
        allocate ( GammaL(Ngrids) )
      END IF
      IF (.not.allocated(GammaS)) THEN
        allocate ( GammaS(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxLpl)) THEN
        allocate ( GRmaxLpl(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxLps)) THEN
        allocate ( GRmaxLps(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxLzs)) THEN
        allocate ( GRmaxLzs(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxPpl)) THEN
        allocate ( GRmaxPpl(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxPzl)) THEN
        allocate ( GRmaxPzl(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxPzs)) THEN
        allocate ( GRmaxPzs(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxSpl)) THEN
        allocate ( GRmaxSpl(Ngrids) )
      END IF
      IF (.not.allocated(GRmaxSps)) THEN
        allocate ( GRmaxSps(Ngrids) )
      END IF
      IF (.not.allocated(KD2N)) THEN
        allocate ( KD2N(Ngrids) )
      END IF
      IF (.not.allocated(KGppL)) THEN
        allocate ( KGppL(Ngrids) )
      END IF
      IF (.not.allocated(KGppS)) THEN
        allocate ( KGppS(Ngrids) )
      END IF
      IF (.not.allocated(KGraL)) THEN
        allocate ( KGraL(Ngrids) )
      END IF
      IF (.not.allocated(KGraP)) THEN
        allocate ( KGraP(Ngrids) )
      END IF
      IF (.not.allocated(KGraS)) THEN
        allocate ( KGraS(Ngrids) )
      END IF
      IF (.not.allocated(KMorPL)) THEN
        allocate ( KMorPL(Ngrids) )
      END IF
      IF (.not.allocated(KMorPS)) THEN
        allocate ( KMorPS(Ngrids) )
      END IF
      IF (.not.allocated(KMorZL)) THEN
        allocate ( KMorZL(Ngrids) )
      END IF
      IF (.not.allocated(KMorZP)) THEN
        allocate ( KMorZP(Ngrids) )
      END IF
      IF (.not.allocated(KMorZS)) THEN
        allocate ( KMorZS(Ngrids) )
      END IF
      IF (.not.allocated(KNH4L)) THEN
        allocate ( KNH4L(Ngrids) )
      END IF
      IF (.not.allocated(KNH4S)) THEN
        allocate ( KNH4S(Ngrids) )
      END IF
      IF (.not.allocated(KNit)) THEN
        allocate ( KNit(Ngrids) )
      END IF
      IF (.not.allocated(KNO3L)) THEN
        allocate ( KNO3L(Ngrids) )
      END IF
      IF (.not.allocated(KNO3S)) THEN
        allocate ( KNO3S(Ngrids) )
      END IF
      IF (.not.allocated(KO2S)) THEN
        allocate ( KO2S(Ngrids) )
      END IF
      IF (.not.allocated(KP2D)) THEN
        allocate ( KP2D(Ngrids) )
      END IF
      IF (.not.allocated(KP2N)) THEN
        allocate ( KP2N(Ngrids) )
      END IF
      IF (.not.allocated(KPL2ZL)) THEN
        allocate ( KPL2ZL(Ngrids) )
      END IF
      IF (.not.allocated(KPL2ZS)) THEN
        allocate ( KPL2ZS(Ngrids) )
      END IF
      IF (.not.allocated(KPS2ZL)) THEN
        allocate ( KPS2ZL(Ngrids) )
      END IF
      IF (.not.allocated(KPS2ZS)) THEN
        allocate ( KPS2ZS(Ngrids) )
      END IF
      IF (.not.allocated(KPL2ZP)) THEN
        allocate ( KPL2ZP(Ngrids) )
      END IF
      IF (.not.allocated(KResPL)) THEN
        allocate ( KResPL(Ngrids) )
      END IF
      IF (.not.allocated(KResPS)) THEN
        allocate ( KResPS(Ngrids) )
      END IF
      IF (.not.allocated(KSiL)) THEN
        allocate ( KSiL(Ngrids) )
      END IF
      IF (.not.allocated(KZL2ZP)) THEN
        allocate ( KZL2ZP(Ngrids) )
      END IF
      IF (.not.allocated(KZS2ZL)) THEN
        allocate ( KZS2ZL(Ngrids) )
      END IF
      IF (.not.allocated(KZS2ZP)) THEN
        allocate ( KZS2ZP(Ngrids) )
      END IF
      IF (.not.allocated(LamL)) THEN
        allocate ( LamL(Ngrids) )
      END IF
      IF (.not.allocated(LamP)) THEN
        allocate ( LamP(Ngrids) )
      END IF
      IF (.not.allocated(LamS)) THEN
        allocate ( LamS(Ngrids) )
      END IF
      IF (.not.allocated(MorPL0)) THEN
        allocate ( MorPL0(Ngrids) )
      END IF
      IF (.not.allocated(MorPS0)) THEN
        allocate ( MorPS0(Ngrids) )
      END IF
      IF (.not.allocated(MorZL0)) THEN
        allocate ( MorZL0(Ngrids) )
      END IF
      IF (.not.allocated(MorZP0)) THEN
        allocate ( MorZP0(Ngrids) )
      END IF
      IF (.not.allocated(MorZS0)) THEN
        allocate ( MorZS0(Ngrids) )
      END IF
      IF (.not.allocated(Nit0)) THEN
        allocate ( Nit0(Ngrids) )
      END IF
      IF (.not.allocated(PARfrac)) THEN
        allocate ( PARfrac(Ngrids) )
      END IF
      IF (.not.allocated(PusaiL)) THEN
        allocate ( PusaiL(Ngrids) )
      END IF
      IF (.not.allocated(PusaiPL)) THEN
        allocate ( PusaiPL(Ngrids) )
      END IF
      IF (.not.allocated(PusaiS)) THEN
        allocate ( PusaiS(Ngrids) )
      END IF
      IF (.not.allocated(PusaiZS)) THEN
        allocate ( PusaiZS(Ngrids) )
      END IF
      IF (.not.allocated(PL2ZLstar)) THEN
        allocate ( PL2ZLstar(Ngrids) )
      END IF
      IF (.not.allocated(PL2ZPstar)) THEN
        allocate ( PL2ZPstar(Ngrids) )
      END IF
      IF (.not.allocated(PL2ZSstar)) THEN
        allocate ( PL2ZSstar(Ngrids) )
      END IF
      IF (.not.allocated(PS2ZLstar)) THEN
        allocate ( PS2ZLstar(Ngrids) )
      END IF
      IF (.not.allocated(PS2ZSstar)) THEN
        allocate ( PS2ZSstar(Ngrids) )
      END IF
      IF (.not.allocated(ResPL0)) THEN
        allocate ( ResPL0(Ngrids) )
      END IF
      IF (.not.allocated(ResPS0)) THEN
        allocate ( ResPS0(Ngrids) )
      END IF
      IF (.not.allocated(RSiN)) THEN
        allocate ( RSiN(Ngrids) )
      END IF
      IF (.not.allocated(setVOpal)) THEN
        allocate ( setVOpal(Ngrids) )
      END IF
      IF (.not.allocated(setVPON)) THEN
        allocate ( setVPON(Ngrids) )
      END IF
      IF (.not.allocated(VD2N0)) THEN
        allocate ( VD2N0(Ngrids) )
      END IF
      IF (.not.allocated(VmaxL)) THEN
        allocate ( VmaxL(Ngrids) )
      END IF
      IF (.not.allocated(VmaxS)) THEN
        allocate ( VmaxS(Ngrids) )
      END IF
      IF (.not.allocated(VO2S0)) THEN
        allocate ( VO2S0(Ngrids) )
      END IF
      IF (.not.allocated(VP2D0)) THEN
        allocate ( VP2D0(Ngrids) )
      END IF
      IF (.not.allocated(VP2N0)) THEN
        allocate ( VP2N0(Ngrids) )
      END IF
      IF (.not.allocated(ZL2ZPstar)) THEN
        allocate ( ZL2ZPstar(Ngrids) )
      END IF
      IF (.not.allocated(ZS2ZLstar)) THEN
        allocate ( ZS2ZLstar(Ngrids) )
      END IF
      IF (.not.allocated(ZS2ZPstar)) THEN
        allocate ( ZS2ZPstar(Ngrids) )
      END IF
#ifdef IRON_LIMIT
      IF (.not.allocated(T_Fe)) THEN
        allocate ( T_Fe(Ngrids) )
      END IF
      IF (.not.allocated(A_Fe)) THEN
        allocate ( A_Fe(Ngrids) )
      END IF
      IF (.not.allocated(B_Fe)) THEN
        allocate ( B_Fe(Ngrids) )
      END IF
      IF (.not.allocated(SK_FeC)) THEN
        allocate ( SK_FeC(Ngrids) )
      END IF
      IF (.not.allocated(LK_FeC)) THEN
        allocate ( LK_FeC(Ngrids) )
      END IF
      IF (.not.allocated(FeRR)) THEN
        allocate ( FeRR(Ngrids) )
      END IF
#endif
#ifdef NEMURO_SAN
      IF (.not.allocated(Fwwt0)) THEN
        allocate( Fwwt0(max_species, Ngrids))
      END IF
      IF (.not.allocated(Fwth0)) THEN
        allocate( Fwth0(max_species, Ngrids))
      END IF
      IF (.not.allocated(Fspstr)) THEN
        allocate( Fspstr(max_species, Ngrids))
      END IF
      IF (.not.allocated(Fspend)) THEN
        allocate( Fspend(max_species, Ngrids))
      END IF
      IF (.not.allocated(FspTmin)) THEN
        allocate( FspTmin(max_species, Ngrids))
      END IF
      IF (.not.allocated(FspTmax)) THEN
        allocate( FspTmax(max_species, Ngrids))
      END IF
      IF (.not.allocated(Hbehave)) THEN
        allocate( Hbehave(max_species, Ngrids))
      END IF
      IF (.not.allocated(Vbehave)) THEN
        allocate( Vbehave(max_species, Ngrids))
      END IF
      IF (.not.allocated(ZSpref)) THEN
        allocate( ZSpref(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(ZSpref_L)) THEN
        allocate( ZSpref_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(ZSpref_J)) THEN
        allocate( ZSpref_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(ZSpref_A)) THEN
        allocate( ZSpref_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(ZLpref)) THEN
        allocate( ZLpref(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(ZLpref_L)) THEN
        allocate( ZLpref_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(ZLpref_J)) THEN
        allocate( ZLpref_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(ZLpref_A)) THEN
        allocate( ZLpref_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(ZPpref)) THEN
        allocate( ZPpref(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(ZPpref_L)) THEN
        allocate( ZPpref_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(ZPpref_J)) THEN
        allocate( ZPpref_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(ZPpref_A)) THEN
        allocate( ZPpref_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(K_ZS)) THEN
        allocate( K_ZS(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(K_ZS_L)) THEN
        allocate( K_ZS_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(K_ZS_J)) THEN
        allocate( K_ZS_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(K_ZS_A)) THEN
        allocate( K_ZS_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(K_ZL)) THEN
        allocate( K_ZL(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(K_ZL_L)) THEN
        allocate( K_ZL_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(K_ZL_J)) THEN
        allocate( K_ZL_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(K_ZL_A)) THEN
        allocate( K_ZL_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(K_ZP)) THEN
        allocate( K_ZP(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(K_ZP_L)) THEN
        allocate( K_ZP_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(K_ZP_J)) THEN
        allocate( K_ZP_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(K_ZP_A)) THEN
        allocate( K_ZP_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(Cal_Z)) THEN
        allocate( Cal_Z(max_species, Ngrids))
      END IF
      IF (.not.allocated(Cal_F)) THEN
        allocate( Cal_F(max_species, Ngrids))
      END IF
      IF (.not.allocated(a_C)) THEN
        allocate( a_C(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(a_C_L)) THEN
        allocate( a_C_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(a_C_J)) THEN
        allocate( a_C_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(a_C_A)) THEN
        allocate( a_C_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(b_C)) THEN
        allocate( b_C(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(b_C_L)) THEN
        allocate( b_C_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(b_C_J)) THEN
        allocate( b_C_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(b_C_A)) THEN
        allocate( b_C_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(pvalmax)) THEN
        allocate( pvalmax(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(pvalmax_L)) THEN
        allocate( pvalmax_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(pvalmax_J)) THEN
        allocate( pvalmax_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(pvalmax_A)) THEN
        allocate( pvalmax_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(a_R)) THEN
        allocate( a_R(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(a_R_L)) THEN
        allocate( a_R_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(a_R_J)) THEN
        allocate( a_R_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(a_R_A)) THEN
        allocate( a_R_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(b_R)) THEN
        allocate( b_R(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(b_R_L)) THEN
        allocate( b_R_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(b_R_J)) THEN
        allocate( b_R_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(b_R_A)) THEN
        allocate( b_R_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(activity)) THEN
        allocate( activity(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(activity_L)) THEN
        allocate( activity_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(activity_J)) THEN
        allocate( activity_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(activity_A)) THEN
        allocate( activity_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(d_R)) THEN
        allocate( d_R(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(d_R_L)) THEN
        allocate( d_R_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(d_R_J)) THEN
        allocate( d_R_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(d_R_A)) THEN
        allocate( d_R_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(Fswim)) THEN
        allocate( Fswim(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(Fswim_L)) THEN
        allocate( Fswim_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(Fswim_J)) THEN
        allocate( Fswim_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(Fswim_A)) THEN
        allocate( Fswim_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(a_AE)) THEN
        allocate( a_AE(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(a_AE_L)) THEN
        allocate( a_AE_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(a_AE_J)) THEN
        allocate( a_AE_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(a_AE_A)) THEN
        allocate( a_AE_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(b_AE)) THEN
        allocate( b_AE(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(b_AE_L)) THEN
        allocate( b_AE_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(b_AE_J)) THEN
        allocate( b_AE_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(b_AE_A)) THEN
        allocate( b_AE_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(AEmax)) THEN
        allocate( AEmax(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(AEmax_L)) THEN
        allocate( AEmax_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(AEmax_J)) THEN
        allocate( AEmax_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(AEmax_A)) THEN
        allocate( AEmax_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(te1)) THEN
        allocate( te1(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(te1_L)) THEN
        allocate( te1_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(te1_J)) THEN
        allocate( te1_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(te1_A)) THEN
        allocate( te1_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(te2)) THEN
        allocate( te2(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(te2_L)) THEN
        allocate( te2_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(te2_J)) THEN
        allocate( te2_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(te2_A)) THEN
        allocate( te2_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(te3)) THEN
        allocate( te3(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(te3_L)) THEN
        allocate( te3_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(te3_J)) THEN
        allocate( te3_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(te3_A)) THEN
        allocate( te3_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(te4)) THEN
        allocate( te4(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(te4_L)) THEN
        allocate( te4_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(te4_J)) THEN
        allocate( te4_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(te4_A)) THEN
        allocate( te4_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(xk1)) THEN
        allocate( xk1(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(xk1_L)) THEN
        allocate( xk1_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(xk1_J)) THEN
        allocate( xk1_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(xk1_A)) THEN
        allocate( xk1_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(xk2)) THEN
        allocate( xk2(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(xk2_L)) THEN
        allocate( xk2_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(xk2_J)) THEN
        allocate( xk2_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(xk2_A)) THEN
        allocate( xk2_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(xk3)) THEN
        allocate( xk3(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(xk3_L)) THEN
        allocate( xk3_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(xk3_J)) THEN
        allocate( xk3_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(xk3_A)) THEN
        allocate( xk3_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(xk4)) THEN
        allocate( xk4(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(xk4_L)) THEN
        allocate( xk4_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(xk4_J)) THEN
        allocate( xk4_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(xk4_A)) THEN
        allocate( xk4_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(cr)) THEN
        allocate( cr(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(cr_L)) THEN
        allocate( cr_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(cr_J)) THEN
        allocate( cr_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(cr_A)) THEN
        allocate( cr_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(tr)) THEN
        allocate( tr(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(tr_L)) THEN
        allocate( tr_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(tr_J)) THEN
        allocate( tr_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(tr_A)) THEN
        allocate( tr_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(Wffeed)) THEN
        allocate( Wffeed(max_species, Ngrids))
      END IF
      IF (.not.allocated(Lffeed)) THEN
        allocate( Lffeed(max_species, Ngrids))
      END IF
      IF (.not.allocated(WeightLJ)) THEN
        allocate( WeightLJ(max_species, Ngrids))
      END IF
      IF (.not.allocated(LengthLJ)) THEN
        allocate( LengthLJ(max_species, Ngrids))
      END IF
      IF (.not.allocated(WeightJA)) THEN
        allocate( WeightJA(max_species, Ngrids))
      END IF
      IF (.not.allocated(LengthJA)) THEN
        allocate( LengthJA(max_species, Ngrids))
      END IF
      IF (.not.allocated(aw2l)) THEN
        allocate( aw2l(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(aw2l_L)) THEN
        allocate( aw2l_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(aw2l_J)) THEN
        allocate( aw2l_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(aw2l_A)) THEN
        allocate( aw2l_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(bw2l)) THEN
        allocate( bw2l(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(bw2l_L)) THEN
        allocate( bw2l_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(bw2l_J)) THEN
        allocate( bw2l_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(bw2l_A)) THEN
        allocate( bw2l_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(dSLk)) THEN
        allocate( dSLk(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(dSLk_L)) THEN
        allocate( dSLk_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(dSLk_J)) THEN
        allocate( dSLk_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(dSLk_A)) THEN
        allocate( dSLk_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(dSLinf)) THEN
        allocate( dSLinf(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(dSLinf_L)) THEN
        allocate( dSLinf_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(dSLinf_J)) THEN
        allocate( dSLinf_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(dSLinf_A)) THEN
        allocate( dSLinf_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(al2w)) THEN
        allocate( al2w(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(al2w_L)) THEN
        allocate( al2w_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(al2w_J)) THEN
        allocate( al2w_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(al2w_A)) THEN
        allocate( al2w_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(bl2w)) THEN
        allocate( bl2w(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(bl2w_L)) THEN
        allocate( bl2w_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(bl2w_J)) THEN
        allocate( bl2w_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(bl2w_A)) THEN
        allocate( bl2w_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(abatch)) THEN
        allocate( abatch(max_species, Ngrids))
      END IF
      IF (.not.allocated(bbatch)) THEN
        allocate( bbatch(max_species, Ngrids))
      END IF
      IF (.not.allocated(T0batch)) THEN
        allocate( T0batch(max_species, Ngrids))
      END IF
      IF (.not.allocated(apof)) THEN
        allocate( apof(max_species, Ngrids))
      END IF
      IF (.not.allocated(bpof)) THEN
        allocate( bpof(max_species, Ngrids))
      END IF
      IF (.not.allocated(T0pof)) THEN
        allocate( T0pof(max_species, Ngrids))
      END IF
      IF (.not.allocated(epg)) THEN
        allocate( epg(max_species, Ngrids))
      END IF
      IF (.not.allocated(eegg)) THEN
        allocate( eegg(max_species, Ngrids))
      END IF
      IF (.not.allocated(megg)) THEN
        allocate( megg(max_species, Ngrids))
      END IF
      IF (.not.allocated(breed)) THEN
        allocate( breed(max_species, Ngrids))
      END IF
      IF (.not.allocated(amature)) THEN
        allocate( amature(max_species, Ngrids))
      END IF
      IF (.not.allocated(bmature)) THEN
        allocate( bmature(max_species, Ngrids))
      END IF
      IF (.not.allocated(pctxwt)) THEN
        allocate( pctxwt(max_species, Ngrids))
      END IF
      IF (.not.allocated(pctgain)) THEN
        allocate( pctgain(max_species, Ngrids))
      END IF
      IF (.not.allocated(Nmort)) THEN
        allocate( Nmort(max_lstages, max_species, Ngrids))
      END IF
      IF (.not.allocated(Nmort_E)) THEN
        allocate( Nmort_E(max_species, Ngrids))
      END IF
      IF (.not.allocated(Nmort_Y)) THEN
        allocate( Nmort_Y(max_species, Ngrids))
      END IF
      IF (.not.allocated(Nmort_L)) THEN
        allocate( Nmort_L(max_species, Ngrids))
      END IF
      IF (.not.allocated(Nmort_J)) THEN
        allocate( Nmort_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(Nmort_A)) THEN
        allocate( Nmort_A(max_species, Ngrids))
      END IF
      IF (.not.allocated(Fymort)) THEN
        allocate( Fymort(max_species, Ngrids))
      END IF
      IF (.not.allocated(DDmort1_J)) THEN
        allocate( DDmort1_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(DDmort2_J)) THEN
        allocate( DDmort2_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(DDmort3_J)) THEN
        allocate( DDmort3_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(DDmort4_J)) THEN
        allocate( DDmort4_J(max_species, Ngrids))
      END IF
      IF (.not.allocated(DDscale_J)) THEN
        allocate( DDscale_J(max_species, Ngrids))
      END IF
# ifdef PREDATOR
      IF (.not.allocated(Pwwt0)) THEN
        allocate( Pwwt0(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(Pwth0)) THEN
        allocate( Pwth0(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(Pswim)) THEN
        allocate( Pswim(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(HIF)) THEN
        allocate( HIF(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(FEUE1)) THEN
        allocate( FEUE1(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(FEUE2)) THEN
        allocate( FEUE2(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(FEUE3)) THEN
        allocate( FEUE3(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(FEUE4)) THEN
        allocate( FEUE4(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(PED1)) THEN
        allocate( PED1(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(PED2)) THEN
        allocate( PED2(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(PED3)) THEN
        allocate( PED3(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(PED4)) THEN
        allocate( PED4(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(Eprot)) THEN
        allocate( Eprot(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(Efat)) THEN
        allocate( Efat(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(RefpFat)) THEN
        allocate( RefpFat(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(RefpPr)) THEN
        allocate( RefpPr(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(Linf)) THEN
        allocate( Linf(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(vbk)) THEN
        allocate( vbk(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(vbt0)) THEN
        allocate( vbt0(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(al2m)) THEN
        allocate( al2m(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(bl2m)) THEN
        allocate( bl2m(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(Pmgstr)) THEN
        allocate( Pmgstr(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(Pmgend)) THEN
        allocate( Pmgend(max_predspecies, Ngrids))
      END IF
      IF (.not.allocated(K_Fish)) THEN
        allocate( K_Fish(max_species, Ngrids))
      END IF
      IF (.not.allocated(Fpref)) THEN
        allocate( Fpref(max_species, Ngrids))
      END IF
      IF (.not.allocated(Pcmax)) THEN
        allocate( Pcmax(max_predspecies, Ngrids))
      END IF
# endif
# ifdef FISHING_FLEET
      IF (.not.allocated(iPort)) THEN
        allocate( iPort(max_ports, Ngrids))
      END IF
      IF (.not.allocated(jPort)) THEN
        allocate( jPort(max_ports, Ngrids))
      END IF
      IF (.not.allocated(EncMax)) THEN
        allocate( EncMax(Ngrids))
      END IF
      IF (.not.allocated(CatchMax)) THEN
        allocate( CatchMax(Ngrids))
      END IF
      IF (.not.allocated(TravCost)) THEN
        allocate( TravCost(Ngrids))
      END IF
      IF (.not.allocated(BoatVel)) THEN
        allocate( BoatVel(Ngrids))
      END IF
      IF (.not.allocated(Qcatch)) THEN
        allocate( Qcatch(Ngrids))
      END IF
      IF (.not.allocated(FishTime)) THEN
        allocate( FishTime(Ngrids))
      END IF
      IF (.not.allocated(EncRate)) THEN
        allocate( EncRate(Ngrids))
      END IF
      IF (.not.allocated(CatchPrice)) THEN
        allocate( CatchPrice(max_ports, Ngrids))
      END IF
# endif
#endif
!-----------------------------------------------------------------------
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
      END IF
#  ifdef NEMURO_SAN
      IF (.not. allocated(idfish))                                      &
     &                   allocate ( idfish(max_species) )
      IF (.not. allocated(idfish_inv))                                  &
     &                   allocate ( idfish_inv(max_species) )
#   ifdef PREDATOR
      IF (.not. allocated(idpred))                                      &
     &                   allocate ( idpred(max_predspecies) )
      IF (.not. allocated(idpred_inv))                                  &
     &                   allocate ( idpred_inv(max_predspecies) )
#   endif
#  endif
#ifdef NEMURO_SAN
!
! NOTE: When adding variables, must update NFishV in mod_param.F
!
      IF (.not. allocated(Nfishperyear))                                &
     &                   allocate(Nfishperyear(Ngrids))
      IF (.not. allocated(Nspecies))                                    &
     &                   allocate(Nspecies(Ngrids))
      IF (.not. allocated(Nyearclass))                                  &
     &                   allocate(Nyearclass(Ngrids))
!   Nfish=Nfishperyear*Nspecies*Nyearclass
      IF (.not. allocated(Nfish))                                       &
     &                   allocate(Nfish(Ngrids))

# ifdef PREDATOR
      IF (.not. allocated(Npredperspecies))                             &
     &                   allocate(Npredperspecies(Ngrids))
      IF (.not. allocated(Npredspecies))                                &
     &                   allocate(Npredspecies(Ngrids))
      IF (.not. allocated(Npred))                                       &
     &                   allocate(Npred(Ngrids))
# endif
# ifdef FISHING_FLEET
      IF (.not. allocated(Nboats))                                      &
     &                   allocate(Nboats(Ngrids))
      IF (.not. allocated(Nports))                                      &
     &                   allocate(Nports(Ngrids))
# endif
#endif
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
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
