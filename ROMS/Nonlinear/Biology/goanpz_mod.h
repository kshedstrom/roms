!=======================================================================
!   Parameters for Gulf of Alaska Biological Model
!      Sarah Hinckley and Elizabeth Dobbins
!=======================================================================
!   PARfrac  Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!=======================================================================
        USE mod_param

        implicit none

        integer, dimension(Ngrids) :: BioIter
        real(r8) :: VertMixIncr
        real(r8), dimension(Ngrids) :: PARfrac       ! nondimensional
!  Bio- conversions
        real(r8) :: xi, ccr, ccrPhL
!  extinction coefficients
        real(r8) :: k_ext, k_chl
!  PhS growth params
        real(r8) :: DiS, DpS
        real(r8) :: alphaPhS
        real(r8) :: psiPhS
        real(r8) :: k1PhS, k2PhS
!  PhL growth params
        real(r8) :: DiL
        real(r8) :: DpL
        real(r8) :: alphaPhL
        real(r8) :: psiPhL
        real(r8) :: k1PhL
        real(r8) :: k2PhL
!  MZS preference params
        real(r8) :: fpPhSMZS, fpPhLMZS
!  MZS growth and feeding params
        real(r8) :: eMZS
        real(r8) :: Q10MZS
        real(r8) :: Q10MZST
        real(r8) :: fMZS
        real(r8) :: kMZS
        real(r8) :: gammaMZS
!  MZL preferences params
        real(r8) :: fpPhSMZL, fpPhLMZL, fpMZSMZL
!  MZL growth and feeding params
        real(r8) :: eMZL
        real(r8) :: Q10MZL
        real(r8) :: Q10MZLT
        real(r8) :: fMZL
        real(r8) :: kMZL
        real(r8) :: gammaMZL
!  Cop preference params
        real(r8) :: fpPhSCop, fpPhLCop, fpMZSCop, fpMZLCop
!  Cop growth and feeding params
        real(r8) :: eCop
        real(r8) :: Q10Cop
        real(r8) :: Q10CopT
        real(r8) :: fCop
        real(r8) :: gammaCop
        real(r8) :: kCop
!  NCa preference params
        real(r8) :: fpPhSNCa, fpPhLNCa, fpMZSNCa, fpMZLNCa
!  NCa growth and feeding params
        real(r8) :: eNCa
        real(r8) :: Q10NCa
        real(r8) :: Q10NCaT
        real(r8) :: fNCa
        real(r8) :: gammaNCa
        real(r8) :: kNCa
!  Eup preference params
        real(r8) :: fpPhSEup, fpPhLEup, fpMZSEup, fpMZLEup, fpCopEup
!  Eup growth and feeding params
        real(r8) :: eEup
        real(r8) :: Q10Eup
        real(r8) :: Q10EupT
        real(r8) :: fEup
        real(r8) :: gammaEup
        real(r8) :: kEup
!  Phytoplankton senescence
        real(r8) :: minmPhS, maxmPhS, NcritPhS
        real(r8) :: minmPhL, maxmPhL, NcritPhL
!  Zoopkankton mortality
        real(r8) :: mMZS, mMZL, mCop, mNCa, mEup
!  predation closure
        real(r8) :: mpredCop, mpredNCa, mpredEup
        real(r8) :: mpredMZS, mpredMZL
!  sinking and regeneration terms
        real(r8) :: regen, dgrad
        real(r8) :: wPhS, wPhL, wDet, terms
!  Terms to define the Iron climatology field
        real(r8) :: Feinlo, Feinhi, Feinh, Feofflo, Feoffhi, Feoffh
!  Terms to define respiration
        real(r8) :: respPhS, respPhL, respMZS, respMZL
        real(r8) :: respCop, respNCa, respEup
!  Iron limitation terms
        real(r8) :: kfePhS, kfePhL, FeC
!  Diapause
        real(r8) :: NCmaxz
        real(r8) :: wNCrise,wNCsink
        real(r8) :: RiseStart, RiseEnd, SinkStart, SinkEnd
#  if defined BIOFLUX
        real(r8), dimension(NAT+NBT,NAT+NBT) :: bflx = 0.0_r8
#  endif
