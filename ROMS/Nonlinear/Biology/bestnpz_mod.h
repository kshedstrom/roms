!=======================================================================
!   Parameters for BEST NPZ Biological Model
!      Georgina Gibson 
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
#ifdef JELLY
!  Jellyfish Parameters
        real(r8) :: eJel, gammaJel
        real(r8) :: respJel,mpredJel
        real(r8) :: fpCopJel, fpNCaJel, fpEupJel
        real(r8) :: Q10Jelr, Q10JelTr,Q10Jele, Q10JelTe
        real(r8) :: bmJ,ktbmJ,TrefJ
#endif

!  Phytoplankton senescence
        real(r8) :: minmPhS, maxmPhS, NcritPhS
        real(r8) :: minmPhL, maxmPhL, NcritPhL
!  Zoopkankton mortality
        real(r8) :: mMZS, mMZL, mCop, mNCa, mEup
!  predation closure
        real(r8) :: mpredCop, mpredNCa, mpredEup
        real(r8) :: mpredMZS, mpredMZL
!  sinking
        real(r8) :: wPhS, wPhL, wDet, wDetF
!  Terms to define the Iron climatology field
        real(r8) :: Feinlo, Feinhi, Feinh, Feofflo, Feoffhi, Feoffh

!  Terms to define respiration
        real(r8) :: respPhS, respPhL, respMZS, respMZL
        real(r8) :: respCop, respNCa, respEup
        real(r8) :: TmaxPhS,TminPhS, Topt_PhS, KtBm_PhS
        real(r8) :: TmaxPhL, TminPhL, Topt_PhL, KtBm_PhL
        real(r8) :: TmaxMZS, KtBm_MZS, TmaxMZL, KtBm_MZL
        real(r8) :: ktbmC,TrefC
        real(r8) :: ktbmN,TrefN
        real(r8) :: ktbmE,TrefE
!  Detrital Remineralization and Nitrification
        real(r8) :: regen, dgrad
        real(r8) :: Pv0, PvT
        real(r8) :: KnT, Nitr0,ToptNtr,ktntr,KNH4Nit
        real(r8) :: tI0,KI

!  Iron limitation terms
        real(r8) :: kfePhS, kfePhL, FeC
!  Diapause
        real(r8) :: NCmaxz
        real(r8) :: wNCrise,wNCsink
        real(r8) :: RiseStart, RiseEnd, SinkStart, SinkEnd
#ifdef BENTHIC
        real(r8) :: bmB,ktbmB,TrefB
        real(r8) :: iremin
        real(r8) :: q10,q10r
        real(r8) :: Rup,KupP,LupD, LupP,KupD
        real(r8) :: Qres,Rres,rmort,eex,eexD,BenPred
        real(r8) :: prefD,prefPS,prefPL,T0ben,T0benr
#endif
#ifdef ICE_BIO
        real(r8) :: alphaIb, betaI,  inhib
        real(r8) :: ksnut1,ksnut2,mu0, R0i
        real(r8) :: rg0,rg,annit
        real(r8) :: aidz
#endif
