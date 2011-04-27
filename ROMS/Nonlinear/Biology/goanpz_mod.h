!=======================================================================
!   Parameters for Gulf of Alaska Biological Model
!      Sarah Hinckley and Elizabeth Dobbins
!=======================================================================
!   PARfrac  Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!=======================================================================
      USE mod_param

      implicit none

!
!  Set biological tracer identification indices.
!
      integer, pointer :: idbio(:)    ! Biological tracers
      integer :: iNO3                 ! Nitrate
      integer :: iNH4                 ! Ammonium
      integer :: iPhS                 ! Small Phytoplankton
      integer :: iPhL                 ! Large Phytoplankton
      integer :: iMZS                 ! Small Microzooplankton
      integer :: iMZL                 ! Large Microzooplankton
      integer :: iCop                 ! Small Coastal Copepods
      integer :: iNCa                 ! Neocalanus
      integer :: iEup                 ! Euphausiids
      integer :: iDet                 ! Detritus
      integer :: iPhSprd              ! Small Phytoplankton Production
      integer :: iPhLprd              ! Large Phytoplankton Production
      integer :: iMZSprd              ! Small Microzooplankton Production
      integer :: iMZLprd              ! Large Microzooplankton Production
      integer :: iCopPrd              ! Copepod production
      integer :: iNCaPrd              ! Neocalanus production
      integer :: iEupPrd              ! Euphausiid production
# ifdef IRON_LIMIT
      integer :: iFe                  ! Iron
# endif
      integer :: NBTS
!
!  Biological parameters.
!
      integer, allocatable :: BioIter(:)
      real(r8) :: VertMixIncr
      real(r8), allocatable :: PARfrac(:)       ! nondimensional
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
#if defined BIOFLUX
      real(r8), dimension(NAT+NBT,NAT+NBT) :: bflx = 0.0_r8
#endif

      integer, allocatable :: idTSvar(:)    ! Stationary production variables
      integer, allocatable :: hisTSid(:,:)  ! history St tracer IDs
      integer, allocatable :: avgTSid(:,:)  ! averages stationary tracer IDs
      integer, allocatable :: avg2TSid(:,:) ! averages stationary tracer IDs

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
!  Determine number of biological tracers.
!-----------------------------------------------------------------------
!
      NBTS = 7
#if defined IRON_LIMIT
      NBT = 11
#else
      NBT = 10
#endif

!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(Parfrac)) THEN
        allocate ( Parfrac(Ngrids) )
      END IF

!-----------------------------------------------------------------------

      allocate ( idTSvar(MST) )
      allocate ( hisTSid(MST,Ngrids) )
      allocate ( avgTSid(MST,Ngrids) )
      allocate ( avg2TSid(MST,Ngrids) )
!
!  Initialize biology diagnostic indices.
!
      ic=0
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
!
!  Set identification indices.
!
      ic=NAT+NPT+NCS+NNS
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      iNO3_=ic+1
      iNH4_=ic+2
      iSiOH=ic+3
      iSphy=ic+4
      iLphy=ic+5
      iSzoo=ic+6
      iLzoo=ic+7
      iSDet=ic+8
      iopal=ic+9
      iPO4_=ic+10
      ic=ic+10
#ifdef OXYGEN
      iOxyg=ic+1
      ic=ic+1
#endif 
#ifdef CARBON
      iTIC_=ic+1
      iTAlk=ic+2
      ic=ic+2
#endif

      !---------------------------------------------
      !Adding stationary production tracers to model
      !---------------------------------------------
      MST = 7
      DO ng=1,Ngrids
        NTS(ng) = 7
      END DO

      RETURN
      END SUBROUTINE initialize_biology
