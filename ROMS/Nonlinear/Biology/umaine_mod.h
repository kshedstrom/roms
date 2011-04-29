!                                                                      !
!  Parameters for UMaine CoSiNE model:                                 !
!                                                                      !
!   reg1     Microzooplankton excretion rate to ammonium [1/day].      !
!   reg2     Mesozooplankton excretion rate to ammonium [1/day].       !
!   gmaxs1   Maximum specific growth rate of small phytoplankton       !
!              [1/day]                                                 !
!   gmaxs2   Maximum specific growth rate of diatom [1/day]            !
!   beta1    Microzooplankton maximum grazing rate [1/day]             !
!   beta2    Mesozooplankton maximum grazing rate [1/day]              !
!   akz1     Half saturation constant for microzooplankton grazing     !
!              [mmol_N/m3]                                             !
!   akz2     Half saturation constant for mesozooplankton grazing      !
!              [mmol_N/m3]                                             !
!   PARfrac  Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!   amaxs2   Initial slope of P-I curve of small phytoplankton         !
!              [1/(Watts/m2)/day]                                      !
!   parsats1 PAR saturation onset parameter of small phytoplankton     !
!              [Watts/m2]                                              !
!   parsats2 PAR saturation onset parameter of diatom [Watts/m2]       !
!              [Watts/m2]                                              !
!   pis1     Ammonium inhibition parameter for small phytoplankton     !
!              [mmol_N/m3]                                             !
!   pis2     Ammonium inhibition parameter for diatom [mmol_N/m3]      !
!   akno3s1  Half saturation concentration for nitrate uptake by       !
!              small phytoplankton [mmol_N/m3].                        !
!   akno3s2  Half saturation concentration for nitrate uptake by       !
!              diatom [mmol_N/m3].                                     !
!   aknh4s1  Half saturation concentration for ammonium uptake by      !
!              small phytoplankton [mmol_N/m3].                        !
!   aknh4s2  Half saturation concentration for ammonium uptake by      !
!              diatom [mmol_N/m3].                                     !
!   akpo4s1  Half saturation concentration for phosphate uptake by     !
!              small phytoplankton [mmol_P/m3].                        !
!   akpo4s2  Half saturation concentration for phosphate uptake by     !
!              diatom [mmol_P/m3].                                     !
!   akco2s1  Half saturation concentration for co2 uptake by           !
!              small phytoplankton [mmol_C/m3].                        !
!   akco2s2  Half saturation concentration for co2 uptake by           !
!              diatom [mmol_C/m3].                                     !
!   aksio4s2 Half saturation constant for silicate uptake by           !
!              diatom [mmol_Si/m3].                                    !
!   ak1      Light attenuation coefficient of water [1/m]              !
!   ak2      Specific light attenuation coefficient for                !
!              phytoplankton [1/m/(mmol_N/m3)].                        !
!   bgamma0   Mesozooplankton specific mortality rate [1/day].         !
!   bgamma1   Grazing efficiency of microzooplankton [nondimensional]. !
!   bgamma2   Grazing efficiency of mesozooplankton [nondimensional].  !
!   bgamma3   Death rate of small phytoplankton [1/day].               !
!   bgamma4   Death rate of large phytoplankton [1/day].               !
!   bgamma5   Decay rate of detritus [1/day].                          !
!   bgamma6                                                            !
!   bgamma7   Nitrafication rate [1/day].                              !
!   wsd      Sinking velocity of detritus [m/day].                     !
!   wsdsi    Sinking velocity of detritus silicate [m/day].            !
!   wsp      Sinking velocity of large phytoplankton [m/day].          !
!   pco2a    Air pCO2 [ppmv].                                          !
!   si2n     Silicate to nitrogen ratio [mol_Si/mol_N].                !
!   p2n      Phosphorus to nitrogen ratio [mol_P/mol_N].               !
!   o2no     Oxygen to nitrate ratio [mol_O2/mol_NO3].                 !
!   o2nh     Oxygen to ammonium ratio [mol_O2/mol_NH4].                !
!   c2n      Carbon to nitrogen ratio [mol_C/mol_N].                   !
!   ro5      Grazing preference for diatom [nondimensional].           !
!   ro6      Grazing preference for mesozooplankton [nondimensional]   !
!   ro7      Grazing preference for detritus [nondimensional].         !
!=======================================================================
!
      USE mod_param

      implicit none
!
!  Set biological tracer identification indices.
!
      integer, allocatable :: idbio(:) ! Biological tracers
      integer :: iNO3_                 ! Nitrate concentration
      integer :: iNH4_                 ! Ammonium concentration
      integer :: iSiOH                 ! Silicate concentration
      integer :: iSphy                 ! Samll phytoplankton
      integer :: iLphy                 ! Diatom concentration
      integer :: iSzoo                 ! Small zooplankotn concentration
      integer :: iLzoo                 ! Mesozooplankotn concentration
      integer :: iSDet                 ! Detritus notrogen concentration
      integer :: iopal                 ! Biogenic silicate concentration
      integer :: iPO4_                 ! Phosphate concentration
#ifdef OXYGEN
      integer :: iOxyg                 ! Dissolved oxygen concentration
#endif
#ifdef CARBON
      integer :: iTIC_                 ! Total inorganic carbon
      integer :: iTAlk                 ! Total alkalinity
#endif

      integer, allocatable :: BioIter(:)

      real(r8), allocatable :: reg1(:)            ! 1/day
      real(r8), allocatable :: reg2(:)            ! 1/day
      real(r8), allocatable :: gmaxs1(:)          ! 1/day
      real(r8), allocatable :: gmaxs2(:)          ! 1/day
      real(r8), allocatable :: beta1(:)           ! 1/day
      real(r8), allocatable :: beta2(:)           ! 1/day
      real(r8), allocatable :: akz1(:)            ! mmol_N/m3
      real(r8), allocatable :: akz2(:)            ! mmol_N/m3
      real(r8), allocatable :: PARfrac(:)         ! nondimensional
      real(r8), allocatable :: amaxs2(:)          ! 1/(Watts/m2)/day
      real(r8), allocatable :: pis1(:)            ! m3/mmol_N
      real(r8), allocatable :: pis2(:)            ! m3/mmol_N
      real(r8), allocatable :: akno3s1(:)         ! mmol_N/m3
      real(r8), allocatable :: aknh4s1(:)         ! mmol_N/m3
      real(r8), allocatable :: akpo4s1(:)         ! mmol_P/m3
      real(r8), allocatable :: akco2s1(:)         ! mmol_C/m3
      real(r8), allocatable :: akno3s2(:)         ! mmol_N/m3
      real(r8), allocatable :: aknh4s2(:)         ! mmol_N/m3
      real(r8), allocatable :: aksio4s2(:)        ! mmol_Si/m3
      real(r8), allocatable :: akpo4s2(:)         ! mmol_P/m3
      real(r8), allocatable :: akco2s2(:)         ! mmol_C/m3
      real(r8), allocatable :: ak1(:)             ! 1/m
      real(r8), allocatable :: ak2(:)             ! 1/m/(mmol_N/m3)
      real(r8), allocatable :: parsats1(:)        ! Watts/m2
      real(r8), allocatable :: parsats2(:)        ! Watts/m2
      real(r8), allocatable :: bgamma0(:)         ! 1/day
      real(r8), allocatable :: bgamma1(:)         ! [nondimensional]
      real(r8), allocatable :: bgamma2(:)         ! [nondimensional]
      real(r8), allocatable :: bgamma3(:)         ! 1/day
      real(r8), allocatable :: bgamma4(:)         ! 1/day
      real(r8), allocatable :: bgamma5(:)         ! 1/day
      real(r8), allocatable :: bgamma6(:)         ! 
      real(r8), allocatable :: bgamma7(:)         ! 1/day
      real(r8), allocatable :: wsd(:)             ! m/day
      real(r8), allocatable :: wsdsi(:)           ! m/day
      real(r8), allocatable :: wsp(:)             ! m/day
      real(r8), allocatable :: si2n(:)            ! mol_Si/mol_N
      real(r8), allocatable :: pco2a(:)           ! ppmv
      real(r8), allocatable :: p2n(:)             ! mol_P/mol_N
      real(r8), allocatable :: o2no(:)            ! mol_O2/mol_NO3
      real(r8), allocatable :: o2nh(:)            ! mol_O2/mol_NH4
      real(r8), allocatable :: c2n(:)             ! mol_C/mol_N
      real(r8), allocatable :: ro5(:)             ! nondimensional
      real(r8), allocatable :: ro6(:)             ! nondimensional
      real(r8), allocatable :: ro7(:)             ! nondimensional
      real(r8), allocatable :: pCO2air(:)         ! ppmv

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
#ifdef CARBON
# ifdef OXYGEN
      NBT=13
# else
      NBT=12
# endif
#else
# ifdef OXYGEN
      NBT=11
# else
      NBT=10
# endif
#endif
!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(reg1)) THEN
        allocate ( reg1(Ngrids) )
      END IF
      IF (.not.allocated(reg2)) THEN
        allocate ( reg2(Ngrids) )
      END IF
      IF (.not.allocated(gmaxs1)) THEN
        allocate ( gmaxs1(Ngrids) )
      END IF
      IF (.not.allocated(gmaxs2)) THEN
        allocate ( gmaxs2(Ngrids) )
      END IF
      IF (.not.allocated(beta1)) THEN
        allocate ( beta1(Ngrids) )
      END IF
      IF (.not.allocated(beta2)) THEN
        allocate ( beta2(Ngrids) )
      END IF
      IF (.not.allocated(akz1)) THEN
        allocate ( akz1(Ngrids) )
      END IF
      IF (.not.allocated(akz2)) THEN
        allocate ( akz2(Ngrids) )
      END IF
      IF (.not.allocated(PARfrac)) THEN
        allocate ( PARfrac(Ngrids) )
      END IF
      IF (.not.allocated(amaxs2)) THEN
        allocate ( amaxs2(Ngrids) )
      END IF
      IF (.not.allocated(pis1)) THEN
        allocate ( pis1(Ngrids) )
      END IF
      IF (.not.allocated(pis2)) THEN
        allocate ( pis2(Ngrids) )
      END IF
      IF (.not.allocated(akno3s1)) THEN
        allocate ( akno3s1(Ngrids) )
      END IF
      IF (.not.allocated(aknh4s1)) THEN
        allocate ( aknh4s1(Ngrids) )
      END IF
      IF (.not.allocated(akpo4s1)) THEN
        allocate ( akpo4s1(Ngrids) )
      END IF
      IF (.not.allocated(akco2s1)) THEN
        allocate ( akco2s1(Ngrids) )
      END IF
      IF (.not.allocated(akno3d2)) THEN
        allocate ( akno3d2(Ngrids) )
      END IF
      IF (.not.allocated(aknh4s2)) THEN
        allocate ( aknh4s2(Ngrids) )
      END IF
      IF (.not.allocated(aksio4s2)) THEN
        allocate ( aksio4s2(Ngrids) )
      END IF
      IF (.not.allocated(akpo4s2)) THEN
        allocate ( akpo4s2(Ngrids) )
      END IF
      IF (.not.allocated(akco2s2)) THEN
        allocate ( akco2s2(Ngrids) )
      END IF
      IF (.not.allocated(ak1)) THEN
        allocate ( ak1(Ngrids) )
      END IF
      IF (.not.allocated(ak2)) THEN
        allocate ( ak2(Ngrids) )
      END IF
      IF (.not.allocated(parsats1)) THEN
        allocate ( parsats1(Ngrids) )
      END IF
      IF (.not.allocated(parsats2)) THEN
        allocate ( parsats2(Ngrids) )
      END IF
      IF (.not.allocated(bgamma0)) THEN
        allocate ( bgamma0(Ngrids) )
      END IF
      IF (.not.allocated(bgamma1)) THEN
        allocate ( bgamma1(Ngrids) )
      END IF
      IF (.not.allocated(bgamma2)) THEN
        allocate ( bgamma2(Ngrids) )
      END IF
      IF (.not.allocated(bgamma3)) THEN
        allocate ( bgamma3(Ngrids) )
      END IF
      IF (.not.allocated(bgamma4)) THEN
        allocate ( bgamma4(Ngrids) )
      END IF
      IF (.not.allocated(bgamma5)) THEN
        allocate ( bgamma5(Ngrids) )
      END IF
      IF (.not.allocated(bgamma6)) THEN
        allocate ( bgamma6(Ngrids) )
      END IF
      IF (.not.allocated(bgamma7)) THEN
        allocate ( bgamma7(Ngrids) )
      END IF
      IF (.not.allocated(wsd)) THEN
        allocate ( wsd(Ngrids) )
      END IF
      IF (.not.allocated(wsdsi)) THEN
        allocate ( wsdsi(Ngrids) )
      END IF
      IF (.not.allocated(wsp)) THEN
        allocate ( wsp(Ngrids) )
      END IF
      IF (.not.allocated(si2n)) THEN
        allocate ( si2n(Ngrids) )
      END IF
      IF (.not.allocated(pco2a)) THEN
        allocate ( pco2a(Ngrids) )
      END IF
      IF (.not.allocated(p2n)) THEN
        allocate ( p2n(Ngrids) )
      END IF
      IF (.not.allocated(o2no)) THEN
        allocate ( o2no(Ngrids) )
      END IF
      IF (.not.allocated(o2nh)) THEN
        allocate ( o2nh(Ngrids) )
      END IF
      IF (.not.allocated(c2n)) THEN
        allocate ( c2n(Ngrids) )
      END IF
      IF (.not.allocated(ro5)) THEN
        allocate ( ro5(Ngrids) )
      END IF
      IF (.not.allocated(ro6)) THEN
        allocate ( ro6(Ngrids) )
      END IF
      IF (.not.allocated(ro7)) THEN
        allocate ( ro7(Ngrids) )
      END IF
      IF (.not.allocated(pCO2air)) THEN
        allocate ( pCO2air(Ngrids) )
      END IF
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
# ifdef OXYGEN
      iOxyg=ic+1
      ic=ic+1
# endif
# ifdef CARBON
      iTIC_=ic+1
      iTAlk=ic+2
      ic=ic+2
# endif

      RETURN
      END SUBROUTINE initialize_biology
