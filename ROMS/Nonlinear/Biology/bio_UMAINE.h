      SUBROUTINE biology (ng,tile)
!
!svn $Id: bio_UMAINE.h 702 2008-08-12 16:44:47Z kate $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  UMaine CoSiNE ecosystem model version 1.0 (March 25, 2009)          !
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  then the global biological fields. The model is based on the        !
!  Carbon, Silicon, Nitrogen Ecosystem (CoSiNE) model (Chai et al.,    !
!  2002). The model state variables are:                               !
!                                                                      !
!    iNO3_      Nitrate                                                !  
!    iSiOH      Silicate                                               !
!    iNH4_      Ammonium                                               !
!    iSphy      Small Phytoplankton                                    !
!    iLphy      Diatoms                                                !
!    iSzoo      Micro Zooplankton                                      !
!    iLphy      Meso Zooplankton                                       !
!    iSDet      Detritus-nitrogen                                      !
!    iopal      Detritus-silicate                                      !
!    iPH4_      Phosphate                                              !
!    iOxyg      Dissolved Oxygen                                       !
!    ITIC_      Total CO2                                              !
!    iTAlk      Total Alkalinity                                       !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Chai, F., R.C. Dugdale, T-H Peng, F.P.Wilkerson, and R.T. Barber  !
!      (2002): One dimensional Ecosystem Model of the Equatorial       !
!      Pacific Upwelling System, Part I: Model Development and Silicon !
!      and Nitorgen Cycle. Deep-Sea Res. II, Vol. 49, No. 13-14,       !
!      2713-2745.                                                      !
!    Dugdale, R.C., R. T. Barber, F. Chai, T.H. Peng, and              !
!      F.P. Wilkerson (2002): One Dimensional Ecosystem Model of the   !
!      Equatorial Pacific Upwelling System, Part II: Sensitivity       !
!      Analysis and Comparison with JGOFS EqPac Data. Deep-Sea Res.    !
!      II, Vol. 49, No. 13-14, 2746-2762.                              !
!                                                                      !
!  Adapted from 1D code develped by Fei Chai of UMaine, and further    !
!  modified and developed by Lei Shi and adding 4 more biological      !
!  statevariables (phospate, dissolved oxygen, total co2 and           !
!  total akalinity). The 3D ROMS-COSINE model code implemented         !
!  by Lei Shi of UMaine.                                               !
!                                                                      !
!  Release Note: current release have no CO2 sea surafce flux.         !
!    Should be there soon. Meanwhile it is recommended                 !
!    to undef CARBON.                                                  !
!                                                                      !
!***********************************************************************
!
!  Additional reference:                                              
!                                                                     
!    Chai, F., M. Jiang, R.T. Barber, R.C. Dugdale, and               
!      Y. Chao (2003): Interdecadal Variation of the Transition Zone  
!      Chlorophyll Front, A Physical-Biological Model Simulation      
!      between 1960 and 1990. Journal of Oceanography, Vol. 59,       
!      461-475.                                                       
!    Chai, F., M-S Jiang, Y. Chao, R.C. Dugdale, F. Chavez, and       
!      R.T. Barber (2007): Modeling Responses of Diatom Productivity  
!      and Biogenic Silica Export to Iron Enrichment in the Equatorial
!      Pacific Ocean. Global Biogeochemical Cycle, Vol. 21, GB3S90,   
!      doi:10.1029/2006GB002804.                                      
!    Jiang, M-S, F. Chai, R.T. Barber, R.C. Dugdale, F. Wilkerson,    
!      and T-H Peng (2003). A nitrate and silicate budget in the      
!      Equatorial Pacific Ocean: A coupled biological-physical model  
!      study. Deep Sea Res. II, Vol. 50 (22-26), 2971-2996.           
!    Liu, G. and F. Chai (2009): Seasonal and interannual variability 
!      of primary and export production in the South China Sea: A     
!      three-dimensional physical-biogeochemical model study. ICES    
!      Journal of Marine Science, 66, 420-431.                        
!    Liu, G. and F. Chai (2009): Seasonal and interannual variation of
!      physical and biological processes during 1994-2001 in the Sea  
!      of Japan/East Sea: a three-dimensional physical-biogeochemical 
!      modeling study. Journal of Marine Systems, in press.           
!    Polovina, J. J., F. Chai, E. A. Howell, D. R. Kobayashi, L. Shi, 
!      and Y. Chao (2008): Ecosystem dynamics at a productivity       
!      gradient: a study of the lower trophic dynamics around the     
!      northern atolls in the Hawaiian Archipelago. Progress in       
!      Oceanography, 77, 217-224.                                     
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
!  Set header file name.
!
#ifdef DISTRIBUTE
      IF (Lbiofile(iNLM)) THEN
#else
      IF (Lbiofile(iNLM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iNLM)=.FALSE.
        BIONAME(iNLM)=__FILE__
      END IF
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
     &                   FORCES(ng) % srflx,                            &
#if defined OXYGEN || defined CARBON
# ifdef BULK_FLUXES
     &                   FORCES(ng) % Uwind,                            &
     &                   FORCES(ng) % Vwind,                            &
# else
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
# endif
#endif
#ifdef CARBON
     &                   OCEAN(ng) % pH,                                &
#endif
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif
      RETURN
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
#ifdef MASKING
     &                         rmask,                                   &
#endif
     &                         Hz, z_r, z_w, srflx,                     &
#if defined OXYGEN || defined CARBON
# ifdef BULK_FLUXES
     &                         Uwind, Vwind,                            &
# else
     &                         sustr, svstr,                            &
# endif
#endif
#ifdef CARBON
     &                         pH,                                      &
#endif
     &                         t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
# if defined OXYGEN || defined CARBON
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
#  else
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
#  endif
# endif
# ifdef CARBON
      real(r8), intent(inout) :: pH(LBi:,LBj:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
# if defined OXYGEN || defined CARBON
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
#  endif
# endif
# ifdef CARBON
      real(r8), intent(inout) :: pH(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
      integer, parameter :: Nsink = 3
#if defined OXYGEN || defined CARBON
      real(r8) :: u10squ, u10spd
#endif

      integer :: Iter, i, indx, isink, ibio, ivar, j, k, ks

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: Minval = 1.0e-6_r8
      real(r8), parameter :: zeptic = 100.0_r8

      real(r8) :: dtdays

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5

      real(r8), dimension(Nsink) :: Wbio

      real(r8), dimension(IminS:ImaxS) :: PARsur

#if defined OXYGEN || defined CARBON
      real(r8), dimension(IminS:ImaxS) :: kw660
#endif

#ifdef OXYGEN
      real(r8), dimension(IminS:ImaxS) :: o2sat
      real(r8), dimension(IminS:ImaxS) :: o2flx
#endif

#ifdef CARBON
      real(r8), dimension(IminS:ImaxS) :: pCO2
      real(r8), dimension(IminS:ImaxS) :: co2flx
#endif

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_bak

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3

      real(r8), dimension(IminS:ImaxS,N(ng)+1) :: PIO
      real(r8), dimension(IminS:ImaxS,N(ng)) :: PAR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ADPT

      real(r8), dimension(N(ng)) :: sinkindx

      real(r8) :: thick
      real(r8) :: alts1, alts2, grows1, grows2, OXR, Q10
      real(r8) :: uno3s2, unh4s2, uno3s1, unh4s1
      real(r8) :: cff0, cff9, cff7, cff8, cff6
      real(r8) :: xco2_in
      real(r8), parameter :: AKOX = 30.0_r8
#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Set time-stepping according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)
!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iLphy
      idsink(2)=iSDet
      idsink(3)=iopal
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wsp(ng)                ! iLphy
      Wbio(2)=wsd(ng)                ! iSDet
      Wbio(3)=wsdsi(ng)              ! iopal
!
!  Compute inverse thickness to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
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
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
        DO ibio=1,NBT
          indx=idbio(ibio)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_bak(i,k,indx)=MAX(t(i,j,k,nstp,indx),0.0_r8)
              Bio(i,k,indx)=Bio_bak(i,k,indx)
            END DO
          END DO
        END DO
#ifdef CARBON
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,iTIC_)=MIN(Bio(i,k,iTIC_),3000.0_r8)
            Bio(i,k,iTIC_)=MAX(Bio(i,k,iTIC_),400.0_r8)
          END DO
        END DO
#endif
!
!  Extract potential temperature and salinity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),35.0_r8)
            Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
!  During the iterative procedure a series of fractional time steps are
!  performed in a chained mode (splitting by different biological
!  conversion processes) in sequence of the main food chain.  In all
!  stages the concentration of the component being consumed is treated
!  in fully implicit manner, so the algorithm guarantees non-negative
!  values, no matter how strong s the concentration of active consuming
!  component (Phytoplankton or Zooplankton).  The overall algorithm,
!  as well as any stage of it, is formulated in conservative form
!  (except explicit sinking) in sense that the sum of concentration of
!  all components is conserved.
!
!
!  In the implicit algorithm, we have for example (N: nitrate,
!                                                  P: phytoplankton),
!
!     N(new) = N(old) - uptake * P(old)     uptake = mu * N / (Kn + N)
!                                                    {Michaelis-Menten}
!  below, we set
!                                           The N in the numerator of
!     cff = mu * P(old) / (Kn + N(old))     uptake is treated implicitly
!                                           as N(new)
!
!  so the time-stepping of the equations becomes:
!
!     N(new) = N(old) / (1 + cff)     (1) when substracting a sink term,
!                                         consuming, divide by (1 + cff)
!  and
!
!     P(new) = P(old) + cff * N(new)  (2) when adding a source term,
!                                         growing, add (cff * source)
!
!  Notice that if you substitute (1) in (2), you will get:
!
!     P(new) = P(old) + cff * N(old) / (1 + cff)    (3)
!
!  If you add (1) and (3), you get
!
!     N(new) + P(new) = N(old) + P(old)
!
!  implying conservation regardless how "cff" is computed. Therefore,
!  this scheme is unconditionally stable regardless of the conversion
!  rate. It does not generate negative values since the constituent
!  to be consumed is always treated implicitly. It is also biased
!  toward damping oscillations.
!
!  The iterative loop below is to iterate toward an universal Backward-
!  Euler treatment of all terms. So if there are oscillations in the
!  system, they are only physical oscillations. These iterations,
!  however, do not improve the accuaracy of the solution.
!
        ITER_LOOP: DO Iter=1,BioIter(ng)
!
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
!
!  Compute attenuation coefficient based on the concentration of
!  microphytoplankotn and diatom concentration within each grid box. 
!  Then, photosynthetically available radiation (PAR) are calculated
!  as the depth-averaged PAR within the vertical grid. Since the light
!  penetreated from surface to bottom, so the PAR are calculated in the 
!  same order. PARsur is surface PAR value. PIO is the PAR value at 
!  surface or bottom of a vertical grid, or at w location vertically.
!
          DO i=Istr,Iend
            PIO(i,N(ng)+1)=PARsur(i)
      IF (PIO(i,N(ng)+1).lt.0) PIO(i,N(ng)+1)=0.0_r8
!            IF (PIO(i,N(ng)+1).gt.0.0_r8) THEN
              DO k=N(ng),1,-1
      cff1=(AK1(ng)+(Bio(i,k,iSphy)+Bio(i,k,iLphy))*AK2(ng))*HZ(i,j,k)
      PIO(i,K)=PIO(i,K+1)*EXP(-cff1)
      PAR(i,K)=(PIO(i,K+1)-PIO(i,K))/cff1
        ADPT(i,K) = 1.0_r8-4.0_r8*z_r(i,j,k)/zeptic
              END DO
!            END IF
          END DO

          DO k=1,N(ng)
            DO i=Istr,Iend
!
!-----------------------------------------------------------------------
!     CALCULATING the temperature dependence of biology processes
!-----------------------------------------------------------------------
!
      Q10=exp(0.069_r8*(Bio(i,k,itemp)-25.0_r8))
!
!-----------------------------------------------------------------------
!     CALCULATING THE OXIDATION RATE OF ORGANIC MATTER
!-----------------------------------------------------------------------
!
!  Any biology processes that consume oxygen will be limited by the
!  availability of dissolved oxygen except the bottom layer.
!
#ifdef OXYGEN
      if(k.lt.N(ng))then
      OXR = Bio(i,k,iOxyg)/(Bio(i,k,iOxyg)+AKOX)
      else
      OXR = 1.0_r8
      endif
#else
      OXR = 1.0_r8
#endif
!
!-----------------------------------------------------------------------
!     CALCULATING THE GROWTH RATE AS NO3,NH4, AND LIGHT;
!     GRAZING, PARTICLE SINKING AND REGENERATION
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!  small phytoplankton nutrient uptake, and growth
!-----------------------------------------------------------------------
!
!  nitrogen limited growth
	cff0=exp(-pis1(ng)*Bio(i,k,iNH4_))
	cff1=Bio(i,k,iNH4_)/aknh4s1(ng)+                                &
     &       cff0*Bio(i,k,iNO3_)/akno3s1(ng)
	cff3=1.0_r8+Bio(i,k,iNH4_)/aknh4s1(ng)+                         &
     &       cff0*Bio(i,k,iNO3_)/akno3s1(ng)
	cff4=cff1/cff3
!  phosphate limited growth
	cff5=Bio(i,k,iPO4_)/(akpo4s1(ng)+Bio(i,k,iPO4_))
!  tco2 limited growth
#ifdef CARBON
	cff6=Bio(i,k,iTIC_)/(akco2s1(ng)+Bio(i,k,iTIC_))
#else
	cff6=1.0_r8
#endif
!  total limitation,
	cff7=min(cff4,cff5,cff6)
	cff8=cff7/(Minval+cff4)
!  only cff0, cff3, and cff8 is used next
!        cff0=exp(-pis1(ng)*Bio(i,k,iNH4_))
	cff1=Bio(i,k,iSphy)/aknh4s1(ng)
	cff2=cff0*Bio(i,k,iSphy)/akno3s1(ng)
!	cff3=1.0_r8+Bio(i,k,iNH4_)/aknh4s1(ng)+                         &
!     &       cff0*Bio(i,k,iNO3_)/akno3s1(ng)
        alts1= 1.0_r8 - exp(-PAR(i,k)*ADPT(i,K)/parsats1(ng))
	cff4=dtdays*Q10*gmaxs1(ng)*alts1*cff8*cff1/cff3
	cff5=dtdays*Q10*gmaxs1(ng)*alts1*cff8*cff2/cff3

	Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff4)
	Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff5)
	cff1=cff4*Bio(i,k,iNH4_)+cff5*Bio(i,k,iNO3_)
	Bio(i,k,iSphy)=Bio(i,k,iSphy)+cff1
	Bio(i,k,iPO4_)=Bio(i,k,iPO4_)-                                  &
     &                p2n(ng)*cff1
#ifdef CARBON
	Bio(i,k,iTIC_)=Bio(i,k,iTIC_)-                                  &
     &                c2n(ng)*cff1
#endif
#ifdef OXYGEN
	Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                                  &
     &                o2nh(ng)*cff4*Bio(i,k,iNH4_)+                     &
     &                o2no(ng)*cff5*Bio(i,k,iNO3_)
#endif
!
!-----------------------------------------------------------------------
!  diatom nutrient uptake, and growth
!-----------------------------------------------------------------------
!
!  nitrogen limited growth
	cff0=exp(-pis2(ng)*Bio(i,k,iNH4_))
	cff1=Bio(i,k,iNH4_)/aknh4s2(ng)+                                &
     &       cff0*Bio(i,k,iNO3_)/akno3s2(ng)
	cff3=1.0_r8+Bio(i,k,iNH4_)/aknh4s2(ng)+                         &
     &       cff0*Bio(i,k,iNO3_)/akno3s2(ng)
	cff4=cff1/cff3
!  silicate limited growth
	cff2=Bio(i,k,iSiOH)/(aksio4s2(ng)+Bio(i,k,iSiOH))
!  phosphate limited growth
	cff5=Bio(i,k,iPO4_)/(akpo4s2(ng)+Bio(i,k,iPO4_))
!  tco2 limited growth
#ifdef CARBON
	cff6=Bio(i,k,iTIC_)/(akco2s2(ng)+Bio(i,k,iTIC_))
#else
	cff6=1.0_r8
#endif
!  total limitation,
	cff7=min(cff4,cff2,cff5,cff6)
	cff8=cff7/(Minval+cff4)
!  only cff0, cff3, and cff8 is used next
!        cff0=exp(-pis2(ng)*Bio(i,k,iNH4_))
	cff1=Bio(i,k,iLphy)/aknh4s2(ng)
	cff2=cff0*Bio(i,k,iLphy)/akno3s2(ng)
!	cff3=1.0_r8+Bio(i,k,iNH4_)/aknh4s2(ng)+                         &
!     &       cff0*Bio(i,k,iNO3_)/akno3s2(ng)
        alts2= 1.0_r8 - exp(-PAR(i,k)*ADPT(i,K)/parsats2(ng))
	cff4=dtdays*Q10*gmaxs2(ng)*alts2*cff8*cff1/cff3
	cff5=dtdays*Q10*gmaxs2(ng)*alts2*cff8*cff2/cff3

	Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff4)
	Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff5)
	cff1=cff4*Bio(i,k,iNH4_)+cff5*Bio(i,k,iNO3_)
	Bio(i,k,iLphy)=Bio(i,k,iLphy)+cff1
	Bio(i,k,iSiOH)=Bio(i,k,iSiOH)-                                  &
     &                 si2n(ng)*cff1
	Bio(i,k,iPO4_)=Bio(i,k,iPO4_)-                                  &
     &                 p2n(ng)*cff1
#ifdef CARBON
	Bio(i,k,iTIC_)=Bio(i,k,iTIC_)-                                  &
     &                 c2n(ng)*cff1
#endif
#ifdef OXYGEN
	Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                                  &
     &                 o2nh(ng)*cff4*Bio(i,k,iNH4_)+                    &
     &                 o2no(ng)*cff5*Bio(i,k,iNO3_)
#endif
!
!-----------------------------------------------------------------------
!  CALAULATING phytoplankotn mortality
!-----------------------------------------------------------------------
!
	cff1=dtdays*Q10*bgamma3(ng)
	Bio(i,k,iSphy)=Bio(i,k,iSphy)/(1.0_r8+cff1)
	cff2=dtdays*Q10*bgamma4(ng)
	Bio(i,k,iLphy)=Bio(i,k,iLphy)/(1.0_r8+cff2)
	Bio(i,k,iSDet)=Bio(i,k,iSDet)+                                  &
      &                cff1*Bio(i,k,iSphy)+cff2*Bio(i,k,iLphy)
	Bio(i,k,iopal)=Bio(i,k,iopal)+si2n(ng)*cff2*Bio(i,k,iLphy)
!
!-----------------------------------------------------------------------
!  CALAULATING zooplankotn excretion
!-----------------------------------------------------------------------
!
	cff1=dtdays*OXR*Q10*reg1(ng)
	Bio(i,k,iSzoo)=Bio(i,k,iSzoo)/(1.0_r8+cff1)
	cff2=dtdays*OXR*Q10*reg2(ng)
	Bio(i,k,iLzoo)=Bio(i,k,iLzoo)/(1.0_r8+cff2)
	cff3=cff1*Bio(i,k,iSzoo)+cff2*Bio(i,k,iLzoo)
	Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+cff3
	Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+p2n(ng)*cff3
#ifdef CARBON
	Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+c2n(ng)*cff3
#endif
#ifdef OXYGEN
# ifdef oxygen_bottom_layer_treatment
	if(k.ne.1)then
	Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-o2nh(ng)*cff3
	endif
# else
	Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-o2nh(ng)*cff3
# endif
#endif
!
!-----------------------------------------------------------------------
!     CALAULATING nitrification
!-----------------------------------------------------------------------
!
	cff1=dtdays*OXR*Q10*bgamma7(ng)
	Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff1)
	Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+cff1*Bio(i,k,iNH4_)
#ifdef OXYGEN
# ifdef oxygen_bottom_layer_treatment
	if(k.ne.1)then
	Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                                  &
     &                 2.0_r8*cff1*Bio(i,k,iNH4_)
	endif
# else
	Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                                  &
     &                 2.0_r8*cff1*Bio(i,k,iNH4_)
# endif
#endif
!
!-----------------------------------------------------------------------
!     CALAULATING reminalization
!-----------------------------------------------------------------------
!
	cff1=dtdays*OXR*Q10*bgamma5(ng)
	Bio(i,k,iSDet)=Bio(i,k,iSDet)/(1.0_r8+cff1)
	Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+cff1*Bio(i,k,iSDet)
	Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+                                  &
     &                 p2n(ng)*cff1*Bio(i,k,iSDet)
#ifdef CARBON
	Bio(i,k,iTIC_)=Bio(i,k,iTIC_)+cff1*Bio(i,k,iSDet)
#endif
#ifdef OXYGEN
# ifdef oxygen_bottom_layer_treatment
	if(k.ne.1)then
	Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                                  &
     &                 o2nh(ng)*cff1*Bio(i,k,iSDet)
	endif
# else
	Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                                  &
     &                 o2nh(ng)*cff1*Bio(i,k,iSDet)
# endif
#endif

	cff2=dtdays*Q10*bgamma5(ng)
	Bio(i,k,iopal)=Bio(i,k,iopal)/(1.0_r8+cff2)
	Bio(i,k,iSiOH)=Bio(i,k,iSiOH)+cff2*Bio(i,k,iopal)
!
!-----------------------------------------------------------------------
!     CALAULATING THE GRAZING RATE 
!-----------------------------------------------------------------------
!
	if(bio(i,k,iSphy).gt.0.005_r8)then
	cff1=dtdays*Q10*beta1(ng)*Bio(i,k,iSzoo)/                       &
     &         (akz1(ng)+Bio(i,k,iSphy))
	bio(i,k,iSphy)=bio(i,k,iSphy)/(1.0_r8+cff1)
	cff2=cff1*bio(i,k,iSphy)
	bio(i,k,iSzoo)=bio(i,k,iSzoo)+bgamma1(ng)*cff2
	bio(i,k,iSDet)=bio(i,k,iSDet)+(1.0-bgamma1(ng))*cff2
	endif

	if(bio(i,k,iLphy).gt.0.005_r8)then
	cff0=ro5(ng)*Bio(i,k,iLphy)+                                    &
     &       ro6(ng)*Bio(i,k,iSzoo)+ro7(ng)*Bio(i,k,iSDet)
	cff0=akz2(ng)*cff0+(ro5(ng)*Bio(i,k,iLphy)*Bio(i,k,iLphy)+      &
     &                  ro6(ng)*Bio(i,k,iSzoo)*Bio(i,k,iSzoo)+          &
     &                  ro7(ng)*Bio(i,k,iSDet)*Bio(i,k,iSDet))
	cff1=dtdays*Q10*beta2(ng)*ro5(ng)*                              &
     &       Bio(i,k,iLphy)*Bio(i,k,iLzoo)/cff0
	cff2=dtdays*Q10*beta2(ng)*ro6(ng)*                              &
     &       Bio(i,k,iSzoo)*Bio(i,k,iLzoo)/cff0
	cff3=dtdays*Q10*beta2(ng)*ro7(ng)*                              &
     &       Bio(i,k,iSDet)*Bio(i,k,iLzoo)/cff0
	bio(i,k,iLphy)=bio(i,k,iLphy)/(1.0_r8+cff1)
	bio(i,k,iSzoo)=bio(i,k,iSzoo)/(1.0_r8+cff2)
	bio(i,k,iSDet)=bio(i,k,iSDet)/(1.0_r8+cff3)
	cff4=cff1*Bio(i,k,iLphy)+                                       &
     &       cff2*Bio(i,k,iSzoo)+cff3*Bio(i,k,iSDet)
	bio(i,k,iLzoo)=bio(i,k,iLzoo)+bgamma2(ng)*cff4
	bio(i,k,iSDet)=bio(i,k,iSDet)+(1.0-bgamma2(ng))*cff4
	bio(i,k,iopal)=bio(i,k,iopal)+cff1*Bio(i,k,iLphy)
	endif
!
!-----------------------------------------------------------------------
!  CALAULATING mesozooplankton removal/mortality
!-----------------------------------------------------------------------
!
	cff1=dtdays*Q10*bgamma0(ng)*Bio(i,k,iLzoo)
	Bio(i,k,iLzoo)=Bio(i,k,iLzoo)/(1.0_r8+cff1)
	Bio(i,k,iSDet)=bio(i,k,iSDet)+                                  &
      &                cff1*bio(i,k,iLzoo)

              END DO
            END DO
	    
#if defined OXYGEN || defined CARBON
!
!-----------------------------------------------------------------------
!     CALCULATING gas transfer velocity at a Schmidt number of 660
!-----------------------------------------------------------------------
!
          k=N(ng)
              DO i=Istr,Iend
!
!  Compute wind speed.
!
# ifdef BULK_FLUXES
            u10squ=Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j)
# else
!
!  drag coefficient is 0.001, and air density is 1.2 kg per cube meter
!
      cff1=rho0/(0.001_r8*1.2_r8)
            u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+     &
     &                       (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
# endif
            u10spd=sqrt(u10squ)
!
!  Compute gas transfer velocity at a Schmidt number of 660.
!
! climatology wind speed (Wanninkhof & Mcgillis, 1999).
!      kw660=1.09*u10spd-0.333*u10squ+0.078*u10spd*u10squ      
!(in units of cm/hr), the one is too pronounced with large speed
! short-term (<1 day) winds (Wanninkhof & Mcgillis, 1999).
!      kw660=0.0283*u10spd*u10squ                     
!(in units of cm/hr)
!
      kw660=0.31*u10squ                     
!
!(in units of cm/hr)
!
!      kw660 (i)=0.39*u10squ
!
!(in units of cm/hr)
!
              END DO
#endif
#ifdef OXYGEN
!
!-----------------------------------------------------------------------
!     CALCULATING THE O2 SURFACE saturation concentration and FLUX  
!-----------------------------------------------------------------------
!
          k=N(ng)
          CALL O2_flux (Istr, Iend, LBi, UBi, LBj, UBj,                 &
     &                     IminS, ImaxS, j,                             &
# ifdef MASKING
     &                     rmask,                                       &
# endif
     &                     Bio(IminS:,k,itemp), Bio(IminS:,k,isalt),    &
     &                     Bio(IminS:,k,iOxyg), kw660,                  &
     &                     1.0_r8, o2sat, o2flx)     
              DO i=Istr,Iend
	bio(i,k,iOxyg)=bio(i,k,iOxyg)+dtdays*o2flx(i)*Hz_inv(i,k)
              END DO
#endif
#ifdef CARBON
!
!-----------------------------------------------------------------------
!  CALCULATING
!  Surface equilibrium partial pressure inorganic carbon (ppmv) at the
!  surface, and CO2 gas exchange.
!  Under current setting, no co2 gas exchange at sea surafce and 
!  subroutine CO2_flux is not included, the code will be updated
!  to include the subroutine late. Lei Shi, 03/25/09
!-----------------------------------------------------------------------
!
          k=N(ng)
!!!  CALL CO2_flux (Istr, Iend, LBi, UBi, LBj, UBj,                    &
!!! &                     IminS, ImaxS, j, DoNewton,                   &
!!!#  ifdef MASKING
!!! &                     rmask,                                       &
!!!#  endif
!!! &                     Bio(IminS:,k,itemp), Bio(IminS:,k,isalt),    &
!!! &                     Bio(IminS:,k,iTIC_), Bio(IminS:,k,iTAlk),    &
!!! &                     kw660, 0_r8,                                 &
!!! &                     0.0_r8, 0.0_r8, pH, xco2_in, pCO2, co2flx)
              DO i=Istr,Iend
      co2flx(i)=0.0_r8
	bio(i,k,iTIC_)=bio(i,k,iTIC_)+dtdays*co2flx(i)*Hz_inv(i,k)
              END DO
#endif
#ifdef CARBON
!
!-----------------------------------------------------------------------
!     adjust the alkalinity
!-----------------------------------------------------------------------
!
              DO i=Istr,Iend
           DO k=1,N(ng)
 	cff0=Bio_bak(i,k,iNO3_)-Bio(i,k,iNO3_)-                         &
     &       (Bio_bak(i,k,iNH4_)-Bio(i,k,iNH4_))
	bio(i,k,iTAlk)=bio(i,k,iTAlk)+cff0
           END DO
              END DO
#endif
!
!-----------------------------------------------------------------------
!     CALCULATING THE SINKING FLUX
!-----------------------------------------------------------------------
!
          SINK_LOOP: DO isink=1,Nsink
            indx=idsink(isink)
              DO i=Istr,Iend
           DO k=1,N(ng)
      thick=HZ(i,j,k)
	cff0=HZ(i,j,k)/dtdays
	cff1=min(0.9_r8*cff0,wbio(isink))
        if(k.eq.N(ng))then
          sinkindx(k) = cff1*Bio(i,k,indx)/thick
        elseif(k.gt.1.and.k.lt.n(ng))then
          sinkindx(k) = cff1*                                           &
     &                 (Bio(i,k,indx)-Bio(i,k+1,indx))/thick
        elseif(k.eq.1)then
          sinkindx(k) = cff1*                                           &
     &                 (-Bio(i,k+1,indx))/thick
        endif
            END DO
           DO k=1,N(ng)
	bio(i,k,indx)=bio(i,k,indx)-dtdays*sinkindx(k)
	bio(i,k,indx)=max(bio(i,k,indx),0.00001_r8)
           END DO
              END DO
        END DO SINK_LOOP
      END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables.
!-----------------------------------------------------------------------
!
#ifdef CARBON
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,iTIC_)=MIN(Bio(i,k,iTIC_),3000.0_r8)
            Bio(i,k,iTIC_)=MAX(Bio(i,k,iTIC_),400.0_r8)
          END DO
        END DO
#endif

        DO ibio=1,NBT
          indx=idbio(ibio)
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,j,k,nnew,indx)=MIN(t(i,j,k,nnew,indx),0.0_r8)+        &
     &                           Hz(i,j,k)*Bio(i,k,indx)
#ifdef TS_MPDATA
              t(i,j,k,3,indx)=t(i,j,k,nnew,indx)*Hz_inv(i,k)
#endif
            END DO
          END DO
        END DO
      END DO J_LOOP

      RETURN
      END SUBROUTINE biology_tile

#ifdef OXYGEN
      SUBROUTINE O2_flux (Istr, Iend,                                   &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, j,                           &
#  ifdef MASKING
     &                       rmask,                                     &
#  endif
     &                       T, S, O2, kw660, ppo, o2sat, O2flx)
!
!***********************************************************************
!                                                                      !
!  This routine computes equilibrium partial pressure of CO2 (pCO2)    !
!  in the surface seawater.                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     LBj        J-dimension lower bound.                              !
!     UBj        J-dimension upper bound.                              !
!     IminS      I-dimension lower bound for private arrays.           !
!     ImaxS      I-dimension upper bound for private arrays.           !
!     j          j-pipelined index.                                    !
!     rmask      Land/Sea masking.                                     !
!     T          Surface temperature (Celsius).                        !
!     S          Surface salinity (PSS).                               !
!     O2         Dissolevd oxygen concentration (micromole O2/m^3)     !
!     kw660      gas transfer velocity at a Schmidt number of 660,     !
!                  accounting for sea ice fraction (cm/hr)             !
!     ppo        surface pressure divided by 1 atm.                    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     o2sat      dissolved oxygen saturation concentration (mmol/m^3)  !
!                  due to air-sea exchange (mmol/m^2/day)              !
!     o2flx      time rate of oxygen O2 flux in the sea surface        !
!                  due to air-sea exchange (mmol/m^2/day)              !
!                                                                      !
!  This subroutine was modified from OCMIP2 code. Lei Shi 03/25/2009.  !
!                                                                      !
!***********************************************************************
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j
!
      real(r8),  intent(in) :: ppo
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: T(IminS:)
      real(r8), intent(in) :: S(IminS:)
      real(r8), intent(in) :: O2(IminS:)
#  else
#   ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: T(IminS:ImaxS)
      real(r8), intent(in) :: S(IminS:ImaxS)
      real(r8), intent(in) :: O2(IminS:ImaxS)
#  endif

      real(r8), intent(in) :: kw660(IminS:ImaxS)

      real(r8), intent(out) :: o2sat(IminS:ImaxS)
      real(r8), intent(out) :: o2flx(IminS:ImaxS)
!
!  Local variable declarations.
!

      integer :: i

      real(r8) :: sco2, kwo2
      real(r8) :: TT, TK, TS, TS2, TS3, TS4, TS5, CO

      real(r8), parameter :: A0 = 2.00907_r8       ! Oxygen
      real(r8), parameter :: A1 = 3.22014_r8       ! saturation
      real(r8), parameter :: A2 = 4.05010_r8       ! coefficients
      real(r8), parameter :: A3 = 4.94457_r8
      real(r8), parameter :: A4 =-0.256847_r8
      real(r8), parameter :: A5 = 3.88767_r8
      real(r8), parameter :: B0 =-0.00624523_r8
      real(r8), parameter :: B1 =-0.00737614_r8
      real(r8), parameter :: B2 =-0.0103410_r8
      real(r8), parameter :: B3 =-0.00817083_r8
      real(r8), parameter :: C0 =-0.000000488682_r8
!
!=======================================================================
!  Determine coefficients.  If land/sea
!  masking, compute only on water points.
!=======================================================================
!
      I_LOOP: DO i=Istr,Iend
#  ifdef MASKING
        IF (rmask(i,j).gt.0.0_r8) THEN
#  endif
!
! ********************************************************************
!                                     
! Computes the oxygen saturation concentration at 1 atm total pressure
! in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
! in permil). 
!
! FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
! THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
!
! o2sato IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
! 0 permil <= S <= 42 permil
! C
! CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil, 
! o2sato = 282.015 mmol/m^3
!
! ********************************************************************
!
      TT  = 298.15_r8-T(i)
      TK  = 273.15_r8+T(i)
      TS  = LOG(TT/TK)
      TS2 = TS**2
      TS3 = TS**3
      TS4 = TS**4
      TS5 = TS**5
      CO  = A0 + A1*TS + A2*TS2 + A3*TS3 + A4*TS4 + A5*TS5              &
     &     + S(i)*(B0 + B1*TS + B2*TS2 + B3*TS3)                        &
     &     + C0*(S(i)*S(i))
      o2sat(i) = EXP(CO)
!
!  Convert from ml/l to mol/m^3
!
      o2sat(i) = (o2sat(i)/22391.6_r8)*1000.0_r8
!
!  Convert from mol/m^3 to mmol/m^3
!
      o2sat(i) = o2sat(i)*1000.0_r8
!
!
!*********************************************************************
!
!  Computes the Schmidt number of oxygen in seawater using the
!  formulation proposed by Keeling et al. (1998, Global Biogeochem.
!  Cycles, 12, 141-163).  Input is temperature in deg C.
!
!*********************************************************************
!
      sco2 = 1638.0_r8 - 81.83_r8*t(i) +                                &
     &       1.483_r8*t(i)**2 - 0.008004_r8*t(i)**3

!
!  Compute the transfer velocity for O2 in m/s
!
!    kwo2 = Kw660 * (sco2(t)/660)**-0.5*0.01/3600.0    !(in  m/sec)
!
      kwo2 = Kw660(i) * sqrt(660.0_r8/sco2)   !(in units of cm/hr)
!
!  Compute the transfer velocity for O2 in m/day
!
      KWO2=KWO2*0.01_r8*24.0_r8
!
!  (in units of m/day)
!
!  Compute the saturation concentrations for O2
!
!      o2sat(i) = o2sato(t,s)*ppo       
!  OCMIP
!      o2sat = dosat(t+273.15,s)     
!  Weiss
!
!  Compute time rate of O2 gas exchange
!
      o2flx(i) = kwo2*(o2sat(i)*ppo-o2(i))
#  ifdef MASKING
      ELSE
        o2sat(i)=0.0_r8
        o2flx(i)=0.0_r8
      END IF
#  endif

      END DO I_LOOP

      RETURN
      END SUBROUTINE O2_flux
# endif


