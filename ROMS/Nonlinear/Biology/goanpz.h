#include "cppdefs.h"
      SUBROUTINE biology (ng,tile)
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  then the global biological fields.                                  !
!                                                                      !
!  Sarah Hinckley''s GOANPZ Code                                        !
!  Implemented by Craig Lewis (CVL)                                    !
!  Modified by Liz Dobbins and Sarah Hinckley                          !
!                                                                      !
!=======================================================================
!
#define        NEWSHADE    /* Use Craig''s formulation for self shading in PAR calc
                       Else use Sarah''s self-shading from original NPZ code */
#undef        KODIAK_IRAD /* Generate irradiance with curve matching Kodiak data
                       Else use shortwave radiation (srflx) as irradiance   */
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping

      integer, intent(in) :: ng, tile

# include "tile.h"
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
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   N(ng), NT(ng), NTS(ng),                        &
     &                   nnew(ng), nstp(ng),                            &
# ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
     &                   FORCES(ng) % srflx,                            &
     &                   OCEAN(ng) % t,                                 &
     &                   OCEAN(ng) % st)
#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15)
#endif
      RETURN
      END SUBROUTINE biology
!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         UBk, UBt, UBst,                          &
     &                         nnew, nstp,                              &
# ifdef MASKING
     &                         rmask,                                   &
# endif
     &                         Hz, z_r, z_w, srflx, t, st)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_scalars
      USE mod_ocean
      USE mod_grid
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: UBk, UBt, UBst
      integer, intent(in) :: nnew, nstp
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: st(LBi:,LBj:,:,:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
      real(r8), intent(inout) :: st(LBi:UBi,LBj:UBj,UBk,3,UBst)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, ibio, itr, itrmx, itrc
      integer :: Iter
      integer :: iday, month, year

      real(r8) :: cff1, cff2, cff3
      real(r8) :: Drate, Pmax, NOup, NHup
      real(r8) :: dtdays
      real(r8) :: LightLim, NOLim, NHLim, IronLim
      real(r8) :: hour, yday, lat, k_phy, Dl, Par1
      real(r8) :: Sal1, Temp1, TmaxPhS, KtBm_PhS, TmaxPhL, KtBm_PhL
      real(r8) :: ParMax,TmaxMZS,KtBm_MZS,TmaxMZL,KtBm_MZL,BasalMet
      real(r8) :: Iron1, kfePh,respPh
      real(r8) :: PON,Pv0,PvT,Dep1,Nitrif,NH4R

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: DBio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_bak
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Prod
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS) :: PARsur
      real(r8), dimension(IminS:ImaxS,N(ng)) :: PAR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Dens
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TestVal
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncPhS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncPhL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncMZS
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncMZL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncCop
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncNeo
      real(r8), dimension(IminS:ImaxS,N(ng)) :: TempFuncEup
      real(r8), dimension(IminS:ImaxS,N(ng)) :: HzL
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: z_wL
#ifdef DIAPAUSE
      logical :: downward = .false., upward = .false.
#endif
!
      real(r8), parameter :: eps  = 1.0E-20_r8
      real(r8), parameter :: minv = 0.0E-20_r8
!----------------------
#include "set_bounds.h"
!
      CALL caldate (r_date, tdays(ng), year, yday, month, iday, hour)
      dtdays = dt(ng)*sec2day/REAL(BioIter(ng),r8)
      k_phy = k_chl / ccr
!
#ifdef DIAPAUSE
!  Based on date, determine if NCa are going upwards or downwards
      downward = .false.
      upward = .false.
      IF ( ( RiseStart.lt.RiseEnd .and.                                 &
     &       yday.ge.RiseStart .and. yday.le.RiseEnd ) .or.             &
     &     ( RiseStart.gt.RiseEnd .and.                                 &
     &      ( yday.ge.RiseStart .or. yday.le.RiseEnd ) ) )  THEN
        upward = .true.
      ELSE IF ( ( SinkStart.lt.SinkEnd .and.                            &
     &       yday.ge.SinkStart .and. yday.le.SinkEnd ) .or.             &
     &     ( SinkStart.gt.SinkEnd .and.                                 &
     &      ( yday.ge.SinkStart .or. yday.le.SinkEnd ) ) )  THEN
        downward = .true.
      END IF
#endif
!
! ----------------------------------------------------------------------
! Begin HORIZONTAL INDEX LOOPING
! ----------------------------------------------------------------------
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
!  Extract biological variables from tracer arrays, and
!  restrict their values to be positive definite.  Removed CVL''s
!  conservation of mass correction because conflicted with SPLINES.
!  For ROMS 2.2+, convert from "flux form" to concentrations by
!  dividing by grid cell thickness.
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_bak(i,k,ibio)=t(i,j,k,nstp,ibio)
              Bio(i,k,ibio)=Bio_bak(i,k,ibio)
              Prod(i,k,ibio)=0.0_r8
              DBio(i,k,ibio)=0.0_r8
            END DO
          END DO
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=t(i,j,k,nstp,itemp)
            Bio(i,k,isalt)=t(i,j,k,nstp,isalt)
            IF (Bio(i,k,itemp) .gt. 35._r8) THEN
              print *, 'Temperature: ', &
     &             Bio(i,k,itemp),i, j, k,ng, yday
              print *,'Tracer: ',t(i,j,k,1,itemp),t(i,j,k,2,itemp), &
     &              t(i,j,k,3,itemp),Hz(i,j,k),nnew
              print *,'Others: ', z_w(i,j,N(ng)),       &
     &                  GRID(ng) % h(i,j)
            END IF
            IF ((grid(ng) % h(i,j) + z_w(i,j,N(ng))) .lt.0.0_r8) THEN
              print *, 'zeta & h: ',z_w(i,j,N(ng)),'  ',grid(ng) % h(i,j)
            END IF
          END DO
        END DO
!
! ----------------------------------------------------------------------
!  Calculate Day Length and Surface PAR
! ----------------------------------------------------------------------
!
!  Calculate Day Length
          DO i=Istr,Iend
#if defined DIURNAL_SRFLUX
!  Day Length is already accounted for in ANA_SWRAD so disable correction
            Dl = 24.0_r8
#else
!  Day Length calculation (orig from R. Davis) from latitude and declination.
!  cff2 is Solar declination from Oberhuber (1988) (COADS documentation)
!            lat = 58.00  orig from C code
!            lat = 45.00_r8  test for EPOC
            lat = GRID(ng) % latr(i,j)
            cff1 = 2.0_r8 * pi * ( yday-1.0_r8 ) / 365.0_r8
            cff2 = 0.006918_r8 - 0.399912_r8*cos(cff1)                  &
     &           + 0.070257_r8*sin(cff1) - 0.006758_r8*cos(2*cff1)      &
     &           + 0.000907_r8*sin(2*cff1) - 0.002697_r8*cos(3*cff1)    &
     &           + 0.00148_r8*sin(3*cff1)
            cff3 = lat * pi /180.0_r8
            IF ( abs( -tan(cff3)*tan(cff2) ) .le. 1.0_r8 ) THEN
              cff1 = acos( -tan(cff3)*tan(cff2) ) * 180.0_r8 / pi
              Dl = 2.0_r8 / 15.0_r8 * cff1
            ELSE
              IF ( yday.gt.90.0_r8 .and. yday.lt.270.0_r8 ) THEN
                Dl = 24.0_r8
              ELSE
                Dl = 0.0_r8
              END IF
            END IF
#endif
!  Calculate PAR at the surface
#ifdef KODIAK_IRAD
!  For PAR, Eyeball fit of data from Hinckley''s ezeroday.dat (E d-1 m-2)
            cff2 = 41.0_r8 - 35.0_r8                                   &
     &           * COS( ( 12.0_r8 + yday) * 2.0_r8 * pi / 365.0_r8 )
#else
!  For PAR, use Shortwave radiation ( = surface solar irradiance)
!  converted from deg C m/s to E/m2/day
            cff2 = srflx(i,j) * rho0 * Cp * 0.394848_r8
#endif
           PAR(i,N(ng)) = PARfrac(ng) * cff2                                &
     &      * exp( k_ext + k_chl*(Bio(i,N(ng),iPhS)/ccr +               &
     &                            Bio(i,N(ng),iPhL)/ccrPhL)**0.428      &
     &      * ( z_r(i,j,N(ng)) -  z_w(i,j,N(ng)) ) )
          END DO
!  Calculate light decay in the water column
#ifdef NEWSHADE
!-----------------------------------------------------------
!  George Blamey''s version after Morel 1988 (in Loukos 1977)
!-----------------------------------------------------------
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              cff1 = k_ext * ( z_r(i,j,k) - z_r(i,j,k+1) )
              cff2 = (k_chl*(Bio(i,k+1,iPhS)/ccr +                      &
     &                       Bio(i,k+1,iPhL)/ccrPhL)**0.428_r8)         &
     &                 * ( z_w(i,j,k) - z_r(i,j,k+1) )
              cff3 = (k_chl*(Bio(i,k,iPhS)/ccr +                        &
     &                       Bio(i,k,iPhL)/ccrPhL)**0.428_r8)           &
                       * ( z_r(i,j,k) - z_w(i,j,k) )
              PAR(i,k) = PAR(i,k+1) * EXP(cff1+cff2+cff3)
            END DO
          END DO
#else
!  Version from Sarah''s old C code (probably wrong)
          DO k=N(ng),1,-1
            DO i=Istr,Iend
              cff3 = z_r(i,j,k)+2.5_r8
              IF ( cff3 .gt. -71.0_r8 ) THEN
                cff1 = k_ext + k_chl *                                  &
     &                  ( Bio(i,k,iPhS) + Bio(i,k,iPhL) ) / ccr
              ELSE
                cff1 = 0.077_r8
              END IF
              PAR(i,k) = PARfrac(ng) * cff2 * exp( cff1 * cff3 )
            END DO
          END DO
#endif
          TmaxPhS = 20_r8
          KtBm_PhS = 0.069_r8
          TmaxPhL = 20
          KtBm_PhL = 0.069
          TmaxMZS = 20
          KtBm_MZS = 0.069
          TmaxMZL = 20
          KtBm_MZL = 0.069
          DO k=1,N(ng)
            DO i=Istr,Iend
              HzL(i,k) = Hz(i,j,k)
              Sal1 = Bio(i,k,isalt)
              Temp1 = Bio(i,k,itemp)
              !------------------------------
              !Compute sigma-t for each depth
              !------------------------------
              Dens(i,k) = ComputeDensity(Temp1,Sal1)
              !---------------------------------
              !Arhonditsis temperature functions
              !---------------------------------
              TempFuncPhS(i,k) = GetPhytoResp2(Temp1,TmaxPhS,           &
     &              KtBm_PhS)
              TempFuncPhL(i,k) = GetPhytoResp2(Temp1,TmaxPhL,           &
     &              KtBm_PhL)
              TempFuncMZS(i,k) = GetPhytoResp2(Temp1,TmaxMZS,           &
     &              KtBm_MZS)
              TempFuncMZL(i,k) = GetPhytoResp2(Temp1,TmaxMZL,           &
     &              KtBm_MZL)
              !-------------------------------
              ! Copepod respiration correction
              !-------------------------------
              TempFuncCop(i,k) = GetCopepodResp(Temp1,respCop)
              !----------------------------------
              ! Neocalanus respiration correction
              !----------------------------------
              TempFuncNeo(i,k) = GetCopepodResp(Temp1,respNCa)
              !----------------------------------
              ! Euphausiid respiration correction
              ! Use the copepod function
              !----------------------------------
              TempFuncEup(i,k) = GetCopepodResp(Temp1,respEup)
            END DO
          END DO
          DO k=0,N(ng)
            DO i=Istr,Iend
              z_wL(i,k) = z_w(i,j,k)
            END DO
          END DO
!
! ----------------------------------------------------------------------
! Begin BIOITER LOOP
! ----------------------------------------------------------------------
!
          ITER_LOOP: DO Iter=1,BioIter(ng)
!               if ( .not. downward) then
                  !------------------------------------
                  !Make Neocalanus go down if temp > 12
                  !------------------------------------
!                 if (Bio(i,k,itemp) .gt. 12._r8 .and. NCa(k) .gt. 0.2_r8) then
!                    downward = .true.
!                   goto 111
!                     end if
!               end if
! 111           Continue
               LightLim = 1.0_r8
               NOLim = 1.0_r8
               NHLim = 1.0_r8
               IronLim = 1.0_r8
!=======================================================================
!  Uptake by Small Phytoplankton
!--------------------------------
            DO k=1,N(ng)
              DO i=Istr,Iend
                !------------------------
                !Growth rate computations
                !------------------------
                Drate = DiS * 10.0_r8 ** (DpS * Bio(i,k,itemp) )
                Pmax = (2.0_r8 ** Drate - 1.0_r8 )
!               old LightLim = TANH( alphaPhS * PAR(k) / Pmax / ccr )
                !------------------
                !Nitrate limitation
                !------------------
                NOLim = Bio(i,k,iNO3) * EXP( -psiPhS * Bio(i,k,iNH4) )  &
     &                    / ( k1PhS + Bio(i,k,iNO3) )
#ifdef IRON_LIMIT
!--------------------------------------------------------------
!  Define Iron Limitation such that it is disabled at high iron
!  concentration (2 micromol Fe m-3)
!--------------------------------------------------------------
                IronLim = eps + Bio(i,k,iFe) / (kfePhS + Bio(i,k,iFe))  &
     &                            * (kfePhS + 2._r8) / 2._r8
#endif
!--------------------------------------------------------------
                !Light limitation function
                !-------------------------
                Par1 = Par(i,k)
                LightLim = GetLightLimIronSml(alphaPhS, PAR1,        &
     &                  Pmax,ccr, IronLim)
                !--------------
                !Nitrate uptake
                !--------------
                NOup = Bio(i,k,iPhS) * Pmax * LightLim * NOLim * IronLim
                !-------------------------------
                !Change in nitrate concentration
                !-------------------------------
                DBio(i,k,iNO3) = DBio(i,k,iNO3) - xi * NOup * dtdays
                !-------------------
                !Ammonium limitation
                !-------------------
                NHLim = Bio(i,k,iNH4) / ( k2PhS + Bio(i,k,iNH4) )
                !-----------------------------
                !Light limitation for ammonium
                !-----------------------------
                LightLim = GetLightLimSml(alphaPhS, PAR1, Pmax,    &
     &                  ccr)
                !---------------
                !Ammonium uptake
                !---------------
                NHup = Bio(i,k,iPhS) * Pmax * LightLim * NHLim
                !--------------------------------
                !Change in ammonium concentration
                !--------------------------------
                DBio(i,k,iNH4) = DBio(i,k,iNH4) - xi * NHup * dtdays
                !----------------------------------------------
                !Change in concentration of small phytoplankton
                !----------------------------------------------
                DBio(i,k,iPhS) = DBio(i,k,iPhS) + ( NOup + NHup ) * dtdays
                !-----------------------------------------
                !Primary production of small phytoplankton
                !-----------------------------------------
                Prod(i,k,iPhS) = Prod(i,k,iPhS) + DBio(i,k,iPhS)
#ifdef IRON_LIMIT
                !----------------------------
                !Change in iron concentration
                !----------------------------
                DBio(i,k,iFe) = DBio(i,k,iFe) - FeC * NOup * dtdays
#endif
#if defined BIOFLUX && defined GAK1D
                IF (i.eq.3.and.j.eq.3) THEN
                  bflx(iNO3,iPhS) = bflx(iNO3,iPhS) + NOup
                  bflx(iNH4,iPhS) = bflx(iNH4,iPhS) + NHup
                ENDIF
#endif
              END DO
            END DO
!=========================================================================
!  Uptake by Large Phytoplankton
!-------------------------------
            DO k=1,N(ng)
              DO i=Istr,Iend
                !------------------------
                !Growth rate computations
                !------------------------

                Drate = DiL * 10.0_r8 ** (DpL * Bio(i,k,itemp) )
                Pmax = (2.0_r8 ** Drate - 1.0_r8 )
                !------------------
                !Nitrate limitation
                !------------------
                NOLim = Bio(i,k,iNO3) * EXP( -psiPhL * Bio(i,k,iNH4) )  &
     &                    / ( k1PhL + Bio(i,k,iNO3) )
#ifdef IRON_LIMIT
!----------------------------------------------------------------
!  Define Iron Limitation such that it is disabled at the highest
!  iron concentration (2 micromol Fe m-3)
!----------------------------------------------------------------
                IronLim = eps + Bio(i,k,iFe) / (kfePhL + Bio(i,k,iFe))  &
     &                           * (kfePhL + 2._r8) / 2._r8
#endif
!----------------------------------------------------------------
                !Light limitation function
                !-------------------------
                Par1 = Par(i,k)
                ParMax = Par(i,N(ng))
                !----------------------------------------------
                !Use composite light curve with iron limitation
                !----------------------------------------------
                LightLim = GetLightLimIron(alphaPhL,PAR1,Pmax,   &
     &                  ccrPhL,IronLim,ParMax)
!                   LightLim = GetLightLimIron2(alphaPhL, PAR1, Pmax, &
!     &                  ccrPhL, IronLim)
                !--------------
                !Nitrate uptake
                !--------------
                NOup = Bio(i,k,iPhL) * Pmax * LightLim * NOLim * IronLim
                !-------------------------------
                !Change in nitrate concentration
                !-------------------------------
                DBio(i,k,iNO3) = DBio(i,k,iNO3) - xi * NOup * dtdays
                !-------------------
                !Ammonium limitation
                !-------------------
                NHLim = Bio(i,k,iNH4) / ( k2PhL + Bio(i,k,iNH4) )
                !-------------------------------------------------
                !Use composite light curve without iron limitation
                !-------------------------------------------------
                LightLim = GetLightLim(alphaPhL,PAR1,Pmax,   &
     &                  ccrPhL,ParMax)
                !----------------------------
                !Use hyperbolic tangent curve
                !----------------------------
!                   LightLim = GetLightLim2(alphaPhL, PAR1, Pmax, ccrPhL)
                !---------------
                !Ammonium uptake
                !---------------
                NHup = Bio(i,k,iPhL) * Pmax * LightLim * NHLim
                !--------------------------------
                !Change in ammonium concentration
                !--------------------------------
                DBio(i,k,iNH4) = DBio(i,k,iNH4) - xi * NHup * dtdays
                !----------------------------------------------
                !Change in concentration of large phytoplankton
                !----------------------------------------------
                DBio(i,k,iPhL) = DBio(i,k,iPhL) + ( NOup + NHup ) * dtdays
                !-----------------------------------------
                !Primary production of large phytoplankton
                !-----------------------------------------
                Prod(i,k,iPhL) = Prod(i,k,iPhL) + DBio(i,k,iPhL)
#ifdef IRON_LIMIT
                !----------------------------
                !Change in iron concentration
                !----------------------------
                DBio(i,k,iFe) = DBio(i,k,iFe) - FeC * NOup * dtdays
#endif
#if defined BIOFLUX && defined GAK1D
                  IF (i.eq.3.and.j.eq.3) THEN
                    bflx(iNO3,iPhL) = bflx(iNO3,iPhL) + NOup
                    bflx(iNH4,iPhL) = bflx(iNH4,iPhL) + NHup
                  END IF
#endif
               END DO
             END DO

!=======================================================================
! Grazing by MZS
!----------------
            DO k=1,N(ng)
              DO i=Istr,Iend
                 !----------------
                 !Food preferences
                 !----------------
                  cff1 = fpPhSMZS * Bio(i,k,iPhS)**2                    &
     &                 + fpPhLMZS * Bio(i,k,iPhL)**2
                 !--------------------------------------
                 !Temperature corrected food consumption
                 !--------------------------------------
                  cff2 = eMZS *                                         &
     &              Q10MZS ** ( (Bio(i,k,itemp)-Q10MZST) / 10.0_r8 )    &
     &                 * Bio(i,k,iMZS) / (fMZS**2 + cff1)
                 !-------------------------------------------------
                 !Change in ammonium concentration due to excretion
                 !Computed from respiration
                 !-------------------------------------------------
!                 DBio(i,k,iNH4) = DBio(i,k,iNH4) +                     &
!      & xi * kMZS * cff1 * cff2 * dtdays
                 !--------------------------------------------------------
                 !Change in small and large phytoplankton due to predation
                 !--------------------------------------------------------
                  DBio(i,k,iPhS) = DBio(i,k,iPhS) - fpPhSMZS *          &
     &                (Bio(i,k,iPhS)**2) * cff2 * dtdays
                  DBio(i,k,iPhL) = DBio(i,k,iPhL) - fpPhLMZS *          &
     &                (Bio(i,k,iPhL)**2) * cff2 * dtdays
                 !--------------------------------
                 !Growth of small microzooplankton
                 !--------------------------------
                  DBio(i,k,iMZS) = DBio(i,k,iMZS) +                     &
     &                         gammaMZS * cff1 * cff2 * dtdays
                 !-------------------------------------
                 !Production for small microzooplankton
                 !-------------------------------------
                  Prod(i,k,iMZS) = Prod(i,k,iMZS) + DBio(i,k,iMZS)
                 !---------------------------
                 ! Additions to detritus pool
                 ! kMZS not needed
                 !---------------------------
                  DBio(i,k,iDet) = DBio(i,k,iDet) +                     &
     &                 (1.0_r8 - gammaMZS) * cff1 * cff2 * dtdays
#if defined BIOFLUX && defined GAK1D
                  IF (i.eq.3.and.j.eq.3) THEN
                    bflx(iMZS,iNH4) = bflx(iMZS,iNH4) + kMZS*cff1*cff2
                    bflx(iPhS,iMZS) = bflx(iPhS,iMZS) +                 &
                      fpPhSMZS*PhS(k)*cff2
                    bflx(iPhL,iMZS) = bflx(iPhL,iMZS) +                 &
                      fpPhSMZL*PhL(k)*cff2
                    bflx(iMZS,iDet) = bflx(iMZS,iDet)  +                &
     &                ( 1.0_r8-kMZS-gammaMZS )*cff1*cff2
                  END IF
#endif
!
               END DO
             END DO

!========================================================================
! Grazing by MZL
!---------------
            DO k=1,N(ng)
              DO i=Istr,Iend
                !----------------
                !Food preferences
                !----------------
                cff1 = fpPhSMZL * Bio(i,k,iPhS)**2                      &
     &               + fpPhLMZL * Bio(i,k,iPhL)**2                      &
     &               + fpMZSMZL * Bio(i,k,iMZS)**2
                 !--------------------------------------
                 !Temperature corrected food consumption
                 !--------------------------------------
                 cff2 = eMZL *                                          &
     &              ( Q10MZL ** ( (Bio(i,k,itemp)-Q10MZLT) / 10.0_r8 ) )&
     &                 * Bio(i,k,iMZL) / (fMZL**2 + cff1)
                 !-------------------------------------------------
                 !Change in ammonium concentration due to excretion
                 !Computed from respiration
                 !-------------------------------------------------
!                 DBio(i,k,iNH4) = DBio(i,k,iNH4) +                     &
!     &                    xi * kMZL * cff1 * cff2 * dtdays
                 !--------------------------------------------------------
                 !Change in small and large phytoplankton due to predation
                 !--------------------------------------------------------
                 DBio(i,k,iPhS) = DBio(i,k,iPhS) - fpPhSMZL *           &
     &                 (Bio(i,k,iPhS)**2) * cff2 * dtdays
                 DBio(i,k,iPhL) = DBio(i,k,iPhL) - fpPhLMZL *           &
     &                 (Bio(i,k,iPhL)**2) * cff2 * dtdays
                 DBio(i,k,iMZS) = DBio(i,k,iMZS) - fpMZSMZL *           &
     &                 (Bio(i,k,iMZS)**2) * cff2 * dtdays
                 !--------------------------------
                 !Growth of large microzooplankton
                 !--------------------------------
                 DBio(i,k,iMZL) = DBio(i,k,iMZL) +                      &
     &                        gammaMZL * cff1 * cff2 * dtdays
                 !------------------------------------
                 !Production of large microzooplankton
                 !------------------------------------
                 Prod(i,k,iMZL) = Prod(i,k,iMZL) + DBio(i,k,iMZL)
                 !---------------------------
                 ! Additions to detritus pool
                        ! kMZL not needed
                 !---------------------------
                 DBio(i,k,iDet) = DBio(i,k,iDet) +                      &
     &                 (1.0_r8 - gammaMZL) * cff1 * cff2 * dtdays
#if defined BIOFLUX && defined GAK1D
                  IF (i.eq.3.and.j.eq.3) THEN
                    bflx(iMZL,iNH4) = bflx(iMZL,iNH4) + kMZL*cff1*cff2
                    bflx(iPhS,iMZL) = bflx(iPhS,iMZL) +                 &
                      fpPhSMZL*PhS(k)*cff2
                    bflx(iPhL,iMZL) = bflx(iPhL,iMZL) +                 &
                      fpPhLMZL*PhL(k)*cff2
                    bflx(iMZL,iDet) = bflx(iMZL,iDet)  +                &
     &                ( 1.0_r8-kMZL-gammaMZL )*cff1*cff2
                  END IF
#endif
               END DO
             END DO

!===========================================================================
! Grazing and Predation by Copepods
!-----------------------------------
            DO k=1,N(ng)
              DO i=Istr,Iend
                !----------------
                !Food preferences
                !----------------
                cff1 = fpPhSCop * Bio(i,k,iPhS)**2                      &
     &               + fpPhLCop * Bio(i,k,iPhL)**2                      &
     &               + fpMZSCop * Bio(i,k,iMZS)**2                      &
     &               + fpMZLCop * Bio(i,k,iMZL)**2
                 !--------------------------------------
                 !Temperature corrected food consumption
                 !--------------------------------------
                 cff2 = eCop *                                          &
     &             (Q10Cop ** ( (Bio(i,k,itemp)-Q10CopT) / 10.0_r8 ))   &
     &                 * Bio(i,k,iCop) / (fCop**2 + cff1)
                 !------------------------
                 !Growth of small copepods
                 !------------------------
                 DBio(i,k,iCop) = DBio(i,k,iCop) +                      &
     &                        gammaCop * cff1 * cff2 * dtdays
                 !------------------
                 !Copepod production
                 !------------------
                 Prod(i,k,iCop) = Prod(i,k,iCop) + DBio(i,k,iCop)
                 !----------------------------------------------
                 !Changes in prey concentration due to predation
                 !----------------------------------------------
                 DBio(i,k,iPhS) = DBio(i,k,iPhS) -                      &
     &                 fpPhSCop * (Bio(i,k,iPhS)**2) * cff2 * dtdays
                 DBio(i,k,iPhL) = DBio(i,k,iPhL) -                      &
     &                 fpPhLCop * (Bio(i,k,iPhL)**2) * cff2 * dtdays
                 DBio(i,k,iMZS) = DBio(i,k,iMZS) -                      &
     &                 fpMZSCop * (Bio(i,k,iMZS)**2) * cff2 * dtdays
                 DBio(i,k,iMZL) = DBio(i,k,iMZL) -                      &
     &                 fpMZLCop * (Bio(i,k,iMZL)**2) * cff2 * dtdays
                 !-------------------------------------------------
                 !Change in ammonium concentration due to excretion
                 !Computed from respiration
                 !-------------------------------------------------
!                 DBio(i,k,iNH4) = DBio(i,k,iNH4) +                     &
!     &                         xi * kCop * cff1 * cff2 * dtdays
                 !---------------------------
                 ! Additions to detritus pool
                 ! kCop not needed
                 !---------------------------
                 DBio(i,k,iDet) = DBio(i,k,iDet) +                      &
     &                 (1.0_r8 - gammaCop) * cff1 * cff2 * dtdays
#if defined BIOFLUX && defined GAK1D
                  IF (i.eq.3.and.j.eq.3) THEN
                    bflx(iPhS,iCop)=bflx(iPhS,iCop)+                    &
     &                  fpPhSCop*Bio(i,k,iPhS)*cff2
                    bflx(iPhL,iCop)=bflx(iPhL,iCop)+                    &
     &                  fpPhLCop*Bio(i,k,iPhL)*cff2
                    bflx(iMZS,iCop)=bflx(iMZS,iCop)+                    &
     &                  fpMZSCop*Bio(i,k,iMZS)*cff2
                    bflx(iMZL,iCop)=bflx(iMZL,iCop)+                    &
     &                  fpMZLCop*Bio(i,k,iMZL)*cff2
                    bflx(iCop,iNH4) = bflx(iCop,iNH4) + kCop*cff1*cff2
                    bflx(iCop,iDet) = bflx(iCop,iDet) +                 &
     &                          ( 1.0_r8-kCop-gammaCop )*cff1*cff2
                  END IF
#endif
             END DO
           END DO

!========================================================================
! Grazing and Predation by NCa
!-----------------------------
            DO k=1,N(ng)
              DO i=Istr,Iend
                 !----------------
                 !Food preferences
                 !----------------
                  cff1 = fpPhSNCa * Bio(i,k,iPhS)**2                    &
     &                 + fpPhLNCa * Bio(i,k,iPhL)**2                    &
     &                 + fpMZSNCa * Bio(i,k,iMZS)**2                    &
     &                 + fpMZLNCa * Bio(i,k,iMZL)**2
                        !--------------------------------------
                 !Temperature corrected food consumption
                 !--------------------------------------
                 cff2 = eNCa *                                          &
     &               (Q10NCa ** ( (Bio(i,k,itemp)-Q10NCaT) / 10.0_r8 ) )&
     &                 * Bio(i,k,iNCa) / (fNCa**2 + cff1)
                 !--------------------
                 !Growth of Neocalanus
                 !--------------------
                 DBio(i,k,iNCa) = DBio(i,k,iNCa) +                      &
     &                        gammaNCa * cff1 * cff2 * dtdays
                 !---------------------
                 !Neocalanus production
                 !---------------------
                 Prod(i,k,iNCa) = Prod(i,k,iNCa) + DBio(i,k,iNCa)
                 !----------------------------------------------
                 !Changes in prey concentration due to predation
                 !----------------------------------------------
                 DBio(i,k,iPhS) = DBio(i,k,iPhS) -                      &
     &                 fpPhSNCa * (Bio(i,k,iPhS)**2) * cff2 * dtdays
                 DBio(i,k,iPhL) = DBio(i,k,iPhL) -                      &
     &                 fpPhLNCa * (Bio(i,k,iPhL)**2) * cff2 * dtdays
                 DBio(i,k,iMZS) = DBio(i,k,iMZS) -                      &
     &                 fpMZSNCa * (Bio(i,k,iMZS)**2) * cff2 * dtdays
                 DBio(i,k,iMZL) = DBio(i,k,iMZL) -                      &
     &                 fpMZLNCa * (Bio(i,k,iMZL)**2) * cff2 * dtdays
                 !-------------------------------------------------
                 !Change in ammonium concentration due to excretion
                 !Computed from respiration
                 !-------------------------------------------------
!                 DBio(i,k,iNH4) = DBio(i,k,iNH4) +                     &
!     &                       xi * kNCa * cff1 * cff2 * dtdays
                 !---------------------------
                 ! Additions to detritus pool
                 ! kNCa not needed
                 !---------------------------
                  DBio(i,k,iDet) = DBio(i,k,iDet) +                     &
     &                 (1.0_r8 - gammaNCa) * cff1 * cff2 * dtdays
#if defined BIOFLUX && defined GAK1D
                  IF (i.eq.3.and.j.eq.3) THEN
                    bflx(iPhS,iNCa)=bflx(iPhS,iNCa)+                    &
     &                       fpPhSNCa*Bio(i,k,iPhS)*cff2
                    bflx(iPhL,iNCa)=bflx(iPhL,iNCa)+                    &
     &                       fpPhLNCa*Bio(i,k,iPhL)*cff2
                    bflx(iMZS,iNCa)=bflx(iMZS,iNCa)+                    &
     &                       fpMZSNCa*Bio(i,k,iMZS)*cff2
                    bflx(iMZL,iNCa)=bflx(iMZL,iNCa)+                    &
     &                       fpMZLNCa*Bio(i,k,iMZL)*cff2
                    bflx(iNCa,iNH4) = bflx(iNCa,iNH4) + kNCa*cff1*cff2
                    bflx(iNCa,iDet) = bflx(iNCa,iDet) +                 &
     &                          ( 1.0_r8-kNCa-gammaNCa )*cff1*cff2
                  END IF
#endif
               END DO
             END DO

!=========================================================================
! Grazing and Predation by Euphuasiids
!-------------------------------------
            DO k=1,N(ng)
              DO i=Istr,Iend
                !----------------
                !Food preferences
                !----------------
                cff1 = fpPhSEup * Bio(i,k,iPhS)**2                      &
     &               + fpPhLEup * Bio(i,k,iPhL)**2                      &
     &               + fpMZSEup * Bio(i,k,iMZS)**2                      &
     &               + fpMZLEup * Bio(i,k,iMZL)**2                      &
     &               + fpCopEup * Bio(i,k,iCop)**2
                cff2 = eEup *                                           &
     &                 Q10Eup ** ( (Bio(i,k,itemp)-Q10EupT) / 10.0_r8 ) &
     &                 * Bio(i,k,iEup) / (fEup**2 + cff1)
                !---------------------
                !Growth of Euphausiids
                !---------------------
                DBio(i,k,iEup) = DBio(i,k,iEup) + gammaEup * cff1 * cff2 * dtdays
                !---------------------
                !Euphausiid production
                !---------------------
                Prod(i,k,iEup) = Prod(i,k,iEup) + DBio(i,k,iEup)
                !----------------------------------------------
                !Changes in prey concentration due to predation
                !----------------------------------------------
                DBio(i,k,iPhS) = DBio(i,k,iPhS) -                       &
     &                 fpPhSEup * (Bio(i,k,iPhS)**2) * cff2 * dtdays
                DBio(i,k,iPhL) = DBio(i,k,iPhL) -                       &
     &                 fpPhLEup * (Bio(i,k,iPhL)**2) * cff2 * dtdays
                DBio(i,k,iMZS) = DBio(i,k,iMZS) -                       &
     &                 fpMZSEup * (Bio(i,k,iMZS)**2) * cff2 * dtdays
                DBio(i,k,iMZL) = DBio(i,k,iMZL) -                       &
     &                 fpMZLEup * (Bio(i,k,iMZL)**2) * cff2 * dtdays
                DBio(i,k,iCop) = DBio(i,k,iCop) -                       &
     &                 fpCopEup * (Bio(i,k,iCop)**2) * cff2 * dtdays
                !-------------------------------------------------
                !Change in ammonium concentration due to excretion
                !Computed from respiration
                !-------------------------------------------------
!                DBio(i,k,iNH4) = DBio(i,k,iNH4) +                       &
!     &                      xi * kEup * cff1 * cff2 * dtdays
                !---------------------------
                ! Additions to detritus pool
                ! kEup not needed
                !---------------------------
                DBio(i,k,iDet) = DBio(i,k,iDet) +                                   &
     &                 (1.0_r8 - kEup - gammaEup) * cff1 * cff2 * dtdays
#if defined BIOFLUX && defined GAK1D
                  IF (i.eq.3.and.j.eq.3) THEN
                    bflx(iPhS,iEup)=bflx(iPhS,iEup)+                    &
     &                        fpPhSEup*Bio(i,k,iPhS)*cff2
                    bflx(iPhL,iEup)=bflx(iPhL,iEup)+                    &
     &                        fpPhLEup*Bio(i,k,iPhL)*cff2
                    bflx(iMZS,iEup)=bflx(iMZS,iEup)+                    &
     &                        fpMZSEup*Bio(i,k,iMZS)*cff2
                    bflx(iMZL,iEup)=bflx(iMZL,iEup)+                    &
     &                        fpMZLEup*Bio(i,k,iMZL)*cff2
                    bflx(iCop,iEup)=bflx(iCop,iEup)+                    &
     &                        fpCopEup*Bio(i,k,iCop)*cff2
                    bflx(iEup,iNH4) = bflx(iEup,iNH4) + kEup*cff1*cff2
                    bflx(iEup,iDet) = bflx(iEup,iDet) +                 &
     &                          ( 1.0_r8-kEup-gammaEup )*cff1*cff2
                  END IF
#endif
               END DO
           END DO

!=======================================================================
! Linear Mortality and Senescence Terms, for phytoplankton only
!--------------------------------------------------------------
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff1 = MAX( minmPhS , maxmPhS -                         &
     &                ( maxmPhS - minmPhS ) * Bio(i,k,iNO3) / NcritPhS )
                cff2 = MAX( minmPhL , maxmPhL -                         &
     &                ( maxmPhL - minmPhL ) * Bio(i,k,iNO3) / NcritPhL )
!
                DBio(i,k,iPhS) = DBio(i,k,iPhS) -                       &
     &                       cff1 * Bio(i,k,iPhS) * dtdays
                DBio(i,k,iPhL) = DBio(i,k,iPhL) -                       &
     &                       cff2 * Bio(i,k,iPhL) * dtdays
                  !------------------------------------------
                  !Using nonlinear mortality for copepods and
                  !euphausiids; therefore linear mortality is
                  !commented out
                  !------------------------------------------
!                  DBio(i,k,iMZS) = DBio(i,k,iMZS) -                     &
!      &                         mMZS * Bio(i,k,iMZS) * dtdays
!                  DBio(i,k,iMZL) = DBio(i,k,iMZL) -                     &
!      &                         mMZL * Bio(i,k,iMZL) * dtdays
!                  DBio(i,k,iCop) = DBio(i,k,iCop) -                     &
!      &                         mCop * Bio(i,k,iCop) * dtdays
!                  DBio(i,k,iNCa) = DBio(i,k,iNCa) -                     &
!      &                         mNCa * Bio(i,k,iNCa) * dtdays
!                  DBio(i,k,iEup) = DBio(i,k,iEup) -                     &
!      &                         mEup * Bio(i,k,iEup) * dtdays
                  !-------------------------------------------------
                  !The contribution of microzooplankton, copepod
                  !and euphausiid mortality is computed by non-
                  !linear mortality.  Therefore the linear mortality
                  !contribution below is commented out.
                  !-------------------------------------------------
                  DBio(i,k,iDet) = DBio(i,k,iDet) +                     &
     &                     ( cff1 * Bio(i,k,iPhS) +                     &
     &                       cff2 * Bio(i,k,iPhL) ) * dtdays
!     &                + mMZS * Bio(i,k,iMZS) + mMZL * Bio(i,k,iMZL)     &
!     &                   + mCop * Bio(i,k,iCop) + mNCa * Bio(i,k,iNCa)  &
!     &                   + mEup * Bio(i,k,iEup) ) * dtdays
#if defined BIOFLUX && defined GAK1D
                  IF (i.eq.3.and.j.eq.3) THEN
                    bflx(iPhS,iDet)= bflx(iPhS,iDet) + cff1*Bio(i,k,iPhS)
                    bflx(iPhL,iDet)= bflx(iPhL,iDet) + cff2*Bio(i,k,iPhL)
                    bflx(iMZS,iDet)= bflx(iMZS,iDet) + mMZS*Bio(i,k,iMZS)
                    bflx(iMZL,iDet)= bflx(iMZL,iDet) + mMZL*Bio(i,k,iMZL)
                    bflx(iCop,iDet)= bflx(iCop,iDet) + mCop*Bio(i,k,iCop)
                    bflx(iNCa,iDet)= bflx(iNCa,iDet) + mNCa*Bio(i,k,iNCa)
                    bflx(iEup,iDet)= bflx(iEup,iDet) + mEup*Bio(i,k,iEup)
                  END IF
#endif
             END DO
           END DO

!==================================================================
! Nonlinear Mortality (closure) Terms
!------------------------------------
            DO k=1,N(ng)
              DO i=Istr,Iend
                DBio(i,k,iMZS) = DBio(i,k,iMZS) -                       &
     &                         mpredMZS*dtdays*Bio(i,k,iMZS)**2
                DBio(i,k,iMZL) = DBio(i,k,iMZL) -                       &
     &                         mpredMZL*dtdays*Bio(i,k,iMZL)**2
                DBio(i,k,iCop) = DBio(i,k,iCop) -                       &
     &                         mpredCop*dtdays*Bio(i,k,iCop)**2
                DBio(i,k,iNCa) = DBio(i,k,iNCa) -                       &
     &                         mpredNCa*dtdays*Bio(i,k,iNCa)**2
                DBio(i,k,iEup) = DBio(i,k,iEup) -                       &
     &                         mpredEup*dtdays*Bio(i,k,iEup)**2
! CVL: This code would be needed for mass conservation, but is not in SH
! CVL: version.  I have included to check for closure.
                !---------------------------------
                !Detritus from nonlinear mortality
                !---------------------------------
                DBio(i,k,iDet) = DBio(i,k,iDet) +                       &
     &                    (mpredMZS * Bio(i,k,iMZS)**2                  &
     &                   + mpredMZL * Bio(i,k,iMZL)**2                  &
     &                   + mpredCop * Bio(i,k,iCop)**2                  &
     &                   + mpredNCa * Bio(i,k,iNCa)**2                  &
     &                   + mpredEup * Bio(i,k,iEup)**2 ) * dtdays
#if defined BIOFLUX && defined GAK1D
                IF (i.eq.3.and.j.eq.3) THEN
                    bflx(iCop,itemp)= bflx(iCop,itemp)+                 &
     &                             mpredCop*Bio(i,k,iCop)**2
                    bflx(iNca,itemp)= bflx(iNCa,itemp)+                 &
     &                             mpredNCa*Bio(i,k,iNCak)**2
                    bflx(iEup,itemp)= bflx(iEup,itemp)+                 &
     &                             mpredEup*Bio(i,k,iEupk)**2
                END IF
#endif
              END DO
            END DO

!=========================================================================
!Respiration, temp functions are temperature corrected respiration factors
!computed above by the GetPhytoResp2 and GetCopepodRes functions.
!-------------------------------------------------------------------------
            DO k=1,N(ng)
              DO i=Istr,Iend
                !---------------------------------------------
                !Corrects basal metabolism for iron limitation
                !Small phytoplankton
                !---------------------------------------------
                BasalMet = respPhL
#ifdef IRON_LIMIT
                Iron1 = Bio(i,k,iFe)
                respPh = respPhS
                kfePh = kfePhS
                BasalMet = GetBasalMetabolism(respPh,kfePh,Iron1)
#endif
                !----------------------------------------------------
                !Subtract respiration losses from small phytoplankton
                !Compute net primary production
                !----------------------------------------------------
                DBio(i,k,iPhS) = DBio(i,k,iPhS) -                       &
     &                  TempFuncPhS(i,k)*BasalMet*dtdays*Bio(i,k,iPhS)
                TestVal(i,k) = DBio(i,k,iPhL)
                Prod(i,k,iPHS) = Prod(i,k,iPHS) -                       &
     &               TempFuncPhS(i,k)*BasalMet*dtdays*Bio(i,k,iPhS)
                !---------------------------------------------
                !Corrects basal metabolism for iron limitation
                !large phytoplankton
                !---------------------------------------------
                BasalMet = respPhL
#ifdef IRON_LIMIT
                respPh = respPhL
                kfePh = kfePhL
                BasalMet = GetBasalMetabolism(respPh,kfePh,Iron1)
#endif
               !----------------------------------------------------
               !Subtract respiration losses from large phytoplankton
               !Compute net primary production
               !----------------------------------------------------
               DBio(i,k,iPhL) = DBio(i,k,iPhL) -                        &
     &               TempFuncPhL(i,k)*BasalMet*dtdays*Bio(i,k,iPhL)
                !  TestVal(i,k) = DBio(i,k,iPhL)
               Prod(i,k,iPHL) = Prod(i,k,iPHL) -                        &
     &               TempFuncPhL(i,k)*BasalMet*dtdays*Bio(i,k,iPhL)
               !-------------------------------------------------------
               !Subtract respiration losses from small microzooplankton
               !-------------------------------------------------------
               BasalMet = respMZS
               DBio(i,k,iMZS) = DBio(i,k,iMZS) -                        &
     &                TempFuncMZS(i,k)*BasalMet*dtdays*Bio(i,k,iMZS)
               !----------------------
               !Compute net production
               !----------------------
               Prod(i,k,iMZS) = Prod(i,k,iMZS) -                        &
     &              TempFuncMZS(i,k)*BasalMet*dtdays*Bio(i,k,iMZS)
               !-----------------------------------------------------------
               !Add ammonium to correct for excretion related to metabolism
               !-----------------------------------------------------------
               DBio(i,k,iNH4) = DBio(i,k,iNH4) +                        &
     &               xi*(TempFuncMZS(i,k)*BasalMet*dtdays*Bio(i,k,iMZS))
               !---------------------------------------------
               !Respiration for large microzooplankton is a
               !proxy for predation within the box and can
               !be included as a loss term but additional
               !ammonium regeneration should be added for
               !mass balance
               !---------------------------------------------
               BasalMet = respMZL
               DBio(i,k,iMZL) = DBio(i,k,iMZL) -                        &
     &               TempFuncMZL(i,k)*BasalMet*dtdays*Bio(i,k,iMZL)
               DBio(i,k,iNH4) = DBio(i,k,iNH4) +                        &
     &               xi*(TempFuncMZL(i,k)*BasalMet*dtdays*Bio(i,k,iMZL))
               !----------------------
               !Compute net production
               !----------------------
               Prod(i,k,iMZL) = Prod(i,k,iMZL) -                        &
     &             TempFuncMZL(i,k)*BasalMet*dtdays*Bio(i,k,iMZL)
               !-----------------------------------------------------
               !Respiration, net production and excretion corrections
               !for copeopds and euphausiids
               !-----------------------------------------------------
               DBio(i,k,iCop) = DBio(i,k,iCop) -                        &
     &                        TempFuncCop(i,k)*Bio(i,k,iCop)*dtdays
               TestVal(i,k) = DBio(i,k,iNCa)
               Prod(i,k,iCop) = Prod(i,k,iCop) -                        &
     &            TempFuncCop(i,k)*Bio(i,k,iCop)*dtdays
               DBio(i,k,iNH4) = DBio(i,k,iNH4) +                        &
     &               xi*(TempFuncCop(i,k)*dtdays*Bio(i,k,iCop))
               DBio(i,k,iNCa) = DBio(i,k,iNCa) -                        &
     &                        TempFuncNeo(i,k)*Bio(i,k,iNCa)*dtdays
               Prod(i,k,iNCa) = Prod(i,k,iNCa) -                        &
     &             TempFuncNeo(i,k)*Bio(i,k,iNCa)*dtdays
               DBio(i,k,iNH4) = DBio(i,k,iNH4) +                        &
     &               xi*(TempFuncNeo(i,k)*dtdays*Bio(i,k,iNCa))
               DBio(i,k,iEup) = DBio(i,k,iEup) -                        &
     &                        TempFuncEup(i,k)*Bio(i,k,iEup)*dtdays
               Prod(i,k,iEup) = Prod(i,k,iEup) -                        &
     &              TempFuncEup(i,k)*Bio(i,k,iEup)*dtdays
               DBio(i,k,iNH4) = DBio(i,k,iNH4) +                        &
     &               xi*(TempFuncEup(i,k)*dtdays*Bio(i,k,iEup))
              END DO
            END DO
!=========================================================================
! Molting: NOTE that it is unclear where this molting equation comes from.
! This is present only for euphausiids, not copepods
!-------------------------------------------------------------------------
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff1 = 0.02_r8 / (10.0_r8 - 0.4_r8 * Bio(i,k,itemp))*   &
     &                       Bio(i,k,iEup)
                DBio(i,k,iDet) = DBio(i,k,iDet) + cff1 * dtdays
                DBio(i,k,iEup) = DBio(i,k,iEup) - cff1 * dtdays
#if defined BIOFLUX && defined GAK1D
                IF (i.eq.3.and.j.eq.3) THEN
                  bflx(iEup,iDet)= bflx(iEup,iDet) + cff1
                END IF
#endif
              END DO
            END DO

!=========================================================================
! Detrital Remineralization
!--------------------------
            DO k=1,N(ng)
              DO i=Istr,Iend
                Temp1 = Bio(i,k,itemp)
                !-------------------------
                !Detrital remineralization
                !From Frost (1993).
                !-------------------------
                !cff1 = regen * dgrad * Bio(i,k,iDet)
                !DBio(i,k,iNH4) = DBio(i,k,iNH4) + xi * cff1 * dtdays
                !DBio(i,k,iDet) = DBio(i,k,iDet) - cff1 * dtdays
                !-------------------
                !From Kawamiya(2000)
                !-------------------
                PON = Bio(i,k,iDet)*xi  !Particulate organic nitrogen
                Pv0 = 0.02   !PON decompositon at 0 deg C (d-1)
                PvT = 0.0693 !Temperature coefficient (deg C-1)
                DBio(i,k,iDet) = DBio(i,k,iDet) -                       &
     &                     ((Pv0*exp(PvT*Temp1)*PON)/xi)*dtdays
                DBio(i,k,iNH4) = DBio(i,k,iNH4) +                       &
     &                     ((Pv0*exp(PvT*Temp1)*PON))*dtdays
                !--------------------
                !Nitrification
               !--------------------
               Dep1 = (z_w(i,j,k)) *(-1._r8)
               NH4R = Bio(i,k,iNH4)
               Par1 = PAR(i,k)
               !--------------------
               !Arondritsis & Denmon
               !--------------------
               !Nitrif = GetNitrif(Temp1,Dep1,NH4R)
               !-----------------
               !Fennel & Kawamiya
               !-----------------
               Nitrif = GetNitrif2(Temp1,Par1,NH4R)
               DBio(i,k,iNH4) = DBio(i,k,iNH4) - Nitrif * dtdays
               DBio(i,k,iNO3) = DBio(i,k,iNO3) + Nitrif * dtdays

#if defined BIOFLUX && defined GAK1D
               IF (i.eq.3.and.j.eq.3) THEN
                 bflx(iDet,iNH4)= bflx(iDet,iNH4) + cff1
               END IF
#endif
            END DO
          END DO
!IminS:ImaxS
! ----------------------------------------------------------------------
! Sinking (code external: CVL)
!
          call BIOSINK(ng, wPhS, Bio(IminS,1,iPhS), DBio(IminS,1,iPhS), &
     &                  HzL, dtdays, z_wL, 1.0_r8, LBi, UBi,            &
     &                  IminS, ImaxS)
          call BIOSINK(ng, wPhL, Bio(IminS,1,iPhL), DBio(IminS,1,iPhL), &
     &                  HzL, dtdays, z_wL, 1.0_r8, LBi, UBi,            &
     &                  IminS, ImaxS)
             !------------------------------------------------------
             !ROMS seems to bomb with shallow water. Use of the
             !DetSINK option seems to have problems in shallow water
             !See if it will run OK for water deeper than 60
	     ! (this check now happens inside DetSINK)
             !------------------------------------------------------
          call DetSINK(ng, wDet, Bio(IminS,1,iDet), DBio(IminS,1,iDet), &
     &                 HzL, dtdays, z_wL, 1.0_r8, Dens, LBi, UBi,       &
     &                 IminS, ImaxS)
#ifdef DIAPAUSE
          IF ( downward ) THEN
             call BIOSINK(ng, wNCsink, Bio(IminS,1,iNCa),               &
     &                DBio(IminS,1,iNCa), HzL,dtdays, z_wL, -1*NCmaxz,  &
     &                LBi, UBi, IminS, ImaxS)
           ELSE IF (upward) THEN
             call BIORISE(ng, wNCrise, Bio(IminS,1,iNCa),               &
     &              DBio(IminS,1,iNCa), HzL,dtdays, z_wL, -1*NCmaxz,    &
     &              LBi, UBi, 6.0_r8, IminS, ImaxS)
           END IF
#endif
!
! ----------------------------------------------------------------------
! Update Bio array
! ----------------------------------------------------------------------
!
          DO itrc=1,NBT
            ibio=idbio(itrc)
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,ibio)=Bio(i,k,ibio)+DBio(i,k,ibio)
              END DO
            END DO
          END DO
        END DO ITER_LOOP
!-----------------------------------------------------------------------
!  Update global tracer variables (m Tunits).
!-----------------------------------------------------------------------
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              t(i,j,k,nnew,ibio)=MAX(t(i,j,k,nnew,ibio)+                &
     &                               (Bio(i,k,ibio)-Bio_bak(i,k,ibio))* &
     &                               Hz(i,j,k),                         &
     &                               0.0_r8)
            END DO
          END DO
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
            st(i,j,k,nnew,iPhSprd) = st(i,j,k,nnew,iPhSprd) +           &
     &                             Prod(i,k,iPhS)
            st(i,j,k,nnew,iPhLprd) = st(i,j,k,nnew,iPhLprd) +           &
     &                             Prod(i,k,iPhL)
            st(i,j,k,nnew,iMZSprd) = st(i,j,k,nnew,iMZSprd) +           &
     &                             Prod(i,k,iMZS)
            st(i,j,k,nnew,iMZLprd) = st(i,j,k,nnew,iMZLprd) +           &
     &                             Prod(i,k,iMZL)
            st(i,j,k,nnew,iCopPrd) = st(i,j,k,nnew,iCopPrd) +           &
     &                             Prod(i,k,iCop)
            st(i,j,k,nnew,iNCaPrd) = st(i,j,k,nnew,iNCaPrd) +           &
     &                             Prod(i,k,iNCa)
            st(i,j,k,nnew,iEupPrd) = st(i,j,k,nnew,iEupPrd) +           &
     &                             Prod(i,k,iEup)
          END DO
        END DO

      END DO J_LOOP

      RETURN
      END SUBROUTINE biology_tile
!
!================================================ Craig V. W. Lewis ==
!                                                                    !
!  This routine extracts the sinking code from biology.F and allows  !
!  it to be used as an external module for any number of fields      !
!  independently; it is a hack, but better than the current setup    !
!                                                                    !
!  Added 1 argument for diapause calculations:                       !
!    zlimit   - lowest depth for sinking.  -1 = sinks out at         !
!               constant rate                                        !
!  Liz Dobbins 4/15/03                                               !
!=====================================================================
!
      subroutine BIOSINK(ng,wBio,Bio,dBioOut,HzL,dtdays,z_wL,zlimit,    &
     &           LBi, UBi, IminS, ImaxS)
!
      USE mod_param
!
      implicit none
!
      integer, intent(in) :: ng, LBi, UBi, IminS, ImaxS
      real(r8), intent(in) :: wBio
      real(r8), intent(in) :: Bio(IminS:ImaxS,N(ng))
      real(r8), intent(inout) :: dBioOut(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: HzL(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: dtdays
      real(r8), intent(in) :: z_wL(IminS:ImaxS,0:N(ng))
      real(r8), intent(in) :: zlimit
!
      integer :: i, k
      real(r8) :: aL, aR, cff1, cff2, cu
      real(r8) :: FC(IminS:ImaxS,0:N(ng))
      real(r8) :: dBio(IminS:ImaxS,0:N(ng))
      real(r8) :: wBiod(IminS:ImaxS,0:N(ng))
!
!
      IF ( zlimit .lt. 0 ) THEN
        DO k=0,N(ng)
          DO i=LBi, UBi
            IF ( z_wL(i,k) .ge. zlimit ) THEN
              wBiod(i,k) = wBio*exp( -1*(z_wL(i,k)-(zlimit/2))**2 /     &
     &          (zlimit/2)**2 )
            ELSE
              wBiod(i,k) = 0.0_r8
            END IF
          END DO
        END DO
      ELSE
        DO k=0,N(ng)
          DO i=LBi, UBi
            wBiod(i,k) = wBio
          END DO
        END DO
      END IF
!
!  Vertical sinking: Vertical advection algorithm based on monotonic,
!  continuous conservative parabolic splines.
!  Construct parabolic splines:  Compute vertical derivatives for
!  "Bio" at W-points.  Neumann boundary conditions are
!  assumed on top and bottom.
!
      DO i=LBi,UBi
        FC(i,0)=0.0_r8
        dBio(i,0)=0.0_r8
      END DO
      do k=1,N(ng)-1
        DO i=LBi,UBi
         cff1=1.0_r8/(2.0_r8*HzL(i,k+1)+                                &
     &        HzL(i,k)*(2.0_r8-FC(i,k-1)))
         FC(i,k)=cff1*HzL(i,k+1)
         dBio(i,k)=cff1*(6.0_r8*(Bio(i,k+1)-Bio(i,k))-                  &
     &        HzL(i,k)*dBio(i,k-1))
        END DO
      END DO

      DO i=LBi,UBi
        dBio(i,N(ng))=0.0_r8
      END DO

      DO k=N(ng)-1,1,-1
        DO i=LBi,UBi
         dBio(i,k)=dBio(i,k)-FC(i,k)*dBio(i,k+1)
        END DO
      END DO
!
!  Convert vertical derivatives "dBio" into field values
!  at grid box interfaces assuming parabolic profiles within each
!  grid box.  Restrict these values to lie between bounds determined
!  from box-averaged values of grid boxes adjscent from above and
!  below. This restriction is part of the PPM-like monotonization
!  procedure.
!
      cff1=1.0_r8/3.0_r8
      DO i=LBi,UBi
        dBio(i,0)=Bio(i,1)          ! -cff*HzL(1)*(dBio(0)+0.5_r8*dBio(1))
        dBio(i,N(ng))=Bio(i,N(ng))  ! +cff*HzL(N(ng))*(dBio(N(ng))+0.5_r8*dBio(N(ng)-1))
      END DO
      DO k=2,N(ng)
        DO i=LBi,UBi
         dBio(i,k-1)=Bio(i,k)-cff1*HzL(i,k)                             &
     &        *(0.5_r8*dBio(i,k)+dBio(i,k-1))
         dBio(i,k-1)=MAX(dBio(i,k-1),MIN(Bio(i,k-1),Bio(i,k)))
         dBio(i,k-1)=MIN(dBio(i,k-1),MAX(Bio(i,k-1),Bio(i,k)))
        END DO
      END DO
!
!  Convert "dBio" into flux-integrated values;  complete
!  PPM flux limiting.  This procedure starts from assigning left and
!  right (aR,aL) values of the interpolating parabolae, then the
!  monotonicity conditions are checked and aL, aR are modified to fit.
!
      DO k=1,N(ng)
        DO i=LBi,UBi
         FC(i,k)=dtdays/HzL(i,k)
         aR=dBio(i,k)
         aL=dBio(i,k-1)
         cff1=(aR-aL)*6.0_r8*(Bio(i,k)-0.5_r8*(aR+aL))
         cff2=(aR-aL)**2
         IF ((aR-Bio(i,k))*(Bio(i,k)-aL).lt.0.0_r8) THEN
            aL=Bio(i,k)
            aR=Bio(i,k)
         ELSE IF (cff1.gt.cff2) THEN
            aL=3.0_r8*Bio(i,k)-2.0_r8*aR
         ELSE IF (cff1.lt.-cff2) THEN
            aR=3.0_r8*Bio(i,k)-2.0_r8*aL
         END IF
         cu=wBio*FC(i,k)
         dBio(i,k-1)=Bio(i,k)-(1.0_r8-cu)*(0.5_r8*(aR-aL)-              &
     &        (0.5_r8*(aR+aL)-Bio(i,k))*(1.0_r8-2.0_r8*cu))
        END DO
      END DO
      DO i=LBi,UBi
        dBio(i,N(ng))=0.0_r8
      END DO
!
! Set change in biological variable.
!
      do k=1,N(ng)
        DO i=LBi,UBi
         dBioOut(i,k) = dBioOut(i,k) + wBiod(i,k)*FC(i,k)*dBio(i,k) -   &
     &        wBiod(i,k-1)*FC(i,k)*dBio(i,k-1)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE BIOSINK
!=====================================================================
!                                                                    !
!   It's Craig's biosink.F turned upside-down, based on an idea of   !
!   Hal Batchelder''s which he used to move NO3 upwards to simulate   !
!   upwelling                                                        !
!                                                                    !
!  Added 1 argument for diapause calculations:                       !
!    zlimit   - lowest depth for sinking.  -1 = constant rate        !
!  Liz Dobbins 4/15/03                                               !
!=====================================================================
!
      subroutine BIORISE(ng, wBio, Bio, dBioOut, HzL, dtdays, z_wL,     &
     &                   zlimit, LBi, UBi, dlimit, IminS, ImaxS)
!
      USE mod_param
!
      implicit none
!
      integer, intent(in) :: ng, LBi,UBi, IminS, ImaxS
      real(r8), intent(in) :: wBio
      real(r8), intent(in) :: Bio(IminS:ImaxS,N(ng))
      real(r8), intent(inout) :: dBioOut(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: HzL(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: dtdays
      real(r8), intent(in) :: z_wL(IminS:ImaxS,0:N(ng))
      real(r8), intent(in) :: zlimit, dlimit
!
      integer :: i, k
      real(r8) :: aL, aR, cff1, cff2, cu
      real(r8) :: FC(IminS:ImaxS,1:N(ng)+1)
      real(r8) :: dBio(IminS:ImaxS,1:N(ng)+1)
      real(r8) :: wBiod(IminS:ImaxS,1:N(ng)+1)
!
      IF ( zlimit .lt. 0 ) THEN
        DO k=1,N(ng)
          DO i=LBi,UBi
            IF ( z_wL(i,k-1).ge.zlimit/2) THEN
              wBiod(i,k) = wBio*exp( -1*(z_wL(i,k-1)-(zlimit/2))**2 /   &
     &          (zlimit/2)**2 )
            ELSE
              wBiod(i,k) = wBio
            END IF
! This check used to be outside the whole function call.
            IF ( Bio(i,N(ng)).ge.dlimit) THEN
              wBiod(i,k) = 0.0_r8
            END IF
          END DO
        END DO
        DO i=LBi,UBi
          wBiod(i,N(ng)+1) = 0.0_r8
        END DO
      ELSE
        DO k=1,N(ng)+1
          DO i=LBi,UBi
            wBiod(i,k) = wBio
          END DO
        END DO
      END IF
!
!  Vertical rising: Vertical advection algorithm based on monotonic,
!  continuous conservative parabolic splines.
!  Construct parabolic splines:  Compute vertical derivatives for
!  "Bio" at W-points.  Neumann boundary conditions are
!  assumed on top and bottom.
!
      DO i=LBi,UBi
        FC(i,N(ng)+1)=0.0_r8
        dBio(i,N(ng)+1)=0.0_r8
      END DO
      DO k=N(ng),2,-1
        DO i=LBi,UBi
         cff1=1.0_r8/(2.0_r8*HzL(i,k-1)+                                  &
     &        HzL(i,k)*(2.0_r8-FC(i,k+1)))
         FC(i,k)=cff1*HzL(i,k-1)
         dBio(i,k)=cff1*(6.0_r8*(Bio(i,k-1)-Bio(i,k))-                        &
     &        HzL(i,k)*dBio(i,k+1))
        END DO
      END DO

      DO i=LBi,UBi
        dBio(i,1)=0.0_r8
      END DO

      DO k=2,N(ng)
        DO i=LBi,UBi
         dBio(i,k)=dBio(i,k)-FC(i,k)*dBio(i,k-1)
        END DO
      END DO
!
!  Convert vertical derivatives "dBio" into field values
!  at grid box interfaces assuming parabolic profiles within each
!  grid box.  Restrict these values to lie between bounds determined
!  from box-averaged values of grid boxes adjscent from above and
!  below. This restriction is part of the PPM-like monotonization
!  procedure.
!
      cff1=1.0_r8/3.0_r8
      DO i=LBi,UBi
        dBio(i,N(ng)+1)=Bio(i,N(ng))  ! -cff*HzL(1)*(dBio(0)+0.5_r8*dBio(1))
        dBio(i,1)=Bio(i,1)            ! +cff*HzL(N(ng))*(dBio(N(ng))+0.5_r8*dBio(N(ng)-1))
      END DO
      DO k=N(ng)-1,1,-1
        DO i=LBi,UBi
         dBio(i,k+1)=Bio(i,k)-cff1*HzL(i,k)                                   &
     &        *(0.5_r8*dBio(i,k)+dBio(i,k+1))
         dBio(i,k+1)=MAX(dBio(i,k+1),MIN(Bio(i,k+1),Bio(i,k)))
         dBio(i,k+1)=MIN(dBio(i,k+1),MAX(Bio(i,k+1),Bio(i,k)))
        END DO
      END DO
!
!  Convert "dBio" into flux-integrated values;  complete
!  PPM flux limiting.  This procedure starts from assigning left and
!  right (aR,aL) values of the interpolating parabolae, then the
!  monotonicity conditions are checked and aL, aR are modified to fit.
!
      DO k=N(ng),1,-1
        DO i=LBi,UBi
         FC(i,k)=dtdays/HzL(i,k)
         aR=dBio(i,k)
         aL=dBio(i,k+1)
         cff1=(aR-aL)*6.0_r8*(Bio(i,k)-0.5_r8*(aR+aL))
         cff2=(aR-aL)**2
         IF ((aR-Bio(i,k))*(Bio(i,k)-aL).lt.0.0_r8) THEN
            aL=Bio(i,k)
            aR=Bio(i,k)
         ELSE IF (cff1.gt.cff2) THEN
            aL=3.0_r8*Bio(i,k)-2.0_r8*aR
         ELSE IF (cff1.lt.-cff2) THEN
            aR=3.0_r8*Bio(i,k)-2.0_r8*aL
         END IF
         cu=wBio*FC(i,k)
         dBio(i,k+1)=Bio(i,k)-(1.0_r8-cu)*(0.5_r8*(aR-aL)-                  &
     &        (0.5_r8*(aR+aL)-Bio(i,k))*(1.0_r8-2.0_r8*cu))
        END DO
      END DO
      DO i=LBi,UBi
        dBio(i,1)=0.0_r8
      END DO
!
! Set change in biological variable.
!
      DO k=1,N(ng)
        DO i=LBi,UBi
         dBioOut(i,k) = dBioOut(i,k) + wBiod(i,k)*FC(i,k)*dBio(i,k) -   &
     &                  wBiod(i,k+1)*FC(i,k)*dBio(i,k+1)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE BIORISE
!====================================================================
        Function ComputeDensity(Temp1,Sal1)
!----------------------------------------------------------------
! Computes the water column density from salinity and temperature
! Returns sigma-t
!----------------------------------------------------------------
        USE mod_kinds

        Real(r8) ComputeDensity
        Real(r8) Temp1, Sal1
        Real(r8) Sig
        Sig = 999.842594 + 0.06793952 * Temp1
        Sig = Sig - 0.00909529 * Temp1 ** 2 +                      &
     &          0.0001001685 * Temp1 ** 3
        Sig = Sig - 0.000001120083 * Temp1 ** 4 +                  &
     &          0.000000006536332 * Temp1 ** 5
        Sig = Sig + 0.824493 * Sal1 - 0.0040899 * Temp1 * Sal1
        Sig = Sig + 0.000076438 * Temp1 ** 2 * Sal1 -              &
     &          0.00000082467 * Temp1 ** 3 * Sal1
        Sig = Sig + 0.0000000053875 * Temp1 ** 4 * Sal1 -          &
     &          0.00572466 * Sal1 ** (3 / 2)
        Sig = Sig + 0.00010227 * Temp1 * Sal1 ** (3 / 2) -         &
     &          0.0000016546 * Temp1 ** 2 * Sal1 ** (3 / 2)
        Sig = Sig + 0.00048314 * Sal1 ** 2
         ComputeDensity = Sig - 1000
        End Function ComputeDensity
!===============================================================
        Function GetPhytoResp2(Temp1, Tref, KbmPh)
!------------------------------------------------------
! Computes the temperature correction for phytoplankton
! respiration according to Arhonditsis 2005.
!------------------------------------------------------
        USE mod_kinds

        Real(r8) GetPhytoResp2
        Real(r8) Temp1      !Temperature, passed
        Real(r8) Tref       !Reference temperature
        Real(r8) KbmPh      !Half saturation, temperature

        Real(r8) Resp       !Returned variable

        Resp = exp(KbmPh * (Temp1 - Tref))
        GetPhytoResp2 = Resp
        Return
        End Function GetPhytoResp2
!=====================================================================
      FUNCTION GetLightLimIronSml(alphaPh, PAR1, Pmax1,          &
     &     CrChlRatio1,IronLim1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLimIronSml
        Real(r8) LightLim,OffSet
        OffSet = 0.0
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)   &
     &             / Pmax1 / CrChlRatio1 / IronLim1)
       GetLightLimIronSml = LightLim
      END FUNCTION GetLightLimIronSml
!=====================================================================
      FUNCTION GetLightLimSml(alphaPh, PAR1, Pmax1, CrChlRatio1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLimSml
        Real(r8) LightLim,OffSet
        OffSet = 0.0
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)   &
     &             / Pmax1 / CrChlRatio1 )
       GetLightLimSml = LightLim
      END FUNCTION GetLightLimSml
!=====================================================================
      FUNCTION GetLightLimIron(alphaPh, PAR1, Pmax1, CrChlRatio1,  &
     &     IronLim1, ParMax)
!------------------------------------------------------------------
! Light lim with varying alpha. Works with iron limitation. Alph is
! a function of the surface light intensity.
!------------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

       Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
       Real(r8) GetLightLimIron
       Real(r8) Alpha,LightLim,OffSet,ParMax

       !Alpha = 1.e-8*EXP(0.48*ParMax) + 0.5
       !if (Alpha .gt. 10) Alpha = 10
       !--------------------------
       !Use a simple step function
       !--------------------------
       if (ParMax.lt.48_r8) then
          Alpha = 1._r8
       else
          Alpha = 4._r8
       end if
       LightLim = TANH( Alpha * Par1/Pmax1/CrChlRatio1/IronLim1)
       GetLightLimIron = LightLim

      END FUNCTION GetLightLimIron
!=======================================================================
      FUNCTION GetLightLim(alphaPh, PAR1, Pmax1, CrChlRatio1, ParMax)
!-----------------------------------------------------------------
! Generates a light lim with varying alphaPh without iron
!-----------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

       Real(r8) :: alphaPh,Pmax1,PAR1,CrChlRatio1
       Real(r8) GetLightLim
       Real(r8) Alpha,LightLim,OffSet,ParMax

!       Alpha = 1.e-8*EXP(0.48*ParMax) + 0.5
!       if (Alpha .gt. 10) Alpha = 10
       !--------------------------
       !Use a simple step function
       !--------------------------
       if (ParMax.lt.48_r8) then
          Alpha = 1._r8
       else
          Alpha = 4._r8
       end if
       LightLim = TANH(alphaPh * PAR1/Pmax1/CrChlRatio1)

       GetLightLim = LightLim
      END FUNCTION GetLightLim
!==============================================================
        Function GetCopepodResp(Temp1,respVal)
!--------------------------------------------------------------
! Computes copepod respiration according to Arhonditsis (2005).
!--------------------------------------------------------------
        USE mod_kinds
        USE mod_param

        Real(r8) GetCopepodResp
        real(r8) :: respVal
        Real(r8) Temp1       !Passed variable
        Real(r8) :: bm  = 0.04   !Basal metabolic rate day**-1
        Real(r8) :: ktbm  = 0.05 !Temperature response degrees C**-1
        Real(r8) :: Tref  = 20   !Reference temperature degrees C
        Real(r8) Resp        !Returned variable

        Resp = respVal * exp(ktbm * (Temp1 - Tref))
        GetCopepodResp = Resp
        Return
        End Function GetCopepodResp
!=================================================================
        Function GetBasalMetabolism(respPh,kfePh,Iron1)
!---------------------------------------------------------
! Computes an iron correction for the basal metabolism for
! the phytoplankton respiration calculation
!---------------------------------------------------------
        USE mod_kinds

        Real(r8) GetBasalMetabolism
        Real(r8) Iron1     !Iron concentration
        Real(r8) kfePh     !Half saturation for iron
        Real(r8) respPh    !Phytoplankton uncorrected basal metabolism
        Real(r8) BaseMet   !Phytoplankton basal metabolism

        BaseMet = Iron1/(kfePh + Iron1)
        BaseMet = BaseMet * ((kfePh +2)/2)
        GetBasalMetabolism = BaseMet * respPh

        Return
        End Function GetBasalMetabolism
!=========================================================
        Function GetNitrif(Temp1,Dep1,NH4R)
!-------------------------------------------------------
! Computes the nitrification with respect to temperature
! according to Arhonditsis (2005).  Generates depth
! correction according to Denman (2003).
!-------------------------------------------------------
        USE mod_kinds
        Real(r8) GetNitrif
                               !--------------------------------
        Real(r8) Temp1, Dep1, NH4R !Passed variables
        Real(r8) NH4conv  /14/     !mg N/mol N
        Real(r8) KNH4Nit  /0.08/   !Half Sat Con mg N/m3/day
        Real(r8) KTNitr   /0.002/  !Temperature responce dec C^2
        Real(r8) ToptNtr  /28/     !Optimum nitrification temp
        Real(r8) Zox      /20/     !50% nitrification depth
        Real(r8) Nexp     /6/      !Exponent to adjust profile shape
        Real(r8) NitrMax  /0.011/  !Maximum nitrification (mM/m3/d
                               !--------------------------------

        Real(r8) Nitr, DepCor

        NH4R = NH4R * NH4conv
        Nitr = NH4R/(KNH4Nit + NH4R)
        Nitr = Nitr * exp(-KTNitr*(Temp1 - ToptNtr)**2)
        DepCor = (Dep1**Nexp)/( (Zox**Nexp) + Dep1**Nexp)
        Nitr = (Nitr * DepCor) * NitrMax
        GetNitrif = Nitr
        Return
        End Function GetNitrif
!========================================================================
        Function GetNitrif2(Temp1,PAR1,NH4R)
        !---------------------------------------------------------
        !Computes nitrificaton from Kawamiya with light correction
        !from Fennel; Kawamiya (2000), Fennel (2006)
        !---------------------------------------------------------
        USE mod_kinds
        USE mod_param
         Real(r8) GetNitrif2
                               !--------------------------------
        Real(r8) :: Temp1, NH4R       !Passed variables
        Real(r8) :: I0 = 0.0095     !Threshold,light inhibition, W m-2
        Real(r8) :: KI = 4.0        !Half Saturation light intensity, W m-2
        Real(r8) :: KN0 = 0.03       !Nitrification at 0 deg C, day-1
        Real(r8) :: KNT = 0.0693     !Temperature coefficient
        Real(r8) :: ParW              !Par in watts
        Real(r8) :: NitrMax           !Maximum nitrification
        Real(r8) :: Nitr              !Nitrification
                               !---------------------------------
        Real(r8) :: cff1, PAR1

        !-----------------------------------
        !Temperature dependent nitrification
        !-----------------------------------
        KI = 1.5
        KN0 = 0.15
        !KNT = 0.07
        NitrMax = (KN0*Exp(KNT*Temp1))*NH4R
        !-----------------------------------
        !Convert PAR in E m-2 d-1 to W day-1
        !-----------------------------------
        ParW = PAR1/0.394848_r8
        !---------------------------------
        !Light correction of nitrification
        !---------------------------------
        cff1 = (ParW-I0)/(KI+ParW-I0)
        Nitr = NitrMax*(1-MAX(0.0_r8,cff1))
        GetNitrif2 = Nitr
        Return
        End Function GetNitrif2
!=====================================================================
      FUNCTION GetLightLimIron2(alphaPh, PAR1, Pmax1, CrChlRatio1, &
     &     IronLim1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLimIron2
        Real(r8) LightLim,OffSet
        OffSet = 0.0_r8
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)   &
     &             / Pmax1 / CrChlRatio1 / IronLim1)
       GetLightLimIron2 = LightLim
      END FUNCTION GetLightLimIron2
!=====================================================================
      FUNCTION GetLightLim2(alphaPh, PAR1, Pmax1, CrChlRatio1)
!---------------------------------------------------------------
! Uses a normal hyperbolic tangent function for light limitation
!---------------------------------------------------------------
      USE mod_kinds
      USE mod_param
!
      implicit none

        Real(r8) :: alphaPh,Pmax1,IronLim1,PAR1,CrChlRatio1
        Real(r8) GetLightLim2
        Real(r8) LightLim,OffSet
        OffSet = 0.0_r8
          LightLim = TANH( alphaPh * MAX((PAR1 - OffSet),0.0_r8)   &
     &             / Pmax1 / CrChlRatio1 )
       GetLightLim2 = LightLim
      END FUNCTION GetLightLim2
!=======================================================================
      subroutine DetSINK(ng, wBio, Bio, dBioOut, HzL, dtdays, z_wL,     &
     &             zlimit, Dens, LBi, UBi, IminS, ImaxS)
!================================================ Craig V. W. Lewis ==
!                                                                    !
!  This routine extracts the sinking code from biology.F and allows  !
!  it to be used as an external module for any number of fields      !
!  independently; it is a hack, but better than the current setup    !
!                                                                    !
!  Added 1 argument for diapause calculations:                       !
!    zlimit   - lowest depth for sinking.  -1 = sinks out at         !
!               constant rate                                        !
!  Liz Dobbins 4/15/03                                               !
!  Modified to use density as proxy for sinking rate                 !
!=====================================================================
!
      USE mod_kinds
      USE mod_param
!
      implicit none
!
      integer, intent(in) :: ng, LBi, UBi, IminS, ImaxS
      real(r8), intent(in) :: wBio
      real(r8), intent(inout) :: Bio(IminS:ImaxS,N(ng))
      real(r8), intent(inout) :: dBioOut(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: Dens(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: HzL(IminS:ImaxS,N(ng))
      real(r8), intent(in) :: dtdays
      real(r8), intent(in) :: z_wL(IminS:ImaxS,0:N(ng))
      real(r8), intent(in) :: zlimit
!
      integer :: i, k
      real(r8) :: aL, aR, cff1, cff2, cu
      real(r8) :: FC(IminS:ImaxS,0:N(ng))
      real(r8) :: dBio(IminS:ImaxS,0:N(ng))
      real(r8) :: wBiod(IminS:ImaxS,0:N(ng))
!---------------------------------------------------------
!     NOTE: wBio is the sinking rate passed from goanpz.in
!---------------------------------------------------------
      DO i=LBi,UBi
        IF (z_wL(i,N(ng)) .gt. -60._r8) THEN
        ! Shallow water case
          DO k=0,N(ng)
            IF ( zlimit .lt. 0 ) THEN
              IF ( z_wL(i,k) .ge. zlimit ) THEN
                wBiod(i,k) = wBio*exp( -1*(z_wL(i,k)-(zlimit/2))**2 /          &
     &          (zlimit/2)**2 )
              ELSE
                wBiod(i,k) = 0.0_r8
              END IF
            ELSE
              wBiod(i,k) = wBio
            END IF
	  END DO
	ELSE
          DO k=1,N(ng)
            wBiod(i,k) = wBio*(26.25 - Dens(i,k))/(0.03 * Dens(i,k))
            if (wBiod(i,k) .lt. 0.0_r8) wBiod(i,k) = 0.0_r8
          END DO
          wBiod(i,0) = wBio*(26.25 - Dens(i,1))/(0.03 * Dens(i,1))
          if (wBiod(i,0) .lt. 0.0_r8) wBiod(i,0) = 0.0_r8
        END IF
      END DO
!
!  Vertical sinking: Vertical advection algorithm based on monotonic,
!  continuous conservative parabolic splines.
!  Construct parabolic splines:  Compute vertical derivatives for
!  "Bio" at W-points.  Neumann boundary conditions are
!  assumed on top and bottom.
!
      DO i=LBi,UBi
        FC(i,0)=0.0_r8
        dBio(i,0)=0.0_r8
      END DO

      DO k=1,N(ng)-1
        DO i=LBi,UBi
         cff1=1.0_r8/(2.0_r8*HzL(i,k+1)+                                  &
     &        HzL(i,k)*(2.0_r8-FC(i,k-1)))
         FC(i,k)=cff1*HzL(i,k+1)
         dBio(i,k)=cff1*(6.0_r8*(Bio(i,k+1)-Bio(i,k))-                        &
     &        HzL(i,k)*dBio(i,k-1))
        END DO
      END DO

      DO i=LBi,UBi
        dBio(i,N(ng))=0.0_r8
      END DO

      DO k=N(ng)-1,1,-1
        DO i=LBi,UBi
         dBio(i,k)=dBio(i,k)-FC(i,k)*dBio(i,k+1)
        END DO
      END DO
!
!  Convert vertical derivatives "dBio" into field values
!  at grid box interfaces assuming parabolic profiles within each
!  grid box.  Restrict these values to lie between bounds determined
!  from box-averaged values of grid boxes adjscent from above and
!  below. This restriction is part of the PPM-like monotonization
!  procedure.
!
      cff1=1.0_r8/3.0_r8
      DO i=LBi,UBi
        dBio(i,0)=Bio(i,1)            ! -cff*HzL(1)*(dBio(0)+0.5_r8*dBio(1))
        dBio(i,N(ng))=Bio(i,N(ng))    ! +cff*HzL(N(ng))*(dBio(N(ng))+0.5_r8*dBio(N(ng)-1))
      END DO
      DO k=2,N(ng)
        DO i=LBi,UBi
         dBio(i,k-1)=Bio(i,k)-cff1*HzL(i,k)                                   &
     &        *(0.5_r8*dBio(i,k)+dBio(i,k-1))
         dBio(i,k-1)=MAX(dBio(i,k-1),MIN(Bio(i,k-1),Bio(i,k)))
         dBio(i,k-1)=MIN(dBio(i,k-1),MAX(Bio(i,k-1),Bio(i,k)))
        END DO
      END DO
!
!  Convert "dBio" into flux-integrated values;  complete
!  PPM flux limiting.  This procedure starts from assigning left and
!  right (aR,aL) values of the interpolating parabolae, then the
!  monotonicity conditions are checked and aL, aR are modified to fit.
!
      DO k=1,N(ng)
        DO i=LBi,UBi
         FC(i,k)=dtdays/HzL(i,k)
         aR=dBio(i,k)
         aL=dBio(i,k-1)
         cff1=(aR-aL)*6.0_r8*(Bio(i,k)-0.5_r8*(aR+aL))
         cff2=(aR-aL)**2
         IF ((aR-Bio(i,k))*(Bio(i,k)-aL).lt.0.0_r8) THEN
            aL=Bio(i,k)
            aR=Bio(i,k)
         ELSE IF (cff1.gt.cff2) THEN
            aL=3.0_r8*Bio(i,k)-2.0_r8*aR
         ELSE IF (cff1.lt.-cff2) THEN
            aR=3.0_r8*Bio(i,k)-2.0_r8*aL
         END IF
         cu=wBio*FC(i,k)
         dBio(i,k-1)=Bio(i,k)-(1.0_r8-cu)*(0.5_r8*(aR-aL)-                  &
     &        (0.5_r8*(aR+aL)-Bio(i,k))*(1.0_r8-2.0_r8*cu))
        END DO
      END DO
      DO i=LBi,UBi
        dBio(i,N(ng))=0.0_r8
      END DO
!
! Set change in biological variable.
!
      DO k=1,N(ng)
        DO i=LBi,UBi
         dBioOut(i,k) = dBioOut(i,k) + wBiod(i,k)*FC(i,k)*dBio(i,k) -   &
     &                 wBiod(i,k-1)*FC(i,k)*dBio(i,k-1)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE DetSINK
