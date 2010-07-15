#include "cppdefs.h"
      MODULE set_avg_mod
#if defined AVERAGES && (!defined ADJOINT && defined NONLINEAR)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  fields.  Due to synchronization, the time-averaged fields are       !
!  computed in delayed mode. All averages are accumulated at the       !
!  beggining of the next time-step.                                    !
!                                                                      !
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!  It computes leas-squares coefficients to detide time-averaged       !
!  fields. Notice that "set_detide" is called last since we need       !
!  the regular time-averages for those fields to detide.               !
!                                                                      !
# endif
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: set_avg

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_avg (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_stepping
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      USE mod_tides
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 5)
# endif
      CALL set_avg_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
# ifdef SOLVE3D
#  ifdef ICE_MODEL
     &                   IOUT, IUOUT, IEOUT,                            &
#  endif
     &                   NOUT,                                          &
# endif
     &                   KOUT)

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      CALL set_detide_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      NTC(ng), KOUT,                              &
#  ifdef SOLVE3D
     &                      NOUT,                                       &
#  endif
     &                      TIDES(ng) % CosOmega,                       &
     &                      TIDES(ng) % SinOmega,                       &
     &                      TIDES(ng) % CosW_avg,                       &
     &                      TIDES(ng) % CosW_sum,                       &
     &                      TIDES(ng) % SinW_avg,                       &
     &                      TIDES(ng) % SinW_sum,                       &
     &                      TIDES(ng) % CosWCosW,                       &
     &                      TIDES(ng) % SinWSinW,                       &
     &                      TIDES(ng) % SinWCosW)
# endif
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 5)
# endif
      RETURN
      END SUBROUTINE set_avg
!
!***********************************************************************
      SUBROUTINE set_avg_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
# ifdef SOLVE3D
#  ifdef ICE_MODEL
     &                         Iout, Iuout, Ieout,                      &
#  endif
     &                         Nout,                                    &
# endif
     &                         Kout)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_average
# if defined FORWARD_WRITE && defined SOLVE3D
      USE mod_coupling
# endif
      USE mod_forces
      USE mod_grid
# ifdef SOLVE3D
      USE mod_mixing
#  ifdef ICE_MODEL
      USE mod_ice
#  endif
# endif
      USE mod_ocean
      USE mod_scalars
# if defined BEST_NPZ || defined BIO_GOANPZ
      USE mod_biology
# endif
!
      USE vorticity_mod, ONLY : vorticity_tile
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kout
# ifdef SOLVE3D
#  ifdef ICE_MODEL
      integer, intent(in) :: Iout, Iuout, Ieout
#  endif
      integer, intent(in) :: Nout
# endif
!
!  Local variable declarations.
!
      integer :: i, it, itrc, j, k
#ifdef  BIOFLUX
       integer ::itrc2
#endif

      real(r8) :: fac
# ifdef SOLVE3D
      real(r8) :: potvor(LBi:UBi,LBj:UBj,N(ng))
      real(r8) :: relvor(LBi:UBi,LBj:UBj,N(ng))
# endif
      real(r8) :: potvor_bar(LBi:UBi,LBj:UBj)
      real(r8) :: relvor_bar(LBi:UBi,LBj:UBj)

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Compute vorticity fields.
!-----------------------------------------------------------------------
!
      IF (Aout(id2dPV,ng).or.Aout(id2dRV,ng).or.                        &
     &    Aout(id3dPV,ng).or.Aout(id3dRV,ng)) THEN
        CALL vorticity_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
# ifdef SOLVE3D
     &                       Kout, Nout,                                &
# else
     &                       Kout,                                      &
# endif
# ifdef MASKING
     &                       GRID(ng) % pmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
# endif
     &                       GRID(ng) % fomn,                           &
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % om_u,                           &
     &                       GRID(ng) % on_v,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
# ifdef SOLVE3D
     &                       GRID(ng) % z_r,                            &
     &                       OCEAN(ng) % pden,                          &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
# endif
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
     &                       OCEAN(ng) % zeta,                          &
# ifdef SOLVE3D
     &                       potvor, relvor,                            &
# endif
     &                       potvor_bar, relvor_bar)
      END IF
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.1)).or.                          &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
!
!  Initialize state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              aVERAGE(ng)%avgu2d(i,j)=OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=OCEAN(ng)%W(i,j,k)*           &
     &                                    GRID(ng)%pm(i,j)*             &
     &                                    GRID(ng)%pn(i,j)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=OCEAN(ng)%wvel(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=OCEAN(ng)%rho(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=OCEAN(ng)%t(i,j,k,Nout,it)
                END DO
              END DO
            END DO
          END IF
        END DO

#  ifdef NEMURO_SED1
        IF (Aout(idPONsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONsed(i,j)=OCEAN(ng)%PONsed(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idOPALsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALsed(i,j)=OCEAN(ng)%OPALsed(i,j)
            END DO
          END DO
        END IF
#  endif

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=MIXING(ng)%Akv(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=MIXING(ng)%Akt(i,j,k,itemp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=MIXING(ng)%Akt(i,j,k,isalt)
              END DO
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_SKPP
        IF (Aout(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsbl(i,j)=MIXING(ng)%hsbl(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (Aout(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghbbl(i,j)=MIXING(ng)%hbbl(i,j)
            END DO
          END DO
        END IF
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Initialize 2D/3D coupling terms.
!
        IF (Aout(idUfx1,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg1(i,j)=COUPLING(ng)%DU_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUfx2,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg2(i,j)=COUPLING(ng)%DU_avg2(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idVfx1,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg1(i,j)=COUPLING(ng)%DV_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVfx2,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg2(i,j)=COUPLING(ng)%DV_avg2(i,j)
            END DO
          END DO
        END IF
# endif
!
!  Initialize surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=FORCES(ng)%sustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=FORCES(ng)%svstr(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=FORCES(ng)%bustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=FORCES(ng)%bvstr(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=FORCES(ng)%stflx(i,j,itemp)
            END DO
          END DO
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=FORCES(ng)%stflx(i,j,isalt)
            END DO
          END DO
        END IF
#  ifdef SSSFLX
        IF (Aout(idSSSf,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsssflx(i,j)=FORCES(ng)%sssflx(i,j,isalt)
            END DO
          END DO
        END IF
#  endif
#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=FORCES(ng)%srflx(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=FORCES(ng)%lhflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=FORCES(ng)%lrflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=FORCES(ng)%shflx(i,j)
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=FORCES(ng)%evap(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgrain(i,j)=FORCES(ng)%rain(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif
# endif

# ifdef ICE_MODEL
        IF (Aout(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avguice(i,j)=ICE(ng)%ui(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout(idUice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgvice(i,j)=ICE(ng)%vi(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgaice(i,j)=ICE(ng)%ai(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghice(i,j)=ICE(ng)%hi(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsno(i,j)=ICE(ng)%hsn(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgtice(i,j)=ICE(ng)%tis(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgtimid(i,j)=ICE(ng)%ti(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idSfwat,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsfwat(i,j)=ICE(ng)%sfwat(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgageice(i,j)=ICE(ng)%ageice(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgiomflx(i,j)=ICE(ng)%io_mflux(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig11(i,j)=ICE(ng)%sig11(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig12(i,j)=ICE(ng)%sig12(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsig22(i,j)=ICE(ng)%sig22(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgT0mk(i,j)=ICE(ng)%t0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgS0mk(i,j)=ICE(ng)%s0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWfr(i,j)=ICE(ng)%wfr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWai(i,j)=ICE(ng)%wai(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWao(i,j)=ICE(ng)%wao(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWio(i,j)=ICE(ng)%wio(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWro(i,j)=ICE(ng)%wro(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgutau_iw(i,j)=ICE(ng)%utau_iw(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgchu_iw(i,j)=ICE(ng)%chu_iw(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef ICE_BIO
#  ifdef BERING_10K
        IF (Aout(idIceNO3,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgINO3(i,j)=ICE(ng)%IceNO3(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIcePhL,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgIPhL(i,j)=ICE(ng)%IcePhL(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIceNH4,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgINH4(i,j)=ICE(ng)%IceNH4(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIceLog,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgILog(i,j)=ICE(ng)%IceLog(i,j,Iout)
            END DO
          END DO
        END IF
#  endif
# endif

# ifdef BEST_NPZ
#  ifdef STATIONARY
        DO itrc=1,NTS(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgst(i,j,k,itrc)=                          &
     &                       OCEAN(ng)%st(i,j,k,Nout,itrc)
              END DO
            END DO
          END DO
        END DO
#  endif
#  ifdef PROD3
        DO itrc=1,NPT3(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgpt3(i,j,k,itrc)=                         &
     &                      OCEAN(ng)%pt3(i,j,k,Nout,itrc)
              END DO
            END DO
          END DO
        END DO
#  endif
#  ifdef PROD2
        DO itrc=1,NPT2(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgpt2(i,j,itrc)=                           &
     &                      OCEAN(ng)%pt2(i,j,Nout,itrc)
              END DO
            END DO
        END DO
#  endif
#  ifdef STATIONARY2
        DO itrc=1,NTS2(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgst2(i,j,itrc)=                           &
     &                      OCEAN(ng)%st2(i,j,Nout,itrc)
              END DO
            END DO
        END DO
#  endif
#  ifdef BENTHIC
           DO itrc=1,NBeT(ng)
             DO k=1,NBL(ng)
               DO j=JstrR,JendR
                 DO i=IstrR,IendR
                   AVERAGE(ng)%avgbt(i,j,k,itrc)=                       &
     &                      OCEAN(ng)%bt(i,j,k,Nout,itrc)
                 END DO
               END DO
             END DO
          END DO
#  endif
#  if defined CLIM_ICE_1D
           DO itrc=1,NIceT(ng)
             DO j=JstrR,JendR
               DO i=IstrR,IendR
                 AVERAGE(ng)%avgit(i,j,itrc)=                           &
     &                      OCEAN(ng)%it(i,j,Nout,itrc)
               END DO
            END DO
          END DO
#  endif
#  ifdef BIOFLUX
        DO itrc=1,NT(ng)
        DO itrc2=1,NT(ng)
                AVERAGE(ng)%avgbflx(itrc,itrc2)=                        &
     &                      OCEAN(ng)%bflx(itrc,itrc2)
          END DO
        END DO
#  endif
#endif

# ifdef BIO_GOANPZ
        DO itrc=1,NTS(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgst(i,j,k,itrc)=                          &
     &                      OCEAN(ng)%st(i,j,k,Nout,itrc)
              END DO
            END DO
          END DO
        END DO
# endif

# ifdef NEARSHORE_MELLOR
!
!  Initialize radiation stresses.
!
        IF (Aout(idU2Sd,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2Sd(i,j)=OCEAN(ng)%ubar_stokes(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2Sd,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2Sd(i,j)=OCEAN(ng)%vbar_stokes(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idU2rs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2RS(i,j)=MIXING(ng)%rustr2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2rs,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2RS(i,j)=MIXING(ng)%rvstr2d(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idW2xx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxx2d(i,j)=MIXING(ng)%Sxx_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idW2xy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxy2d(i,j)=MIXING(ng)%Sxy_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idW2yy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSyy2d(i,j)=MIXING(ng)%Syy_bar(i,j)
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idU3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3Sd(i,j,k)=OCEAN(ng)%u_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3Sd(i,j,k)=OCEAN(ng)%v_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idU3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3RS(i,j,k)=MIXING(ng)%rustr3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3RS(i,j,k)=MIXING(ng)%rvstr3d(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idW3xx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxx3d(i,j,k)=MIXING(ng)%Sxx(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3xy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxy3d(i,j,k)=MIXING(ng)%Sxy(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3yy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSyy3d(i,j,k)=MIXING(ng)%Syy(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3zx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzx3d(i,j,k)=MIXING(ng)%Szx(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3zy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzy3d(i,j,k)=MIXING(ng)%Szy(i,j,k)
              END DO
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Initialize vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=potvor_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=relvor_bar(i,j)
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=potvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=relvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
# endif
!
!  Initialize quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=OCEAN(ng)%zeta(i,j,Kout)*          &
     &                               OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=OCEAN(ng)%ubar(i,j,Kout)*          &
     &                               OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=OCEAN(ng)%vbar(i,j,Kout)*          &
     &                               OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=OCEAN(ng)%u(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=OCEAN(ng)%v(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=0.25_r8*                       &
     &                                   (OCEAN(ng)%u(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%u(i+1,j  ,k,Nout))* &
     &                                   (OCEAN(ng)%v(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%v(i  ,j+1,k,Nout))
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=GRID(ng)%Huon(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=GRID(ng)%Hvom(i,j,k)
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)*     &
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=0.5_r8*                   &
     &                                        OCEAN(ng)%u(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i-1,j,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i  ,j,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=0.5_r8*                   &
     &                                        OCEAN(ng)%v(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i,j-1,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i,j  ,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF

          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=0.5_r8*                &
     &                                           GRID(ng)%Huon(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i-1,j,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i  ,j,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=0.5_r8*                &
     &                                           GRID(ng)%Hvom(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i,j-1,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i,j  ,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG(ng)) THEN
!
!  Accumulate state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=AVERAGE(ng)%avgzeta(i,j)+        &
     &                                 OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=AVERAGE(ng)%avgu2d(i,j)+          &
     &                                OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=AVERAGE(ng)%avgv2d(i,j)+          &
     &                                OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=AVERAGE(ng)%avgu3d(i,j,k)+    &
     &                                    OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=AVERAGE(ng)%avgv3d(i,j,k)+    &
     &                                    OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=AVERAGE(ng)%avgw3d(i,j,k)+    &
     &                                    OCEAN(ng)%W(i,j,k)*           &
     &                                    GRID(ng)%pm(i,j)*             &
     &                                    GRID(ng)%pn(i,j)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=AVERAGE(ng)%avgwvel(i,j,k)+  &
     &                                     OCEAN(ng)%wvel(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=AVERAGE(ng)%avgrho(i,j,k)+    &
     &                                    OCEAN(ng)%rho(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=AVERAGE(ng)%avgt(i,j,k,it)+&
     &                                       OCEAN(ng)%t(i,j,k,Nout,it)
                END DO
              END DO
            END DO
          END IF
        END DO

#  ifdef NEMURO_SED1
        IF (Aout(idPONsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONsed(i,j)=AVERAGE(ng)%avgPONsed(i,j)+    &
     &                                 OCEAN(ng)%PONsed(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idOPALsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALsed(i,j)=AVERAGE(ng)%avgOPALsed(i,j)+  &
     &                                 OCEAN(ng)%OPALsed(i,j)
            END DO
          END DO
        END IF
#  endif

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=AVERAGE(ng)%avgAKv(i,j,k)+    &
     &                                    MIXING(ng)%Akv(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=AVERAGE(ng)%avgAKt(i,j,k)+    &
     &                                    MIXING(ng)%Akt(i,j,k,itemp)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=AVERAGE(ng)%avgAKs(i,j,k)+    &
     &                                    MIXING(ng)%Akt(i,j,k,isalt)
              END DO
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_SKPP
        IF (Aout(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsbl(i,j)=AVERAGE(ng)%avghsbl(i,j)+        &
     &                                 MIXING(ng)%hsbl(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (Aout(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghbbl(i,j)=AVERAGE(ng)%avghbbl(i,j)+        &
     &                                 MIXING(ng)%hbbl(i,j)
            END DO
          END DO
        END IF
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Accumulate 2D/3D coupling terms.
!
        IF (Aout(idUfx1,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg1(i,j)=AVERAGE(ng)%avgDU_avg1(i,j)+  &
     &                                    COUPLING(ng)%DU_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUfx2,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg2(i,j)=AVERAGE(ng)%avgDU_avg2(i,j)+  &
     &                                    COUPLING(ng)%DU_avg2(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idVfx1,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg1(i,j)=AVERAGE(ng)%avgDV_avg1(i,j)+  &
     &                                    COUPLING(ng)%DV_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVfx2,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg2(i,j)=AVERAGE(ng)%avgDV_avg2(i,j)+  &
     &                                    COUPLING(ng)%DV_avg2(i,j)
            END DO
          END DO
        END IF
# endif
!
!  Accumulate surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=AVERAGE(ng)%avgsus(i,j)+          &
     &                                FORCES(ng)%sustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=AVERAGE(ng)%avgsvs(i,j)+          &
     &                                FORCES(ng)%svstr(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=AVERAGE(ng)%avgbus(i,j)+          &
     &                                FORCES(ng)%bustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=AVERAGE(ng)%avgbvs(i,j)+          &
     &                                FORCES(ng)%bvstr(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=AVERAGE(ng)%avgstf(i,j)+          &
     &                                FORCES(ng)%stflx(i,j,itemp)
            END DO
          END DO
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=AVERAGE(ng)%avgswf(i,j)+          &
     &                                FORCES(ng)%stflx(i,j,isalt)
            END DO
          END DO
        END IF
#  ifdef SSSFLX
        IF (Aout(idSSSf,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsssflx(i,j)=AVERAGE(ng)%avgsssflx(i,j)+    &
     &                                FORCES(ng)%sssflx(i,j,isalt)
            END DO
          END DO
        END IF
#  endif
#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=AVERAGE(ng)%avgsrf(i,j)+          &
     &                                FORCES(ng)%srflx(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=AVERAGE(ng)%avglhf(i,j)+          &
     &                                FORCES(ng)%lhflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=AVERAGE(ng)%avgshf(i,j)+          &
     &                                FORCES(ng)%shflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=AVERAGE(ng)%avglrf(i,j)+          &
     &                                FORCES(ng)%lrflx(i,j)
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=AVERAGE(ng)%avgevap(i,j)+        &
     &                                 FORCES(ng)%evap(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgrain(i,j)=AVERAGE(ng)%avgrain(i,j)+        &
     &                                 FORCES(ng)%rain(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif
# endif

# ifdef ICE_MODEL
        IF (Aout(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avguice(i,j)=AVERAGE(ng)%avguice(i,j)+        &
     &                                 ICE(ng)%ui(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout(idVice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgvice(i,j)=AVERAGE(ng)%avgvice(i,j)+        &
     &                                 ICE(ng)%vi(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgaice(i,j)=AVERAGE(ng)%avgaice(i,j)+        &
     &                                 ICE(ng)%ai(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avghice(i,j)=AVERAGE(ng)%avghice(i,j)+        &
     &                                 ICE(ng)%hi(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avghsno(i,j)=AVERAGE(ng)%avghsno(i,j)+        &
     &                                 ICE(ng)%hsn(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgtice(i,j)=AVERAGE(ng)%avgtice(i,j)+        &
     &                                 ICE(ng)%tis(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgtimid(i,j)=AVERAGE(ng)%avgtimid(i,j)+      &
     &                                 ICE(ng)%ti(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idSfwat,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsfwat(i,j)=AVERAGE(ng)%avgsfwat(i,j)+      &
     &                                 ICE(ng)%sfwat(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgageice(i,j)=AVERAGE(ng)%avgageice(i,j)+    &
     &                                 ICE(ng)%ageice(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgiomflx(i,j)=AVERAGE(ng)%avgiomflx(i,j)+    &
     &                                 ICE(ng)%io_mflux(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsig11(i,j)=AVERAGE(ng)%avgsig11(i,j)+      &
     &                                 ICE(ng)%sig11(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsig12(i,j)=AVERAGE(ng)%avgsig12(i,j)+      &
     &                                 ICE(ng)%sig12(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsig22(i,j)=AVERAGE(ng)%avgsig22(i,j)+      &
     &                                 ICE(ng)%sig22(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgT0mk(i,j)=AVERAGE(ng)%avgT0mk(i,j)+        &
     &                                 ICE(ng)%t0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgS0mk(i,j)=AVERAGE(ng)%avgS0mk(i,j)+        &
     &                                 ICE(ng)%s0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWfr(i,j)=AVERAGE(ng)%avgWfr(i,j)+          &
     &                                 ICE(ng)%wfr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWai(i,j)=AVERAGE(ng)%avgWai(i,j)+          &
     &                                 ICE(ng)%wai(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWao(i,j)=AVERAGE(ng)%avgWao(i,j)+          &
     &                                 ICE(ng)%wao(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWio(i,j)=AVERAGE(ng)%avgWio(i,j)+          &
     &                                 ICE(ng)%wio(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWro(i,j)=AVERAGE(ng)%avgWro(i,j)+          &
     &                                 ICE(ng)%wro(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgutau_iw(i,j)=AVERAGE(ng)%avgutau_iw(i,j)+  &
     &                                 ICE(ng)%utau_iw(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgchu_iw(i,j)=AVERAGE(ng)%avgchu_iw(i,j)+    &
     &                                 ICE(ng)%chu_iw(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef ICE_BIO
#  if defined BERING_10K
        IF (Aout(idIcePhL,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgIPhL(i,j)=AVERAGE(ng)%avgIPhL(i,j)+        &
     &                                 ICE(ng)%IcePhL(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIceNO3,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgINO3(i,j)=AVERAGE(ng)%avgINO3(i,j)+        &
     &                                 ICE(ng)%IceNO3(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIceNH4,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgINH4(i,j)=AVERAGE(ng)%avgINH4(i,j)+        &
     &                                 ICE(ng)%IceNH4(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout(idIceLog,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgILog(i,j)=AVERAGE(ng)%avgILog(i,j)+        &
     &                                 ICE(ng)%IceLog(i,j,Iout)
            END DO
          END DO
        END IF
#  endif
# endif

# ifdef BEST_NPZ
!--------------------------------------------
! values are accumulated in the bestnpz.h file
! so do not sum them here
!--------------------------------------------
#  ifdef STATIONARY
         DO itrc=1,NTS(ng)
           DO k=1,N(ng)
             DO j=JstrR,JendR
               DO i=IstrR,IendR
                 AVERAGE(ng)%avgst(i,j,k,itrc)=                         &
     &                     OCEAN(ng)%st(i,j,k,Nout,itrc)
               END DO
             END DO
           END DO
         END DO
#endif
#ifdef PROD3
         DO itrc=1,NPT3(ng)
           DO k=1,N(ng)
             DO j=JstrR,JendR
               DO i=IstrR,IendR
                 AVERAGE(ng)%avgpt3(i,j,k,itrc)=                        &
     &                     OCEAN(ng)%pt3(i,j,k,Nout,itrc)
               END DO
             END DO
           END DO
         END DO
#endif
#ifdef PROD2
         DO itrc=1,NPT2(ng)
             DO j=JstrR,JendR
               DO i=IstrR,IendR
                 AVERAGE(ng)%avgpt2(i,j,itrc)=                          &
     &                     OCEAN(ng)%pt2(i,j,Nout,itrc)
               END DO
             END DO
         END DO
#endif
#ifdef STATIONARY2
         DO itrc=1,NTS2(ng)
!          DO k=1,N(ng)
             DO j=JstrR,JendR
               DO i=IstrR,IendR
                 AVERAGE(ng)%avgst2(i,j,itrc)=                          &
     &                     OCEAN(ng)%st2(i,j,Nout,itrc)
               END DO
             END DO
!          END DO
         END DO
#endif
#ifdef BENTHIC
          DO itrc=1,NBeT(ng)
            DO k=1,NBL(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgbt(i,j,k,itrc)=                        &
     &                   AVERAGE(ng)%avgbt(i,j,k,itrc)+                 &
     &                   OCEAN(ng)%bt(i,j,k,Nout,itrc)
                END DO
              END DO
            END DO
          END DO

#endif
#if defined CLIM_ICE_1D
          DO itrc=1,NIceT(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgit(i,j,itrc)=avgit(i,j,itrc)+            &
     &                           OCEAN(ng)%it(i,j,Nout,itrc)

              END DO
            END DO
          END DO
#  endif
#  ifdef BIOFLUX
         DO itrc=1,NT(ng)
           DO itrc2=1,NT(ng)
            AVERAGE(ng)%avgbflx(itrc,itrc2)=OCEAN(ng)%bflx(itrc,itrc2)
           END DO
         END DO
#  endif
# endif

# ifdef BIO_GOANPZ
!--------------------------------------------
! values are accumulated in the goanpz.h file
! so do not sum them here
!--------------------------------------------
         DO itrc=1,NTS(ng)
           DO k=1,N(ng)
             DO j=JstrR,JendR
               DO i=IstrR,IendR
                 AVERAGE(ng)%avgst(i,j,k,itrc)=                         &
     &                       OCEAN(ng)%st(i,j,k,Nout,itrc)
               END DO
             END DO
           END DO
         END DO
# endif

# ifdef NEARSHORE_MELLOR
!
!  Accumulate radiation stresses.
!
        IF (Aout(idU2Sd,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2Sd(i,j)=AVERAGE(ng)%avgu2Sd(i,j)+        &
     &                                 OCEAN(ng)%ubar_stokes(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2Sd,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2Sd(i,j)=AVERAGE(ng)%avgv2Sd(i,j)+        &
     &                                 OCEAN(ng)%vbar_stokes(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idU2rs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2RS(i,j)=AVERAGE(ng)%avgu2RS(i,j)+        &
     &                                 MIXING(ng)%rustr2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2rs,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2RS(i,j)=AVERAGE(ng)%avgv2RS(i,j)+        &
     &                                 MIXING(ng)%rvstr2d(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idW2xx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxx2d(i,j)=AVERAGE(ng)%avgSxx2d(i,j)+      &
     &                                  MIXING(ng)%Sxx_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idW2xy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxy2d(i,j)=AVERAGE(ng)%avgSxy2d(i,j)+      &
     &                                  MIXING(ng)%Sxy_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idW2yy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSyy2d(i,j)=AVERAGE(ng)%avgSyy2d(i,j)+      &
     &                                  MIXING(ng)%Syy_bar(i,j)
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idU3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3Sd(i,j,k)=AVERAGE(ng)%avgu3Sd(i,j,k)+  &
     &                                     OCEAN(ng)%u_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3Sd(i,j,k)=AVERAGE(ng)%avgv3Sd(i,j,k)+  &
     &                                     OCEAN(ng)%v_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idU3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3RS(i,j,k)=AVERAGE(ng)%avgu3RS(i,j,k)+  &
     &                                     MIXING(ng)%rustr3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3RS(i,j,k)=AVERAGE(ng)%avgv3RS(i,j,k)+  &
     &                                     MIXING(ng)%rvstr3d(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idW3xx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxx3d(i,j,k)=AVERAGE(ng)%avgSxx3d(i,j,k)+&
     &                                      MIXING(ng)%Sxx(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3xy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxy3d(i,j,k)=AVERAGE(ng)%avgSxy3d(i,j,k)+&
     &                                      MIXING(ng)%Sxy(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3yy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSyy3d(i,j,k)=AVERAGE(ng)%avgSyy3d(i,j,k)+&
     &                                      MIXING(ng)%Syy(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3zx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzx3d(i,j,k)=AVERAGE(ng)%avgSzx3d(i,j,k)+&
     &                                      MIXING(ng)%Szx(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3zy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzy3d(i,j,k)=AVERAGE(ng)%avgSzy3d(i,j,k)+&
     &                                      MIXING(ng)%Szy(i,j,k)
              END DO
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Accumulate vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=AVERAGE(ng)%avgpvor2d(i,j)+    &
     &                                   potvor_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=AVERAGE(ng)%avgrvor2d(i,j)+    &
     &                                   relvor_bar(i,j)
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=AVERAGE(ng)%avgpvor3d(i,j, &
     &                                                             k)+  &
     &                                       potvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=AVERAGE(ng)%avgrvor3d(i,j, &
     &                                                             k)+  &
     &                                       relvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
# endif
!
!  Accumulate quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=AVERAGE(ng)%avgZZ(i,j)+            &
     &                               OCEAN(ng)%zeta(i,j,Kout)*          &
     &                               OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=AVERAGE(ng)%avgU2(i,j)+            &
     &                               OCEAN(ng)%ubar(i,j,Kout)*          &
     &                               OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=AVERAGE(ng)%avgV2(i,j)+            &
     &                               OCEAN(ng)%vbar(i,j,Kout)*          &
     &                               OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=AVERAGE(ng)%avgUU(i,j,k)+      &
     &                                   OCEAN(ng)%u(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=AVERAGE(ng)%avgVV(i,j,k)+      &
     &                                   OCEAN(ng)%v(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=AVERAGE(ng)%avgUV(i,j,k)+      &
     &                                   0.25_r8*                       &
     &                                   (OCEAN(ng)%u(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%u(i+1,j  ,k,Nout))* &
     &                                   (OCEAN(ng)%v(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%v(i  ,j+1,k,Nout))
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=AVERAGE(ng)%avgHuon(i,j,k)+  &
     &                                     GRID(ng)%Huon(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=AVERAGE(ng)%avgHvom(i,j,k)+  &
     &                                     GRID(ng)%Hvom(i,j,k)
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=AVERAGE(ng)%avgTT(i,j,k,  &
     &                                                          it)+    &
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)*     &
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=AVERAGE(ng)%avgUT(i,j,k,  &
     &                                                          it)+    &
     &                                        0.5_r8*                   &
     &                                        OCEAN(ng)%u(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i-1,j,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i  ,j,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=AVERAGE(ng)%avgVT(i,j,k,  &
     &                                                          it)+    &
     &                                        0.5_r8*                   &
     &                                        OCEAN(ng)%v(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i,j-1,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i,j  ,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF

          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=AVERAGE(ng)%avgHuonT(i,&
     &                                                       j,k,it)+   &
     &                                           0.5_r8*                &
     &                                           GRID(ng)%Huon(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i-1,j,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i  ,j,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=AVERAGE(ng)%avgHvomT(i,&
     &                                                       j,k,it)+   &
     &                                           0.5_r8*                &
     &                                           GRID(ng)%Hvom(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i,j-1,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i,j  ,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsAVG(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
        fac=1.0_r8/REAL(nAVG(ng),r8)
        IF (SOUTH_WEST_TEST) THEN
          AVGtime(ng)=AVGtime(ng)+REAL(nAVG(ng),r8)*dt(ng)
        END IF
!
!  Process state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=fac*AVERAGE(ng)%avgzeta(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=fac*AVERAGE(ng)%avgu2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=fac*AVERAGE(ng)%avgv2d(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=fac*AVERAGE(ng)%avgu3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=fac*AVERAGE(ng)%avgv3d(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=fac*AVERAGE(ng)%avgw3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=fac*                         &
     &                                     AVERAGE(ng)%avgwvel(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=fac*AVERAGE(ng)%avgrho(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=fac*                       &
     &                                       AVERAGE(ng)%avgt(i,j,k,it)
                END DO
              END DO
            END DO
          END IF
        END DO

#  ifdef NEMURO_SED1
        IF (Aout(idPONsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgPONsed(i,j)=fac*AVERAGE(ng)%avgPONsed(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idOPALsed,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgOPALsed(i,j)=fac*                          &
     &                                    AVERAGE(ng)%avgOPALsed(i,j)
            END DO
          END DO
        END IF
#  endif

#  if defined LMD_MIXING || defined MY25_MIXING || defined GLS_MIXING
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=fac*AVERAGE(ng)%avgAKv(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=fac*AVERAGE(ng)%avgAKt(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idSdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKs(i,j,k)=fac*AVERAGE(ng)%avgAKs(i,j,k)
              END DO
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_SKPP
        IF (Aout(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghsbl(i,j)=fac*AVERAGE(ng)%avghsbl(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (Aout(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avghbbl(i,j)=fac*AVERAGE(ng)%avghbbl(i,j)
            END DO
          END DO
        END IF
#  endif
# endif

# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Process 2D/3D coupling terms.
!
        IF (Aout(idUfx1,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg1(i,j)=fac*                          &
     &                                    AVERAGE(ng)%avgDU_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUfx2,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgDU_avg2(i,j)=fac*                          &
     &                                    AVERAGE(ng)%avgDU_avg2(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idVfx1,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg1(i,j)=fac*                          &
     &                                    AVERAGE(ng)%avgDV_avg1(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVfx2,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgDV_avg2(i,j)=fac*                          &
     &                                    AVERAGE(ng)%avgDV_avg2(i,j)
            END DO
          END DO
        END IF
# endif
!
!  Process surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=fac*AVERAGE(ng)%avgsus(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=fac*AVERAGE(ng)%avgsvs(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=fac*AVERAGE(ng)%avgbus(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=fac*AVERAGE(ng)%avgbvs(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=fac*AVERAGE(ng)%avgstf(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgswf(i,j)=fac*AVERAGE(ng)%avgswf(i,j)
            END DO
          END DO
        END IF
#  ifdef SSSFLX
        IF (Aout(idSSSf,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsssflx(i,j)=fac*AVERAGE(ng)%avgsssflx(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef SHORTWAVE
        IF (Aout(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsrf(i,j)=fac*AVERAGE(ng)%avgsrf(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglhf(i,j)=fac*AVERAGE(ng)%avglhf(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgshf(i,j)=fac*AVERAGE(ng)%avgshf(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avglrf(i,j)=fac*AVERAGE(ng)%avglrf(i,j)
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgevap(i,j)=fac*AVERAGE(ng)%avgevap(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgrain(i,j)=fac*AVERAGE(ng)%avgrain(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif
# endif

#  ifdef ICE_MODEL
        IF (Aout(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avguice(i,j)=fac*AVERAGE(ng)%avguice(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgvice(i,j)=fac*AVERAGE(ng)%avgvice(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgaice(i,j)=fac*AVERAGE(ng)%avgaice(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avghice(i,j)=fac*AVERAGE(ng)%avghice(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avghsno(i,j)=fac*AVERAGE(ng)%avghsno(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgtice(i,j)=fac*AVERAGE(ng)%avgtice(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgtimid(i,j)=fac*AVERAGE(ng)%avgtimid(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idSfwat,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsfwat(i,j)=fac*AVERAGE(ng)%avgsfwat(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgageice(i,j)=fac*AVERAGE(ng)%avgageice(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgiomflx(i,j)=fac*AVERAGE(ng)%avgiomflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsig11(i,j)=fac*AVERAGE(ng)%avgsig11(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsig12(i,j)=fac*AVERAGE(ng)%avgsig12(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsig22(i,j)=fac*AVERAGE(ng)%avgsig22(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgT0mk(i,j)=fac*AVERAGE(ng)%avgT0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgS0mk(i,j)=fac*AVERAGE(ng)%avgS0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWfr(i,j)=fac*AVERAGE(ng)%avgWfr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWai(i,j)=fac*AVERAGE(ng)%avgWai(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWao(i,j)=fac*AVERAGE(ng)%avgWao(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWio(i,j)=fac*AVERAGE(ng)%avgWio(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgWro(i,j)=fac*AVERAGE(ng)%avgWro(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgutau_iw(i,j)=                              &
     &                      fac*AVERAGE(ng)%avgutau_iw(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgchu_iw(i,j)=fac*AVERAGE(ng)%avgchu_iw(i,j)
            END DO
          END DO
        END IF
# endif
# ifdef ICE_BIO
#  ifdef BERING_10K
        IF (Aout(idIcePhL,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgIPhL(i,j)=fac*AVERAGE(ng)%avgIPhL(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idIceNO3,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgINO3(i,j)=fac*AVERAGE(ng)%avgINO3(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idIceNH4,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgINH4(i,j)=fac*AVERAGE(ng)%avgINH4(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idIceLog,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgILog(i,j)=fac*AVERAGE(ng)%avgILog(i,j)
            END DO
          END DO
        END IF
#  endif
# endif

# ifdef BEST_NPZ
#  ifdef BENTHIC
        DO itrc=1,NBeT(ng)
          DO k=1,NBL(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgbt(i,j,k,itrc)=                          &
     &                     fac*AVERAGE(ng)%avgbt(i,j,k,itrc)
              END DO
            END DO
          END DO
        END DO
#  endif
#  ifdef CLIM_ICE_1D
        DO itrc=1,NIceT(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgit(i,j,itrc)=                              &
     &                     fac*AVERAGE(ng)%avgit(i,j,itrc)
            END DO
          END DO
        END DO
#  endif
#  ifdef STATIONARY
        DO itrc=1,NTS(ng)
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgst(i,j,k,itrc)=                          &
     &                     fac*AVERAGE(ng)%avgst(i,j,k,itrc)
              END DO
            END DO
          END DO
        END DO
#  endif
#  ifdef STATIONARY2
        DO itrc=1,NTS2(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgst2(i,j,itrc)=                             &
     &                     fac*AVERAGE(ng)%avgst2(i,j,itrc)
            END DO
          END DO
        END DO
#  endif
# endif

# ifdef NEARSHORE_MELLOR
!
!  Process radiation stresses.
!
        IF (Aout(idU2Sd,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2Sd(i,j)=fac*AVERAGE(ng)%avgu2Sd(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2Sd,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2Sd(i,j)=fac*AVERAGE(ng)%avgv2Sd(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idU2rs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2RS(i,j)=fac*AVERAGE(ng)%avgu2RS(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2rs,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2RS(i,j)=fac*AVERAGE(ng)%avgv2RS(i,j)
            END DO
          END DO
        END IF

        IF (Aout(idW2xx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxx2d(i,j)=fac*AVERAGE(ng)%avgSxx2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idW2xy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSxy2d(i,j)=fac*AVERAGE(ng)%avgSxy2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idW2yy,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgSyy2d(i,j)=fac*AVERAGE(ng)%avgSyy2d(i,j)
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idU3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3Sd(i,j,k)=fac*                         &
     &                                     AVERAGE(ng)%avgu3Sd(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3Sd(i,j,k)=fac*                         &
     &                                     AVERAGE(ng)%avgv3Sd(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idU3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3RS(i,j,k)=fac*                         &
     &                                     AVERAGE(ng)%avgu3RS(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3RS(i,j,k)=fac*                         &
     &                                     AVERAGE(ng)%avgv3RS(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idW3xx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxx3d(i,j,k)=fac*                        &
     &                                      AVERAGE(ng)%avgSxx3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3xy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSxy3d(i,j,k)=fac*                        &
     &                                      AVERAGE(ng)%avgSxy3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3yy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSyy3d(i,j,k)=fac*                        &
     &                                      AVERAGE(ng)%avgSyy3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3zx,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzx3d(i,j,k)=fac*                        &
     &                                      AVERAGE(ng)%avgSzx3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3zy,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgSzy3d(i,j,k)=fac*                        &
     &                                      AVERAGE(ng)%avgSzy3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Process vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=fac*                           &
     &                                   AVERAGE(ng)%avgpvor2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=fac*                           &
     &                                   AVERAGE(ng)%avgrvor2d(i,j)
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=fac*                       &
     &                                      AVERAGE(ng)%avgpvor3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=fac*                       &
     &                                      AVERAGE(ng)%avgrvor3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
# endif
!
!  Process quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=fac*AVERAGE(ng)%avgZZ(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=fac*AVERAGE(ng)%avgU2(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=fac*AVERAGE(ng)%avgV2(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=fac*AVERAGE(ng)%avgUU(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=fac*AVERAGE(ng)%avgVV(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=fac*AVERAGE(ng)%avgUV(i,j,k)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=fac*                         &
     &                                     AVERAGE(ng)%avgHuon(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=fac*                         &
     &                                     AVERAGE(ng)%avgHvom(i,j,k)
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=fac*                      &
     &                                       AVERAGE(ng)%avgTT(i,j,k,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=fac*                      &
     &                                       AVERAGE(ng)%avgUT(i,j,k,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=fac*                      &
     &                                       AVERAGE(ng)%avgVT(i,j,k,it)
                END DO
              END DO
            END DO
          END IF

          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=fac*                   &
     &                                    AVERAGE(ng)%avgHuonT(i,j,k,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=fac*                   &
     &                                    AVERAGE(ng)%avgHvomT(i,j,k,it)
                END DO
              END DO
            END DO
          END IF
        END DO
# endif
      END IF

      RETURN
      END SUBROUTINE set_avg_tile

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!***********************************************************************
      SUBROUTINE set_detide_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            NTC, Kout,                            &
#  ifdef SOLVE3D
     &                            Nout,                                 &
#  endif
     &                            CosOmega, SinOmega,                   &
     &                            CosW_avg, CosW_sum,                   &
     &                            SinW_avg, SinW_sum,                   &
     &                            CosWCosW, SinWSinW, SinWCosW)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_average
      USE mod_ocean
      USE mod_scalars
      USE mod_tides
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kout
#  ifdef SOLVE3D
      integer, intent(in) :: Nout
#  endif
      integer, intent(in) :: NTC
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: CosOmega(:)
      real(r8), intent(in) :: SinOmega(:)
      real(r8), intent(inout) :: CosW_avg(:)
      real(r8), intent(inout) :: CosW_sum(:)
      real(r8), intent(inout) :: SinW_avg(:)
      real(r8), intent(inout) :: SinW_sum(:)
      real(r8), intent(inout) :: CosWCosW(:,:)
      real(r8), intent(inout) :: SinWSinW(:,:)
      real(r8), intent(inout) :: SinWCosW(:,:)
#  else
      real(r8), intent(in) :: CosOmega(NTC)
      real(r8), intent(in) :: SinOmega(NTC)
      real(r8), intent(inout) :: CosW_avg(NTC)
      real(r8), intent(inout) :: CosW_sum(NTC)
      real(r8), intent(inout) :: SinW_avg(NTC)
      real(r8), intent(inout) :: SinW_sum(NTC)
      real(r8), intent(inout) :: CosWCosW(NTC,NTC)
      real(r8), intent(inout) :: SinWSinW(NTC,NTC)
      real(r8), intent(inout) :: SinWCosW(NTC,NTC)
#  endif
!
!  Local variable declarations.
!
      integer :: i, it, j, k
      integer :: NTC2, mk, nk

      integer, dimension(2*NTC+1) :: indx

      real(r8) :: fac, fac1
      real(r8) :: Hsum, d

      real(r8), dimension(0:2*NTC) :: Ak
      real(r8), dimension(0:2*NTC) :: tide_harmonics
      real(r8), dimension(0:2*NTC,0:2*NTC) :: C, Y

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      NTC2=2*NTC

      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.1)).or.                          &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
!
!  Compute least-squares coefficients to detide time-averaged fields.
!  Notice that the coefficients are always accumulated and not
!  re-initialized.  This allows better tidal fit as the simulation
!  progresses.
!
        IF (SOUTH_WEST_TEST) THEN
          Hcount(ng)=Hcount(ng)+1
          DO nk=1,NTC
            SinW_avg(nk)=SinOmega(nk)
            CosW_avg(nk)=CosOmega(nk)
            SinW_sum(nk)=SinW_sum(nk)+SinOmega(nk)
            CosW_sum(nk)=CosW_sum(nk)+CosOmega(nk)
            DO mk=1,NTC
              SinWSinW(mk,nk)=SinWSinW(mk,nk)+SinOmega(mk)*SinOmega(nk)
              CosWCosW(mk,nk)=CosWCosW(mk,nk)+CosOmega(mk)*CosOmega(nk)
              SinWCosW(mk,nk)=SinWCosW(mk,nk)+SinOmega(mk)*CosOmega(nk)
            END DO
          END DO
          tide_harmonics(0)=1.0_r8
          DO nk=1,NTC
            tide_harmonics(nk    )=SinOmega(nk)
            tide_harmonics(nk+NTC)=CosOmega(nk)
          END DO
        END IF
!
!  Initialize.
!
        IF (Aout(idFsuD,ng)) THEN
          DO nk=0,NTC2
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                TIDES(ng)%zeta_tide(i,j,nk)=TIDES(ng)%zeta_tide(i,j,    &
     &                                                          nk)+    &
     &                                      OCEAN(ng)%zeta(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idu2dD,ng)) THEN
          DO nk=0,NTC2
            DO j=JstrR,JendR
              DO i=Istr,IendR
                TIDES(ng)%ubar_tide(i,j,nk)=TIDES(ng)%ubar_tide(i,j,    &
     &                                                          nk)+    &
     &                                      OCEAN(ng)%ubar(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv2dD,ng)) THEN
          DO nk=0,NTC2
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                TIDES(ng)%vbar_tide(i,j,nk)=TIDES(ng)%vbar_tide(i,j,    &
     &                                                          nk)+    &
     &                                      OCEAN(ng)%zeta(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idu3dD,ng)) THEN
          DO nk=0,NTC2
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  TIDES(ng)%u_tide(i,j,k,nk)=TIDES(ng)%u_tide(i,j,k,    &
     &                                                        nk)+      &
     &                                       OCEAN(ng)%u(i,j,k,Nout)*   &
     &                                       tide_harmonics(nk)
                END DO
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv3dD,ng)) THEN
          DO nk=0,NTC2
            DO k=1,N(ng)
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  TIDES(ng)%v_tide(i,j,k,nk)=TIDES(ng)%v_tide(i,j,k,    &
     &                                                        nk)+      &
     &                                       OCEAN(ng)%v(i,j,k,Nout)*   &
     &                                       tide_harmonics(nk)
                END DO
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTrcD(it),ng)) THEN
            DO nk=0,NTC2
              DO k=1,N(ng)
                DO j=Jstr,JendR
                  DO i=IstrR,IendR
                    TIDES(ng)%t_tide(i,j,k,nk,it)=TIDES(ng)%t_tide(i,j, &
     &                                                        k,nk,it)+ &
     &                                            OCEAN(ng)%t(i,j,k,    &
     &                                                        Nout,it)* &
     &                                            tide_harmonics(nk)
                  END DO
                END DO
              END DO
            END DO
          END IF
        END DO
#  endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG(ng)) THEN
!
!  Accumukate Detide least-squares coefficients. They only vary in time
!  since omega (as computed in set_tides) uses model time coordinate.
!
        IF (SOUTH_WEST_TEST) THEN
          Hcount(ng)=Hcount(ng)+1
          DO nk=1,NTC
            SinW_avg(nk)=SinW_avg(nk)+SinOmega(nk)
            CosW_avg(nk)=CosW_avg(nk)+CosOmega(nk)
            SinW_sum(nk)=SinW_sum(nk)+SinOmega(nk)
            CosW_sum(nk)=CosW_sum(nk)+CosOmega(nk)
            DO mk=1,NTC
              SinWSinW(mk,nk)=SinWSinW(mk,nk)+SinOmega(mk)*SinOmega(nk)
              CosWCosW(mk,nk)=CosWCosW(mk,nk)+CosOmega(mk)*CosOmega(nk)
              SinWCosW(mk,nk)=SinWCosW(mk,nk)+SinOmega(mk)*CosOmega(nk)
            END DO
          END DO
          tide_harmonics(0)=1.0_r8
          DO nk=1,NTC
            tide_harmonics(nk    )=SinOmega(nk)
            tide_harmonics(nk+NTC)=CosOmega(nk)
          END DO
        END IF
!
!  Accumulate.
!
        IF (Aout(idFsuD,ng)) THEN
          DO nk=0,NTC2
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                TIDES(ng)%zeta_tide(i,j,nk)=TIDES(ng)%zeta_tide(i,j,    &
     &                                                          nk)+    &
     &                                      OCEAN(ng)%zeta(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idu2dD,ng)) THEN
          DO nk=0,NTC2
            DO j=JstrR,JendR
              DO i=Istr,IendR
               TIDES(ng)%ubar_tide(i,j,nk)=TIDES(ng)%ubar_tide(i,j,    &
     &                                                          nk)+    &
     &                                      OCEAN(ng)%ubar(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv2dD,ng)) THEN
          DO nk=0,NTC2
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                TIDES(ng)%vbar_tide(i,j,nk)=TIDES(ng)%vbar_tide(i,j,    &
     &                                                          nk)+    &
     &                                      OCEAN(ng)%zeta(i,j,Kout)*   &
     &                                      tide_harmonics(nk)
              END DO
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idu3dD,ng)) THEN
          DO nk=0,NTC2
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,IendR
                  TIDES(ng)%u_tide(i,j,k,nk)=TIDES(ng)%u_tide(i,j,k,    &
     &                                                        nk)+      &
     &                                       OCEAN(ng)%u(i,j,k,Nout)*   &
     &                                       tide_harmonics(nk)
                END DO
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv3dD,ng)) THEN
          DO nk=0,NTC2
            DO k=1,N(ng)
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  TIDES(ng)%v_tide(i,j,k,nk)=TIDES(ng)%v_tide(i,j,k,    &
     &                                                        nk)+      &
     &                                       OCEAN(ng)%v(i,j,k,Nout)*   &
     &                                       tide_harmonics(nk)
                END DO
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTrcD(it),ng)) THEN
            DO nk=0,NTC2
              DO k=1,N(ng)
                DO j=Jstr,JendR
                  DO i=IstrR,IendR
                    TIDES(ng)%t_tide(i,j,k,nk,it)=TIDES(ng)%t_tide(i,j, &
     &                                                        k,nk,it)+ &
     &                                            OCEAN(ng)%t(i,j,k,    &
     &                                                        Nout,it)* &
     &                                            tide_harmonics(nk)
                  END DO
                END DO
              END DO
            END DO
          END IF
        END DO
#  endif
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsAVG(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
        fac=1.0_r8/REAL(nAVG(ng),r8)
!
!  Compute detide least-squares coefficients.  Build coefficient squared
!  matrix C(0:2*NTC,0:2*NTC) to invert. It is 2*NTC because we are
!  solving for real components Ak and Bk. The zero rows and column is
!  for the coefficients associated with the time mean.
!
!        F(t) = Fmean + SUM [ Ak sin(omega(k)*t) ]
!                     + SUM [ Bk cos(omega(k)*t) ]   for k=1:NTC
!
!  In the code below, all the arrays are collapsed into a single
!  dimension index such that:
!
!            k=0               mean term
!            k=1:NTC           sine terms
!            k=NTC+1:2*NTC     cosine terms
!
        IF (SOUTH_WEST_TEST) THEN
          C(0,0)=1.0_r8                     ! time-averaged coefficient
          fac1=1.0_r8/REAL(Hcount(ng),r8)   ! global summation factor
          DO nk=1,NTC
            C(0,nk    )=fac1*SinW_sum(nk)
            C(0,nk+NTC)=fac1*CosW_sum(nk)
            C(nk,0    )=C(0,nk)             ! symmetric
            C(nk+NTC,0)=C(0,nk+NTC)         ! symmetric
            DO mk=1,NTC
              C(mk,nk)=fac1*SinWSinW(mk,nk)
              C(mk,nk+NTC)=fac1*SinWCosW(mk,nk)
              C(mk+NTC,nk)=fac1*SinWCosW(nk,mk)
              C(mk+NTC,nk+NTC)=fac1*CosWCosW(mk,nk)
            END DO
          END DO
          DO nk=0,NTC2
            DO mk=0,NTC2
              C(nk,mk)=C(mk,nk)
            END DO
          END DO
!
!  Invert least-squares coefficient matrix by LU decomposition.
!
          DO mk=0,NTC2
            DO nk=0,NTC2
              Y(mk,nk)=0.0_r8
            END DO
            Y(mk,mk)=1.0_r8                ! identity matrix
          END DO
          CALL ludcmp (C(0,0), NTC2+1, NTC2+1, indx, d)
!
!  Find inverse by columns. The matrix Y will now contain the inverse
!  of the least-squares coefficient matrix C, which will have been
!  destroyed.
!
          DO nk=0,NTC2
            CALL lubksb (C(0,0), NTC2+1, NTC2+1, indx, Y(0,nk))
          END DO
!
!  Compute time-averaged harmonics for current field average window.
!
          tide_harmonics(0)=1.0_r8
          DO nk=1,NTC
            tide_harmonics(nk    )=fac*SinW_avg(nk)
            tide_harmonics(nk+NTC)=fac*cosW_avg(nk)
          END DO
!
!  Scale inverse by the global summation factor.
!
          DO nk=0,NTC2
            DO mk=0,NTC2
              Y(nk,mk)=fac1*Y(nk,mk)
            END DO
          END DO
        END IF
!
!  Process accumulated detided averages. Notice that the regular
!  time-averaged field values used here are the ones computed in
!  routine "set_avg".
!
        IF (Aout(idFsuD,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              Hsum=0.0_r8
              DO nk=0,NTC2
                Ak(nk)=0.0_r8
                DO mk=0,NTC2
!!                Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%zeta_tide(i,j,mk)
                  Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%zeta_tide(i,j,mk)
                END DO
                Hsum=Hsum+Ak(nk)*tide_harmonics(nk)
              END DO
              TIDES(ng)%zeta_detided(i,j)=AVERAGE(ng)%avgzeta(i,j)-     &
     &                                    Hsum
            END DO
          END DO
        END IF

        IF (Aout(idu2dD,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              Hsum=0.0_r8
              DO nk=0,NTC2
                Ak(nk)=0.0_r8
                DO mk=0,NTC2
!!                Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%ubar_tide(i,j,mk)
                  Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%ubar_tide(i,j,mk)
                END DO
                Hsum=Hsum+Ak(nk)*tide_harmonics(nk)
              END DO
              TIDES(ng)%ubar_detided(i,j)=AVERAGE(ng)%avgu2d(i,j)-      &
     &                                    Hsum
            END DO
          END DO
        END IF

        IF (Aout(idv2dD,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              Hsum=0.0_r8
              DO nk=0,NTC2
                Ak(nk)=0.0_r8
                DO mk=0,NTC2
!!                Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%vbar_tide(i,j,mk)
                  Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%vbar_tide(i,j,mk)
                END DO
                Hsum=Hsum+Ak(nk)*tide_harmonics(nk)
              END DO
              TIDES(ng)%vbar_detided(i,j)=AVERAGE(ng)%avgv2d(i,j)-      &
     &                                    Hsum
            END DO
          END DO
        END IF

#  ifdef SOLVE3D
        IF (Aout(idu3dD,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                Hsum=0.0_r8
                DO nk=0,NTC2
                  Ak(nk)=0.0_r8
                  DO mk=0,NTC2
!!                  Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%u_tide(i,j,k,mk)
                    Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%u_tide(i,j,k,mk)
                  END DO
                  Hsum=Hsum+Ak(nk)*tide_harmonics(nk)
                END DO
                TIDES(ng)%u_detided(i,j,k)=AVERAGE(ng)%avgu3d(i,j,k)-   &
     &                                     Hsum
              END DO
            END DO
          END DO
        END IF

        IF (Aout(idv3dD,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                Hsum=0.0_r8
                DO nk=0,NTC2
                  Ak(nk)=0.0_r8
                  DO mk=0,NTC2
!!                  Ak(nk)=Ak(nk)+Y(nk,mk)*TIDES(ng)%v_tide(i,j,k,mk)
                    Ak(nk)=Ak(nk)+Y(mk,nk)*TIDES(ng)%v_tide(i,j,k,mk)
                  END DO
                  Hsum=Hsum+Ak(nk)*tide_harmonics(nk)
                END DO
                TIDES(ng)%v_detided(i,j,k)=AVERAGE(ng)%avgv3d(i,j,k)-   &
     &                                     Hsum
              END DO
            END DO
          END DO
        END IF

        DO it=1,NAT
          IF (Aout(idTrcD(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,JendR
                DO i=IstrR,IendR
                  Hsum=0.0_r8
                  DO nk=0,NTC2
                    Ak(nk)=0.0_r8
                    DO mk=0,NTC2
!!                    Ak(nk)=Ak(nk)+                                    &
!!   &                       Y(nk,mk)*TIDES(ng)%t_tide(i,j,k,mk,it)
                      Ak(nk)=Ak(nk)+                                    &
     &                       Y(mk,nk)*TIDES(ng)%t_tide(i,j,k,mk,it)
                    END DO
                    Hsum=Hsum+Ak(nk)*tide_harmonics(nk)
                  END DO
                  TIDES(ng)%t_detided(i,j,k,it)=AVERAGE(ng)%avgt(i,j,k, &
     &                                                           it)-   &
     &                                          Hsum
                END DO
              END DO
            END DO
          END IF
        END DO
      END IF
#  endif
      RETURN
      END SUBROUTINE set_detide_tile
# endif
#endif
      END MODULE set_avg_mod
