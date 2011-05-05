#include "cppdefs.h"
      MODULE set_avg2_mod
#if defined AVERAGES2 && (!defined ADJOINT && defined NONLINEAR)
!
!svn $Id: set_avg.F 702 2008-08-12 16:44:47Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  fields.  Due to synchronization, the time-averaged fields are       !
!  computed in delayed mode. All averages are accumulated at the       !
!  beggining of the next time-step.                                    !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: set_avg2

      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_avg2 (ng, tile)
!***********************************************************************
!
      USE mod_param
# if defined FORWARD_WRITE && defined SOLVE3D
      USE mod_coupling
# endif
      USE mod_stepping
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
      CALL set_avg2_tile (ng, tile,                                     &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   KOUT,                                          &
# ifdef ICE_MODEL
     &                   IOUT, IUOUT, IEOUT,                            &
# endif
     &                   NOUT)

# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 5)
# endif
      RETURN
      END SUBROUTINE set_avg2
!
!***********************************************************************
      SUBROUTINE set_avg2_tile (ng, tile,                               &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Kout,                                    &
# ifdef ICE_MODEL
     &                         Iout, Iuout, Ieout,                      &
# endif
     &                         Nout)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_average2
      USE mod_forces
# ifdef SOLVE3D
      USE mod_grid
# endif
# ifdef ICE_MODEL
      USE mod_ice
# endif
      USE mod_mixing
      USE mod_ocean
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
      integer :: i, itrc, j

      real(r8) :: fac, fac1

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG2(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG2(ng)).and.                                &
     &     (MOD(iic(ng)-1,nAVG2(ng)).eq.1)).or.                         &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
!
!  Initialize 2D fields.
!
        IF (Aout2(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgzeta(i,j)=OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout2(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgu2d (i,j)=OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout2(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgv2d (i,j)=OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout2(idUvel,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgu3d(i,j)=OCEAN(ng)%u(i,j,N(ng),Nout)
            END DO
          END DO
        END IF
        IF (Aout2(idVvel,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgv3d(i,j)=OCEAN(ng)%v(i,j,N(ng),Nout)
            END DO
          END DO
        END IF

        IF (Aout2(idDano,ng)) THEN
          DO j=JstrR,JendR 
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgrho(i,j)=OCEAN(ng)%rho(i,j,N(ng))
            END DO
          END DO
        END IF
        DO itrc=1,NT(ng)
          IF (Aout2(idTvar(itrc),ng)) THEN
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE2(ng)%avgt(i,j,itrc)=                            &
     &                    OCEAN(ng)%t(i,j,N(ng),Nout,itrc)
              END DO
            END DO
          END IF
        END DO
#  ifdef LMD_SKPP
        IF (Aout2(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghsbl(i,j)=MIXING(ng)%hsbl(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (Aout2(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghbbl(i,j)=MIXING(ng)%hbbl(i,j)
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Initialize surface and bottom fluxes.
!
        IF (Aout2(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsus(i,j)=FORCES(ng)%sustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsvs(i,j)=FORCES(ng)%svstr(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgbus(i,j)=FORCES(ng)%bustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgbvs(i,j)=FORCES(ng)%bvstr(i,j)
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout2(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgstf(i,j)=FORCES(ng)%stflx(i,j,itemp)
            END DO
          END DO
        END IF
        IF (Aout2(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgswf(i,j)=FORCES(ng)%stflx(i,j,isalt)
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout2(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsrf(i,j)=FORCES(ng)%srflx(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout2(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avglhf(i,j)=FORCES(ng)%lhflx(i,j)
            END DO
          END DO
        END IF 
        IF (Aout2(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avglrf(i,j)=FORCES(ng)%lrflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgshf(i,j)=FORCES(ng)%shflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avguwind(i,j)=FORCES(ng)%Uwind(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgvwind(i,j)=FORCES(ng)%Vwind(i,j)
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout2(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgevap(i,j)=FORCES(ng)%evap(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgrain(i,j)=FORCES(ng)%rain(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif
#  ifdef ICE_MODEL
        IF (Aout2(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avguice(i,j)=ICE(ng)%ui(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout2(idUice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgvice(i,j)=ICE(ng)%vi(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout2(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgaice(i,j)=ICE(ng)%ai(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout2(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghice(i,j)=ICE(ng)%hi(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout2(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghsno(i,j)=ICE(ng)%hsn(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout2(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgtice(i,j)=ICE(ng)%tis(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgtimid(i,j)=ICE(ng)%ti(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout2(idSfwat,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsfwat(i,j)=ICE(ng)%sfwat(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout2(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgageice(i,j)=ICE(ng)%ageice(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout2(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgiomflx(i,j)=ICE(ng)%io_mflux(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsig11(i,j)=ICE(ng)%sig11(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout2(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsig12(i,j)=ICE(ng)%sig12(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout2(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsig22(i,j)=ICE(ng)%sig22(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout2(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgT0mk(i,j)=ICE(ng)%t0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgS0mk(i,j)=ICE(ng)%s0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWfr(i,j)=ICE(ng)%wfr(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWai(i,j)=ICE(ng)%wai(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWao(i,j)=ICE(ng)%wao(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWio(i,j)=ICE(ng)%wio(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgWro(i,j)=ICE(ng)%wro(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgutau_iw(i,j)=ICE(ng)%utau_iw(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgchu_iw(i,j)=ICE(ng)%chu_iw(i,j)
            END DO
          END DO
        END IF
#  endif
# endif

# ifdef SOLVE3D
!
!  Initialized fields associated with tracers.
!
#  ifdef BIO_GOANPZ
        DO itrc=1,NTS(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgst(i,j,itrc)=                             &
     &                       OCEAN(ng)%st(i,j,N(ng),Nout,itrc)
            END DO
          END DO
        END DO
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG2(ng)) THEN
!
!  Accumulate 2D fields.
!
        IF (Aout2(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgzeta(i,j)=AVERAGE2(ng)%avgzeta(i,j)+      &
     &                                  OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout2(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgu2d (i,j)=AVERAGE2(ng)%avgu2d (i,j)+      &
     &                                  OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout2(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgv2d (i,j)=AVERAGE2(ng)%avgv2d (i,j)+      &
     &                                  OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout2(idUvel,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgu3d(i,j)=AVERAGE2(ng)%avgu3d(i,j)+        &
     &                                  OCEAN(ng)%u(i,j,N(ng),Nout)
            END DO
          END DO
        END IF
        IF (Aout2(idVvel,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgv3d(i,j)=AVERAGE2(ng)%avgv3d(i,j)+        &
     &                                  OCEAN(ng)%v(i,j,N(ng),Nout)
            END DO
          END DO
        END IF

        IF (Aout2(idDano,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgrho(i,j)=AVERAGE2(ng)%avgrho(i,j)+        &
     &                                  OCEAN(ng)%rho(i,j,N(ng))
            END DO
          END DO
        END IF
        DO itrc=1,NT(ng)
          IF (Aout2(idTvar(itrc),ng)) THEN
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE2(ng)%avgt(i,j,itrc)=AVERAGE2(ng)%avgt(i,j,itrc)+&
     &                                  OCEAN(ng)%t(i,j,N(ng),Nout,itrc)
              END DO
            END DO
          END IF
        END DO
#  ifdef LMD_SKPP
        IF (Aout2(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghsbl(i,j)=AVERAGE2(ng)%avghsbl(i,j)+      &
     &                                  MIXING(ng)%hsbl(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (Aout2(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghbbl(i,j)=AVERAGE2(ng)%avghbbl(i,j)+      &
     &                                  MIXING(ng)%hbbl(i,j)
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Accumulate surface and bottom fluxes.
!
        IF (Aout2(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsus(i,j)=AVERAGE2(ng)%avgsus(i,j)+        &
     &                                 FORCES(ng)%sustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsvs(i,j)=AVERAGE2(ng)%avgsvs(i,j)+        &
     &                                 FORCES(ng)%svstr(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgbus(i,j)=AVERAGE2(ng)%avgbus(i,j)+        &
     &                                 FORCES(ng)%bustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgbvs(i,j)=AVERAGE2(ng)%avgbvs(i,j)+        &
     &                                 FORCES(ng)%bvstr(i,j)
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout2(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgstf(i,j)=AVERAGE2(ng)%avgstf(i,j)+        &
     &                                  FORCES(ng)%stflx(i,j,itemp)
            END DO
          END DO
        END IF
        IF (Aout2(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgswf(i,j)=AVERAGE2(ng)%avgswf(i,j)+        &
     &                                  FORCES(ng)%stflx(i,j,isalt)
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout2(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsrf(i,j)=AVERAGE2(ng)%avgsrf(i,j)+        &
     &                                  FORCES(ng)%srflx(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout2(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avglhf(i,j)=AVERAGE2(ng)%avglhf(i,j)+        &
     &                                  FORCES(ng)%lhflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgshf(i,j)=AVERAGE2(ng)%avgshf(i,j)+        &
     &                                  FORCES(ng)%shflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avglrf(i,j)=AVERAGE2(ng)%avglrf(i,j)+        &
     &                                  FORCES(ng)%lrflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avguwind(i,j)=AVERAGE2(ng)%avguwind(i,j)+    &
     &                                  FORCES(ng)%Uwind(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgvwind(i,j)=AVERAGE2(ng)%avgvwind(i,j)+    &
     &                                  FORCES(ng)%Vwind(i,j)
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout2(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgevap(i,j)=AVERAGE2(ng)%avgevap(i,j)+      &
     &                                  FORCES(ng)%evap(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgrain(i,j)=AVERAGE2(ng)%avgrain(i,j)+      &
     &                                  FORCES(ng)%rain(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif
#  ifdef ICE_MODEL
        IF (Aout2(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avguice(i,j)=AVERAGE2(ng)%avguice(i,j)+      &
     &                                  ICE(ng)%ui(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout2(idVice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgvice(i,j)=AVERAGE2(ng)%avgvice(i,j)+      &
     &                                  ICE(ng)%vi(i,j,Iuout)
            END DO
          END DO
        END IF
        IF (Aout2(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgaice(i,j)=AVERAGE2(ng)%avgaice(i,j)+      &
     &                                  ICE(ng)%ai(i,j,Iout)
            END DO
          END DO
        END IF 
        IF (Aout2(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avghice(i,j)=AVERAGE2(ng)%avghice(i,j)+      &
     &                                  ICE(ng)%hi(i,j,Iout)
            END DO
          END DO
        END IF 
        IF (Aout2(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avghsno(i,j)=AVERAGE2(ng)%avghsno(i,j)+      &
     &                                  ICE(ng)%hsn(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout2(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgtice(i,j)=AVERAGE2(ng)%avgtice(i,j)+      &
     &                                  ICE(ng)%tis(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgtimid(i,j)=AVERAGE2(ng)%avgtimid(i,j)+    &
     &                                  ICE(ng)%ti(i,j,Iout)
            END DO
          END DO
        END IF 
        IF (Aout2(idSfwat,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsfwat(i,j)=AVERAGE2(ng)%avgsfwat(i,j)+    &
     &                                  ICE(ng)%sfwat(i,j,Iout)
            END DO
          END DO
        END IF 
        IF (Aout2(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgageice(i,j)=AVERAGE2(ng)%avgageice(i,j)+  &
     &                                  ICE(ng)%ageice(i,j,Iout)
            END DO
          END DO
        END IF
        IF (Aout2(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgiomflx(i,j)=AVERAGE2(ng)%avgiomflx(i,j)+  &
     &                                  ICE(ng)%io_mflux(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsig11(i,j)=AVERAGE2(ng)%avgsig11(i,j)+    &
     &                                  ICE(ng)%sig11(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout2(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsig12(i,j)=AVERAGE2(ng)%avgsig12(i,j)+    &
     &                                  ICE(ng)%sig12(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout2(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsig22(i,j)=AVERAGE2(ng)%avgsig22(i,j)+    &
     &                                  ICE(ng)%sig22(i,j,Ieout)
            END DO
          END DO
        END IF
        IF (Aout2(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgT0mk(i,j)=AVERAGE2(ng)%avgT0mk(i,j)+      &
     &                                  ICE(ng)%t0mk(i,j)
            END DO
          END DO
        END IF 
        IF (Aout2(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgS0mk(i,j)=AVERAGE2(ng)%avgS0mk(i,j)+      &
     &                                  ICE(ng)%s0mk(i,j)
            END DO
          END DO
        END IF 
        IF (Aout2(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWfr(i,j)=AVERAGE2(ng)%avgWfr(i,j)+        &
     &                                 ICE(ng)%wfr(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWai(i,j)=AVERAGE2(ng)%avgWai(i,j)+        &
     &                                 ICE(ng)%wai(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWao(i,j)=AVERAGE2(ng)%avgWao(i,j)+        &
     &                                 ICE(ng)%wao(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWio(i,j)=AVERAGE2(ng)%avgWio(i,j)+        &
     &                                 ICE(ng)%wio(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWro(i,j)=AVERAGE2(ng)%avgWro(i,j)+        &
     &                                 ICE(ng)%wro(i,j)
            END DO
          END DO
        END IF 
        IF (Aout2(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgutau_iw(i,j)=AVERAGE2(ng)%avgutau_iw(i,j)+&
     &                                  ICE(ng)%utau_iw(i,j)
            END DO
          END DO
        END IF 
        IF (Aout2(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgchu_iw(i,j)=AVERAGE2(ng)%avgchu_iw(i,j)+  &
     &                                  ICE(ng)%chu_iw(i,j)
            END DO
          END DO
        END IF
#  endif
# endif

# ifdef BIO_GOANPZ
!--------------------------------------------
! values are accumulated in the goanpz.h file
! so do not sum them here
!--------------------------------------------
        DO itrc=1,NTS(ng)
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgst(i,j,itrc)=                             &
     &                       OCEAN(ng)%st(i,j,N(ng),Nout,itrc)
            END DO
          END DO
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsAVG2(ng)).and.                                 &
     &    (MOD(iic(ng)-1,nAVG2(ng)).eq.0).and.                          &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
        fac=1.0_r8/REAL(nAVG2(ng),r8)
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          AVG2time(ng)=AVG2time(ng)+REAL(nAVG2(ng),r8)*dt(ng)
        END IF
!
!  Process state variables.
!
        IF (Aout2(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgzeta(i,j)=fac*AVERAGE2(ng)%avgzeta(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgu2d (i,j)=fac*AVERAGE2(ng)%avgu2d (i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgv2d (i,j)=fac*AVERAGE2(ng)%avgv2d (i,j)
            END DO
          END DO
        END IF
# ifdef SOLVE3D
        IF (Aout2(idUvel,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgu3d(i,j)=fac*AVERAGE2(ng)%avgu3d(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVvel,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgv3d(i,j)=fac*AVERAGE2(ng)%avgv3d(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idDano,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgrho(i,j)=fac*AVERAGE2(ng)%avgrho(i,j)
            END DO
          END DO
        END IF
        DO itrc=1,NT(ng)
          IF (Aout2(idTvar(itrc),ng)) THEN
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE2(ng)%avgt(i,j,itrc)=fac*                        &
     &                      AVERAGE2(ng)%avgt(i,j,itrc)
              END DO
            END DO
          END IF
        END DO
#  ifdef LMD_SKPP
        IF (Aout2(idHsbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghsbl(i,j)=fac*AVERAGE2(ng)%avghsbl(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef LMD_BKPP
        IF (Aout2(idHbbl,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avghbbl(i,j)=fac*AVERAGE2(ng)%avghbbl(i,j)
            END DO
          END DO
        END IF
#  endif
# endif
!
!  Process surface and bottom fluxes.
!
        IF (Aout2(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsus(i,j)=fac*AVERAGE2(ng)%avgsus(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsvs(i,j)=fac*AVERAGE2(ng)%avgsvs(i,j)
            END DO
          END DO
        END IF

        IF (Aout2(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgbus(i,j)=fac*AVERAGE2(ng)%avgbus(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgbvs(i,j)=fac*AVERAGE2(ng)%avgbvs(i,j)
            END DO
          END DO
        END IF

# ifdef SOLVE3D
        IF (Aout2(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgstf(i,j)=fac*AVERAGE2(ng)%avgstf(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idTsur(isalt),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgswf(i,j)=fac*AVERAGE2(ng)%avgswf(i,j)
            END DO
          END DO
        END IF
#  ifdef SHORTWAVE
        IF (Aout2(idSrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgsrf(i,j)=fac*AVERAGE2(ng)%avgsrf(i,j)
            END DO
          END DO
        END IF
#  endif
#  ifdef BULK_FLUXES
        IF (Aout2(idLhea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avglhf(i,j)=fac*AVERAGE2(ng)%avglhf(i,j)
            END DO
          END DO
        END IF 
        IF (Aout2(idShea,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgshf(i,j)=fac*AVERAGE2(ng)%avgshf(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idLrad,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avglrf(i,j)=fac*AVERAGE2(ng)%avglrf(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idUair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avguwind(i,j)=fac*AVERAGE2(ng)%avguwind(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVair,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgvwind(i,j)=fac*AVERAGE2(ng)%avgvwind(i,j)
            END DO
          END DO
        END IF
#   ifdef EMINUSP
        IF (Aout2(idevap,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgevap(i,j)=fac*AVERAGE2(ng)%avgevap(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idrain,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgrain(i,j)=fac*AVERAGE2(ng)%avgrain(i,j)
            END DO
          END DO
        END IF
#   endif
#  endif

#  ifdef ICE_MODEL
        IF (Aout2(idUice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avguice(i,j)=fac*AVERAGE2(ng)%avguice(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idVice,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE2(ng)%avgvice(i,j)=fac*AVERAGE2(ng)%avgvice(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idAice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgaice(i,j)=fac*AVERAGE2(ng)%avgaice(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idHice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avghice(i,j)=fac*AVERAGE2(ng)%avghice(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idHsno,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avghsno(i,j)=fac*AVERAGE2(ng)%avghsno(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idTice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgtice(i,j)=fac*AVERAGE2(ng)%avgtice(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idTimid,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgtimid(i,j)=fac*                           &
     &                      AVERAGE2(ng)%avgtimid(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idSfwat,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsfwat(i,j)=fac*                           &
     &                      AVERAGE2(ng)%avgsfwat(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idAgeice,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgageice(i,j)=fac*                          &
     &                      AVERAGE2(ng)%avgageice(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idIomflx,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgiomflx(i,j)=fac*                          &
     &                      AVERAGE2(ng)%avgiomflx(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idSig11,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsig11(i,j)=fac*                           &
     &                      AVERAGE2(ng)%avgsig11(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idSig12,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsig12(i,j)=fac*                           &
     &                      AVERAGE2(ng)%avgsig12(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idSig22,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgsig22(i,j)=fac*                           &
     &                      AVERAGE2(ng)%avgsig22(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idT0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgT0mk(i,j)=fac*AVERAGE2(ng)%avgT0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idS0mk,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgS0mk(i,j)=fac*AVERAGE2(ng)%avgS0mk(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWfr,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWfr(i,j)=fac*AVERAGE2(ng)%avgWfr(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWai,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWai(i,j)=fac*AVERAGE2(ng)%avgWai(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWao,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWao(i,j)=fac*AVERAGE2(ng)%avgWao(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWio,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWio(i,j)=fac*AVERAGE2(ng)%avgWio(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idWro,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgWro(i,j)=fac*AVERAGE2(ng)%avgWro(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idTauiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgutau_iw(i,j)=fac*                         &
     &                      AVERAGE2(ng)%avgutau_iw(i,j)
            END DO
          END DO
        END IF
        IF (Aout2(idChuiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE2(ng)%avgchu_iw(i,j)=fac*                          &
     &                      AVERAGE2(ng)%avgchu_iw(i,j)
            END DO
          END DO
        END IF
#  endif
# endif

# ifdef SOLVE3D
!#  ifdef BIO_GOANPZ
!-----------------------------------------------------------
! Computes average of time series
! For production, totals are needed so this is commented out
!-----------------------------------------------------------
!        DO itrc=1,NTS(ng)
!            DO j-JstrR,JendR
!              DO i-IstrR,IendR
!                AVERAGE2(ng)%anvst(i,j,itrc)=fac*AVERAGE2(ng)%avgt(i,j,itrc)
!            END DO
!          END DO
!        END DO
!#  endif
# endif
!
      END IF

      RETURN
      END SUBROUTINE set_avg2_tile
#endif
      END MODULE set_avg2_mod
