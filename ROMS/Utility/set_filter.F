#include "cppdefs.h"
#ifdef FILTERED
      SUBROUTINE set_filter(ng, tile)
!
!=======================================================================
!  Copyright (c) 2002-2009 ROMS/TOMS Group                             !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  fields.  Due to synchronization, the time-averaged fields are       !
!  computed in delayed mode. All averages are accumulated at the       !
!  beginning of the next time-step.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
!
      implicit none

      integer, intent(in) :: ng, tile

# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 5)
# endif
      CALL set_filter_tile (ng, tile,                                   &
# ifdef ICE_MODEL
     &                   IOUT, IUOUT, IEOUT,                            &
# endif
     &                   KOUT, NOUT)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 5)
# endif
      RETURN
      END SUBROUTINE set_filter
!
!***********************************************************************
      SUBROUTINE set_filter_tile (ng, tile,                             &
# ifdef ICE_MODEL
     &                         Iout, Iuout, Ieout,                      &
# endif
     &                         Kout, Nout)
!***********************************************************************
!
      USE mod_param
      USE mod_average
      USE mod_filter
      USE mod_forces
      USE mod_grid
# ifdef ICE_MODEL
      USE mod_ice
# endif
      USE mod_mixing
      USE mod_ocean
      USE mod_scalars
      USE mod_parallel
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Kout, Nout
# ifdef ICE_MODEL
      integer, intent(in) :: Iout, Iuout, Ieout
# endif
!
!  Local variable declarations.
!
      integer  :: i, itrc, j, k, iw, f
      real(r8) :: fac

#include "set_bounds.h"

!
!-----------------------------------------------------------------------
! Initialize time-filtered arrays when appropriate.  They aren't all
! initialized at beginning of run - skip right to summing, which is
! OK because it was initialized in init_arrays
!-----------------------------------------------------------------------
!
!  handle tracers other than T,S without filtering
!
!     IF (Master) print *, 'start of SET_FILTER'
      IF (NT(ng) .gt. NAT) THEN
        IF (((iic(ng).gt.ntsAVG(ng)).and.                               &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.1)).or.                          &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
          DO itrc=NAT+1,NT(ng)
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,itrc)=                       &
     &                    OCEAN(ng)%t(i,j,k,Nout,itrc)
                END DO
              END DO
            END DO
          END DO
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
        ELSE IF (iic(ng).gt.ntsAVG(ng)) THEN
          DO itrc=NAT+1,NT(ng)
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,itrc)=                         &
     &                        AVERAGE(ng)%avgt(i,j,k,itrc)+             &
     &                        OCEAN(ng)%t(i,j,k,Nout,itrc)
                END DO
              END DO
            END DO
          END DO
        END IF
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!-----------------------------------------------------------------------
!
        IF ((iic(ng).gt.ntsAVG(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
        fac=1.0_r8/REAL(nAVG(ng),r8)
          DO itrc=NAT+1,NT(ng)
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,itrc)=                       &
     &                        fac*AVERAGE(ng)%avgt(i,j,k,itrc)
                END DO
              END DO
            END DO
          END DO
        END IF
      END IF
!     IF (Master) THEN
!       print *, 'SET_FILTER  nw:',nw,'dtfil',dtfil,'dt',INT(dt(ng))
!       print *, 'SET_FILTER  nfil:',nfil,'navg',navg,'nfil/nw',nfil/nw
!       print *, 'SET_FILTER  iic-1:',iic-1
!       print *, 'SET_FILTER  ntstart:',ntstart,'nrrec',nrrec
!     END IF
      IF ((iic(ng).gt.ntsavg(ng)).and.                                  &
     &     (MOD(iic(ng)-ntstart(ng),nfil/nw).eq.0).and.                 &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
        DO f = 1,FIL
          IF (fcount(f).eq.0) THEN
            iw=1
             IF (Master) THEN
               print *, 'SET_FILTER initialize sum', f, '    iw:', iw
             END IF
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                FILTER(ng)%filzeta(i,j,f)=                              &
     &                                   OCEAN(ng)%zeta(i,j,Kout)*fw(iw)
                FILTER(ng)%filu2d (i,j,f)=                              &
     &                                   OCEAN(ng)%ubar(i,j,Kout)*fw(iw)
                FILTER(ng)%filv2d (i,j,f)=                              &
     &                                   OCEAN(ng)%vbar(i,j,Kout)*fw(iw)
              END DO
            END DO
# ifdef SOLVE3D
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filu3d(i,j,k,f)=                           &
     &                                   OCEAN(ng)%u(i,j,k,Nout)*fw(iw)
                  FILTER(ng)%filv3d(i,j,k,f)=                           &
     &                                   OCEAN(ng)%v(i,j,k,Nout)*fw(iw)
                  FILTER(ng)%filrho(i,j,k,f)=                           &
     &                                   OCEAN(ng)%rho(i,j,k)*fw(iw)
                END DO
              END DO
              DO itrc=1,NAT
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filt(i,j,k,itrc,f)=                      &
     &                         OCEAN(ng)%t(i,j,k,Nout,itrc)*fw(iw)
                  END DO
                END DO
              END DO
            END DO
            DO k=0,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filw3d(i,j,k,f)=                           &
     &                                 OCEAN(ng)%W(i,j,k)*fw(iw)*       &
     &                                 GRID(ng)%pm(i,j)*GRID(ng)%pn(i,j)
                END DO
              END DO
            END DO
#  ifdef LMD_SKPP
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                FILTER(ng)%filhsbl(i,j,f)=MIXING(ng)%hsbl(i,j)*fw(iw)
              END DO
            END DO
#  endif
#  ifdef LMD_BKPP
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                FILTER(ng)%filhbbl(i,j,f)=MIXING(ng)%hbbl(i,j)*fw(iw)
              END DO
            END DO
#  endif
#  ifdef AVERAGES_AKV
            DO k=0,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filAKv(i,j,k,f)=                           &
     &                                      MIXING(ng)%Akv(i,j,k)*fw(iw)
                END DO
              END DO
            END DO
#  endif
#  ifdef AVERAGES_AKT
            DO k=0,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filAKt(i,j,k,f)=                           &
     &                                MIXING(ng)%Akt(i,j,k,itemp)*fw(iw)
                END DO
              END DO
            END DO
#  endif
#  ifdef AVERAGES_AKS
            DO k=0,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filAKs(i,j,k,f)=                           &
     &                                MIXING(ng)%Akt(i,j,k,isalt)*fw(iw)
                END DO
              END DO
            END DO
#  endif
#  ifdef AVERAGES_FLUXES
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                FILTER(ng)%filstf(i,j,f)=                               &
     &                                FORCES(ng)%stflx(i,j,itemp)*fw(iw)
                FILTER(ng)%filswf(i,j,f)=                               &
     &                                FORCES(ng)%stflx(i,j,isalt)*fw(iw)
#   ifdef BULK_FLUXES
                FILTER(ng)%fillhf(i,j,f)= FORCES(ng)%lhflx(i,j)*fw(iw)
                FILTER(ng)%fillrf(i,j,f)= FORCES(ng)%shflx(i,j)*fw(iw)
                FILTER(ng)%filshf(i,j,f)= FORCES(ng)%lrflx(i,j)*fw(iw)
#   endif
#   ifdef SHORTWAVE
                FILTER(ng)%filsrf(i,j,f)= FORCES(ng)%srflx(i,j)*fw(iw)
#   endif
              END DO
            END DO
#  endif
# endif
# ifdef AVERAGES_FLUXES
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                FILTER(ng)%filsus(i,j,f) = FORCES(ng)%sustr(i,j)*fw(iw)
                FILTER(ng)%filsvs(i,j,f) = FORCES(ng)%svstr(i,j)*fw(iw)
              END DO
            END DO
# endif
#  ifdef ICE_MODEL
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                FILTER(ng)%filui(i,j,f) = ICE(ng)%ui(i,j,Iuout)*fw(iw)
                FILTER(ng)%filvi(i,j,f) = ICE(ng)%vi(i,j,Iuout)*fw(iw)
                FILTER(ng)%filai(i,j,f) = ICE(ng)%ai(i,j,Iout)*fw(iw)
                FILTER(ng)%filhi(i,j,f) = ICE(ng)%hi(i,j,Iout)*fw(iw)
                FILTER(ng)%filhsn(i,j,f) = ICE(ng)%hsn(i,j,Iout)*fw(iw)
                FILTER(ng)%filsfwat(i,j,f) = ICE(ng)%sfwat(i,j,Iout)*   &
     &                            fw(iw)
                FILTER(ng)%filtis(i,j,f) = ICE(ng)%tis(i,j)*fw(iw)
                FILTER(ng)%filti(i,j,f) = ICE(ng)%ti(i,j,Iout)*fw(iw)
                FILTER(ng)%filAgeice(i,j,f) = ICE(ng)%ageice(i,j,Iout)* &
     &                            fw(iw)
                FILTER(ng)%fils0mk(i,j,f) = ICE(ng)%s0mk(i,j)*fw(iw)
                FILTER(ng)%filt0mk(i,j,f) = ICE(ng)%t0mk(i,j)*fw(iw)
                FILTER(ng)%filWfr(i,j,f) = ICE(ng)%wfr(i,j)*fw(iw)
                FILTER(ng)%filWai(i,j,f) = ICE(ng)%wai(i,j)*fw(iw)
                FILTER(ng)%filWao(i,j,f) = ICE(ng)%wao(i,j)*fw(iw)
                FILTER(ng)%filWio(i,j,f) = ICE(ng)%wio(i,j)*fw(iw)
                FILTER(ng)%filWro(i,j,f) = ICE(ng)%wro(i,j)*fw(iw)
                FILTER(ng)%filsig11(i,j,f) = ICE(ng)%sig11(i,j,Ieout)*  &
     &                            fw(iw)
                FILTER(ng)%filsig12(i,j,f) = ICE(ng)%sig12(i,j,Ieout)*  &
     &                            fw(iw)
                FILTER(ng)%filsig22(i,j,f) = ICE(ng)%sig22(i,j,Ieout)*  &
     &                            fw(iw)
                FILTER(ng)%filio_mflux(i,j,f) = ICE(ng)%io_mflux(i,j)*  &
     &                            fw(iw)
                FILTER(ng)%filutau_iw(i,j,f) = ICE(ng)%utau_iw(i,j)*    &
     &                            fw(iw)
                FILTER(ng)%filchu_iw(i,j,f) = ICE(ng)%chu_iw(i,j)*fw(iw)
              END DO
            END DO
#  endif
            fcount(f)=iw
!
!-----------------------------------------------------------------------
!  Accumulate time-filtered fields.
!-----------------------------------------------------------------------
!
          ELSE
            iw=fcount(f)+1
             IF (Master) THEN
               print *, 'SET_FILTER add to sum', f, '    iw:', iw
             END IF
            IF (fcount(f).gt.0) THEN
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filzeta(i,j,f)=FILTER(ng)%filzeta(i,j,f)+  &
     &                               OCEAN(ng)%zeta(i,j,Kout)*fw(iw)
                  FILTER(ng)%filu2d (i,j,f)=FILTER(ng)%filu2d (i,j,f)+  &
     &                               OCEAN(ng)%ubar(i,j,Kout)*fw(iw)
                  FILTER(ng)%filv2d (i,j,f)=FILTER(ng)%filv2d (i,j,f)+  &
     &                               OCEAN(ng)%vbar(i,j,Kout)*fw(iw)
                END DO
              END DO
# ifdef SOLVE3D
              DO k=1,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filu3d(i,j,k,f)=                         &
     &                                    FILTER(ng)%filu3d(i,j,k,f)+   &
     &                                    OCEAN(ng)%u(i,j,k,Nout)*fw(iw)
                    FILTER(ng)%filv3d(i,j,k,f)=                         &
     &                                    FILTER(ng)%filv3d(i,j,k,f)+   &
     &                                    OCEAN(ng)%v(i,j,k,Nout)*fw(iw)
                    FILTER(ng)%filrho(i,j,k,f)=                         &
     &                                    FILTER(ng)%filrho(i,j,k,f)+   &
     &                                       OCEAN(ng)%rho(i,j,k)*fw(iw)
                  END DO
                END DO
                DO itrc=1,NAT
                  DO j=JstrR,JendR
                    DO i=IstrR,IendR
                      FILTER(ng)%filt(i,j,k,itrc,f)=                    &
     &                           FILTER(ng)%filt(i,j,k,itrc,f)+   &
     &                           OCEAN(ng)%t(i,j,k,Nout,itrc)*fw(iw)
                    END DO
                  END DO
                END DO
              END DO
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filw3d(i,j,k,f)=                         &
     &                                  FILTER(ng)%filw3d(i,j,k,f)+     &
     &                                  OCEAN(ng)%W(i,j,k)*fw(iw)*      &
     &                                  GRID(ng)%pm(i,j)*               &
     &                                  GRID(ng)%pn(i,j)
                  END DO
                END DO
              END DO
#  ifdef LMD_SKPP
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filhsbl(i,j,f)=FILTER(ng)%filhsbl(i,j,f)+  &
     &                               MIXING(ng)%hsbl(i,j)*fw(iw)
                END DO
              END DO
#  endif
#  ifdef LMD_BKPP
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filhbbl(i,j,f)=FILTER(ng)%filhbbl(i,j,f)+  &
     &                               MIXING(ng)%hbbl(i,j)*fw(iw)
                END DO
              END DO
#  endif
#  ifdef AVERAGES_AKV
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filAKv(i,j,k,f)=                         &
     &                                      FILTER(ng)%filAKv(i,j,k,f)+ &
     &                                      MIXING(ng)%Akv(i,j,k)*fw(iw)
                  END DO
                END DO
              END DO
#  endif
#  ifdef AVERAGES_AKT
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filAKt(i,j,k,f)=                         &
     &                                FILTER(ng)%filAKt(i,j,k,f)+       &
     &                                MIXING(ng)%Akt(i,j,k,itemp)*fw(iw)
                  END DO
                END DO
              END DO
#  endif
#  ifdef AVERAGES_AKS
              DO k=0,N(ng)
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    FILTER(ng)%filAKs(i,j,k,f)=                         &
     &                                FILTER(ng)%filAKs(i,j,k,f)+       &
     &                                MIXING(ng)%Akt(i,j,k,isalt)*fw(iw)
                  END DO
                END DO
              END DO
#  endif
#  ifdef AVERAGES_FLUXES
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filstf(i,j,f)=FILTER(ng)%filstf(i,j,f)+    &
     &                                FORCES(ng)%stflx(i,j,itemp)*fw(iw)
                  FILTER(ng)%filswf(i,j,f)=FILTER(ng)%filswf(i,j,f)+    &
     &                                FORCES(ng)%stflx(i,j,isalt)*fw(iw)
#   ifdef BULK_FLUXES
                  FILTER(ng)%fillhf(i,j,f)=FILTER(ng)%fillhf(i,j,f)+    &
     &                                      FORCES(ng)%lhflx(i,j)*fw(iw)
                  FILTER(ng)%fillrf(i,j,f)=FILTER(ng)%fillrf(i,j,f)+    &
     &                                      FORCES(ng)%shflx(i,j)*fw(iw)
                  FILTER(ng)%filshf(i,j,f)=FILTER(ng)%filshf(i,j,f)+    &
     &                                      FORCES(ng)%lrflx(i,j)*fw(iw)
#   endif
#   ifdef SHORTWAVE
                  FILTER(ng)%filsrf(i,j,f)=FILTER(ng)%filsrf(i,j,f)+    &
     &                                      FORCES(ng)%srflx(i,j)*fw(iw)
#   endif
                END DO
              END DO
#  endif
# endif
# ifdef AVERAGES_FLUXES
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filsus(i,j,f)=FILTER(ng)%filsus(i,j,f)+    &
     &                                      FORCES(ng)%sustr(i,j)*fw(iw)
                  FILTER(ng)%filsvs(i,j,f)=FILTER(ng)%filsvs(i,j,f)+    &
     &                                      FORCES(ng)%svstr(i,j)*fw(iw)
                END DO
              END DO
# endif
# ifdef ICE_MODEL
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  FILTER(ng)%filui(i,j,f)=FILTER(ng)%filui(i,j,f)+      &
     &                                    ICE(ng)%ui(i,j,Iuout)*fw(iw)
                  FILTER(ng)%filvi(i,j,f)=FILTER(ng)%filvi(i,j,f)+      &
     &                                    ICE(ng)%vi(i,j,Iuout)*fw(iw)
                  FILTER(ng)%filai(i,j,f)=FILTER(ng)%filai(i,j,f)+      &
     &                                    ICE(ng)%ai(i,j,Iout)*fw(iw)
                  FILTER(ng)%filhi(i,j,f)=FILTER(ng)%filhi(i,j,f)+      &
     &                                    ICE(ng)%hi(i,j,Iout)*fw(iw)
                  FILTER(ng)%filhsn(i,j,f)=FILTER(ng)%filhsn(i,j,f)+    &
     &                                    ICE(ng)%hsn(i,j,Iout)*fw(iw)
                  FILTER(ng)%filsfwat(i,j,f)=FILTER(ng)%filsfwat(i,j,f)+&
     &                                    ICE(ng)%sfwat(i,j,Iout)*fw(iw)
                  FILTER(ng)%filtis(i,j,f)=FILTER(ng)%filtis(i,j,f)+    &
     &                                    ICE(ng)%tis(i,j)*fw(iw)
                  FILTER(ng)%filti(i,j,f)=FILTER(ng)%filAgeice(i,j,f)+  &
     &                                ICE(ng)%ageice(i,j,Iout)*fw(iw)
                  FILTER(ng)%filti(i,j,f)=FILTER(ng)%filti(i,j,f)+      &
     &                                    ICE(ng)%ti(i,j,Iout)*fw(iw)
                  FILTER(ng)%fils0mk(i,j,f)=FILTER(ng)%fils0mk(i,j,f)+  &
     &                                    ICE(ng)%s0mk(i,j)*fw(iw)
                  FILTER(ng)%filt0mk(i,j,f)=FILTER(ng)%filt0mk(i,j,f)+  &
     &                                    ICE(ng)%t0mk(i,j)*fw(iw)
                  FILTER(ng)%filWfr(i,j,f)=FILTER(ng)%filWfr(i,j,f)+  &
     &                                    ICE(ng)%wfr(i,j)*fw(iw)
                  FILTER(ng)%filWai(i,j,f)=FILTER(ng)%filWai(i,j,f)+  &
     &                                    ICE(ng)%wai(i,j)*fw(iw)
                  FILTER(ng)%filWao(i,j,f)=FILTER(ng)%filWao(i,j,f)+  &
     &                                    ICE(ng)%wao(i,j)*fw(iw)
                  FILTER(ng)%filWio(i,j,f)=FILTER(ng)%filWio(i,j,f)+  &
     &                                    ICE(ng)%wio(i,j)*fw(iw)
                  FILTER(ng)%filWro(i,j,f)=FILTER(ng)%filWro(i,j,f)+  &
     &                                    ICE(ng)%wro(i,j)*fw(iw)
                  FILTER(ng)%filsig11(i,j,f)=FILTER(ng)%filsig11(i,j,f)+&
     &                               ICE(ng)%sig11(i,j,Ieout)*fw(iw)
                  FILTER(ng)%filsig12(i,j,f)=FILTER(ng)%filsig12(i,j,f)+&
     &                               ICE(ng)%sig12(i,j,Ieout)*fw(iw)
                  FILTER(ng)%filsig22(i,j,f)=FILTER(ng)%filsig22(i,j,f)+&
     &                               ICE(ng)%sig22(i,j,Ieout)*fw(iw)
                  FILTER(ng)%filio_mflux(i,j,f)=                        &
     &                                FILTER(ng)%filio_mflux(i,j,f)+    &
     &                                    ICE(ng)%io_mflux(i,j)*fw(iw)
                  FILTER(ng)%filutau_iw(i,j,f)=                         &
     &                                FILTER(ng)%filutau_iw(i,j,f)+     &
     &                                    ICE(ng)%utau_iw(i,j)*fw(iw)
                  FILTER(ng)%filchu_iw(i,j,f)=                          &
     &                                FILTER(ng)%filchu_iw(i,j,f)+      &
     &                                    ICE(ng)%chu_iw(i,j)*fw(iw)
                END DO
              END DO
# endif
            END IF
            fcount(f)=iw
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Transfer finished sums into averages arrays for writing out.
!-----------------------------------------------------------------------
!
      IF ((MOD(iic(ng)-1,navg(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN
        DO f=1,FIL
          IF (fcount(f).eq.nw) THEN
             IF (Master) THEN
               print *,'SET_FILTER transferring sums', f
             END IF
            IF (SOUTH_WEST_TEST) THEN
              avgtime(ng)=time(ng)-0.5_r8*REAL(nfil,r8)*dt(ng)
            END IF
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgzeta(i,j)=FILTER(ng)%filzeta(i,j,f)
                AVERAGE(ng)%avgu2d (i,j)=FILTER(ng)%filu2d (i,j,f)
                AVERAGE(ng)%avgv2d (i,j)=FILTER(ng)%filv2d (i,j,f)
              END DO
            END DO
# ifdef SOLVE3D
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgu3d(i,j,k)=FILTER(ng)%filu3d(i,j,k,f)
                  AVERAGE(ng)%avgv3d(i,j,k)=FILTER(ng)%filv3d(i,j,k,f)
                  AVERAGE(ng)%avgrho(i,j,k)=FILTER(ng)%filrho(i,j,k,f)
                END DO
              END DO
              DO itrc=1,NAT
                DO j=JstrR,JendR
                  DO i=IstrR,IendR
                    AVERAGE(ng)%avgt(i,j,k,itrc)=                       &
     &                                FILTER(ng)%filt(i,j,k,itrc,f)
                  END DO
                END DO
              END DO
            END DO
            DO k=0,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgw3d(i,j,k)=FILTER(ng)%filw3d(i,j,k,f)
                END DO
              END DO
            END DO
#  ifdef LMD_SKPP
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avghsbl(i,j)=FILTER(ng)%filhsbl(i,j,f)
              END DO
            END DO
#  endif
#  ifdef LMD_BKPP
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avghbbl(i,j)=FILTER(ng)%filhbbl(i,j,f)
              END DO
            END DO
#  endif
#  ifdef AVERAGES_AKV
            DO k=0,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgAKv(i,j,k)=FILTER(ng)%filAKv(i,j,k,f)
                END DO
              END DO
            END DO
#  endif
#  ifdef AVERAGES_AKT
            DO k=0,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgAKt(i,j,k)=FILTER(ng)%filAKt(i,j,k,f)
                END DO
              END DO
            END DO
#  endif
#  ifdef AVERAGES_AKS
            DO k=0,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgAKs(i,j,k)=FILTER(ng)%filAKs(i,j,k,f)
                END DO
              END DO
            END DO
#  endif
#  ifdef AVERAGES_FLUXES
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgstf(i,j)=FILTER(ng)%filstf(i,j,f)
                AVERAGE(ng)%avgswf(i,j)=FILTER(ng)%filswf(i,j,f)
#   ifdef BULK_FLUXES
                AVERAGE(ng)%avglhf(i,j)=FILTER(ng)%fillhf(i,j,f)
                AVERAGE(ng)%avglrf(i,j)=FILTER(ng)%fillrf(i,j,f)
                AVERAGE(ng)%avgshf(i,j)=FILTER(ng)%filshf(i,j,f)
#   endif
#   ifdef SHORTWAVE
                AVERAGE(ng)%avgsrf(i,j)=FILTER(ng)%filsrf(i,j,f)
#   endif
              END DO
            END DO
#  endif
# endif
# ifdef AVERAGES_FLUXES
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgsus(i,j)=FILTER(ng)%filsus(i,j,f)
                AVERAGE(ng)%avgsvs(i,j)=FILTER(ng)%filsvs(i,j,f)
              END DO
            END DO
# endif
# ifdef ICE_MODEL
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avguice(i,j)=FILTER(ng)%filui(i,j,f)
                AVERAGE(ng)%avgvice(i,j)=FILTER(ng)%filvi(i,j,f)
                AVERAGE(ng)%avgaice(i,j)=FILTER(ng)%filai(i,j,f)
                AVERAGE(ng)%avghice(i,j)=FILTER(ng)%filhi(i,j,f)
                AVERAGE(ng)%avghsno(i,j)=FILTER(ng)%filhsn(i,j,f)
                AVERAGE(ng)%avgsfwat(i,j)=FILTER(ng)%filsfwat(i,j,f)
                AVERAGE(ng)%avgtice(i,j)=FILTER(ng)%filtis(i,j,f)
                AVERAGE(ng)%avgtimid(i,j)=FILTER(ng)%filti(i,j,f)
                AVERAGE(ng)%avgAgeice(i,j)=FILTER(ng)%filAgeice(i,j,f)
                AVERAGE(ng)%avgS0mk(i,j)=FILTER(ng)%fils0mk(i,j,f)
                AVERAGE(ng)%avgT0mk(i,j)=FILTER(ng)%filt0mk(i,j,f)
                AVERAGE(ng)%avgWfr(i,j)=FILTER(ng)%filWfr(i,j,f)
                AVERAGE(ng)%avgWai(i,j)=FILTER(ng)%filWai(i,j,f)
                AVERAGE(ng)%avgWao(i,j)=FILTER(ng)%filWao(i,j,f)
                AVERAGE(ng)%avgWio(i,j)=FILTER(ng)%filWio(i,j,f)
                AVERAGE(ng)%avgWro(i,j)=FILTER(ng)%filWro(i,j,f)
                AVERAGE(ng)%avgsig11(i,j)=FILTER(ng)%filsig11(i,j,f)
                AVERAGE(ng)%avgsig12(i,j)=FILTER(ng)%filsig12(i,j,f)
                AVERAGE(ng)%avgsig22(i,j)=FILTER(ng)%filsig22(i,j,f)
                AVERAGE(ng)%avgiomflx(i,j)=                             &
     &                          FILTER(ng)%filio_mflux(i,j,f)
                AVERAGE(ng)%avgutau_iw(i,j)=FILTER(ng)%filutau_iw(i,j,f)
                AVERAGE(ng)%avgchu_iw(i,j)=FILTER(ng)%filchu_iw(i,j,f)
              END DO
            END DO
# endif
            fcount(f)=0
          END IF
        END DO
      END IF
!     IF (Master) print *, 'end of SET_FILTER'
      RETURN
      END SUBROUTINE set_filter_tile
#else
      SUBROUTINE set_filter
      RETURN
      END SUBROUTINE set_filter
#endif
