#include "cppdefs.h"
      MODULE set_vbc_mod
#ifdef NONLINEAR
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module sets vertical boundary conditons for momentum and       !
!  tracers.                                                            !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: set_vbc

      CONTAINS

# ifdef SOLVE3D
!
!***********************************************************************
      SUBROUTINE set_vbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iNLM, 6)
#  endif
      CALL set_vbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                      &
     &                   GRID(ng) % Hz,                                 &
#  if !defined BBL_MODEL || defined ICESHELF
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % z_w,                                &
#  endif
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
#  if defined ICESHELF
     &                   GRID(ng) % zice,                               &
#  endif
     &                   OCEAN(ng) % t,                                 &
#  if !defined BBL_MODEL || defined ICESHELF
     &                   OCEAN(ng) % u,                                 &
     &                   OCEAN(ng) % v,                                 &
#  endif
#  ifdef QCORRECTION
     &                   FORCES(ng) % dqdt,                             &
     &                   FORCES(ng) % sst,                              &
#  endif
#  ifdef RDRG_GRID
     &                   FORCES(ng) % rdrg_grid,                        &
#  endif
#  if defined SCORRECTION || defined SRELAXATION
     &                   FORCES(ng) % sss,                              &
#  endif
#  if defined SSSFLX || defined WRT_SSSFLX
     &                   FORCES(ng) % sssflx,                           &
#  endif
#  if defined ICESHELF
#   ifdef SHORTWAVE
     &                   FORCES(ng) % srflx,                            &
#   endif
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
#  endif
#  ifndef BBL_MODEL
     &                   FORCES(ng) % bustr,                            &
     &                   FORCES(ng) % bvstr,                            &
#  endif
     &                   FORCES(ng) % stflx,                            &
     &                   FORCES(ng) % btflx)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 6)
#  endif
      RETURN
      END SUBROUTINE set_vbc
!
!***********************************************************************
      SUBROUTINE set_vbc_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nrhs,                                    &
     &                         Hz,                                      &
#  if !defined BBL_MODEL || defined ICESHELF
     &                         z_r, z_w,                                &
#  endif
#  if defined MASKING
     &                         rmask,                                   &
#  endif
#  if defined ICESHELF
     &                         zice,                                    &
#  endif
     &                         t,                                       &
#  if !defined BBL_MODEL || defined ICESHELF
     &                         u, v,                                    &
#  endif
#  ifdef QCORRECTION
     &                         dqdt, sst,                               &
#  endif
#  ifdef RDRG_GRID
     &                         rdrg_grid,                               &
#  endif
#  if defined SCORRECTION || defined SRELAXATION
     &                         sss,                                     &
#  endif
#  if defined SSSFLX || defined WRT_SSSFLX
     &                         sssflx,                                  &
#  endif
#  if defined ICESHELF
#   ifdef SHORTWAVE
     &                         srflx,                                   &
#   endif
     &                         sustr, svstr,                            &
#  endif
#  ifndef BBL_MODEL
     &                         bustr, bvstr,                            &
#  endif
     &                         stflx, btflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
#   if !defined BBL_MODEL || defined ICESHELF
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
#   endif
#   if defined MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#   endif
#   if defined ICESHELF
      real(r8), intent(in) :: zice(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#   if !defined BBL_MODEL || defined ICESHELF
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
#   endif
#   ifdef QCORRECTION
      real(r8), intent(in) :: dqdt(LBi:,LBj:)
      real(r8), intent(in) :: sst(LBi:,LBj:)
#   endif
#   ifdef RDRG_GRID
      real(r8), intent(in) :: rdrg_grid(LBi:,LBj:)
#   endif
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: sss(LBi:,LBj:)
#   endif
#   if defined SSSFLX || defined WRT_SSSFLX
      real(r8), intent(inout) :: sssflx(LBi:,LBj:)
#   endif
#   if defined ICESHELF
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:,LBj:)
#    endif
      real(r8), intent(inout) :: sustr(LBi:,LBj:)
      real(r8), intent(inout) :: svstr(LBi:,LBj:)
#   endif
#   ifndef BBL_MODEL
      real(r8), intent(inout) :: bustr(LBi:,LBj:)
      real(r8), intent(inout) :: bvstr(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
      real(r8), intent(inout) :: btflx(LBi:,LBj:,:)
#  else
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
#   if !defined BBL_MODEL || defined ICESHELF
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#   endif
#   if defined MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#   endif
#   if defined ICESHELF
      real(r8), intent(in) :: zice(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
#   if !defined BBL_MODEL || defined ICESHELF
      real(r8), intent(in) :: u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(in) :: v(LBi:UBi,LBj:UBj,N(ng),2)
#   endif
#   ifdef QCORRECTION
      real(r8), intent(in) :: dqdt(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: sst(LBi:UBi,LBj:UBj)
#   endif
#   ifdef RDRG_GRID
      real(r8), intent(in) :: rdrg_grid(LBi:UBi,LBj:UBj)
#   endif
#   if defined SCORRECTION || defined SRELAXATION
      real(r8), intent(in) :: sss(LBi:UBi,LBj:UBj)
#   endif
#   if defined SSSFLX || defined WRT_SSSFLX
      real(r8), intent(in) :: sssflx(LBi:UBi,LBj:UBj)
#   endif
#   if defined ICESHELF
#    ifdef SHORTWAVE
      real(r8), intent(inout) :: srflx(LBi:UBi,LBj:UBj)
#    endif
      real(r8), intent(inout) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: svstr(LBi:UBi,LBj:UBj)
#   endif
#   ifndef BBL_MODEL
      real(r8), intent(inout) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bvstr(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: stflx(LBi:UBi,LBj:UBj,NT(ng))
      real(r8), intent(inout) :: btflx(LBi:UBi,LBj:UBj,NT(ng))
#  endif

!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, j, itrc

# ifdef ICESHELF
      real(r8), parameter :: gamma = 0.0001_r8
      real(r8), parameter :: hlfreeze = 334000.0_r8
      real(r8), parameter :: refSalt = 34.4_r8
      real(r8), parameter :: trelax = 30.0_r8 * 86400.0_r8
      real(r8), parameter :: sfcTemp = -1.9_r8
      real(r8), parameter :: sfcSalt = 34.4_r8

      real(r8) :: temp_f
# endif

#  if !defined BBL_MODEL || defined ICESHELF
      real(r8) :: cff1, cff2, cff3
#  endif

#  if (!defined BBL_MODEL || defined ICESHELF) && defined UV_LOGDRAG
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
#  endif

#  include "set_bounds.h"

#  ifdef QCORRECTION
!
!-----------------------------------------------------------------------
!  Add in flux correction to surface net heat flux (degC m/s).
!-----------------------------------------------------------------------
!
! Add in net heat flux correction.
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          stflx(i,j,itemp)=stflx(i,j,itemp)+                            &
     &                     dqdt(i,j)*(t(i,j,N(ng),nrhs,itemp)-sst(i,j))
        END DO
      END DO
#  endif
#  ifdef SALINITY
!
!-----------------------------------------------------------------------
!  Multiply fresh water flux with surface salinity. If appropriate,
!  apply correction.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
#   if defined SCORRECTION
#     ifdef ICE_THERMO
          stflx(i,j,isalt)=stflx(i,j,isalt) -                           &
     &                     Tnudg(isalt,ng)*Hz(i,j,N(ng))*               &
     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
#     ifdef MASKING
          stflx(i,j,isalt) = rmask(i,j)*stflx(i,j,isalt)
#     endif
#    else
          stflx(i,j,isalt)=stflx(i,j,isalt)*t(i,j,N(ng),nrhs,isalt)-    &
     &                     Tnudg(isalt,ng)*Hz(i,j,N(ng))*               &
     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
#     endif
#   elif defined SRELAXATION
          stflx(i,j,isalt)=-Tnudg(isalt,ng)*Hz(i,j,N(ng))*              &
     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
#     ifdef MASKING
          stflx(i,j,isalt) = rmask(i,j)*stflx(i,j,isalt)
#     endif
#   else
#     ifdef ICE_THERMO
          stflx(i,j,isalt)=stflx(i,j,isalt)
#     else
          stflx(i,j,isalt)=stflx(i,j,isalt)*t(i,j,N(ng),nrhs,isalt)
#     endif
#     ifdef MASKING
          stflx(i,j,isalt) = rmask(i,j)*stflx(i,j,isalt)
#     endif
#   endif
#   ifdef SSSFLX
          stflx(i,j,isalt)=stflx(i,j,isalt)+sssflx(i,j)
#     ifdef MASKING
          stflx(i,j,isalt) = rmask(i,j)*stflx(i,j,isalt)
#     endif

#   endif
#   if (defined SCORRECTION || defined SRELAXATION) && \
        defined WRT_SSSFLX
          sssflx(i,j) = -Tnudg(isalt,ng)*Hz(i,j,N(ng))*                 &
     &                     (t(i,j,N(ng),nrhs,isalt)-sss(i,j))
#     ifdef MASKING
          sssflx(i,j) = rmask(i,j)*sssflx(i,j)
#     endif
#   endif
          btflx(i,j,isalt)=btflx(i,j,isalt)*t(i,j,1,nrhs,isalt)
        END DO
      END DO
#  endif
#  ifdef ICESHELF
!
!-----------------------------------------------------------------------
!  If ice shelf cavities, zero out for now the surface tracer flux
!  over the ice.
!  Modification for the ICETEST (4/24/06):  Now compute the surface
!  tracer fluxes underneath ice in the manner defined by the test
!  specification.  Note that the units for stflx for heat are degC m/s 
!  (which means the rho*C_p part is removed here) and for salt are 
!  psu m/s.
!  Further modification for the ICETEST 2.01 case (4/24/06):  Compute
!  heat and salt fluxes for the open water/atmosphere interface to match
!  the relaxation calculation in sections 3.2.6 and 3.2.7 of John 
!  Hunter's document.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (zice(i,j).ne.0.0_r8) THEN
            temp_f=0.0939_r8-0.057_r8*t(i,j,N(ng),nrhs,isalt)+          &
     &             7.6410E-4_r8*zice(i,j)
            stflx(i,j,itemp)=gamma*(temp_f-t(i,j,N(ng),nrhs,itemp))
            stflx(i,j,isalt)=Cp*stflx(i,j,itemp)*refSalt/hlfreeze
          ELSE
            stflx(i,j,itemp)=Hz(i,j,N(ng))*                             &
     &                       (sfcTemp-t(i,j,N(ng),nrhs,itemp))/trelax
            stflx(i,j,isalt)=Hz(i,j,N(ng))*                             &
     &                       (sfcSalt-t(i,j,N(ng),nrhs,isalt))/trelax
          END IF
        END DO
      END DO
#   ifdef SHORTWAVE
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (zice(i,j).ne.0.0_r8) THEN
            srflx(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif
!
!-----------------------------------------------------------------------
!  If ice shelf cavities, replace surface wind stress with ice shelf
!  cavity stress (m2/s2).
!-----------------------------------------------------------------------

#   if defined UV_LOGDRAG
!
!  Set logarithmic ice shelf cavity stress.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=1.0_r8/LOG((z_w(i,j,N(ng))-z_r(i,j,N(ng)))/Zob(ng))
          cff2=vonKar*vonKar*cff1*cff1
          wrk(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff2))
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            cff1=0.25_r8*(v(i  ,j  ,N(ng),nrhs)+                        &
     &                    v(i  ,j+1,N(ng),nrhs)+                        &
     &                    v(i-1,j  ,N(ng),nrhs)+                        &
     &                    v(i-1,j+1,N(ng),nrhs))
            cff2=SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+cff1*cff1)
            sustr(i,j)=-0.5_r8*(wrk(i-1,j)+wrk(i,j))*                   &
     &                 u(i,j,N(ng),nrhs)*cff2
          END IF
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            cff1=0.25_r8*(u(i  ,j  ,N(ng),nrhs)+                        &
     &                    u(i+1,j  ,N(ng),nrhs)+                        &
     &                    u(i  ,j-1,N(ng),nrhs)+                        &
     &                    u(i+1,j-1,N(ng),nrhs))
            cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
            svstr(i,j)=-0.5_r8*(wrk(i,j-1)+wrk(i,j))*                   &
     &                 v(i,j,N(ng),nrhs)*cff2
          END IF
        END DO
      END DO
#   elif defined UV_QDRAG
!
!  Set quadratic ice shelf cavity stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            cff1=0.25_r8*(v(i  ,j  ,N(ng),nrhs)+                        &
     &                    v(i  ,j+1,N(ng),nrhs)+                        &
     &                    v(i-1,j  ,N(ng),nrhs)+                        &
     &                    v(i-1,j+1,N(ng),nrhs))
            cff2=SQRT(u(i,j,N(ng),nrhs)*u(i,j,N(ng),nrhs)+cff2*cff2)
#ifdef RDRG_GRID
            sustr(i,j)=-0.5_r8*(rdrg_grid(i,j)+rdrg_grid(i-1,j))*       &
     &                  u(i,j,N(ng),nrhs)*cff2
#else
            sustr(i,j)=-rdrg2(ng)*u(i,j,N(ng),nrhs)*cff2
#endif
          END IF
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            cff1=0.25_r8*(u(i  ,j  ,N(ng),nrhs)+                        &
     &                    u(i+1,j  ,N(ng),nrhs)+                        &
     &                    u(i  ,j-1,N(ng),nrhs)+                        &
     &                    u(i+1,j-1,N(ng),nrhs))
            cff2=SQRT(cff1*cff1+v(i,j,N(ng),nrhs)*v(i,j,N(ng),nrhs))
#ifdef RDRG_GRID
            svstr(i,j)=-0.5_r8*(rdrg_grid(i,j)+rdrg_grid(i,j-1))*       &
     &                  v(i,j,N(ng),nrhs)*cff2
#else
            svstr(i,j)=-rdrg2(ng)*v(i,j,N(ng),nrhs)*cff2
#endif
          END IF
        END DO
      END DO
#   elif defined UV_LDRAG
!
!  Set linear ice shelf cavity stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
#ifdef RDRG_GRID
            sustr(i,j)=-0.5_r8*(rdrg_grid(i,j)+rdrg_grid(i-1,j))*       &
     &                  u(i,j,N(ng),nrhs)
#else
            sustr(i,j)=-rdrg(ng)*u(i,j,N(ng),nrhs)
#endif
          END IF
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
#ifdef RDRG_GRID
            svstr(i,j)=-0.5_r8*(rdrg_grid(i,j)+rdrg_grid(i,j-1))*       &
     &                  v(i,j,N(ng),nrhs)
#else
            svstr(i,j)=-rdrg(ng)*v(i,j,N(ng),nrhs)
#endif
          END IF
        END DO
      END DO
#   else
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          IF (zice(i,j)*zice(i-1,j).ne.0.0_r8) THEN
            sustr(i,j)=0.0_r8
          END IF
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          IF (zice(i,j)*zice(i,j-1).ne.0.0_r8) THEN
            svstr(i,j)=0.0_r8
          END IF
        END DO
      END DO
#   endif
!
!  Apply periodic or gradient boundary conditions for output
!  purposes only.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  sustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  svstr)
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    sustr, svstr)
#   endif
#  endif
#  ifndef BBL_MODEL
!
!-----------------------------------------------------------------------
!  Set kinematic bottom momentum flux (m2/s2).
!-----------------------------------------------------------------------

#   if defined UV_LOGDRAG
!
!  Set logarithmic bottom stress.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          cff1=1.0_r8/LOG((z_r(i,j,1)-z_w(i,j,0))/Zob(ng))
          cff2=vonKar*vonKar*cff1*cff1
          wrk(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff2))
        END DO
      END DO
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(v(i  ,j  ,1,nrhs)+                              &
     &                  v(i  ,j+1,1,nrhs)+                              &
     &                  v(i-1,j  ,1,nrhs)+                              &
     &                  v(i-1,j+1,1,nrhs))
          cff2=SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff1*cff1)
          bustr(i,j)=0.5_r8*(wrk(i-1,j)+wrk(i,j))*                      &
     &               u(i,j,1,nrhs)*cff2
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(u(i  ,j  ,1,nrhs)+                              &
     &                  u(i+1,j  ,1,nrhs)+                              &
     &                  u(i  ,j-1,1,nrhs)+                              &
     &                  u(i+1,j-1,1,nrhs))
          cff2=SQRT(cff1*cff1+v(i,j,1,nrhs)*v(i,j,1,nrhs))
          bvstr(i,j)=0.5_r8*(wrk(i,j-1)+wrk(i,j))*                      &
     &               v(i,j,1,nrhs)*cff2
        END DO
      END DO
#   elif defined UV_QDRAG
!
!  Set quadratic bottom stress.
!
#ifdef DRAG_LIMITER
!
! Here, 0.6 is a limiting nondimensional number for how
! much of the bottom velocity can be taken out in one timestep
!
      cff3 = 0.6
#endif
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(v(i  ,j  ,1,nrhs)+                              &
     &                  v(i  ,j+1,1,nrhs)+                              &
     &                  v(i-1,j  ,1,nrhs)+                              &
     &                  v(i-1,j+1,1,nrhs))
          cff2=SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff1*cff1)
#ifdef RDRG_GRID
          cff1=0.5_r8*(rdrg_grid(i,j)+rdrg_grid(i-1,j))
#else
          cff1=rdrg2(ng)
#endif
#ifdef DRAG_LIMITER
          cff1 = min((cff3*0.5*(Hz(i,j,1)+Hz(i-1,j,1))/(cff2*dt(ng))),  &
     &               cff1)
#endif
          bustr(i,j)=cff1*u(i,j,1,nrhs)*cff2
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(u(i  ,j  ,1,nrhs)+                              &
     &                  u(i+1,j  ,1,nrhs)+                              &
     &                  u(i  ,j-1,1,nrhs)+                              &
     &                  u(i+1,j-1,1,nrhs))
          cff2=SQRT(cff1*cff1+v(i,j,1,nrhs)*v(i,j,1,nrhs))
#ifdef RDRG_GRID
          cff1=0.5_r8*(rdrg_grid(i,j)+rdrg_grid(i,j-1))
#else
          cff1=rdrg2(ng)
#endif
#ifdef DRAG_LIMITER
          cff1 = min((cff3*0.5*(Hz(i,j,1)+Hz(i,j-1,1))/(cff2*dt(ng))),  &
     &               cff1)
#endif
          bvstr(i,j)=cff1*v(i,j,1,nrhs)*cff2
        END DO
      END DO
#    elif defined UV_LDRAG
!
!  Set linear bottom stress.
!
#ifdef DRAG_LIMITER
!
! Here, 0.6 is a limiting nondimensional number for how
! much of the bottom velocity can be taken out in one timestep
!
      cff2 = 0.6
#endif
      DO j=Jstr,Jend
        DO i=IstrU,Iend
#ifdef RDRG_GRID
          cff1=0.5_r8*(rdrg_grid(i,j)+rdrg_grid(i-1,j))
#else
          cff1=rdrg(ng)
#endif
#ifdef DRAG_LIMITER
          cff1 = min((cff2*0.5*(Hz(i,j,1)+Hz(i-1,j,1))/dt(ng)),cff1)
#endif
          bustr(i,j)=cff1*u(i,j,1,nrhs)
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
#ifdef RDRG_GRID
          cff1=0.5_r8*(rdrg_grid(i,j)+rdrg_grid(i,j-1))
#else
          cff1=rdrg(ng)
#endif
#ifdef DRAG_LIMITER
          cff1 = min((cff2*0.5*(Hz(i,j,1)+Hz(i,j-1,1))/dt(ng)),cff1)
#endif
          bvstr(i,j)=cff1*v(i,j,1,nrhs)
        END DO
      END DO
#    endif
!
!  Apply boundary conditions.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstr)
#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bustr, bvstr)
#   endif
#  endif
      RETURN
      END SUBROUTINE set_vbc_tile

# else

!
!***********************************************************************
      SUBROUTINE set_vbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iNLM, 6)
#  endif
      CALL set_vbc_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs(ng), kstp(ng), knew(ng),                  &
     &                   GRID(ng) % om_u,                               &
     &                   GRID(ng) % om_v,                               &
     &                   GRID(ng) % on_u,                               &
     &                   GRID(ng) % on_v,                               &
     &                   OCEAN(ng) % ubar,                              &
     &                   OCEAN(ng) % vbar,                              &
     &                   FORCES(ng) % bustr,                            &
     &                   FORCES(ng) % bvstr)
#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 6)
#  endif
      RETURN
      END SUBROUTINE set_vbc
!
!***********************************************************************
      SUBROUTINE set_vbc_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         krhs, kstp, knew,                        &
     &                         om_u, om_v, on_u, on_v,                  &
     &                         ubar, vbar, bustr, bvstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: krhs, kstp, knew
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: ubar(LBi:,LBj:,:)
      real(r8), intent(in) :: vbar(LBi:,LBj:,:)
      real(r8), intent(inout) :: bustr(LBi:,LBj:)
      real(r8), intent(inout) :: bvstr(LBi:,LBj:)
#  else
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(in) :: vbar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: bustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: bvstr(LBi:UBi,LBj:UBj)
#  endif
!
!  Local variable declarations.
!
# ifdef DISTRIBUTE
#  ifdef EW_PERIODIC
      logical :: EWperiodic=.TRUE.
#  else
      logical :: EWperiodic=.FALSE.
#  endif
#  ifdef NS_PERIODIC
      logical :: NSperiodic=.TRUE.
#  else
      logical :: NSperiodic=.FALSE.
#  endif
# endif
      integer :: i, j

      real(r8) :: cff1, cff2

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set kinematic barotropic bottom momentum stress (m2/s2).
!-----------------------------------------------------------------------

#  if defined UV_LDRAG
!
!  Set linear bottom stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
#ifdef RDRG_GRID
          bustr(i,j)=0.5_r8*(rdrg_grid(i,j)+rdrg_grid(i-1,j))*          &
     &                  ubar(i,j,krhs)
#else
          bustr(i,j)=rdrg(ng)*ubar(i,j,krhs)
#endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
#ifdef RDRG_GRID
          bvstr(i,j)=0.5_r8*(rdrg_grid(i,j)+rdrg_grid(i,j-1))*          &
     &                  vbar(i,j,krhs)
#else
          bvstr(i,j)=rdrg(ng)*vbar(i,j,krhs)
#endif
        END DO
      END DO
#  elif defined UV_QDRAG
!
!  Set quadratic bottom stress.
!
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=0.25_r8*(vbar(i  ,j  ,krhs)+                             &
     &                  vbar(i  ,j+1,krhs)+                             &
     &                  vbar(i-1,j  ,krhs)+                             &
     &                  vbar(i-1,j+1,krhs))
          cff2=SQRT(ubar(i,j,krhs)*ubar(i,j,krhs)+cff1*cff1)
#ifdef RDRG_GRID
          bustr(i,j)=0.5_r8*(rdrg_grid(i,j)+rdrg_grid(i-1,j))*          &
     &                  u(i,j,krhs)*cff2
#else
          bustr(i,j)=rdrg2(ng)*ubar(i,j,krhs)*cff2
#endif
        END DO
      END DO
      DO j=JstrV,Jend
        DO i=Istr,Iend
          cff1=0.25_r8*(ubar(i  ,j  ,krhs)+                             &
     &                  ubar(i+1,j  ,krhs)+                             &
     &                  ubar(i  ,j-1,krhs)+                             &
     &                  ubar(i+1,j-1,krhs))
          cff2=SQRT(cff1*cff1+vbar(i,j,krhs)*vbar(i,j,krhs))
#ifdef RDRG_GRID
          bvstr(i,j)=0.5_r8*(rdrg_grid(i,j)+rdrg_grid(i,j-1))*          &
     &                  vbar(i,j,krhs)*cff2
#else
          bvstr(i,j)=rdrg2(ng)*vbar(i,j,krhs)*cff2
#endif
        END DO
      END DO
#  endif
!
!  Apply boundary conditions.
!
      CALL bc_u2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bustr)
      CALL bc_v2d_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  bvstr)
#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, EWperiodic, NSperiodic,         &
     &                    bustr, bvstr)
#  endif
      RETURN
      END SUBROUTINE set_vbc_tile
# endif
#endif
      END MODULE set_vbc_mod
