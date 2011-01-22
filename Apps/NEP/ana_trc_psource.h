      SUBROUTINE ana_trc_psource (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets analytical tracer and mass point Sources       !
!  and/or Sinks.  River runoff can be consider as a point source.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_trc_sources
      USE mod_stepping
!
      integer, intent(in) :: ng, tile, model

#include "tile.h"
!
      CALL ana_trc_psource_tile (ng, tile, model,                       &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Msrcpt(ng),Nsrcpt(ng),                     &
     &                       TRC_SOURCES(ng) % Isrcpt,                  &
     &                       TRC_SOURCES(ng) % Jsrcpt,                  &
     &                       TRC_SOURCES(ng) % Lsrcpt,                  &
     &                       TRC_SOURCES(ng) % Tsrcpt   )
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(49)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_trc_psource
!
!***********************************************************************
      SUBROUTINE ana_trc_psource_tile (ng, tile, model, 	        &
     &                             LBi, UBi, LBj, UBj,		        &
     &                             IminS, ImaxS, JminS, JmaxS,	        &
     &                             Msrcpt,Nsrcpt,		        &
     &                             Isrcpt, Jsrcpt, Lsrcpt,	        &
     &                             Tsrcpt  )
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
#ifndef ASSUMED_SHAPE
      USE mod_trc_sources, ONLY : Msrcpt
#endif
#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti, mp_bcastl
      USE distribute_mod, ONLY : mp_collect, mp_reduce
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

      integer, intent(out) :: Msrcpt,Nsrcpt
!
#ifdef ASSUMED_SHAPE
      logical, intent(out) :: Lsrcpt(:,:)

      integer, intent(out) :: Isrcpt(:)
      integer, intent(inout) :: Jsrcpt(:)
      real(r8), intent(out) :: Tsrcpt(:,:,:)
#else
      logical, intent(out) :: Lsrcpt(Msrcpt,NT(ng))

      integer, intent(out) :: Isrcpt(Msrcpt)
      integer, intent(out) :: Jsrcpt(Msrcpt)

      real(r8), intent(out) :: Tsrcpt(Msrcpt,N(ng),NT(ng))
#endif
!
!  Local variable declarations.
!
      integer :: is, i, j, k

#if defined DISTRIBUTE
      real(r8), dimension(2) :: buffer

      character (len=3), dimension(2) :: io_handle
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set tracer and/or mass point sources and/or sink.
!-----------------------------------------------------------------------
!
      IF (iic(ng).eq.ntstart(ng)) THEN
!
!  Set-up point Sources/Sink number (Nsrc), direction (Dsrc), I- and
!  J-grid locations (Isrc,Jsrc), and logical switch for type of tracer
!  to apply (Lsrc). Currently, the direction can be along XI-direction
!  (Dsrc = 0) or along ETA-direction (Dsrc > 0).  The mass sources are
!  located at U- or V-points so the grid locations should range from
!  1 =< Isrc =< L  and  1 =< Jsrc =< M.
!
          Lsrcpt(:,:)=.FALSE.
          Nsrcpt=33
          DO is=1,Nsrcpt
            Isrcpt(is)=73+is
            Jsrcpt(is)=INT(-0.7*Isrcpt(is)+473)
            Lsrcpt(is,inert(1))=.TRUE.
          END DO
          Nsrcpt=33+136
          DO is=34,Nsrcpt
            Jsrcpt(is)=408+is-34
            Isrcpt(is)=INT((Jsrcpt(is)-124)/3.2)
            Lsrcpt(is,inert(2))=.TRUE.
          END DO
#ifdef DISTRIBUTE
!
!  Broadcast point sources/sinks information to all nodes.
!
        CALL mp_bcasti (ng, iNLM, Nsrcpt)
        CALL mp_bcasti (ng, iNLM, Isrcpt)
        CALL mp_bcasti (ng, iNLM, Jsrcpt)
        CALL mp_bcastl (ng, iNLM, Lsrcpt)
#endif
!
!  Set-up vertically integrated mass transport (m3/s) of point
!  Sources/Sinks (positive in the positive U- or V-direction and
!  viceversa).
!
!  Set-up tracer (tracer units) point Sources/Sinks.
!
        DO k=1,N(ng)
          DO is=1,33
            Tsrcpt(is,k,inert(1))=1.
          END DO
        END DO
        DO k=1,N(ng)
          DO is=34,Nsrcpt
            Tsrcpt(is,k,inert(2))=1.
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE ana_trc_psource_tile
