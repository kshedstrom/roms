      SUBROUTINE ana_trc_psource (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
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
      SUBROUTINE ana_trc_psource_tile (ng, tile, model,                 &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Msrcpt,Nsrcpt,                       &
     &                             Isrcpt, Jsrcpt, Lsrcpt,              &
     &                             Tsrcpt  )
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
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

      integer, intent(in) :: Msrcpt
      integer, intent(out) :: Nsrcpt
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
!  to apply (Lsrc). No direction for these sources. The mass sources are 
!  located at rho-points so the grid locations should range from
!  1 =< Isrc =< Lm  and  1 =< Jsrc =< Mm.
!
          Lsrcpt(:,:)=.FALSE.
!  62.6 N from 179.27 to 180.51 E
          Nsrcpt=13
          DO is=1,Nsrcpt
            Isrcpt(is)=286 + is
            Jsrcpt(is)=279
            Lsrcpt(is,inert(1))=.TRUE.
          END DO
!  62 N to 65 N at 194.5 E
          Nsrcpt=Nsrcpt+61
          DO is=1,61
            Isrcpt(is+13)=435.333333 - is*.333333
            Jsrcpt(is+13)=291 + is
            Lsrcpt(is+13,inert(3))=.TRUE.
          END DO
!  69.5 N from 15 to 18 E (roughly)
          Nsrcpt=Nsrcpt+18
          DO is=1,18
            Isrcpt(is+74)=179+is
            Jsrcpt(is+74)=1024
            Lsrcpt(is+74,inert(5))=.TRUE.
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
          DO is=1,13
            Tsrcpt(is,k,inert(1))=1.
          END DO
        END DO
        DO k=1,N(ng)
          DO is=1,61
            Tsrcpt(is+13,k,inert(3))=1.
          END DO
        END DO
        DO k=1,N(ng)
          DO is=1,18
            Tsrcpt(is+74,k,inert(5))=1.
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE ana_trc_psource_tile
