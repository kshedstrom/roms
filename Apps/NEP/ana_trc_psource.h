      SUBROUTINE ana_trc_psource (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
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
!  to apply (Lsrc). No direction for these sources. The mass sources are 
!  located at rho-points so the grid locations should range from
!  1 =< Isrc =< Lm  and  1 =< Jsrc =< Mm.
!
          Lsrcpt(:,:)=.FALSE.
          Nsrcpt=5
          DO is=1,Nsrcpt
            Isrcpt(is)=94+is
            Jsrcpt(is)=402
            Lsrcpt(is,inert(1))=.TRUE.
          END DO
          Nsrcpt=6+231
          Jsrcpt(6:Nsrcpt)=(/ &
     &  407, 407, 406, 405, 405, 405, 405, 406, 406, 407, 407, 408, 409, &
     &  410, 410, 411, 412, 413, 414, 415, 416, 416, 417, 417, 418, 419, &
     &  419, 420, 420, 421, 422, 422, 423, 424, 424, 425, 425, 426, 427, &
     &  427, 428, 429, 430, 430, 431, 431, 432, 433, 433, 433, 434, 434, &
     &  435, 435, 436, 437, 438, 439, 440, 441, 442, 443, 443, 444, 445, &
     &  446, 447, 447, 448, 449, 450, 451, 452, 452, 453, 454, 454, 454, &
     &  455, 455, 456, 456, 456, 457, 457, 458, 458, 459, 459, 460, 460, &
     &  461, 461, 462, 462, 463, 463, 464, 465, 465, 466, 467, 468, 468, &
     &  469, 470, 470, 471, 471, 472, 473, 474, 475, 475, 476, 477, 478, &
     &  479, 479, 480, 480, 480, 481, 481, 482, 483, 484, 484, 485, 485, &
     &  486, 486, 486, 487, 487, 487, 486, 486, 485, 486, 486, 487, 488, &
     &  488, 489, 489, 489, 490, 490, 490, 491, 491, 492, 492, 493, 494, &
     &  495, 495, 496, 496, 497, 497, 498, 498, 499, 500, 501, 501, 502, &
     &  502, 503, 504, 505, 505, 506, 507, 507, 508, 509, 510, 510, 511, &
     &  512, 513, 514, 515, 515, 516, 517, 517, 518, 518, 519, 519, 520, &
     &  520, 521, 521, 522, 522, 523, 523, 524, 524, 525, 525, 526, 526, &
     &  527, 527, 528, 528, 529, 530, 531, 531, 532, 532, 533, 534, 534, &
     &  534, 535, 535, 536, 536, 537, 537, 538, 538, 539 /)

          Jsrcpt(6:Nsrcpt)=(/  &
     &   95,  95,  95,  96,  97,  98,  99,  99, 100, 100, 101, 101, 101, &
     &  101, 101, 100, 100, 100, 100, 100, 100, 101, 101, 102, 102, 102, &
     &  103, 103, 104, 104, 104, 105, 105, 105, 106, 106, 107, 107, 107, &
     &  108, 108, 108, 108, 109, 109, 109, 108, 108, 108, 107, 106, 106, &
     &  105, 106, 106, 106, 106, 106, 106, 106, 106, 106, 107, 107, 107, &
     &  107, 107, 107, 106, 106, 106, 106, 106, 107, 107, 107, 108, 109, &
     &  109, 110, 110, 111, 112, 112, 113, 113, 114, 114, 115, 115, 116, &
     &  116, 117, 117, 118, 118, 119, 119, 119, 120, 120, 120, 120, 120, &
     &  119, 119, 119, 118, 119, 119, 119, 119, 119, 119, 118, 118, 118, &
     &  118, 119, 119, 120, 121, 121, 122, 122, 122, 122, 123, 123, 124, &
     &  124, 125, 126, 126, 127, 127, 128, 128, 129, 129, 130, 130, 130, &
     &  131, 131, 132, 133, 133, 134, 135, 135, 136, 136, 137, 137, 137, &
     &  137, 137, 136, 136, 135, 135, 134, 134, 133, 133, 133, 134, 134, &
     &  135, 135, 135, 135, 136, 136, 136, 137, 137, 137, 137, 138, 138, &
     &  138, 138, 138, 138, 138, 137, 137, 137, 136, 136, 135, 135, 134, &
     &  134, 133, 133, 132, 132, 131, 131, 130, 130, 129, 129, 128, 128, &
     &  127, 127, 126, 126, 125, 125, 125, 125, 124, 124, 123, 123, 123, &
     &  122, 121, 121, 120, 120, 119, 119, 118, 118, 117 /)
          DO is=6,Nsrcpt
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
