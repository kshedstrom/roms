#include "cppdefs.h"
      MODULE fish_land_mod
#if defined NONLINEAR && defined NEMURO_SAN
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group         Mark Hadfield   !
!    Licensed under a MIT/X style license             John M. Klinck   !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Reflect or relocate fish due to land.                               !
!                                                                      !
!=======================================================================

      implicit none

      PRIVATE
      PUBLIC  :: fish_land

      CONTAINS
!
!***********************************************************************
      SUBROUTINE fish_land (ng, tile, LBi, UBi, LBj, UBj,               &
     &                          itimem1, itime, itimep1,                &
# ifdef MASKING
     &                          rmask,                                  &
# endif
     &                          fishthread, bounded, track, alive,      &
     &                          fish_count, fish_list, fishnodes)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_biology
      USE mod_grid
      USE mod_types
      USE nrutil
      USE mod_fish
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, LBi, UBi, LBj, UBj
      integer, intent(in) :: itimem1, itime, itimep1
      integer, intent(in) :: fish_count(LBi:UBi,LBj:UBj)

      type(fishnode), intent(in) :: fish_list(LBi:UBi,LBj:UBj)
      type(fishnode), target, intent(in) :: fishnodes(Nfish(ng))

      logical, intent(in) :: fishthread(Nfish(ng))
      logical, intent(in) :: bounded(Nfish(ng))
      logical, intent(in) :: alive(Nfish(ng))

# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfish(ng))
!
!  Local variable declarations.
!
      integer :: i, j, ii, jj, Ir, Jr
      integer :: ifish, fid

      type(fishnode), pointer :: thisfish

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Reflect or relocate fish due to land boundary
!-----------------------------------------------------------------------
!
      DO i=MAX(Istr-1,1),MIN(Iend+1,Lm(ng))
        DO j=MAX(Jstr-1,1),MIN(Jend+1,Mm(ng))
          IF (fish_count(i,j).gt.0) THEN
            thisfish => fish_list(i,j) % next
            DO ifish=1,fish_count(i,j)
              fid = thisfish % fish
              IF (fishthread(fid).and.bounded(fid).and.                 &
     &            alive(fid)) THEN
! Approach 1: turn around
                Ir=NINT(track(ixgrd,itimep1,fid))
                Jr=NINT(track(iygrd,itimep1,fid))
# ifdef MASKING
                IF (rmask(Ir,Jr).lt.0.5_r8) THEN
                  track(ixgrd,itimep1,fid)=track(ixgrd,itimem1,fid)
                  track(iygrd,itimep1,fid)=track(iygrd,itimem1,fid)
                  track(ixrhs,itime,fid)=-track(ixrhs,itime,fid)
                  track(iyrhs,itime,fid)=-track(iyrhs,itime,fid)
                END IF
# endif
! Approach 2: relocate
!                Ir=NINT(track(ixgrd,itime,fid))
!                Jr=NINT(track(iygrd,itime,fid))
! West and east cells
!                DO jj=Jr-1,Jr+1
!                  IF ((rmask(Ir-1,jj).lt.0.5_r8).and.                   &
!     &                (track(ixgrd,itime,fid).lt.REAL(Ir,r8)))          &
!     &              track(ixgrd,itimep1,fid)=REAL(Ir,r8)
!                  IF ((rmask(Ir+1,jj).lt.0.5_r8).and.                   &
!     &                (track(ixgrd,itime,fid).gt.REAL(Ir,r8)))          &
!     &              track(ixgrd,itimep1,fid)=REAL(Ir,r8)
!                END DO
! South and north cells
!                DO ii=Ir-1,Ir+1
!                  IF ((rmask(ii,Jr-1).lt.0.5_r8).and.                   &
!     &                (track(iygrd,itime,fid).lt.REAL(Jr,r8)))          &
!     &              track(iygrd,itimep1,fid)=REAL(Jr,r8)
!                  IF ((rmask(ii,Jr+1).lt.0.5_r8).and.                   &
!     &                (track(iygrd,itime,fid).gt.REAL(Jr,r8)))          &
!     &              track(iygrd,itimep1,fid)=REAL(Jr,r8)
!                END DO
              END IF
              thisfish => thisfish % next
            END DO
          END IF
        END DO
      END DO
!
      RETURN
      END SUBROUTINE fish_land
#endif
      END MODULE fish_land_mod
