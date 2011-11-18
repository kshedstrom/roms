#include "cppdefs.h"
      MODULE icenh4bc_mod
#ifdef ICE_MODEL
#ifdef ICE_BIO
!
!***********************************************************************
!  Compute lateral boundary conditions for snow thickness.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC iceNH4bc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE icenh4bc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL iceNH4bc_tile (ng, tile,                                     &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 liold(ng), linew(ng),                            &
     &                 ICE(ng) % ui,                                    &
     &                 ICE(ng) % vi,                                    &
     &                 ICE(ng) % IceNH4)
      RETURN
      END SUBROUTINE iceNH4bc

!
!***********************************************************************
      SUBROUTINE iceNH4bc_tile (ng, tile,                               &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           liold, linew,                          &
     &                           ui, vi, iceNH4)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_scalars

      implicit none

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: liold, linew

# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: iceNH4(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: iceNH4(LBi:UBi,LBj:UBj,2)
# endif

!
!  Local variable declarations.
!
      integer :: i, j, know

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
        know=liold

!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, clamped boundary condition.
!
        IF (LBC(iwest,isINH4,ng)%clamped)) THEN
          DO j=Jstr,Jend
            IF(ui(1,j,linew).ge.0._r8) THEN
               iceNH4(0,j,linew)=BOUNDARY(ng)%iceNH4_west(j)
#ifdef MASKING
               iceNH4(0,j,linew)=iceNH4(0,j,linew)*                     &
     &                     GRID(ng)%rmask(0,j)
#endif
#ifdef WET_DRY
               iceNH4(0,j,linew)=iceNH4(0,j,linew)*                     &
     &                     GRID(ng)%rmask_wet(0,j)
#endif
            ELSE
               iceNH4(0,j,linew)=iceNH4(1,j,liold)
#ifdef MASKING
               iceNH4(0,j,linew)=iceNH4(0,j,linew)*                     &
     &                     GRID(ng)%rmask(0,j)
#endif
#ifdef WET_DRY
               iceNH4(0,j,linew)=iceNH4(0,j,linew)*                     &
     &                     GRID(ng)%rmask_wet(0,j)
#endif
            ENDIF
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isINH4,ng)%gradient)) THEN
          DO j=Jstr,Jend
            iceNH4(0,j,linew)=iceNH4(1,j,linew)
#ifdef MASKING
            iceNH4(0,j,linew)=iceNH4(0,j,linew)*                        &
     &                     GRID(ng)%rmask(0,j)
#endif
#ifdef WET_DRY
            iceNH4(0,j,linew)=iceNH4(0,j,linew)*                        &
     &                     GRID(ng)%rmask_wet(0,j)
#endif
          END DO
        ELSE
!
!  Western edge, closed boundary condition.
!
          DO j=Jstr,Jend
            iceNH4(0,j,linew)=iceNH4(1,j,linew)
#ifdef MASKING
            iceNH4(0,j,linew)=iceNH4(0,j,linew)*                        &
     &                     GRID(ng)%rmask(0,j)
#endif
#ifdef WET_DRY
            iceNH4(0,j,linew)=iceNH4(0,j,linew)*                        &
     &                     GRID(ng)%rmask_wet(0,j)
#endif
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, clamped boundary condition.
!
        IF (LBC(ieast,isINH4,ng)%clamped)) THEN
          DO j=Jstr,Jend
            IF(ui(Lm(ng)+1,j,linew).le.0._r8) THEN
               iceNH4(Lm(ng)+1,j,linew)=BOUNDARY(ng)%iceNH4_east(j)
#ifdef MASKING
               iceNH4(Lm(ng)+1,j,linew)=iceNH4(Lm(ng)+1,j,linew)*       &
     &                            GRID(ng)%rmask(Lm(ng)+1,j)
#endif
#ifdef WET_DRY
               iceNH4(Lm(ng)+1,j,linew)=iceNH4(Lm(ng)+1,j,linew)*       &
     &                            GRID(ng)%rmask_wet(Lm(ng)+1,j)
#endif
            ELSE
              iceNH4(Lm(ng)+1,j,linew)=iceNH4(Lm(ng),j,liold)
#ifdef MASKING
              iceNH4(Lm(ng)+1,j,linew)=iceNH4(Lm(ng)+1,j,linew)*        &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#endif
#ifdef WET_DRY
              iceNH4(Lm(ng)+1,j,linew)=iceNH4(Lm(ng)+1,j,linew)*        &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#endif
            ENDIF
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isINH4,ng)%gradient)) THEN
          DO j=Jstr,Jend
            iceNH4(Lm(ng)+1,j,linew)=iceNH4(Lm(ng),j,linew)
#ifdef MASKING
            iceNH4(Lm(ng)+1,j,linew)=iceNH4(Lm(ng)+1,j,linew)*          &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#endif
#ifdef WET_DRY
            iceNH4(Lm(ng)+1,j,linew)=iceNH4(Lm(ng)+1,j,linew)*          &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#endif
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE
          DO j=Jstr,Jend
            iceNH4(Lm(ng)+1,j,linew)=iceNH4(Lm(ng),j,linew)
#ifdef MASKING
            iceNH4(Lm(ng)+1,j,linew)=iceNH4(Lm(ng)+1,j,linew)*          &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#endif
#ifdef WET_DRY
            iceNH4(Lm(ng)+1,j,linew)=iceNH4(Lm(ng)+1,j,linew)*          &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#endif
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, clamped boundary condition.
!
        IF (LBC(isouth,isINH4,ng)%clamped)) THEN
          DO i=Istr,Iend
            IF (vi(i,1,linew).ge.0._r8) THEN
              iceNH4(i,0,linew)=BOUNDARY(ng)%iceNH4_south(i)
#ifdef MASKING
              iceNH4(i,0,linew)=iceNH4(i,0,linew)*                      &
     &                   GRID(ng)%rmask(i,0)
#endif
#ifdef WET_DRY
              iceNH4(i,0,linew)=iceNH4(i,0,linew)*                      &
     &                   GRID(ng)%rmask_wet(i,0)
#endif
            ELSE
              iceNH4(i,0,linew)=iceNH4(i,1,liold)
#ifdef MASKING
              iceNH4(i,0,linew)=iceNH4(i,0,linew)*                      &
     &                   GRID(ng)%rmask(i,0)
#endif
#ifdef WET_DRY
              iceNH4(i,0,linew)=iceNH4(i,0,linew)*                      &
     &                   GRID(ng)%rmask_wet(i,0)
#endif
            ENDIF
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isINH4,ng)%gradient)) THEN
          DO i=Istr,Iend
            iceNH4(i,0,linew)=iceNH4(i,1,linew)
#ifdef MASKING
            iceNH4(i,0,linew)=iceNH4(i,0,linew)*                        &
     &                   GRID(ng)%rmask(i,0)
#endif
#ifdef WET_DRY
            iceNH4(i,0,linew)=iceNH4(i,0,linew)*                        &
     &                   GRID(ng)%rmask_wet(i,0)
#endif
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE
          DO i=Istr,Iend
            iceNH4(i,0,linew)=iceNH4(i,1,linew)
#ifdef MASKING
            iceNH4(i,0,linew)=iceNH4(i,0,linew)*                        &
     &                   GRID(ng)%rmask(i,0)
#endif
#ifdef WET_DRY
            iceNH4(i,0,linew)=iceNH4(i,0,linew)*                        &
     &                   GRID(ng)%rmask_wet(i,0)
#endif
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, clamped boundary condition.
!
        IF (LBC(inorth,isINH4,ng)%clamped)) THEN
          DO i=Istr,Iend
            IF(vi(i,Mm(ng)+1,linew).le.0._r8) THEN
              iceNH4(i,Mm(ng)+1,linew)=BOUNDARY(ng)%iceNH4_north(i)
#ifdef MASKING
              iceNH4(i,Mm(ng)+1,linew)=iceNH4(i,Mm(ng)+1,linew)*        &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#endif
#ifdef WET_DRY
              iceNH4(i,Mm(ng)+1,linew)=iceNH4(i,Mm(ng)+1,linew)*        &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#endif
            ELSE
              iceNH4(i,Mm(ng)+1,linew)=iceNH4(i,Mm(ng),liold)
#ifdef MASKING
              iceNH4(i,Mm(ng)+1,linew)=iceNH4(i,Mm(ng)+1,linew)*        &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#endif
#ifdef WET_DRY
              iceNH4(i,Mm(ng)+1,linew)=iceNH4(i,Mm(ng)+1,linew)*        &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#endif
            ENDIF
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isINH4,ng)%gradient)) THEN
          DO i=Istr,Iend
            iceNH4(i,Mm(ng)+1,linew)=iceNH4(i,Mm(ng),linew)
#ifdef MASKING
            iceNH4(i,Mm(ng)+1,linew)=iceNH4(i,Mm(ng)+1,linew)*          &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#endif
#ifdef WET_DRY
            iceNH4(i,Mm(ng)+1,linew)=iceNH4(i,Mm(ng)+1,linew)*          &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#endif
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE
          DO i=Istr,Iend
            iceNH4(i,Mm(ng)+1,linew)=iceNH4(i,Mm(ng),linew)
#ifdef MASKING
            iceNH4(i,Mm(ng)+1,linew)=iceNH4(i,Mm(ng)+1,linew)*          &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#endif
#ifdef WET_DRY
            iceNH4(i,Mm(ng)+1,linew)=iceNH4(i,Mm(ng)+1,linew)*          &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#endif
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.EWperiodic(ng).and. .not.NSperiodic(ng)) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          iceNH4(0,0,linew)=0.5_r8*(iceNH4(1,0,linew)+                  &
     &                         iceNH4(0,1,linew))
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          iceNH4(Lm(ng)+1,0,linew)=0.5_r8*(iceNH4(Lm(ng)+1,1,linew)+    &
     &                                iceNH4(Lm(ng)  ,0,linew))
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          iceNH4(0,Mm(ng)+1,linew)=0.5_r8*(iceNH4(0,Mm(ng)  ,linew)+    &
     &                                iceNH4(1,Mm(ng)+1,linew))
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          iceNH4(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                       &
     &             (iceNH4(Lm(ng)+1,Mm(ng)  ,linew)+                    &
     &              iceNH4(Lm(ng)  ,Mm(ng)+1,linew))
        END IF
      END IF
      RETURN
      END SUBROUTINE icenh4bc_tile
#endif
#endif
      END MODULE icenh4bc_mod
