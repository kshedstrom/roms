#include "cppdefs.h"
      MODULE uibc_mod
#ifdef ICE_MODEL
!***********************************************************************
!  Compute the lateral boundary conditions on the ice U-velocity.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC uibc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE uibc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

#include "tile.h"
!
      CALL  uibc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 liuol(ng), liunw(ng),                            &
     &                 ICE(ng) % ui)
      RETURN
      END SUBROUTINE uibc
!
!***********************************************************************
      SUBROUTINE uibc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liuol, liunw,                               &
     &                      ui)
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
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: liuol, liunw

# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: ui(LBi:,LBj:,:)
# else
      real(r8), intent(inout) :: ui(LBi:UBi,LBj:UBj,2)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, know
      real(r8), parameter :: eps =1.0E-20_r8
      real(r8) :: Ce, Cx
      real(r8) :: cff, dUde, dUdt, dUdx, tau

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
        know=liuol
!
#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN

# if defined WEST_M2RADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
          grad(Istr  ,j)=ui(Istr  ,j  ,know)-                           &
     &                   ui(Istr  ,j-1,know)
          grad(Istr+1,j)=ui(Istr+1,j  ,know)-                           &
     &                   ui(Istr+1,j-1,know)
        END DO
        DO j=Jstr,Jend
          dUdt=ui(Istr+1,j,know)-ui(Istr+1,j,liunw)
          dUdx=ui(Istr+1,j,liunw)-ui(Istr+2,j,liunw)
#  ifdef WEST_M2NUDGING
          IF ((dUdt*dUdx).lt.0.0_r8) THEN
            tau=M2obc_in(ng,iwest)
          ELSE
            tau=M2obc_out(ng,iwest)
          END IF
          tau=tau*dt(ng)
#  endif
          IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(Istr+1,j)+grad(Istr+1,j+1))).gt.0.0_r8)
THEN
            dUde=grad(Istr+1,j  )
          ELSE
            dUde=grad(Istr+1,j+1)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
          Cx=dUdt*dUdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dUdt*dUde,-cff))
#  else
          Ce=0.0_r8
#  endif
          ui(Istr,j,liunw)=(cff*ui(Istr  ,j,know)+                      &
     &                       Cx *ui(Istr+1,j,liunw)-                    &
     &                       MAX(Ce,0.0_r8)*grad(Istr,j  )-             &
     &                       MIN(Ce,0.0_r8)*grad(Istr,j+1))/            &
     &                      (cff+Cx)
#  ifdef WEST_M2NUDGING
         ui(Istr,j,liunw)=ui(Istr,j,liunw)+                             &
     &                     tau*(BOUNDARY(ng)%ui_west(j)-                &
     &                          ui(Istr,j,know))
#  endif
#  ifdef MASKING
          ui(Istr,j,liunw)=ui(Istr,j,liunw)*                            &
     &                      GRID(ng)%umask(Istr,j)
#  endif
        END DO

# elif defined WEST_MICLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ui(1,j,liunw)=BOUNDARY(ng)%ui_west(j)
#  ifdef MASKING
          ui(1,j,liunw)=ui(1,j,liunw)*                                  &
     &                   GRID(ng)%umask(1,j)
#  endif
#  ifdef WET_DRY
          ui(1,j,liunw)=ui(1,j,liunw)*                                  &
     &                   GRID(ng)%umask_wet(1,j)
#  endif
        END DO
# elif defined WEST_MIGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ui(1,j,liunw)=ui(2,j,liunw)
#  ifdef MASKING
          ui(1,j,liunw)=ui(1,j,liunw)*                                  &
     &                   GRID(ng)%umask(1,j)
#  endif
#  ifdef WET_DRY
          ui(1,j,liunw)=ui(1,j,liunw)*                                  &
     &                   GRID(ng)%umask_wet(1,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ui(1,j,liunw)=0.0_r8
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN

# if defined EAST_M2RADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
          grad(Iend  ,j)=ui(Iend  ,j  ,know)-                           &
     &                   ui(Iend  ,j-1,know)
          grad(Iend+1,j)=ui(Iend+1,j  ,know)-                           &
     &                   ui(Iend+1,j-1,know)
        END DO
        DO j=Jstr,Jend
          dUdt=ui(Iend,j,know)-ui(Iend  ,j,liunw)
          dUdx=ui(Iend,j,liunw)-ui(Iend-1,j,liunw)
#  ifdef EAST_M2NUDGING
          IF ((dUdt*dUdx).lt.0.0_r8) THEN
            tau=M2obc_in(ng,ieast)
          ELSE
            tau=M2obc_out(ng,ieast)
          END IF
          tau=tau*dt(ng)
#  endif
          IF ((dUdt*dUdx).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
            dUde=grad(Iend,j)
          ELSE
            dUde=grad(Iend,j+1)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
          Cx=dUdt*dUdx
#  ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dUdt*dUde,-cff))
#  else
          Ce=0.0_r8
#  endif
          ui(Iend+1,j,liunw)=(cff*ui(Iend+1,j,know)+                    &
     &                         Cx *ui(Iend  ,j,liunw)-                  &
     &                         MAX(Ce,0.0_r8)*grad(Iend+1,j  )-         &
     &                         MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/        &
     &                        (cff+Cx)
#  ifdef EAST_M2NUDGING
          ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)+ &
     &                        tau*(BOUNDARY(ng)%ui_east(j)-             &
     &                             ui(Iend+1,j,know))
#  endif
#  ifdef MASKING
          ui(Iend+1,j,liunw)=ui(Iend+1,j,liunw)*                        &
     &                        GRID(ng)%umask(Iend+1,j)
#  endif
        END DO

# elif defined EAST_MICLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ui(Lm(ng)+1,j,liunw)=BOUNDARY(ng)%ui_east(j)
#  ifdef MASKING
          ui(Lm(ng)+1,j,liunw)=ui(Lm(ng)+1,j,liunw)*                    &
     &                          GRID(ng)%umask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          ui(Lm(ng)+1,j,liunw)=ui(Lm(ng)+1,j,liunw)*                    &
     &                          GRID(ng)%umask_wet(Lm(ng)+1,j)
#  endif
        END DO
# elif defined EAST_MIGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ui(Lm(ng)+1,j,liunw)=ui(Lm(ng),j,liunw)
#  ifdef MASKING
          ui(Lm(ng)+1,j,liunw)=ui(Lm(ng)+1,j,liunw)*                    &
     &                          GRID(ng)%umask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          ui(Lm(ng)+1,j,liunw)=ui(Lm(ng)+1,j,liunw)*                    &
     &                          GRID(ng)%umask_wet(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ui(Lm(ng)+1,j,liunw)=0.0_r8
        END DO
# endif
      END IF
#endif
#ifndef NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN

# if defined SOUTH_M2RADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO i=IstrU-1,Iend
          grad(i,Jstr-1)=ui(i+1,Jstr-1,know)-                           &
     &                   ui(i  ,Jstr-1,know)
          grad(i,Jstr  )=ui(i+1,Jstr  ,know)-                           &
     &                   ui(i  ,Jstr  ,know)
        END DO
        DO i=IstrU,Iend
          dUdt=ui(i,Jstr,know)-ui(i,Jstr  ,liunw)
          dUde=ui(i,Jstr,liunw)-ui(i,Jstr+1,liunw)
#  ifdef SOUTH_M2NUDGING
          IF ((dUdt*dUde).lt.0.0_r8) THEN
            tau=M2obc_in(ng,isouth)
          ELSE
            tau=M2obc_out(ng,isouth)
          END IF
          tau=tau*dt(ng)
#  endif
          IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(i-1,Jstr)+grad(i,Jstr))).gt.0.0_r8) THEN
            dUdx=grad(i-1,Jstr)
          ELSE
            dUdx=grad(i  ,Jstr)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dUdt*dUde
          ui(i,Jstr-1,liunw)=(cff*ui(i,Jstr-1,know)+                    &
     &                         Ce *ui(i,Jstr  ,liunw)-                  &
     &                         MAX(Cx,0.0_r8)*grad(i-1,Jstr-1)-         &
     &                         MIN(Cx,0.0_r8)*grad(i  ,Jstr-1))/        &
     &                        (cff+Ce)
#  ifdef SOUTH_M2NUDGING
          ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)+                        &
     &                        tau*(BOUNDARY(ng)%ui_south(i)-            &
     &                             ui(i,Jstr-1,know))
#  endif
#  ifdef MASKING
          ui(i,Jstr-1,liunw)=ui(i,Jstr-1,liunw)*                        &
     &                        GRID(ng)%umask(i,Jstr-1)
#  endif
        END DO

# elif defined SOUTH_MICLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=IstrU,Iend
          ui(i,0,liunw)=BOUNDARY(ng)%ui_south(i)
#  ifdef MASKING
          ui(i,0,liunw)=ui(i,0,liunw)*                                  &
     &                   GRID(ng)%umask(i,0)
#  endif
#  ifdef WET_DRY
          ui(i,0,liunw)=ui(i,0,liunw)*                                  &
     &                   GRID(ng)%umask_wet(i,0)
#  endif
        END DO
# elif defined SOUTH_MIGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=IstrU,Iend
          ui(i,0,liunw)=ui(i,1,liunw)
#  ifdef MASKING
          ui(i,0,liunw)=ui(i,0,liunw)*                                  &
     &                   GRID(ng)%umask(i,0)
#  endif
#  ifdef WET_DRY
          ui(i,0,liunw)=ui(i,0,liunw)*                                  &
     &                   GRID(ng)%umask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          ui(i,0,liunw)=gamma2(ng)*ui(i,1,liunw)
#  ifdef MASKING
          ui(i,0,liunw)=ui(i,0,liunw)*                                  &
     &                   GRID(ng)%umask(i,0)
#  endif
#  ifdef WET_DRY
          ui(i,0,liunw)=ui(i,0,liunw)*                                  &
     &                   GRID(ng)%umask_wet(i,0)
#  endif
        END DO
#  undef I_RANGE
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN

# if defined NORTH_M2RADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO i=IstrU-1,Iend
          grad(i,Jend  )=ui(i+1,Jend  ,know)-                           &
     &                   ui(i  ,Jend  ,know)
          grad(i,Jend+1)=ui(i+1,Jend+1,know)-                           &
     &                   ui(i  ,Jend+1,know)
        END DO
        DO i=IstrU,Iend
          dUdt=ui(i,Jend,know)-ui(i,Jend  ,liunw)
          dUde=ui(i,Jend,liunw)-ui(i,Jend-1,liunw)
#  ifdef NORTH_M2NUDGING
          IF ((dUdt*dUde).lt.0.0_r8) THEN
            tau=M2obc_in(ng,inorth)
          ELSE
            tau=M2obc_out(ng,inorth)
          END IF
          tau=tau*dt(ng)
#  endif
          IF ((dUdt*dUde).lt.0.0_r8) dUdt=0.0_r8
          IF ((dUdt*(grad(i-1,Jend)+grad(i,Jend))).gt.0.0_r8) THEN
            dUdx=grad(i-1,Jend)
          ELSE
            dUdx=grad(i  ,Jend)
          END IF
          cff=MAX(dUdx*dUdx+dUde*dUde,eps)
#  ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dUdt*dUdx,-cff))
#  else
          Cx=0.0_r8
#  endif
          Ce=dUdt*dUde
          ui(i,Jend+1,liunw)=(cff*ui(i,Jend+1,know)+                    &
     &                         Ce *ui(i,Jend  ,liunw)-                  &
     &                         MAX(Cx,0.0_r8)*grad(i-1,Jend+1)-         &
     &                         MIN(Cx,0.0_r8)*grad(i  ,Jend+1))/        &
     &                        (cff+Ce)
#  ifdef NORTH_M2NUDGING
          ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)+                        &
     &                        tau*(BOUNDARY(ng)%ui_north(i)-            &
     &                             ui(i,Jend+1,know))
#  endif
#  ifdef MASKING
          ui(i,Jend+1,liunw)=ui(i,Jend+1,liunw)*                        &
     &                        GRID(ng)%umask(i,Jend+1)
#  endif
        END DO

# elif defined NORTH_MICLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=IstrU,Iend
          ui(i,Mm(ng)+1,liunw)=BOUNDARY(ng)%ui_north(i)
#  ifdef MASKING
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%umask_wet(i,Mm(ng)+1)
#  endif
        END DO
# elif defined NORTH_MIGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=IstrU,Iend
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng),liunw)
#  ifdef MASKING
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%umask_wet(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition: free slip (gamma2=1)  or
!                                            no   slip (gamma2=-1).
!
#  ifdef EW_PERIODIC
#   define I_RANGE IstrU,Iend
#  else
#   define I_RANGE Istr,IendR
#  endif
        DO i=I_RANGE
          ui(i,Mm(ng)+1,liunw)=gamma2(ng)*ui(i,Mm(ng),liunw)
#  ifdef MASKING
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%umask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          ui(i,Mm(ng)+1,liunw)=ui(i,Mm(ng)+1,liunw)*                    &
     &                          GRID(ng)%umask_wet(i,Mm(ng)+1)
#  endif
        END DO
#  undef I_RANGE
# endif
      END IF
#endif
#if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
        ui(1,0,liunw)=0.5_r8*(ui(2,0,liunw)+                            &
     &                         ui(1,1,liunw))
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
        ui(Lm(ng)+1,0,liunw)=0.5_r8*(ui(Lm(ng)  ,0,liunw)+              &
     &                                ui(Lm(ng)+1,1,liunw))
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
        ui(1,Mm(ng)+1,liunw)=0.5_r8*(ui(2,Mm(ng)+1,liunw)+              &
     &                                ui(1,Mm(ng)  ,liunw))
      END IF
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        ui(Lm(ng)+1,Mm(ng)+1,liunw)=0.5_r8*                             &
     &                               (ui(Lm(ng)  ,Mm(ng)+1,liunw)+      &
     &                                ui(Lm(ng)+1,Mm(ng)  ,liunw))
      END IF
#endif
      RETURN
      END SUBROUTINE uibc_tile
#endif

      END MODULE uibc_mod
