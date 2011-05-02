#include "cppdefs.h"
      MODULE aibc_mod
#ifdef ICE_MODEL
! 
!***********************************************************************
!  Compute lateral boundary conditions for the ice concentration.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC aibc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE aibc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL aibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                liold(ng), linew(ng),                             &
     &                ICE(ng) % ui,                                     &
     &                ICE(ng) % vi,                                     &
     &                ICE(ng) % ai)
      RETURN
      END SUBROUTINE aibc

!
!***********************************************************************
      SUBROUTINE aibc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      ui, vi, ai)
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
      integer, intent(in) :: liold, linew

# ifdef ASSUMED_SHAPE
      real(r8), intent(in)    :: ui(LBi:,LBj:,:)
      real(r8), intent(in)    :: vi(LBi:,LBj:,:)
      real(r8), intent(inout) :: ai(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: ai(LBi:UBi,LBj:UBj,2)
# endif

!
!  Local variable declarations.
!
      integer :: i, j, know
      real(r8), parameter :: eps =1.0E-20_r8
      real(r8) :: Ce, Cx, cff, dTde, dTdt, dTdx, tau 

      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Set time-indices
!-----------------------------------------------------------------------
!
        know=liold

# ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN

#  if defined WEST_AIRADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
          grad(Istr-1,j)=ai(Istr-1,j,know)-ai(Istr-1,j-1,know)
#   ifdef MASKING
          grad(Istr-1,j)=grad(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
#   endif
          grad(Istr,j)=ai(Istr,j,know)-ai(Istr,j-1,know)
#   ifdef MASKING
          grad(Istr,j)=grad(Istr,j)*GRID(ng)%vmask(Istr,j)
#   endif
        END DO
        DO j=Jstr,Jend
          dTdt=ai(Istr,j,know)-ai(Istr  ,j,linew)
          dTdx=ai(Istr,j,linew)-ai(Istr+1,j,linew)
#   ifdef WEST_AINUDGING
          tau=Tobc_out(1,ng,iwest)
          IF ((dTdt*dTdx).lt.0.0_r8) tau=Tobc_in(1,ng,iwest)
          tau=tau*dt(ng)
#   endif
          IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
          IF ((dTdt*(grad(Istr,j)+grad(Istr,j+1))).gt.0.0_r8) THEN
            dTde=grad(Istr,j  )
          ELSE
            dTde=grad(Istr,j+1)
          END IF
          cff=MAX(dTdx*dTdx+dTde*dTde,eps)
          Cx=dTdt*dTdx
#   ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dTdt*dTde,-cff))
#   else
          Ce=0.0_r8
#   endif
          ai(Istr-1,j,linew)=(cff*ai(Istr-1,j,know)+                    &
     &                        Cx *ai(Istr  ,j,linew)-                   &
     &                        MAX(Ce,0.0_r8)*grad(Istr-1,j  )-          &
     &                        MIN(Ce,0.0_r8)*grad(Istr-1,j+1))/         &
     &                              (cff+Cx)
#   ifdef WEST_AINUDGING
          ai(Istr-1,j,linew)=ai(Istr-1,j,linew)+                        &
     &                       tau*(BOUNDARY(ng)%ai_west(j)-              &
     &                                   ai(Istr-1,j,know))
#   endif
#   ifdef MASKING
          ai(Istr-1,j,linew)=ai(Istr-1,j,linew)*                        &
     &                              GRID(ng)%rmask(Istr-1,j)
#   endif
        END DO

# elif defined WEST_AICLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ai(0,j,linew)=BOUNDARY(ng)%ai_west(j)
#  ifdef MASKING
          ai(0,j,linew)=ai(0,j,linew)*                                  &
     &                 GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          ai(0,j,linew)=ai(0,j,linew)*                                  &
     &                 GRID(ng)%rmask_wet(0,j)
#  endif
        END DO

# elif defined WEST_AIMIXED
!
!  Western edge, clamped on inflow, gradient on outflow.
!
        DO j=Jstr,Jend
          IF (ui(1,j,linew).ge.0._r8) THEN
            ai(0,j,linew)=BOUNDARY(ng)%ai_west(j)
#  ifdef MASKING
            ai(0,j,linew)=ai(0,j,linew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
            ai(0,j,linew)=ai(0,j,linew)*                                &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ELSE
             ai(0,j,linew)=ai(1,j,liold)
#  ifdef MASKING
             ai(0,j,linew)=ai(0,j,linew)*                               &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             ai(0,j,linew)=ai(0,j,linew)*                               &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          END IF
        END DO
# elif defined WEST_AIGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ai(0,j,linew)=ai(1,j,linew)
#  ifdef MASKING
          ai(0,j,linew)=ai(0,j,linew)*                                  &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          ai(0,j,linew)=ai(0,j,linew)*                                  &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ai(0,j,linew)=ai(1,j,linew)
#  ifdef MASKING
          ai(0,j,linew)=ai(0,j,linew)*                                  &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          ai(0,j,linew)=ai(0,j,linew)*                                  &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# endif

# ifdef CHUKCHI
        DO j=Jstr,Jend
          ai(1,j,linew)=BOUNDARY(ng)%ai_west(j)
#  ifdef MASKING
          ai(1,j,linew)=ai(1,j,linew)*                                  &
     &                 GRID(ng)%rmask(1,j)
#  endif
#  ifdef WET_DRY
          ai(1,j,linew)=ai(1,j,linew)*                                  &
     &                 GRID(ng)%rmask_wet(1,j)
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN

#  if defined EAST_AIRADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
         grad(Iend,j)=ai(Iend,j,know)-ai(Iend,j-1,know)
#   ifdef MASKING
         grad(Iend,j)=grad(Iend,j)*GRID(ng)%vmask(Iend  ,j)
#   endif
         grad(Iend+1,j)=ai(Iend+1,j,know)-ai(Iend+1,j-1,know)
#   ifdef MASKING
         grad(Iend+1,j)=grad(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
#   endif
        END DO
        DO j=Jstr,Jend
          dTdt=ai(Iend,j,know)-ai(Iend  ,j,linew)
          dTdx=ai(Iend,j,linew)-ai(Iend-1,j,linew)
#   ifdef EAST_AINUDGING
          tau=Tobc_out(1,ng,ieast)
          IF ((dTdt*dTdx).lt.0.0_r8) tau=Tobc_in(1,ng,ieast)
          tau=tau*dt(ng)
#   endif
          IF ((dTdt*dTdx).lt.0.0_r8) dTdt=0.0_r8
          IF ((dTdt*(grad(Iend,j)+grad(Iend,j+1))).gt.0.0_r8) THEN
            dTde=grad(Iend,j  )
          ELSE
            dTde=grad(Iend,j+1)
          END IF
          cff=MAX(dTdx*dTdx+dTde*dTde,eps)
          Cx=dTdt*dTdx
#   ifdef RADIATION_2D
          Ce=MIN(cff,MAX(dTdt*dTde,-cff))
#   else
          Ce=0.0_r8
#   endif
          ai(Iend+1,j,linew)=(cff*ai(Iend+1,j,know)+                    &
     &                        Cx *ai(Iend  ,j,linew)-                   &
     &                        MAX(Ce,0.0_r8)*grad(Iend+1,j  )-          &
     &                        MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/         &
     &                              (cff+Cx)
#   ifdef EAST_AINUDGING
          ai(Iend+1,j,linew)=ai(Iend+1,j,linew)+                        &
     &             tau*(BOUNDARY(ng)%ai_east(j)-ai(Iend+1,j,know))
#   endif
#   ifdef MASKING
          ai(Iend+1,j,linew)=ai(Iend+1,j,linew)*                        &
     &                              GRID(ng)%rmask(Iend+1,j)
#   endif
        END DO

# elif defined EAST_AICLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          ai(Lm(ng)+1,j,linew)=BOUNDARY(ng)%ai_east(j)
#  ifdef MASKING
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                    &
     &                        GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                    &
     &                        GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO

# elif defined EAST_AIMIXED
!
!  Eastern edge, clamped on inflow, gradient on outflow.
!
        DO j=Jstr,Jend
          IF (ui(Lm(ng)+1,j,linew).le.0._r8) THEN
             ai(Lm(ng)+1,j,linew)=BOUNDARY(ng)%ai_east(j)
#  ifdef MASKING
             ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ELSE
             ai(Lm(ng)+1,j,linew)=ai(Lm(ng),j,liold)
#  ifdef MASKING
             ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          END IF
        END DO
# elif defined EAST_AIGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng),j,linew)
#  ifdef MASKING
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng),j,linew)
#  ifdef MASKING
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          ai(Lm(ng)+1,j,linew)=ai(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# endif

# ifdef CHUKCHI
        DO j=Jstr,Jend
          ai(Lm(ng),j,linew)=BOUNDARY(ng)%ai_east(j)
#  ifdef MASKING
          ai(Lm(ng),j,linew)=ai(Lm(ng),j,linew)*                        &
     &                        GRID(ng)%rmask(Lm(ng),j)
#  endif
#  ifdef WET_DRY
          ai(Lm(ng),j,linew)=ai(Lm(ng),j,linew)*                        &
     &                        GRID(ng)%rmask_wet(Lm(ng),j)
#  endif
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

#  if defined SOUTH_AIRADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO i=Istr,Iend+1
          grad(i,Jstr)=ai(i,Jstr,know)-ai(i-1,Jstr,know)
#   ifdef MASKING
          grad(i,Jstr)=grad(i,Jstr)*GRID(ng)%umask(i,Jstr)
#   endif
          grad(i,Jstr-1)=ai(i,Jstr-1,know)-ai(i-1,Jstr-1,know)
#   ifdef MASKING
          grad(i,Jstr-1)=grad(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
#   endif
        END DO
        DO i=Istr,Iend
          dTdt=ai(i,Jstr,know)-ai(i,Jstr  ,linew)
          dTde=ai(i,Jstr,linew)-ai(i,Jstr+1,linew)
#   ifdef SOUTH_AINUDGING
          tau=Tobc_out(1,ng,isouth)
          IF ((dTdt*dTde).lt.0.0_r8) tau=Tobc_in(1,ng,isouth)
          tau=tau*dt(ng)
#   endif
          IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
          IF ((dTdt*(grad(i,Jstr)+grad(i+1,Jstr))).gt.0.0_r8) THEN
            dTdx=grad(i  ,Jstr)
          ELSE
            dTdx=grad(i+1,Jstr)
          END IF
          cff=MAX(dTdx*dTdx+dTde*dTde,eps)
#   ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
#   else
          Cx=0.0_r8
#   endif
          Ce=dTdt*dTde
          ai(i,Jstr-1,linew)=(cff*ai(i,Jstr-1,know)+                    &
     &                        Ce *ai(i,Jstr  ,linew)-                   &
     &                        MAX(Cx,0.0_r8)*grad(i  ,Jstr-1)-          &
     &                        MIN(Cx,0.0_r8)*grad(i+1,Jstr-1))/         &
     &                              (cff+Ce)
#   ifdef SOUTH_AINUDGING
          ai(i,Jstr-1,linew)=ai(i,Jstr-1,linew)+                        &
     &           tau*(BOUNDARY(ng)%ai_south(i)-ai(i,Jstr-1,know))
#   endif
#   ifdef MASKING
          ai(i,Jstr-1,linew)=ai(i,Jstr-1,linew)*                        &
     &                              GRID(ng)%rmask(i,Jstr-1)
#   endif
        END DO

# elif defined SOUTH_AICLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          ai(i,0,linew)=BOUNDARY(ng)%ai_south(i)
#  ifdef MASKING
          ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO

# elif defined SOUTH_AIMIXED
!
!  Southern edge, clamped on inflow, gradient on outflow.
!
        DO i=Istr,Iend
          IF (vi(i,1,linew).ge.0._r8) THEN
             ai(i,0,linew)=BOUNDARY(ng)%ai_south(i)
#  ifdef MASKING
             ai(i,0,linew)=ai(i,0,linew)*                               &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             ai(i,0,linew)=ai(i,0,linew)*                               &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ELSE
          ai(i,0,linew)=ai(i,1,liold)
#  ifdef MASKING
          ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          END IF
        END DO
# elif defined SOUTH_AIGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          ai(i,0,linew)=ai(i,1,linew)
#  ifdef MASKING
          ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          ai(i,0,linew)=ai(i,1,linew)
#  ifdef MASKING
          ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          ai(i,0,linew)=ai(i,0,linew)*                                  &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# endif

# ifdef CHUKCHI
        DO i=Istr,Iend
          ai(i,1,linew)=BOUNDARY(ng)%ai_south(i)
#  ifdef MASKING
          ai(i,1,linew)=ai(i,1,linew)*                                  &
     &                   GRID(ng)%rmask(i,1)
#  endif
#  ifdef WET_DRY
          ai(i,1,linew)=ai(i,1,linew)*                                  &
     &                   GRID(ng)%rmask_wet(i,1)
#  endif
        END DO
# endif
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN

#  if defined NORTH_AIRADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO i=Istr,Iend+1
          grad(i,Jend)=ai(i,Jend,know)-ai(i-1,Jend,know)
#   ifdef MASKING
          grad(i,Jend)=grad(i,Jend)*GRID(ng)%umask(i,Jend)
#   endif
          grad(i,Jend+1)=ai(i,Jend+1,know)-ai(i-1,Jend+1,know)
#   ifdef MASKING
          grad(i,Jend+1)=grad(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
#   endif
        END DO
        DO i=Istr,Iend
          dTdt=ai(i,Jend,know)-ai(i,Jend  ,linew)
          dTde=ai(i,Jend,linew)-ai(i,Jend-1,linew)
#   ifdef NORTH_AINUDGING
          tau=Tobc_out(1,ng,inorth)
          IF ((dTdt*dTde).lt.0.0_r8) tau=Tobc_in(1,ng,inorth)
          tau=tau*dt(ng)
#   endif
          IF ((dTdt*dTde).lt.0.0_r8) dTdt=0.0_r8
          IF ((dTdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
            dTdx=grad(i  ,Jend)
          ELSE
            dTdx=grad(i+1,Jend)
          END IF
          cff=MAX(dTdx*dTdx+dTde*dTde,eps)
#   ifdef RADIATION_2D
          Cx=MIN(cff,MAX(dTdt*dTdx,-cff))
#   else
          Cx=0.0_r8
#   endif
          Ce=dTdt*dTde
          ai(i,Jend+1,linew)=(cff*ai(i,Jend+1,know)+                    &
     &                        Ce *ai(i,Jend  ,linew)-                   &
     &                        MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-          &
     &                        MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/         &
     &                              (cff+Ce)
#   ifdef NORTH_AINUDGING
          ai(i,Jend+1,linew)=ai(i,Jend+1,linew)+                        &
     &              tau*(BOUNDARY(ng)%ai_north(i)-ai(i,Jend+1,know))
#   endif
#   ifdef MASKING
          ai(i,Jend+1,linew)=ai(i,Jend+1,linew)*                        &
     &                              GRID(ng)%rmask(i,Jend+1)
#   endif
        END DO

# elif defined NORTH_AICLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          ai(i,Mm(ng)+1,linew)=BOUNDARY(ng)%ai_north(i)
#  ifdef MASKING
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO

# elif defined NORTH_AIMIXED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          IF (vi(i,Mm(ng)+1,linew).le.0._r8) THEN
             ai(i,Mm(ng)+1,linew)=BOUNDARY(ng)%ai_north(i)
#  ifdef MASKING
             ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                 &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                 &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ELSE
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng),liold)
#  ifdef MASKING
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          END IF
        END DO
# elif defined NORTH_AIGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng),linew)
#  ifdef MASKING
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng),linew)
#  ifdef MASKING
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          ai(i,Mm(ng)+1,linew)=ai(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# endif

# ifdef CHUKCHI
        DO i=Istr,Iend
          ai(i,Mm(ng),linew)=BOUNDARY(ng)%ai_north(i)
#  ifdef MASKING
          ai(i,Mm(ng),linew)=ai(i,Mm(ng),linew)*                        &
     &                          GRID(ng)%rmask(i,Mm(ng))
#  endif
#  ifdef WET_DRY
          ai(i,Mm(ng),linew)=ai(i,Mm(ng),linew)*                        &
     &                          GRID(ng)%rmask_wet(i,Mm(ng))
#  endif
        END DO
# endif
      END IF
#endif
# if !defined EW_PERIODIC && !defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
        ai(0,0,linew)=0.5_r8*(ai(1,0,linew)+                            &
     &                         ai(0,1,linew))
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
        ai(Lm(ng)+1,0,linew)=0.5_r8*(ai(Lm(ng)+1,1,linew)+              &
     &                                ai(Lm(ng)  ,0,linew))
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
        ai(0,Mm(ng)+1,linew)=0.5_r8*(ai(0,Mm(ng)  ,linew)+              &
     &                                ai(1,Mm(ng)+1,linew))
      END IF
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        ai(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                             &
     &             (ai(Lm(ng)+1,Mm(ng)  ,linew)+                        &
     &              ai(Lm(ng)  ,Mm(ng)+1,linew))
      END IF
# endif
      RETURN
      END SUBROUTINE aibc_tile
#endif

      END MODULE aibc_mod
