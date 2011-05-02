#include "cppdefs.h"
      MODULE hibc_mod
#ifdef ICE_MODEL
!
!***********************************************************************
!  Compute lateral boundary conditions for the ice thickness.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC hibc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE hibc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL hibc_tile (ng, tile,                                         &
     &                LBi, UBi, LBj, UBj,                               &
     &                IminS, ImaxS, JminS, JmaxS,                       &
     &                liold(ng), linew(ng),                             &
     &                ICE(ng) % ui,                                     &
     &                ICE(ng) % vi,                                     &
     &                ICE(ng) % hi)
      RETURN
      END SUBROUTINE hibc

!
!***********************************************************************
      SUBROUTINE hibc_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      liold, linew,                               &
     &                      ui, vi, hi)
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
      real(r8), intent(inout) :: hi(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hi(LBi:UBi,LBj:UBj,2)
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

#ifndef EW_PERIODIC
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN

#  if defined WEST_HIRADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
          grad(Istr-1,j)=hi(Istr-1,j,know)-hi(Istr-1,j-1,know)
#   ifdef MASKING
          grad(Istr-1,j)=grad(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
#   endif
          grad(Istr,j)=hi(Istr,j,know)-hi(Istr,j-1,know)
#   ifdef MASKING
          grad(Istr,j)=grad(Istr,j)*GRID(ng)%vmask(Istr,j)
#   endif
        END DO
        DO j=Jstr,Jend
          dTdt=hi(Istr,j,know)-hi(Istr  ,j,linew)
          dTdx=hi(Istr,j,linew)-hi(Istr+1,j,linew)
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
          hi(Istr-1,j,linew)=(cff*hi(Istr-1,j,know)+                    &
     &                        Cx *hi(Istr  ,j,linew)-                   &
     &                        MAX(Ce,0.0_r8)*grad(Istr-1,j  )-          &
     &                        MIN(Ce,0.0_r8)*grad(Istr-1,j+1))/         &
     &                              (cff+Cx)
#   ifdef WEST_AINUDGING
          hi(Istr-1,j,linew)=hi(Istr-1,j,linew)+                        &
     &                       tau*(BOUNDARY(ng)%hi_west(j)-              &
     &                                   hi(Istr-1,j,know))
#   endif
#   ifdef MASKING
          hi(Istr-1,j,linew)=hi(Istr-1,j,linew)*                        &
     &                              GRID(ng)%rmask(Istr-1,j)
#   endif
        END DO

# elif defined WEST_HICLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          hi(0,j,linew)=BOUNDARY(ng)%hi_west(j)
#  ifdef MASKING
          hi(0,j,linew)=hi(0,j,linew)*                                  &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          hi(0,j,linew)=hi(0,j,linew)*                                  &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO

# elif defined WEST_HIMIXED
!
!  Western edge, clamped on inflow, gradient on outflow.
!
        DO j=Jstr,Jend
          IF (ui(1,j,linew).ge.0._r8) THEN
            hi(0,j,linew)=BOUNDARY(ng)%hi_west(j)
#  ifdef MASKING
            hi(0,j,linew)=hi(0,j,linew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
            hi(0,j,linew)=hi(0,j,linew)*                                &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ELSE
            hi(0,j,linew)=hi(1,j,liold)
#  ifdef MASKING
            hi(0,j,linew)=hi(0,j,linew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
            hi(0,j,linew)=hi(0,j,linew)*                                &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          END IF
        END DO
# elif defined WEST_HIGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          hi(0,j,linew)=hi(1,j,linew)
#  ifdef MASKING
          hi(0,j,linew)=hi(0,j,linew)*                                  &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          hi(0,j,linew)=hi(0,j,linew)*                                  &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          hi(0,j,linew)=hi(1,j,linew)
#  ifdef MASKING
          hi(0,j,linew)=hi(0,j,linew)*                                  &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          hi(0,j,linew)=hi(0,j,linew)*                                  &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# endif

# ifdef CHUKCHI
        DO j=Jstr,Jend
          hi(1,j,linew)=BOUNDARY(ng)%hi_west(j)
#  ifdef MASKING
          hi(1,j,linew)=hi(1,j,linew)*                                  &
     &                   GRID(ng)%rmask(1,j)
#  endif
#  ifdef WET_DRY
          hi(1,j,linew)=hi(1,j,linew)*                                  &
     &                   GRID(ng)%rmask_wet(1,j)
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

#  if defined EAST_HIRADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
         grad(Iend,j)=hi(Iend,j,know)-hi(Iend,j-1,know)
#   ifdef MASKING
         grad(Iend,j)=grad(Iend,j)*GRID(ng)%vmask(Iend  ,j)
#   endif
         grad(Iend+1,j)=hi(Iend+1,j,know)-hi(Iend+1,j-1,know)
#   ifdef MASKING
         grad(Iend+1,j)=grad(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
#   endif
        END DO
        DO j=Jstr,Jend
          dTdt=hi(Iend,j,know)-hi(Iend  ,j,linew)
          dTdx=hi(Iend,j,linew)-hi(Iend-1,j,linew)
#   ifdef EAST_HINUDGING
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
          hi(Iend+1,j,linew)=(cff*hi(Iend+1,j,know)+                    &
     &                        Cx *hi(Iend  ,j,linew)-                   &
     &                        MAX(Ce,0.0_r8)*grad(Iend+1,j  )-          &
     &                        MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/         &
     &                              (cff+Cx)
#   ifdef EAST_HINUDGING
          hi(Iend+1,j,linew)=hi(Iend+1,j,linew)+                        &
     &             tau*(BOUNDARY(ng)%hi_east(j)-hi(Iend+1,j,know))
#   endif
#   ifdef MASKING
          hi(Iend+1,j,linew)=hi(Iend+1,j,linew)*                        &
     &                              GRID(ng)%rmask(Iend+1,j)
#   endif
        END DO

# elif defined EAST_HICLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          hi(Lm(ng)+1,j,linew)=BOUNDARY(ng)%hi_east(j)
#  ifdef MASKING
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO

# elif defined EAST_HIMIXED
!
!  Eastern edge, clamped on inflow, gradient on outflow.
!
        DO j=Jstr,Jend
          IF (ui(Lm(ng)+1,j,linew).le.0._r8) THEN
             hi(Lm(ng)+1,j,linew)=BOUNDARY(ng)%hi_east(j)
#  ifdef MASKING
             hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ELSE
             hi(Lm(ng)+1,j,linew)=hi(Lm(ng),j,liold)
#  ifdef MASKING
             hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                 &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          END IF
        END DO
# elif defined EAST_HIGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng),j,linew)
#  ifdef MASKING
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask(Lm(ng),j)
#  endif
#  ifdef WET_DRY
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask_wet(Lm(ng),j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng),j,linew)
#  ifdef MASKING
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask(Lm(ng),j)
#  endif
#  ifdef WET_DRY
          hi(Lm(ng)+1,j,linew)=hi(Lm(ng)+1,j,linew)*                    &
     &                          GRID(ng)%rmask_wet(Lm(ng),j)
#  endif
        END DO
# endif

# ifdef CHUKCHI
        DO j=Jstr,Jend
          hi(Lm(ng),j,linew)=BOUNDARY(ng)%hi_east(j)
#  ifdef MASKING
          hi(Lm(ng),j,linew)=hi(Lm(ng),j,linew)*                        &
     &                          GRID(ng)%rmask(Lm(ng),j)
#  endif
#  ifdef WET_DRY
          hi(Lm(ng),j,linew)=hi(Lm(ng),j,linew)*                        &
     &                          GRID(ng)%rmask_wet(Lm(ng),j)
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

#  if defined SOUTH_HIRADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO i=Istr,Iend+1
          grad(i,Jstr)=hi(i,Jstr,know)-hi(i-1,Jstr,know)
#   ifdef MASKING
          grad(i,Jstr)=grad(i,Jstr)*GRID(ng)%umask(i,Jstr)
#   endif
          grad(i,Jstr-1)=hi(i,Jstr-1,know)-hi(i-1,Jstr-1,know)
#   ifdef MASKING
          grad(i,Jstr-1)=grad(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
#   endif
        END DO
        DO i=Istr,Iend
          dTdt=hi(i,Jstr,know)-hi(i,Jstr  ,linew)
          dTde=hi(i,Jstr,linew)-hi(i,Jstr+1,linew)
#   ifdef SOUTH_HINUDGING
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
          hi(i,Jstr-1,linew)=(cff*hi(i,Jstr-1,know)+                    &
     &                        Ce *hi(i,Jstr  ,linew)-                   &
     &                        MAX(Cx,0.0_r8)*grad(i  ,Jstr-1)-          &
     &                        MIN(Cx,0.0_r8)*grad(i+1,Jstr-1))/         &
     &                              (cff+Ce)
#   ifdef SOUTH_AINUDGING
          hi(i,Jstr-1,linew)=hi(i,Jstr-1,linew)+                        &
     &           tau*(BOUNDARY(ng)%hi_south(i)-hi(i,Jstr-1,know))
#   endif
#   ifdef MASKING
          hi(i,Jstr-1,linew)=hi(i,Jstr-1,linew)*                        &
     &                              GRID(ng)%rmask(i,Jstr-1)
#   endif
        END DO

# elif defined SOUTH_HICLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          hi(i,0,linew)=BOUNDARY(ng)%hi_south(i)
#  ifdef MASKING
          hi(i,0,linew)=hi(i,0,linew)*                                  &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          hi(i,0,linew)=hi(i,0,linew)*                                  &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO

# elif defined SOUTH_HIMIXED
!
!  Southern edge, clamped on inflow, gradient on outflow.
!
        DO i=Istr,Iend
          IF (vi(i,1,linew).ge.0._r8) THEN
             hi(i,0,linew)=BOUNDARY(ng)%hi_south(i)
#  ifdef MASKING
             hi(i,0,linew)=hi(i,0,linew)*                               &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             hi(i,0,linew)=hi(i,0,linew)*                               &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ELSE
             hi(i,0,linew)=hi(i,1,linew)
#  ifdef MASKING
             hi(i,0,linew)=hi(i,0,linew)*                               &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             hi(i,0,linew)=hi(i,0,linew)*                               &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          END IF
        END DO
# elif defined SOUTH_HIGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          hi(i,0,linew)=hi(i,1,linew)
#  ifdef MASKING
          hi(i,0,linew)=hi(i,0,linew)*                                  &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          hi(i,0,linew)=hi(i,0,linew)*                                  &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          hi(i,0,linew)=hi(i,1,linew)
#  ifdef MASKING
          hi(i,0,linew)=hi(i,0,linew)*                                  &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          hi(i,0,linew)=hi(i,0,linew)*                                  &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# endif

# ifdef CHUKCHI
        DO i=Istr,Iend
          hi(i,1,linew)=BOUNDARY(ng)%hi_south(i)
#  ifdef MASKING
          hi(i,1,linew)=hi(i,1,linew)*                                  &
     &                   GRID(ng)%rmask(i,1)
#  endif
#  ifdef WET_DRY
          hi(i,1,linew)=hi(i,1,linew)*                                  &
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

#  if defined NORTH_HIRADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO i=Istr,Iend+1
          grad(i,Jend)=hi(i,Jend,know)-hi(i-1,Jend,know)
#   ifdef MASKING
          grad(i,Jend)=grad(i,Jend)*GRID(ng)%umask(i,Jend)
#   endif
          grad(i,Jend+1)=hi(i,Jend+1,know)-hi(i-1,Jend+1,know)
#   ifdef MASKING
          grad(i,Jend+1)=grad(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
#   endif
        END DO
        DO i=Istr,Iend
          dTdt=hi(i,Jend,know)-hi(i,Jend  ,linew)
          dTde=hi(i,Jend,linew)-hi(i,Jend-1,linew)
#   ifdef NORTH_HINUDGING
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
          hi(i,Jend+1,linew)=(cff*hi(i,Jend+1,know)+                    &
     &                        Ce *hi(i,Jend  ,linew)-                   &
     &                        MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-          &
     &                        MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/         &
     &                              (cff+Ce)
#   ifdef NORTH_HINUDGING
          hi(i,Jend+1,linew)=hi(i,Jend+1,linew)+                        &
     &              tau*(BOUNDARY(ng)%hi_north(i)-hi(i,Jend+1,know))
#   endif
#   ifdef MASKING
          hi(i,Jend+1,linew)=hi(i,Jend+1,linew)*                        &
     &                              GRID(ng)%rmask(i,Jend+1)
#   endif
        END DO

# elif defined NORTH_HICLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          hi(i,Mm(ng)+1,linew)=BOUNDARY(ng)%hi_north(i)
#  ifdef MASKING
          hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO

# elif defined NORTH_HIMIXED
!
!  Northern edge, clamped on inflow, gradient on outflow.
!
        DO i=Istr,Iend
          IF (vi(i,Mm(ng)+1,linew).le.0._r8) THEN
             hi(i,Mm(ng)+1,linew)=BOUNDARY(ng)%hi_north(i)
#  ifdef MASKING
             hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                 &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                 &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
           ELSE
             hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng),liold)
#  ifdef MASKING
             hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                 &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
             hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                 &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
           END IF
        END DO
# elif defined NORTH_HIGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng),linew)
#  ifdef WET_DRY
          hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng),linew)
#  ifdef MASKING
          hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          hi(i,Mm(ng)+1,linew)=hi(i,Mm(ng)+1,linew)*                    &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# endif

# ifdef CHUKCHI
        DO i=Istr,Iend
          hi(i,Mm(ng),linew)=BOUNDARY(ng)%hi_north(i)
#  ifdef MASKING
          hi(i,Mm(ng),linew)=hi(i,Mm(ng),linew)*                        &
     &                          GRID(ng)%rmask(i,Mm(ng))
#  endif
#  ifdef WET_DRY
          hi(i,Mm(ng),linew)=hi(i,Mm(ng),linew)*                        &
     &                          GRID(ng)%rmask_wet(i,Mm(ng))
#  endif
        END DO
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
        hi(0,0,linew)=0.5_r8*(hi(1,0,linew)+                            &
     &                         hi(0,1,linew))
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
        hi(Lm(ng)+1,0,linew)=0.5_r8*(hi(Lm(ng)+1,1,linew)+              &
     &                                hi(Lm(ng)  ,0,linew))
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
        hi(0,Mm(ng)+1,linew)=0.5_r8*(hi(0,Mm(ng)  ,linew)+              &
     &                                hi(1,Mm(ng)+1,linew))
      END IF
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        hi(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                             &
     &     (hi(Lm(ng)+1,Mm(ng)  ,linew)+                                &
     &      hi(Lm(ng)  ,Mm(ng)+1,linew))
      END IF
#endif
      RETURN
      END SUBROUTINE hibc_tile
#endif

      END MODULE hibc_mod
