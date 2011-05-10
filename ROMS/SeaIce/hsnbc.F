#include "cppdefs.h"
      MODULE hsnbc_mod
#ifdef ICE_MODEL
!
!***********************************************************************
!  Compute lateral boundary conditions for snow thickness.
!***********************************************************************

      implicit none

      PRIVATE
      PUBLIC hsnbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE hsnbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ice
      USE mod_stepping
!
      integer, intent(in) :: ng, tile

# include "tile.h"
!
      CALL hsnbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 liold(ng), linew(ng),                            &
     &                 ICE(ng) % ui,                                    &
     &                 ICE(ng) % vi,                                    &
     &                 ICE(ng) % hsn)
      RETURN
      END SUBROUTINE hsnbc

!
!***********************************************************************
      SUBROUTINE hsnbc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       liold, linew,                              &
     &                       ui, vi, hsn)
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
      real(r8), intent(inout) :: hsn(LBi:,LBj:,:)
# else
      real(r8), intent(in)    :: ui(LBi:UBi,LBj:UBj,2)
      real(r8), intent(in)    :: vi(LBi:UBi,LBj:UBj,2)
      real(r8), intent(inout) :: hsn(LBi:UBi,LBj:UBj,2)
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

#  if defined WEST_HSNRADIATION
!
!  Western edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
          grad(Istr-1,j)=hsn(Istr-1,j,know)-hsn(Istr-1,j-1,know)
#   ifdef MASKING
          grad(Istr-1,j)=grad(Istr-1,j)*GRID(ng)%vmask(Istr-1,j)
#   endif
          grad(Istr,j)=hsn(Istr,j,know)-hsn(Istr,j-1,know)
#   ifdef MASKING
          grad(Istr,j)=grad(Istr,j)*GRID(ng)%vmask(Istr,j)
#   endif
        END DO
        DO j=Jstr,Jend
          dTdt=hsn(Istr,j,know)-hsn(Istr  ,j,linew)
          dTdx=hsn(Istr,j,linew)-hsn(Istr+1,j,linew)
#   ifdef WEST_HSNNUDGING
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
          hsn(Istr-1,j,linew)=(cff*hsn(Istr-1,j,know)+                  &
     &                        Cx *hsn(Istr  ,j,linew)-                  &
     &                        MAX(Ce,0.0_r8)*grad(Istr-1,j  )-          &
     &                        MIN(Ce,0.0_r8)*grad(Istr-1,j+1))/         &
     &                              (cff+Cx)
#   ifdef WEST_HSNNUDGING
          hsn(Istr-1,j,linew)=hsn(Istr-1,j,linew)+                      &
     &                       tau*(BOUNDARY(ng)%hsn_west(j)-             &
     &                                   hsn(Istr-1,j,know))
#   endif
#   ifdef MASKING
          hsn(Istr-1,j,linew)=hsn(Istr-1,j,linew)*                      &
     &                              GRID(ng)%rmask(Istr-1,j)
#   endif
        END DO

# elif defined WEST_HSNCLAMPED
!
!  Western edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          hsn(0,j,linew)=BOUNDARY(ng)%hsn_west(j)
#  ifdef MASKING
          hsn(0,j,linew)=hsn(0,j,linew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          hsn(0,j,linew)=hsn(0,j,linew)*                                &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO

# elif defined WEST_HSNMIXED
!
!  Western edge, clamped on inflow, gradient on outflow.
!
        DO j=Jstr,Jend
          IF (ui(1,j,linew).ge.0._r8) THEN
             hsn(0,j,linew)=BOUNDARY(ng)%hsn_west(j)
#  ifdef MASKING
             hsn(0,j,linew)=hsn(0,j,linew)*                             &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             hsn(0,j,linew)=hsn(0,j,linew)*                             &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          ELSE
             hsn(0,j,linew)=hsn(1,j,liold)
#  ifdef MASKING
             hsn(0,j,linew)=hsn(0,j,linew)*                             &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
             hsn(0,j,linew)=hsn(0,j,linew)*                             &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
          END IF
        END DO
# elif defined WEST_HSNGRADIENT
!
!  Western edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          hsn(0,j,linew)=hsn(1,j,linew)
#  ifdef MASKING
          hsn(0,j,linew)=hsn(0,j,linew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          hsn(0,j,linew)=hsn(0,j,linew)*                                &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# else
!
!  Western edge, closed boundary condition.
!
        DO j=Jstr,Jend
          hsn(0,j,linew)=hsn(1,j,linew)
#  ifdef MASKING
          hsn(0,j,linew)=hsn(0,j,linew)*                                &
     &                   GRID(ng)%rmask(0,j)
#  endif
#  ifdef WET_DRY
          hsn(0,j,linew)=hsn(0,j,linew)*                                &
     &                   GRID(ng)%rmask_wet(0,j)
#  endif
        END DO
# endif

# ifdef CHUKCHI
        DO j=Jstr,Jend
          hsn(1,j,linew)=BOUNDARY(ng)%hsn_west(j)
#  ifdef MASKING
          hsn(1,j,linew)=hsn(1,j,linew)*                                &
     &                   GRID(ng)%rmask(1,j)
#  endif
#  ifdef WET_DRY
          hsn(1,j,linew)=hsn(1,j,linew)*                                &
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

#  if defined EAST_HSNRADIATION
!
!  Eastern edge, implicit upstream radiation condition.
!
        DO j=Jstr,Jend+1
         grad(Iend,j)=hsn(Iend,j,know)-hsn(Iend,j-1,know)
#   ifdef MASKING
         grad(Iend,j)=grad(Iend,j)*GRID(ng)%vmask(Iend  ,j)
#   endif
         grad(Iend+1,j)=hsn(Iend+1,j,know)-hsn(Iend+1,j-1,know)
#   ifdef MASKING
         grad(Iend+1,j)=grad(Iend+1,j)*GRID(ng)%vmask(Iend+1,j)
#   endif
        END DO
        DO j=Jstr,Jend
          dTdt=hsn(Iend,j,know)-hsn(Iend  ,j,linew)
          dTdx=hsn(Iend,j,linew)-hsn(Iend-1,j,linew)
#   ifdef EAST_HSNNUDGING
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
          hsn(Iend+1,j,linew)=(cff*hsn(Iend+1,j,know)+                  &
     &                        Cx *hsn(Iend  ,j,linew)-                  &
     &                        MAX(Ce,0.0_r8)*grad(Iend+1,j  )-          &
     &                        MIN(Ce,0.0_r8)*grad(Iend+1,j+1))/         &
     &                              (cff+Cx)
#   ifdef EAST_HSNNUDGING
          hsn(Iend+1,j,linew)=hsn(Iend+1,j,linew)+                      &
     &             tau*(BOUNDARY(ng)%hsn_east(j)-hsn(Iend+1,j,know))
#   endif
#   ifdef MASKING
          hsn(Iend+1,j,linew)=hsn(Iend+1,j,linew)*                      &
     &                              GRID(ng)%rmask(Iend+1,j)
#   endif
        END DO

# elif defined EAST_HSNCLAMPED
!
!  Eastern edge, clamped boundary condition.
!
        DO j=Jstr,Jend
          hsn(Lm(ng)+1,j,linew)=BOUNDARY(ng)%hsn_east(j)
#  ifdef MASKING
          hsn(Lm(ng)+1,j,linew)=hsn(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          hsn(Lm(ng)+1,j,linew)=hsn(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO

# elif defined EAST_HSNMIXED
!
!  Eastern edge, clamped on inflow, gradient on outflow.
!
        DO j=Jstr,Jend
          IF (ui(Lm(ng)+1,j,linew).le.0._r8) THEN
             hsn(Lm(ng)+1,j,linew)=BOUNDARY(ng)%hsn_east(j)
#  ifdef MASKING
             hsn(Lm(ng)+1,j,linew)=hsn(Lm(ng)+1,j,linew)*               &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             hsn(Lm(ng)+1,j,linew)=hsn(Lm(ng)+1,j,linew)*               &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          ELSE
             hsn(Lm(ng)+1,j,linew)=hsn(Lm(ng),j,liold)
#  ifdef MASKING
             hsn(Lm(ng)+1,j,linew)=hsn(Lm(ng)+1,j,linew)*               &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
             hsn(Lm(ng)+1,j,linew)=hsn(Lm(ng)+1,j,linew)*               &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
          END IF
        END DO
# elif defined EAST_HSNGRADIENT
!
!  Eastern edge, gradient boundary condition.
!
        DO j=Jstr,Jend
          hsn(Lm(ng)+1,j,linew)=hsn(Lm(ng),j,linew)
#  ifdef MASKING
          hsn(Lm(ng)+1,j,linew)=hsn(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          hsn(Lm(ng)+1,j,linew)=hsn(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# else
!
!  Eastern edge, closed boundary condition.
!
        DO j=Jstr,Jend
          hsn(Lm(ng)+1,j,linew)=hsn(Lm(ng),j,linew)
#  ifdef MASKING
          hsn(Lm(ng)+1,j,linew)=hsn(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask(Lm(ng)+1,j)
#  endif
#  ifdef WET_DRY
          hsn(Lm(ng)+1,j,linew)=hsn(Lm(ng)+1,j,linew)*                  &
     &                          GRID(ng)%rmask_wet(Lm(ng)+1,j)
#  endif
        END DO
# endif

# ifdef CHUKCHI
        DO j=Jstr,Jend
          hsn(Lm(ng),j,linew)=BOUNDARY(ng)%hsn_east(j)
#  ifdef MASKING
          hsn(Lm(ng),j,linew)=hsn(Lm(ng),j,linew)*                      &
     &                          GRID(ng)%rmask(Lm(ng),j)
#  endif
#  ifdef WET_DRY
          hsn(Lm(ng),j,linew)=hsn(Lm(ng),j,linew)*                      &
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

#  if defined SOUTH_HSNRADIATION
!
!  Southern edge, implicit upstream radiation condition.
!
        DO i=Istr,Iend+1
          grad(i,Jstr)=hsn(i,Jstr,know)-hsn(i-1,Jstr,know)
#   ifdef MASKING
          grad(i,Jstr)=grad(i,Jstr)*GRID(ng)%umask(i,Jstr)
#   endif
          grad(i,Jstr-1)=hsn(i,Jstr-1,know)-hsn(i-1,Jstr-1,know)
#   ifdef MASKING
          grad(i,Jstr-1)=grad(i,Jstr-1)*GRID(ng)%umask(i,Jstr-1)
#   endif
        END DO
        DO i=Istr,Iend
          dTdt=hsn(i,Jstr,know)-hsn(i,Jstr  ,linew)
          dTde=hsn(i,Jstr,linew)-hsn(i,Jstr+1,linew)
#   ifdef SOUTH_HSNNUDGING
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
          hsn(i,Jstr-1,linew)=(cff*hsn(i,Jstr-1,know)+                  &
     &                        Ce *hsn(i,Jstr  ,linew)-                  &
     &                        MAX(Cx,0.0_r8)*grad(i  ,Jstr-1)-          &
     &                        MIN(Cx,0.0_r8)*grad(i+1,Jstr-1))/         &
     &                              (cff+Ce)
#   ifdef SOUTH_HSNNUDGING
          hsn(i,Jstr-1,linew)=hsn(i,Jstr-1,linew)+                      &
     &           tau*(BOUNDARY(ng)%hsn_south(i)-hsn(i,Jstr-1,know))
#   endif
#   ifdef MASKING
          hsn(i,Jstr-1,linew)=hsn(i,Jstr-1,linew)*                      &
     &                              GRID(ng)%rmask(i,Jstr-1)
#   endif
        END DO

# elif defined SOUTH_HSNCLAMPED
!
!  Southern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          hsn(i,0,linew)=BOUNDARY(ng)%hsn_south(i)
#  ifdef MASKING
          hsn(i,0,linew)=hsn(i,0,linew)*                                &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          hsn(i,0,linew)=hsn(i,0,linew)*                                &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO

# elif defined SOUTH_HSNMIXED
!
!  Southern edge, clamped on inflow, gradient on outflow.
!
        DO i=Istr,Iend
          IF (vi(i,1,linew).ge.0._r8) THEN
             hsn(i,0,linew)=BOUNDARY(ng)%hsn_south(i)
#  ifdef MASKING
             hsn(i,0,linew)=hsn(i,0,linew)*                             &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             hsn(i,0,linew)=hsn(i,0,linew)*                             &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          ELSE
             hsn(i,0,linew)=hsn(i,1,liold)
#  ifdef MASKING
             hsn(i,0,linew)=hsn(i,0,linew)*                             &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
             hsn(i,0,linew)=hsn(i,0,linew)*                             &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
          END IF
        END DO
# elif defined SOUTH_HSNGRADIENT
!
!  Southern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          hsn(i,0,linew)=hsn(i,1,linew)
#  ifdef MASKING
          hsn(i,0,linew)=hsn(i,0,linew)*                                &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          hsn(i,0,linew)=hsn(i,0,linew)*                                &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# else
!
!  Southern edge, closed boundary condition.
!
        DO i=Istr,Iend
          hsn(i,0,linew)=hsn(i,1,linew)
#  ifdef MASKING
          hsn(i,0,linew)=hsn(i,0,linew)*                                &
     &                   GRID(ng)%rmask(i,0)
#  endif
#  ifdef WET_DRY
          hsn(i,0,linew)=hsn(i,0,linew)*                                &
     &                   GRID(ng)%rmask_wet(i,0)
#  endif
        END DO
# endif

# ifdef CHUKCHI
        DO i=Istr,Iend
          hsn(i,1,linew)=BOUNDARY(ng)%hsn_south(i)
#  ifdef MASKING
          hsn(i,1,linew)=hsn(i,1,linew)*                                &
     &                   GRID(ng)%rmask(i,1)
#  endif
#  ifdef WET_DRY
          hsn(i,1,linew)=hsn(i,1,linew)*                                &
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

#  if defined NORTH_HSNRADIATION
!
!  Northern edge, implicit upstream radiation condition.
!
        DO i=Istr,Iend+1
          grad(i,Jend)=hsn(i,Jend,know)-hsn(i-1,Jend,know)
#   ifdef MASKING
          grad(i,Jend)=grad(i,Jend)*GRID(ng)%umask(i,Jend)
#   endif
          grad(i,Jend+1)=hsn(i,Jend+1,know)-hsn(i-1,Jend+1,know)
#   ifdef MASKING
          grad(i,Jend+1)=grad(i,Jend+1)*GRID(ng)%umask(i,Jend+1)
#   endif
        END DO
        DO i=Istr,Iend
          dTdt=hsn(i,Jend,know)-hsn(i,Jend  ,linew)
          dTde=hsn(i,Jend,linew)-hsn(i,Jend-1,linew)
#   ifdef NORTH_HSNNUDGING
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
          hsn(i,Jend+1,linew)=(cff*hsn(i,Jend+1,know)+                  &
     &                        Ce *hsn(i,Jend  ,linew)-                  &
     &                        MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-          &
     &                        MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/         &
     &                              (cff+Ce)
#   ifdef NORTH_HSNNUDGING
          hsn(i,Jend+1,linew)=hsn(i,Jend+1,linew)+                      &
     &              tau*(BOUNDARY(ng)%hsn_north(i)-hsn(i,Jend+1,know))
#   endif
#   ifdef MASKING
          hsn(i,Jend+1,linew)=hsn(i,Jend+1,linew)*                      &
     &                              GRID(ng)%rmask(i,Jend+1)
#   endif
        END DO

# elif defined NORTH_HSNCLAMPED
!
!  Northern edge, clamped boundary condition.
!
        DO i=Istr,Iend
          hsn(i,Mm(ng)+1,linew)=BOUNDARY(ng)%hsn_north(i)
#  ifdef MASKING
          hsn(i,Mm(ng)+1,linew)=hsn(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          hsn(i,Mm(ng)+1,linew)=hsn(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO

# elif defined NORTH_HSNMIXED
!
!  Northern edge, clamped on inflow, gradient on outflow.
!
        DO i=Istr,Iend
          IF (vi(i,Mm(ng)+1,linew).le.0._r8) THEN
            hsn(i,Mm(ng)+1,linew)=BOUNDARY(ng)%hsn_north(i)
#  ifdef MASKING
            hsn(i,Mm(ng)+1,linew)=hsn(i,Mm(ng)+1,linew)*                &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
            hsn(i,Mm(ng)+1,linew)=hsn(i,Mm(ng)+1,linew)*                &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          ELSE
            hsn(i,Mm(ng)+1,linew)=hsn(i,Mm(ng),liold)
#  ifdef MASKING
            hsn(i,Mm(ng)+1,linew)=hsn(i,Mm(ng)+1,linew)*                &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
            hsn(i,Mm(ng)+1,linew)=hsn(i,Mm(ng)+1,linew)*                &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
          END IF
        END DO
# elif defined NORTH_HSNGRADIENT
!
!  Northern edge, gradient boundary condition.
!
        DO i=Istr,Iend
          hsn(i,Mm(ng)+1,linew)=hsn(i,Mm(ng),linew)
#  ifdef MASKING
          hsn(i,Mm(ng)+1,linew)=hsn(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          hsn(i,Mm(ng)+1,linew)=hsn(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# else
!
!  Northern edge, closed boundary condition.
!
        DO i=Istr,Iend
          hsn(i,Mm(ng)+1,linew)=hsn(i,Mm(ng),linew)
#  ifdef MASKING
          hsn(i,Mm(ng)+1,linew)=hsn(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask(i,Mm(ng)+1)
#  endif
#  ifdef WET_DRY
          hsn(i,Mm(ng)+1,linew)=hsn(i,Mm(ng)+1,linew)*                  &
     &                          GRID(ng)%rmask_wet(i,Mm(ng)+1)
#  endif
        END DO
# endif

# ifdef CHUKCHI
        DO i=Istr,Iend
          hsn(i,Mm(ng),linew)=BOUNDARY(ng)%hsn_north(i)
#  ifdef MASKING
          hsn(i,Mm(ng),linew)=hsn(i,Mm(ng),linew)*                      &
     &                          GRID(ng)%rmask(i,Mm(ng))
#  endif
#  ifdef WET_DRY
          hsn(i,Mm(ng),linew)=hsn(i,Mm(ng),linew)*                      &
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
        hsn(0,0,linew)=0.5_r8*(hsn(1,0,linew)+                          &
     &                         hsn(0,1,linew))
      END IF
      IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
        hsn(Lm(ng)+1,0,linew)=0.5_r8*(hsn(Lm(ng)+1,1,linew)+            &
     &                                hsn(Lm(ng)  ,0,linew))
      END IF
      IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
        hsn(0,Mm(ng)+1,linew)=0.5_r8*(hsn(0,Mm(ng)  ,linew)+            &
     &                                hsn(1,Mm(ng)+1,linew))
      END IF
      IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        hsn(Lm(ng)+1,Mm(ng)+1,linew)=0.5_r8*                            &
     &             (hsn(Lm(ng)+1,Mm(ng)  ,linew)+                       &
     &              hsn(Lm(ng)  ,Mm(ng)+1,linew))
      END IF
#endif
      RETURN
      END SUBROUTINE hsnbc_tile
#endif

      END MODULE hsnbc_mod
