#include "cppdefs.h"
      MODULE akbc_mod
#if defined SOLVE3D && defined TKE_MIXING
!
!svn $Id$
!================================================== Robin Robertson ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets lateral boundary conditions for       !
!  vertical diffusivities associated                        !
!  with the Mellor and Yamada or Nakanishi-Nino closures.                        !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: akbc_tile

      CONTAINS
!
!***********************************************************************
      SUBROUTINE akbc (ng, tile, nout)
!***********************************************************************
!
      USE mod_param
      USE mod_mixing
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nout
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL akbc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 MIXING(ng)% Akv,                                 &
     &                 MIXING(ng)% Akt)
      RETURN
      END SUBROUTINE akbc
!
!***********************************************************************
      SUBROUTINE akbc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Akv, Akt)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
# else
      real(r8), intent(inout) :: Akv(LBi:UBi,LBj:UBj,0:N(ng))
      real(r8), intent(inout) :: Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT)
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, it

# include "set_bounds.h"

!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, gradient boundary condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend
            Akv(Istr-1,j,k)=Akv(Istr,j,k)
            DO it=1,NAT
              Akt(Istr-1,j,k,it)=Akt(Istr,j,k,it)
            END DO
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, gradient boundary condition.
!
        DO k=0,N(ng)
          DO j=Jstr,Jend
            Akv(Iend+1,j,k)=Akv(Iend,j,k)
            DO it=1,NAT
              Akt(Iend+1,j,k,it)=Akt(Iend,j,k,it)
            END DO
          END DO
        END DO
      END IF

!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, gradient boundary condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend
            Akv(i,Jstr-1,k)=Akv(i,Jstr,k)
            DO it=1,NAT
              Akt(i,Jstr-1,k,it)=Akt(i,Jstr,k,it)
            END DO
          END DO
        END DO
      END IF

!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, gradient boundary condition.
!
        DO k=0,N(ng)
          DO i=Istr,Iend
            Akv(i,Jend+1,k)=Akv(i,Jend,k)
            DO it=1,NAT
              Akt(i,Jend+1,k,it)=Akt(i,Jend,k,it)
            END DO
          END DO
        END DO
      END IF

!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            DO k=0,N(ng)
              Akv(Istr-1,Jstr-1,k)=0.5_r8*(Akv(Istr  ,Jstr-1,k)+  &
     &                                 Akv(Istr-1,Jstr  ,k))
              DO it=1,NAT
                Akt(Istr-1,Jstr-1,k,it)=                          &
     &                  0.5_r8*(Akt(Istr  ,Jstr-1,k,it)+          &
     &                          Akt(Istr-1,Jstr  ,k,it))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=0,N(ng)
              Akv(Iend+1,Jstr-1,k)=0.5_r8*(Akv(Iend  ,Jstr-1,k)+  &
     &                                     Akv(Iend+1,Jstr  ,k))
              DO it=1,NAT
                Akt(Iend+1,Jstr-1,k,it)=                          &
     &                  0.5_r8*(Akt(Iend  ,Jstr-1,k,it)+          &
     &                          Akt(Iend+1,Jstr  ,k,it))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=0,N(ng)
              Akv(Istr-1,Jend+1,k)=0.5_r8*(Akv(Istr  ,Jend+1,k)+  &
     &                                     Akv(Istr-1,Jend  ,k))
              DO it=1,NAT
                Akt(Istr-1,Jend+1,k,it)=                          &
     &                  0.5_r8*(Akt(Istr  ,Jend+1,k,it)+          &
     &                          Akt(Istr-1,Jend  ,k,it))
              END DO
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=0,N(ng)
              Akv(Iend+1,Jend+1,k)=0.5_r8*(Akv(Iend  ,Jend+1,k)+  &
     &                                     Akv(Iend+1,Jend  ,k))
              DO it=1,NAT
                Akt(Iend+1,Jend+1,k,it)=                          &
     &                  0.5_r8*(Akt(Iend  ,Jend+1,k,it)+          &
     &                          Akt(Iend+1,Jend  ,k,it))
              END DO
            END DO
          END IF
        END IF
      END IF

      RETURN
      END SUBROUTINE akbc_tile
#endif
      END MODULE akbc_mod
