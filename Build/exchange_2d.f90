      MODULE exchange_2d_mod
!
!svn $Id: exchange_2d.F 975 2009-05-05 22:51:13Z kate $
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This package contains periodic boundary conditions routines for 2D  !
!  variables.                                                          !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    exchange_p2d_tile    periodic conditions/exchange at PSI-points   !
!    exchange_r2d_tile    periodic conditions/exchange at RHO-points   !
!    exchange_u2d_tile    periodic conditions/exchange at U-points     !
!    exchange_v2d_tile    periodic conditions/exchange at V-points     !
!                                                                      !
!=======================================================================
!
      implicit none
      CONTAINS
!
!***********************************************************************
      SUBROUTINE exchange_p2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(r8), intent(inout) :: A(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrR, IstrT, IstrU, Iend, IendR, IendT
      integer :: Jstr, JstrR, JstrT, JstrV, Jend, JendR, JendT
!
      Istr =BOUNDS(ng)%Istr (tile)
      IstrR=BOUNDS(ng)%IstrR(tile)
      IstrT=BOUNDS(ng)%IstrT(tile)
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      IendR=BOUNDS(ng)%IendR(tile)
      IendT=BOUNDS(ng)%IendT(tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JstrT=BOUNDS(ng)%JstrT(tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
      JendR=BOUNDS(ng)%JendR(tile)
      JendT=BOUNDS(ng)%JendT(tile)
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Istr.eq.1) THEN
          DO j=Jstr,JendR
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
          END DO
        END IF
        IF (Iend.eq.Lm(ng)) THEN
          DO j=Jstr,JendR
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
      RETURN
      END SUBROUTINE exchange_p2d_tile
!
!***********************************************************************
      SUBROUTINE exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(r8), intent(inout) :: A(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrR, IstrT, IstrU, Iend, IendR, IendT
      integer :: Jstr, JstrR, JstrT, JstrV, Jend, JendR, JendT
!
      Istr =BOUNDS(ng)%Istr (tile)
      IstrR=BOUNDS(ng)%IstrR(tile)
      IstrT=BOUNDS(ng)%IstrT(tile)
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      IendR=BOUNDS(ng)%IendR(tile)
      IendT=BOUNDS(ng)%IendT(tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JstrT=BOUNDS(ng)%JstrT(tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
      JendR=BOUNDS(ng)%JendR(tile)
      JendT=BOUNDS(ng)%JendT(tile)
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Istr.eq.1) THEN
          DO j=JstrR,JendR
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
          END DO
        END IF
        IF (Iend.eq.Lm(ng)) THEN
          DO j=JstrR,JendR
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
      RETURN
      END SUBROUTINE exchange_r2d_tile
!
!***********************************************************************
      SUBROUTINE exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(r8), intent(inout) :: A(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrR, IstrT, IstrU, Iend, IendR, IendT
      integer :: Jstr, JstrR, JstrT, JstrV, Jend, JendR, JendT
!
      Istr =BOUNDS(ng)%Istr (tile)
      IstrR=BOUNDS(ng)%IstrR(tile)
      IstrT=BOUNDS(ng)%IstrT(tile)
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      IendR=BOUNDS(ng)%IendR(tile)
      IendT=BOUNDS(ng)%IendT(tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JstrT=BOUNDS(ng)%JstrT(tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
      JendR=BOUNDS(ng)%JendR(tile)
      JendT=BOUNDS(ng)%JendT(tile)
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Istr.eq.1) THEN
          DO j=JstrR,JendR
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
          END DO
        END IF
        IF (Iend.eq.Lm(ng)) THEN
          DO j=JstrR,JendR
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
      RETURN
      END SUBROUTINE exchange_u2d_tile
!
!***********************************************************************
      SUBROUTINE exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              A)
!***********************************************************************
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      real(r8), intent(inout) :: A(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrR, IstrT, IstrU, Iend, IendR, IendT
      integer :: Jstr, JstrR, JstrT, JstrV, Jend, JendR, JendT
!
      Istr =BOUNDS(ng)%Istr (tile)
      IstrR=BOUNDS(ng)%IstrR(tile)
      IstrT=BOUNDS(ng)%IstrT(tile)
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      IendR=BOUNDS(ng)%IendR(tile)
      IendT=BOUNDS(ng)%IendT(tile)
      Jstr =BOUNDS(ng)%Jstr (tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JstrT=BOUNDS(ng)%JstrT(tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
      JendR=BOUNDS(ng)%JendR(tile)
      JendT=BOUNDS(ng)%JendT(tile)
!
!-----------------------------------------------------------------------
!  East-West periodic boundary conditions.
!-----------------------------------------------------------------------
!
        IF (Istr.eq.1) THEN
          DO j=Jstr,JendR
            A(Lm(ng)+1,j)=A(1,j)
            A(Lm(ng)+2,j)=A(2,j)
          END DO
        END IF
        IF (Iend.eq.Lm(ng)) THEN
          DO j=Jstr,JendR
            A(-2,j)=A(Lm(ng)-2,j)
            A(-1,j)=A(Lm(ng)-1,j)
            A( 0,j)=A(Lm(ng)  ,j)
          END DO
        END IF
      RETURN
      END SUBROUTINE exchange_v2d_tile
      END MODULE exchange_2d_mod
