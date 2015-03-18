#include "cppdefs.h"
      MODULE mod_fishing
#if defined NEMURO_SAN && defined FISHING_FLEET
!
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!  Module for fishing fleet.                                           !
!=======================================================================
!
        USE mod_kinds
        USE mod_param
        USE mod_types
        USE mod_ocean
        USE mod_grid
        USE mod_fish
        USE mod_fleet
        USE mod_biology
        USE mod_scalars
        USE mod_stepping
        USE mod_parallel
        USE nrutil
        USE ran_state, ONLY: ran_seed
# ifdef DISTRIBUTE
          USE mod_param
          USE mod_ncparam, ONLY : r2dvar, r3dvar
# endif

        implicit none

        integer, allocatable :: fish_count_all(:,:)
        integer, allocatable :: nearport_all(:,:)
        integer, allocatable :: cpuecount_all(:,:)
        real(r8), allocatable :: xr_all(:,:)
        real(r8), allocatable :: yr_all(:,:)
        real(r8), allocatable :: cpue_all(:,:,:)
        real(r8), allocatable :: avgcpue_all(:,:)
        real(r8), allocatable :: percpue_all(:,:)
        real(r8), allocatable :: new_climcpue_all(:,:,:)
        real(r8), allocatable :: climcpue_all(:,:,:)
        real(r8), allocatable :: distport_all(:,:)

        PUBLIC  :: initports, initcpue, bestlocation, catch, updatecpue

        CONTAINS

        SUBROUTINE gather_array(ng, model, LBi, UBi, LBj, UBj,          &
     &                          xr, yr, fish_count, distport, nearport, &
     &                          cpue, avgcpue, percpue, new_climcpue,   &
     &                          climcpue)
          USE mod_kinds
# ifdef DISTRIBUTE
          USE distribute_mod, ONLY : mp_gather2d, mp_gather3d
# endif
          integer, intent(in) :: ng, model
          integer, intent(in) :: LBi, UBi, LBj, UBj
# ifdef ASSUMED_SHAPE
          integer, intent(in) :: fish_count(LBi:,LBj:)
          integer, intent(in) :: nearport(LBi:,LBj:)
          real(r8), intent(in) :: xr(LBi:,LBj:)
          real(r8), intent(in) :: yr(LBi:,LBj:)
          real(r8), intent(in) :: distport(LBi:,LBj:)
          real(r8), intent(in) :: cpue(LBi:,LBj:,:)
          real(r8), intent(in) :: avgcpue(LBi:,LBj:)
          real(r8), intent(in) :: percpue(LBi:,LBj:)
          real(r8), intent(in) :: new_climcpue(LBi:,LBj:,:)
          real(r8), intent(in) :: climcpue(LBi:,LBj:,:)
# else
          integer, intent(in) :: fish_count((LBi:UBi,LBj:UBj)
          integer, intent(in) :: nearport((LBi:UBi,LBj:UBj)
          real(r8), intent(in) :: xr((LBi:UBi,LBj:UBj)
          real(r8), intent(in) :: yr((LBi:UBi,LBj:UBj)
          real(r8), intent(in) :: distport((LBi:UBi,LBj:UBj)
          real(r8), intent(in) :: cpue((LBi:UBi,LBj:UBj,1:10)
          real(r8), intent(in) :: avgcpue((LBi:UBi,LBj:UBj)
          real(r8), intent(in) :: percpue((LBi:UBi,LBj:UBj)
          real(r8), intent(in) :: new_climcpue((LBi:UBi,LBj:UBj,1:12)
          real(r8), intent(in) :: climcpue((LBi:UBi,LBj:UBj,1:12)
# endif
          integer :: i, j, l
# ifdef DISTRIBUTE
          integer :: Npts
          real(r8), dimension(LBi:UBi,LBj:UBj) :: r_fshcnt
          real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: r_fshcnt_all
          real(r8), dimension(LBi:UBi,LBj:UBj) :: r_nrport
          real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: r_nrport_all
          real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)) :: line2d
          real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)*10) :: line3d_10
          real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)*12) :: line3d_12
# endif
!
# ifdef DISTRIBUTE
          line2d = 0.0_r8
          line3d_10 = 0.0_r8
          line3d_12 = 0.0_r8
# endif
          IF (.not. ALLOCATED(fish_count_all)) THEN
            allocate(fish_count_all(0:(Lm(ng)+1),0:(Mm(ng)+1)))
          END IF
          IF (.not. ALLOCATED(xr_all)) THEN
            allocate(xr_all(0:(Lm(ng)+1),0:(Mm(ng)+1)))
          END IF
          IF (.not. ALLOCATED(yr_all)) THEN
            allocate(yr_all(0:(Lm(ng)+1),0:(Mm(ng)+1)))
          END IF
          IF (.not. ALLOCATED(cpuecount_all)) THEN
            allocate(cpuecount_all(0:(Lm(ng)+1),0:(Mm(ng)+1)))
          END IF
          IF (.not. ALLOCATED(nearport_all)) THEN
            allocate(nearport_all(0:(Lm(ng)+1),0:(Mm(ng)+1)))
          END IF
          IF (.not. ALLOCATED(distport_all)) THEN
            allocate(distport_all(0:(Lm(ng)+1),0:(Mm(ng)+1)))
          END IF
          IF (.not. ALLOCATED(cpue_all)) THEN
            allocate(cpue_all(0:(Lm(ng)+1),0:(Mm(ng)+1),1:10))
          END IF
          IF (.not. ALLOCATED(avgcpue_all)) THEN
            allocate(avgcpue_all(0:(Lm(ng)+1),0:(Mm(ng)+1)))
          END IF
          IF (.not. ALLOCATED(percpue_all)) THEN
            allocate(percpue_all(0:(Lm(ng)+1),0:(Mm(ng)+1)))
          END IF
          IF (.not. ALLOCATED(new_climcpue_all)) THEN
            allocate(new_climcpue_all(0:(Lm(ng)+1),0:(Mm(ng)+1),1:12))
          END IF
          IF (.not. ALLOCATED(climcpue_all)) THEN
            allocate(climcpue_all(0:(Lm(ng)+1),0:(Mm(ng)+1),1:12))
          END IF
! Gather 2D integer variables
# ifdef DISTRIBUTE
          Npts=IOBOUNDS(ng)%xi_rho*IOBOUNDS(ng)%eta_rho
          DO i=LBi,UBi
            DO j=LBj,UBj
              r_fshcnt(i,j)=REAL(fish_count(i,j),r8)
              r_nrport(i,j)=REAL(nearport(i,j),r8)
            END DO
          END DO
          CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,              &
     &                      1, r2dvar, 1.0_r8,                          &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
     &                      r_fshcnt, Npts, line2d, .false.)
          r_fshcnt_all = reshape(line2d,(/ Lm(ng)+2, Mm(ng)+2 /))
          CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,              &
     &                      1, r2dvar, 1.0_r8,                          &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
     &                      r_nrport, Npts, line2d, .false.)
          r_nrport_all = reshape(line2d,(/ Lm(ng)+2, Mm(ng)+2 /))
          DO i=0,Lm(ng)+1
            DO j=0,Mm(ng)+1
              fish_count_all(i,j)=INT(r_fshcnt_all(i,j))
              nearport_all(i,j)=INT(r_nrport_all(i,j))
! reset cpuecount each time
              cpuecount_all(i,j)=0
            END DO
          END DO
! Gather 2D real variables
          CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,              &
     &                      1, r2dvar, 1.0_r8,                          &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
     &                      xr, Npts, line2d, .false.)
          xr_all = reshape(line2d,(/ Lm(ng)+2, Mm(ng)+2 /))
          CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,              &
     &                      1, r2dvar, 1.0_r8,                          &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
     &                      yr, Npts, line2d, .false.)
          yr_all = reshape(line2d,(/ Lm(ng)+2, Mm(ng)+2 /))
          CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,              &
     &                      1, r2dvar, 1.0_r8,                          &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
     &                      distport, Npts, line2d, .false.)
          distport_all = reshape(line2d,(/ Lm(ng)+2, Mm(ng)+2 /))
          CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,              &
     &                      1, r2dvar, 1.0_r8,                          &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
     &                      avgcpue, Npts, line2d, .false.)
          avgcpue_all = reshape(line2d,(/ Lm(ng)+2, Mm(ng)+2 /))
          CALL mp_gather2d (ng, model, LBi, UBi, LBj, UBj,              &
     &                      1, r2dvar, 1.0_r8,                          &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
     &                      percpue, Npts, line2d, .false.)
          percpue_all = reshape(line2d,(/ Lm(ng)+2, Mm(ng)+2 /))
! Gather 3D real variables
          Npts=IOBOUNDS(ng)%xi_rho*IOBOUNDS(ng)%eta_rho*10
          CALL mp_gather3d (ng, model, LBi, UBi, LBj, UBj,              &
     &                      1, 10, 1, r3dvar, 1.0_r8,                   &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
     &                      cpue, Npts, line3d_10, .false.)
          cpue_all = reshape(line3d_10,(/ Lm(ng)+2, Mm(ng)+2, 10 /))
          Npts=IOBOUNDS(ng)%xi_rho*IOBOUNDS(ng)%eta_rho*12
          CALL mp_gather3d (ng, model, LBi, UBi, LBj, UBj,              &
     &                      1, 12, 1, r3dvar, 1.0_r8,                   &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
     &                      new_climcpue, Npts, line3d_12, .false.)
          new_climcpue_all = reshape(line3d_12,                         &
     &                           (/ Lm(ng)+2, Mm(ng)+2, 12 /))
          CALL mp_gather3d (ng, model, LBi, UBi, LBj, UBj,              &
     &                      1, 12, 1, r3dvar, 1.0_r8,                   &
#  ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#  endif
     &                      climcpue, Npts, line3d_12, .false.)
          climcpue_all = reshape(line3d_12,                             &
     &                           (/ Lm(ng)+2, Mm(ng)+2, 12 /))
# else
          DO i=0,Lm(ng)+1
            DO j=0,Mm(ng)+1
              nearport_all(i,j)=nearport(i,j)
              xr_all(i,j)=xr(i,j)
              yr_all(i,j)=yr(i,j)
              distport_all(i,j)=distport(i,j)
              avgcpue_all(i,j)=avgcpue(i,j)
              percpue_all(i,j)=percpue(i,j)
! reset cpuecount each time
              cpuecount_all(i,j)=0
              DO l=1,10
                cpue_all(i,j,l)=cpue(i,j,l)
              END DO
              DO l=1,12
                new_climcpue_all(i,j,l)=new_climcpue(i,j,l)
                climcpue_all(i,j,l)=climcpue(i,j,l)
              END DO
            END DO
          END DO
# endif
!        print*, 'END gather'
        END SUBROUTINE gather_array
!
        SUBROUTINE scatter_array(ng, model, LBi, UBi, LBj, UBj,         &
     &                          distport, nearport,                     &
     &                          cpue, avgcpue, percpue, new_climcpue,   &
     &                          climcpue)
# ifdef DISTRIBUTE
          USE distribute_mod, ONLY : mp_scatter2d, mp_scatter3d
# endif
          integer, intent(in) :: ng, model
          integer, intent(in) :: LBi, UBi, LBj, UBj
# ifdef ASSUMED_SHAPE
          integer, intent(out) :: nearport(LBi:,LBj:)
          real(r8), intent(out) :: distport(LBi:,LBj:)
          real(r8), intent(out) :: cpue(LBi:,LBj:,:)
          real(r8), intent(out) :: avgcpue(LBi:,LBj:)
          real(r8), intent(out) :: percpue(LBi:,LBj:)
          real(r8), intent(out) :: new_climcpue(LBi:,LBj:,:)
          real(r8), intent(out) :: climcpue(LBi:,LBj:,:)
# else
          integer, intent(out) :: nearport((LBi:UBi,LBj:UBj)
          real(r8), intent(out) :: distport((LBi:UBi,LBj:UBj)
          real(r8), intent(out) :: cpue((LBi:UBi,LBj:UBj,1:10)
          real(r8), intent(out) :: avgcpue((LBi:UBi,LBj:UBj)
          real(r8), intent(out) :: percpue((LBi:UBi,LBj:UBj)
          real(r8), intent(out) :: new_climcpue((LBi:UBi,LBj:UBj,1:12)
          real(r8), intent(out) :: climcpue((LBi:UBi,LBj:UBj,1:12)
# endif
          integer :: i, j, l, np
# ifdef DISTRIBUTE
          integer :: Npts
          real(r8) :: Vmin, Vmax
          real(r8), dimension(LBi:UBi,LBj:UBj) :: r_nrport
          real(r8), dimension(0:(Lm(ng)+1),0:(Mm(ng)+1)) :: r_nrport_all
          real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)+2) :: line2d
          real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)*10+2) :: line3d_10
          real(r8), dimension((Lm(ng)+2)*(Mm(ng)+2)*12+2) :: line3d_12
# endif
!
# ifdef DISTRIBUTE
! Gather 2D integer variables
          Npts=(Lm(ng)+2)*(Mm(ng)+2)+2
          DO i=0,Lm(ng)+1
            DO j=0,Mm(ng)+1
              r_nrport_all(i,j)=REAL(nearport_all(i,j),r8)
            END DO
          END DO
          line2d = 0.0_r8
          np=0
          DO j=0,Mm(ng)+1
            DO i=0,Lm(ng)+1
              np=np+1
              line2d(np)=r_nrport_all(i,j)
            END DO
          END DO
          IF (Master) THEN
            Vmin=line2d(1)
            Vmax=line2d(1)
            DO np=1,Npts-2
              Vmin=MIN(Vmin,line2d(np))
              Vmax=MAX(Vmax,line2d(np))
            END DO
!            print*, 'SCAT,NPRT', Vmin,Vmax
          END IF
          CALL mp_scatter2d (ng, model, LBi, UBi, LBj, UBj,             &
     &                       NghostPoints, r2dvar, Vmin, Vmax,          &
     &                       Npts, line2d, r_nrport)
          DO i=LBi,UBi
            DO j=LBj,UBj
              nearport(i,j)=INT(r_nrport(i,j))
            END DO
          END DO
! Gather 2D real variables
          line2d = 0.0_r8
          np=0
          DO j=0,Mm(ng)+1
            DO i=0,Lm(ng)+1
              np=np+1
              line2d(np)=distport_all(i,j)
            END DO
          END DO
          IF (Master) THEN
            Vmin=line2d(1)
            Vmax=line2d(1)
            DO np=1,Npts-2
              Vmin=MIN(Vmin,line2d(np))
              Vmax=MAX(Vmax,line2d(np))
            END DO
!            print*, 'SCAT,DIST', Vmin,Vmax
          END IF
          CALL mp_scatter2d (ng, model, LBi, UBi, LBj, UBj,             &
     &                       NghostPoints, r2dvar, Vmin, Vmax,          &
     &                       Npts, line2d, distport)
          line2d = 0.0_r8
          np=0
          DO j=0,Mm(ng)+1
            DO i=0,Lm(ng)+1
              np=np+1
              line2d(np)=avgcpue_all(i,j)
            END DO
          END DO
          IF (Master) THEN
            Vmin=line2d(1)
            Vmax=line2d(1)
            DO np=1,Npts-2
              Vmin=MIN(Vmin,line2d(np))
              Vmax=MAX(Vmax,line2d(np))
            END DO
!            print*, 'SCAT,ACPUE', Vmin,Vmax
          END IF
          CALL mp_scatter2d (ng, model, LBi, UBi, LBj, UBj,             &
     &                       NghostPoints, r2dvar, Vmin, Vmax,          &
     &                       Npts, line2d, avgcpue)
          line2d = 0.0_r8
          np=0
          DO j=0,Mm(ng)+1
            DO i=0,Lm(ng)+1
              np=np+1
              line2d(np)=percpue_all(i,j)
            END DO
          END DO
          IF (Master) THEN
            Vmin=line2d(1)
            Vmax=line2d(1)
            DO np=1,Npts-2
              Vmin=MIN(Vmin,line2d(np))
              Vmax=MAX(Vmax,line2d(np))
            END DO
!            print*, 'SCAT,PCPUE', Vmin,Vmax
          END IF
          CALL mp_scatter2d (ng, model, LBi, UBi, LBj, UBj,             &
     &                       NghostPoints, r2dvar, Vmin, Vmax,          &
     &                       Npts, line2d, percpue)
! Gather 3D real variables
          Npts=(Lm(ng)+2)*(Mm(ng)+2)*10+2
          line3d_10 = 0.0_r8
          np=0
          DO l=1,10
            DO j=0,Mm(ng)+1
              DO i=0,Lm(ng)+1
                np=np+1
                line3d_10(np)=cpue_all(i,j,l)
              END DO
            END DO
          END DO
          IF (Master) THEN
            Vmin=line3d_10(1)
            Vmax=line3d_10(1)
            DO np=1,Npts-2
              Vmin=MIN(Vmin,line3d_10(np))
              Vmax=MAX(Vmax,line3d_10(np))
            END DO
!            print*, 'SCAT,CPUE', Vmin,Vmax
          END IF
          CALL mp_scatter3d (ng, model, LBi, UBi, LBj, UBj, 1, 10,      &
     &                       NghostPoints, r3dvar, Vmin, Vmax,          &
     &                       Npts, line3d_10, cpue)
          Npts=(Lm(ng)+2)*(Mm(ng)+2)*12+2
          line3d_12 = 0.0_r8
          np=0
          DO l=1,12
            DO j=0,Mm(ng)+1
              DO i=0,Lm(ng)+1
                np=np+1
                line3d_12(np)=new_climcpue_all(i,j,l)
              END DO
            END DO
          END DO
          IF (Master) THEN
            Vmin=line3d_12(1)
            Vmax=line3d_12(1)
            DO np=1,Npts-2
              Vmin=MIN(Vmin,line3d_12(np))
              Vmax=MAX(Vmax,line3d_12(np))
            END DO
!            print*, 'SCAT,YCPUE', Vmin,Vmax
          END IF
          CALL mp_scatter3d (ng, model, LBi, UBi, LBj, UBj, 1, 12,      &
     &                       NghostPoints, r3dvar, Vmin, Vmax,          &
     &                       Npts, line3d_12, new_climcpue)
          line3d_12 = 0.0_r8
          np=0
          DO l=1,12
            DO j=0,Mm(ng)+1
              DO i=0,Lm(ng)+1
                np=np+1
                line3d_12(np)=climcpue_all(i,j,l)
              END DO
            END DO
          END DO
          IF (Master) THEN
            Vmin=line3d_12(1)
            Vmax=line3d_12(1)
            DO np=1,Npts-2
              Vmin=MIN(Vmin,line3d_12(np))
              Vmax=MAX(Vmax,line3d_12(np))
            END DO
!            print*, 'SCAT,CCPUE', Vmin,Vmax
          END IF
          CALL mp_scatter3d (ng, model, LBi, UBi, LBj, UBj, 1, 12,      &
     &                       NghostPoints, r3dvar, Vmin, Vmax,          &
     &                       Npts, line3d_12, climcpue)
# else
          DO i=0,Lm(ng)+1
            DO j=0,Mm(ng)+1
              nearport(i,j)=nearport_all(i,j)
!             xr(i,j)=xr_all(i,j)
!             yr(i,j)=yr_all(i,j)
              distport(i,j)=distport_all(i,j)
              avgcpue(i,j)=avgcpue_all(i,j)
              percpue(i,j)=percpue_all(i,j)
              DO l=1,10
                cpue(i,j,l)=cpue_all(i,j,l)
              END DO
              DO l=1,12
                new_climcpue(i,j,l)=new_climcpue_all(i,j,l)
                climcpue(i,j,l)=climcpue_all(i,j,l)
              END DO
            END DO
          END DO
# endif
          IF (Master) THEN
            IF (ALLOCATED(fish_count_all)) THEN
              deallocate(fish_count_all)
            END IF
            IF (ALLOCATED(xr_all)) THEN
              deallocate(xr_all)
            END IF
            IF (ALLOCATED(yr_all)) THEN
              deallocate(yr_all)
            END IF
            IF (ALLOCATED(cpuecount_all)) THEN
              deallocate(cpuecount_all)
            END IF
            IF (ALLOCATED(nearport_all)) THEN
              deallocate(nearport_all)
            END IF
            IF (ALLOCATED(distport_all)) THEN
              deallocate(distport_all)
            END IF
            IF (ALLOCATED(cpue_all)) THEN
              deallocate(cpue_all)
            END IF
            IF (ALLOCATED(avgcpue_all)) THEN
              deallocate(avgcpue_all)
            END IF
            IF (ALLOCATED(percpue_all)) THEN
              deallocate(percpue_all)
            END IF
            IF (ALLOCATED(new_climcpue_all)) THEN
              deallocate(new_climcpue_all)
            END IF
            IF (ALLOCATED(climcpue_all)) THEN
              deallocate(climcpue_all)
            END IF
          END IF
!        print*, 'END scatter'
        END SUBROUTINE scatter_array
!
        SUBROUTINE initports(ng)
          integer, intent(in) :: ng
          integer :: i, j, ip, ipc, jpc
          real(r8) :: dist2port, dxm, dym, dist
!
! Find closest port and distance from a given i,j location on grid
          DO i=0,Lm(ng)+1
            DO j=0,Mm(ng)+1
              dist2port=1.0e30_r8
              DO ip=1,Nports(ng)
                ipc=iPort(ip,ng)
                jpc=jPort(ip,ng)
                dxm=xr_all(i,j)-xr_all(ipc,jpc)
                dym=yr_all(i,j)-yr_all(ipc,jpc)
                dist=0.001_r8*(dxm**2+dym**2)**0.5_r8
                IF (dist.lt.dist2port) THEN
                  dist2port=dist
                  distport_all(i,j)=dist
                  nearport_all(i,j)=ip
                END IF
              END DO
!              print*, 'INIP', i, j, distport_all(i,j), nearport_all(i,j)
            END DO
          END DO
!        print*, 'END iniport'
        END SUBROUTINE initports
!
        SUBROUTINE initcpue(ng)
          integer, intent(in) :: ng
          integer :: i, j, l, ijcell, ifid, mo, l1, l2
          logical :: foundsard
          real(r8) :: sumcpue, daymo
!
! initialize CPUE values based on fish_by_cell information
          DO i=0,Lm(ng)+1
            DO j=0,Mm(ng)+1
              sumcpue=0.0_r8
              IF (fish_count_all(i,j).gt.0) THEN
!                foundsard=.TRUE.
                foundsard=.FALSE.
                ijcell=i+(Lm(ng)+2)*j
                DO ifid=1,Nfish(ng)
                  IF ((FISHES(ng)%species(ifid).eq.if_sardine).and.     &
     &                (FISHES(ng)%alive(ifid)).and.                     &
     &                (FISHES(ng)%lifestage(ifid).ge.if_subadult)) THEN
                    IF (FISHES(ng)%cellid(ifid).eq.ijcell) THEN
                      foundsard=.TRUE.
                    END IF
                  END IF
                END DO
! CPUE for last 10 days
                DO l=1,10
                  IF (foundsard) THEN
                    cpue_all(i,j,l)=CatchMax(ng)
                  ELSE
                    cpue_all(i,j,l)=0.0_r8
                  END IF
                  sumcpue=sumcpue+cpue_all(i,j,l)
                END DO
              END IF
! 10-day average CPUE (actual and perceived)
              avgcpue_all(i,j)=sumcpue/10.0_r8
              percpue_all(i,j)=avgcpue_all(i,j)
! CPUE climatology based on averaging period of Ndays
              daymo=days_year/12.0_r8
              DO mo=1,12
                l1=INT(REAL(mo-1,r8)*daymo)+1
                l2=INT(REAL(mo,r8)*daymo)
                sumcpue=0.0_r8
                DO l=l1,l2
                  sumcpue=sumcpue+avgcpue_all(i,j)
                END DO
                climcpue_all(i,j,mo)=sumcpue/REAL(l2-l1+1,r8)
              END DO
!              IF (avgcpue_all(i,j).gt.0.0_r8) THEN
!                print*, 'INIC', i, j, avgcpue_all(i,j)
!              END IF
            END DO
          END DO
!        print*, 'END inicpue'
        END SUBROUTINE initcpue
!
        SUBROUTINE bestlocation(ng,ib)
          integer, intent(in) :: ng, ib
          integer :: i, j, ip, idport, ibc, jbc
          real(r8) :: dxm, dym, dist
          real(r8) :: dist2loc, time2loc, dist2port, time2port
          real(r8) :: time2fish, exprev, t1
          real(r8) :: exprevOpt, ev1mu, ev1beta, snudg
!
          ev1mu=0.0_r8
          ev1beta=500.0_r8
!
!          do i=0,Lm(ng)+1
!          do j=0,Mm(ng)+1
!            print*, 'CPUE', i, j, avgcpue_all(i,j) 
!          end do
!          end do
          BOATS(ng)%blocflag(ib)=0
          BOATS(ng)%t2LocOpt(ib)=0.0_r8
          BOATS(ng)%t2PortOpt(ib)=0.0_r8
          exprevOpt=-99.0_r8
          ibc=BOATS(ng)%boat(ibiloc,ib)
          jbc=BOATS(ng)%boat(ibjloc,ib)
          DO i=0,Lm(ng)+1
            DO j=0,Mm(ng)+1
              dxm=xr_all(i,j)-xr_all(ibc,jbc)
              dym=yr_all(i,j)-yr_all(ibc,jbc)
              dist2loc=0.001_r8*(dxm**2+dym**2)**0.5_r8
              time2loc=dist2loc/BoatVel(ng)
              dist2port=distport_all(i,j)
              idport=nearport_all(i,j)
              time2port=dist2port/BoatVel(ng)
              time2fish=BOATS(ng)%tAvail(ib)-time2loc-time2port
              IF (time2fish.gt.FishTime(ng)) THEN
! NOTE: percpue from cpue routine in KR code
                exprev=percpue_all(i,j)*CatchPrice(idport,ng)-          &
     &                       (time2loc+time2port)*TravCost(ng)
                CALL ran1 (snudg)
                t1=ev1mu-ev1beta*LOG(-LOG(snudg))
                exprev=exprev+t1
              ELSE
                exprev=-99.0_r8
              END IF
              IF ((exprev.gt.0.0_r8).and.(exprev.gt.exprevOpt)) THEN
                exprevOpt=exprev
                BOATS(ng)%t2LocOpt(ib)=time2loc
                BOATS(ng)%t2PortOpt(ib)=time2port
                BOATS(ng)%boat(ibiloc,ib)=i
                BOATS(ng)%boat(ibjloc,ib)=j
                BOATS(ng)%boat(ibport,ib)=idport
                BOATS(ng)%blocflag(ib)=1
              END IF
!              IF (BOATS(ng)%blocflag(ib).ne.0) THEN
!                print*, 'BLOC1', ib, BOATS(ng)%boat(ibiloc,ib),         &
!     &                               BOATS(ng)%boat(ibjloc,ib)
!              ELSE
!                print*, 'BLOC0', BOATS(ng)%blocflag(ib)
!              END IF
            END DO
          END DO
!        print*, 'END bestloc'
        END SUBROUTINE bestlocation
!
        SUBROUTINE catch(ng, ib)
          integer, intent(in) :: ng, ib
          integer :: ibc, jbc, ifid, ienc, isMax, ijcell
          integer :: nsard, nsardenc
          type(fishnode), pointer :: thisfish
          logical :: fishing
          real(r8) :: snudg, time2fish, overcatch
          real(r8) :: sbio, sbioMax, Fworth, fac, Fmort
!
          BOATS(ng)%bfishflag(ib)=0
          BOATS(ng)%tAvail(ib)=BOATS(ng)%tAvail(ib)-                    &
     &                  BOATS(ng)%t2LocOpt(ib)-BOATS(ng)%t2PortOpt(ib)
          ibc=BOATS(ng)%boat(ibiloc,ib)
          jbc=BOATS(ng)%boat(ibjloc,ib)
          ijcell=ibc+(Lm(ng)+2)*jbc
!          print*, 'CATCH', ibc, jbc, BOATS(ng)%tAvail(ib)
          IF (fish_count_all(ibc,jbc).gt.0) THEN
            fishing=.TRUE.
            nsard=0
            nsardenc=0
            DO ienc=1,EncMax(ng)
              IF (fishing) THEN
! Find sardine with largest biomass
                sbioMax=0.0_r8
                isMax=0.0_r8
                DO ifid=1,Nfish(ng)
                  IF ((FISHES(ng)%species(ifid).eq.if_sardine).and.     &
     &                (FISHES(ng)%alive(ifid)).and.                     &
     &                (FISHES(ng)%lifestage(ifid).ge.if_subadult).and.  &
                      (FISHES(ng)%cellid(ifid).eq.ijcell)) THEN
! Impose encounter rate
                    nsard=nsard+1
                    CALL ran1 (snudg)
                    IF (snudg.lt.EncRate(ng)) THEN
                      nsardenc=nsardenc+1
                      sbio=0.001_r8*FISHES(ng)%bioenergy(ifwwt,ifid)*   &
     &                              FISHES(ng)%bioenergy(ifworth,ifid)
                      IF (sbio.gt.sbioMax) THEN
                        sbioMax=sbio
                        isMax=ifid
                      END IF
                    END IF
                  END IF
                END DO
                IF (sbioMax.gt.0.0_r8) THEN
                  Fworth=FISHES(ng)%bioenergy(ifworth,isMax)
! Found at least one sardine, so adjust fishing time (only once)
                  IF (nsardenc.eq.1) BOATS(ng)%tAvail(ib)=              &
                          BOATS(ng)%tAvail(ib)-FishTime(ng)
! Fish sardine with largest biomass
                  BOATS(ng)%catch(ib)=BOATS(ng)%catch(ib)+              &
     &                                  Qcatch(ng)*sbioMax
                  FISHES(ng)%bioenergy(ifworth,isMax)=                  &
     &                                         (1.0_r8-Qcatch(ng))*     &
     &                          FISHES(ng)%bioenergy(ifworth,isMax)
                  IF (FISHES(ng)%bioenergy(ifworth,isMax).lt.           &
     &                                           1.0_r8) THEN
                    FISHES(ng)%alive(isMax)=.FALSE.
                    IF (FISHES(ng)%deathby(isMax).eq.0)                 &
     &                             FISHES(ng)%deathby(isMax)=4
                  END IF
! Stop fishing if full or out of time
                  IF ((BOATS(ng)%catch(ib).gt.CatchMax(ng)).or.         &
     &                (BOATS(ng)%tAvail(ib).le.0.0_r8)) THEN
                    fishing=.FALSE.
                    BOATS(ng)%bfishflag(ib)=2
! Put overcatch back into fish
                    overcatch=BOATS(ng)%catch(ib)-CatchMax(ng)
                    IF (overcatch.gt.0.0_r8) THEN
                      BOATS(ng)%catch(ib)=CatchMax(ng)
                      IF (FISHES(ng)%bioenergy(ifwwt,isMax).gt.         &
     &                       0.01_r8*Fwwt0(if_sardine,ng)) THEN
                        FISHES(ng)%bioenergy(ifworth,isMax)= &
     &                  FISHES(ng)%bioenergy(ifworth,isMax)+overcatch*  &
     &                     1000.0_r8/FISHES(ng)%bioenergy(ifwwt,isMax)
                      END IF
                      IF (FISHES(ng)%bioenergy(ifworth,isMax).ge.       &
     &                                               1.0_r8) THEN
                        FISHES(ng)%alive(isMax)=.TRUE.
                        FISHES(ng)%deathby(isMax)=0
                      END IF
                    END IF
                  END IF
                  fac=FISHES(ng)%bioenergy(ifworth,isMax)/Fworth
                  Fmort=-LOG(MAX(fac,0.01_r8))
                  FISHES(ng)%fmortF(isMax)=FISHES(ng)%fmortF(isMax)+    &
     &                                                         Fmort
                END IF
              END IF
            END DO
          END IF
! Set boat flags
          IF (fishing) BOATS(ng)%bfishflag(ib)=3
          IF (nsard.eq.0) THEN
            BOATS(ng)%bfishflag(ib)=0
          ELSE
            IF (nsardenc.eq.0) BOATS(ng)%bfishflag(ib)=1
          END IF
          cpuecount_all(ibc,jbc)=cpuecount_all(ibc,jbc)+1
          cpue_all(ibc,jbc,1)=cpue_all(ibc,jbc,1)+BOATS(ng)%catch(ib)
!        print*, 'END catch'
        END SUBROUTINE catch
!
        SUBROUTINE updatecpue(ng)
          integer, intent(in) :: ng
          integer :: i, j, l, iday, imo, mo, l1, l2
          real(r8) :: nyears, sumcpue, percpue1, daymo
!
! Shift CPUE
          DO i=0,Lm(ng)+1
            DO j=0,Mm(ng)+1
              DO l=9,1,-1
                cpue_all(i,j,l+1)=cpue_all(i,j,l)
              END DO
              cpue_all(i,j,1)=0.0_r8
            END DO
          END DO

! Update CPUE based on today's catch
          nyears=REAL(INT(time(ng)/86400.0_r8/days_year),r8)
          iday=INT(time(ng)/86400.0_r8-days_year*nyears)+1
          iday=MIN(MAX(iday,1),INT(days_year))
          daymo=days_year/12.0_r8
          DO mo=1,12
            l1=INT(REAL(mo-1,r8)*daymo)+1
            l2=INT(REAL(mo,r8)*daymo)
            IF ((iday.ge.l1).and.(iday.le.l2)) imo=mo
          END DO
          DO i=0,Lm(ng)+1
            DO j=0,Mm(ng)+1
              IF (cpuecount_all(i,j).ge.1) THEN
                cpue_all(i,j,1)=cpue_all(i,j,1)/REAL(cpuecount_all(i,j))
                sumcpue=0.0_r8
                DO l=1,10
                  sumcpue=sumcpue+cpue_all(i,j,l)
                END DO
                percpue_all(i,j)=sumcpue/10.0_r8
              ELSE
                cpue_all(i,j,1)=0.0_r8
                percpue1=climcpue_all(i,j,imo)
                DO l=2,10
                  IF (cpue_all(i,j,l).gt.percpue1)                      &
     &              percpue1=cpue_all(i,j,l)
                END DO
                sumcpue=0.0_r8
                DO l=1,10
                  IF (l.eq.1) THEN
                    sumcpue=sumcpue+percpue1
                  ELSE
                    sumcpue=sumcpue+cpue_all(i,j,l)
                  END IF
                END DO
                percpue_all(i,j)=sumcpue/10.0_r8
              END IF
              sumcpue=0.0_r8
              DO l=1,10
                sumcpue=sumcpue+cpue_all(i,j,l)
              END DO
              avgcpue_all(i,j)=sumcpue/10.0_r8
              new_climcpue_all(i,j,imo)=new_climcpue_all(i,j,imo)+      &
     &                                             cpue_all(i,j,1)
            END DO
          END DO
!        print*, 'END updcpue'
        END SUBROUTINE updatecpue
!
        SUBROUTINE updatecpue_clim(ng)
          integer, intent(in) :: ng
          integer :: i, j, mo, l1, l2
          real(r8) :: sumcpue, daymo
!
! Yearly update of climatological CPUE
          daymo=days_year/12.0_r8
          DO i=0,Lm(ng)+1
            DO j=0,Mm(ng)+1
              DO mo=1,12
                l1=INT(REAL(mo-1,r8)*daymo)+1
                l2=INT(REAL(mo,r8)*daymo)
                climcpue_all(i,j,mo)=new_climcpue_all(i,j,mo)/          &
     &                                        REAL(l2-l1+1,r8)
                new_climcpue_all(i,j,mo)=0.0_r8
              END DO
            END DO
          END DO
!        print*, 'END updcpue_clim'
        END SUBROUTINE updatecpue_clim
#endif
      END MODULE mod_fishing
