#include "cppdefs.h"
      MODULE new_fish_mod
#if defined NEMURO_SAN
!
!svn $Id: step_floats.F 997 2009-06-04 17:42:20Z kate $
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Check to see if any fish spawned and make new superindividuals      !
!  from them.
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: new_fish

      CONTAINS
!
!***********************************************************************
      SUBROUTINE new_fish
!***********************************************************************
!
      USE mod_param
      USE mod_fish
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend
      integer :: i
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10)
# endif

      IF (Master) THEN
        call tree_init

        DO i=1,Nfish(ng)
          IF (bioenergy(i, ifeggs) .ne. 0) THEN
            CALL tree_insert(caviar(i), dist(i), i)
        END DO

        call tree_traverse
        CALL tree_collect(Nsuper, Nfound, eggs, moms)
        print *, 'after collect: ', Nfound
        do i=1,Nfound
          print *, i, eggs(i), moms(i)
        end do
        CALL tree_destroy
      END IF

      CALL new_fish_tile (ng, Lstr, Lend,                            &
     &                       knew(ng), nnew(ng), nfm3(ng), nfm2(ng),    &
     &                       nfm1(ng), nf(ng), nfp1(ng),                &
     &                       FLT(ng) % bounded,                         &
     &                       FLT(ng) % Ftype,                           &
     &                       FLT(ng) % Tinfo,                           &
     &                       FLT(ng) % Fz0,                             &
     &                       FLT(ng) % track)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10)
# endif
      RETURN
      END SUBROUTINE step_floats
!
!***********************************************************************
      SUBROUTINE step_floats_tile (ng, Lstr, Lend,                      &
     &                             knew, nnew,                          &
     &                             nfm3, nfm2, nfm1, nf, nfp1,          &
     &                             bounded, Ftype, Tinfo, Fz0, track)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
# ifdef FLOAT_VWALK
      USE mod_mixing
# endif
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
      USE interp_floats_mod
# ifdef DIAPAUSE
      USE interp_floats_diapW_mod
# endif
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_collect
# endif
      USE utility_mod, ONLY : nrng
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, Lstr, Lend
      integer, intent(in) :: knew, nnew, nfm3, nfm2, nfm1, nf, nfp1
!
# ifdef ASSUMED_SHAPE
      integer, intent(in) :: Ftype(:)
      real(r8), intent(inout) :: Tinfo(0:,:)
      real(r8), intent(inout) :: Fz0(:)

      logical, intent(inout) :: bounded(:)
      real(r8), intent(inout) :: track(:,0:,:)
# else
      integer, intent(in) :: Ftype(Nfloats(ng))
      real(r8), intent(inout) :: Tinfo(0:izrhs,Nfloats(ng))
      real(r8), intent(inout) :: Fz0(Nfloats(ng))

      logical, intent(inout) :: bounded(Nfloats(ng))
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfloats(ng))
# endif
!
!  Local variable declarations.
!
      logical, parameter :: Gmask = .FALSE.
# ifdef MASKING
      logical, parameter :: Lmask = .TRUE.
# else
      logical, parameter :: Lmask = .FALSE.
# endif
      logical, dimension(Lstr:Lend) :: MyThread

      integer :: LBi, UBi, LBj, UBj
      integer :: Ir, Jr, Npts, i, i1, i2, j, j1, j2, itrc, l, k

      real(r8), parameter :: Fspv = 0.0_r8

      real(r8) :: cff1, cff2, cff3, cff4, cff5, cff6, cff7, cff8, cff9
      real(r8) :: p1, p2, q1, q2, xrhs, yrhs, zrhs, zfloat

# ifdef FLOAT_VWALK
      integer :: ierr
      integer :: iseed = 149876
# endif

      real(r8), dimension(Lstr:Lend) :: nudg

# ifdef DISTRIBUTE
      real(r8) :: Xstr, Xend, Ystr, Yend
      real(r8), dimension(Nfloats(ng)*NFV(ng)*(NFT+1)) :: Fwrk
# endif
!
! Set tile array bounds.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)

# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
! In distributed-memory configuration, determine which node bounds the
! current location of the floats. Assign non-bounded floats to the
! master node.
!-----------------------------------------------------------------------
!
! The strategy here is to build a switch that processes only the floats
! contained within the node bounds. The trajectory data for the new
! time-level (nfp1) is initialized to Fspv. These values are used during
! recombining step at the end of the routine.  Since a SUM reduction is
! carried-out, setting Fspv to zero means the floats only contribute in
! their own tile.
!
      Npts=NFV(ng)*(NFT+1)*Nfloats(ng)

      Xstr=REAL(BOUNDS(ng)%Istr(MyRank),r8)-0.5_r8
      Xend=REAL(BOUNDS(ng)%Iend(MyRank),r8)+0.5_r8
      Ystr=REAL(BOUNDS(ng)%Jstr(MyRank),r8)-0.5_r8
      Yend=REAL(BOUNDS(ng)%Jend(MyRank),r8)+0.5_r8
      DO l=Lstr,Lend
        MyThread(l)=.FALSE.
        IF ((Xstr.le.track(ixgrd,nf,l)).and.                            &
     &      (track(ixgrd,nf,l).lt.Xend).and.                            &
     &      (Ystr.le.track(iygrd,nf,l)).and.                            &
     &      (track(iygrd,nf,l).lt.Yend)) THEN
          MyThread(l)=.TRUE.
        ELSE IF (Master.and.(.not.bounded(l))) THEN
          MyThread(l)=.TRUE.
        ELSE
          DO j=0,NFT
            DO i=1,NFV(ng)
              track(i,j,l)=Fspv
            END DO
          END DO
        END IF
      END DO
# else
      DO l=Lstr,Lend
        MyThread(l)=.TRUE.
      END DO
# endif

   
!
!-----------------------------------------------------------------------
!  Determine floats status.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!  If appropriate, activate the release of new floats and set initial
!  positions for all time levels.
!-----------------------------------------------------------------------
!
      DO l=Lstr,Lend
!      print*,'time=',time(ng),'Tinfo=',Tinfo(itstr,l),'bounded=',bounded
        IF (.not.bounded(l).and.                                        &
     &      (time(ng)-dt(ng).le.Tinfo(itstr,l).and.                     &
     &       time(ng)+dt(ng).gt.Tinfo(itstr,l))) THEN
          bounded(l)=.TRUE.
# ifdef DISTRIBUTE
          IF ((Xstr.le.Tinfo(ixgrd,l)).and.                             &
     &        (Tinfo(ixgrd,l).lt.Xend).and.                             &
     &        (Ystr.le.Tinfo(iygrd,l)).and.                             &
     &        (Tinfo(iygrd,l).lt.Yend)) THEN
            DO j=0,NFT
              track(ixgrd,j,l)=Tinfo(ixgrd,l)
              track(iygrd,j,l)=Tinfo(iygrd,l)
#  ifdef SOLVE3D
              track(izgrd,j,l)=Tinfo(izgrd,l)
#  endif

            END DO
            MyThread(l)=.TRUE.
          ELSE
            MyThread(l)=.FALSE.
            DO j=0,NFT
              DO i=1,NFV(ng)
                track(i,j,l)=Fspv
              END DO
            END DO
          END IF
# else
          MyThread(l)=.TRUE.
          DO j=0,NFT
            track(ixgrd,j,l)=Tinfo(ixgrd,l)
            track(iygrd,j,l)=Tinfo(iygrd,l)
#  ifdef SOLVE3D
            track(izgrd,j,l)=Tinfo(izgrd,l)
#  endif
          END DO
# endif
        END IF
      END DO
!
!  If newly released floats, initialize slopes at all time levels.
!
      DO l=Lstr,Lend
        IF (MyThread(l).and.bounded(l).and.                             &
     &      (time(ng)-dt(ng).le.Tinfo(itstr,l).and.                     &
     &       time(ng)+dt(ng).gt.Tinfo(itstr,l))) THEN
          xrhs=track(ixrhs,nfp1,l)
          yrhs=track(iyrhs,nfp1,l)
# ifdef SOLVE3D
          zrhs=track(izrhs,nfp1,l)
# endif
          DO i=0,NFT
            track(ixrhs,i,l)=xrhs
            track(iyrhs,i,l)=yrhs
# ifdef SOLVE3D
            track(izrhs,i,l)=zrhs
# endif
          END DO
        END IF
      END DO
!
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect fish on all nodes.
!-----------------------------------------------------------------------
!
      Fwrk=RESHAPE(track,(/Npts/))
      CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)
      track=RESHAPE(Fwrk,(/NFV(ng),NFT+1,Nfloats(ng)/))
# endif
!
      RETURN
      END SUBROUTINE step_floats_tile
#endif
      END MODULE new_fish_mod
