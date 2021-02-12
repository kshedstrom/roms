#include "cppdefs.h"
      MODULE step_fish_mod
#if defined NONLINEAR && defined NEMURO_SAN
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2020 The ROMS/TOMS Group        John M. Klinck   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine time-steps  simulated  fishes  trajectories using a    !
!  fourth-order Milne predictor and fourth-order Hamming corrector.    !
!                                                                      !
!  Vertical diffusion is optionally represented by a random walk,      !
!  in which case a forward scheme is used for vertical displacement.   !
!  The probability distribution for the vertical displacement is       !
!  Gaussian and includes a correction for the vertical gradient in     !
!  diffusion coefficient                                               !
!                                                                      !
! Reference:                                                           !
!                                                                      !
!  Hunter, J.R, P.D. Craig, and H.E. Philips, 1993: On the use of      !
!    random walk models with spatially variable diffusivity,           !
!    Journal of Computational Physics, 106, 366-376.                   !
!
!  Rose, K.A., Fiechter, J., Curchitser, E.N., Hedstrom, K., Bernal, M.,
!    Creekmore, S., Haynie, A., Ito, S., Lluch-Cota, S., Megrey, B.A.,
!    Edwards, C.A., Checkley, D., Koslow, T., McClatchie, S., Werner, F.,
!    MacCall, A. and Agostini, V.. (2015). Demonstration of a
!    fully-coupled end-to-end model for small pelagic fish using sardine
!    and anchovy in the California Current Progress in Oceanography.
!    doi: 10.1016/j.pocean.2015.01.012
!
!  Fiechter, J., D. D. Huff, B. T. Martin, D. W. Jackson, C. A. Edwards, K.
!    A. Rose, E. N. Curchitser, K. S. Hedstrom, S. T. Lindley, and B. K.
!    Wells. (2015). "Environmental conditions impacting juvenile Chinook
!    salmon growth off central California: An ecosystem model analysis".
!    Geophys. Res. Lett.. 42, doi: 10.1002/2015GL063046 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: step_fish

      CONTAINS
!
!***********************************************************************
      SUBROUTINE step_fish (ng, tile, Lstr, Lend)
!***********************************************************************
!
      USE mod_param
      USE mod_fish
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lstr, Lend
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 10, __LINE__, __FILE__)
# endif
      CALL step_fish_tile (ng, tile, Lstr, Lend,                        &
     &                       knew(ng), nnew(ng), nfm3(ng), nfm2(ng),    &
     &                       nfm1(ng), nf(ng), nfp1(ng),                &
     &                       FISHES(ng) % bounded,                      &
     &                       FISHES(ng) % Ftype,                        &
     &                       FISHES(ng) % Tinfo,                        &
     &                       FISHES(ng) % Fz0,                          &
     &                       FISHES(ng) % rwalk,                        &
     &                       FISHES(ng) % r2walk,                       &
     &                       FISHES(ng) % feedback,                     &
     &                       FISHES(ng) % bioenergy,                    &
     &                       FISHES(ng) % species,                      &
     &                       FISHES(ng) % alive,                        &
     &                       FISHES(ng) % lifestage,                    &
     &                       FISHES(ng) % swimtype,                     &
     &                       FISHES(ng) % deathby,                      &
     &                       FISHES(ng) % cellid,                       &
     &                       FISHES(ng) % egg_dur,                      &
     &                       FISHES(ng) % egg_num,                      &
     &                       FISHES(ng) % ysac_dur,                     &
     &                       FISHES(ng) % ysac_num,                     &
     &                       FISHES(ng) % larv_dur,                     &
     &                       FISHES(ng) % larv_num,                     &
     &                       FISHES(ng) % juv_dur,                      &
     &                       FISHES(ng) % juv_num,                      &
     &                       FISHES(ng) % suba_num,                     &
     &                       FISHES(ng) % fmortN,                       &
     &                       FISHES(ng) % fmortS,                       &
     &                       FISHES(ng) % fmortP,                       &
     &                       FISHES(ng) % fmortPsum,                    &
     &                       FISHES(ng) % num_free,                     &
     &                       FISHES(ng) % num_super,                    &
     &                       FISHES(ng) % track)
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 10, __LINE__, __FILE__)
# endif
      RETURN
      END SUBROUTINE step_fish
!
!***********************************************************************
      SUBROUTINE step_fish_tile (ng, tile, Lstr, Lend,                  &
     &                             knew, nnew,                          &
     &                             nfm3, nfm2, nfm1, nf, nfp1,          &
     &                             bounded, Ftype, Tinfo, Fz0, rwalk,   &
     &                             r2walk, feedback, bioenergy,         &
     &                             species, alive, lifestage,           &
     &                             swimtype, deathby, cellid,           &
     &                             egg_dur, egg_num,                    &
     &                             ysac_dur, ysac_num,                  &
     &                             larv_dur, larv_num,                  &
     &                             juv_dur, juv_num, suba_num,          &
     &                             fmortN, fmortS, fmortP, fmortPsum,   &
     &                             num_free, num_super, track)
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
      USE mod_fish
# ifdef PREDATOR
      USE mod_pred
# endif
      USE mod_scalars
      USE mod_biology
      USE fish_growth_mod
      USE fish_mort_mod
      USE fish_elh_mod
      USE fish_update_mod
      USE fish_swim_mod
      USE fish_land_mod
      USE interp_floats_mod
      USE mod_floats, ONLY : flt_Lagran, flt_Isobar, flt_Geopot
!
# ifdef DISTRIBUTE
      USE distribute_mod
# endif
      USE utility_mod, ONLY : nrng
      USE nrutil
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, Lstr, Lend
      integer, intent(in) :: knew, nnew, nfm3, nfm2, nfm1, nf, nfp1
!
# ifdef ASSUMED_SHAPE
      integer, intent(in)     :: Ftype(:)
      real(r8), intent(inout) :: Tinfo(0:,:)
      real(r8), intent(inout) :: Fz0(:)
      real(r8), intent(inout) :: rwalk(:)
      real(r8), intent(inout) :: r2walk(:)
      real(r8), intent(inout) :: feedback(:,:)

      logical, intent(inout)  :: bounded(:)
      real(r8), intent(inout) :: bioenergy(:,:)
      integer, intent(inout)  :: species(:)
      logical, intent(inout)  :: alive(:)
      integer, intent(inout)  :: lifestage(:)
      integer, intent(inout)  :: swimtype(:,:)
      integer, intent(inout)  :: deathby(:)
      integer, intent(inout)  :: cellid(:)
      real(r8), intent(inout) :: egg_dur(:)
      real(r8), intent(inout) :: egg_num(:)
      real(r8), intent(inout) :: ysac_dur(:)
      real(r8), intent(inout) :: ysac_num(:)
      real(r8), intent(inout) :: larv_dur(:)
      real(r8), intent(inout) :: larv_num(:)
      real(r8), intent(inout) :: juv_dur(:)
      real(r8), intent(inout) :: juv_num(:)
      real(r8), intent(inout) :: suba_num(:)
      real(r8), intent(inout) :: fmortN(:)
      real(r8), intent(inout) :: fmortS(:)
      real(r8), intent(inout) :: fmortP(:)
      real(r8), intent(inout) :: fmortPsum(:)
      integer, intent(in)     :: num_free(:)
      integer, intent(out)    :: num_super(:)
      real(r8), intent(inout) :: track(:,0:,:)
# else
      integer, intent(in)     :: Ftype(Nfish(ng))
      real(r8), intent(inout) :: Tinfo(0:izrhs,Nfish(ng))
      real(r8), intent(inout) :: Fz0(Nfish(ng))
      real(r8), intent(inout) :: rwalk(Nfish(ng)*3)
      real(r8), intent(inout) :: r2walk(Nfish(ng)*6)
      real(r8), intent(inout) :: feedback(NT(ng),Nfish(ng))

      logical, intent(inout)  :: bounded(Nfish(ng))
      real(r8), intent(inout) :: bioenergy(NFishV(ng),Nfish(ng))
      integer, intent(inout)  :: species(Nfish(ng))
      logical, intent(inout)  :: alive(Nfish(ng))
      integer, intent(inout)  :: lifestage(Nfish(ng))
      integer, intent(inout)  :: swimtype(2,Nfish(ng))
      integer, intent(inout)  :: deathby(Nfish(ng))
      integer, intent(inout)  :: cellid(Nfish(ng))
      real(r8), intent(inout) :: egg_dur(Nfish(ng))
      real(r8), intent(inout) :: egg_num(Nfish(ng))
      real(r8), intent(inout) :: ysac_dur(Nfish(ng))
      real(r8), intent(inout) :: ysac_num(Nfish(ng))
      real(r8), intent(inout) :: larv_dur(Nfish(ng))
      real(r8), intent(inout) :: larv_num(Nfish(ng))
      real(r8), intent(inout) :: juv_dur(Nfish(ng))
      real(r8), intent(inout) :: juv_num(Nfish(ng))
      real(r8), intent(inout) :: suba_num(Nfish(ng))
      real(r8), intent(inout) :: fmortN(Nfish(ng))
      real(r8), intent(inout) :: fmortS(Nfish(ng))
      real(r8), intent(inout) :: fmortP(Nfish(ng))
      real(r8), intent(inout) :: fmortPsum(Nfish(ng))
      real(r8), intent(inout) :: track(NFV(ng),0:NFT,Nfish(ng))
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
      logical, dimension(Lstr:Lend) :: MyFishThread

      integer :: Ir, Jr, Npts, i, i1, i2, j, j1, j2, itrc, l, k, isp
      integer :: NptsF, NptsL, NptsM, NptsI

      integer :: ct

      real(r8), parameter :: Fspv = 0.0_r8
      integer, parameter :: iFspv = 0
      logical, parameter :: lFspv = .false.

      real(r8) :: cff1, cff2, cff3, cff4, cff5, cff6, cff7, cff8, cff9
      real(r8) :: p1, p2, q1, q2, xrhs, yrhs, zrhs, zfloat
      real(r8) :: dtfish
!      real(r8) :: tspawn, Fworth0
!      real(r8) :: Nymort

# ifdef FLOAT_VWALK
      integer :: ierr
      integer :: iseed = 149876
# endif

      real(r8), dimension(Lstr:Lend) :: nudg

# ifdef DISTRIBUTE
      real(r8) :: Xstr, Xend, Ystr, Yend
      real(r8), dimension(Nfish(ng)*NFV(ng)*(NFT+1)) :: Fwrk
      real(r8), dimension(Nfish(ng)*NFishV(ng)) :: FwrkF
      real(r8), dimension(Nfish(ng)*Nspecies(ng)) :: FwrkM
      logical,  dimension(Nfish(ng)) :: FwrkL
      integer,  dimension(Nfish(ng)*2) :: FwrkI
# endif
!
! Set tile array bounds.
!
# include "tile.h"

      dtfish=3600.0_r8  ! one hour
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
! In distributed-memory configuration, determine which node bounds the
! current location of the fish. Assign non-bounded fish to the
! master node.
!-----------------------------------------------------------------------
!
! The strategy here is to build a switch that processes only the fish
! contained within the node bounds. The trajectory data for the new
! time-level (nfp1) is initialized to Fspv. These values are used during
! recombining step at the end of the routine.  Since a SUM reduction is
! carried-out, setting Fspv to zero means the fish only contribute in
! their own tile.
!
      Npts=NFV(ng)*(NFT+1)*Nfish(ng)
      NptsF=NFishV(ng)*Nfish(ng)
      NptsL=Nfish(ng)
      NptsM=Nspecies(ng)*Nfish(ng)
      NptsI=2*Nfish(ng)

      Xstr=REAL(BOUNDS(ng)%Istr(MyRank),r8)-0.5_r8
      Xend=REAL(BOUNDS(ng)%Iend(MyRank),r8)+0.5_r8
      Ystr=REAL(BOUNDS(ng)%Jstr(MyRank),r8)-0.5_r8
      Yend=REAL(BOUNDS(ng)%Jend(MyRank),r8)+0.5_r8
      DO l=Lstr,Lend
        MyFishThread(l)=.FALSE.
        IF ((Xstr.le.track(ixgrd,nf,l)).and.                            &
     &      (track(ixgrd,nf,l).lt.Xend).and.                            &
     &      (Ystr.le.track(iygrd,nf,l)).and.                            &
     &      (track(iygrd,nf,l).lt.Yend)) THEN
          MyFishThread(l)=.TRUE.
        ELSE IF (Master.and.(.not.bounded(l))) THEN
          MyFishThread(l)=.TRUE.
        ELSE
          DO j=0,NFT
            DO i=1,NFV(ng)
              track(i,j,l)=Fspv
            END DO
          END DO
          DO i=1,NFishV(ng)
            bioenergy(i,l)=Fspv
          END DO
          egg_dur(l)=Fspv
          egg_num(l)=Fspv
          ysac_dur(l)=Fspv
          ysac_num(l)=Fspv
          larv_dur(l)=Fspv
          larv_num(l)=Fspv
          juv_dur(l)=Fspv
          juv_num(l)=Fspv
          suba_num(l)=Fspv
          fmortN(l)=Fspv
          fmortS(l)=Fspv
          fmortP(l)=Fspv
          fmortPsum(l)=Fspv
          lifestage(l)=iFspv
          swimtype(1,l)=iFspv
          swimtype(2,l)=iFspv
          deathby(l)=iFspv
          cellid(l)=iFspv
          alive(l)=lFspv
        END IF
      END DO
# else
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      DO l=Lstr,Lend
        MyFishThread(l)=.TRUE.
      END DO
# endif
# ifdef DISTRIBUTE
      IF (Master) THEN
        CALL gasdev (rwalk)
      END IF
      CALL mp_bcastf (ng, iNLM, rwalk)
# elif defined _OPENMP
!$OMP SINGLE
      CALL gasdev (rwalk)
!$OMP END SINGLE
# else
      IF (Lstr.eq.1) THEN
        CALL gasdev (rwalk)
      END IF
# endif
!
!-----------------------------------------------------------------------
! Compute nudging velocities for vertical random walk.
!-----------------------------------------------------------------------
!
# if defined SOLVE3D && defined FLOAT_VWALK

      DO l=Lstr,Lend
        nudg(l)=0.0_r8
      END DO

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nf, ifakt, Nfish(ng),             &
     &                    isBw3d, w3dvar, Lmask, spval, nudg,           &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    MIXING(ng) % Akt(:,:,:,itemp),                &
     &                    MyFishThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nf, ifdak, Nfish(ng),             &
     &                    isBr3d, r3dvar, Lmask, spval, nudg,           &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    MIXING(ng) % dAktdz,                          &
     &                    MyFishThread, bounded, track)

      CALL nrng (iseed, nudg(Lstr:Lend), Lend-Lstr+1, ierr)
      DO l=Lstr,Lend
!       nudg(l) = rwalk(l)
        IF (MyFishThread(l).and.bounded(l)) THEN
          nudg(l)=SQRT(2.0_r8*MAX(track(ifakt,nf,l),0.0_r8)/dtfish)*    &
     &            nudg(l)+track(ifdak,nf,l)
        ELSE
          nudg(l)=0.0_r8
        END IF
      END DO

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nf, izrhs, Nfish(ng),             &
     &                    isBw3d, -w3dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % W,                                &
     &                    MyFishThread, bounded, track)
!
!  For now, use a forward time step for vertical position.
!
      DO l=Lstr,Lend
        IF (MyFishThread(l).and.bounded(l)) THEN
          track(izgrd,nfp1,l)=track(izgrd,nf,l)+                        &
     &                        dtfish*track(izrhs,nf,l)
        END IF
      END DO
# endif

      DO l=Lstr,Lend
        nudg(l)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Predictor step: compute first guess fish locations using a
!                  4th-order Milne time-stepping scheme.
!-----------------------------------------------------------------------
!
      cff1=8.0_r8/3.0_r8
      cff2=4.0_r8/3.0_r8
!
      DO l=Lstr,Lend
        IF (MyFishThread(l).and.bounded(l).and.                         &
     &      alive(l).and.(lifestage(l).le.if_larva)) THEN
!
          track(ixgrd,nfp1,l)=track(ixgrd,nfm3,l)+                      &
     &                        dtfish*(cff1*track(ixrhs,nf  ,l)-         &
     &                                cff2*track(ixrhs,nfm1,l)+         &
     &                                cff1*track(ixrhs,nfm2,l))
          track(iygrd,nfp1,l)=track(iygrd,nfm3,l)+                      &
     &                        dtfish*(cff1*track(iyrhs,nf  ,l)-         &
     &                                cff2*track(iyrhs,nfm1,l)+         &
     &                                cff1*track(iyrhs,nfm2,l))
        END IF
!
# if defined SOLVE3D && !defined FLOAT_VWALK
!
!  Compute vertical position (grid units) 3D Lagrangian fish.
!
        IF (MyFishThread(l).and.bounded(l).and.                         &
     &      alive(l).and.(lifestage(l).le.if_yolksac)) THEN
          IF (Ftype(l).eq.flt_Lagran) THEN
            track(izgrd,nfp1,l)=track(izgrd,nfm3,l)+                    &
     &                          dtfish*(cff1*track(izrhs,nf  ,l)-       &
     &                                  cff2*track(izrhs,nfm1,l)+       &
     &                                  cff1*track(izrhs,nfm2,l))
!
!  Compute vertical position (grid units) for isobaric fish
!  (p=g*(z+zeta)=constant) or geopotential fish (constant depth).
!  Use bilinear interpolation to determine vertical position.
!
          ELSE IF ((Ftype(l).eq.flt_Isobar).or.                         &
     &             (Ftype(l).eq.flt_Geopot)) THEN
            Ir=NINT(track(ixgrd,nfp1,l))
            Jr=NINT(track(iygrd,nfp1,l))
!
            i1=MIN(MAX(Ir  ,0),Lm(ng)+1)
            i2=MIN(MAX(Ir+1,1),Lm(ng)+1)
            j1=MIN(MAX(Jr  ,0),Mm(ng)+1)
            j2=MIN(MAX(Jr+1,0),Mm(ng)+1)
!
            p2=REAL(i2-i1,r8)*(track(ixgrd,nfp1,l)-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(track(iygrd,nfp1,l)-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
#  ifdef MASKING
            cff7=p1*q1*GRID(ng)%z_w(i1,j1,N(ng))*GRID(ng)%rmask(i1,j1)+ &
     &           p2*q1*GRID(ng)%z_w(i2,j1,N(ng))*GRID(ng)%rmask(i2,j1)+ &
     &           p1*q2*GRID(ng)%z_w(i1,j2,N(ng))*GRID(ng)%rmask(i1,j2)+ &
     &           p2*q2*GRID(ng)%z_w(i2,j2,N(ng))*GRID(ng)%rmask(i2,j2)
            cff8=p1*q1*GRID(ng)%rmask(i1,j1)+                           &
     &           p2*q1*GRID(ng)%rmask(i2,j1)+                           &
     &           p1*q2*GRID(ng)%rmask(i1,j2)+                           &
     &           p2*q2*GRID(ng)%rmask(i2,j2)
            cff9=0.0_r8
            IF (cff8.gt.0.0_r8) cff9=cff7/cff8
#  else
            cff9=p1*q1*GRID(ng)%z_w(i1,j1,N(ng))+                       &
     &           p2*q1*GRID(ng)%z_w(i2,j1,N(ng))+                       &
     &           p1*q2*GRID(ng)%z_w(i1,j2,N(ng))+                       &
     &           p2*q2*GRID(ng)%z_w(i2,j2,N(ng))
#  endif
            cff6=cff9
!
            IF (Ftype(l).eq.flt_Geopot) THEN
              zfloat=Fz0(l)
            ELSE IF (Ftype(l).eq.flt_Isobar) THEN
              zfloat=Fz0(l)+cff9
            END IF
!
            DO k=N(ng),0,-1
!            DO k=N(ng)-1,1,-1
#  ifdef MASKING
              cff7=p1*q1*GRID(ng)%z_w(i1,j1,k)*GRID(ng)%rmask(i1,j1)+   &
     &             p2*q1*GRID(ng)%z_w(i2,j1,k)*GRID(ng)%rmask(i2,j1)+   &
     &             p1*q2*GRID(ng)%z_w(i1,j2,k)*GRID(ng)%rmask(i1,j2)+   &
     &             p2*q2*GRID(ng)%z_w(i2,j2,k)*GRID(ng)%rmask(i2,j2)
              cff8=p1*q1*GRID(ng)%rmask(i1,j1)+                         &
     &             p2*q1*GRID(ng)%rmask(i2,j1)+                         &
     &             p1*q2*GRID(ng)%rmask(i1,j2)+                         &
     &             p2*q2*GRID(ng)%rmask(i2,j2)
              IF (cff8.gt.0.0_r8) THEN
                cff5=cff7/cff8
              ELSE
                cff5=0.0_r8
              END IF
#  else
              cff5=p1*q1*GRID(ng)%z_w(i1,j1,k)+                         &
     &             p2*q1*GRID(ng)%z_w(i2,j1,k)+                         &
     &             p1*q2*GRID(ng)%z_w(i1,j2,k)+                         &
     &             p2*q2*GRID(ng)%z_w(i2,j2,k)
#  endif
              IF ((zfloat-cff5)*(cff6-zfloat).ge.0.0_r8) THEN
                track(izgrd,nfp1,l)=REAL(k,r8)+(zfloat-cff5)/(cff6-cff5)
              END IF
              cff6=cff5
            END DO
! Check for shallowing water
            IF (zfloat < cff5) THEN
              print *, 'STEP_FLOAT getting too shallow', zfloat, cff5, l
              Fz0(l) = cff5 + 5.0
            END IF
          ELSE
            track(izgrd,nfp1,l)=track(izgrd,nfm3,l)+                    &
     &                          dtfish*(cff1*track(izrhs,nf  ,l)-       &
     &                                  cff2*track(izrhs,nfm1,l)+       &
     &                                  cff1*track(izrhs,nfm2,l))
          END IF
# endif
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Calculate slopes at new time-step.
!-----------------------------------------------------------------------
!
# ifdef SOLVE3D
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, ixrhs, Nfish(ng),           &
     &                    isUvel, -u3dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % u(:,:,:,nnew),                    &
     &                    MyFishThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, iyrhs, Nfish(ng),           &
     &                    isVvel, -v3dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % v(:,:,:,nnew),                    &
     &                    MyFishThread, bounded, track)

#  if !defined FLOAT_VWALK
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nfp1, izrhs, Nfish(ng),           &
     &                    isBw3d, -w3dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % W,                                &
     &                    MyFishThread, bounded, track)
#  endif
# else
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfp1, ixrhs, Nfish(ng),           &
     &                    isUbar, -u2dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % ubar(:,:,knew),                   &
     &                    MyFishThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfp1, iyrhs, Nfish(ng),           &
     &                    isVbar, -v2dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % vbar(:,:,knew),                   &
     &                    MyFishThread, bounded, track)
# endif
!
!-----------------------------------------------------------------------
!  Corrector step: correct fish locations using a 4th-order
!                  Hamming time-stepping scheme.
!-----------------------------------------------------------------------
!
      cff1=9.0_r8/8.0_r8
      cff2=1.0_r8/8.0_r8
      cff3=3.0_r8/8.0_r8
      cff4=6.0_r8/8.0_r8
!
      DO l=Lstr,Lend
        IF (MyFishThread(l).and.bounded(l).and.                         &
     &      alive(l).and.(lifestage(l).le.if_larva)) THEN
!
          track(ixgrd,nfp1,l)=cff1*track(ixgrd,nf  ,l)-                 &
     &                        cff2*track(ixgrd,nfm2,l)+                 &
     &                        dtfish*(cff3*track(ixrhs,nfp1,l)+         &
     &                                cff4*track(ixrhs,nf  ,l)-         &
     &                                cff3*track(ixrhs,nfm1,l))
          track(iygrd,nfp1,l)=cff1*track(iygrd,nf  ,l)-                 &
     &                        cff2*track(iygrd,nfm2,l)+                 &
     &                        dtfish*(cff3*track(iyrhs,nfp1,l)+         &
     &                                cff4*track(iyrhs,nf  ,l)-         &
     &                                cff3*track(iyrhs,nfm1,l))
        END IF
!
# if defined SOLVE3D && !defined FLOAT_VWALK
!
!  Compute vertical position (grid units) 3D Lagrangian fish.
!
        IF (MyFishThread(l).and.bounded(l).and.                         &
     &      alive(l).and.(lifestage(l).le.if_yolksac)) THEN
          IF (Ftype(l).eq.flt_Lagran) THEN
!            track(izgrd,nfp1,l)=cff1*track(izgrd,nf  ,l)-               &
!     &                          cff2*track(izgrd,nfm2,l)+               &
!     &                          dtfish*(cff3*track(izrhs,nfp1,l)+       &
!     &                                  cff4*track(izrhs,nf  ,l)-       &
!     &                                  cff3*track(izrhs,nfm1,l))
            track(izgrd,nfp1,l)=track(izgrd,nf,l)+                      &
     &                          dtfish*track(izrhs,nf,l)
!
!  Compute vertical position (grid units) for isobaric fish
!  (p=g*(z+zeta)=constant) or geopotential fish (constant depth).
!  Use bilinear interpolation to determine vertical position.
!
          ELSE IF ((Ftype(l).eq.flt_Isobar).or.                         &
     &             (Ftype(l).eq.flt_Geopot)) THEN
            Ir=NINT(track(ixgrd,nfp1,l))
            Jr=NINT(track(iygrd,nfp1,l))
!
            i1=MIN(MAX(Ir  ,0),Lm(ng)+1)
            i2=MIN(MAX(Ir+1,1),Lm(ng)+1)
            j1=MIN(MAX(Jr  ,0),Mm(ng)+1)
            j2=MIN(MAX(Jr+1,0),Mm(ng)+1)
!
            p2=REAL(i2-i1,r8)*(track(ixgrd,nfp1,l)-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(track(iygrd,nfp1,l)-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
#  ifdef MASKING
            cff7=p1*q1*GRID(ng)%z_w(i1,j1,N(ng))*GRID(ng)%rmask(i1,j1)+ &
     &           p2*q1*GRID(ng)%z_w(i2,j1,N(ng))*GRID(ng)%rmask(i2,j1)+ &
     &           p1*q2*GRID(ng)%z_w(i1,j2,N(ng))*GRID(ng)%rmask(i1,j2)+ &
     &           p2*q2*GRID(ng)%z_w(i2,j2,N(ng))*GRID(ng)%rmask(i2,j2)
            cff8=p1*q1*GRID(ng)%rmask(i1,j1)+                           &
     &           p2*q1*GRID(ng)%rmask(i2,j1)+                           &
     &           p1*q2*GRID(ng)%rmask(i1,j2)+                           &
     &           p2*q2*GRID(ng)%rmask(i2,j2)
            IF (cff8.gt.0.0_r8) THEN
              cff9=cff7/cff8
            ELSE
              cff9=0.0_r8
            END IF
#  else
            cff9=p1*q1*GRID(ng)%z_w(i1,j1,N(ng))+                       &
     &           p2*q1*GRID(ng)%z_w(i2,j1,N(ng))+                       &
     &           p1*q2*GRID(ng)%z_w(i1,j2,N(ng))+                       &
     &           p2*q2*GRID(ng)%z_w(i2,j2,N(ng))
#  endif
            cff6=cff9
!
            IF (Ftype(l).eq.flt_Geopot) THEN
              zfloat=Fz0(l)
            ELSE IF (Ftype(l).eq.flt_Isobar) THEN
              zfloat=Fz0(l)+cff9
            END IF
!
            DO k=N(ng)-1,1,-1
#  ifdef MASKING
              cff7=p1*q1*GRID(ng)%z_w(i1,j1,k)*GRID(ng)%rmask(i1,j1)+   &
     &             p2*q1*GRID(ng)%z_w(i2,j1,k)*GRID(ng)%rmask(i2,j1)+   &
     &             p1*q2*GRID(ng)%z_w(i1,j2,k)*GRID(ng)%rmask(i1,j2)+   &
     &             p2*q2*GRID(ng)%z_w(i2,j2,k)*GRID(ng)%rmask(i2,j2)
              cff8=p1*q1*GRID(ng)%rmask(i1,j1)+                         &
     &             p2*q1*GRID(ng)%rmask(i2,j1)+                         &
     &             p1*q2*GRID(ng)%rmask(i1,j2)+                         &
     &             p2*q2*GRID(ng)%rmask(i2,j2)
              cff5=0.0_r8
              IF (cff8.gt.0.0_r8) cff5=cff7/cff8
#  else
              cff5=p1*q1*GRID(ng)%z_w(i1,j1,k)+                         &
     &             p2*q1*GRID(ng)%z_w(i2,j1,k)+                         &
     &             p1*q2*GRID(ng)%z_w(i1,j2,k)+                         &
     &             p2*q2*GRID(ng)%z_w(i2,j2,k)
#  endif
              IF ((zfloat-cff5)*(cff6-zfloat).ge.0.0_r8) THEN
                track(izgrd,nfp1,l)=REAL(k,r8)+(zfloat-cff5)/(cff6-cff5)
              END IF
              cff6=cff5
            END DO
! Check for shallowing water
            IF (zfloat < cff5) THEN
              print *, 'STEP_FLOAT getting too shallow', zfloat,        &
     &            cff5, l
              Fz0(l) = cff5 + 5.0
            END IF
          ELSE
            track(izgrd,nfp1,l)=cff1*track(izgrd,nf  ,l)-               &
     &                          cff2*track(izgrd,nfm2,l)+               &
     &                          dtfish*(cff3*track(izrhs,nfp1,l)+       &
     &                                  cff4*track(izrhs,nf  ,l)-       &
     &                                  cff3*track(izrhs,nfm1,l))
          END IF
# endif
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Fish swimming behavior for juveniles and adults
!-----------------------------------------------------------------------
!
      CALL fish_swim (ng, tile, LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                    nfm2, nfm1, nf, nfp1, nnew,                   &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % h,                                 &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
# ifdef SOLVE3D
     &                    GRID(ng) % Hz,                                &
# endif
     &                    MyFishThread, bounded, rwalk, r2walk,         &
     &                    track, bioenergy, alive,                      &
     &                    species, lifestage, swimtype,                 &
     &                    OCEAN(ng) % fish_count,                       &
     &                    OCEAN(ng) % fish_list,                        &
     &                    FISHES(ng) % fishnodes)
!
!-----------------------------------------------------------------------
!  Check for land
!-----------------------------------------------------------------------
!
      CALL fish_land (ng, tile, LBi, UBi, LBj, UBj,                     &
     &                    nfm1, nf, nfp1,                               &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    MyFishThread, bounded, track, alive,          &
     &                    OCEAN(ng) % fish_count,                       &
     &                    OCEAN(ng) % fish_list,                        &
     &                    FISHES(ng) % fishnodes)
!
!-----------------------------------------------------------------------
! Calculate "cell-by-fish" for fishing fleet
! (don't know how to dealwith fish_list on entire grid)
!-----------------------------------------------------------------------
!
      DO l=Lstr,Lend
        IF (MyFishThread(l).and.bounded(l)) THEN
          Ir=NINT(track(ixgrd,nf,l))
          Jr=NINT(track(iygrd,nf,l))
          cellid(l)=Ir+(Lm(ng)+2)*Jr
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Determine fish status.
!-----------------------------------------------------------------------
!
# ifdef EW_PERIODIC
      cff1=REAL(Lm(ng),r8)
      DO l=Lstr,Lend
        IF (MyFishThread(l).and.bounded(l)) THEN
          IF (track(ixgrd,nfp1,l).ge.REAL(Lm(ng)+1,r8)-0.5_r8) THEN
            track(ixgrd,nfp1,l)=track(ixgrd,nfp1,l)-cff1
            track(ixgrd,nf  ,l)=track(ixgrd,nf  ,l)-cff1
            track(ixgrd,nfm1,l)=track(ixgrd,nfm1,l)-cff1
            track(ixgrd,nfm2,l)=track(ixgrd,nfm2,l)-cff1
            track(ixgrd,nfm3,l)=track(ixgrd,nfm3,l)-cff1
          ELSE IF (track(ixgrd,nfp1,l).lt.0.5_r8) THEN
            track(ixgrd,nfp1,l)=cff1+track(ixgrd,nfp1,l)
            track(ixgrd,nf  ,l)=cff1+track(ixgrd,nf  ,l)
            track(ixgrd,nfm1,l)=cff1+track(ixgrd,nfm1,l)
            track(ixgrd,nfm2,l)=cff1+track(ixgrd,nfm2,l)
            track(ixgrd,nfm3,l)=cff1+track(ixgrd,nfm3,l)
          END IF
        END IF
      END DO
#  ifdef DISTRIBUTE
      IF (NtileI(ng).gt.1) THEN
        Fwrk=RESHAPE(track,(/Npts/))
        CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)
        track=RESHAPE(Fwrk,(/NFV(ng),NFT+1,Nfish(ng)/))
        DO l=Lstr,Lend
          IF ((Xstr.le.track(ixgrd,nfp1,l)).and.                        &
     &        (track(ixgrd,nfp1,l).lt.Xend).and.                        &
     &        (Ystr.le.track(iygrd,nfp1,l)).and.                        &
     &        (track(iygrd,nfp1,l).lt.Yend)) THEN
            MyFishThread(l)=.TRUE.
          ELSE IF (Master.and.(.not.bounded(l))) THEN
            MyFishThread(l)=.TRUE.
          ELSE
            MyFishThread(l)=.FALSE.
            DO j=0,NFT
              DO i=1,NFV(ng)
                track(i,j,l)=Fspv
              END DO
            END DO
            DO i=1,NFishV(ng)
              bioenergy(i,l)=Fspv
            END DO
            egg_dur(l)=Fspv
            egg_num(l)=Fspv
            ysac_dur(l)=Fspv
            ysac_num(l)=Fspv
            larv_dur(l)=Fspv
            larv_num(l)=Fspv
            juv_dur(l)=Fspv
            juv_num(l)=Fspv
            suba_num(l)=Fspv
            fmortN(l)=Fspv
            fmortS(l)=Fspv
            fmortP(l)=Fspv
            fmortPsum(l)=Fspv
            lifestage(l)=iFspv
            swimtype(1,l)=iFspv
            swimtype(2,l)=iFspv
            deathby(l)=iFspv
            cellid(l)=iFspv
            alive(l)=lFspv
          END IF
        END DO
      END IF
#  endif
# else
      DO l=Lstr,Lend
        IF (MyFishThread(l).and.bounded(l)) THEN
! JF: Relocation at starting position for non-periodic case
!          IF ((track(ixgrd,nfp1,l).ge.REAL(Lm(ng)+1,r8)-0.5_r8).or.     &
!     &        (track(ixgrd,nfp1,l).lt.0.5_r8)) THEN
!            DO j=0,NFT
!              track(ixgrd,j,l)=Tinfo(ixgrd,l)
!              track(iygrd,j,l)=Tinfo(iygrd,l)
!              track(izgrd,j,l)=Tinfo(izgrd,l)
!              track(ixrhs,j,l)=0.0_r8
!              track(iyrhs,j,l)=0.0_r8
!              track(izrhs,j,l)=0.0_r8
!            END DO
!          END IF
! JF: Crude reflective boundary condition for non-periodic case
          IF (track(ixgrd,nfp1,l).ge.REAL(Lm(ng)+1,r8)-1.5_r8) THEN
            DO j=0,NFT
              track(ixgrd,j,l)=REAL(Lm(ng)+1-2)
              track(ixrhs,j,l)=-track(ixrhs,nfp1,l)
            END DO
          END IF
          IF (track(ixgrd,nfp1,l).lt.1.5_r8) THEN
            DO j=0,NFT
              track(ixgrd,j,l)=2.0_r8
              track(ixrhs,j,l)=-track(ixrhs,nfp1,l)
            END DO
          END IF
! Original ROMS code
!          IF ((track(ixgrd,nfp1,l).ge.REAL(Lm(ng)+1,r8)-0.5_r8).or.     &
!     &        (track(ixgrd,nfp1,l).lt.0.5_r8)) THEN
!            bounded(l)=.FALSE.
!          END IF
        END IF
      END DO
# endif
# ifdef NS_PERIODIC
      cff1=REAL(Mm(ng),r8)
      DO l=Lstr,Lend
        IF (MyFishThread(l).and.bounded(l)) THEN
          IF (track(iygrd,nfp1,l).ge.REAL(Mm(ng)+1,r8)-0.5_r8) THEN
            track(iygrd,nfp1,l)=track(iygrd,nfp1,l)-cff1
            track(iygrd,nf  ,l)=track(iygrd,nf  ,l)-cff1
            track(iygrd,nfm1,l)=track(iygrd,nfm1,l)-cff1
            track(iygrd,nfm2,l)=track(iygrd,nfm2,l)-cff1
            track(iygrd,nfm3,l)=track(iygrd,nfm3,l)-cff1
          ELSE IF (track(iygrd,nfp1,l).lt.0.5_r8) THEN
            track(iygrd,nfp1,l)=cff1+track(iygrd,nfp1,l)
            track(iygrd,nf  ,l)=cff1+track(iygrd,nf  ,l)
            track(iygrd,nfm1,l)=cff1+track(iygrd,nfm1,l)
            track(iygrd,nfm2,l)=cff1+track(iygrd,nfm2,l)
            track(iygrd,nfm3,l)=cff1+track(iygrd,nfm3,l)
          END IF
        END IF
      END DO
#  ifdef DISTRIBUTE
      IF (NtileJ(ng).gt.1) THEN
        Fwrk=RESHAPE(track,(/Npts/))
        CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)
        track=RESHAPE(Fwrk,(/NFV(ng),NFT+1,Nfish(ng)/))
        DO l=Lstr,Lend
          IF ((Xstr.le.track(ixgrd,nfp1,l)).and.                        &
     &        (track(ixgrd,nfp1,l).lt.Xend).and.                        &
     &        (Ystr.le.track(iygrd,nfp1,l)).and.                        &
     &        (track(iygrd,nfp1,l).lt.Yend)) THEN
            MyFishThread(l)=.TRUE.
          ELSE IF (Master.and.(.not.bounded(l))) THEN
            MyFishThread(l)=.TRUE.
          ELSE
            MyFishThread(l)=.FALSE.
            DO j=0,NFT
              DO i=1,NFV(ng)
                track(i,j,l)=Fspv
              END DO
            END DO
            DO i=1,NFishV(ng)
              bioenergy(i,l)=Fspv
            END DO
            egg_dur(l)=Fspv
            egg_num(l)=Fspv
            ysac_dur(l)=Fspv
            ysac_num(l)=Fspv
            larv_dur(l)=Fspv
            larv_num(l)=Fspv
            juv_dur(l)=Fspv
            juv_num(l)=Fspv
            suba_num(l)=Fspv
            fmortN(l)=Fspv
            fmortS(l)=Fspv
            fmortP(l)=Fspv
            fmortPsum(l)=Fspv
            lifestage(l)=iFspv
            swimtype(1,l)=iFspv
            swimtype(2,l)=iFspv
            deathby(l)=iFspv
            cellid(l)=iFspv
            alive(l)=lFspv
          END IF
        END DO
      END IF
#  endif
# else
      DO l=Lstr,Lend
        IF (MyFishThread(l).and.bounded(l)) THEN
! JF: Relocation at starting position for non-periodic case
!          IF ((track(iygrd,nfp1,l).ge.REAL(Mm(ng)+1,r8)-0.5_r8).or.     &
!     &        (track(iygrd,nfp1,l).lt.0.5_r8)) THEN
!            DO j=0,NFT
!              track(ixgrd,j,l)=Tinfo(ixgrd,l)
!              track(iygrd,j,l)=Tinfo(iygrd,l)
!              track(izgrd,j,l)=Tinfo(izgrd,l)
!              track(ixrhs,j,l)=0.0_r8
!              track(iyrhs,j,l)=0.0_r8
!              track(izrhs,j,l)=0.0_r8
!            END DO
!          END IF
! JF: Crude reflective boundary condition for non-periodic case
          IF (track(iygrd,nfp1,l).ge.REAL(Mm(ng)+1,r8)-1.5_r8) THEN
            DO j=0,NFT
              track(iygrd,j,l)=REAL(Mm(ng)+1-2)
              track(iyrhs,j,l)=-track(iyrhs,nfp1,l)
            END DO
          END IF
          IF (track(iygrd,nfp1,l).lt.1.5_r8) THEN
            DO j=0,NFT
              track(iygrd,j,l)=2.0_r8
              track(iyrhs,j,l)=-track(iyrhs,nfp1,l)
            END DO
          END IF
! Original ROMS code
!          IF ((track(iygrd,nfp1,l).ge.REAL(Mm(ng)+1,r8)-0.5_r8).or.     &
!     &        (track(iygrd,nfp1,l).lt.0.5_r8)) THEN
!            bounded(l)=.FALSE.
!          END IF
        END IF
      END DO
# endif
# ifdef SOLVE3D
!
!  Reflect fish at surface or bottom.
!
      DO l=Lstr,Lend
        IF (MyFishThread(l).and.bounded(l)) THEN
          IF (track(izgrd,nfp1,l).gt.REAL(N(ng),r8))                    &
     &      track(izgrd,nfp1,l)=REAL(N(ng),r8)-0.5_r8
          IF (track(izgrd,nfp1,l).lt.0.0_r8)                            &
     &      track(izgrd,nfp1,l)=0.5_r8
        END IF
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Calculate slopes with corrected locations.
!-----------------------------------------------------------------------
!
# ifdef SOLVE3D
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, ixrhs, Nfish(ng),           &
     &                    isUvel, -u3dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % u(:,:,:,nnew),                    &
     &                    MyFishThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, iyrhs, Nfish(ng),           &
     &                    isVvel, -v3dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % v(:,:,:,nnew),                    &
     &                    MyFishThread, bounded, track)

#  if !defined FLOAT_VWALK
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nfp1, izrhs, Nfish(ng),           &
     &                    isBw3d, -w3dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % W,                                &
     &                    MyFishThread, bounded, track)
#  endif
# else
      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfp1, ixrhs, Nfish(ng),           &
     &                    isUbar, -u2dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % ubar(:,:,knew),                   &
     &                    MyFishThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,                 &
     &                    Lstr, Lend, nfp1, iyrhs, Nfish(ng),           &
     &                    isVbar, -v2dvar, Lmask, spval, nudg,          &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % vbar(:,:,knew),                   &
     &                    MyFishThread, bounded, track)
# endif
!
!  If newly released fish, initialize slopes at all time levels.
!
!      DO l=Lstr,Lend
!        IF (MyFishThread(l).and.bounded(l).and.                             &
!     &      (time(ng)-dtfish.le.Tinfo(itstr,l).and.                     &
!     &       time(ng)+dtfish.gt.Tinfo(itstr,l))) THEN
!          xrhs=track(ixrhs,nfp1,l)
!          yrhs=track(iyrhs,nfp1,l)
!# ifdef SOLVE3D
!          zrhs=track(izrhs,nfp1,l)
!# endif
!          DO i=0,NFT
!            track(ixrhs,i,l)=xrhs
!            track(iyrhs,i,l)=yrhs
!# ifdef SOLVE3D
!            track(izrhs,i,l)=zrhs
!# endif
!          END DO
!        END IF
!      END DO
!
!-----------------------------------------------------------------------
!  Interpolate various output variables at the corrected locations.
!-----------------------------------------------------------------------
!
      IF (spherical) THEN
        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,               &
     &                      Lstr, Lend, nfp1, iflon, Nfish(ng),         &
     &                      isBr2d, r2dvar, Gmask, spval, nudg,         &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      GRID(ng) % lonr,                            &
     &                      MyFishThread, bounded, track)

        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,               &
     &                      Lstr, Lend, nfp1, iflat, Nfish(ng),         &
     &                      isBr2d, r2dvar, Gmask, spval, nudg,         &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      GRID(ng) % latr,                            &
     &                      MyFishThread, bounded, track)
      ELSE
        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,               &
     &                      Lstr, Lend, nfp1, iflon, Nfish(ng),         &
     &                      isBr2d, r2dvar, Gmask, spval, nudg,         &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      GRID(ng) % xr,                              &
     &                      MyFishThread, bounded, track)

        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, 1,               &
     &                      Lstr, Lend, nfp1, iflat, Nfish(ng),         &
     &                      isBr2d, r2dvar, Gmask, spval, nudg,         &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
# endif
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      GRID(ng) % yr,                              &
     &                      MyFishThread, bounded, track)
      END IF
# ifdef SOLVE3D

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                    Lstr, Lend, nfp1, idpth, Nfish(ng),           &
     &                    isBw3d, w3dvar, Lmask, spval, nudg,           &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    GRID(ng) % z_w,                               &
     &                    MyFishThread, bounded, track)

      CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                    Lstr, Lend, nfp1, ifden, Nfish(ng),           &
     &                    isBr3d, r3dvar, Lmask, spval, nudg,           &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % Hz,                                &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    OCEAN(ng) % rho,                              &
     &                    MyFishThread, bounded, track)

      DO itrc=1,NT(ng)
        CALL interp_floats (ng, LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                      Lstr, Lend, nfp1, itrc+NFV(ng)-NT(ng),      &
     &                      Nfish(ng), isTvar(itrc),                    &
     &                      r3dvar, Lmask, spval, nudg,                 &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      GRID(ng) % Hz,                              &
# ifdef MASKING
     &                      GRID(ng) % rmask,                           &
# endif
     &                      OCEAN(ng) % t(:,:,:,nnew,itrc),             &
     &                      MyFishThread, bounded, track)
      END DO
!
      CALL fish_growth (ng, tile, LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                      nfp1, nnew,                                 &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
# ifdef SOLVE3D
     &                      GRID(ng) % Hz,                              &
# endif
     &                      MyFishThread, bounded, track(:,:,:),        &
     &                      feedback, bioenergy, alive, species,        &
     &                      lifestage,                                  &
     &                      OCEAN(ng) % fish_count,                     &
     &                      OCEAN(ng) % fish_list,                      &
     &                      FISHES(ng) % fishnodes)
!
      CALL fish_mort (ng, tile, LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                    nfp1, nnew,                                   &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % pn,                                &
     &                    MyFishThread, bounded, track(:,:,:),          &
     &                    bioenergy, alive, species,                    &
     &                    lifestage, deathby,                           &
     &                    fmortN, fmortS, fmortP, fmortPsum,            &
# ifdef PREDATOR
     &                    OCEAN(ng) % pred_count,                       &
     &                    OCEAN(ng) % pred_list,                        &
     &                    PREDS(ng) % prednodes,                        &
     &                    PREDS(ng) % bioenergy,                        &
     &                    PREDS(ng) % species,                          &
     &                    PREDS(ng) % track,                            &
# endif
     &                    OCEAN(ng) % fish_count,                       &
     &                    OCEAN(ng) % fish_list,                        &
     &                    FISHES(ng) % fishnodes)
!
      CALL fish_elh (ng, tile, LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nfp1, nnew,                                    &
     &                   MyFishThread, bounded, track(:,:,:),           &
     &                   bioenergy, alive, species,                     &
     &                   lifestage, deathby, fmortN,                    &
     &                   egg_dur, egg_num,                              &
     &                   ysac_dur, ysac_num,                            &
     &                   larv_dur, larv_num,                            &
     &                   juv_dur, juv_num, suba_num,                    &
     &                   OCEAN(ng) % fish_count,                        &
     &                   OCEAN(ng) % fish_list,                         &
     &                   FISHES(ng) % fishnodes)
!
      CALL fish_update (ng, tile, LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                      nfp1, nnew,                                 &
     &                      MyFishThread, bounded, track(:,:,:),        &
     &                      bioenergy, alive, species,                  &
     &                      lifestage, deathby,                         &
     &                      egg_dur, egg_num,                           &
     &                      ysac_dur, ysac_num,                         &
     &                      larv_dur, larv_num,                         &
     &                      juv_dur, juv_num, suba_num,                 &
     &                      OCEAN(ng) % fish_count,                     &
     &                      OCEAN(ng) % fish_list,                      &
     &                      FISHES(ng) % fishnodes)

!
# endif
# ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Collect fish on all nodes.
!-----------------------------------------------------------------------
!

      !Fwrk=RESHAPE(track,(/Npts/))
      !CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)
      !track=RESHAPE(Fwrk,(/NFV(ng),NFT+1,Nfish(ng)/))
      ct=1
      DO l=1,Nfish(ng)
        DO j=0,NFT
          DO i=1,NFV(ng)
               Fwrk(ct) = track(i,j,l)
               ct=ct+1
          ENDDO
        ENDDO
      ENDDO

      CALL mp_collect (ng, iNLM, Npts, Fspv, Fwrk)

      ct=1
      DO l=1,Nfish(ng)
        DO j=0,NFT
          DO i=1,NFV(ng)
               track(i,j,l) = Fwrk(ct)
               ct=ct+1
          ENDDO
        ENDDO
      ENDDO

      !FwrkF=RESHAPE(bioenergy,(/NptsF/))
      !CALL mp_collect (ng, iNLM, NptsF, Fspv, FwrkF)
      !bioenergy=RESHAPE(FwrkF,(/NFishV(ng),Nfish(ng)/))

      ct=1
      DO l=1,Nfish(ng)
        DO j=1,NFishV(ng)
            FwrkF(ct) = bioenergy(j,l)
            ct=ct+1
        ENDDO
      ENDDO

      CALL mp_collect (ng, iNLM, NptsF, Fspv, FwrkF)

      ct=1
      DO l=1,Nfish(ng)
        DO j=1,NFishV(ng)
            bioenergy(j,l) = FwrkF(ct)
            ct=ct+1
        ENDDO
      ENDDO

      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, egg_dur)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, egg_num)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, ysac_dur)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, ysac_num)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, larv_dur)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, larv_num)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, juv_dur)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, juv_num)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, suba_num)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, fmortN)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, fmortS)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, fmortP)
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, fmortPsum)
      CALL mp_collect_i (ng, iNLM, Nfish(ng), iFspv, lifestage)

      !FwrkI=RESHAPE(swimtype,(/NptsI/))
      !CALL mp_collect_i (ng, iNLM, NptsI, iFspv, FwrkI)
      !swimtype=RESHAPE(FwrkI,(/2,Nfish(ng)/))

      ct=1
      DO l=1,Nfish(ng)
        DO j=1,2
            FwrkI(ct) = swimtype(j,l)
            ct=ct+1
        ENDDO
      ENDDO

      CALL mp_collect_i (ng, iNLM, NptsI, iFspv, FwrkI)

      ct=1
      DO l=1,Nfish(ng)
        DO j=1,2
            swimtype(j,l) = FwrkI(ct)
            ct=ct+1
        ENDDO
      ENDDO

      CALL mp_collect_i (ng, iNLM, NptsI, iFspv,swimtype)

      CALL mp_collect_i (ng, iNLM, Nfish(ng), iFspv, deathby)
      CALL mp_collect_i (ng, iNLM, Nfish(ng), iFspv, cellid)
      CALL mp_collect_l (ng, iNLM, Nfish(ng), alive)
!
#ifdef FOO
!  Collect the alive status.
      Fwrk=Fspv
      DO l=1,Nfish(ng)
        IF (bounded(l)) THEN
          Fwrk(l)=1.0_r8
        END IF
      END DO
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, Fwrk)
      DO l=1,Nfish(ng)
        IF (Fwrk(l).ne.Fspv) THEN
          alive(l)=.TRUE.
        ELSE
          alive(l)=.FALSE.
        END IF
      END DO
#endif
!
!  Collect the bounded status switch.
      Fwrk=Fspv
      DO l=1,Nfish(ng)
        IF (bounded(l)) THEN
          Fwrk(l)=1.0_r8
        END IF
      END DO
      CALL mp_collect (ng, iNLM, Nfish(ng), Fspv, Fwrk)
      DO l=1,Nfish(ng)
        IF (Fwrk(l).ne.Fspv) THEN
          bounded(l)=.TRUE.
        ELSE
          bounded(l)=.FALSE.
        END IF
      END DO
!
# endif
      RETURN
      END SUBROUTINE step_fish_tile
#endif
      END MODULE step_fish_mod
