#include "cppdefs.h"
#ifdef NEMURO_SAN
      SUBROUTINE ini_fish
!
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes the fish.                                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_fish
      USE mod_biology
      USE mod_grid
      USE nrutil
# ifdef DISTRIBUTE
      USE mod_strings
# endif
!
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcastf, mp_bcasts
# endif
      USE ran_state, ONLY: ran_seed
!
      implicit none
!
!  Local variable declarations.
!
      logical :: Lwrite
      logical :: lflag
      logical :: landcell

      integer :: Itile, Jtile, Nghost, Ntiles, tile
      integer :: Imin, Imax, Jmin, Jmax
# ifdef DISTRIBUTE
      integer :: MaxHaloLenI, MaxHaloLenJ
# endif
      integer :: inp, out, itrc, ng, npts, sequence
      integer :: i, isp, ils, l, fid, j

      real(r8) :: cff, Fwwght0, Fworth0, snudg
      real(r8) :: LfromW, matchk
      real(r8) :: xini, yini, zini
      real(r8), parameter :: spv = 0.0_r8
!
!-----------------------------------------------------------------------
!  Read in and report input model parameters.
!-----------------------------------------------------------------------
!
!  Set input units.
!
# if defined DISTRIBUTE || defined MODEL_COUPLING
      Lwrite=Master
      inp=1
      out=stdout
# else
      Lwrite=Master
      inp=stdinp
      out=stdout
# endif
!
!-----------------------------------------------------------------------
!  Read in fish input parameters.
!-----------------------------------------------------------------------
!
      OPEN (46, FILE=TRIM(fishposnam), FORM='formatted', STATUS='old')

      CALL read_FishPar (46, out, Lwrite)

!  Don't have these things for each grid... yet.
      ng = 1
!
!  Fill in parameters that are lifestage dependent.
      DO isp=1,Nspecies(ng)
        DO ils=if_egg,if_adult
          IF (ils.eq.if_larva) THEN
            ZSpref(ils,isp,ng)=ZSpref_L(isp,ng)
            ZLpref(ils,isp,ng)=ZLpref_L(isp,ng)
            ZPpref(ils,isp,ng)=ZPpref_L(isp,ng)
            K_ZS(ils,isp,ng)=K_ZS_L(isp,ng)
            K_ZL(ils,isp,ng)=K_ZL_L(isp,ng)
            K_ZP(ils,isp,ng)=K_ZP_L(isp,ng)
            a_C(ils,isp,ng)=a_C_L(isp,ng)
            b_C(ils,isp,ng)=b_C_L(isp,ng)
            pvalmax(ils,isp,ng)=pvalmax_L(isp,ng)
            a_R(ils,isp,ng)=a_R_L(isp,ng)
            b_R(ils,isp,ng)=b_R_L(isp,ng)
            activity(ils,isp,ng)=activity_L(isp,ng)
            d_R(ils,isp,ng)=d_R_L(isp,ng)
            Fswim(ils,isp,ng)=Fswim_L(isp,ng)
            a_AE(ils,isp,ng)=a_AE_L(isp,ng)
            b_AE(ils,isp,ng)=b_AE_L(isp,ng)
            AEmax(ils,isp,ng)=AEmax_L(isp,ng)
            te1(ils,isp,ng)=te1_L(isp,ng)
            te2(ils,isp,ng)=te2_L(isp,ng)
            te3(ils,isp,ng)=te3_L(isp,ng)
            te4(ils,isp,ng)=te4_L(isp,ng)
            xk1(ils,isp,ng)=xk1_L(isp,ng)
            xk2(ils,isp,ng)=xk2_L(isp,ng)
            xk3(ils,isp,ng)=xk3_L(isp,ng)
            xk4(ils,isp,ng)=xk4_L(isp,ng)
            cr(ils,isp,ng)=cr_L(isp,ng)
            tr(ils,isp,ng)=tr_L(isp,ng)
            aw2l(ils,isp,ng)=aw2l_L(isp,ng)
            bw2l(ils,isp,ng)=bw2l_L(isp,ng)
            dSLk(ils,isp,ng)=dSLk_L(isp,ng)
            dSLinf(ils,isp,ng)=dSLinf_L(isp,ng)
            al2w(ils,isp,ng)=al2w_L(isp,ng)
            bl2w(ils,isp,ng)=bl2w_L(isp,ng)
          ELSE IF (ils.eq.if_juvenile) THEN
            ZSpref(ils,isp,ng)=ZSpref_J(isp,ng)
            ZLpref(ils,isp,ng)=ZLpref_J(isp,ng)
            ZPpref(ils,isp,ng)=ZPpref_J(isp,ng)
            K_ZS(ils,isp,ng)=K_ZS_J(isp,ng)
            K_ZL(ils,isp,ng)=K_ZL_J(isp,ng)
            K_ZP(ils,isp,ng)=K_ZP_J(isp,ng)
            a_C(ils,isp,ng)=a_C_J(isp,ng)
            b_C(ils,isp,ng)=b_C_J(isp,ng)
            pvalmax(ils,isp,ng)=pvalmax_J(isp,ng)
            a_R(ils,isp,ng)=a_R_J(isp,ng)
            b_R(ils,isp,ng)=b_R_J(isp,ng)
            activity(ils,isp,ng)=activity_J(isp,ng)
            d_R(ils,isp,ng)=d_R_J(isp,ng)
            Fswim(ils,isp,ng)=Fswim_J(isp,ng)
            a_AE(ils,isp,ng)=a_AE_J(isp,ng)
            b_AE(ils,isp,ng)=b_AE_J(isp,ng)
            AEmax(ils,isp,ng)=AEmax_J(isp,ng)
            te1(ils,isp,ng)=te1_J(isp,ng)
            te2(ils,isp,ng)=te2_J(isp,ng)
            te3(ils,isp,ng)=te3_J(isp,ng)
            te4(ils,isp,ng)=te4_J(isp,ng)
            xk1(ils,isp,ng)=xk1_J(isp,ng)
            xk2(ils,isp,ng)=xk2_J(isp,ng)
            xk3(ils,isp,ng)=xk3_J(isp,ng)
            xk4(ils,isp,ng)=xk4_J(isp,ng)
            cr(ils,isp,ng)=cr_J(isp,ng)
            tr(ils,isp,ng)=tr_J(isp,ng)
            aw2l(ils,isp,ng)=aw2l_J(isp,ng)
            bw2l(ils,isp,ng)=bw2l_J(isp,ng)
            dSLk(ils,isp,ng)=dSLk_J(isp,ng)
            dSLinf(ils,isp,ng)=dSLinf_J(isp,ng)
            al2w(ils,isp,ng)=al2w_J(isp,ng)
            bl2w(ils,isp,ng)=bl2w_J(isp,ng)
          ELSE IF ((ils.eq.if_subadult).or.(ils.eq.if_adult)) THEN
            ZSpref(ils,isp,ng)=ZSpref_A(isp,ng)
            ZLpref(ils,isp,ng)=ZLpref_A(isp,ng)
            ZPpref(ils,isp,ng)=ZPpref_A(isp,ng)
            K_ZS(ils,isp,ng)=K_ZS_A(isp,ng)
            K_ZL(ils,isp,ng)=K_ZL_A(isp,ng)
            K_ZP(ils,isp,ng)=K_ZP_A(isp,ng)
            a_C(ils,isp,ng)=a_C_A(isp,ng)
            b_C(ils,isp,ng)=b_C_A(isp,ng)
            pvalmax(ils,isp,ng)=pvalmax_A(isp,ng)
            a_R(ils,isp,ng)=a_R_A(isp,ng)
            b_R(ils,isp,ng)=b_R_A(isp,ng)
            activity(ils,isp,ng)=activity_A(isp,ng)
            d_R(ils,isp,ng)=d_R_A(isp,ng)
            Fswim(ils,isp,ng)=Fswim_A(isp,ng)
            a_AE(ils,isp,ng)=a_AE_A(isp,ng)
            b_AE(ils,isp,ng)=b_AE_A(isp,ng)
            AEmax(ils,isp,ng)=AEmax_A(isp,ng)
            te1(ils,isp,ng)=te1_A(isp,ng)
            te2(ils,isp,ng)=te2_A(isp,ng)
            te3(ils,isp,ng)=te3_A(isp,ng)
            te4(ils,isp,ng)=te4_A(isp,ng)
            xk1(ils,isp,ng)=xk1_A(isp,ng)
            xk2(ils,isp,ng)=xk2_A(isp,ng)
            xk3(ils,isp,ng)=xk3_A(isp,ng)
            xk4(ils,isp,ng)=xk4_A(isp,ng)
            cr(ils,isp,ng)=cr_A(isp,ng)
            tr(ils,isp,ng)=tr_A(isp,ng)
            aw2l(ils,isp,ng)=aw2l_A(isp,ng)
            bw2l(ils,isp,ng)=bw2l_A(isp,ng)
            dSLk(ils,isp,ng)=dSLk_A(isp,ng)
            dSLinf(ils,isp,ng)=dSLinf_A(isp,ng)
            al2w(ils,isp,ng)=al2w_A(isp,ng)
            bl2w(ils,isp,ng)=bl2w_A(isp,ng)
          ELSE
            ZSpref(ils,isp,ng)=0.0_r8
            ZLpref(ils,isp,ng)=0.0_r8
            ZPpref(ils,isp,ng)=0.0_r8
            K_ZS(ils,isp,ng)=0.0_r8
            K_ZL(ils,isp,ng)=0.0_r8
            K_ZP(ils,isp,ng)=0.0_r8
            a_C(ils,isp,ng)=0.0_r8
            b_C(ils,isp,ng)=0.0_r8
            pvalmax(ils,isp,ng)=0.0_r8
            a_R(ils,isp,ng)=0.0_r8
            b_R(ils,isp,ng)=0.0_r8
            activity(ils,isp,ng)=0.0_r8
            d_R(ils,isp,ng)=0.0_r8
            Fswim(ils,isp,ng)=0.0_r8
            a_AE(ils,isp,ng)=0.0_r8
            b_AE(ils,isp,ng)=0.0_r8
            AEmax(ils,isp,ng)=0.0_r8
            te1(ils,isp,ng)=0.0_r8
            te2(ils,isp,ng)=0.0_r8
            te3(ils,isp,ng)=0.0_r8
            te4(ils,isp,ng)=0.0_r8
            xk1(ils,isp,ng)=0.0_r8
            xk2(ils,isp,ng)=0.0_r8
            xk3(ils,isp,ng)=0.0_r8
            xk4(ils,isp,ng)=0.0_r8
            cr(ils,isp,ng)=0.0_r8
            tr(ils,isp,ng)=0.0_r8
            aw2l(ils,isp,ng)=0.0_r8
            bw2l(ils,isp,ng)=0.0_r8
            dSLk(ils,isp,ng)=0.0_r8
            dSLinf(ils,isp,ng)=0.0_r8
            al2w(ils,isp,ng)=0.0_r8
            bl2w(ils,isp,ng)=0.0_r8
          END IF
        END DO
      END DO
!
! Fill out stage dependent natural mortality
      DO isp=1,Nspecies(ng)
        DO ils=if_egg,if_adult
          IF (ils.eq.if_egg) THEN
            Nmort(ils,isp,ng)=Nmort_E(isp,ng)
          ELSE IF (ils.eq.if_yolksac) THEN
            Nmort(ils,isp,ng)=Nmort_Y(isp,ng)
          ELSE IF (ils.eq.if_larva) THEN
            Nmort(ils,isp,ng)=Nmort_L(isp,ng)
          ELSE IF (ils.eq.if_juvenile) THEN
            Nmort(ils,isp,ng)=Nmort_J(isp,ng)
          ELSE IF ((ils.eq.if_subadult).or.(ils.eq.if_adult)) THEN
            Nmort(ils,isp,ng)=Nmort_A(isp,ng)
          END IF
        END DO
      END DO

# ifdef DISTRIBUTE
      IF (Master) THEN
        CALL ran1 (FISHES(ng)%r2walk)
      END IF
      CALL mp_bcastf (ng, iNLM, FISHES(ng)%r2walk)
# elif defined _OPENMP
!$OMP SINGLE
      CALL ran1 (FISHES(ng)%r2walk)
!$OMP END SINGLE
# else
!     IF (Lstr.eq.1) THEN
        CALL ran1 (FISHES(ng)%r2walk)
!     END IF
# endif
!
!  KLUDGE WARNING!
!      idfish(if_herring) = 1
!      idfish_inv(1) = if_herring
      idfish(if_anchovy) = 1
      idfish_inv(1) = if_anchovy
      idfish(if_sardine) = 2
      idfish_inv(2) = if_sardine

      DO isp=1,Nspecies(ng)
        FISHES(ng)%species_base(isp)=(isp-1)*Nyearclass(ng)*            &
     &                                       Nfishperyear(ng)
        FISHES(ng)%age_base(isp)=FISHES(ng)%species_base(isp)
      END DO

      DO isp=1,Nspecies(ng)
        FISHES(ng)%num_free(isp)=Nfishperyear(ng)
        FISHES(ng)%num_super(isp) = FISHES(ng)%num_free(isp) /          &
     &        (Fspend(isp,ng)-Fspstr(isp,ng) + 1)
        FISHES(ng)%next_free(isp)=1+FISHES(ng)%age_base(isp)
      END DO
      IF (Master) THEN
        print *, 'First Year for Fish'
        print *, 'Free fish ', FISHES(ng) % num_free
        print *, 'Fish ages ', FISHES(ng) % age_base
      END IF

!  Put the new fish in here, all except age 0.
      DO isp=1,Nspecies(ng)
        DO i=1,Nyearclass(ng)-1
! NB: initialization for fish wet weight loosely adapted from KR anchovy
!          Fwwght0=MIN(35.0_r8,Fwwt0(isp,ng)+REAL(i-1,r8)*5.0_r8)
!          Fwwght0=MIN(3.5_r8*Fwwt0(isp,ng),                            &
!     &                Fwwt0(isp,ng)*(1.0_r8+0.5_r8*REAL(i-1,r8)))
          IF (isp.eq.if_anchovy) THEN
            Fwwght0=MIN(2.5_r8*Fwwt0(isp,ng),Fwwt0(isp,ng)*            &
     &                  (1.0_r8+0.4_r8*REAL(Nyearclass(ng)-i-1,r8)))
          ELSE IF (isp.eq.if_sardine) THEN
            Fwwght0=MIN(4.0_r8*Fwwt0(isp,ng),Fwwt0(isp,ng)*            &
     &                  (1.0_r8+0.6_r8*REAL(Nyearclass(ng)-i-1,r8)))
          END IF
          Fworth0=Fwth0(isp,ng)*1.0e6_r8/REAL(Nfishperyear(ng),r8)*    &
     &              EXP(-Nmort_A(isp,ng)*REAL(Nyearclass(ng)-i-1,r8))
          DO l=1,Nfishperyear(ng)
            fid=FISHES(ng)%species_base(isp)+i*Nfishperyear(ng)+l
# ifdef CCS30
! ---- FOR CCS CONFIGURATION ----
            CALL ran1 (snudg)
            IF (snudg.lt.0.5_r8) THEN
! southern domain
              yini=20.0_r8+2.0_r8*snudg*25_r8
              IF (isp.eq.if_anchovy) THEN
! for anchovy (inshore)
                xini=31.0_r8
              ELSE IF (isp.eq.if_sardine) THEN
! for sardine (offshore)
                xini=30.0_r8
              END IF
            ELSE
              yini=55.0_r8+2.0_r8*(snudg-0.5_r8)*25_r8
              IF (isp.eq.if_anchovy) THEN
! for anchovy (inshore)
                xini=27.0_r8
              ELSE IF (isp.eq.if_sardine) THEN
! for sardine (offshore)
                xini=26.0_r8
              END IF
            END IF
            CALL ran1 (snudg)
            xini=xini+2.0_r8*(snudg-0.5_r8)
            zini=REAL(N(ng),r8)
# endif
# ifdef WC13
! ---- FOR WC CONFIGURATION ----
            CALL ran1 (snudg)
! old            yini=MAX(5.0_r8,snudg*REAL(Mm(ng)-4,r8))
! old            IF (yini.lt.33.0_r8) THEN
            yini=10.0_r8+snudg*REAL(Mm(ng)-20,r8)
            IF (yini.lt.29.0_r8) THEN
              IF (isp.eq.if_anchovy) THEN
! for anchovy (inshore)
! old                xini=47.0_r8-(yini-3.0_r8)*25.0_r8/30.0_r8
                xini=40.0_r8-(yini-10.0_r8)*15.0_r8/19.0_r8
              ELSE IF (isp.eq.if_sardine) THEN
! for sardine (offshore)
! old                xini=44.0_r8-(yini-3.0_r8)*25.0_r8/30.0_r8
                xini=37.0_r8-(yini-10.0_r8)*15.0_r8/19.0_r8
              END IF
            ELSE
              IF (isp.eq.if_anchovy) THEN
! for anchovy (inshore)
! old                xini=22.0_r8
                xini=25.0_r8
              ELSE IF (isp.eq.if_sardine) THEN
! for sardine (offshore)
! old                xini=19.0_r8
                xini=22.0_r8
              END IF
            END IF
            CALL ran1 (snudg)
            xini=xini+2.0_r8*(snudg-0.5_r8)
            zini=REAL(N(ng),r8)
# endif
# ifdef WC12
! ---- FOR WC CONFIGURATION ----
            CALL ran1 (snudg)
            yini=40.0_r8+snudg*REAL(Mm(ng)-70,r8)
            IF (yini.lt.100.0_r8) THEN
              IF (isp.eq.if_anchovy) THEN
! for anchovy (inshore)
                xini=92.0_r8+(100-yini)*0.7_r8
              ELSE IF (isp.eq.if_sardine) THEN
! for sardine (offshore)
                xini=85.0_r8+(100-yini)*0.7_r8
              END IF
            ELSE
              IF (isp.eq.if_anchovy) THEN
! for anchovy (inshore)
                xini=92.0_r8
              ELSE IF (isp.eq.if_sardine) THEN
! for sardine (offshore)
! old                xini=19.0_r8
                xini=85.0_r8
              END IF
            END IF
            CALL ran1 (snudg)
            xini=xini-4.0_r8*snudg
            zini=REAL(N(ng),r8)
# endif
# ifdef CCS1
! ----  FOR CCS1 CONFIGURATION : may need fixes !!! RD ----

            CALL ran1 (snudg)
            yini=220.0_r8+snudg*REAL(177,r8)
            IF (yini.lt.330.0_r8) THEN
              IF (isp.eq.if_anchovy) THEN
! for anchovy (inshore)
                xini=115.0_r8
              ELSE IF (isp.eq.if_sardine) THEN
! for sardine (offshore)
                xini=107.0_r8
              END IF
            ELSE
              IF (isp.eq.if_anchovy) THEN
! for anchovy (inshore)
                xini=115.0_r8+(yini-330)*0.57_r8
              ELSE IF (isp.eq.if_sardine) THEN
! for sardine (offshore)
                xini=107.0_r8+(yini-330)*0.57_r8
              END IF
            END IF
            CALL ran1 (snudg)
            xini=xini-4.0_r8*snudg
            zini=REAL(N(ng),r8)
# endif
# ifdef NWPACIFIC
!--- Shinichi custom fish initial position ---
            CALL ran1 (snudg)
            IF (isp.eq.if_sardine) THEN
! sardines (shinichi) :
! lowest_point_j=171.
! lowest_point_i=154.
! extension_j=100.
! mid_point1_j=258.
! mid_point1_i=197.
! upper_point_j=279
! upper_point_i=183
!
!               yini=171.0_r8+snudg*REAL(100,r8)
!               IF (yini .lt. 258.) THEN
!                  xini = 154.0_r8 + (yini - 171.0_r8) * (197.0_r8 - 154.0_r8) / (258.0_r8 - 171.0_r8)
!               ELSE
!                  xini = 197.0_r8 + (yini - 258.0_r8) * (183.0_r8 - 197.0_r8) / (279.0_r8 - 258.0_r8)
!               ENDIF

! sardines (takaaki) :
! make sure sardines are offshore, anchovies inshore
! lowest_point_j=177.
! lowest_point_i=162.
! extension_j=100.
! mid_point1_j=261.
! mid_point1_i=203.
! upper_point_j=282
! upper_point_i=190

               yini=177.0_r8+snudg*REAL(100,r8)
               IF (yini .lt. 261.) THEN
                  xini = 162.0_r8 + (yini - 177.0_r8) * (203.0_r8 - 162.0_r8) / (261.0_r8 - 177.0_r8)
               ELSE
                  xini = 203.0_r8 + (yini - 261.0_r8) * (190.0_r8 - 203.0_r8) / (282.0_r8 - 261.0_r8)
               ENDIF

            ELSEIF (isp.eq.if_anchovy) THEN

! anchovies (shinichi) : 
! lowest_point_j=166.
! lowest_point_i=152.
! extension_j=105.
! mid_point1_j=258.
! mid_point1_i=197.
! upper_point_j=279
! upper_point_i=183
!               yini=166.0_r8+snudg*REAL(105,r8)
!               IF (yini .lt. 258.) THEN
!                  xini = 152.0_r8 + (yini - 166.0_r8) * (197.0_r8 - 152.0_r8) / (258.0_r8 - 166.0_r8)
!               ELSE
!                  xini = 197.0_r8 + (yini - 258.0_r8) * (183.0_r8 - 197.0_r8) / (279.0_r8 - 258.0_r8)
!               ENDIF

! anchovies (takaaki) : 
! west_point_j=172.
! lowest_point_i=160.
! extension_j=105.
! mid_point1_j=261.
! mid_point1_i=203.
! upper_point_j=282
! upper_point_i=190

               yini=172.0_r8+snudg*REAL(105,r8)
               IF (yini .lt. 261.) THEN
                  xini = 160.0_r8 + (yini - 172.0_r8) * (203.0_r8 - 160.0_r8) / (261.0_r8 - 172.0_r8)
               ELSE
                  xini = 203.0_r8 + (yini - 261.0_r8) * (190.0_r8 - 203.0_r8) / (282.0_r8 - 261.0_r8)
               ENDIF

            ELSE
              STOP ! no other species
            ENDIF
            CALL ran1 (snudg)
            xini=xini-4.0_r8*snudg
            zini=REAL(N(ng),r8)

#endif

# ifdef DOUBLE_GYRE
! ---- FOR TEST CONFIGURATION ----
            CALL ran1 (snudg)
!            xini=4.0_r8+(snudg-0.5_r8)
!            CALL ran1 (snudg)
!            yini=7.0_r8+(snudg-0.5_r8)
!            xini=REAL(i+7,r8)+0.1_r8*(snudg-0.5_r8)
!            yini=REAL(i+7,r8)+0.1_r8*(snudg-0.5_r8)
            xini=13.0_r8
            yini=27.0_r8
            zini=REAL(N(ng),r8)
# endif
            FISHES(ng)%Tinfo(ixgrd,fid)=xini
            FISHES(ng)%Tinfo(iygrd,fid)=yini
            FISHES(ng)%Tinfo(izgrd,fid)=zini
            DO j=0,NFT
              FISHES(ng)%track(ixgrd,j,fid)=xini
              FISHES(ng)%track(iygrd,j,fid)=yini
#  ifdef SOLVE3D
              FISHES(ng)%track(izgrd,j,fid)=zini
              FISHES(ng)%track(ixrhs,j,fid)=0.0_r8
              FISHES(ng)%track(iyrhs,j,fid)=0.0_r8
              FISHES(ng)%track(izrhs,j,fid)=0.0_r8
#  endif
            END DO
!            FISHES(ng)%bioenergy(ifwwt,fid)=Fwwght0*(FISHES(ng)%        &
!     &               r2walk(fid+2*Nfish(ng))+0.5_r8)
            CALL ran1 (snudg)
            cff=1.0_r8+0.2_r8*(snudg-0.5_r8)
            FISHES(ng)%bioenergy(ifwwt,fid)=cff*Fwwght0
            FISHES(ng)%bioenergy(ifworth,fid)=Fworth0
            FISHES(ng)%bioenergy(ifage,fid)=REAL(Nyearclass(ng)-i,r8)
            FISHES(ng)%bioenergy(ifbday,fid)=REAL(Nyearclass(ng)-i,r8)
            FISHES(ng)%bioenergy(ifripe,fid)=0.0_r8
            FISHES(ng)%bioenergy(ifrest,fid)=0.0_r8
            FISHES(ng)%bioenergy(ifdevl,fid)=0.0_r8
            FISHES(ng)%bioenergy(ifatre,fid)=0.0_r8
            FISHES(ng)%bioenergy(ifgonj,fid)=0.0_r8
            FISHES(ng)%bioenergy(ifrcvr,fid)=0.0_r8
            FISHES(ng)%species(fid)=idfish_inv(isp)
            FISHES(ng)%alive(fid)=.TRUE.
            FISHES(ng)%swimtype(1,fid)=Hbehave(isp,ng)
            FISHES(ng)%swimtype(2,fid)=Vbehave(isp,ng)
            FISHES(ng)%species(fid)=idfish_inv(isp)
            FISHES(ng)%bounded(fid)=.TRUE.
            FISHES(ng)%cellid(fid)=nint(xini)+(Lm(ng)+2)*nint(yini)
            FISHES(ng)%deathby(fid)=0
! Set lifestage and fish length
            LfromW=(1.0_r8/aw2l(if_adult,isp,ng)*                       &
     &                     FISHES(ng)%bioenergy(ifwwt,fid))**           &
     &                         (1.0_r8/bw2l(if_adult,isp,ng))
            FISHES(ng)%bioenergy(iflngth,fid)=MIN(LfromW,               &
     &                           dSLinf(if_adult,isp,ng))
! Option 1: let spanwing routine decide if mature
!            FISHES(ng)%lifestage(fid)=if_subadult
! Option 2: calculate maturity from length
            matchk=1.0_r8/(1.0_r8+exp(-(amature(isp,ng)+                &
     &                            bmature(isp,ng)*LfromW)))
            IF (matchk.ge.0.5_r8) THEN
              FISHES(ng)%lifestage(fid)=if_adult
            ELSE
              FISHES(ng)%lifestage(fid)=if_subadult
            ENDIF 
          END DO
        END DO
      END DO

      END SUBROUTINE ini_fish

      SUBROUTINE read_FishPar (inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in input fish parameters.                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_fish
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval
      integer :: i, j, igrid, mc, nc, ng, status, isp

      integer, dimension(Ngrids) :: ncount, nentry
      integer, dimension(Ngrids) :: Npos

      integer, allocatable :: Fcoor(:,:), Fcount(:,:), Ftype(:,:)

      integer :: decode_line, load_i, load_l, load_r

      real(r8) :: xfloat, yfloat, zfloat

      real(r8), dimension(100) :: Rval

      real(r8), allocatable :: Ft0(:,:), Fx0(:,:), Fy0(:,:), Fz0(:,:)
      real(r8), allocatable :: Fdt(:,:), Fdx(:,:), Fdy(:,:), Fdz(:,:)

      character (len=35) :: frmt
      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in initial fish locations.
!-----------------------------------------------------------------------
!
!  Notice I added one when allocating local scratch arrays to avoid
!  out of bounds in some compilers when reading the last blank line
!  which signal termination of input data.
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=20,END=30) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lfish') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lfish)
          ELSE IF (TRIM(KeyWord).eq.'Nspecies') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nspecies)
          ELSE IF (TRIM(KeyWord).eq.'Nyearclass') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nyearclass)
          ELSE IF (TRIM(KeyWord).eq.'Nfishperyear') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nfishperyear)
          ELSE IF (TRIM(KeyWord).eq.'Fwwt0') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Fwwt0)
          ELSE IF (TRIM(KeyWord).eq.'Fwth0') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Fwth0)
          ELSE IF (TRIM(KeyWord).eq.'Fspstr') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Fspstr)
          ELSE IF (TRIM(KeyWord).eq.'Fspend') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Fspend)
          ELSE IF (TRIM(KeyWord).eq.'FspTmin') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, FspTmin)
          ELSE IF (TRIM(KeyWord).eq.'FspTmax') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, FspTmax)
          ELSE IF (TRIM(KeyWord).eq.'Hbehave') THEN
            Npts=load_i(Nval, Rval, max_species*Ngrids, Hbehave)
          ELSE IF (TRIM(KeyWord).eq.'Vbehave') THEN
            Npts=load_i(Nval, Rval, max_species*Ngrids, Vbehave)
          ELSE IF (TRIM(KeyWord).eq.'ZSpref_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, ZSpref_L)
          ELSE IF (TRIM(KeyWord).eq.'ZSpref_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, ZSpref_J)
          ELSE IF (TRIM(KeyWord).eq.'ZSpref_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, ZSpref_A)
          ELSE IF (TRIM(KeyWord).eq.'ZLpref_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, ZLpref_L)
          ELSE IF (TRIM(KeyWord).eq.'ZLpref_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, ZLpref_J)
          ELSE IF (TRIM(KeyWord).eq.'ZLpref_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, ZLpref_A)
          ELSE IF (TRIM(KeyWord).eq.'ZPpref_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, ZPpref_L)
          ELSE IF (TRIM(KeyWord).eq.'ZPpref_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, ZPpref_J)
          ELSE IF (TRIM(KeyWord).eq.'ZPpref_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, ZPpref_A)
          ELSE IF (TRIM(KeyWord).eq.'K_ZS_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, K_ZS_L)
          ELSE IF (TRIM(KeyWord).eq.'K_ZS_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, K_ZS_J)
          ELSE IF (TRIM(KeyWord).eq.'K_ZS_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, K_ZS_A)
          ELSE IF (TRIM(KeyWord).eq.'K_ZL_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, K_ZL_L)
          ELSE IF (TRIM(KeyWord).eq.'K_ZL_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, K_ZL_J)
          ELSE IF (TRIM(KeyWord).eq.'K_ZL_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, K_ZL_A)
          ELSE IF (TRIM(KeyWord).eq.'K_ZP_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, K_ZP_L)
          ELSE IF (TRIM(KeyWord).eq.'K_ZP_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, K_ZP_J)
          ELSE IF (TRIM(KeyWord).eq.'K_ZP_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, K_ZP_A)
          ELSE IF (TRIM(KeyWord).eq.'Cal_Z') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Cal_Z)
          ELSE IF (TRIM(KeyWord).eq.'Cal_F') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Cal_F)
          ELSE IF (TRIM(KeyWord).eq.'a_C_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, a_C_L)
          ELSE IF (TRIM(KeyWord).eq.'a_C_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, a_C_J)
          ELSE IF (TRIM(KeyWord).eq.'a_C_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, a_C_A)
          ELSE IF (TRIM(KeyWord).eq.'b_C_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, b_C_L)
          ELSE IF (TRIM(KeyWord).eq.'b_C_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, b_C_J)
          ELSE IF (TRIM(KeyWord).eq.'b_C_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, b_C_A)
          ELSE IF (TRIM(KeyWord).eq.'pvalmax_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, pvalmax_L)
          ELSE IF (TRIM(KeyWord).eq.'pvalmax_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, pvalmax_J)
          ELSE IF (TRIM(KeyWord).eq.'pvalmax_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, pvalmax_A)
          ELSE IF (TRIM(KeyWord).eq.'a_R_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, a_R_L)
          ELSE IF (TRIM(KeyWord).eq.'a_R_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, a_R_J)
          ELSE IF (TRIM(KeyWord).eq.'a_R_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, a_R_A)
          ELSE IF (TRIM(KeyWord).eq.'b_R_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, b_R_L)
          ELSE IF (TRIM(KeyWord).eq.'b_R_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, b_R_J)
          ELSE IF (TRIM(KeyWord).eq.'b_R_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, b_R_A)
          ELSE IF (TRIM(KeyWord).eq.'activity_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, activity_L)
          ELSE IF (TRIM(KeyWord).eq.'activity_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, activity_J)
          ELSE IF (TRIM(KeyWord).eq.'activity_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, activity_A)
          ELSE IF (TRIM(KeyWord).eq.'d_R_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, d_R_L)
          ELSE IF (TRIM(KeyWord).eq.'d_R_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, d_R_J)
          ELSE IF (TRIM(KeyWord).eq.'d_R_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, d_R_A)
          ELSE IF (TRIM(KeyWord).eq.'Fswim_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Fswim_L)
          ELSE IF (TRIM(KeyWord).eq.'Fswim_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Fswim_J)
          ELSE IF (TRIM(KeyWord).eq.'Fswim_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Fswim_A)
          ELSE IF (TRIM(KeyWord).eq.'a_AE_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, a_AE_L)
          ELSE IF (TRIM(KeyWord).eq.'a_AE_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, a_AE_J)
          ELSE IF (TRIM(KeyWord).eq.'a_AE_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, a_AE_A)
          ELSE IF (TRIM(KeyWord).eq.'b_AE_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, b_AE_L)
          ELSE IF (TRIM(KeyWord).eq.'b_AE_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, b_AE_J)
          ELSE IF (TRIM(KeyWord).eq.'b_AE_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, b_AE_A)
          ELSE IF (TRIM(KeyWord).eq.'AEmax_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, AEmax_L)
          ELSE IF (TRIM(KeyWord).eq.'AEmax_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, AEmax_J)
          ELSE IF (TRIM(KeyWord).eq.'AEmax_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, AEmax_A)
          ELSE IF (TRIM(KeyWord).eq.'te1_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, te1_L)
          ELSE IF (TRIM(KeyWord).eq.'te1_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, te1_J)
          ELSE IF (TRIM(KeyWord).eq.'te1_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, te1_A)
          ELSE IF (TRIM(KeyWord).eq.'te2_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, te2_L)
          ELSE IF (TRIM(KeyWord).eq.'te2_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, te2_J)
          ELSE IF (TRIM(KeyWord).eq.'te2_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, te2_A)
          ELSE IF (TRIM(KeyWord).eq.'te3_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, te3_L)
          ELSE IF (TRIM(KeyWord).eq.'te3_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, te3_J)
          ELSE IF (TRIM(KeyWord).eq.'te3_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, te3_A)
          ELSE IF (TRIM(KeyWord).eq.'te4_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, te4_L)
          ELSE IF (TRIM(KeyWord).eq.'te4_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, te4_J)
          ELSE IF (TRIM(KeyWord).eq.'te4_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, te4_A)
          ELSE IF (TRIM(KeyWord).eq.'xk1_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, xk1_L)
          ELSE IF (TRIM(KeyWord).eq.'xk1_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, xk1_J)
          ELSE IF (TRIM(KeyWord).eq.'xk1_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, xk1_A)
          ELSE IF (TRIM(KeyWord).eq.'xk2_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, xk2_L)
          ELSE IF (TRIM(KeyWord).eq.'xk2_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, xk2_J)
          ELSE IF (TRIM(KeyWord).eq.'xk2_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, xk2_A)
          ELSE IF (TRIM(KeyWord).eq.'xk3_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, xk3_L)
          ELSE IF (TRIM(KeyWord).eq.'xk3_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, xk3_J)
          ELSE IF (TRIM(KeyWord).eq.'xk3_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, xk3_A)
          ELSE IF (TRIM(KeyWord).eq.'xk4_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, xk4_L)
          ELSE IF (TRIM(KeyWord).eq.'xk4_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, xk4_J)
          ELSE IF (TRIM(KeyWord).eq.'xk4_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, xk4_A)
          ELSE IF (TRIM(KeyWord).eq.'cr_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, cr_L)
          ELSE IF (TRIM(KeyWord).eq.'cr_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, cr_J)
          ELSE IF (TRIM(KeyWord).eq.'cr_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, cr_A)
          ELSE IF (TRIM(KeyWord).eq.'tr_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, tr_L)
          ELSE IF (TRIM(KeyWord).eq.'tr_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, tr_J)
          ELSE IF (TRIM(KeyWord).eq.'tr_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, tr_A)
          ELSE IF (TRIM(KeyWord).eq.'Wffeed') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Wffeed)
          ELSE IF (TRIM(KeyWord).eq.'Lffeed') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Lffeed)
          ELSE IF (TRIM(KeyWord).eq.'WeightLJ') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, WeightLJ)
          ELSE IF (TRIM(KeyWord).eq.'LengthLJ') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, LengthLJ)
          ELSE IF (TRIM(KeyWord).eq.'WeightJA') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, WeightJA)
          ELSE IF (TRIM(KeyWord).eq.'LengthJA') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, LengthJA)
          ELSE IF (TRIM(KeyWord).eq.'aw2l_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, aw2l_L)
          ELSE IF (TRIM(KeyWord).eq.'aw2l_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, aw2l_J)
          ELSE IF (TRIM(KeyWord).eq.'aw2l_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, aw2l_A)
          ELSE IF (TRIM(KeyWord).eq.'bw2l_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, bw2l_L)
          ELSE IF (TRIM(KeyWord).eq.'bw2l_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, bw2l_J)
          ELSE IF (TRIM(KeyWord).eq.'bw2l_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, bw2l_A)
          ELSE IF (TRIM(KeyWord).eq.'dSLk_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, dSLk_L)
          ELSE IF (TRIM(KeyWord).eq.'dSLk_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, dSLk_J)
          ELSE IF (TRIM(KeyWord).eq.'dSLk_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, dSLk_A)
          ELSE IF (TRIM(KeyWord).eq.'dSLinf_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, dSLinf_L)
          ELSE IF (TRIM(KeyWord).eq.'dSLinf_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, dSLinf_J)
          ELSE IF (TRIM(KeyWord).eq.'dSLinf_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, dSLinf_A)
          ELSE IF (TRIM(KeyWord).eq.'al2w_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, al2w_L)
          ELSE IF (TRIM(KeyWord).eq.'al2w_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, al2w_J)
          ELSE IF (TRIM(KeyWord).eq.'al2w_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, al2w_A)
          ELSE IF (TRIM(KeyWord).eq.'bl2w_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, bl2w_L)
          ELSE IF (TRIM(KeyWord).eq.'bl2w_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, bl2w_J)
          ELSE IF (TRIM(KeyWord).eq.'bl2w_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, bl2w_A)
          ELSE IF (TRIM(KeyWord).eq.'abatch') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, abatch)
          ELSE IF (TRIM(KeyWord).eq.'bbatch') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, bbatch)
          ELSE IF (TRIM(KeyWord).eq.'T0batch') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, T0batch)
          ELSE IF (TRIM(KeyWord).eq.'apof') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, apof)
          ELSE IF (TRIM(KeyWord).eq.'bpof') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, bpof)
          ELSE IF (TRIM(KeyWord).eq.'T0pof') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, T0pof)
          ELSE IF (TRIM(KeyWord).eq.'epg') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, epg)
          ELSE IF (TRIM(KeyWord).eq.'eegg') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, eegg)
          ELSE IF (TRIM(KeyWord).eq.'megg') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, megg)
          ELSE IF (TRIM(KeyWord).eq.'breed') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, breed)
          ELSE IF (TRIM(KeyWord).eq.'amature') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, amature)
          ELSE IF (TRIM(KeyWord).eq.'bmature') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, bmature)
          ELSE IF (TRIM(KeyWord).eq.'pctxwt') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, pctxwt)
          ELSE IF (TRIM(KeyWord).eq.'pctgain') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, pctgain)
          ELSE IF (TRIM(KeyWord).eq.'Nmort_E') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Nmort_E)
          ELSE IF (TRIM(KeyWord).eq.'Nmort_Y') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Nmort_Y)
          ELSE IF (TRIM(KeyWord).eq.'Nmort_L') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Nmort_L)
          ELSE IF (TRIM(KeyWord).eq.'Nmort_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Nmort_J)
          ELSE IF (TRIM(KeyWord).eq.'Nmort_A') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Nmort_A)
!          ELSE IF (TRIM(KeyWord).eq.'Nymort') THEN
!            Npts=load_r(Nval, Rval, max_species*Ngrids, Nymort)
          ELSE IF (TRIM(KeyWord).eq.'Fymort') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Fymort)
          ELSE IF (TRIM(KeyWord).eq.'DDmort1_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, DDmort1_J)
          ELSE IF (TRIM(KeyWord).eq.'DDmort2_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, DDmort2_J)
          ELSE IF (TRIM(KeyWord).eq.'DDmort3_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, DDmort3_J)
          ELSE IF (TRIM(KeyWord).eq.'DDmort4_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, DDmort4_J)
          ELSE IF (TRIM(KeyWord).eq.'DDscale_J') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, DDscale_J)
          END IF
        END IF
      END DO
  20  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  30  CONTINUE
!
!-----------------------------------------------------------------------
!  Set initial fish locations
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        Nfish(ng)=Nfishperyear(ng)*Nspecies(ng)*Nyearclass(ng)
        Npos(ng)=Nfish(ng)
      END DO
      Npts=Nfish(1)
      IF (Ngrids.gt.1) Npts=MAXVAL(Nfish)
      allocate ( Fcoor (Npts,Ngrids) )
      allocate ( Fcount(Npts,Ngrids) )
      allocate ( Ftype (Npts,Ngrids) )
      allocate ( Ft0(Npts,Ngrids) )
      allocate ( Fx0(Npts,Ngrids) )
      allocate ( Fy0(Npts,Ngrids) )
      allocate ( Fz0(Npts,Ngrids) )
      allocate ( Fdt(Npts,Ngrids) )
      allocate ( Fdx(Npts,Ngrids) )
      allocate ( Fdy(Npts,Ngrids) )
      allocate ( Fdz(Npts,Ngrids) )
      DO ng=1,Ngrids
        CALL allocate_fish (ng)
      END DO
      ncount(1:Ngrids)=0
      nentry(1:Ngrids)=0
      DO i=1,Npts
        igrid=1
        Fcoor (nentry(igrid)+1,igrid)=0
        Ftype (nentry(igrid)+1,igrid)=1         ! 3-D Largrangian fish (neutrally buoyant)
!        Ftype (nentry(igrid)+1,igrid)=3        ! 3-D Geopotential fish (fixed depth)
        Fcount(nentry(igrid)+1,igrid)=1
        Ft0(nentry(igrid)+1,igrid)=0.0_r8
        Fx0(nentry(igrid)+1,igrid)=0.0_r8
        Fy0(nentry(igrid)+1,igrid)=0.0_r8
        Fz0(nentry(igrid)+1,igrid)=-5.0_r8
        Fdt(nentry(igrid)+1,igrid)=0.0_r8
        Fdx(nentry(igrid)+1,igrid)=0.0_r8
        Fdy(nentry(igrid)+1,igrid)=0.0_r8
        Fdz(nentry(igrid)+1,igrid)=0.0_r8
        ncount(igrid)=ncount(igrid)+Fcount(nentry(igrid)+1,igrid)
        nentry(igrid)=nentry(igrid)+1        
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (ncount(ng).ne.Npos(ng)) THEN
            IF (Master) WRITE (stdout,60) ncount(ng), Npos(ng)
            exit_flag=4
            RETURN
          END IF
          IF (Lfish(ng)) THEN
            WRITE (out,70) ng
            DO i=1,nentry(ng)
              IF (.not.spherical.and.(Fcoor(i,ng).eq.0)) THEN
                frmt='(i1,i2,i5,f10.4,2f8.2,f8.2,4f9.3)'
              ELSE
                frmt='(i1,i2,i5,f10.4,3f8.2,4f9.3)'
              END IF
              WRITE (out,frmt) Fcoor(i,ng), Ftype(i,ng), Fcount(i,ng),  &
     &                         Ft0(i,ng), Fx0(i,ng), Fy0(i,ng),         &
     &                         Fz0(i,ng), Fdt(i,ng), Fdx(i,ng),         &
     &                         Fdy(i,ng), Fdz(i,ng)
            END DO
            WRITE (out,80) Nfish(ng),                                   &
     &            'Nfish',                                              &
     &            'Number of fish trajectories to compute.'
            WRITE (out,80) Nfishperyear(ng),                            &
     &            'Nfishperyear',                                       &
     &            'Number of fish per year per species.'
            WRITE (out,80) Nspecies(ng),                                &
     &            'Nspecies',                                           &
     &            'Number of fish species.'
            WRITE (out,80) Nyearclass(ng),                              &
     &            'Nyearclass',                                         &
     &            'Number of year classes per species.'
            WRITE (out,80) Npos(ng),                                    &
     &            'Npos',                                               &
     &            'Number of fish initial positions.'
            DO isp=1,Nspecies(ng)
              WRITE (out,*)
              WRITE (out,*) 'PARAMETERS FOR SPECIES ',isp
              WRITE (out,*) '--------------------------------------'
              WRITE (out,110) Fwwt0(isp,ng), 'Fwwt0',                   &
     &                'Fish initial wet weight (g).'
              WRITE (out,110) Fwth0(isp,ng), 'Fwth0',                   &
     &                'Fish initial worth (millions of fishes).'
              WRITE (out,110) Fspstr(isp,ng), 'Fspstr',                 &
     &              'Fish start of spawning season (yearday).'
              WRITE (out,110) Fspend(isp,ng), 'Fspend',                 &
     &              'Fish end of spawning season (yearday).'
              WRITE (out,110) FspTmin(isp,ng), 'FspTmin',               &
     &              'Min temperature for pawning season (C).'
              WRITE (out,110) FspTmax(isp,ng), 'FspTmax',               &
     &              'Max temperature for pawning season (C).'
              WRITE (out,80) Hbehave(isp,ng),                           &
     &              'Hbehave',                                          &
     &              'Fish horizontal behavior'
              WRITE (out,80) Vbehave(isp,ng),                           &
     &              'Vbehave',                                          &
     &              'Fish vertical behavior'
              WRITE (out,110) ZSpref_L(isp,ng), 'ZSpref_L',             &
     &              'Larva ZSpref'
              WRITE (out,110) ZSpref_J(isp,ng), 'ZSpref_J',             &
     &              'Juven ZSpref'
              WRITE (out,110) ZSpref_A(isp,ng), 'ZSpref_A',             &
     &              'Adult ZSpref'
              WRITE (out,110) ZLpref_L(isp,ng), 'ZLpref_L',             &
     &              'Larva ZLpref'
              WRITE (out,110) ZLpref_J(isp,ng), 'ZLpref_J',             &
     &              'Juven ZLpref'
              WRITE (out,110) ZLpref_A(isp,ng), 'ZLpref_A',             &
     &              'Adult ZLpref'
              WRITE (out,110) ZPpref_L(isp,ng), 'ZPpref_L',             &
     &              'Larva ZPpref'
              WRITE (out,110) ZPpref_J(isp,ng), 'ZPpref_J',             &
     &              'Juven ZPpref'
              WRITE (out,110) ZPpref_A(isp,ng), 'ZPpref_A',             &
     &              'Adult ZPpref'
              WRITE (out,110) K_ZS_L(isp,ng), 'K_ZS_L',                 &
     &              'Larva K_ZS'
              WRITE (out,110) K_ZS_J(isp,ng), 'K_ZS_J',                 &
     &              'Juven K_ZS'
              WRITE (out,110) K_ZS_A(isp,ng), 'K_ZS_A',                 &
     &              'Adult K_ZS'
              WRITE (out,110) K_ZL_L(isp,ng), 'K_ZL_L',                 &
     &              'Larva K_ZL'
              WRITE (out,110) K_ZL_J(isp,ng), 'K_ZL_J',                 &
     &              'Juven K_ZL'
              WRITE (out,110) K_ZL_A(isp,ng), 'K_ZL_A',                 &
     &              'Adult K_ZL'
              WRITE (out,110) K_ZP_L(isp,ng), 'K_ZP_L',                 &
     &              'Larva K_ZP'
              WRITE (out,110) K_ZP_J(isp,ng), 'K_ZP_J',                 &
     &              'Juven K_ZP'
              WRITE (out,110) K_ZP_A(isp,ng), 'K_ZP_A',                 &
     &              'Adult K_ZP'
              WRITE (out,110) Cal_Z(isp,ng), 'Cal_Z',                   &
     &              'Zoop energy content'
              WRITE (out,110) Cal_F(isp,ng), 'Cal_F',                   &
     &              'Fish energy content'
              WRITE (out,110) a_C_L(isp,ng), 'a_C_L',                   &
     &              'Larva consumption'
              WRITE (out,110) a_C_J(isp,ng), 'a_C_J',                   &
     &              'Juven consumption'
              WRITE (out,110) a_C_A(isp,ng), 'a_C_A',                   &
     &              'Adult consumption'
              WRITE (out,110) b_C_L(isp,ng), 'b_C_L',                   &
     &              'Larva consumption'
              WRITE (out,110) b_C_J(isp,ng), 'b_C_J',                   &
     &              'Juven consumption'
              WRITE (out,110) b_C_A(isp,ng), 'b_C_A',                   &
     &              'Adult consumption'
              WRITE (out,110) pvalmax_L(isp,ng), 'pvalmax_L',           &
     &              'Larva Max Pvalue'
              WRITE (out,110) pvalmax_J(isp,ng), 'pvalmax_J',           &
     &              'Juven Max Pvalue'
              WRITE (out,110) pvalmax_A(isp,ng), 'pvalmax_A',           &
     &              'Adult Max Pvalue'
              WRITE (out,110) a_R_L(isp,ng), 'a_R_L',                   &
     &              'Larva respiration'
              WRITE (out,110) a_R_J(isp,ng), 'a_R_J',                   &
     &              'Juven respiration'
              WRITE (out,110) a_R_A(isp,ng), 'a_R_A',                   &
     &              'Adult respiration'
              WRITE (out,110) b_R_L(isp,ng), 'b_R_L',                   &
     &              'Larva respiration'
              WRITE (out,110) b_R_J(isp,ng), 'b_R_J',                   &
     &              'Juven respiration'
              WRITE (out,110) b_R_A(isp,ng), 'b_R_A',                   &
     &              'Adult respiration'
              WRITE (out,110) activity_L(isp,ng), 'activity_L',         &
     &              'Larva activity'
              WRITE (out,110) activity_J(isp,ng), 'activity_J',         &
     &              'Juven activity'
              WRITE (out,110) activity_A(isp,ng), 'activity_A',         &
     &              'Adult activity'
              WRITE (out,110) d_R_L(isp,ng), 'd_R_L',                   &
     &              'Larva swimming'
              WRITE (out,110) d_R_J(isp,ng), 'd_R_J',                   &
     &              'Juven swimming'
              WRITE (out,110) d_R_A(isp,ng), 'd_R_A',                   &
     &              'Adult swimming'
              WRITE (out,110) Fswim_L(isp,ng), 'Fswim_L',               &
     &              'Larva swimming'
              WRITE (out,110) Fswim_J(isp,ng), 'Fswim_J',               &
     &              'Juven swimming'
              WRITE (out,110) Fswim_A(isp,ng), 'Fswim_A',               &
     &              'Adult swimming'
              WRITE (out,110) a_AE_L(isp,ng), 'a_AE_L',                 &
     &              'Larva assimilation'
              WRITE (out,110) a_AE_J(isp,ng), 'a_AE_J',                 &
     &              'Juven assimilation'
              WRITE (out,110) a_AE_A(isp,ng), 'a_AE_A',                 &
     &              'Adult assimilation'
              WRITE (out,110) b_AE_L(isp,ng), 'b_AE_L',                 &
     &              'Larva assimilation'
              WRITE (out,110) b_AE_J(isp,ng), 'b_AE_J',                 &
     &              'Juven assimilation'
              WRITE (out,110) b_AE_A(isp,ng), 'b_AE_A',                 &
     &              'Adult assimilation'
              WRITE (out,110) AEmax_L(isp,ng), 'AEmax_L',               &
     &              'Larva assimilation'
              WRITE (out,110) AEmax_J(isp,ng), 'AEmax_J',               &
     &              'Juven assimilation'
              WRITE (out,110) AEmax_A(isp,ng), 'AEmax_A',               &
     &              'Adult assimilation'
              WRITE (out,110) te1_L(isp,ng), 'te1_L',                   &
     &              'Larva temperature dependence'
              WRITE (out,110) te1_J(isp,ng), 'te1_J',                   &
     &              'Juven temperature dependence'
              WRITE (out,110) te1_A(isp,ng), 'te1_A',                   &
     &              'Adult temperature dependence'
              WRITE (out,110) te2_L(isp,ng), 'te2_L',                   &
     &              'Larva temperature dependence'
              WRITE (out,110) te2_J(isp,ng), 'te2_J',                   &
     &              'Juven temperature dependence'
              WRITE (out,110) te2_A(isp,ng), 'te2_A',                   &
     &              'Adult temperature dependence'
              WRITE (out,110) te3_L(isp,ng), 'te3_L',                   &
     &              'Larva temperature dependence'
              WRITE (out,110) te3_J(isp,ng), 'te3_J',                   &
     &              'Juven temperature dependence'
              WRITE (out,110) te3_A(isp,ng), 'te3_A',                   &
     &              'Adult temperature dependence'
              WRITE (out,110) te4_L(isp,ng), 'te4_L',                   &
     &              'Larva temperature dependence'
              WRITE (out,110) te4_J(isp,ng), 'te4_J',                   &
     &              'Juven temperature dependence'
              WRITE (out,110) te4_A(isp,ng), 'te4_A',                   &
     &              'Adult temperature dependence'
              WRITE (out,110) xk1_L(isp,ng), 'xk1_L',                   &
     &              'Larva temperature dependence'
              WRITE (out,110) xk1_J(isp,ng), 'xk1_J',                   &
     &              'Juven temperature dependence'
              WRITE (out,110) xk1_A(isp,ng), 'xk1_A',                   &
     &              'Adult temperature dependence'
              WRITE (out,110) xk2_L(isp,ng), 'xk2_L',                   &
     &              'Larva temperature dependence'
              WRITE (out,110) xk2_J(isp,ng), 'xk2_J',                   &
     &              'Juven temperature dependence'
              WRITE (out,110) xk2_A(isp,ng), 'xk2_A',                   &
     &              'Adult temperature dependence'
              WRITE (out,110) xk3_L(isp,ng), 'xk3_L',                   &
     &              'Larva temperature dependence'
              WRITE (out,110) xk3_J(isp,ng), 'xk3_J',                   &
     &              'Juven temperature dependence'
              WRITE (out,110) xk3_A(isp,ng), 'xk3_A',                   &
     &              'Adult temperature dependence'
              WRITE (out,110) xk4_L(isp,ng), 'xk4_L',                   &
     &              'Larva temperature dependence'
              WRITE (out,110) xk4_J(isp,ng), 'xk4_J',                   &
     &              'Juven temperature dependence'
              WRITE (out,110) xk4_A(isp,ng), 'xk4_A',                   &
     &              'Adult temperature dependence'
              WRITE (out,110) cr_L(isp,ng), 'cr_L',                     &
     &              'Adult temperature dependence'
              WRITE (out,110) cr_J(isp,ng), 'cr_J',                     &
     &              'Adult temperature dependence'
              WRITE (out,110) cr_A(isp,ng), 'cr_A',                     &
     &              'Adult temperature dependence'
              WRITE (out,110) tr_L(isp,ng), 'tr_L',                     &
     &              'Adult temperature dependence'
              WRITE (out,110) tr_J(isp,ng), 'tr_J',                     &
     &              'Adult temperature dependence'
              WRITE (out,110) tr_A(isp,ng), 'tr_A',                     &
     &              'Adult temperature dependence'
              WRITE (out,110) Wffeed(isp,ng), 'Wffeed',                 &
     &              'Weight at first feeding'
              WRITE (out,110) Lffeed(isp,ng), 'Lffeed',                 &
     &              'Length at first feeding'
              WRITE (out,110) WeightLJ(isp,ng), 'WeightLJ',             &
     &              'Weight for L to J transition' 
              WRITE (out,110) LengthLJ(isp,ng), 'LengthLJ',             &
     &              'Length for L to J transition' 
              WRITE (out,110) WeightJA(isp,ng), 'WeightJA',             &
     &              'Weight for J to A transition' 
              WRITE (out,110) LengthJA(isp,ng), 'LengthJA',             &
     &              'Length for J to A transition' 
              WRITE (out,110) aw2l_L(isp,ng), 'aw2l_L',                 &
     &              'Larva weight-length conversion'
              WRITE (out,110) aw2l_J(isp,ng), 'aw2l_J',                 &
     &              'Juven weight-length conversion'
              WRITE (out,110) aw2l_A(isp,ng), 'aw2l_A',                 &
     &              'Adult weight-length conversion'
              WRITE (out,110) bw2l_L(isp,ng), 'bw2l_L',                 &
     &              'Larva weight-length conversion'
              WRITE (out,110) bw2l_J(isp,ng), 'bw2l_J',                 &
     &              'Juven weight-length conversion'
              WRITE (out,110) bw2l_A(isp,ng), 'bw2l_A',                 &
     &              'Adult weight-length conversion'
              WRITE (out,110) dSLk_L(isp,ng), 'dSLk_L',                 &
     &              'Larva weight-length conversion'
              WRITE (out,110) dSLk_J(isp,ng), 'dSLk_J',                 &
     &              'Juven weight-length conversion'
              WRITE (out,110) dSLk_A(isp,ng), 'dSLk_A',                 &
     &              'Adult weight-length conversion'
              WRITE (out,110) dSLinf_L(isp,ng), 'dSLinf_L',             &
     &              'Larva weight-length conversion'
              WRITE (out,110) dSLinf_J(isp,ng), 'dSLinf_J',             &
     &              'Juven weight-length conversion'
              WRITE (out,110) dSLinf_A(isp,ng), 'dSLinf_A',             &
     &              'Adult weight-length conversion'
              WRITE (out,110) al2w_L(isp,ng), 'al2w_L',                 &
     &              'Larva weight-length conversion'
              WRITE (out,110) al2w_J(isp,ng), 'al2w_J',                 &
     &              'Juven weight-length conversion'
              WRITE (out,110) al2w_A(isp,ng), 'al2w_A',                 &
     &              'Adult weight-length conversion'
              WRITE (out,110) bl2w_L(isp,ng), 'bl2w_L',                 &
     &              'Larva weight-length conversion'
              WRITE (out,110) bl2w_J(isp,ng), 'bl2w_J',                 &
     &              'Juven weight-length conversion'
              WRITE (out,110) bl2w_A(isp,ng), 'bl2w_A',                 &
     &              'Adult weight-length conversion'
              WRITE (out,110) abatch(isp,ng), 'abatch',                 &
     &              'Batch duration'
              WRITE (out,110) bbatch(isp,ng), 'bbatch',                 &
     &              'Batch duration'
              WRITE (out,110) T0batch(isp,ng), 'T0batch',               &
     &              'Batch duration'
              WRITE (out,110) apof(isp,ng), 'apof',                     &
     &              'Batch duration'
              WRITE (out,110) bpof(isp,ng), 'bpof',                     &
     &              'Batch duration'
              WRITE (out,110) T0pof(isp,ng), 'T0pof',                   &
     &              'Batch duration'
              WRITE (out,110) epg(isp,ng), 'epg',                       &
     &              'Batch size'
              WRITE (out,110) eegg(isp,ng), 'eegg',                     &
     &              'Batch size'
              WRITE (out,110) megg(isp,ng), 'megg',                     &
     &              'Batch size'
              WRITE (out,110) breed(isp,ng), 'breed',                   &
     &              'Breeding strategy'
              WRITE (out,110) amature(isp,ng), 'amature',               &
     &              'Maturity'
              WRITE (out,110) bmature(isp,ng), 'bmature',               &
     &              'Maturity'
              WRITE (out,110) pctxwt(isp,ng), 'pctxwt',                 &
     &              'Weight percentage'
              WRITE (out,110) pctgain(isp,ng), 'pctgain',               &
     &              'Gain percentage'
              WRITE (out,110) Nmort_E(isp,ng), 'Nmort_E',               &
     &              'Egg daily natural mortality'
              WRITE (out,110) Nmort_Y(isp,ng), 'Nmort_Y',               &
     &              'Yolksac daily natural mortality'
              WRITE (out,110) Nmort_L(isp,ng), 'Nmort_L',               &
     &              'Larva daily natural mortality'
              WRITE (out,110) Nmort_J(isp,ng), 'Nmort_J',               &
     &              'Juven daily natural mortality'
              WRITE (out,110) Nmort_A(isp,ng), 'Nmort_A',               &
     &              'Adult yearly natural mortality'
!              WRITE (out,110) Nymort(isp,ng), 'Nymort',                 &
!     &              'Natural mortality'
              WRITE (out,110) Fymort(isp,ng), 'Fymort',                 &
     &              'Fishing mortality'
              WRITE (out,110) DDmort1_J(isp,ng), 'DDmort1_J',         &
     &              'Juvenile density-dependent mortality parm 1'
              WRITE (out,110) DDmort2_J(isp,ng), 'DDmort2_J',         &
     &              'Juvenile density-dependent mortality parm 2'
              WRITE (out,110) DDmort3_J(isp,ng), 'DDmort3_J',         &
     &              'Juvenile density-dependent mortality parm 3'
              WRITE (out,110) DDmort4_J(isp,ng), 'DDmort4_J',         &
     &              'Juvenile density-dependent mortality parm 4'
              WRITE (out,110) DDscale_J(isp,ng), 'DDscale_J',           &
     &              'Juvenile density-dependent scaling factor'
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Process initial fish locations.
!-----------------------------------------------------------------------
!
!  Set time of float release (seconds after model initialization) and
!  initial float horizontal positions (grid units).  Fill the initial
!  vertical level or depth position.
!
!  For fish, we need to save some spaces at the beginning for spawned
!  eggs. Put the fish from fish.in into the end of the array so they will
!  be booted last when making new generations.
!
      DO ng=1,Ngrids
        mc=Nfish(ng)-nentry(ng)
        nc=Nfish(ng)-nentry(ng)
        DO i=1,nentry(ng)
          IF (Fcount(i,ng).eq.1) THEN
            nc=nc+1
            FISHES(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng))*day2sec
            FISHES(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
            FISHES(ng)%Ftype(nc)=Ftype(i,ng)
            IF (Fcoor(i,ng).eq.0) THEN
              FISHES(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,Fx0(i,ng)),     &
     &                                REAL(Lm(ng),r8)+0.5_r8)
              FISHES(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,Fy0(i,ng)),     &
     &                                REAL(Mm(ng),r8)+0.5_r8)
            ELSE
              mc=mc+1
              FISHES(ng)%Flon(mc)=Fx0(i,ng)
              FISHES(ng)%Flat(mc)=Fy0(i,ng)
              FISHES(ng)%Findex(mc)=nc
            END IF
            IF (Fz0(i,ng).le.0.0_r8) THEN
              FISHES(ng)%Ftype(nc)=Ftype(i,ng)
              FISHES(ng)%Fz0(nc)=Fz0(i,ng)
            END IF
          ELSE IF (Fcount(i,ng).gt.1) THEN
            DO j=1,Fcount(i,ng)
              nc=nc+1
              IF (Fdt(i,ng).gt.0.0_r8) THEN
                FISHES(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng)+           &
     &                                   REAL(j-1,r8)*Fdt(i,ng))*       &
     &                                  day2sec
                FISHES(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
                FISHES(ng)%Ftype(nc)=Ftype(i,ng)
                IF (Fcoor(i,ng).eq.0) THEN
                  FISHES(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,Fx0(i,ng)), &
     &                                    REAL(Lm(ng),r8)+0.5_r8)
                  FISHES(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,Fy0(i,ng)), &
     &                                    REAL(Mm(ng),r8)+0.5_r8)
                ELSE
                  mc=mc+1
                  FISHES(ng)%Flon(mc)=Fx0(i,ng)
                  FISHES(ng)%Flat(mc)=Fy0(i,ng)
                  FISHES(ng)%Findex(mc)=nc
                END IF
                IF (Fz0(i,ng).le.0.0_r8) THEN
                  FISHES(ng)%Ftype(nc)=Ftype(i,ng)
                  FISHES(ng)%Fz0(nc)=FISHES(ng)%Tinfo(izgrd,nc)
                END IF
              ELSE
                FISHES(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng))*day2sec
                IF (Fdz(i,ng).le.0.1_r8) THEN
                  FISHES(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
                ELSE
                  IF (Fz0(i,ng).gt.0.0_r8) THEN
                    zfloat=Fz0(i,ng)+REAL(j-1,r8)*Fdz(i,ng)
                    FISHES(ng)%Tinfo(izgrd,nc)=MIN(MAX(0.0_r8,zfloat),  &
     &                                      REAL(N(ng),r8))
                  ELSE
                    FISHES(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)+               &
     &                                      REAL(j-1,r8)*Fdz(i,ng)
                  END IF
                END IF
                IF (Fcoor(i,ng).eq.0) THEN
                  xfloat=Fx0(i,ng)+REAL(j-1,r8)*Fdx(i,ng)
                  FISHES(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,xfloat),    &
     &                                    REAL(Lm(ng),r8)+0.5_r8)
                  yfloat=Fy0(i,ng)+REAL(j-1,r8)*Fdy(i,ng)
                  FISHES(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,yfloat),    &
     &                                    REAL(Mm(ng),r8)+0.5_r8)
                ELSE
                  mc=mc+1
                  FISHES(ng)%Flon(mc)=Fx0(i,ng)+REAL(j-1,r8)*Fdx(i,ng)
                  FISHES(ng)%Flat(mc)=Fy0(i,ng)+REAL(j-1,r8)*Fdy(i,ng)
                  FISHES(ng)%Findex(mc)=nc
                END IF
              END IF
            END DO
          END IF
        END DO
        FISHES(ng)%Findex(0)=mc
      END DO
!
!  Deallocate local arrays.
!
      deallocate ( Fcoor )
      deallocate ( Fcount )
      deallocate ( Ftype )
      deallocate ( Ft0 )
      deallocate ( Fx0 )
      deallocate ( Fy0 )
      deallocate ( Fz0 )
      deallocate ( Fdt )
      deallocate ( Fdx )
      deallocate ( Fdy )
      deallocate ( Fdz )
!
  40  FORMAT (/,' READ_FishPar - Error while reading floats',           &
     &          ' locations in input script: ',a)
  50  FORMAT (/,' READ_FishPar - Error while processing line: ',/,a)
  60  FORMAT (/,' READ_FishPar - Inconsistent number of floats to',     &
     &          ' process: ', 2i6,/,18x,'change input script.')
  70  FORMAT (/,/,' Fish Initial Locations, Grid: ',i2.2,               &
     &        /,  ' ==================================',/,/,            &
     &        15x,'Ft0',5x,'Fx0',5x,'Fy0',5x,'Fz0',                     &
     &        6x,'Fdt',6x,'Fdx',6x,'Fdy',6x,'Fdz',/)
 80   FORMAT (/,1x,i10,2x,a,t30,a)
 110  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)

      RETURN
      END SUBROUTINE read_FishPar
#endif
