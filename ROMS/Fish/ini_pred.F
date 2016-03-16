#include "cppdefs.h"
#if defined NEMURO_SAN && defined PREDATOR
      SUBROUTINE ini_pred
!
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes the predators.                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_pred
      USE mod_biology
      USE mod_fish
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
      integer :: i, isp, l, pid, j

      real(r8) :: cff, Pwwght0, Pworth0, snudg
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
!  Read in predator input parameters.
!-----------------------------------------------------------------------
!
      OPEN (47, FILE=TRIM(predposnam), FORM='formatted', STATUS='old')

      CALL read_PredPar (47, out, Lwrite)

!  Don't have these things for each grid... yet.
      ng = 1

# ifdef DISTRIBUTE
      IF (Master) THEN
        CALL ran1 (PREDS(ng)%rwalk)
      END IF
      CALL mp_bcastf (ng, iNLM, PREDS(ng)%rwalk)
# elif defined _OPENMP
!$OMP SINGLE
      CALL ran1 (PREDS(ng)%rwalk)
!$OMP END SINGLE
# else
!     IF (Lstr.eq.1) THEN
        CALL ran1 (PREDS(ng)%rwalk)
!     END IF
# endif
!
!  KLUDGE WARNING!
      idpred(ip_albacore) = 1
      idpred_inv(1) = ip_albacore

!  Put the new predators in here
      pid=0
      DO isp=1,Npredspecies(ng)
        DO i=1,Npredperspecies(ng)
          pid=pid+1
          Pwwght0=Pwwt0(isp,ng)
          Pworth0=Pwth0(isp,ng)
# ifdef CCS30
! ---- FOR CCS CONFIGURATION ----
!          IF (PREDS(ng)%rwalk(pid).lt.0.5_r8) THEN
          CALL ran1 (snudg)
          IF (snudg.lt.0.5_r8) THEN
! southern domain
            yini=20.0_r8+2.0_r8*snudg*25_r8
            IF (isp.eq.ip_albacore) THEN
              xini=29.0_r8
            END IF
          ELSE
! northern domain
!           yini=55.0_r8+2.0_r8*(PREDS(ng)%rwalk(pid)-0.5_r8)*25_r8
            yini=55.0_r8+2.0_r8*(snudg-0.5_r8)*25_r8
            IF (isp.eq.ip_albacore) THEN
              xini=25.0_r8
            END IF
          END IF
!          xini=xini+2.0_r8*(PREDS(ng)%rwalk(pid+Npred(ng))-0.5_r8)
          CALL ran1 (snudg)
          xini=xini+2.0_r8*(snudg-0.5_r8)
          zini=REAL(N(ng),r8)
# endif
# ifdef WC13
! ---- FOR WC CONFIGURATION ----
!          yini=MAX(5.0_r8,PREDS(ng)%rwalk(pid)*REAL(Mm(ng)-4,r8))
          CALL ran1 (snudg)
          yini=10.0_r8+snudg*REAL(Mm(ng)-20,r8)
          IF (yini.lt.33.0_r8) THEN
            IF (isp.eq.ip_albacore) THEN
! old              xini=44.0_r8-(yini-3.0_r8)*25.0_r8/30.0_r8
              xini=36.0_r8-(yini-10.0_r8)*15.0_r8/19.0_r8
            END IF
          ELSE
            IF (isp.eq.ip_albacore) THEN
! old              xini=18.0_r8
              xini=21.0_r8
            END IF
          END IF
!          xini=xini+2.0_r8*(PREDS(ng)%rwalk(pid+Npred(ng))-0.5_r8)
          CALL ran1 (snudg)
          xini=xini+2.0_r8*(snudg-0.5_r8)
          zini=REAL(N(ng),r8)
# endif
# ifdef WC12
! ---- FOR WC CONFIGURATION ----
!          yini=MAX(5.0_r8,PREDS(ng)%rwalk(pid)*REAL(Mm(ng)-4,r8))
          CALL ran1 (snudg)
          yini=30.0_r8+snudg*REAL(Mm(ng)-60,r8)
          IF (yini.lt.99.0_r8) THEN
            IF (isp.eq.ip_albacore) THEN
! old              xini=44.0_r8-(yini-3.0_r8)*25.0_r8/30.0_r8
              xini=108.0_r8-(yini-30.0_r8)*15.0_r8/19.0_r8
            END IF
          ELSE
            IF (isp.eq.ip_albacore) THEN
! old              xini=18.0_r8
              xini=53.0_r8
            END IF
          END IF
!          xini=xini+2.0_r8*(PREDS(ng)%rwalk(pid+Npred(ng))-0.5_r8)
          CALL ran1 (snudg)
          xini=xini+6.0_r8*(snudg-0.5_r8)
          zini=REAL(N(ng),r8)
# endif
# ifdef CCS1
! ---- FOR CCS1 CONFIGURATION : may need fixes !!! RD ----
          CALL ran1 (snudg)
          yini=219.0_r8+snudg*REAL(198,r8)
          IF (yini.lt.349.0_r8) THEN
            IF (isp.eq.ip_albacore) THEN
              xini=76.0_r8
            END IF
          ELSE
            IF (isp.eq.ip_albacore) THEN
              xini=76.0_r8-(yini-349.0_r8)*0.58_r8
            END IF
          END IF
          CALL ran1 (snudg)
          xini=xini+6.0_r8*(snudg-0.5_r8)
          zini=REAL(N(ng),r8)
# endif
# ifdef NWPACIFIC
! Shinichi custom values for predator
! lowest_pt_j=90.
! lowest_pt_i=175.
! extension_j=45.
! extension_i=30.
!          CALL ran1 (snudg)
!          yini=90.0_r8+snudg*REAL(45,r8)
!          xini=175.0_r8
!          CALL ran1 (snudg)
!          xini=xini+30.0_r8*snudg
!          zini=REAL(N(ng),r8)

! Takaaki custom values for predator
! lowest_pt_j=99.
! lowest_pt_i=181.
! extension_j=45.
! extension_i=30.
          CALL ran1 (snudg)
          yini=99.0_r8+snudg*REAL(45,r8)
          xini=181.0_r8
          CALL ran1 (snudg)
          xini=xini+30.0_r8*snudg
          zini=REAL(N(ng),r8)


# endif
# ifdef NEMSAN_TEST
! ---- FOR TEST CONFIGURATION ----
!          xini=REAL(12,r8)+10.0_r8*(PREDS(ng)%rwalk(pid)-0.5_r8)
!          yini=REAL(12,r8)+10.0_r8*(PREDS(ng)%rwalk(pid)-0.5_r8)
          CALL ran1 (snudg)
!          xini=4.0_r8+(snudg-0.5_r8)
!          CALL ran1 (snudg)
!          yini=7.0_r8+(snudg-0.5_r8)
          xini=REAL(12,r8)+10.0_r8*(snudg-0.5_r8)
          yini=REAL(12,r8)+10.0_r8*(snudg-0.5_r8)
          zini=REAL(N(ng),r8)
# endif
          PREDS(ng)%Tinfo(ixgrd,pid)=xini
          PREDS(ng)%Tinfo(iygrd,pid)=yini
          PREDS(ng)%Tinfo(izgrd,pid)=zini
          DO j=0,NFT
            PREDS(ng)%track(ixgrd,j,pid)=xini
            PREDS(ng)%track(iygrd,j,pid)=yini
            PREDS(ng)%track(ixrhs,j,pid)=0.0_r8
            PREDS(ng)%track(iyrhs,j,pid)=0.0_r8
#  ifdef SOLVE3D
            PREDS(ng)%track(izgrd,j,pid)=zini
            PREDS(ng)%track(izrhs,j,pid)=0.0_r8
#  endif
          END DO
!          PREDS(ng)%bioenergy(ipwwt,pid)=Pwwght0* &
!     &              (PREDS(ng)%rwalk(pid+2*Npred(ng))+0.5_r8)
          CALL ran1 (snudg)
          cff=1.0_r8+0.2_r8*(snudg-0.5_r8)
          PREDS(ng)%bioenergy(ipwwt,pid)=cff*Pwwght0
          PREDS(ng)%bioenergy(ipworth,pid)=Pworth0
          PREDS(ng)%species(pid)=idpred_inv(isp)
          PREDS(ng)%alive(pid)=.TRUE.
          PREDS(ng)%bounded(pid)=.TRUE.
        END DO
      END DO

      END SUBROUTINE ini_pred

      SUBROUTINE read_PredPar (inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in input predators parameters.                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_fish
      USE mod_pred
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
      integer :: i, j, igrid, mc, nc, ng, status, isp, ifsp

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
!  Read in initial predator locations.
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
          IF (TRIM(KeyWord).eq.'Lpred') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lpred)
          ELSE IF (TRIM(KeyWord).eq.'Npredspecies') THEN
            Npts=load_i(Nval, Rval, Ngrids, Npredspecies)
          ELSE IF (TRIM(KeyWord).eq.'Npredperspecies') THEN
            Npts=load_i(Nval, Rval, Ngrids, Npredperspecies)
          ELSE IF (TRIM(KeyWord).eq.'Pwwt0') THEN
            Npts=load_r(Nval, Rval, max_predspecies*Ngrids, Pwwt0)
          ELSE IF (TRIM(KeyWord).eq.'Pwth0') THEN
            Npts=load_r(Nval, Rval, max_predspecies*Ngrids, Pwth0)
          ELSE IF (TRIM(KeyWord).eq.'Pmgstr') THEN
            Npts=load_r(Nval, Rval, max_predspecies*Ngrids, Pmgstr)
          ELSE IF (TRIM(KeyWord).eq.'Pmgend') THEN
            Npts=load_r(Nval, Rval, max_predspecies*Ngrids, Pmgend)
          ELSE IF (TRIM(KeyWord).eq.'Fpref') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, Fpref)
          ELSE IF (TRIM(KeyWord).eq.'K_Fish') THEN
            Npts=load_r(Nval, Rval, max_species*Ngrids, K_Fish)
          ELSE IF (TRIM(KeyWord).eq.'Pcmax') THEN
            Npts=load_r(Nval, Rval, max_predspecies*Ngrids, Pcmax)
          ELSE IF (TRIM(KeyWord).eq.'Pswim') THEN
            Npts=load_r(Nval, Rval, max_predspecies*Ngrids, Pswim)
          END IF
        END IF
      END DO
  20  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  30  CONTINUE
!
!-----------------------------------------------------------------------
!  Set initial predator locations
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        Npred(ng)=Npredperspecies(ng)*Npredspecies(ng)
        Npos(ng)=Npred(ng)
      END DO
      Npts=Npred(1)
      IF (Ngrids.gt.1) Npts=MAXVAL(Npred)
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
        CALL allocate_pred (ng)
      END DO
      ncount(1:Ngrids)=0
      nentry(1:Ngrids)=0
      DO i=1,Npts
        igrid=1
        Fcoor (nentry(igrid)+1,igrid)=0
        Ftype (nentry(igrid)+1,igrid)=1         ! 3-D Largrangian pred
        Fcount(nentry(igrid)+1,igrid)=1
        Ft0(nentry(igrid)+1,igrid)=0.0_r8
        Fx0(nentry(igrid)+1,igrid)=0.0_r8
        Fy0(nentry(igrid)+1,igrid)=0.0_r8
        Fz0(nentry(igrid)+1,igrid)=30.0_r8
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
          IF (Lpred(ng)) THEN
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
            WRITE (out,80) Npred(ng),                                   &
     &            'Npred',                                              &
     &            'Number of predator trajectories to compute.'
            WRITE (out,80) Npredperspecies(ng),                         &
     &            'Npredperspecies',                                    &
     &            'Number of predators per species.'
            WRITE (out,80) Npredspecies(ng),                            &
     &            'Npredspecies',                                       &
     &            'Number of predator species.'
            WRITE (out,80) Npos(ng),                                    &
     &            'Npos',                                               &
     &            'Number of predator initial positions.'
            DO isp=1,Npredspecies(ng)
              WRITE (out,*)
              WRITE (out,*) 'PARAMETERS FOR PREDATOR SPECIES ',isp
              WRITE (out,*) '--------------------------------------'
              WRITE (out,110) Pwwt0(isp,ng), 'Pwwt0',                   &
     &                'Pred initial wet weight (g).'
              WRITE (out,110) Pwth0(isp,ng), 'Pwth0',                   &
     &                'Pred initial worth.'
              WRITE (out,110) Pmgstr(isp,ng), 'Pmgstr',                 &
     &                'Pred migration start.'
              WRITE (out,110) Pmgend(isp,ng), 'Pmgend',                 &
     &                'Pred migration end.'
              DO ifsp=1,Nspecies(ng)
                WRITE (out,110) Fpref(ifsp,ng), 'Fpref',                &
     &                'Pred preference for Fish.'
              END DO
              DO ifsp=1,Nspecies(ng)
                WRITE (out,110) K_Fish(ifsp,ng), 'K_Fish',              &
     &                'Pred half-saturation constant for Fish.'
              END DO
              WRITE (out,110) Pcmax(isp,ng), 'Pcmax',                   &
     &              'Pred maxmimum comsumption rate.'
              WRITE (out,110) Pswim(isp,ng), 'Pswim',                   &
     &              'Pred swimming speed.'
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Process initial predator locations.
!-----------------------------------------------------------------------
!
!  Set time of float release (seconds after model initialization) and
!  initial float horizontal positions (grid units).  Fill the initial
!  vertical level or depth position.
!
!  For pred, we need to save some spaces at the beginning for spawned
!  eggs. Put the pred from pred.in into the end of the array so they will
!  be booted last when making new generations.
!
      DO ng=1,Ngrids
        mc=Npred(ng)-nentry(ng)
        nc=Npred(ng)-nentry(ng)
        DO i=1,nentry(ng)
          IF (Fcount(i,ng).eq.1) THEN
            nc=nc+1
            PREDS(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng))*day2sec
            PREDS(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
            PREDS(ng)%Ftype(nc)=Ftype(i,ng)
            IF (Fcoor(i,ng).eq.0) THEN
              PREDS(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,Fx0(i,ng)),      &
     &                                REAL(Lm(ng),r8)+0.5_r8)
              PREDS(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,Fy0(i,ng)),      &
     &                                REAL(Mm(ng),r8)+0.5_r8)
            ELSE
              mc=mc+1
              PREDS(ng)%Flon(mc)=Fx0(i,ng)
              PREDS(ng)%Flat(mc)=Fy0(i,ng)
              PREDS(ng)%Findex(mc)=nc
            END IF
            IF (Fz0(i,ng).le.0.0_r8) THEN
              PREDS(ng)%Ftype(nc)=Ftype(i,ng)
              PREDS(ng)%Fz0(nc)=Fz0(i,ng)
            END IF
          ELSE IF (Fcount(i,ng).gt.1) THEN
            DO j=1,Fcount(i,ng)
              nc=nc+1
              IF (Fdt(i,ng).gt.0.0_r8) THEN
                PREDS(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng)+            &
     &                                   REAL(j-1,r8)*Fdt(i,ng))*       &
     &                                  day2sec
                PREDS(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
                PREDS(ng)%Ftype(nc)=Ftype(i,ng)
                IF (Fcoor(i,ng).eq.0) THEN
                  PREDS(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,Fx0(i,ng)),  &
     &                                    REAL(Lm(ng),r8)+0.5_r8)
                  PREDS(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,Fy0(i,ng)),  &
     &                                    REAL(Mm(ng),r8)+0.5_r8)
                ELSE
                  mc=mc+1
                  PREDS(ng)%Flon(mc)=Fx0(i,ng)
                  PREDS(ng)%Flat(mc)=Fy0(i,ng)
                  PREDS(ng)%Findex(mc)=nc
                END IF
                IF (Fz0(i,ng).le.0.0_r8) THEN
                  PREDS(ng)%Ftype(nc)=Ftype(i,ng)
                  PREDS(ng)%Fz0(nc)=PREDS(ng)%Tinfo(izgrd,nc)
                END IF
              ELSE
                PREDS(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng))*day2sec
                IF (Fdz(i,ng).le.0.1_r8) THEN
                  PREDS(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
                ELSE
                  IF (Fz0(i,ng).gt.0.0_r8) THEN
                    zfloat=Fz0(i,ng)+REAL(j-1,r8)*Fdz(i,ng)
                    PREDS(ng)%Tinfo(izgrd,nc)=MIN(MAX(0.0_r8,zfloat),   &
     &                                      REAL(N(ng),r8))
                  ELSE
                    PREDS(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)+                &
     &                                      REAL(j-1,r8)*Fdz(i,ng)
                  END IF
                END IF
                IF (Fcoor(i,ng).eq.0) THEN
                  xfloat=Fx0(i,ng)+REAL(j-1,r8)*Fdx(i,ng)
                  PREDS(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,xfloat),     &
     &                                    REAL(Lm(ng),r8)+0.5_r8)
                  yfloat=Fy0(i,ng)+REAL(j-1,r8)*Fdy(i,ng)
                  PREDS(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,yfloat),     &
     &                                    REAL(Mm(ng),r8)+0.5_r8)
                ELSE
                  mc=mc+1
                  PREDS(ng)%Flon(mc)=Fx0(i,ng)+REAL(j-1,r8)*Fdx(i,ng)
                  PREDS(ng)%Flat(mc)=Fy0(i,ng)+REAL(j-1,r8)*Fdy(i,ng)
                  PREDS(ng)%Findex(mc)=nc
                END IF
              END IF
            END DO
          END IF
        END DO
        PREDS(ng)%Findex(0)=mc
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
  40  FORMAT (/,' READ_PredPar - Error while reading floats',           &
     &          ' locations in input script: ',a)
  50  FORMAT (/,' READ_PredPar - Error while processing line: ',/,a)
  60  FORMAT (/,' READ_PredPar - Inconsistent number of floats to',     &
     &          ' process: ', 2i6,/,18x,'change input script.')
  70  FORMAT (/,/,' Pred Initial Locations, Grid: ',i2.2,               &
     &        /,  ' ==================================',/,/,            &
     &        15x,'Ft0',5x,'Fx0',5x,'Fy0',5x,'Fz0',                     &
     &        6x,'Fdt',6x,'Fdx',6x,'Fdy',6x,'Fdz',/)
 80   FORMAT (/,1x,i10,2x,a,t30,a)
 110  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)

      RETURN
      END SUBROUTINE read_PredPar
#endif
