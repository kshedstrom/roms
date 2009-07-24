#include "cppdefs.h"
      SUBROUTINE ini_fish (model)
!
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
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
#ifdef DISTRIBUTE
      USE mod_strings
#endif
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti, mp_bcasts
#endif
      USE ran_state, ONLY: ran_seed
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: model
!
!  Local variable declarations.
!
      logical :: Lwrite

      integer :: Itile, Jtile, Nghost, Ntiles, tile
      integer :: Imin, Imax, Jmin, Jmax
#ifdef DISTRIBUTE
      integer :: MaxHaloLenI, MaxHaloLenJ
#endif
      integer :: inp, out, itrc, ng, npts, sequence

      real(r8) :: cff
      real(r8), parameter :: spv = 0.0_r8
!
!-----------------------------------------------------------------------
!  Read in and report input model parameters.
!-----------------------------------------------------------------------
!
!  Set input units.
!
#if defined DISTRIBUTE || defined MODEL_COUPLING
      Lwrite=Master
      inp=1
      out=stdout
#else
      Lwrite=Master
      inp=stdinp
      out=stdout
#endif
!

#ifdef NEMURO_SAN
!
!-----------------------------------------------------------------------
!  Read in fish and ship input parameters.
!-----------------------------------------------------------------------
!
      OPEN (46, FILE=TRIM(fishposnam), FORM='formatted', STATUS='old')

      CALL read_FishPar (model, 46, out, Lwrite)
!
# ifdef FLEET
      OPEN (47, FILE=TRIM(shipposnam), FORM='formatted', STATUS='old')

      CALL read_ShipsPar (model, 47, out, Lwrite)
# endif
#endif

      END SUBROUTINE ini_fish

#ifdef NEMURO_SAN
      SUBROUTINE read_FishPar (model, inp, out, Lwrite)
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
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval
      integer :: i, j, igrid, mc, nc, ng, status

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
          ELSE IF (TRIM(KeyWord).eq.'FISHRREC') THEN
            Npts=load_i(Nval, Rval, Ngrids, fishrrec)
          ELSE IF (TRIM(KeyWord).eq.'Fwwt0') THEN
            Npts=load_r(Nval, Rval, Ngrids, Fwwt0)
          ELSE IF (TRIM(KeyWord).eq.'Fage0') THEN
            Npts=load_r(Nval, Rval, Ngrids, Fage0)
          ELSE IF (TRIM(KeyWord).eq.'ZSpref') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZSpref)
          ELSE IF (TRIM(KeyWord).eq.'ZLpref') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZLpref)
          ELSE IF (TRIM(KeyWord).eq.'ZPpref') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZPpref)
          ELSE IF (TRIM(KeyWord).eq.'Nfishperyear') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nfishperyear)
          ELSE IF (TRIM(KeyWord).eq.'Nspecies') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nspecies)
          ELSE IF (TRIM(KeyWord).eq.'Nyearclass') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nyearclass)
          ELSE IF (TRIM(KeyWord).eq.'Nbox') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nbox)
          ELSE IF (TRIM(KeyWord).eq.'Nbatch') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nbatch)
          ELSE IF (TRIM(KeyWord).eq.'Npos') THEN
            Npts=load_i(Nval, Rval, Ngrids, Npos)
          ELSE IF (TRIM(KeyWord).eq.'POS') THEN
            DO ng=1,Ngrids
              Nfish(ng) = Nfishperyear(ng)*Nspecies(ng)*Nyearclass(ng)
            END DO
            Npts=Npos(1)+1
            IF (Ngrids.gt.1) Npts=MAXVAL(Npos)+1
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
            DO WHILE (.TRUE.)
              READ (inp,*,ERR=30,END=30) igrid,                         &
     &                                   Fcoor (nentry(igrid)+1,igrid), &
     &                                   Ftype (nentry(igrid)+1,igrid), &
     &                                   Fcount(nentry(igrid)+1,igrid), &
     &                                   Ft0(nentry(igrid)+1,igrid),    &
     &                                   Fx0(nentry(igrid)+1,igrid),    &
     &                                   Fy0(nentry(igrid)+1,igrid),    &
     &                                   Fz0(nentry(igrid)+1,igrid),    &
     &                                   Fdt(nentry(igrid)+1,igrid),    &
     &                                   Fdx(nentry(igrid)+1,igrid),    &
     &                                   Fdy(nentry(igrid)+1,igrid),    &
     &                                   Fdz(nentry(igrid)+1,igrid)
              IF (igrid.gt.Ngrids) THEN
                IF (Lwrite) WRITE (out,40) fposnam
                exit_flag=4
                RETURN
              END IF
              ncount(igrid)=ncount(igrid)+Fcount(nentry(igrid)+1,igrid)
              nentry(igrid)=nentry(igrid)+1        
            END DO
          END IF
        END IF
      END DO
  20  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  30  CONTINUE
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
                IF (Fdz(i,ng).eq.0.0_r8) THEN
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

! was elsewhere
!            WRITE (out,100) Fwwt0(ng), 'Fwwt0',                         &
!     &            'Fish initial wet weight (g).'
!            WRITE (out,100) ZSpref(ng), 'ZSpref',                       &
!     &            'Fish preference for ZS (-).'
!            WRITE (out,100) ZLpref(ng), 'ZLpref',                       &
!     &            'Fish preference for ZL (-).'
!            WRITE (out,100) ZPpref(ng), 'ZPpref',                       &
!     &            'Fish preference for ZP (-).'
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

      RETURN
      END SUBROUTINE read_FishPar
# ifdef FLEET
      SUBROUTINE read_ShipsPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in input ship parameters.                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ships
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval
      integer :: i, j, igrid, mc, nc, ng, status

      integer, dimension(Ngrids) :: ncount, nentry

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
!  Read in initial ship locations.
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
          IF (TRIM(KeyWord).eq.'Lships') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lships)
          ELSE IF (TRIM(KeyWord).eq.'SHIPRREC') THEN
            Npts=load_i(Nval, Rval, Ngrids, shiprrec)
          ELSE IF (TRIM(KeyWord).eq.'NSHIPS') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nships)
          ELSE IF (TRIM(KeyWord).eq.'POS') THEN
            Npts=Nships(1)+1
            IF (Ngrids.gt.1) Npts=MAXVAL(Nships)+1
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
              CALL allocate_ships (ng)
            END DO
            ncount(1:Ngrids)=0
            nentry(1:Ngrids)=0
            DO WHILE (.TRUE.)
              READ (inp,*,ERR=30,END=30) igrid,                         &
     &                                   Fcoor (nentry(igrid)+1,igrid), &
     &                                   Ftype (nentry(igrid)+1,igrid), &
     &                                   Fcount(nentry(igrid)+1,igrid), &
     &                                   Ft0(nentry(igrid)+1,igrid),    &
     &                                   Fx0(nentry(igrid)+1,igrid),    &
     &                                   Fy0(nentry(igrid)+1,igrid),    &
     &                                   Fz0(nentry(igrid)+1,igrid),    &
     &                                   Fdt(nentry(igrid)+1,igrid),    &
     &                                   Fdx(nentry(igrid)+1,igrid),    &
     &                                   Fdy(nentry(igrid)+1,igrid),    &
     &                                   Fdz(nentry(igrid)+1,igrid)
              IF (igrid.gt.Ngrids) THEN
                IF (Lwrite) WRITE (out,40) fposnam
                exit_flag=4
                RETURN
              END IF
              ncount(igrid)=ncount(igrid)+Fcount(nentry(igrid)+1,igrid)
              nentry(igrid)=nentry(igrid)+1        
            END DO
          END IF
        END IF
      END DO
  20  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  30  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (ncount(ng).ne.Nships(ng)) THEN
            IF (Master) WRITE (stdout,60) ncount(ng), Nships(ng)
            exit_flag=4
            RETURN
          END IF
          IF (Lships(ng)) THEN
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
            WRITE (out,80) Nships(ng),                                  &
     &            'Nships',                                             &
     &            'Number of ship trajectories to compute.'
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Process initial ship locations.
!-----------------------------------------------------------------------
!
!  Set time of ship release (seconds after model initialization) and
!  initial ship horizontal positions (grid units).  Fill the initial
!  vertical level or depth position. FOO - SHIPS ARE AT SURFACE!
!
      DO ng=1,Ngrids
        mc=0
        nc=0
        DO i=1,nentry(ng)
          IF (Fcount(i,ng).eq.1) THEN
            nc=nc+1
            SHIPS(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng))*day2sec
            SHIPS(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
            SHIPS(ng)%Ftype(nc)=Ftype(i,ng)
            IF (Fcoor(i,ng).eq.0) THEN
              SHIPS(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,Fx0(i,ng)),      &
     &                                REAL(Lm(ng),r8)+0.5_r8)
              SHIPS(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,Fy0(i,ng)),      &
     &                                REAL(Mm(ng),r8)+0.5_r8)
            ELSE
              mc=mc+1
              SHIPS(ng)%Flon(mc)=Fx0(i,ng)
              SHIPS(ng)%Flat(mc)=Fy0(i,ng)
              SHIPS(ng)%Findex(mc)=nc
            END IF
            IF (Fz0(i,ng).le.0.0_r8) THEN
              SHIPS(ng)%Ftype(nc)=Ftype(i,ng)
              SHIPS(ng)%Fz0(nc)=Fz0(i,ng)
            END IF
          ELSE IF (Fcount(i,ng).gt.1) THEN
            DO j=1,Fcount(i,ng)
              nc=nc+1
              IF (Fdt(i,ng).gt.0.0_r8) THEN
                SHIPS(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng)+            &
     &                                   REAL(j-1,r8)*Fdt(i,ng))*       &
     &                                  day2sec
                SHIPS(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
                SHIPS(ng)%Ftype(nc)=Ftype(i,ng)
                IF (Fcoor(i,ng).eq.0) THEN
                  SHIPS(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,Fx0(i,ng)), &
     &                                    REAL(Lm(ng),r8)+0.5_r8)
                  SHIPS(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,Fy0(i,ng)), &
     &                                    REAL(Mm(ng),r8)+0.5_r8)
                ELSE
                  mc=mc+1
                  SHIPS(ng)%Flon(mc)=Fx0(i,ng)
                  SHIPS(ng)%Flat(mc)=Fy0(i,ng)
                  SHIPS(ng)%Findex(mc)=nc
                END IF
                IF (Fz0(i,ng).le.0.0_r8) THEN
                  SHIPS(ng)%Ftype(nc)=Ftype(i,ng)
                  SHIPS(ng)%Fz0(nc)=SHIPS(ng)%Tinfo(izgrd,nc)
                END IF
              ELSE
                SHIPS(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng))*day2sec
                IF (Fdz(i,ng).eq.0.0_r8) THEN
                  SHIPS(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
                ELSE
                  IF (Fz0(i,ng).gt.0.0_r8) THEN
                    zfloat=Fz0(i,ng)+REAL(j-1,r8)*Fdz(i,ng)
                    SHIPS(ng)%Tinfo(izgrd,nc)=MIN(MAX(0.0_r8,zfloat),   &
     &                                      REAL(N(ng),r8))
                  ELSE
                    SHIPS(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)+                &
     &                                      REAL(j-1,r8)*Fdz(i,ng)
                  END IF
                END IF
                IF (Fcoor(i,ng).eq.0) THEN
                  xfloat=Fx0(i,ng)+REAL(j-1,r8)*Fdx(i,ng)
                  SHIPS(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,xfloat),     &
     &                                    REAL(Lm(ng),r8)+0.5_r8)
                  yfloat=Fy0(i,ng)+REAL(j-1,r8)*Fdy(i,ng)
                  SHIPS(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,yfloat),     &
     &                                    REAL(Mm(ng),r8)+0.5_r8)
                ELSE
                  mc=mc+1
                  SHIPS(ng)%Flon(mc)=Fx0(i,ng)+REAL(j-1,r8)*Fdx(i,ng)
                  SHIPS(ng)%Flat(mc)=Fy0(i,ng)+REAL(j-1,r8)*Fdy(i,ng)
                  SHIPS(ng)%Findex(mc)=nc
                END IF
              END IF
            END DO
          END IF
        END DO
        SHIPS(ng)%Findex(0)=mc
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
  40  FORMAT (/,' READ_ShipsPar - Error while reading floats',          &
     &          ' locations in input script: ',a)
  50  FORMAT (/,' READ_ShipsPar - Error while processing line: ',/,a)
  60  FORMAT (/,' READ_ShipsPar - Inconsistent number of floats to',    &
     &          ' process: ', 2i6,/,18x,'change input script.')
  70  FORMAT (/,/,' Ships Initial Locations, Grid: ',i2.2,              &
     &        /,  ' ==================================',/,/,            &
     &        15x,'Ft0',5x,'Fx0',5x,'Fy0',5x,'Fz0',                     &
     &        6x,'Fdt',6x,'Fdx',6x,'Fdy',6x,'Fdz',/)
 80   FORMAT (/,1x,i10,2x,a,t30,a)

      RETURN
      END SUBROUTINE read_ShipsPar
# endif
#endif
