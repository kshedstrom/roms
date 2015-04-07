#include "cppdefs.h"
#if defined NEMURO_SAN && defined FISHING_FLEET
      SUBROUTINE ini_fleet
!
!================================================== Kate Hedstrom ======
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine initializes the fishing fleet.                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_fleet
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
      integer :: i, ib, ip, icell, jcell

      real(r8) :: snudg
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
!  Read in fishing fleet input parameters.
!-----------------------------------------------------------------------
!
      OPEN (48, FILE=TRIM(fleetposnam), FORM='formatted', STATUS='old')

      CALL read_FleetPar (48, out, Lwrite)

!  Don't have these things for each grid... yet.
      ng = 1

# ifdef DISTRIBUTE
      IF (Master) THEN
        CALL ran1 (BOATS(ng)%rwalk)
      END IF
      CALL mp_bcastf (ng, iNLM, BOATS(ng)%rwalk)
# elif defined _OPENMP
!$OMP SINGLE
      CALL ran1 (BOATS(ng)%rwalk)
!$OMP END SINGLE
# else
!     IF (Lstr.eq.1) THEN
        CALL ran1 (BOATS(ng)%rwalk)
!     END IF
# endif

!
!  Put the new boats in here
      DO ib=1,Nboats(ng)
!        ip=INT(REAL(Nports(ng),r8)*BOATS(ng)%rwalk(ib))+1
        CALL ran1 (snudg)
        ip=INT(REAL(Nports(ng),r8)*snudg)+1
        BOATS(ng)%boat(ibport,ib)=ip
        BOATS(ng)%boat(ibiloc,ib)=iPort(ip,ng)
        BOATS(ng)%boat(ibjloc,ib)=jPort(ip,ng)
      END DO
      BOATS(ng)%initcpue=.TRUE.
      BOATS(ng)%initports=.TRUE.

      END SUBROUTINE ini_fleet

      SUBROUTINE read_FleetPar (inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in input fleet parameters.                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_fleet
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
      integer :: i, j, igrid, mc, nc, ng, status, ip

      integer :: decode_line, load_i, load_l, load_r

      real(r8), dimension(100) :: Rval

      character (len=35) :: frmt
      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in fishing fleet parameters.
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
          IF (TRIM(KeyWord).eq.'Lfleet') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lfleet)
          ELSE IF (TRIM(KeyWord).eq.'Nboats') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nboats)
          ELSE IF (TRIM(KeyWord).eq.'Nports') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nports)
          ELSE IF (TRIM(KeyWord).eq.'CatchMax') THEN
            Npts=load_r(Nval, Rval, Ngrids, CatchMax)
          ELSE IF (TRIM(KeyWord).eq.'TravCost') THEN
            Npts=load_r(Nval, Rval, Ngrids, TravCost)
          ELSE IF (TRIM(KeyWord).eq.'BoatVel') THEN
            Npts=load_r(Nval, Rval, Ngrids, BoatVel)
          ELSE IF (TRIM(KeyWord).eq.'Qcatch') THEN
            Npts=load_r(Nval, Rval, Ngrids, Qcatch)
          ELSE IF (TRIM(KeyWord).eq.'FishTime') THEN
            Npts=load_r(Nval, Rval, Ngrids, FishTime)
          ELSE IF (TRIM(KeyWord).eq.'EncMax') THEN
            Npts=load_i(Nval, Rval, Ngrids, EncMax)
          ELSE IF (TRIM(KeyWord).eq.'EncRate') THEN
            Npts=load_r(Nval, Rval, Ngrids, EncRate)
          ELSE IF (TRIM(KeyWord).eq.'CatchPrice') THEN
            Npts=load_r(Nval, Rval, max_ports*Ngrids, CatchPrice)
          ELSE IF (TRIM(KeyWord).eq.'iPort') THEN
            Npts=load_i(Nval, Rval, max_ports*Ngrids, iPort)
          ELSE IF (TRIM(KeyWord).eq.'jPort') THEN
            Npts=load_i(Nval, Rval, max_ports*Ngrids, jPort)
          END IF
        END IF
      END DO
  20  IF (Master) WRITE (out,50) line
      exit_flag=4
      RETURN
  30  CONTINUE
!
!-----------------------------------------------------------------------
!  Allocate fleet variables.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        CALL allocate_fleet (ng)
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lfleet(ng)) THEN
            WRITE (out,70) ng
            WRITE (out,80) Nboats(ng),                                  &
     &            'Nboats',                                             &
     &            'Number of boats in fishing fleet.'
            WRITE (out,80) Nports(ng),                                  &
     &            'Nports',                                             &
     &            'Number of ports for fishing fleet.'
            WRITE (out,*)
            WRITE (out,*) 'PARAMETERS FOR FISHING FLEET'
            WRITE (out,*) '----------------------------'
            WRITE (out,110) CatchMax(ng), 'CatchMax',                   &
     &              'Maximum daily catch (kg).'
            WRITE (out,110) TravCost(ng), 'TravCost',                   &
     &              'Cost per hour of traveling ($).'
            WRITE (out,110) BoatVel(ng), 'BoatVel',                     &
     &              'Boat motoring speed (km/h).'
            WRITE (out,110) Qcatch(ng), 'Qcatch',                       &
     &              'Catchability (-).'
            WRITE (out,110) FishTime(ng), 'FishTime',                   &
     &              'Time to fish at a location (hr).'
            WRITE (out,80) EncMax(ng),                                  &
     &              'EncMax',                                           &
     &              'Maximum numer of encounters with fish (-).'
            WRITE (out,110) EncRate(ng), 'EncRate',                     &
     &              'Mean encounter rate (-).'
            DO ip=1,Nports(ng)
              WRITE (out,110) CatchPrice(ip,ng), 'CatchPrice',          &
     &              'Price for catch by port ($/kg).'
            END DO
            DO ip=1,Nports(ng)
              WRITE (out,80) iPort(ip,ng),                              &
     &              'iPort',                                            &
     &              'Port grid i location.'
            END DO
            DO ip=1,Nports(ng)
              WRITE (out,80) jPort(ip,ng),                              &
     &              'jPort',                                            &
     &              'Port grid j location.'
            END DO
          END IF
        END DO
      END IF
!
!
  40  FORMAT (/,' READ_FleetPar - Error while reading floats',          &
     &          ' locations in input script: ',a)
  50  FORMAT (/,' READ_FleetPar - Error while processing line: ',/,a)
  60  FORMAT (/,' READ_FleetPar - Inconsistent number of floats to',    &
     &          ' process: ', 2i6,/,18x,'change input script.')
  70  FORMAT (/,/,' Fleet Initial Locations, Grid: ',i2.2,              &
     &        /,  ' ==================================',/,/,            &
     &        15x,'Ft0',5x,'Fx0',5x,'Fy0',5x,'Fz0',                     &
     &        6x,'Fdt',6x,'Fdx',6x,'Fdy',6x,'Fdz',/)
 80   FORMAT (/,1x,i10,2x,a,t30,a)
 110  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)

      RETURN
      END SUBROUTINE read_FleetPar
#endif
