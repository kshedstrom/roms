#include "cppdefs.h"
      SUBROUTINE close_io
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
! This subroutine flushes and closes all IO files.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
#if defined BIOFLUX
      USE mod_biology
#endif
!
      implicit none
!
!  Local variable declarations.
!
      logical :: First

      integer :: MyError, i, ng, status
#if defined FILTERED_RST
      integer :: ifile
#endif

#if defined BIOFLUX
      integer :: j
      real(r8) :: dtdays
!
!---------------------------------------------------------------------
!  Write biological fluxes to standard out.
!---------------------------------------------------------------------
!
      
      write(stdout,*) ' '
      write(stdout,*) ' '
      write(stdout,134) ' Cumulative biological fluxes'
      DO ng=1,Ngrids
        IF (Master) THEN
# ifdef BIO_GOANPZ
!         Scale fluxes by maximum
          xi = 0.0_r8
          DO i = itemp,iDet
            DO j = itemp,iDet
              xi = MAX(xi, bflx(i,j))
            END DO
          END DO
          xi = 10.0_r8/xi
          write(stdout,134) ' To  :  From:'
          write(stdout,131) 'NO3','NH4','PhS','PhL','MZS','MZL','Cop',  &
     &                      'NCa','Eup','Det'
          write(stdout,132) 'NO3 : ', (xi*bflx(i,iNO3), i=iNO3,IDet)
          write(stdout,132) 'NH4 : ', (xi*bflx(i,iNH4), i=iNO3,IDet)
          write(stdout,132) 'PhS : ', (xi*bflx(i,iPhS), i=iNO3,IDet)
          write(stdout,132) 'PhL : ', (xi*bflx(i,iPhL), i=iNO3,IDet)
          write(stdout,132) 'MZS : ', (xi*bflx(i,iMZS), i=iNO3,IDet)
          write(stdout,132) 'MZL : ', (xi*bflx(i,iMZL), i=iNO3,IDet)
          write(stdout,132) 'Cop : ', (xi*bflx(i,iCop), i=iNO3,IDet)
          write(stdout,132) 'NCa : ', (xi*bflx(i,iNCa), i=iNO3,IDet)
          write(stdout,132) 'Eup : ', (xi*bflx(i,iEup), i=iNO3,IDet)
          write(stdout,132) 'Det : ', (xi*bflx(i,iDet), i=iNO3,IDet)
          write(stdout,133) 'pred: ', (xi*bflx(i,itemp), i=iCop,iEup)
          write(stdout,*) ' '
 131  format('BIOFLUX  ',7x,10a6)      
 132  format('BIOFLUX  ',a7,10f6.1)      
 133  format('BIOFLUX  ',a7,36x,3f6.1)      
 134  format('BIOFLUX  ',a)      
# endif
        END IF
      END DO
#endif
!
      SourceFile='close_io.F'
!
!-----------------------------------------------------------------------
!  Close output NetCDF files. Set file indices to closed state.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        IF (ncRSTid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncRSTid(ng))
        END IF
#if defined FORWARD_READ || defined FORWARD_WRITE
        IF ((ncFWDid(ng).ne.-1).and.(ncFWDid(ng).eq.ncHISid(ng))) THEN
          ncFWDid(ng)=-1
        END IF
        IF (ncFWDid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncFWDid(ng))
        END IF
#endif
        IF (ncHISid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncHISid(ng))
        END IF
#ifdef ADJOINT
        IF (ncADJid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iADM, ncADJid(ng))
        END IF
#endif
#ifdef TANGENT
        IF (ncTLMid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iTLM, ncTLMid(ng))
        END IF
#endif
#ifdef AVERAGES
        IF (ncAVGid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncAVGid(ng))
        END IF
#endif
#ifdef AVERAGES2
        IF (ncAVG2id(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncAVG2id(ng))
        END IF
#endif
#ifdef FILTERED_RST
        DO ifile=1,nfile
          IF (ncFILid(ifile,ng).ne.-1) THEN
            CALL netcdf_close (ng, iNLM, ncFILid(ifile,ng))
          END IF
        END DO
#endif
#ifdef DIAGNOSTICS
        IF (ncDIAid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncDIAid(ng))
        END IF
#endif
#ifdef FLOATS
        IF (ncFLTid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncFLTid(ng))
        END IF
#endif
#ifdef NEMURO_SAN
        IF (ncFISHid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncFISHid(ng))
        END IF
# ifdef FLEET
        IF (ncSHIPSid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncSHIPSid(ng))
        END IF
# endif
#endif
#ifdef STATIONS
        IF (ncSTAid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iNLM, ncSTAid(ng))
        END IF
#endif
#if defined WEAK_CONSTRAINT   && \
   (defined POSTERIOR_ERROR_F || defined POSTERIOR_ERROR_I)
        IF (ncERRid(ng).ne.-1) THEN
          CALL netcdf_close (ng, iTLM, ncERRid(ng))
        END IF
#endif
!
!  Report number of time records written.
!
        IF (Master) THEN
          WRITE (stdout,10) ng
          IF (NrecHIS(ng).gt.0) THEN
            WRITE (stdout,20) 'HISTORY', NrecHIS(ng)
          END IF
          IF (NrecRST(ng).gt.0) THEN
            IF (LcycleRST(ng)) THEN
              IF (NrecRST(ng).gt.1) THEN
                NrecRST(ng)=2
              ELSE
                NrecRST(ng)=1
              END IF
            END IF
            WRITE (stdout,20) 'RESTART', NrecRST(ng)
          END IF
#ifdef ADJOINT
          IF (NrecADJ(ng).gt.0) THEN
            WRITE (stdout,20) 'ADJOINT', NrecADJ(ng)
          END IF
#endif
#ifdef TANGENT
          IF (NrecTLM(ng).gt.0) THEN
            WRITE (stdout,20) 'TANGENT', NrecTLM(ng)
          END IF
#endif
#ifdef AVERAGES
          IF (NrecAVG(ng).gt.0) THEN
            WRITE (stdout,20) 'AVERAGE', NrecAVG(ng)
          END IF
#endif
#ifdef AVERAGES2
          IF (NrecAVG2(ng).gt.0) THEN
            WRITE (stdout,20) 'AVERAGE', NrecAVG2(ng)
          END IF
#endif
#ifdef FILTERED_RST
          IF (NrecFIL(ng).gt.0) THEN
            WRITE (stdout,20) 'FILTER RESTART' , NrecFIL(ng)
          END IF
#endif
#ifdef STATIONS
          IF (NrecSTA(ng).gt.0) THEN
            WRITE (stdout,20) 'STATION', NrecSTA(ng)
          END IF
#endif
#if defined WEAK_CONSTRAINT   && \
   (defined POSTERIOR_ERROR_F || defined POSTERIOR_ERROR_I)
          IF (NrecERR(ng).gt.0) THEN
            WRITE (stdout,20) 'ERROR  ', NrecERR(ng)
          END IF
#endif
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Report analytical header files used.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        First=.TRUE.
        DO i=1,48
          IF (LEN_TRIM(ANANAME(i)).gt.0) THEN
            IF (First) THEN
              First=.FALSE.
              WRITE (stdout,30) ' Analytical header files used:'
            END IF
            WRITE (stdout,'(5x,a)') TRIM(ADJUSTL(ANANAME(i)))
          END IF
        END DO
      END IF

#ifdef BIOLOGY
!
!-----------------------------------------------------------------------
!  Report biology model header files used.
!-----------------------------------------------------------------------
!
      IF (Master) THEN
        First=.TRUE.
        DO i=1,4
          IF (LEN_TRIM(BIONAME(i)).gt.0) THEN
            IF (First) THEN
              First=.FALSE.
              WRITE (stdout,30) ' Biology model header files used:'
            END IF
            WRITE (stdout,'(5x,a)') TRIM(ADJUSTL(BIONAME(i)))
          END IF
        END DO
      END IF
#endif
!
!-----------------------------------------------------------------------
!  If applicable, report internal exit errors.
!-----------------------------------------------------------------------
!
      IF (Master.and.(exit_flag.ne.NoError)) THEN
        WRITE (stdout,40) Rerror(exit_flag), exit_flag
      END IF
      IF (exit_flag.eq.NoError) THEN
#ifndef NO_DATE_STAMP
        CALL get_date (date_str)
#endif
        IF (Master) WRITE (stdout,50) TRIM(date_str)
      ELSE IF ((exit_flag.eq.1).or.(blowup.ne.0)) THEN
        IF (Master) WRITE (stdout,60)
      ELSE IF (exit_flag.eq.2) THEN
        IF (Master) WRITE (stdout,70) nf90_strerror(ioerror)
      ELSE IF (exit_flag.eq.3) THEN
        IF (Master) WRITE (stdout,80) nf90_strerror(ioerror)
      ELSE IF (exit_flag.eq.4) THEN
        IF (Master) WRITE (stdout,90)
      ELSE IF (exit_flag.eq.5) THEN
        IF (Master) WRITE (stdout,100)
      ELSE IF (exit_flag.eq.6) THEN
        IF (Master) WRITE (stdout,110)
      ELSE IF (exit_flag.eq.7) THEN
        IF (Master) WRITE (stdout,120)
      ELSE IF (exit_flag.eq.8) THEN
        IF (Master) WRITE (stdout,130)
      END IF
!
 10   FORMAT (/,' ROMS/TOMS - Output NetCDF summary for Grid ',         &
     &        i2.2,':')
 20   FORMAT (13x,'number of time records written in ',                 &
     &        a,' file = ',i8.8)
 30   FORMAT (/,a,/)
 40   FORMAT (/,a,i3,/)
 50   FORMAT (/,' ROMS/TOMS: DONE... ',a)
 60   FORMAT (/,' MAIN: Abnormal termination: BLOWUP.')
 70   FORMAT (/,' ERROR: Abnormal termination: NetCDF INPUT.',/,        &
     &          ' REASON: ',a)
 80   FORMAT (/,' ERROR: Abnormal termination: NetCDF OUTPUT.',/,       &
     &          ' REASON: ',a)
 90   FORMAT (/,' ERROR: I/O related problem.')
100   FORMAT (/,' ERROR: Illegal model configuration.')
110   FORMAT (/,' ERROR: Illegal domain partition.')
120   FORMAT (/,' ERROR: Illegal input parameter.')
130   FORMAT (/,' ERROR: Fatal algorithm result.')

      RETURN
      END SUBROUTINE close_io
