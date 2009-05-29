#include "cppdefs.h"
#ifdef NONLINEAR
      SUBROUTINE output (ng)
!
!svn $Id: output.F 982 2009-05-22 20:04:01Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine manages nonlinear model output. It creates output   !
!  NetCDF files and writes out data into NetCDF files. If requested,   !
!  it can create several history and/or time-averaged files to avoid   !
!  generating too large files during a single model run.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
# ifdef FILTERED
      USE mod_filter, ONLY: nfil, FIL
# endif
# ifdef FLOATS
      USE mod_floats
# endif
# ifdef NEMURO_SAN
      USE mod_fish
# endif
# if defined FOUR_DVAR || defined VERIFICATION
      USE mod_fourdvar
# endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasts
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: Ldefine, NewFile

      integer :: ifile, lstr, status, tile
# ifdef FILTERED_RST
      integer :: i
# endif
!
      SourceFile='output.F'

# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 8)
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, process nonlinear history NetCDF file.
!-----------------------------------------------------------------------
!
!  Turn off checking for analytical header files.
!
      IF (Lanafile) THEN
        Lanafile=.FALSE.
      END IF
!
!  Create output history NetCDF file or prepare existing file to
!  append new data to it.  Also,  notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefHIS(ng)) THEN
        IF (ndefHIS(ng).gt.0) THEN
          IF (iic(ng).eq.ntstart(ng)) THEN
            NewFile=.TRUE.
            IF (idefHIS(ng).lt.0) THEN
              idefHIS(ng)=((ntstart(ng)-1)/ndefHIS(ng))*ndefHIS(ng)
              IF (idefHIS(ng).lt.iic(ng)-1) THEN
                idefHIS(ng)=idefHIS(ng)+ndefHIS(ng)
              END IF
            END IF
            IF (MOD(iic(ng)-1,ndefHIS(ng)).gt.0) THEN
              Ldefine=ldefout(ng)
            ELSE
              Ldefine=.TRUE.
            END IF
          ELSE
            NewFile=.FALSE.
          END IF
          IF ((iic(ng)-1).eq.idefHIS(ng)) THEN
            idefHIS(ng)=idefHIS(ng)+ndefHIS(ng)
            Ldefine=.TRUE.
            NewFile=.TRUE.
! think about this one...
            IF (nHIS(ng).ne.ndefHIS(ng).and.iic(ng).eq.ntstart(ng)) THEN
              idefHIS(ng)=idefHIS(ng)+nHIS(ng)
            END IF
          END IF
          IF (NewFile) THEN
            NrecHIS(ng)=0
            ifile=(iic(ng)-1)/ndefHIS(ng)+1
            IF (Master) THEN
              lstr=LEN_TRIM(HISbase(ng))
              WRITE (HISname(ng),10) HISbase(ng)(1:lstr-3),ifile
#ifdef NEP5
  10          FORMAT (a,'_',i5.5,'.nc')
#else
  10          FORMAT (a,'_',i4.4,'.nc')
#endif
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, HISname(ng))
# endif
            IF (ncHISid(ng).ne.-1) THEN
              CALL netcdf_close (ng, iNLM, ncHISid(ng))
            END IF
            CALL def_his (ng, Ldefine)
            IF (exit_flag.ne.NoError) RETURN
            LwrtHIS(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_his (ng, ldefout(ng))
            IF (exit_flag.ne.NoError) RETURN
            LwrtHIS(ng)=.TRUE.
            LdefHIS(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into history NetCDF file.  Avoid writing initial
!  conditions in perturbation mode computations.
!
      IF (LwrtHIS(ng)) THEN
        IF (LwrtPER(ng)) THEN
          IF ((iic(ng).gt.ntstart(ng)).and.                             &
     &        (MOD(iic(ng)-1,nHIS(ng)).eq.0)) THEN
            IF (nrrec(ng).eq.0.or.iic(ng).ne.ntstart(ng)) THEN
              CALL wrt_his (ng)
            END IF
            IF (exit_flag.ne.NoError) RETURN
          END IF
        ELSE
          IF (MOD(iic(ng)-1,nHIS(ng)).eq.0) THEN
            CALL wrt_his (ng)
            IF (exit_flag.ne.NoError) RETURN
          END IF
        END IF
      END IF

# if defined AVERAGES && !defined ADJOINT
!
!-----------------------------------------------------------------------
!  If appropriate, process time-averaged NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output time-averaged NetCDF file or prepare existing file
!  to append new data to it. Also, notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefAVG(ng)) THEN
        IF (ndefAVG(ng).gt.0) THEN
          IF (iic(ng).eq.ntstart(ng)) THEN
            NewFile=.TRUE.
            IF (idefAVG(ng).lt.0) THEN
              idefAVG(ng)=((ntstart(ng)-1)/ndefAVG(ng))*ndefAVG(ng)
              IF (idefAVG(ng).lt.iic(ng)-1) THEN
                idefAVG(ng)=idefAVG(ng)+ndefAVG(ng)
              END IF
            END IF
            IF (MOD(iic(ng)-1,ndefAVG(ng)).gt.0) THEN
              Ldefine=ldefout(ng)
            ELSE
              Ldefine=.TRUE.
            END IF
          ELSE
            NewFile=.FALSE.
          END IF
          IF ((iic(ng)-1).eq.idefAVG(ng)) THEN
            idefAVG(ng)=idefAVG(ng)+ndefAVG(ng)
            Ldefine=.TRUE.
            NewFile=.TRUE.
! think about this one...
            IF (nAVG(ng).ne.ndefAVG(ng).and.iic(ng).eq.ntstart(ng)) THEN
              idefAVG(ng)=idefAVG(ng)+nAVG(ng)
            END IF
          END IF
          IF (NewFile) THEN
            ifile=(iic(ng)-1)/ndefAVG(ng)+1
            IF (Master) THEN
              lstr=LEN_TRIM(AVGbase(ng))
              WRITE (AVGname(ng),20) AVGbase(ng)(1:lstr-3),ifile
#ifdef NEP5
  20          FORMAT (a,'_',i5.5,'.nc')
#else
  20          FORMAT (a,'_',i4.4,'.nc')
#endif
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, AVGname(ng))
# endif
            IF (ncAVGid(ng).ne.-1) THEN
              CALL netcdf_close (ng, iNLM, ncAVGid(ng))
            END IF
            CALL def_avg (ng, Ldefine)
            IF (exit_flag.ne.NoError) RETURN
            LwrtAVG(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_avg (ng, ldefout(ng))
            IF (exit_flag.ne.NoError) RETURN
            LwrtAVG(ng)=.TRUE.
            LdefAVG(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into time-averaged NetCDF file.
!
      IF (LwrtAVG(ng)) THEN
        IF ((iic(ng).gt.ntstart(ng)).and.                               &
     &      (MOD(iic(ng)-1,nAVG(ng)).eq.0)) THEN
          CALL wrt_avg (ng)
          IF (exit_flag.ne.NoError) RETURN
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
          CALL wrt_tides (ng)
          IF (exit_flag.ne.NoError) RETURN
# endif
        END IF
      END IF
# endif
# if defined AVERAGES2
!
!-----------------------------------------------------------------------
!  If appropriate, process secondary time-averaged NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output time-averaged NetCDF file or prepare existing file
!  to append new data to it. Also, notice that it is possible to
!  create several files during a single model run.
!
      IF (LdefAVG2(ng)) THEN
        IF (ndefAVG2(ng).gt.0) THEN
          IF (iic(ng).eq.ntstart(ng)) THEN
            NewFile=.TRUE.
            IF (idefAVG2(ng).lt.0) THEN
              idefAVG2(ng)=((ntstart(ng)-1)/ndefAVG2(ng))*ndefAVG2(ng)
              IF (idefAVG2(ng).lt.iic(ng)-1) THEN
                idefAVG2(ng)=idefAVG2(ng)+ndefAVG2(ng)
              END IF
            END IF
            IF (MOD(iic(ng)-1,ndefAVG2(ng)).gt.0) THEN
              Ldefine=ldefout(ng)
            ELSE
              Ldefine=.TRUE.
            END IF
          ELSE
            NewFile=.FALSE.
          END IF
          IF ((iic(ng)-1).eq.idefAVG2(ng)) THEN
            idefAVG2(ng)=idefAVG2(ng)+ndefAVG2(ng)
            Ldefine=.TRUE.
            NewFile=.TRUE.
! think about this one...
            IF (nAVG2(ng).ne.ndefAVG2(ng).and.iic(ng).eq.ntstart(ng))   &
     &              THEN
              idefAVG2(ng)=idefAVG2(ng)+nAVG2(ng)
            END IF
          END IF
          IF (NewFile) THEN
            ifile=(iic(ng)-1)/ndefAVG2(ng)+1
            IF (Master) THEN
              lstr=LEN_TRIM(AVG2base(ng))
              WRITE (AVG2name(ng),20) AVG2base(ng)(1:lstr-3),ifile
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, AVG2name(ng))
# endif
            IF (ncAVG2id(ng).ne.-1) THEN
              CALL netcdf_close(ng, iNLM, ncAVG2id(ng))
            END IF
            CALL def_avg2 (ng, Ldefine)
            IF (exit_flag.ne.NoError) RETURN
            LwrtAVG2(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_avg2 (ng, ldefout(ng))
            IF (exit_flag.ne.NoError) RETURN
            LwrtAVG2(ng)=.TRUE.
            LdefAVG2(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into time-averaged NetCDF file.
!
      IF (LwrtAVG2(ng)) THEN
        IF ((iic(ng).gt.ntstart(ng)).and.                               &
     &      (MOD(iic(ng)-1,nAVG2(ng)).eq.0)) THEN
          CALL wrt_avg2 (ng)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif
# ifdef DIAGNOSTICS
!
!-----------------------------------------------------------------------
!  If appropriate, process time-averaged diagnostics NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output time-averaged diagnostics NetCDF file or prepare
!  existing file to append new data to it. Also, notice that it is
!  possible to create several files during a single model run.
!
      IF (LdefDIA(ng)) THEN
        IF (ndefDIA(ng).gt.0) THEN
          IF (iic(ng).eq.ntstart(ng)) THEN
            NewFile=.TRUE.
            IF (idefDIA(ng).lt.0) THEN
              idefDIA(ng)=((ntstart(ng)-1)/ndefDIA(ng))*ndefDIA(ng)
              IF (idefDIA(ng).lt.iic(ng)-1) THEN
                idefDIA(ng)=idefDIA(ng)+ndefDIA(ng)
              END IF
            END IF
            IF (MOD(iic(ng)-1,ndefDIA(ng)).gt.0) THEN
              Ldefine=ldefout(ng)
            ELSE
              Ldefine=.TRUE.
            END IF
          ELSE
            NewFile=.FALSE.
          END IF
          IF ((iic(ng)-1).eq.idefDIA(ng)) THEN
            idefDIA(ng)=idefDIA(ng)+ndefDIA(ng)
            Ldefine=.TRUE.
            NewFile=.TRUE.
! think about this one...
            IF (nDIA(ng).ne.ndefDIA(ng).and.iic(ng).eq.ntstart(ng)) THEN
              idefDIA(ng)=idefDIA(ng)+nDIA(ng)
            END IF
          END IF
          IF (NewFile) THEN
            ifile=(iic(ng)-1)/ndefDIA(ng)+1
            IF (Master) THEN
              lstr=LEN_TRIM(DIAbase(ng))
              WRITE (DIAname(ng),30) DIAbase(ng)(1:lstr-3),ifile
#ifdef NEP5
  30          FORMAT (a,'_',i5.5,'.nc')
#else
  30          FORMAT (a,'_',i4.4,'.nc')
#endif
            END IF
# ifdef DISTRIBUTE
            CALL mp_bcasts (ng, iNLM, DIAname(ng))
# endif
            IF (ncDIAid(ng).ne.-1) THEN
              CALL netcdf_close (ng, iNLM, ncDIAid(ng))
            END IF
            CALL def_diags (ng, Ldefine)
            IF (exit_flag.ne.NoError) RETURN
            LwrtDIA(ng)=.TRUE.
          END IF
        ELSE
          IF (iic(ng).eq.ntstart(ng)) THEN
            CALL def_diags (ng, ldefout(ng))
            IF (exit_flag.ne.NoError) RETURN
            LwrtDIA(ng)=.TRUE.
            LdefDIA(ng)=.FALSE.
          END IF
        END IF
      END IF
!
!  Write out data into time-averaged diagnostics NetCDF file.
!
      IF (LwrtDIA(ng)) THEN
        IF ((iic(ng).gt.ntstart(ng)).and.                               &
     &      (MOD(iic(ng)-1,nDIA(ng)).eq.0)) THEN
          CALL wrt_diags (ng)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif
# ifdef STATIONS
!
!-----------------------------------------------------------------------
!  If appropriate, process stations NetCDF file.
!-----------------------------------------------------------------------
!
      IF (Nstation(ng).gt.0.and.(nSTA(ng).gt.0)) THEN
!
!  Create output station NetCDF file or prepare existing file to
!  append new data to it.
!
        IF (LdefSTA(ng).and.(iic(ng).eq.ntstart(ng))) THEN
          CALL def_station (ng, ldefout(ng))
          IF (exit_flag.ne.NoError) RETURN
          LdefSTA(ng)=.FALSE.
        END IF
!
!  Write out data into stations NetCDF file.
!
        IF (MOD(iic(ng)-1,nSTA(ng)).eq.0) THEN
          CALL wrt_station (ng)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif
# ifdef FLOATS
!
!-----------------------------------------------------------------------
!  If appropriate, process floats NetCDF file.
!-----------------------------------------------------------------------
!
      IF (Nfloats(ng).gt.0.and.(nFLT(ng).gt.0)) THEN
!
!  Create output floats NetCDF file or prepare existing file to
!  append new data to it.
!
        IF (LdefFLT(ng)) THEN
          IF (frrec(ng).eq.0) THEN
            NewFile=.TRUE.
          ELSE
            NewFile=.FALSE.
          END IF
          CALL def_floats (ng, NewFile)
          IF (exit_flag.ne.NoError) RETURN
          LdefFLT(ng)=.FALSE.
        END IF
!
!  Write out data into floats NetCDF file.
!
        IF ((MOD(iic(ng)-1,nFLT(ng)).eq.0).and.                         &
     &      ((frrec(ng).eq.0).or.(iic(ng).ne.ntstart(ng)))) THEN
          CALL wrt_floats (ng)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif
# ifdef NEMURO_SAN
!
!-----------------------------------------------------------------------
!  If appropriate, process fish NetCDF file.
!-----------------------------------------------------------------------
!
      IF (Nfish(ng).gt.0.and.(nFISHES(ng).gt.0)) THEN
!
!  Create output fish NetCDF file or prepare existing file to
!  append new data to it.
!
        IF (LdefFISH(ng)) THEN
          IF (fishrrec(ng).eq.0) THEN
            NewFile=.TRUE.
          ELSE
            NewFile=.FALSE.
          END IF
          CALL def_fish (ng, NewFile)
          IF (exit_flag.ne.NoError) RETURN
          LdefFISH(ng)=.FALSE.
        END IF
!
!  Write out data into fish NetCDF file.
!
        IF ((MOD(iic(ng)-1,nFISHES(ng)).eq.0).and.                      &
     &      ((fishrrec(ng).eq.0).or.(iic(ng).ne.ntstart(ng)))) THEN
          CALL wrt_fish (ng)
          IF (exit_flag.ne.NoError) RETURN
        END IF
      END IF
# endif
!
!-----------------------------------------------------------------------
!  If appropriate, process restart NetCDF file.
!-----------------------------------------------------------------------
!
!  Create output restart NetCDF file or prepare existing file to
!  append new data to it.
!
      IF (LdefRST(ng)) THEN
        CALL def_rst (ng)
        IF (exit_flag.ne.NoError) RETURN
        LwrtRST(ng)=.TRUE.
        LdefRST(ng)=.FALSE.
#ifdef FILTERED_RST
        CALL def_filt (ng)
#endif
      END IF
!
!  Write out data into restart NetCDF file.
!
      IF (LwrtRST(ng)) THEN
        IF ((iic(ng).gt.ntstart(ng)).and.                               &
     &      (MOD(iic(ng)-1,nRST(ng)).eq.0)) THEN
          CALL wrt_rst (ng)
          IF (exit_flag.ne.NoError) RETURN
# ifdef FILTERED_RST
          DO i=1,FIL
            CALL wrt_filt (ng, i)
            IF (exit_flag.ne.NoError) RETURN
          END DO
# endif
        END IF
      END IF
# if (defined FOUR_DVAR || defined VERIFICATION) && \
     !defined OBS_SENSITIVITY
#  ifdef OBSERVATIONS
!
!-----------------------------------------------------------------------
!  If appropriate, process and write model state at observation
!  locations. Compute misfit (model-observations) cost function.
!-----------------------------------------------------------------------
!
#   ifdef DISTRIBUTE
        tile=MyRank
#   else
        tile=-1
#   endif
      IF (((time(ng)-0.5_r8*dt(ng)).le.ObsTime(ng)).and.                &
     &    (ObsTime(ng).lt.(time(ng)+0.5_r8*dt(ng)))) THEN
        ProcessObs=.TRUE.
        CALL obs_read (ng, iNLM, .FALSE.)
        CALL obs_write (ng, tile, iNLM)
#   if !(defined W4DVAR || defined IOM)
        CALL obs_cost (ng, iNLM)
#   endif
      ELSE
        ProcessObs=.FALSE.
      END IF
#  endif
# endif
# ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off output data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 8)
# endif
      RETURN
      END SUBROUTINE output
#else
      SUBROUTINE output
      RETURN
      END SUBROUTINE output
#endif
