#include "cppdefs.h"
#ifdef TIMESERIES
      SUBROUTINE def_tms (ng, ldef)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates timeserie data NetCDF file, it defines its       !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE def_var_mod, ONLY : def_var
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng

      logical, intent(in) :: ldef
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 25

      logical :: got_var(NV)

      integer :: i, j, recdim, stadim
      integer :: status

      integer :: DimIDs(32), pgrd(2)
      integer :: Vsize(4)

      integer :: def_dim

# ifdef SOLVE3D
      integer :: itrc

      !integer :: bgrd(3), rgrd(3), wgrd(3)
      integer :: rgrd(1)
# endif

      real(r8) :: Aval(6)

      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      SourceFile='def_tms.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=TMS(ng)%name
!
      IF (Master) THEN
        IF (ldef) THEN
          WRITE (stdout,10) TRIM(ncname)
        ELSE
          WRITE (stdout,20) TRIM(ncname)
        END IF
      END IF
!
!=======================================================================
!  Create a new timeserie data file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), TMS(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          IF (Master) WRITE (stdout,30) TRIM(ncname)
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Define file dimensions.
!-----------------------------------------------------------------------
!
        DimIDs=0
!
# ifdef SOLVE3D
        status=def_dim(ng, iNLM, TMS(ng)%ncid, ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TMS(ng)%ncid, ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TMS(ng)%ncid, ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

# endif

        status=def_dim(ng, iNLM, TMS(ng)%ncid, ncname, 'tms' ,          &
     &                 1, DimIDs(13))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TMS(ng)%ncid, ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, TMS(ng)%ncid, ncname,                  &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

        recdim=DimIDs(12)
        stadim=DimIDs(13)
!
!  Define dimension vector for point variables.
!
        pgrd(1)=DimIDs(13)
        pgrd(2)=DimIDs(12)
# ifdef SOLVE3D
!
!  Define dimension vector for cast variables at vertical RHO-points.
!
!        rgrd(1)=DimIDs( 9)
!        rgrd(2)=DimIDs(13)
!        rgrd(3)=DimIDs(12)
        rgrd(1)=DimIDs(12)
!
!  Define dimension vector for cast variables at vertical W-points.
!
!        wgrd(1)=DimIDs(10)
!        wgrd(2)=DimIDs(13)
!        wgrd(3)=DimIDs(12)
!
!  Define dimension vector for sediment bed layer type variables.
!
!        bgrd(1)=DimIDs(16)
!        bgrd(2)=DimIDs(13)
!        bgrd(3)=DimIDs(12)
# endif
!
!  Initialize unlimited time record dimension.
!
        TMS(ng)%Rindex=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define time-recordless information variables.
!-----------------------------------------------------------------------
!
        CALL def_info (ng, iNLM, TMS(ng)%ncid, ncname, DimIDs)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Define variables and their attributes.
!-----------------------------------------------------------------------
!
!  Define model time.
!
        Vinfo( 1)=Vname(1,idtime)
        Vinfo( 2)=Vname(2,idtime)
        IF (INT(time_ref).eq.-2) THEN
          Vinfo( 3)='seconds since 1968-05-23 00:00:00 GMT'
          Vinfo( 4)='gregorian'
        ELSE IF (INT(time_ref).eq.-1) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='360_day'
        ELSE IF (INT(time_ref).eq.0) THEN
          Vinfo( 3)='seconds since 0001-01-01 00:00:00'
          Vinfo( 4)='julian'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='gregorian'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status=def_var(ng, iNLM, TMS(ng)%ncid, TMS(ng)%Vid(idtime),     &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN




# ifdef SOLVE3D
! RD KEEP THAT !
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
!          IF (Sout(idTvar(itrc),ng)) THEN
            Vinfo( 1)=Vname(1,idTvar(itrc))
            Vinfo( 2)=Vname(2,idTvar(itrc))
            Vinfo( 3)=Vname(3,idTvar(itrc))
            Vinfo(14)=Vname(4,idTvar(itrc))
            Vinfo(16)=Vname(1,idtime)
            status=def_var(ng, iNLM, TMS(ng)%ncid, TMS(ng)%Tid(itrc),   &
     &                     NF_FOUT, 1, rgrd, Aval, Vinfo, ncname,       &
     &                     SetFillVal = .TRUE.,                         &
     &                     SetParAccess = .FALSE.)
            IF (exit_flag.ne.NoError) RETURN
!          END IF
        END DO

# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, TMS(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, TMS(ng)%ncid, ncname)
        IF (exit_flag.ne.NoError) RETURN

      END IF DEFINE
!
!=======================================================================
!  Open an existing timeseries file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=TMS(ng)%name
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open timeseries file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, TMS(ng)%ncid)
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,50) TRIM(ncname)
          RETURN
        END IF
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  timeseries variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            TMS(ng)%Vid(idtime)=var_id(i)
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              TMS(ng)%Tid(itrc)=var_id(i)
            END IF
          END DO
# endif
        END DO
!
!  Check if timeserie variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
! RD keep that
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc)).and.Sout(idTvar(itrc),ng)) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
!
!  Set unlimited time record dimension to the appropriate value.
!
!!      TMS(ng)%Rindex=rec_size
        !TMS(ng)%Rindex=(ntstart(ng)-1)/nTMS(ng)
        TMS(ng)%Rindex=(ntstart(ng)-1)/1
      END IF QUERY
!
  10  FORMAT (6x,'DEF_TMS - creating timeseries file: ',a)
  20  FORMAT (6x,'DEF_TMS - inquiring timeseries file: ',a)
  30  FORMAT (/,' DEF_TMS - unable to create timeseries NetCDF ', &
     &        'file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_TMS - cannot open timeseries NetCDF file:', &
     &        a)
  60  FORMAT (/,' DEF_TMS - unable to find variable: ',a,2x,        &
     &        ' in timeseries NetCDF file: ',a)

      RETURN
      END SUBROUTINE def_tms
#else
      SUBROUTINE def_tms
      RETURN
      END SUBROUTINE def_tms
#endif
