#include "cppdefs.h"
#if defined NEMURO_SAN && defined FISHING_FLEET
      SUBROUTINE def_fleet (ng, ldef)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates FLEET NetCDF file, it defines dimensions,      !
!  attributes, and variables.                                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fleet
      USE mod_fish
      USE mod_biology
# ifdef FOUR_DVAR
      USE mod_fourdvar
# endif
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# ifdef SEDIMENT
      USE mod_sediment
# endif
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

      logical :: got_var(-6:NV)

      integer :: fleetdim, i, itrc, j, l
      integer :: recdim, status

      integer :: DimIDs(33), fgrd(2), t2dgrd(3), start(2), total(2)
      integer :: Vsize(4)

      integer :: def_dim, varid

      real(r8) :: Aval(6), Tinp(Nboats(ng))

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: ncname
!
      SourceFile='def_fleet.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=FLEET(ng)%name
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
!  Create a new fleet data file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), FLEET(ng)%ncid)
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
        status=def_dim(ng, iNLM, FLEET(ng)%ncid, ncname, 'xi_rho',      &
     &                 IOBOUNDS(ng)%xi_rho, DimIDs( 1))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FLEET(ng)%ncid, ncname, 'xi_u',        &
     &                 IOBOUNDS(ng)%xi_u, DimIDs( 2))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FLEET(ng)%ncid, ncname, 'xi_v',        &
     &                 IOBOUNDS(ng)%xi_v, DimIDs( 3))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FLEET(ng)%ncid, ncname, 'xi_psi',      &
     &                 IOBOUNDS(ng)%xi_psi, DimIDs( 4))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FLEET(ng)%ncid, ncname, 'eta_rho',     &
     &                 IOBOUNDS(ng)%eta_rho, DimIDs( 5))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FLEET(ng)%ncid, ncname, 'eta_u',       &
     &                 IOBOUNDS(ng)%eta_u, DimIDs( 6))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FLEET(ng)%ncid, ncname, 'eta_v',       &
     &                 IOBOUNDS(ng)%eta_v, DimIDs( 7))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FLEET(ng)%ncid, ncname, 'eta_psi',     &
     &                 IOBOUNDS(ng)%eta_psi, DimIDs( 8))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FLEET(ng)%ncid, ncname, 'tracer',      &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FLEET(ng)%ncid, ncname, 'drifter' ,    &
     &                 Nboats(ng), DimIDs(15))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FLEET(ng)%ncid, ncname, 'boundary',    &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iNLM, FLEET(ng)%ncid, ncname, 'Nstate',      &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iNLM, FLEET(ng)%ncid, ncname,                &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

        recdim=DimIDs(12)
        fleetdim=DimIDs(15)
!
!  Define dimension vectors for point variables.
!
        fgrd(1)=DimIDs(15)
        fgrd(2)=DimIDs(12)
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        FLEET(ng)%Rindex=0
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
        CALL def_info (ng, iNLM, FLEET(ng)%ncid, ncname, DimIDs)
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
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idtime),                           &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fleet X-grid locations.
!
        Vinfo( 1)='Xport'
        Vinfo( 2)='port x location'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='xport, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idBxprt),                          &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fleet Y-grid locations.
!
        Vinfo( 1)='Yport'
        Vinfo( 2)='port y location'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='yport, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idByprt),                          &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fleet port locations.
!
        Vinfo( 1)='Nport'
        Vinfo( 2)='port number'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='Nport, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idBnprt),                          &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define number of fishing locations.
!
        Vinfo( 1)='Nfloc'
        Vinfo( 2)='numer of fishing locations'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='Nfloc, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idBnfloc),                         &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fleet X-grid 1st fishing location.
!
        Vinfo( 1)='Xfish1'
        Vinfo( 2)='fishing x location 1'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='xfish1, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idBxfsh1),                         &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fleet Y-grid 1st fishing location.
!
        Vinfo( 1)='Yfish1'
        Vinfo( 2)='fishing y location 1'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='yfish1, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idByfsh1),                         &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fleet X-grid 2st fishing location.
!
        Vinfo( 1)='Xfish2'
        Vinfo( 2)='fishing x location 2'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='xfish2, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idBxfsh2),                         &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fleet Y-grid 1st fishing location.
!
        Vinfo( 1)='Yfish2'
        Vinfo( 2)='fishing y location 2'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='yfish2, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idByfsh2),                         &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fleet catch
!
        Vinfo( 1)='Catch'
        Vinfo( 2)='total daily catch'
        Vinfo( 3)='kg'
        Vinfo(14)='catch, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idBcatch),                         &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fleet cummulative catch
!
        Vinfo( 1)='CumCatch'
        Vinfo( 2)='cummulative annual catch'
        Vinfo( 3)='kg'
        Vinfo(14)='cumcatch, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idBcmctch),                        &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define average catch per unit effort
!
        Vinfo( 1)='AvgCPUE'
        Vinfo( 2)='average catch per unit effort'
        Vinfo( 3)='kg'
        Vinfo(14)='avgcpue, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idBacpue),                         &
     &                 NF_FOUT, 3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define perceived catch per unit effort
!
        Vinfo( 1)='PerCPUE'
        Vinfo( 2)='perceived catch per unit effort'
        Vinfo( 3)='kg'
        Vinfo(14)='percpue, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idBpcpue),                         &
     &                 NF_FOUT, 3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define distance to nearest port
!
        Vinfo( 1)='DistPort'
        Vinfo( 2)='distance to nearest port'
        Vinfo( 3)='km'
        Vinfo(14)='distport, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idBpdist),                         &
     &                 NF_FOUT, 3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define nearest port number
!
        Vinfo( 1)='NearPort'
        Vinfo( 2)='nearest port number'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='nearport, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FLEET(ng)%ncid,                        &
     &                 FLEET(ng)%Vid(idBpnear),                         &
     &                 NF_FOUT, 3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Initialize unlimited time record dimension.
!
        FLEET(ng)%Rindex=0
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, FLEET(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, FLEET(ng)%ncid, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing fleet file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=FLEET(ng)%name
!
!  Inquire about the dimensions and check for consistency.
!
        CALL netcdf_check_dim (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Get the size of the drifter dimension.
!
!       DO i=1,n_dim
!         IF (TRIM(dim_name(i)).eq.'drifter') THEN
!           Nboats(ng)=dim_size(i)
!           EXIT
!         END IF
!       END DO
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open fleet file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, FLEET(ng)%ncid)
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
!  fleet variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            FLEET(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Xport') THEN
            got_var(idBxprt)=.TRUE.
            FLEET(ng)%Vid(idBxprt)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Yport') THEN
            got_var(idByprt)=.TRUE.
            FLEET(ng)%Vid(idByprt)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Nport') THEN
            got_var(idBnprt)=.TRUE.
            FLEET(ng)%Vid(idBnprt)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Nfloc') THEN
            got_var(idBnfloc)=.TRUE.
            FLEET(ng)%Vid(idBnfloc)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Xfish1') THEN
            got_var(idBxfsh1)=.TRUE.
            FLEET(ng)%Vid(idBxfsh1)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Yfish1') THEN
            got_var(idByfsh1)=.TRUE.
            FLEET(ng)%Vid(idByfsh1)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Xfish2') THEN
            got_var(idBxfsh2)=.TRUE.
            FLEET(ng)%Vid(idBxfsh2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Yfish2') THEN
            got_var(idByfsh2)=.TRUE.
            FLEET(ng)%Vid(idByfsh2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Catch') THEN
            got_var(idBcatch)=.TRUE.
            FLEET(ng)%Vid(idBcatch)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'CumCatch') THEN
            got_var(idBcmctch)=.TRUE.
            FLEET(ng)%Vid(idBcmctch)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'AvgCPUE') THEN
            got_var(idBacpue)=.TRUE.
            FLEET(ng)%Vid(idBacpue)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'PerCPUE') THEN
            got_var(idBpcpue)=.TRUE.
            FLEET(ng)%Vid(idBpcpue)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'DistPort') THEN
            got_var(idBpdist)=.TRUE.
            FLEET(ng)%Vid(idBpdist)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'NearPort') THEN
            got_var(idBpnear)=.TRUE.
            FLEET(ng)%Vid(idBpnear)=var_id(i)
          END IF
        END DO
!
!  Check if fleet variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idBxprt)) THEN
          IF (Master) WRITE (stdout,60) 'Xport', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idByprt)) THEN
          IF (Master) WRITE (stdout,60) 'Yport', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idBnprt)) THEN
          IF (Master) WRITE (stdout,60) 'Nport', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Initialize fleet positions to the appropriate values.
!-----------------------------------------------------------------------
!
!  Set-up fleet time record.
!
        IF (fleetrrec(ng).lt.0) THEN
          FLEET(ng)%Rindex=rec_size
        ELSE
         FLEET(ng)%Rindex=ABS(fleetrrec(ng))
        END IF
        FLEET(ng)%Nrec=FLEET(ng)%Rindex
!
!  Read in fleet nondimentional horizontal positions.
!
        CALL netcdf_get_fvar (ng, iNLM, ncname, 'Xport',                &
     &                        Tinp,                                     &
     &                        ncid = FLEET(ng)%ncid,                    &
     &                        start = (/1,FLEET(ng)%Rindex/),           &
     &                        total = (/Nboats(ng),1/))
        IF (exit_flag.ne.NoError) RETURN

        DO l=1,Nboats(ng)
          BOATS(ng)%boat(ibiloc,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname, 'Yport',                &
     &                        Tinp,                                     &
     &                        ncid = FLEET(ng)%ncid,                    &
     &                        start = (/1,FLEET(ng)%Rindex/),           &
     &                        total = (/Nboats(ng),1/))
        IF (exit_flag.ne.NoError) RETURN

        DO l=1,Nboats(ng)
          BOATS(ng)%boat(ibjloc,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname, 'Nport',                &
     &                        Tinp,                                     &
     &                        ncid = FLEET(ng)%ncid,                    &
     &                        start = (/1,FLEET(ng)%Rindex/),           &
     &                        total = (/Nboats(ng),1/))
        IF (exit_flag.ne.NoError) RETURN

        DO l=1,Nboats(ng)
          BOATS(ng)%boat(ibport,l)=Tinp(l)
        END DO

        IF (exit_flag.ne.NoError) RETURN

      END IF QUERY
!
  10  FORMAT (6x,'DEF_FLEET  - creating fleet file: ',a)
  20  FORMAT (6x,'DEF_FLEET  - inquiring history file: ',a)
  30  FORMAT (/,' DEF_FLEET - unable to create fleet NetCDF',           &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_FLEET - unable to open fleet NetCDF file: ',a)
  60  FORMAT (/,' DEF_FLEET - unable to find variable: ',a,2x,          &
     &        ' in fleet NetCDF file: ',a)

      RETURN
      END SUBROUTINE def_fleet
#else
      SUBROUTINE def_fleet
      RETURN
      END SUBROUTINE def_fleet
#endif
