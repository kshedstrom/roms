#include "cppdefs.h"
#ifdef NEMURO_SAN
      SUBROUTINE def_ships (ng, ldef)
!
!svn $Id: def_ships.F 975 2009-05-05 22:51:13Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates SHIPS NetCDF file, it defines dimensions,      !
!  attributes, and variables.                                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ships
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
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_collect, mp_bcasti
# endif
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
      integer, parameter :: Natt = 24

      logical :: got_var(-6:NV)

      integer :: shipsdim, i, itrc, j, l
      integer :: recdim, status

      integer :: DimIDs(31), fgrd(2), start(2), total(2)
      integer :: Vsize(4)

      integer :: def_dim, varid

# ifdef DISTRIBUTE
      integer :: Npts

      real(r8), parameter :: Fspv = 0.0_r8
      real(r8), dimension(Nships(ng)*NFV(ng)) :: Fwrk
# endif

      real(r8) :: Aval(6), Tinp(Nships(ng))

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: ncname
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=SHIPname(ng)
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
!  Create a new ships data file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), ncSHIPSid(ng))
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
        status=def_dim(ng, iNLM, ncSHIPSid(ng), ncname, 's_rho',        &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncSHIPSid(ng), ncname, 's_w',          &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncSHIPSid(ng), ncname, 'tracer',       &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iNLM, ncSHIPSid(ng), ncname, 'Nbed',         &
     &                 Nbed, DimIDs(16))
        IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iNLM, ncSHIPSid(ng), ncname, 'Nphy',         &
     &                 Nphy, DimIDs(25))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncSHIPSid(ng), ncname, 'Nbac',         &
     &                 Nbac, DimIDs(26))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncSHIPSid(ng), ncname, 'Ndom',         &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncSHIPSid(ng), ncname, 'Nfec',         &
     &                 Nfec, DimIDs(28))
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

        status=def_dim(ng, iNLM, ncSHIPSid(ng), ncname, 'drifter' ,     &
     &                 Nships(ng), DimIDs(15))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncSHIPSid(ng), ncname, 'boundary',     &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iNLM, ncSHIPSid(ng), ncname, 'Nstate',       &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iNLM, ncSHIPSid(ng), ncname,                 &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

        recdim=DimIDs(12)
        shipsdim=DimIDs(15)
!
!  Define dimension vectors for point variables.
!
        fgrd(1)=DimIDs(15)
        fgrd(2)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        tSHIPSindx(ng)=0
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
        CALL def_info (ng, iNLM, ncSHIPSid(ng), ncname, DimIDs)
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
        status=def_var(ng, iNLM, ncSHIPSid(ng), shipsVid(idtime,ng),    &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ship types
!
        Vinfo( 1)='Ftype'
        Vinfo( 2)='Ship type'
        Vinfo( 9)='one day'
        Vinfo(10)='ten day'
        status=def_var(ng, iNLM, ncSHIPSid(ng), varid, nf90_int,        &
     &                 1, (/shipsdim/), Aval, Vinfo, ncname)
!
!  Define ships X-grid locations.
!
        Vinfo( 1)='Xgrid'
        Vinfo( 2)='x-grid ships locations'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(Lm(ng)+1,r8)
        Vinfo(14)='Xgrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, ncSHIPSid(ng), shipsVid(idXgrd,ng),    &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ship Y-grid locations.
!
        Vinfo( 1)='Ygrid'
        Vinfo( 2)='Y-grid ship locations'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(Mm(ng)+1,r8)
        Vinfo(14)='Ygrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, ncSHIPSid(ng), shipsVid(idYgrd,ng),    &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Define ship Z-grid locations.
!
        Vinfo( 1)='Zgrid'
        Vinfo( 2)='Z-grid ship locations'
        Vinfo( 5)='valid_min'
        Vinfo( 6)='valid_max'
        Aval(2)=0.0_r8
        Aval(3)=REAL(N(ng),r8)
        Vinfo(14)='Zgrid, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, ncSHIPSid(ng), shipsVid(idZgrd,ng),    &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Define ship (lon,lat) or (x,y) locations.
!
        IF (spherical) THEN
          Vinfo( 1)='lon'
          Vinfo( 2)='longitude of ship trajectories'
          Vinfo( 3)='degree_east'
          Vinfo( 5)='valid_min'
          Vinfo( 6)='valid_max'
          Vinfo(14)='lon, scalar, series'
          Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
# endif
          Aval(2)=-180.0_r8
          Aval(3)=180.0_r8
          status=def_var(ng, iNLM, ncSHIPSid(ng), shipsVid(idglon,ng),  &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN

          Vinfo( 1)='lat'
          Vinfo( 2)='latitude of ship trajectories'
          Vinfo( 3)='degree_north'
          Vinfo( 5)='valid_min'
          Vinfo( 6)='valid_max'
          Vinfo(14)='lat, scalar, series'
          Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
# endif
          Aval(2)=-90.0_r8
          Aval(3)=90.0_r8
          status=def_var(ng, iNLM, ncSHIPSid(ng), shipsVid(idglat,ng),  &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
       ELSE
          Vinfo( 1)='x'
          Vinfo( 2)='x-location of ship trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='x, scalar, series'
          Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
# endif
          status=def_var(ng, iNLM, ncSHIPSid(ng), shipsVid(idglon,ng),  &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN

          Vinfo( 1)='y'
          Vinfo( 2)='y-location of ship trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='y, scalar, series'
          Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
# endif
          status=def_var(ng, iNLM, ncSHIPSid(ng), shipsVid(idglat,ng),  &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

# ifdef SOLVE3D
!
!  Define ship depths. (FOO - doesn't make sense)
!
        Vinfo( 1)='depth'
        Vinfo( 2)='depth of ship trajectories'
        Vinfo( 3)='meter'
        Vinfo(14)='depth, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, ncSHIPSid(ng), shipsVid(iddpth,ng),    &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define density anomaly.
!
        Vinfo( 1)=Vname(1,idDano)
        Vinfo( 2)=Vname(2,idDano)
        Vinfo( 3)=Vname(3,idDano)
        Vinfo(14)=Vname(4,idDano)
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, ncSHIPSid(ng), shipsVid(idDano,ng),    &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define tracer type variables.
!
        DO itrc=1,NT(ng)
          Vinfo( 1)=Vname(1,idTvar(itrc))
          Vinfo( 2)=Vname(2,idTvar(itrc))
          Vinfo( 3)=Vname(3,idTvar(itrc))
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
#  endif
#  ifdef SEDIMENT
          DO i=1,NST
            IF (itrc.eq.idsed(i)) THEN
              WRITE (Vinfo(19),40) 1000.0_r8*Sd50(i,ng)
            END IF
          END DO
#  endif
          status=def_var(ng, iNLM, ncSHIPSid(ng), shipsTid(itrc,ng),    &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
#  ifdef NEMURO_SAN
!  Define ships wet weight
        Vinfo( 1)='fwwt'
        Vinfo( 2)='ships wet weight'
        Vinfo( 3)='g'
        Vinfo(14)='fwwt, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, ncSHIPSid(ng),shipsVid(ifwwt,ng),      &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
!
!  Initialize unlimited time record dimension.
!
        tSHIPSindx(ng)=0
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        IF (OutThread) THEN
          status=nf90_enddef(ncSHIPSid(ng))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,50) TRIM(ncname)
            exit_flag=3
            ioerror=status
          END IF
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, ncSHIPSid(ng), ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
      IF (LdefSHIPS(ng).and.OutThread) THEN
        status=nf90_inq_varid(ncSHIPSid(ng),'Ftype',varid)
        status=nf90_put_var(ncSHIPSid(ng), varid, SHIPS(ng)%Ftype)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) 'Ftype', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
      END IF
!
!=======================================================================
!  Open an existing ships file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=SHIPname(ng)
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
!           Nships(ng)=dim_size(i)
!           EXIT
!         END IF
!       END DO
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open ships file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, ncSHIPSid(ng))
        IF (exit_flag.ne.NoError) THEN
          WRITE (stdout,60) TRIM(ncname)
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
!  ship variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            shipsVid(idtime,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Xgrid') THEN
            got_var(idXgrd)=.TRUE.
            shipsVid(idXgrd,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Ygrid') THEN
            got_var(idYgrd)=.TRUE.
            shipsVid(idYgrd,ng)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.'Zgrid') THEN
            got_var(idZgrd)=.TRUE.
            shipsVid(idZgrd,ng)=var_id(i)
# endif
          ELSE IF (spherical.and.TRIM(var_name(i)).eq.'lon') THEN
            got_var(idglon)=.TRUE.
            shipsVid(idglon,ng)=var_id(i)
          ELSE IF (spherical.and.TRIM(var_name(i)).eq.'lat') THEN
            got_var(idglat)=.TRUE.
            shipsVid(idglat,ng)=var_id(i)
          ELSE IF (.not.spherical.and.TRIM(var_name(i)).eq.'x') THEN
            got_var(idglon)=.TRUE.
            shipsVid(idglon,ng)=var_id(i)
          ELSE IF (.not.spherical.and.TRIM(var_name(i)).eq.'y') THEN
            got_var(idglat)=.TRUE.
            shipsVid(idglat,ng)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.'depth') THEN
            got_var(iddpth)=.TRUE.
            shipsVid(iddpth,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            shipsVid(idDano,ng)=var_id(i)
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              shipsTid(itrc,ng)=var_id(i)
            END IF
          END DO
# endif
        END DO
!
!  Check if ships variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idXgrd)) THEN
          IF (Master) WRITE (stdout,70) 'Xgrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idYgrd)) THEN
          IF (Master) WRITE (stdout,70) 'Ygrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idZgrd)) THEN
          IF (Master) WRITE (stdout,70) 'Zgrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idglon)) THEN
          IF (spherical) THEN
            IF (Master) WRITE (stdout,70) 'lon', TRIM(ncname)
          ELSE
            IF (Master) WRITE (stdout,70) 'x', TRIM(ncname)
          END IF
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idglat)) THEN
          IF (spherical) THEN
            IF (Master) WRITE (stdout,70) 'lat', TRIM(ncname)
          ELSE
            IF (Master) WRITE (stdout,70) 'y', TRIM(ncname)
          END IF
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(iddpth)) THEN
          IF (Master) WRITE (stdout,70) 'depth', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano)) THEN
          IF (Master) WRITE (stdout,70) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
# ifdef SOLVE3D
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,70) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
!
!-----------------------------------------------------------------------
!  Initialize ships positions to the appropriate values.
!-----------------------------------------------------------------------
!
!  Set-up ship time record.
!
        IF (shipsrrec(ng).lt.0) THEN
          tSHIPSindx(ng)=rec_size
        ELSE
         tSHIPSindx(ng)=ABS(shipsrrec(ng))
        END IF
        NrecSHIPS(ng)=tSHIPSindx(ng)
!
!  Read in ship nondimentional horizontal positions.
!
        IF (OutThread) THEN
          start(1)=1
          total(1)=Nships(ng)
          start(2)=tSHIPSindx(ng)
          total(2)=1
          IF (exit_flag.eq.NoError) THEN
            status=nf90_get_var(ncSHIPSid(ng), shipsVid(idXgrd,ng),     &
     &                          Tinp, start, total)
            IF (status.eq.nf90_noerr) THEN
              DO l=1,Nships(ng)
                IF ((Tinp(l).gt.REAL(Lm(ng)+1,r8)-0.5_r8).or.           &
     &              (Tinp(l).lt.0.5_r8)) THEN
                  SHIPS(ng)%bounded(l)=.FALSE.
                ELSE
                  SHIPS(ng)%bounded(l)=.TRUE.
                  DO i=0,NFT
                    SHIPS(ng)%track(ixgrd,i,l)=Tinp(l)
                    SHIPS(ng)%track(ixrhs,i,l)=0.0_r8
                  END DO
                END IF
              END DO
            ELSE
              IF (Master) WRITE (stdout,80) 'Xgrid', tSHIPSindx(ng),    &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

          IF (exit_flag.eq.NoError) THEN
            status=nf90_get_var(ncSHIPSid(ng), shipsVid(idYgrd,ng),     &
     &                          Tinp, start, total)
            IF (status.eq.nf90_noerr) THEN
              DO l=1,Nships(ng)
                IF ((Tinp(l).gt.REAL(Mm(ng)+1,r8)-0.5_r8).or.           &
     &              (Tinp(l).lt.0.5_r8)) THEN
                  SHIPS(ng)%bounded(l)=.FALSE.
                ELSE
                  SHIPS(ng)%bounded(l)=.TRUE.
                  DO i=0,NFT
                    SHIPS(ng)%track(iygrd,i,l)=Tinp(l)
                    SHIPS(ng)%track(iyrhs,i,l)=0.0_r8
                  END DO
                END IF
              END DO
            ELSE
              IF (Master) WRITE (stdout,80) 'Ygrid', tSHIPSindx(ng),    &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF

# ifdef SOLVE3D
          IF (exit_flag.eq.NoError) THEN
            status=nf90_get_var(ncSHIPSid(ng), shipsVid(idZgrd,ng),     &
     &                          Tinp, start, total)
            IF (status.eq.nf90_noerr) THEN
              DO l=1,Nships(ng)
                IF ((Tinp(l).gt.REAL(N(ng),r8)).or.                     &
     &              (Tinp(l).lt.0.0_r8)) THEN
                  SHIPS(ng)%bounded(l)=.FALSE.
                ELSE
                  SHIPS(ng)%bounded(l)=.TRUE.
                  DO i=0,NFT
                    SHIPS(ng)%track(izgrd,i,l)=Tinp(l)
                    SHIPS(ng)%track(izrhs,i,l)=0.0_r8
                  END DO
                END IF
              END DO
            ELSE
              IF (Master) WRITE (stdout,80) 'Zgrid', tSHIPSindx(ng),    &
     &                                      TRIM(ncname)
              exit_flag=3
              ioerror=status
            END IF
          END IF
# endif
        END IF

# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
        IF (exit_flag.ne.NoError) RETURN

# ifdef DISTRIBUTE
!
!  Broadcast restart ship information to all nodes.
!
        Npts=NFV(ng)*(NFT+1)*Nships(ng)
        CALL mp_collect (ng, iNLM, Npts, Fspv, SHIPS(ng)%track)
!
!  Collect the bounded status switch.
!
        Fwrk=Fspv
        DO l=1,Nships(ng)
          IF (SHIPS(ng)%bounded(l)) THEN
            Fwrk(l)=1.0_r8
          END IF
        END DO
        CALL mp_collect (ng, iNLM, Nships(ng), Fspv, Fwrk)
        DO l=1,Nships(ng)
          IF (Fwrk(l).ne.Fspv) THEN
            SHIPS(ng)%bounded(l)=.TRUE.
          ELSE
            SHIPS(ng)%bounded(l)=.FALSE.
          END IF
        END DO
# endif
      END IF QUERY
!
  10  FORMAT (6x,'DEF_SHIPS  - creating ships file: ',a)
  20  FORMAT (6x,'DEF_SHIPS  - inquiring history file: ',a)
  30  FORMAT (/,' DEF_SHIPS - unable to create ship NetCDF',            &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_SHIPS - unable to leave definition mode for',     &
     &        ' file: ', a)
  60  FORMAT (/,' DEF_SHIPS - unable to open ships NetCDF file: ',a)
  70  FORMAT (/,' DEF_SHIPS - unable to find variable: ',a,2x,          &
     &        ' in ships NetCDF file: ',a)
  80  FORMAT (/,' DEF_SHIPS - error while reading variable: ',a,2x,     &
     &        ' at time record = ', i6.6,/,13x,                         &
     &        ' in ships NetCDF file: ',a)

      RETURN
      END SUBROUTINE def_ships
#else
      SUBROUTINE def_ships
      RETURN
      END SUBROUTINE def_ships
#endif
