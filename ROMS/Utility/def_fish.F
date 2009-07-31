#include "cppdefs.h"
#ifdef NEMURO_SAN
      SUBROUTINE def_fish (ng, ldef)
!
!svn $Id: def_fish.F 975 2009-05-05 22:51:13Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates FISH NetCDF file, it defines dimensions,       !
!  attributes, and variables.                                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fish
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

      integer :: fishdim, i, itrc, j, l
      integer :: recdim, status

      integer :: DimIDs(31), fgrd(2), start(2), total(2)
      integer :: Vsize(4)

      integer :: def_dim, varid

      real(r8) :: Aval(6), Tinp(Nfish(ng))

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: ncname
!
      SourceFile='def_fish.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=FISHname(ng)
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
!  Create a new fish data file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), ncFISHid(ng))
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
        status=def_dim(ng, iNLM, ncFISHid(ng), ncname, 's_rho',         &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncFISHid(ng), ncname, 's_w',           &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncFISHid(ng), ncname, 'tracer',        &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iNLM, ncFISHid(ng), ncname, 'Nbed',          &
     &                 Nbed, DimIDs(16))
        IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iNLM, ncFISHid(ng), ncname, 'Nphy',          &
     &                 Nphy, DimIDs(25))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncFISHid(ng), ncname, 'Nbac',          &
     &                 Nbac, DimIDs(26))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncFISHid(ng), ncname, 'Ndom',          &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncFISHid(ng), ncname, 'Nfec',          &
     &                 Nfec, DimIDs(28))
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

        status=def_dim(ng, iNLM, ncFISHid(ng), ncname, 'drifter' ,      &
     &                 Nfish(ng), DimIDs(15))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, ncFISHid(ng), ncname, 'boundary',      &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iNLM, ncFISHid(ng), ncname, 'Nstate',        &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iNLM, ncFISHid(ng), ncname,                  &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

        recdim=DimIDs(12)
        fishdim=DimIDs(15)
!
!  Define dimension vectors for point variables.
!
        fgrd(1)=DimIDs(15)
        fgrd(2)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        tFISHindx(ng)=0
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
        CALL def_info (ng, iNLM, ncFISHid(ng), ncname, DimIDs)
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
        status=def_var(ng, iNLM, ncFISHid(ng), fishVid(idtime,ng),      &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish types
!
        Vinfo( 1)='Ftype'
        Vinfo( 2)='Float type'
        Vinfo( 9)='3-D Lagrangian'
        Vinfo(10)='constant depth'
        status=def_var(ng, iNLM, ncFISHid(ng), varid, nf90_int,         &
     &                 1, (/fishdim/), Aval, Vinfo, ncname)
!
!  Define fish X-grid locations.
!
        Vinfo( 1)='Xgrid'
        Vinfo( 2)='x-grid fish locations'
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
        status=def_var(ng, iNLM, ncFISHid(ng), fishVid(idXgrd,ng),      &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish Y-grid locations.
!
        Vinfo( 1)='Ygrid'
        Vinfo( 2)='Y-grid fish locations'
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
        status=def_var(ng, iNLM, ncFISHid(ng), fishVid(idYgrd,ng),      &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Define fish Z-grid locations.
!
        Vinfo( 1)='Zgrid'
        Vinfo( 2)='Z-grid fish locations'
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
        status=def_var(ng, iNLM, ncFISHid(ng), fishVid(idZgrd,ng),      &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Define fish (lon,lat) or (x,y) locations.
!
        IF (spherical) THEN
          Vinfo( 1)='lon'
          Vinfo( 2)='longitude of fish trajectories'
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
          status=def_var(ng, iNLM, ncFISHid(ng), fishVid(idglon,ng),    &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN

          Vinfo( 1)='lat'
          Vinfo( 2)='latitude of fish trajectories'
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
          status=def_var(ng, iNLM, ncFISHid(ng), fishVid(idglat,ng),    &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
       ELSE
          Vinfo( 1)='x'
          Vinfo( 2)='x-location of fish trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='x, scalar, series'
          Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
# endif
          status=def_var(ng, iNLM, ncFISHid(ng), fishVid(idglon,ng),    &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN

          Vinfo( 1)='y'
          Vinfo( 2)='y-location of fish trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='y, scalar, series'
          Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
# endif
          status=def_var(ng, iNLM, ncFISHid(ng), fishVid(idglat,ng),    &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

# ifdef SOLVE3D
!
!  Define fish depths.
!
        Vinfo( 1)='depth'
        Vinfo( 2)='depth of fish trajectories'
        Vinfo( 3)='meter'
        Vinfo(14)='depth, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, ncFISHid(ng), fishVid(iddpth,ng),      &
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
        status=def_var(ng, iNLM, ncFISHid(ng), fishVid(idDano,ng),      &
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
          status=def_var(ng, iNLM, ncFISHid(ng), fishTid(itrc,ng),      &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
# endif
# ifdef NEMURO_SAN
!  Define fish wet weight
        Vinfo( 1)='fwwt'
        Vinfo( 2)='fish wet weight'
        Vinfo( 3)='g'
        Vinfo(14)='fwwt, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, ncFISHid(ng),fishVid(idfwwt,ng),       &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!  Define fish worth
        Vinfo( 1)='fworth'
        Vinfo( 2)='fish worth'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='fworth, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, ncFISHid(ng),fishVid(idfworth,ng),     &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!  Define fish P-value
        Vinfo( 1)='fpval'
        Vinfo( 2)='fish p-value'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='fpval, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, ncFISHid(ng),fishVid(idfpval,ng),      &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!  Define fish spawning location
        Vinfo( 1)='fspwnloc'
        Vinfo( 2)='fish spawning location'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='fspwnloc, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, ncFISHid(ng),fishVid(idspwnloc,ng),    &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Initialize unlimited time record dimension.
!
        tFISHindx(ng)=0
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, ncFISHid(ng))
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, ncFISHid(ng), ncname)
        IF (exit_flag.ne.NoError) RETURN
        IF (OutThread) THEN
          status=nf90_inq_varid(ncFISHid(ng),'Ftype',varid)
          status=nf90_put_var(ncFISHid(ng), varid, FISHES(ng)%Ftype)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,70) 'Ftype', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END IF
      END IF DEFINE
!
!=======================================================================
!  Open an existing fish file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=FISHname(ng)
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
!           Nfish(ng)=dim_size(i)
!           EXIT
!         END IF
!       END DO
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open fish file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, ncFISHid(ng))
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
!  fish variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            fishVid(idtime,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Xgrid') THEN
            got_var(idXgrd)=.TRUE.
            fishVid(idXgrd,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Ygrid') THEN
            got_var(idYgrd)=.TRUE.
            fishVid(idYgrd,ng)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.'Zgrid') THEN
            got_var(idZgrd)=.TRUE.
            fishVid(idZgrd,ng)=var_id(i)
# endif
          ELSE IF (spherical.and.TRIM(var_name(i)).eq.'lon') THEN
            got_var(idglon)=.TRUE.
            fishVid(idglon,ng)=var_id(i)
          ELSE IF (spherical.and.TRIM(var_name(i)).eq.'lat') THEN
            got_var(idglat)=.TRUE.
            fishVid(idglat,ng)=var_id(i)
          ELSE IF (.not.spherical.and.TRIM(var_name(i)).eq.'x') THEN
            got_var(idglon)=.TRUE.
            fishVid(idglon,ng)=var_id(i)
          ELSE IF (.not.spherical.and.TRIM(var_name(i)).eq.'y') THEN
            got_var(idglat)=.TRUE.
            fishVid(idglat,ng)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.'depth') THEN
            got_var(iddpth)=.TRUE.
            fishVid(iddpth,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            fishVid(idDano,ng)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idfwwt))) THEN
            got_var(idfwwt)=.TRUE.
            fishVid(idfwwt,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idfworth))) THEN
            got_var(idfworth)=.TRUE.
            fishVid(idfworth,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idfpval))) THEN
            got_var(idfpval)=.TRUE.
            fishVid(idfpval,ng)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idspwnloc))) THEN
            got_var(idspwnloc)=.TRUE.
            fishVid(idspwnloc,ng)=var_id(i)
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              fishTid(itrc,ng)=var_id(i)
            END IF
          END DO
# endif
        END DO
!
!  Check if fish variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idXgrd)) THEN
          IF (Master) WRITE (stdout,60) 'Xgrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idYgrd)) THEN
          IF (Master) WRITE (stdout,60) 'Ygrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(idZgrd)) THEN
          IF (Master) WRITE (stdout,60) 'Zgrid', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
        IF (.not.got_var(idglon)) THEN
          IF (spherical) THEN
            IF (Master) WRITE (stdout,60) 'lon', TRIM(ncname)
          ELSE
            IF (Master) WRITE (stdout,60) 'x', TRIM(ncname)
          END IF
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idglat)) THEN
          IF (spherical) THEN
            IF (Master) WRITE (stdout,60) 'lat', TRIM(ncname)
          ELSE
            IF (Master) WRITE (stdout,60) 'y', TRIM(ncname)
          END IF
          exit_flag=3
          RETURN
        END IF
# ifdef SOLVE3D
        IF (.not.got_var(iddpth)) THEN
          IF (Master) WRITE (stdout,60) 'depth', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idDano)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,60) TRIM(Vname(1,idTvar(itrc))),  &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
# endif
! Fish wet weight
        IF (.not.got_var(idfwwt)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idfwwt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
! Fish worth
        IF (.not.got_var(idfworth)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idfworth)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
!-----------------------------------------------------------------------
!  Initialize fish positions to the appropriate values.
!-----------------------------------------------------------------------
!
!  Set-up fish time record.
!
        IF (fishrrec(ng).lt.0) THEN
          tFISHindx(ng)=rec_size
        ELSE
         tFISHindx(ng)=ABS(fishrrec(ng))
        END IF
        NrecFISH(ng)=tFISHindx(ng)
!
!  Read in fish nondimentional horizontal positions.
!
        CALL netcdf_get_fvar (ng, iNLM, ncname, 'Xgrid',                &
     &                        Tinp,                                     &
     &                        ncid = ncFISHid(ng),                      &
     &                        start = (/1,tFISHindx(ng)/),              &
     &                        total = (/Nfish(ng),1/))
        IF (exit_flag.ne.NoError) RETURN

        DO l=1,Nfish(ng)
          IF ((Tinp(l).gt.REAL(Lm(ng)+1,r8)-0.5_r8).or.                 &
     &        (Tinp(l).lt.0.5_r8)) THEN
            FISHES(ng)%bounded(l)=.FALSE.
          ELSE
            FISHES(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              FISHES(ng)%track(ixgrd,i,l)=Tinp(l)
              FISHES(ng)%track(ixrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname, 'Ygrid',                &
     &                        Tinp,                                     &
     &                        ncid = ncFISHid(ng),                      &
     &                        start = (/1,tFISHindx(ng)/),              &
     &                        total = (/Nfish(ng),1/))
        IF (exit_flag.ne.NoError) RETURN

        DO l=1,Nfish(ng)
          IF ((Tinp(l).gt.REAL(Mm(ng)+1,r8)-0.5_r8).or.                 &
     &        (Tinp(l).lt.0.5_r8)) THEN
            FISHES(ng)%bounded(l)=.FALSE.
          ELSE
            FISHES(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              FISHES(ng)%track(iygrd,i,l)=Tinp(l)
              FISHES(ng)%track(iyrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO

# ifdef SOLVE3D
        CALL netcdf_get_fvar (ng, iNLM, ncname, 'Zgrid',                &
     &                        Tinp,                                     &
     &                        ncid = ncFISHid(ng),                      &
     &                        start = (/1,tFISHindx(ng)/),              &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          IF ((Tinp(l).gt.REAL(N(ng),r8)).or.                           &
     &        (Tinp(l).lt.0.0_r8)) THEN
            FISHES(ng)%bounded(l)=.FALSE.
          ELSE
            FISHES(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              FISHES(ng)%track(izgrd,i,l)=Tinp(l)
              FISHES(ng)%track(izrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO
# endif
!
!  Read fish variables
!
        CALL netcdf_get_fvar (ng, iNLM, ncname, 'fwwt',                 &
     &                        Tinp,                                     &
     &                        ncid = ncFISHid(ng),                      &
     &                        start = (/1,tFISHindx(ng)/),              &
     &                        total = (/Nfish(ng),1/))
        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifwwt,l)=Tinp(l)
        END DO
        CALL netcdf_get_fvar (ng, iNLM, ncname, 'fworth',               &
     &                        Tinp,                                     &
     &                        ncid = ncFISHid(ng),                      &
     &                        start = (/1,tFISHindx(ng)/),              &
     &                        total = (/Nfish(ng),1/))
        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifworth,l)=Tinp(l)
        END DO
        CALL netcdf_get_fvar (ng, iNLM, ncname, 'fbatch',               &
     &                        Tinp,                                     &
     &                        ncid = ncFISHid(ng),                      &
     &                        start = (/1,tFISHindx(ng)/),              &
     &                        total = (/Nfish(ng),1/))
        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifbatch,l)=Tinp(l)
        END DO

        IF (exit_flag.ne.NoError) RETURN

      END IF QUERY
!
  10  FORMAT (6x,'DEF_FISH  - creating fish file: ',a)
  20  FORMAT (6x,'DEF_FISH  - inquiring history file: ',a)
  30  FORMAT (/,' DEF_FISH - unable to create fish NetCDF',             &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_FISH - unable to open fish NetCDF file: ',a)
  60  FORMAT (/,' DEF_FISH - unable to find variable: ',a,2x,           &
     &        ' in fish NetCDF file: ',a)
  70  FORMAT (/,' DEF_FISH - unable to write variable: ',a,2x,          &
     &        ' in fish NetCDF file: ',a)

      RETURN
      END SUBROUTINE def_fish
#else
      SUBROUTINE def_fish
      RETURN
      END SUBROUTINE def_fish
#endif
