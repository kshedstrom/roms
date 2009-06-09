#include "cppdefs.h"
#ifdef FILTERED_RST
      SUBROUTINE def_filt (ng)
!
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine creates filter NetCDF file, it defines its             !
!  dimensions, attributes, and variables.                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_filter
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
# endif
      USE def_var_mod, ONLY : def_var
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 25

      logical :: Ldefine, got_var(NV)

      integer :: i, j, ifield, nrec, nvd3, nvd4
      integer :: recdim, status, ifile, varid

      integer :: DimIDs(32), t2dgrd(3), u2dgrd(3), v2dgrd(3)

      integer :: def_info, wrt_info

# ifdef SOLVE3D
      integer :: itrc

      integer :: b3dgrd(4), t3dgrd(4), u3dgrd(4), v3dgrd(4), w3dgrd(4)
# endif

      integer :: def_dim

      real(r8) :: Aval(6)

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: fname, ncname, type
!
!=======================================================================
!  Create a new filter NetCDF file.
!=======================================================================
!
!  Activate creation of filter NetCDF file.  Create a new filter
!  file if during a restart run, the filter filename "FILname"
!  is different than the restart filename "RSTname".  I would prefer to
!  not create a new filter file during a restart run, but simply write to 
!  one assumed to have been created before, but this doesn't work yet.
!
      IF (exit_flag.ne.NoError) RETURN
      Ldefine=.FALSE.
      IF (((nrrec(ng).eq.0).and.(iic(ng).eq.ntstart(ng))).or.           &
     &      ((nrrec(ng).ne.0))) THEN
        Ldefine=.TRUE.
      END IF
!
      DO ifile = 1,nfile
        ncname=FILname(ifile,ng)
        IF (Ldefine) THEN
          CALL netcdf_create(ng, iNLM, TRIM(ncname), ncFILid(ifile,ng))
          IF (exit_flag.ne.NoError) THEN
            IF (Master) WRITE (stdout,30) TRIM(ncname)
            RETURN
          END IF
!
!-----------------------------------------------------------------------
!  Define the dimensions of staggered fields.
!-----------------------------------------------------------------------
!
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'xi_rho', IOBOUNDS(ng)%xi_rho, DimIDs( 1))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'xi_u', IOBOUNDS(ng)%xi_u, DimIDs( 2))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'xi_v', IOBOUNDS(ng)%xi_v, DimIDs( 3))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'xi_psi', IOBOUNDS(ng)%xi_psi, DimIDs( 4))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'eta_rho', IOBOUNDS(ng)%eta_rho, DimIDs( 5))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'eta_u', IOBOUNDS(ng)%eta_u, DimIDs( 6))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'eta_v', IOBOUNDS(ng)%eta_v, DimIDs( 7))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'eta_psi', IOBOUNDS(ng)%eta_psi, DimIDs( 8))
          IF (exit_flag.ne.NoError) RETURN
# if defined WRITE_WATER && defined MASKING
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'xy_rho', IOBOUNDS(ng)%xy_rho, DimIDs(17))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'xy_u', IOBOUNDS(ng)%xy_u, DimIDs(18))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'xy_v', IOBOUNDS(ng)%xy_v, DimIDs(19))
          IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef SOLVE3D
#  if defined WRITE_WATER && defined MASKING
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'xyz_rho', IOBOUNDS(ng)%xy_rho*N(ng), DimIDs(20))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'xyz_u', IOBOUNDS(ng)%xy_u*N(ng), DimIDs(21))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'xyz_v', IOBOUNDS(ng)%xy_v*N(ng), DimIDs(22))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'xyz_w', IOBOUNDS(ng)%xy_rho*(N(ng)+1), DimIDs(23))
          IF (exit_flag.ne.NoError) RETURN
#  endif
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'N', N(ng), DimIDs( 9))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             's_rho', N(ng), DimIDs( 9))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             's_w', N(ng)+1, DimIDs(10))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'tracer', NT(ng), DimIDs(11))
          IF (exit_flag.ne.NoError) RETURN
# endif
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'boundary', 4, DimIDs(14))
          IF (exit_flag.ne.NoError) RETURN
          status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
     &             'filterLevel', nf90_unlimited, DimIDs(32))
          IF (exit_flag.ne.NoError) RETURN
!         status = def_dim(ng, iNLM, ncFILid(ifile,ng), ncname,         &
!     &            TRIM(ADJUSTL(Vname(5,idtime))),                      &
!     &            nf90_unlimited, DimIDs(12))
!         IF (exit_flag.ne.NoError) RETURN
        recdim=DimIDs(32)
!
!  Set number of dimensions for output variables.
!
# if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
# else
        nvd3=3
        nvd4=4
# endif
!
!  Define dimension vectors for staggered tracer type variables.
!
# if defined WRITE_WATER && defined MASKING
        t2dgrd(1)=DimIDs(17)
        t2dgrd(2)=DimIDs(32)
#  ifdef SOLVE3D
        t3dgrd(1)=DimIDs(20)
        t3dgrd(2)=DimIDs(32)
#  endif
# else
        t2dgrd(1)=DimIDs( 1)
        t2dgrd(2)=DimIDs( 5)
        t2dgrd(3)=DimIDs(32)
#  ifdef SOLVE3D
        t3dgrd(1)=DimIDs( 1)
        t3dgrd(2)=DimIDs( 5)
        t3dgrd(3)=DimIDs( 9)
        t3dgrd(4)=DimIDs(32)
#  endif
# endif
!
!  Define dimension vectors for staggered u-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
        u2dgrd(1)=DimIDs(18)
        u2dgrd(2)=DimIDs(32)
#  ifdef SOLVE3D
        u3dgrd(1)=DimIDs(21)
        u3dgrd(2)=DimIDs(32)
#  endif
# else
        u2dgrd(1)=DimIDs( 2)
        u2dgrd(2)=DimIDs( 6)
        u2dgrd(3)=DimIDs(32)
#  ifdef SOLVE3D
        u3dgrd(1)=DimIDs( 2)
        u3dgrd(2)=DimIDs( 6)
        u3dgrd(3)=DimIDs( 9)
        u3dgrd(4)=DimIDs(32)
#  endif
# endif
!
!  Define dimension vectors for staggered v-momemtum type variables.
!
# if defined WRITE_WATER && defined MASKING
        v2dgrd(1)=DimIDs(19)
        v2dgrd(2)=DimIDs(32)
#  ifdef SOLVE3D
        v3dgrd(1)=DimIDs(22)
        v3dgrd(2)=DimIDs(32)
#  endif
# else
        v2dgrd(1)=DimIDs( 3)
        v2dgrd(2)=DimIDs( 7)
        v2dgrd(3)=DimIDs(32)
#  ifdef SOLVE3D
        v3dgrd(1)=DimIDs( 3)
        v3dgrd(2)=DimIDs( 7)
        v3dgrd(3)=DimIDs( 9)
        v3dgrd(4)=DimIDs(32)
#  endif
# endif
# ifdef SOLVE3D
!
!  Define dimension vector for staggered w-momemtum type variables.
!
#  if defined WRITE_WATER && defined MASKING
        w3dgrd(1)=DimIDs(23)
        w3dgrd(2)=DimIDs(32)
#  else
        w3dgrd(1)=DimIDs( 1)
        w3dgrd(2)=DimIDs( 5)
        w3dgrd(3)=DimIDs(10)
        w3dgrd(4)=DimIDs(32)
#  endif
!
# endif
!
!  Initialize unlimited time record dimension.
!
        tFILindx(ng)=0
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!-----------------------------------------------------------------------
!  Define common variables and their attributes.
!-----------------------------------------------------------------------
          IF (Master) WRITE (stdout,10) TRIM(ncname)
!
!  Define model time as filter levels.
!
          Vinfo( 1)='filterLevel'
          Vinfo( 2)='Number of days needed in the time filtering'
          status = def_var(ng, iNLM, ncFILid(ifile,ng),                 &
     &           varid, NF90_INT, 1, (/recdim/), Aval,                  &
     &           Vinfo, ncname, SetParAccess = .FALSE.)
          IF (exit_flag.ne.NoError) RETURN
!
!  Define file type, in lieu of whole def_info.
!
          type='ROMS/TOMS filter status file for restarts'
          IF (Master) THEN
            status = nf90_put_att(ncFILid(ifile,ng), nf90_global,       &
     &               'type', LEN_TRIM(type))
          END IF
!
!  Define ocean_time as scalar.
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
          Vinfo( 4)='365.25_day'
        ELSE IF (time_ref.gt.0.0_r8) THEN
          WRITE (Vinfo( 3),'(a,1x,a)') 'seconds since', TRIM(r_text)
          Vinfo( 4)='standard'
        END IF
        Vinfo(14)=Vname(4,idtime)
        status = def_var(ng, iNLM, ncFILid(ifile,ng),                   &
     &           filVid(idtime,ng), NF_TYPE, 0, (/ recdim /),           &
     &           Aval, Vinfo, ncname, SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define filter counts.
!
        Vinfo( 1)='fcount'
        Vinfo( 2)=Vname(2,idtime)
        status = def_var(ng, iNLM, ncFILid(ifile,ng),                   &
     &           varid, NF90_INT,                                       &
     &                 1, (/recdim/), Aval, Vinfo, ncname,              &
     &                 SetParAccess = .FALSE.)
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Define variables and their attributes to each filter file separately.
!-----------------------------------------------------------------------
      IF (Ldefine) THEN
!
!  Define free-surface.
!
        ncname=FILname(1,ng)
        Vinfo( 1)=Vname(1,idFsur)
        Vinfo( 2)=Vname(2,idFsur)
        Vinfo( 3)=Vname(3,idFsur)
        Vinfo(14)=Vname(4,idFsur)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idFsur,ng),r8)
        status=def_var(ng, iNLM, ncFILid(1,ng), filVid(idFsur,ng),      &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 2D momentum in the XI-direction.
!
        Vinfo( 1)=Vname(1,idUbar)
        Vinfo( 2)=Vname(2,idUbar)
        Vinfo( 3)=Vname(3,idUbar)
        Vinfo(14)=Vname(4,idUbar)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUbar,ng),r8)
        status=def_var(ng, iNLM, ncFILid(1,ng), filVid(idUbar,ng),      &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 2D momentum in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVbar)
        Vinfo( 2)=Vname(2,idVbar)
        Vinfo( 3)=Vname(3,idVbar)
        Vinfo(14)=Vname(4,idVbar)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVbar,ng),r8)
        status=def_var(ng, iNLM, ncFILid(1,ng), filVid(idVbar,ng),      &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# ifdef SOLVE3D
!
!  Define 3D momentum component in the XI-direction.
!
        ncname=FILname(2,ng)
        Vinfo( 1)=Vname(1,idUvel)
        Vinfo( 2)=Vname(2,idUvel)
        Vinfo( 3)=Vname(3,idUvel)
        Vinfo(14)=Vname(4,idUvel)
        Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUvel,ng),r8)
        status=def_var(ng, iNLM, ncFILid(2,ng), filVid(idUvel,ng),      &
     &                 NF_FOUT, nvd4, u3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 3D momentum component in the ETA-direction.
!
        ncname=FILname(3,ng)
        Vinfo( 1)=Vname(1,idVvel)
        Vinfo( 2)=Vname(2,idVvel)
        Vinfo( 3)=Vname(3,idVvel)
        Vinfo(14)=Vname(4,idVvel)
        Vinfo(16)='filter level'
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvel,ng),r8)
        status=def_var(ng, iNLM, ncFILid(3,ng), filVid(idVvel,ng),      &
     &                 NF_FOUT, nvd4, v3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  ifndef FILTRIM
!
!  Define density anomaly.
!
        ncname=FILname(7,ng)
        Vinfo( 1)=Vname(1,idDano)
        Vinfo( 2)=Vname(2,idDano)
        Vinfo( 3)=Vname(3,idDano)
        Vinfo(14)=Vname(4,idDano)
        Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idDano,ng),r8)
        status=def_var(ng, iNLM, ncFILid(7,ng), filVid(idDano,ng),      &
     &                 NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define tracer type variables.
!
        ncname=FILname(8,ng)
        DO itrc=1,NAT
          Vinfo( 1)=Vname(1,idTvar(itrc))
          Vinfo( 2)=Vname(2,idTvar(itrc))
          Vinfo( 3)=Vname(3,idTvar(itrc))
          Vinfo(14)=Vname(4,idTvar(itrc))
          Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
          Vinfo(20)='mask_rho'
#   endif
          Vinfo(22)='coordinates'
          Aval(5)=REAL(r3dvar,r8)
          status=def_var(ng, iNLM, ncFILid(8,ng), filTid(itrc,ng),      &
     &                 NF_FOUT, nvd4, t3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
        END DO
!
!  Define vertical (omega) velocity.
!
        ncname=FILname(9,ng)
        Vinfo( 1)=Vname(1,idOvel)
        Vinfo( 2)=Vname(2,idOvel)
        Vinfo( 3)=Vname(3,idOvel)
        Vinfo(14)=Vname(4,idOvel)
        Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idOvel,ng),r8)
        status=def_var(ng, iNLM, ncFILid(9,ng), filVid(idOvel,ng),      &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  ifdef LMD_SKPP
!
!  Define depth of surface boundary layer.
!
        ncname=FILname(1,ng)
        Vinfo( 1)=Vname(1,idHsbl)
        Vinfo( 2)=Vname(2,idHsbl)
        Vinfo( 3)=Vname(3,idHsbl)
        Vinfo(14)=Vname(4,idHsbl)
        Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idHsbl,ng),r8)
        status=def_var(ng, iNLM, ncFILid(1,ng), filVid(idHsbl,ng),      &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  ifdef LMD_BKPP
!
!  Define depth of bottom boundary layer.
!
        Vinfo( 1)=Vname(1,idHbbl)
        Vinfo( 2)=Vname(2,idHbbl)
        Vinfo( 3)=Vname(3,idHbbl)
        Vinfo(14)=Vname(4,idHbbl)
        Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idHbbl,ng),r8)
        status=def_var(ng, iNLM, ncFILid(1,ng), filVid(idHbbl,ng),      &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  ifdef AVERAGES_AKV
!
!  Define vertical viscosity coefficient.
!
        ncname=FILname(6,ng)
        Vinfo( 1)=Vname(1,idVvis)
        Vinfo( 2)=Vname(2,idVvis)
        Vinfo( 3)=Vname(3,idVvis)
        Vinfo(14)=Vname(4,idVvis)
        Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVvis,ng),r8)
        status=def_var(ng, iNLM, ncFILid(6,ng), filVid(idVvis,ng),      &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  ifdef AVERAGES_AKT
!
!  Define vertical diffusion coefficient for potential temperature.
!
        ncname=FILname(4,ng)
        Vinfo( 1)=Vname(1,idTdif)
        Vinfo( 2)=Vname(2,idTdif)
        Vinfo( 3)=Vname(3,idTdif)
        Vinfo(14)=Vname(4,idTdif)
        Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTdif,ng),r8)
        status=def_var(ng, iNLM, ncFILid(4,ng), filVid(idTdif,ng),      &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  ifdef AVERAGES_AKS
!
!  Define vertical diffusion coefficient for salinity.
!
        ncname=FILname(5,ng)
        Vinfo( 1)=Vname(1,idSdif)
        Vinfo( 2)=Vname(2,idSdif)
        Vinfo( 3)=Vname(3,idSdif)
        Vinfo(14)=Vname(4,idSdif)
        Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSdif,ng),r8)
        status=def_var(ng, iNLM, ncFILid(5,ng), filVid(idSdif,ng),      &
     &                 NF_FOUT, nvd4, w3dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#  endif
#  ifdef AVERAGES_FLUXES
!
!  Define surface net heat flux.
!
        ncname=FILname(1,ng)
        Vinfo( 1)=Vname(1,idTsur(itemp))
        Vinfo( 2)=Vname(2,idTsur(itemp))
        Vinfo( 3)=Vname(3,idTsur(itemp))
        Vinfo(11)='upward flux, cooling'
        Vinfo(12)='downward flux, heating'
        Vinfo(14)=Vname(4,idTsur(itemp))
        Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTsur(itemp),ng),r8)
        status=def_var(ng, iNLM, ncFILid(1,ng),                         &
     &                 filVid(idTsur(itemp),ng),                        &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define surface net salt flux.
!
        ncname=FILname(1,ng)
        Vinfo( 1)=Vname(1,idTsur(isalt))
        Vinfo( 2)=Vname(2,idTsur(isalt))
        Vinfo( 3)=Vname(3,idTsur(isalt))
        Vinfo(11)='upward flux, rain'
        Vinfo(12)='downward flux, evap'
        Vinfo(14)=Vname(4,idTsur(isalt))
        Vinfo(16)='filter level'
#   if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#   endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTsur(isalt),ng),r8)
        status=def_var(ng, iNLM, ncFILid(1,ng),                         &
     &                 filVid(idTsur(isalt),ng),                        &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#   ifdef BULK_FLUXES
!
!  Define latent heat flux.
!
        Vinfo( 1)=Vname(1,idLhea)
        Vinfo( 2)=Vname(2,idLhea)
        Vinfo( 3)=Vname(3,idLhea)
        Vinfo(11)='upward flux, cooling'
        Vinfo(12)='downward flux, heating'
        Vinfo(14)=Vname(4,idLhea)
        Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#    endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idLhea,ng),r8)
        status=def_var(ng, iNLM, ncFILid(1,ng), filVid(idLhea,ng),      &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define sensible heat flux.
!
        Vinfo( 1)=Vname(1,idShea)
        Vinfo( 2)=Vname(2,idShea)
        Vinfo( 3)=Vname(3,idShea)
        Vinfo(11)='upward flux, cooling'
        Vinfo(12)='downward flux, heating'
        Vinfo(14)=Vname(4,idShea)
        Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#    endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idShea,ng),r8)
        status=def_var(ng, iNLM, ncFILid(1,ng), filVid(idShea,ng),      &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define longwave radiation flux.
!
        Vinfo( 1)=Vname(1,idLrad)
        Vinfo( 2)=Vname(2,idLrad)
        Vinfo( 3)=Vname(3,idLrad)
        Vinfo(11)='upward flux, cooling'
        Vinfo(12)='downward flux, heating'
        Vinfo(14)=Vname(4,idLrad)
        Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#    endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idLrad,ng),r8)
        status=def_var(ng, iNLM, ncFILid(1,ng), filVid(idLrad,ng),      &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#   endif
#   ifdef SHORTWAVE
!
!  Define shortwave radiation flux.
!
        Vinfo( 1)=Vname(1,idSrad)
        Vinfo( 2)=Vname(2,idSrad)
        Vinfo( 3)=Vname(3,idSrad)
        Vinfo(11)='upward flux, cooling'
        Vinfo(12)='downward flux, heating'
        Vinfo(14)=Vname(4,idSrad)
        Vinfo(16)=Vname(1,idtime)
#    if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
#    endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSrad,ng),r8)
        status=def_var(ng, iNLM, ncFILid(1,ng), filVid(idSrad,ng),      &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
!
!  Define surface u-momentum stress.
!
        Vinfo( 1)=Vname(1,idUsms)
        Vinfo( 2)=Vname(2,idUsms)
        Vinfo( 3)=Vname(3,idUsms)
        Vinfo(14)=Vname(4,idUsms)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUsms,ng),r8)
        status=def_var(ng, iNLM, ncFILid(1,ng), filVid(idUsms,ng),      &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define surface v-momentum stress.
!
        Vinfo( 1)=Vname(1,idVsms)
        Vinfo( 2)=Vname(2,idVsms)
        Vinfo( 3)=Vname(3,idVsms)
        Vinfo(14)=Vname(4,idVsms)
        Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
#  endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVsms,ng),r8)
        status=def_var(ng, iNLM, ncFILid(1,ng), filVid(idVsms,ng),      &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
# ifdef ICE_MODEL
!
!  Define 2D ice momentum in the XI-direction.
!
        ncname=FILname(10,ng)
        Vinfo( 1)=Vname(1,idUice)
        Vinfo( 2)=Vname(2,idUice)
        Vinfo( 3)=Vname(3,idUice)
        Vinfo(14)=Vname(4,idUice)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_u'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idUice,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idUice,ng),     &
     &                 NF_FOUT, nvd3, u2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define 2D ice momentum in the ETA-direction.
!
        Vinfo( 1)=Vname(1,idVice)
        Vinfo( 2)=Vname(2,idVice)
        Vinfo( 3)=Vname(3,idVice)
        Vinfo(14)=Vname(4,idVice)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_v'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idVice,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idVice,ng),     &
     &                 NF_FOUT, nvd3, v2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice concentration
!
        Vinfo( 1)=Vname(1,idAice)
        Vinfo( 2)=Vname(2,idAice)
        Vinfo( 3)=Vname(3,idAice)
        Vinfo(14)=Vname(4,idAice)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idAice,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idAice,ng),     &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice thickness
!
        Vinfo( 1)=Vname(1,idHice)
        Vinfo( 2)=Vname(2,idHice)
        Vinfo( 3)=Vname(3,idHice)
        Vinfo(14)=Vname(4,idHice)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idHice,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idHice,ng),     &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define snow thickness
!
        Vinfo( 1)=Vname(1,idHsno)
        Vinfo( 2)=Vname(2,idHsno)
        Vinfo( 3)=Vname(3,idHsno)
        Vinfo(14)=Vname(4,idHsno)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idHsno,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idHsno,ng),     &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define surface melt ponds
!
        Vinfo( 1)=Vname(1,idSfwat)
        Vinfo( 2)=Vname(2,idSfwat)
        Vinfo( 3)=Vname(3,idSfwat)
        Vinfo(14)=Vname(4,idSfwat)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSfwat,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idSfwat,ng),    &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice surface temperature
!
        Vinfo( 1)=Vname(1,idTice)
        Vinfo( 2)=Vname(2,idTice)
        Vinfo( 3)=Vname(3,idTice)
        Vinfo(14)=Vname(4,idTice)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTice,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idTice,ng),     &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define internal ice temperature
!
        Vinfo( 1)=Vname(1,idTimid)
        Vinfo( 2)=Vname(2,idTimid)
        Vinfo( 3)=Vname(3,idTimid)
        Vinfo(14)=Vname(4,idTimid)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTimid,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idTimid,ng),    &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define under-ice salinity
!
        Vinfo( 1)=Vname(1,idS0mk)
        Vinfo( 2)=Vname(2,idS0mk)
        Vinfo( 3)=Vname(3,idS0mk)
        Vinfo(14)=Vname(4,idS0mk)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idS0mk,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idS0mk,ng),     &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define under-ice temperature
!
        Vinfo( 1)=Vname(1,idT0mk)
        Vinfo( 2)=Vname(2,idT0mk)
        Vinfo( 3)=Vname(3,idT0mk)
        Vinfo(14)=Vname(4,idT0mk)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idT0mk,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idT0mk,ng),     &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define age of ice
!
        Vinfo( 1)=Vname(1,idAgeice)
        Vinfo( 2)=Vname(2,idAgeice)
        Vinfo( 3)=Vname(3,idAgeice)
        Vinfo(14)=Vname(4,idAgeice)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idAgeice,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idAgeice,ng),   &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define sig11 component of internal ice stress
!
        Vinfo( 1)=Vname(1,idSig11)
        Vinfo( 2)=Vname(2,idSig11)
        Vinfo( 3)=Vname(3,idSig11)
        Vinfo(14)=Vname(4,idSig11)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSig11,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idSig11,ng),    &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define sig12 component of internal ice stress
!
        Vinfo( 1)=Vname(1,idSig12)
        Vinfo( 2)=Vname(2,idSig12)
        Vinfo( 3)=Vname(3,idSig12)
        Vinfo(14)=Vname(4,idSig12)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSig12,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idSig12,ng),    &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define sig22 component of internal ice stress
!
        Vinfo( 1)=Vname(1,idSig22)
        Vinfo( 2)=Vname(2,idSig22)
        Vinfo( 3)=Vname(3,idSig22)
        Vinfo(14)=Vname(4,idSig22)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idSig22,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idSig22,ng),    &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice-ocean mass flux
!
        Vinfo( 1)=Vname(1,idIomflx)
        Vinfo( 2)=Vname(2,idIomflx)
        Vinfo( 3)=Vname(3,idIomflx)
        Vinfo(14)=Vname(4,idIomflx)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idIomflx,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idIomflx,ng),   &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice-water friction velocity
!
        Vinfo( 1)=Vname(1,idTauiw)
        Vinfo( 2)=Vname(2,idTauiw)
        Vinfo( 3)=Vname(3,idTauiw)
        Vinfo(14)=Vname(4,idTauiw)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idTauiw,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idTauiw,ng),    &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice-water momentum transfer coefficient
!
        Vinfo( 1)=Vname(1,idChuiw)
        Vinfo( 2)=Vname(2,idChuiw)
        Vinfo( 3)=Vname(3,idChuiw)
        Vinfo(14)=Vname(4,idChuiw)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idChuiw,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idChuiw,ng),    &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define frazil ice growth rate.
!
        Vinfo( 1)=Vname(1,idWfr)
        Vinfo( 2)=Vname(2,idWfr)
        Vinfo( 3)=Vname(3,idWfr)
        Vinfo(14)=Vname(4,idWfr)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idWfr,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idWfr,ng),      &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice-air melt rate.
!
        Vinfo( 1)=Vname(1,idWai)
        Vinfo( 2)=Vname(2,idWai)
        Vinfo( 3)=Vname(3,idWai)
        Vinfo(14)=Vname(4,idWai)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idWai,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idWai,ng),      &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ocean-air ice growth rate.
!
        Vinfo( 1)=Vname(1,idWao)
        Vinfo( 2)=Vname(2,idWao)
        Vinfo( 3)=Vname(3,idWao)
        Vinfo(14)=Vname(4,idWao)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idWao,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idWao,ng),      &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice-ocean ice melt/growth rate.
!
        Vinfo( 1)=Vname(1,idWio)
        Vinfo( 2)=Vname(2,idWio)
        Vinfo( 3)=Vname(3,idWio)
        Vinfo(14)=Vname(4,idWio)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idWio,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idWio,ng),      &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ice melt runoff rate.
!
        Vinfo( 1)=Vname(1,idWro)
        Vinfo( 2)=Vname(2,idWro)
        Vinfo( 3)=Vname(3,idWro)
        Vinfo(14)=Vname(4,idWro)
        Vinfo(16)='filter level'
# if defined WRITE_WATER && defined MASKING
        Vinfo(20)='mask_rho'
# endif
        Vinfo(22)='coordinates'
        Aval(5)=REAL(Iinfo(1,idWro,ng),r8)
        status=def_var(ng, iNLM, ncFILid(10,ng), filVid(idWro,ng),      &
     &                 NF_FOUT, nvd3, t2dgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        DO ifile=1,nfile
          CALL netcdf_enddef (ng, iNLM, FILname(ifile,ng),              &
     &                        ncFILid(ifile,ng))
          IF (exit_flag.ne.NoError) RETURN
        END DO
      END IF
!
!=======================================================================
!  Open an existing filter files, and prepare for appending data.
!=======================================================================
!
      IF (.not.Ldefine) THEN
!
!  Inquire about the contents of filter NetCDF files:  Don't inquire
!  about the dimensions and variables.  Check for consistency.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO

        DO ifile=1,nfile
          ncname=FILname(ifile,ng)
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

          IF (Master) WRITE (stdout, 20) TRIM(ncname)
          ifield=0
!
!  Open filter files for read/write.
!
          CALL netcdf_open (ng, iNLM, ncname, 1, ncFILid(ifile,ng))
          IF (exit_flag.ne.NoError) THEN
            IF (Master) WRITE (stdout,50) TRIM(ncname)
            RETURN
          END IF
#ifdef FOOO
          DO i=1,n_var
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFltLvl))) THEN
              got_var(idFltLvl)=.TRUE.
              filVid(idFltLvl,ng) = var_id(i)
            END IF
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idFltLvl))) THEN
              got_var(idUvel)=.TRUE.
              filVid(idUvel,ng) = var_id(i)
            END IF
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
              got_var(idUvel)=.TRUE.
              filVid(idUvel,ng) = var_id(i)
            END IF
#ifdef SOLVE3D
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idUvel))) THEN
              got_var(idUvel)=.TRUE.
              filVid(idUvel,ng) = var_id(i)
            END IF
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvel))) THEN
              got_var(idVvel)=.TRUE.
              filVid(idVvel,ng) = var_id(i)
            END IF
            DO itrc=1,NAT
              IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
                got_var(idTvar(itrc))=.TRUE.
                filTid(itrc,ng) = var_id(i)
              END IF
            END DO
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
              got_var(idDano)=.TRUE.
              filVid(idDano,ng) = var_id(i)
            END IF
# ifdef LMD_SKPP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHsbl))) THEN
              got_var(idHsbl)=.TRUE.
              filVid(idHsbl,ng) = var_id(i)
            END IF
# endif
# ifdef LMD_BKPP
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idHbbl))) THEN
              got_var(idHbbl)=.TRUE.
              filVid(idHbbl,ng) = var_id(i)
            END IF
# endif
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idVvis))) THEN
              got_var(idVvis)=.TRUE.
              filVid(idVvis,ng) = var_id(i)
            END IF
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTdif))) THEN
              got_var(idTdif)=.TRUE.
              filVid(idTdif,ng) = var_id(i)
            END IF
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idSdif))) THEN
              got_var(idSdif)=.TRUE.
              filVid(idSdif,ng) = var_id(i)
            END IF
#endif
          END DO
        END DO
!
!  Check if initialization variables are available in input NetCDF
!  file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFsur)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idFsur)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idUbar)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVbar)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVbar)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef SOLVE3D
        IF (.not.got_var(idUvel)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idUvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idVvel)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idVvel)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# ifndef FILTRIM
        DO itrc=1,NAT
          IF (.not.got_var(idTvar(itrc))) THEN
            IF (Master) WRITE (stdout,40) TRIM(Vname(1,idTvar(itrc))),  &
     &                                  TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        IF (.not.got_var(idDano)) THEN
          IF (Master) WRITE (stdout,40) TRIM(Vname(1,idDano)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
!
!  Set unlimited time record dimension to current value.
!
!     tFILindx(ng)=tRSTindx(ng)
!
#endif /* FOO */
        END DO
      END IF

  10  FORMAT (6x,'DEF_FILT  - creating filter file:   ',a)
  20  FORMAT (6x,'DEF_FILT  - inquiring filter file: ',a)
  30  FORMAT (/,' DEF_FILT - unable to create filter NetCDF file: ',a)
  40  FORMAT (/,' DEF_FILT - unable to find variable: ',a,2x,           &
     &          ' in filter NetCDF file: ',a)
  50  FORMAT (/,' DEF_FILT - unable to open filter NetCDF file: ',a)
      RETURN
      END SUBROUTINE def_filt
#else
      SUBROUTINE def_filt
      RETURN
      END SUBROUTINE def_filt
#endif
