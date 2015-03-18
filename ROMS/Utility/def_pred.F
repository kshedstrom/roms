#include "cppdefs.h"
#if defined NEMURO_SAN && defined PREDATOR
      SUBROUTINE def_pred (ng, ldef)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine creates PRED NetCDF file, it defines dimensions,       !
!  attributes, and variables.                                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_pred
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

      integer :: preddim, i, itrc, j, l
      integer :: recdim, status, specdim, ycdim

      integer :: DimIDs(33), fgrd(2), start(2), total(2)
      integer :: Vsize(4)

      integer :: def_dim, varid

      real(r8) :: Aval(6), Tinp(Npred(ng))

      character (len=80) :: Vinfo(Natt)
      character (len=80) :: ncname
!
      SourceFile='def_pred.F'
!
!-----------------------------------------------------------------------
!  Set and report file name.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
      ncname=PRED(ng)%name
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
!  Create a new pred data file.
!=======================================================================
!
      DEFINE : IF (ldef) THEN
        CALL netcdf_create (ng, iNLM, TRIM(ncname), PRED(ng)%ncid)
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
        status=def_dim(ng, iNLM, PRED(ng)%ncid, ncname, 's_rho',        &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, PRED(ng)%ncid, ncname, 's_w',          &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, PRED(ng)%ncid, ncname, 'tracer',       &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iNLM, PRED(ng)%ncid, ncname, 'Nbed',         &
     &                 Nbed, DimIDs(16))
        IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iNLM, PRED(ng)%ncid, ncname, 'Nphy',         &
     &                 Nphy, DimIDs(25))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, PRED(ng)%ncid, ncname, 'Nbac',         &
     &                 Nbac, DimIDs(26))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, PRED(ng)%ncid, ncname, 'Ndom',         &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, PRED(ng)%ncid, ncname, 'Nfec',         &
     &                 Nfec, DimIDs(28))
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

        status=def_dim(ng, iNLM, PRED(ng)%ncid, ncname, 'drifter' ,     &
     &                 Npred(ng), DimIDs(15))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, PRED(ng)%ncid, ncname, 'species' ,     &
     &                 Npredspecies(ng), DimIDs(32))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, PRED(ng)%ncid, ncname, 'boundary',     &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iNLM, PRED(ng)%ncid, ncname, 'Nstate',       &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iNLM, PRED(ng)%ncid, ncname,                 &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

        recdim=DimIDs(12)
        preddim=DimIDs(15)
        specdim=DimIDs(32)
        ycdim=DimIDs(33)
!
!  Define dimension vectors for point variables.
!
        fgrd(1)=DimIDs(15)
        fgrd(2)=DimIDs(12)
!
!  Initialize unlimited time record dimension.
!
        PRED(ng)%Rindex=0
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
        CALL def_info (ng, iNLM, PRED(ng)%ncid, ncname, DimIDs)
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
        status=def_var(ng, iNLM, PRED(ng)%ncid, PRED(ng)%Vid(idtime),   &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!*********************************************************************
! This next batch gets written out, but not as a timeseries, they get
! overwritten with the current lifestage, etc., of the pred for restart
! purposes.
!*********************************************************************
!
!  Define pred species
!
        Vinfo( 1)='Species'
        Vinfo( 2)='Pred species'
        status=def_var(ng, iNLM, PRED(ng)%ncid, PRED(ng)%Vid(idPspc),   &
     &                 nf90_int, 1, (/preddim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define pred alive
!
        Vinfo( 1)='Alive'
        Vinfo( 2)='Pred alive'
        Vinfo( 9)='Pred alive'
        Vinfo(10)='Pred dead'
        status=def_var(ng, iNLM, PRED(ng)%ncid, PRED(ng)%Vid(idPalv),   &
     &                 nf90_int, 1, (/preddim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define pred X-grid locations.
!
        Vinfo( 1)='Xgrid'
        Vinfo( 2)='x-grid pred locations'
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
        status=def_var(ng, iNLM, PRED(ng)%ncid, PRED(ng)%Vid(idXgrd),   &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define pred Y-grid locations.
!
        Vinfo( 1)='Ygrid'
        Vinfo( 2)='Y-grid pred locations'
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
        status=def_var(ng, iNLM, PRED(ng)%ncid, PRED(ng)%Vid(idYgrd),   &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Define pred Z-grid locations.
!
        Vinfo( 1)='Zgrid'
        Vinfo( 2)='Z-grid pred locations'
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
        status=def_var(ng, iNLM, PRED(ng)%ncid, PRED(ng)%Vid(idZgrd),   &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Define pred (lon,lat) or (x,y) locations.
!
        IF (spherical) THEN
          Vinfo( 1)='lon'
          Vinfo( 2)='longitude of pred trajectories'
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
          status=def_var(ng, iNLM, PRED(ng)%ncid, PRED(ng)%Vid(idglon), &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN

          Vinfo( 1)='lat'
          Vinfo( 2)='latitude of pred trajectories'
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
          status=def_var(ng, iNLM, PRED(ng)%ncid, PRED(ng)%Vid(idglat), &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
       ELSE
          Vinfo( 1)='x'
          Vinfo( 2)='x-location of pred trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='x, scalar, series'
          Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
# endif
          status=def_var(ng, iNLM, PRED(ng)%ncid, PRED(ng)%Vid(idglon), &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN

          Vinfo( 1)='y'
          Vinfo( 2)='y-location of pred trajectories'
          Vinfo( 3)='meter'
          Vinfo(14)='y, scalar, series'
          Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
          Vinfo(24)='_FillValue'
          Aval(6)=spval
# endif
          status=def_var(ng, iNLM, PRED(ng)%ncid, PRED(ng)%Vid(idglat), &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END IF

# ifdef SOLVE3D
!
!  Define pred depths.
!
        Vinfo( 1)='depth'
        Vinfo( 2)='depth of pred trajectories'
        Vinfo( 3)='meter'
        Vinfo(14)='depth, scalar, series'
        Vinfo(16)=Vname(1,idtime)
#  ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
#  endif
        status=def_var(ng, iNLM, PRED(ng)%ncid, PRED(ng)%Vid(iddpth),   &
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
        status=def_var(ng, iNLM, PRED(ng)%ncid, PRED(ng)%Vid(idDano),   &
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
          status=def_var(ng, iNLM, PRED(ng)%ncid, PRED(ng)%Tid(itrc),   &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
# endif
!
!  Define pred wet weight
!
        Vinfo( 1)='Pwwt'
        Vinfo( 2)='pred wet weight'
        Vinfo( 3)='g'
        Vinfo(14)='pwwt, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, PRED(ng)%ncid,PRED(ng)%Vid(idPwwt),    &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define pred worth
!
        Vinfo( 1)='Pworth'
        Vinfo( 2)='pred worth'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='pworth, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, PRED(ng)%ncid,PRED(ng)%Vid(idPworth),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define pred P-value
!
        Vinfo( 1)='Ppval'
        Vinfo( 2)='pred p-value'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='ppval, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, PRED(ng)%ncid,PRED(ng)%Vid(idPpval),   &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define Pred consumption on species 1
!
        Vinfo( 1)='PcsmF1'
        Vinfo( 2)='pred csmp on sp1'
        Vinfo( 3)='gPrey / gFish / hour'
        Vinfo(14)='pcsmf1, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, PRED(ng)%ncid,PRED(ng)%Vid(idPcsmF1),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define Pred consumption on species 2
!
        Vinfo( 1)='PcsmF2'
        Vinfo( 2)='pred csmp on sp2'
        Vinfo( 3)='gPrey / gFish / hour'
        Vinfo(14)='pcsmf2, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, PRED(ng)%ncid,PRED(ng)%Vid(idPcsmF2),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define Pred consumption on species 3
!
        Vinfo( 1)='PcsmF3'
        Vinfo( 2)='pred csmp on sp3'
        Vinfo( 3)='gPrey / gFish / hour'
        Vinfo(14)='pcsmf3, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, PRED(ng)%ncid,PRED(ng)%Vid(idPcsmF3),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define Pred consumption on species 4
!
        Vinfo( 1)='PcsmF4'
        Vinfo( 2)='pred csmp on sp4'
        Vinfo( 3)='gPrey / gFish / hour'
        Vinfo(14)='pcsmf4, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, PRED(ng)%ncid,PRED(ng)%Vid(idPcsmF4),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define Pred consumption on species 5
!
        Vinfo( 1)='PcsmF5'
        Vinfo( 2)='pred csmp on sp5'
        Vinfo( 3)='gPrey / gFish / hour'
        Vinfo(14)='pcsmf5, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, PRED(ng)%ncid,PRED(ng)%Vid(idPcsmF5),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Initialize unlimited time record dimension.
!
        PRED(ng)%rindex=0
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, PRED(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, PRED(ng)%ncid, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing pred file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=PRED(ng)%name
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
!           Npred(ng)=dim_size(i)
!           EXIT
!         END IF
!       END DO
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, iNLM, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Open pred file for read/write.
!
        CALL netcdf_open (ng, iNLM, ncname, 1, PRED(ng)%ncid)
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
!  pred variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            PRED(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Pwwt') THEN
            got_var(idPwwt)=.TRUE.
            PRED(ng)%Vid(idPwwt)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Pworth') THEN
            got_var(idPworth)=.TRUE.
            PRED(ng)%Vid(idPworth)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Ppval') THEN
            got_var(idPpval)=.TRUE.
            PRED(ng)%Vid(idPpval)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Species') THEN
            got_var(idPspc)=.TRUE.
            PRED(ng)%Vid(idPspc)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Alive') THEN
            got_var(idPalv)=.TRUE.
            PRED(ng)%Vid(idPalv)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Xgrid') THEN
            got_var(idXgrd)=.TRUE.
            PRED(ng)%Vid(idXgrd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Ygrid') THEN
            got_var(idYgrd)=.TRUE.
            PRED(ng)%Vid(idYgrd)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.'Zgrid') THEN
            got_var(idZgrd)=.TRUE.
            PRED(ng)%Vid(idZgrd)=var_id(i)
# endif
          ELSE IF (spherical.and.TRIM(var_name(i)).eq.'lon') THEN
            got_var(idglon)=.TRUE.
            PRED(ng)%Vid(idglon)=var_id(i)
          ELSE IF (spherical.and.TRIM(var_name(i)).eq.'lat') THEN
            got_var(idglat)=.TRUE.
            PRED(ng)%Vid(idglat)=var_id(i)
          ELSE IF (.not.spherical.and.TRIM(var_name(i)).eq.'x') THEN
            got_var(idglon)=.TRUE.
            PRED(ng)%Vid(idglon)=var_id(i)
          ELSE IF (.not.spherical.and.TRIM(var_name(i)).eq.'y') THEN
            got_var(idglat)=.TRUE.
            PRED(ng)%Vid(idglat)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.'depth') THEN
            got_var(iddpth)=.TRUE.
            PRED(ng)%Vid(iddpth)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            PRED(ng)%Vid(idDano)=var_id(i)
# endif
          ELSE IF (TRIM(var_name(i)).eq.'PcsmF1') THEN
            got_var(idPcsmF1)=.TRUE.
            PRED(ng)%Vid(idPcsmF1)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'PcsmF2') THEN
            got_var(idPcsmF2)=.TRUE.
            PRED(ng)%Vid(idPcsmF2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'PcsmF3') THEN
            got_var(idPcsmF3)=.TRUE.
            PRED(ng)%Vid(idPcsmF3)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'PcsmF4') THEN
            got_var(idPcsmF4)=.TRUE.
            PRED(ng)%Vid(idPcsmF4)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'PcsmF5') THEN
            got_var(idPcsmF5)=.TRUE.
            PRED(ng)%Vid(idPcsmF5)=var_id(i)
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              PRED(ng)%Tid(itrc)=var_id(i)
            END IF
          END DO
# endif
        END DO
!
!  Check if pred variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPwwt)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPwwt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPworth)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPworth)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPspc)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPspc)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idPalv)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idPalv)),          &
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
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,iddpth)),          &
     &                                  TRIM(ncname)
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
!
!-----------------------------------------------------------------------
!  Initialize pred positions to the appropriate values.
!-----------------------------------------------------------------------
!
!  Set-up pred time record.
!
        IF (predrrec(ng).lt.0) THEN
          PRED(ng)%rindex=rec_size
        ELSE
         PRED(ng)%rindex=ABS(predrrec(ng))
        END IF
        PRED(ng)%Nrec=PRED(ng)%rindex
!
!  Read in pred nondimentional horizontal positions.
!
        CALL netcdf_get_fvar (ng, iNLM, ncname, 'Xgrid',                &
     &                        Tinp,                                     &
     &                        ncid = PRED(ng)%ncid,                     &
     &                        start = (/1,PRED(ng)%rindex/),            &
     &                        total = (/Npred(ng),1/))
        IF (exit_flag.ne.NoError) RETURN

        DO l=1,Npred(ng)
          IF ((Tinp(l).gt.REAL(Lm(ng)+1,r8)-0.5_r8).or.                 &
     &        (Tinp(l).lt.0.5_r8)) THEN
            PREDS(ng)%bounded(l)=.FALSE.
          ELSE
            PREDS(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              PREDS(ng)%track(ixgrd,i,l)=Tinp(l)
              PREDS(ng)%track(ixrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname, 'Ygrid',                &
     &                        Tinp,                                     &
     &                        ncid = PRED(ng)%ncid,                     &
     &                        start = (/1,PRED(ng)%rindex/),            &
     &                        total = (/Npred(ng),1/))
        IF (exit_flag.ne.NoError) RETURN

        DO l=1,Npred(ng)
          IF ((Tinp(l).gt.REAL(Mm(ng)+1,r8)-0.5_r8).or.                 &
     &        (Tinp(l).lt.0.5_r8)) THEN
            PREDS(ng)%bounded(l)=.FALSE.
          ELSE
            PREDS(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              PREDS(ng)%track(iygrd,i,l)=Tinp(l)
              PREDS(ng)%track(iyrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO

# ifdef SOLVE3D
        CALL netcdf_get_fvar (ng, iNLM, ncname, 'Zgrid',                &
     &                        Tinp,                                     &
     &                        ncid = PRED(ng)%ncid,                     &
     &                        start = (/1,PRED(ng)%rindex/),            &
     &                        total = (/Npred(ng),1/))

        DO l=1,Npred(ng)
          IF ((Tinp(l).gt.REAL(N(ng),r8)).or.                           &
     &        (Tinp(l).lt.0.0_r8)) THEN
            PREDS(ng)%bounded(l)=.FALSE.
          ELSE
            PREDS(ng)%bounded(l)=.TRUE.
            DO i=0,NFT
              PREDS(ng)%track(izgrd,i,l)=Tinp(l)
              PREDS(ng)%track(izrhs,i,l)=0.0_r8
            END DO
          END IF
        END DO
# endif
!
!  Read pred variables
!
        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idPwwt)),                    &
     &                        Tinp,                                     &
     &                        ncid = PRED(ng)%ncid,                     &
     &                        start = (/1,PRED(ng)%rindex/),            &
     &                        total = (/Npred(ng),1/))
        DO l=1,Npred(ng)
          PREDS(ng)%bioenergy(ipwwt,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idPworth)),                  &
     &                        Tinp,                                     &
     &                        ncid = PRED(ng)%ncid,                     &
     &                        start = (/1,PRED(ng)%rindex/),            &
     &                        total = (/Npred(ng),1/))
        DO l=1,Npred(ng)
          PREDS(ng)%bioenergy(ipworth,l)=Tinp(l)
        END DO

        CALL netcdf_get_ivar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idPspc)),                    &
     &                        PREDS(ng)%species,                        &
     &                        ncid = PRED(ng)%ncid,                     &
     &                        start = (/1,PRED(ng)%rindex/),            &
     &                        total = (/Npred(ng),1/))

        CALL netcdf_get_lvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idPalv)),                    &
     &                        PREDS(ng)%alive,                          &
     &                        ncid = PRED(ng)%ncid,                     &
     &                        start = (/1,PRED(ng)%rindex/),            &
     &                        total = (/Npred(ng),1/))

        IF (exit_flag.ne.NoError) RETURN

      END IF QUERY
!
  10  FORMAT (6x,'DEF_PRED  - creating pred file: ',a)
  20  FORMAT (6x,'DEF_PRED  - inquiring pred file: ',a)
  30  FORMAT (/,' DEF_PRED - unable to create pred NetCDF',             &
     &        ' file: ',a)
  40  FORMAT (1pe11.4,1x,'millimeter')
  50  FORMAT (/,' DEF_PRED - unable to open pred NetCDF file: ',a)
  60  FORMAT (/,' DEF_PRED - unable to find variable: ',a,2x,           &
     &        ' in pred NetCDF file: ',a)
  70  FORMAT (/,' DEF_PRED - unable to write variable: ',a,2x,          &
     &        ' in pred NetCDF file: ',a)

      RETURN
      END SUBROUTINE def_pred
#else
      SUBROUTINE def_pred
      RETURN
      END SUBROUTINE def_pred
#endif
