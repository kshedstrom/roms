#include "cppdefs.h"
#ifdef NEMURO_SAN
      SUBROUTINE def_fish (ng, ldef)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
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

      integer :: fishdim, i, itrc, j, l
      integer :: recdim, status, specdim, ycdim

      integer :: DimIDs(33), fgrd(2), start(2), total(2)
      integer :: Vsize(4)

      integer :: def_dim, varid

      integer :: TinpI(Nfish(ng))
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
      ncname=FISH(ng)%name
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
        CALL netcdf_create (ng, iNLM, TRIM(ncname), FISH(ng)%ncid)
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
        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname, 's_rho',        &
     &                 N(ng), DimIDs( 9))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname, 's_w',          &
     &                 N(ng)+1, DimIDs(10))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname, 'tracer',       &
     &                 NT(ng), DimIDs(11))
        IF (exit_flag.ne.NoError) RETURN

#  ifdef SEDIMENT
        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname, 'Nbed',         &
     &                 Nbed, DimIDs(16))
        IF (exit_flag.ne.NoError) RETURN
#  endif

#  ifdef ECOSIM
        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname, 'Nphy',         &
     &                 Nphy, DimIDs(25))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname, 'Nbac',         &
     &                 Nbac, DimIDs(26))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname, 'Ndom',         &
     &                 Ndom, DimIDs(27))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname, 'Nfec',         &
     &                 Nfec, DimIDs(28))
        IF (exit_flag.ne.NoError) RETURN
#  endif
# endif

        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname, 'drifter' ,     &
     &                 Nfish(ng), DimIDs(15))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname, 'species' ,     &
     &                 Nspecies(ng), DimIDs(32))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname, 'yearclass' ,   &
     &                 Nyearclass(ng), DimIDs(33))
        IF (exit_flag.ne.NoError) RETURN

        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname, 'boundary',     &
     &                 4, DimIDs(14))
        IF (exit_flag.ne.NoError) RETURN

# ifdef FOUR_DVAR
        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname, 'Nstate',       &
     &                 NstateVar(ng), DimIDs(29))
        IF (exit_flag.ne.NoError) RETURN
# endif

        status=def_dim(ng, iNLM, FISH(ng)%ncid, ncname,                 &
     &                 TRIM(ADJUSTL(Vname(5,idtime))),                  &
     &                 nf90_unlimited, DimIDs(12))
        IF (exit_flag.ne.NoError) RETURN

        recdim=DimIDs(12)
        fishdim=DimIDs(15)
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
        FISH(ng)%Rindex=0
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
        CALL def_info (ng, iNLM, FISH(ng)%ncid, ncname, DimIDs)
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
        status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idtime),   &
     &                 NF_TYPE, 1, (/recdim/), Aval, Vinfo, ncname,     &
     &                 SetParAccess = .FALSE.)
        IF (exit_flag.ne.NoError) RETURN
!
!*********************************************************************
! This next batch gets written out, but not as a timeseries, they get
! overwritten with the current lifestage, etc., of the fish for restart
! purposes.
!*********************************************************************
!
!  Define fish species
!
        Vinfo( 1)='Species'
        Vinfo( 2)='Fish species'
        status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idFspc),   &
     &                 nf90_int, 1, (/fishdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish behavior
!
        Vinfo( 1)='Hbehave'
        Vinfo( 2)='Fish horizontal behavior'
        status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idFhbhv),  &
     &                 nf90_int, 1, (/fishdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
        Vinfo( 1)='Vbehave'
        Vinfo( 2)='Fish vertical behavior'
        status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idFvbhv),  &
     &                 nf90_int, 1, (/fishdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish alive
!
        Vinfo( 1)='Alive'
        Vinfo( 2)='Fish alive'
        Vinfo( 9)='Fish alive'
        Vinfo(10)='Fish dead'
        status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idFalv),   &
     &                 nf90_int, 1, (/fishdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define time since last spawn (needed on restart)
!
        Vinfo( 1)='Tspawn'
        Vinfo( 2)='Time since last spawn'
        status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idFTspwn), &
     &                 NF_FOUT, 1, (/fishdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define initial worth (needed on restart)
!
        Vinfo( 1)='Iniworth'
        Vinfo( 2)='Initial worth'
        status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idFiniw),  &
     &                 NF_FOUT, 1, (/fishdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish age base
!
        Vinfo( 1)='Fish_age_base'
        Vinfo( 2)='Fish age base'
        status=def_var(ng, iNLM, FISH(ng)%ncid, varid, nf90_int,        &
     &           1, (/specdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish species base
!
        Vinfo( 1)='Fish_spec_base'
        Vinfo( 2)='Fish species base'
        status=def_var(ng, iNLM, FISH(ng)%ncid, varid, nf90_int,        &
     &           1, (/specdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish num free
!
        Vinfo( 1)='Fish_num_free'
        Vinfo( 2)='Number of free superindividuals'
        status=def_var(ng, iNLM, FISH(ng)%ncid, varid, nf90_int,        &
     &           1, (/specdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish next free
!
        Vinfo( 1)='Fish_next_free'
        Vinfo( 2)='Next free superindividual'
        status=def_var(ng, iNLM, FISH(ng)%ncid, varid, nf90_int,        &
     &           1, (/specdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish number per spawn
!
        Vinfo( 1)='Fish_num_super'
        Vinfo( 2)='Number superindividuals per spawn'
        status=def_var(ng, iNLM, FISH(ng)%ncid, varid, nf90_int,        &
     &           1, (/specdim/), Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
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
        status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idXgrd),   &
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
        status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idYgrd),   &
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
        status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idZgrd),   &
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
          status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idglon), &
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
          status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idglat), &
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
          status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idglon), &
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
          status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idglat), &
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
        status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(iddpth),   &
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
        status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Vid(idDano),   &
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
          status=def_var(ng, iNLM, FISH(ng)%ncid, FISH(ng)%Tid(itrc),   &
     &                   NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
          IF (exit_flag.ne.NoError) RETURN
        END DO
# endif
!
!  Define fish wet weight
!
        Vinfo( 1)='Fwwt'
        Vinfo( 2)='fish wet weight'
        Vinfo( 3)='g'
        Vinfo(14)='fwwt, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFwwt),    &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish worth
!
        Vinfo( 1)='Fworth'
        Vinfo( 2)='fish worth'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='fworth, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFworth),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish length
!
        Vinfo( 1)='Flength'
        Vinfo( 2)='fish length'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='flength, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFlngth),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!
!  Define fish P-value
!
        Vinfo( 1)='Fpval'
        Vinfo( 2)='fish p-value'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='fpval, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFpval),   &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish age
!
        Vinfo( 1)='Fage'
        Vinfo( 2)='fish age'
        Vinfo( 3)='year'
        Vinfo(14)='fage, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFage),    &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish birthday
!
        Vinfo( 1)='Fbday'
        Vinfo( 2)='fish birthday'
        Vinfo( 3)='year'
        Vinfo(14)='fbday, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFbday),   &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define PS consumption
!
        Vinfo( 1)='FcsmPS'
        Vinfo( 2)='PS consumption'
        Vinfo( 3)='gPrey / gFish / hour'
        Vinfo(14)='fcsmps, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFcsmPS),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define PL consumption
!
        Vinfo( 1)='FcsmPL'
        Vinfo( 2)='PL consumption'
        Vinfo( 3)='gPrey / gFish / hour'
        Vinfo(14)='fcsmpl, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFcsmPL),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ZS consumption
!
        Vinfo( 1)='FcsmZS'
        Vinfo( 2)='ZS consumption'
        Vinfo( 3)='gPrey / gFish / hour'
        Vinfo(14)='fcsmzs, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFcsmZS),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ZL consumption
!
        Vinfo( 1)='FcsmZL'
        Vinfo( 2)='ZL consumption'
        Vinfo( 3)='gPrey / gFish / hour'
        Vinfo(14)='fcsmzl, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFcsmZL),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define ZP consumption
!
        Vinfo( 1)='FcsmZP'
        Vinfo( 2)='ZP consumption'
        Vinfo( 3)='gPrey / gFish / hour'
        Vinfo(14)='fcsmzp, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFcsmZP),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define Fish consumption on species 1
!
!        Vinfo( 1)='FcsmF1'
!        Vinfo( 2)='fish csmp on sp1'
!        Vinfo( 3)='gPrey / gFish / hour'
!        Vinfo(14)='fcsmf1, scalar, series'
!        Vinfo(16)=Vname(1,idtime)
!# ifndef NO_4BYTE_REALS
!        Vinfo(24)='_FillValue'
!        Aval(6)=spval
!# endif
!        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFcsmF1),  &
!     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
!        IF (exit_flag.ne.NoError) RETURN
!
!  Define Fish consumption on species 2
!
!        Vinfo( 1)='FcsmF2'
!        Vinfo( 2)='fish csmp on sp2'
!        Vinfo( 3)='gPrey / gFish / hour'
!        Vinfo(14)='fcsmf2, scalar, series'
!        Vinfo(16)=Vname(1,idtime)
!# ifndef NO_4BYTE_REALS
!        Vinfo(24)='_FillValue'
!        Aval(6)=spval
!# endif
!        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFcsmF2),  &
!     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
!        IF (exit_flag.ne.NoError) RETURN
!
!  Define Fish consumption on species 3
!
!        Vinfo( 1)='FcsmF3'
!        Vinfo( 2)='fish csmp on sp3'
!        Vinfo( 3)='gPrey / gFish / hour'
!        Vinfo(14)='fcsmf3, scalar, series'
!        Vinfo(16)=Vname(1,idtime)
!# ifndef NO_4BYTE_REALS
!        Vinfo(24)='_FillValue'
!        Aval(6)=spval
!# endif
!        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFcsmF3),  &
!     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
!        IF (exit_flag.ne.NoError) RETURN
!
!  Define Fish consumption on species 4
!
!        Vinfo( 1)='FcsmF4'
!        Vinfo( 2)='fish csmp on sp4'
!        Vinfo( 3)='gPrey / gFish / hour'
!        Vinfo(14)='fcsmf4, scalar, series'
!        Vinfo(16)=Vname(1,idtime)
!# ifndef NO_4BYTE_REALS
!        Vinfo(24)='_FillValue'
!        Aval(6)=spval
!# endif
!        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFcsmF4),  &
!     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
!        IF (exit_flag.ne.NoError) RETURN
!
!  Define Fish consumption on species 5
!
!        Vinfo( 1)='FcsmF5'
!        Vinfo( 2)='fish csmp on sp5'
!        Vinfo( 3)='gPrey / gFish / hour'
!        Vinfo(14)='fcsmf5, scalar, series'
!        Vinfo(16)=Vname(1,idtime)
!# ifndef NO_4BYTE_REALS
!        Vinfo(24)='_FillValue'
!        Aval(6)=spval
!# endif
!        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFcsmF5),  &
!     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
!        IF (exit_flag.ne.NoError) RETURN
!
!  Define respiration
!
        Vinfo( 1)='Fresp'
        Vinfo( 2)='Respiration'
        Vinfo( 3)='gPrey / gFish / hour'
        Vinfo(14)='fresp, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFresp),   &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define egg duration
!
        Vinfo( 1)='EggDur'
        Vinfo( 2)='egg duration'
        Vinfo( 3)='days'
        Vinfo(14)='eggdur, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFdegg),   &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define numbers of eggs
!
        Vinfo( 1)='EggNum'
        Vinfo( 2)='number of eggs'
        Vinfo( 3)='individual'
        Vinfo(14)='eggnum, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFnegg),   &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define yolk sac duration
!
        Vinfo( 1)='YsacDur'
        Vinfo( 2)='yolk sac duration'
        Vinfo( 3)='days'
        Vinfo(14)='ysacdur, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFdysac),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define numbers of yolk sacs
!
        Vinfo( 1)='YsacNum'
        Vinfo( 2)='number of yolk sacs'
        Vinfo( 3)='individual'
        Vinfo(14)='ysacnum, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFnysac),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define larval duration
!
        Vinfo( 1)='LarvDur'
        Vinfo( 2)='larval duration'
        Vinfo( 3)='days'
        Vinfo(14)='larvdur, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFdlarv),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define numbers of larvae
!
        Vinfo( 1)='LarvNum'
        Vinfo( 2)='number of larvae'
        Vinfo( 3)='individual'
        Vinfo(14)='larvnum, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFnlarv),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define juvenile duration
!
        Vinfo( 1)='JuvDur'
        Vinfo( 2)='juvenile duration'
        Vinfo( 3)='days'
        Vinfo(14)='juvdur, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFdjuv),   &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define numbers of juveniles
!
        Vinfo( 1)='JuvNum'
        Vinfo( 2)='number of juveniles'
        Vinfo( 3)='individual'
        Vinfo(14)='juvnum, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFnjuv),   &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define numbers of subadults
!
        Vinfo( 1)='SubaNum'
        Vinfo( 2)='number of subadults'
        Vinfo( 3)='individual'
        Vinfo(14)='subanum, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFnsuba),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define natural mortality
!
        Vinfo( 1)='FmortN'
        Vinfo( 2)='natural mortality'
        Vinfo( 3)='day-1'
        Vinfo(14)='fmortn, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFmortN),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define starvation mortality
!
        Vinfo( 1)='FmortS'
        Vinfo( 2)='starvation mortality'
        Vinfo( 3)='day-1'
        Vinfo(14)='fmorts, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFmortS),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fishing mortality
!
        Vinfo( 1)='FmortF'
        Vinfo( 2)='fishing mortality'
        Vinfo( 3)='day-1'
        Vinfo(14)='fmortf, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFmortF),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define predation mortality
!
        Vinfo( 1)='FmortP'
        Vinfo( 2)='predation mortality'
        Vinfo( 3)='day-1'
        Vinfo(14)='fmortp, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFmortP),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define summed predation mortality
!
        Vinfo( 1)='FmortPsum'
        Vinfo( 2)='summed predation mortality'
        Vinfo( 3)='day-1'
        Vinfo(14)='fmrtps, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFmrtPs),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish lifestage
!
        Vinfo( 1)='Lifestage'
        Vinfo( 2)='fish lifestage'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='lifestage, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=0
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFlfsg),   &
     &                 nf90_int, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish cause of death (1=starve, 2=Nmort, 3=Pmort, 4=Fmort)
!
        Vinfo( 1)='Deathby'
        Vinfo( 2)='fish cause of death'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='deathby, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=0
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFdeath),  &
     &                 nf90_int, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish batch count
!
        Vinfo( 1)='Fbatch'
        Vinfo( 2)='fish batch count'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='fbatch, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFbatch),  &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Define fish egg count
!
        Vinfo( 1)='Feggs'
        Vinfo( 2)='fish egg count'
        Vinfo( 3)='nondimensional'
        Vinfo(14)='feggs, scalar, series'
        Vinfo(16)=Vname(1,idtime)
# ifndef NO_4BYTE_REALS
        Vinfo(24)='_FillValue'
        Aval(6)=spval
# endif
        status=def_var(ng, iNLM, FISH(ng)%ncid,FISH(ng)%Vid(idFeggs),   &
     &                 NF_FOUT, 2, fgrd, Aval, Vinfo, ncname)
        IF (exit_flag.ne.NoError) RETURN
!
!  Initialize unlimited time record dimension.
!
        FISH(ng)%Rindex=0
!
!-----------------------------------------------------------------------
!  Leave definition mode.
!-----------------------------------------------------------------------
!
        CALL netcdf_enddef (ng, iNLM, ncname, FISH(ng)%ncid)
        IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Write out time-recordless, information variables.
!-----------------------------------------------------------------------
!
        CALL wrt_info (ng, iNLM, FISH(ng)%ncid, ncname)
        IF (exit_flag.ne.NoError) RETURN
      END IF DEFINE
!
!=======================================================================
!  Open an existing fish file, check its contents, and prepare for
!  appending data.
!=======================================================================
!
      QUERY : IF (.not.ldef) THEN
        ncname=FISH(ng)%name
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
        CALL netcdf_open (ng, iNLM, ncname, 1, FISH(ng)%ncid)
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
            FISH(ng)%Vid(idtime)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Fwwt') THEN
            got_var(idFwwt)=.TRUE.
            FISH(ng)%Vid(idFwwt)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Fworth') THEN
            got_var(idFworth)=.TRUE.
            FISH(ng)%Vid(idFworth)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Flength') THEN
            got_var(idFlngth)=.TRUE.
            FISH(ng)%Vid(idFlngth)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Fpval') THEN
            got_var(idFpval)=.TRUE.
            FISH(ng)%Vid(idFpval)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Fage') THEN
            got_var(idFage)=.TRUE.
            FISH(ng)%Vid(idFage)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Fbday') THEN
            got_var(idFbday)=.TRUE.
            FISH(ng)%Vid(idFbday)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Fyield') THEN
            got_var(idFyield)=.TRUE.
            FISH(ng)%Vid(idFyield)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FcsmPS') THEN
            got_var(idFcsmPS)=.TRUE.
            FISH(ng)%Vid(idFcsmPS)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FcsmPL') THEN
            got_var(idFcsmPL)=.TRUE.
            FISH(ng)%Vid(idFcsmPL)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FcsmZS') THEN
            got_var(idFcsmZS)=.TRUE.
            FISH(ng)%Vid(idFcsmZS)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FcsmZL') THEN
            got_var(idFcsmZL)=.TRUE.
            FISH(ng)%Vid(idFcsmZL)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FcsmZP') THEN
            got_var(idFcsmZP)=.TRUE.
            FISH(ng)%Vid(idFcsmZP)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FcsmF1') THEN
            got_var(idFcsmF1)=.TRUE.
            FISH(ng)%Vid(idFcsmF1)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FcsmF2') THEN
            got_var(idFcsmF2)=.TRUE.
            FISH(ng)%Vid(idFcsmF2)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FcsmF3') THEN
            got_var(idFcsmF3)=.TRUE.
            FISH(ng)%Vid(idFcsmF3)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FcsmF4') THEN
            got_var(idFcsmF4)=.TRUE.
            FISH(ng)%Vid(idFcsmF4)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FcsmF5') THEN
            got_var(idFcsmF5)=.TRUE.
            FISH(ng)%Vid(idFcsmF5)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Fresp') THEN
            got_var(idFresp)=.TRUE.
            FISH(ng)%Vid(idFresp)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'EggDur') THEN
            got_var(idFdegg)=.TRUE.
            FISH(ng)%Vid(idFdegg)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'EggNum') THEN
            got_var(idFnegg)=.TRUE.
            FISH(ng)%Vid(idFnegg)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'YsacDur') THEN
            got_var(idFdysac)=.TRUE.
            FISH(ng)%Vid(idFdysac)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'YsacNum') THEN
            got_var(idFnysac)=.TRUE.
            FISH(ng)%Vid(idFnysac)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'LarvDur') THEN
            got_var(idFdlarv)=.TRUE.
            FISH(ng)%Vid(idFdlarv)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'LarvNum') THEN
            got_var(idFnlarv)=.TRUE.
            FISH(ng)%Vid(idFnlarv)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'JuvDur') THEN
            got_var(idFdjuv)=.TRUE.
            FISH(ng)%Vid(idFdjuv)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'JuvNum') THEN
            got_var(idFnjuv)=.TRUE.
            FISH(ng)%Vid(idFnjuv)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'SubaNum') THEN
            got_var(idFnsuba)=.TRUE.
            FISH(ng)%Vid(idFnsuba)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FmortN') THEN
            got_var(idFmortN)=.TRUE.
            FISH(ng)%Vid(idFmortN)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FmortS') THEN
            got_var(idFmortS)=.TRUE.
            FISH(ng)%Vid(idFmortS)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FmortF') THEN
            got_var(idFmortF)=.TRUE.
            FISH(ng)%Vid(idFmortF)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FmortP') THEN
            got_var(idFmortP)=.TRUE.
            FISH(ng)%Vid(idFmortP)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'FmortPsum') THEN
            got_var(idFmrtPs)=.TRUE.
            FISH(ng)%Vid(idFmrtPs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Fbatch') THEN
            got_var(idFbatch)=.TRUE.
            FISH(ng)%Vid(idFbatch)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Feggs') THEN
            got_var(idFeggs)=.TRUE.
            FISH(ng)%Vid(idFeggs)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Species') THEN
            got_var(idFspc)=.TRUE.
            FISH(ng)%Vid(idFspc)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Lifestage') THEN
            got_var(idFlfsg)=.TRUE.
            FISH(ng)%Vid(idFlfsg)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Deathby') THEN
            got_var(idFdeath)=.TRUE.
            FISH(ng)%Vid(idFdeath)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Hbehave') THEN
            got_var(idFhbhv)=.TRUE.
            FISH(ng)%Vid(idFhbhv)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Vbehave') THEN
            got_var(idFvbhv)=.TRUE.
            FISH(ng)%Vid(idFvbhv)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Alive') THEN
            got_var(idFalv)=.TRUE.
            FISH(ng)%Vid(idFalv)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Tspawn') THEN
            got_var(idFTspwn)=.TRUE.
            FISH(ng)%Vid(idFTspwn)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Iniworth') THEN
            got_var(idFiniw)=.TRUE.
            FISH(ng)%Vid(idFiniw)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Xgrid') THEN
            got_var(idXgrd)=.TRUE.
            FISH(ng)%Vid(idXgrd)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.'Ygrid') THEN
            got_var(idYgrd)=.TRUE.
            FISH(ng)%Vid(idYgrd)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.'Zgrid') THEN
            got_var(idZgrd)=.TRUE.
            FISH(ng)%Vid(idZgrd)=var_id(i)
# endif
          ELSE IF (spherical.and.TRIM(var_name(i)).eq.'lon') THEN
            got_var(idglon)=.TRUE.
            FISH(ng)%Vid(idglon)=var_id(i)
          ELSE IF (spherical.and.TRIM(var_name(i)).eq.'lat') THEN
            got_var(idglat)=.TRUE.
            FISH(ng)%Vid(idglat)=var_id(i)
          ELSE IF (.not.spherical.and.TRIM(var_name(i)).eq.'x') THEN
            got_var(idglon)=.TRUE.
            FISH(ng)%Vid(idglon)=var_id(i)
          ELSE IF (.not.spherical.and.TRIM(var_name(i)).eq.'y') THEN
            got_var(idglat)=.TRUE.
            FISH(ng)%Vid(idglat)=var_id(i)
# ifdef SOLVE3D
          ELSE IF (TRIM(var_name(i)).eq.'depth') THEN
            got_var(iddpth)=.TRUE.
            FISH(ng)%Vid(iddpth)=var_id(i)
          ELSE IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idDano))) THEN
            got_var(idDano)=.TRUE.
            FISH(ng)%Vid(idDano)=var_id(i)
# endif
          END IF
# ifdef SOLVE3D
          DO itrc=1,NT(ng)
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              got_var(idTvar(itrc))=.TRUE.
              FISH(ng)%Tid(itrc)=var_id(i)
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
        IF (.not.got_var(idFwwt)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFwwt)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFworth)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFworth)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFlngth)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFlngth)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFage)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFage)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFbday)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFbday)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFspc)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFspc)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFlfsg)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFlfsg)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFdeath)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFdeath)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFdeath)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFdeath)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFhbhv)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFhbhv)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFvbhv)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFvbhv)),         &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFalv)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFalv)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFTspwn)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFTspwn)),        &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        IF (.not.got_var(idFiniw)) THEN
          IF (Master) WRITE (stdout,60) TRIM(Vname(1,idFiniw)),         &
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
!  Initialize fish positions to the appropriate values.
!-----------------------------------------------------------------------
!
!  Set-up fish time record.
!
        IF (fishrrec(ng).lt.0) THEN
          FISH(ng)%Rindex=rec_size
        ELSE
          FISH(ng)%Rindex=ABS(fishrrec(ng))
        END IF
        FISH(ng)%Nrec=FISH(ng)%Rindex
!
!  Read in fish nondimentional horizontal positions.
!
        CALL netcdf_get_fvar (ng, iNLM, ncname, 'Xgrid',                &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
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
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
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
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
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
        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFwwt)),                    &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))
        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifwwt,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFworth)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))
        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifworth,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFlngth)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))
        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(iflngth,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFage)),                    &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))
        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifage,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFbday)),                   &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))
        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifbday,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFbatch)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))
        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifbatch,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFeggs)),                   &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))
        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifeggs,l)=Tinp(l)
        END DO

        CALL netcdf_get_ivar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFspc)),                    &
     &                        FISHES(ng)%species,                       &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        CALL netcdf_get_ivar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFlfsg)),                   &
     &                        FISHES(ng)%lifestage,                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        CALL netcdf_get_ivar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFdeath)),                  &
     &                        FISHES(ng)%deathby,                       &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        CALL netcdf_get_ivar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFhbhv)),                   &
     &                        TinpI,                                    &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))
        DO l=1,Nfish(ng)
          FISHES(ng)%swimtype(1,l)=TinpI(l)
        END DO
        CALL netcdf_get_ivar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFvbhv)),                   &
     &                        TinpI,                                    &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))
        DO l=1,Nfish(ng)
          FISHES(ng)%swimtype(2,l)=TinpI(l)
        END DO

        CALL netcdf_get_lvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFalv)),                    &
     &                        FISHES(ng)%alive,                         &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFTspwn)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1/),                            &
     &                        total = (/Nfish(ng)/))
        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(iftspwn,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFiniw)),                   &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1/),                            &
     &                        total = (/Nfish(ng)/))
        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifiniwth,l)=Tinp(l)
        END DO

     !-----------------------------------------------------------------
     ! RD : this whole section was missing :
        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFdegg)),                   &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%egg_dur(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFnegg)),                   &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%egg_num(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFnysac)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%ysac_num(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFcsmPL)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifcsmPL,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFcsmPS)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifcsmPS,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFcsmZL)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifcsmZL,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFcsmZP)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifcsmZP,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFcsmZS)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifcsmZS,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFmortF)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%fmortF(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFmortN)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%fmortS(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFmortP)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%fmortP(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFmortS)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%fmortS(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFmrtPs)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%fmortPsum(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFpval)),                   &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifpval,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFresp)),                   &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%bioenergy(ifresp,l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFdjuv)),                   &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%juv_dur(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFnjuv)),                   &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%juv_num(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFdlarv)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%larv_dur(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFnlarv)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%larv_num(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFnsuba)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%suba_num(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFdysac)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%ysac_dur(l)=Tinp(l)
        END DO

        CALL netcdf_get_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(1,idFnysac)),                  &
     &                        Tinp,                                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1,FISH(ng)%Rindex/),            &
     &                        total = (/Nfish(ng),1/))

        DO l=1,Nfish(ng)
          FISHES(ng)%ysac_num(l)=Tinp(l)
        END DO


     !-----------------------------------------------------------------

        CALL netcdf_get_ivar (ng, iNLM, ncname,                         &
     &                        'Fish_age_base',                          &
     &                        FISHES(ng)%age_base,                      &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1/),                            &
     &                        total = (/Nspecies(ng)/))

        CALL netcdf_get_ivar (ng, iNLM, ncname,                         &
     &                        'Fish_spec_base',                         &
     &                        FISHES(ng)%species_base,                  &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1/),                            &
     &                        total = (/Nspecies(ng)/))

        CALL netcdf_get_ivar (ng, iNLM, ncname,                         &
     &                        'Fish_num_free',                          &
     &                        FISHES(ng)%num_free,                      &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1/),                            &
     &                        total = (/Nspecies(ng)/))

        CALL netcdf_get_ivar (ng, iNLM, ncname,                         &
     &                        'Fish_next_free',                         &
     &                        FISHES(ng)%next_free,                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1/),                            &
     &                        total = (/Nspecies(ng)/))

        CALL netcdf_get_ivar (ng, iNLM, ncname,                         &
     &                        'Fish_num_super',                         &
     &                        FISHES(ng)%num_super,                     &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        start = (/1/),                            &
     &                        total = (/Nspecies(ng)/))

        IF (exit_flag.ne.NoError) RETURN

        IF (Master) THEN
          print *, 'First Year for Fish'
          print *, 'Free fish ', FISHES(ng) % num_free
          print *, 'Fish ages ', FISHES(ng) % age_base
        END IF
      END IF QUERY
!
  10  FORMAT (6x,'DEF_FISH  - creating fish file: ',a)
  20  FORMAT (6x,'DEF_FISH  - inquiring fish file: ',a)
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
