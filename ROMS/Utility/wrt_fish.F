#include "cppdefs.h"
#ifdef NEMURO_SAN
      SUBROUTINE wrt_fish (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes simulated drifter trajectories into fish     !
!  NetCDF file.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fish
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: itrc, l, status

      real(r8), dimension(Nfish(ng)) :: Tout
!
      SourceFile='wrt_fish.F'
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set time record index.
!
      tFISHindx(ng)=tFISHindx(ng)+1
      NrecFISH(ng)=NrecFISH(ng)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/tFISHindx(ng)/), (/1/),                   &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish species.
!
      CALL netcdf_put_ivar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFspc)),                      &
     &                      FISHES(ng) % species,                       &
     &                      (/1/), (/Nfish(ng)/),                       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFspc,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish lifestage.
!
!      CALL netcdf_put_ivar (ng, iNLM, FISHname(ng),                     &
!     &                      TRIM(Vname(1,idFlfsg)),                     &
!     &                      FISHES(ng) % lifestage,                     &
!     &                      (/1/), (/Nfish(ng)/),                       &
!     &                      ncid = ncFISHid(ng),                        &
!     &                      varid = fishVid(idFlfsg,ng))
!      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish swimtype.
!
      CALL netcdf_put_ivar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFswmt)),                     &
     &                      FISHES(ng) % swimtype,                      &
     &                      (/1/), (/Nfish(ng)/),                       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFswmt,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish aliveness.
!
      CALL netcdf_put_lvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFalv)),                      &
     &                      FISHES(ng) % alive,                         &
     &                      (/1/), (/Nfish(ng)/),                       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFalv,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish eat fish logic.
!
      CALL netcdf_put_lvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFeatf)),                     &
     &                      FISHES(ng) % eatfish,                       &
     &                      (/1/), (/Nfish(ng)/),                       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFeatf,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish age base.
!
      CALL netcdf_put_ivar (ng, iNLM, FISHname(ng),                     &
     &                      'Fish_age_base',                            &
     &                      FISHES(ng) % age_base,                      &
     &                      (/1/), (/Nspecies(ng)/),                    &
     &                      ncid = ncFISHid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish species base.
!
      CALL netcdf_put_ivar (ng, iNLM, FISHname(ng),                     &
     &                      'Fish_spec_base',                           &
     &                      FISHES(ng) % species_base,                  &
     &                      (/1/), (/Nspecies(ng)/),                    &
     &                      ncid = ncFISHid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish num free.
!
      CALL netcdf_put_ivar (ng, iNLM, FISHname(ng),                     &
     &                      'Fish_num_free',                            &
     &                      FISHES(ng) % num_free,                      &
     &                      (/1/), (/Nspecies(ng)/),                    &
     &                      ncid = ncFISHid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish next free.
!
      CALL netcdf_put_ivar (ng, iNLM, FISHname(ng),                     &
     &                      'Fish_next_free',                           &
     &                      FISHES(ng) % next_free,                     &
     &                      (/1/), (/Nspecies(ng)/),                    &
     &                      ncid = ncFISHid(ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish X-grid locations.
!
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%track(ixgrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      'Xgrid', Tout,                              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idXgrd,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish Y-grid locations.
!
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%track(iygrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      'Ygrid', Tout,                              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idYgrd,ng))
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Write out fish Z-grid locations.
!
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%track(izgrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      'Zgrid', Tout,                              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idZgrd,ng))
      IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Write out fish (lon,lat) or (x,y) locations.
!
      DO l=1,Nfish(ng)
        Tout(l)=FISHES(ng)%track(iflon,nf(ng),l)
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                   &
     &                        'lon', Tout,                              &
     &                        (/1,tFISHindx(ng)/), (/Nfish(ng),1/),     &
     &                        ncid = ncFISHid(ng),                      &
     &                        varid = fishVid(idglon,ng))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                   &
     &                        'x', Tout,                                &
     &                        (/1,tFISHindx(ng)/), (/Nfish(ng),1/),     &
     &                        ncid = ncFISHid(ng),                      &
     &                        varid = fishVid(idglon,ng))
      END IF
      IF (exit_flag.ne.NoError) RETURN
!
      DO l=1,Nfish(ng)
        Tout(l)=FISHES(ng)%track(iflat,nf(ng),l)
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                   &
     &                        'lat', Tout,                              &
     &                        (/1,tFISHindx(ng)/), (/Nfish(ng),1/),     &
     &                        ncid = ncFISHid(ng),                      &
     &                        varid = fishVid(idglat,ng))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                   &
     &                        'y', Tout,                                &
     &                        (/1,tFISHindx(ng)/), (/Nfish(ng),1/),     &
     &                        ncid = ncFISHid(ng),                      &
     &                        varid = fishVid(idglat,ng))
      END IF
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Write out fish depths.
!
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%track(idpth,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      'depth', Tout,                              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(iddpth,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out density anomaly.
!
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%track(ifden,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idDano)), Tout,                &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idDano,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        DO l=1,Nfish(ng)
          IF (FISHES(ng)%bounded(l)) THEN
            Tout(l)=FISHES(ng)%track(itrc+10,nf(ng),l)
          ELSE
            Tout(l)=spval
          END IF
        END DO
        CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                   &
     &                        TRIM(Vname(1,idTvar(itrc))), Tout,        &
     &                        (/1,tFISHindx(ng)/), (/Nfish(ng),1/),     &
     &                        ncid = ncFISHid(ng),                      &
     &                        varid = fishTid(itrc,ng))
        IF (exit_flag.ne.NoError) RETURN
      END DO
#  ifdef NEMURO_SAN
!  Write fish wet weight.
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifwwt,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFwwt)), Tout,                &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFwwt,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish worth.
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifworth,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFworth)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFworth,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish P-value.
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifpval,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFpval)), Tout,               &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFpval,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish age
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifage,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFage)), Tout,                &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFage,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish birthday
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifbday,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFbday)), Tout,               &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFbday,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish yield
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifyield,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFyield)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFyield,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write PS consumption
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmPS,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFcsmPS)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFcsmPS,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write PL consumption
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmPL,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFcsmPL)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFcsmPL,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write ZS consumption
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmZS,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFcsmZS)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFcsmZS,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write ZL consumption
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmZL,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFcsmZL)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFcsmZL,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write ZP consumption
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmZP,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFcsmZP)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFcsmZP,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write Fish consumption on species 1
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmF1,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFcsmF1)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFcsmF1,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write Fish consumption on species 2
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmF2,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFcsmF2)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFcsmF2,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write Fish consumption on species 3
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmF3,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFcsmF3)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFcsmF3,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write Fish consumption on species 4
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmF4,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFcsmF4)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFcsmF4,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write Fish consumption on species 5
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmF5,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFcsmF5)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFcsmF5,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write respiration
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifresp,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFresp)), Tout,               &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFresp,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write egg duration
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%egg_dur(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFdegg)), Tout,               &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFdegg,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write number of eggs
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%egg_num(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFnegg)), Tout,               &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFnegg,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write yolk sac duration
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%ysac_dur(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFdysac)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFdysac,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write number of yolk sacs
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%ysac_num(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFnysac)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFnysac,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write larval duration
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%larv_dur(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFdlarv)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFdlarv,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write number of larvae
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%larv_num(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFnlarv)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFnlarv,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write juvenile duration
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%juv_dur(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFdjuv)), Tout,               &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFdjuv,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write number of juveniles
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%juv_num(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFnjuv)), Tout,               &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFnjuv,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write number of subadults
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%suba_num(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFnsuba)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFnsuba,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write natural mortality
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%fmortN(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFmortN)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFmortN,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write starvation mortality
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%fmortS(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFmortS)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFmortS,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write fishing mortality
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%fmortF(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFmortF)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFmortF,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish lifestage
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%lifestage(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFlfsg)), Tout,               &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFlfsg,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish batch.
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifbatch,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idFbatch)), Tout,              &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idFbatch,ng))
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize fish NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, FISHname(ng), ncFISHid(ng))
!
#else
      SUBROUTINE wrt_fish
#endif
      RETURN
      END SUBROUTINE wrt_fish
