#include "cppdefs.h"
#ifdef NEMURO_SAN
      SUBROUTINE wrt_fish (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
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
      USE mod_biology
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

      integer, dimension(Nfish(ng)) :: ToutI
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
      FISH(ng)%Rindex=FISH(ng)%Rindex+1
      FISH(ng)%Nrec=FISH(ng)%Nrec+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/FISH(ng)%Rindex/), (/1/),                 &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish species.
!
      CALL netcdf_put_ivar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFspc)),                      &
     &                      FISHES(ng) % species,                       &
     &                      (/1/), (/Nfish(ng)/),                       &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFspc))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish lifestage.
!
      CALL netcdf_put_ivar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFlfsg)),                     &
     &                      FISHES(ng) % lifestage,                     &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFlfsg))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish cause of death
!
      CALL netcdf_put_ivar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFdeath)),                    &
     &                      FISHES(ng) % deathby,                       &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFdeath))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish behavior.
!
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          ToutI(l)=FISHES(ng)%swimtype(1,l)
        ELSE
          ToutI(l)=0
        END IF
      END DO
      CALL netcdf_put_ivar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFhbhv)),                     &
     &                      ToutI,                                      &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFhbhv))
      IF (exit_flag.ne.NoError) RETURN
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          ToutI(l)=FISHES(ng)%swimtype(2,l)
        ELSE
          ToutI(l)=0
        END IF
      END DO
      CALL netcdf_put_ivar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFvbhv)),                     &
     &                      ToutI,                                      &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFvbhv))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish aliveness.
!
      CALL netcdf_put_lvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFalv)),                      &
     &                      FISHES(ng) % alive,                         &
     &                      (/1/), (/Nfish(ng)/),                       &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFalv))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out time since last spawn.
!
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(iftspwn,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFTspwn)),                    &
     &                      Tout, (/1/), (/Nfish(ng)/),                 &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFTspwn))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out initial worth.
!
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifiniwth,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFiniw)),                     &
     &                      Tout, (/1/), (/Nfish(ng)/),                 &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFiniw))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish age base.
!
      CALL netcdf_put_ivar (ng, iNLM, FISH(ng)%name,                    &
     &                      'Fish_age_base',                            &
     &                      FISHES(ng) % age_base,                      &
     &                      (/1/), (/Nspecies(ng)/),                    &
     &                      ncid = FISH(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish species base.
!
      CALL netcdf_put_ivar (ng, iNLM, FISH(ng)%name,                    &
     &                      'Fish_spec_base',                           &
     &                      FISHES(ng) % species_base,                  &
     &                      (/1/), (/Nspecies(ng)/),                    &
     &                      ncid = FISH(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish num free.
!
      CALL netcdf_put_ivar (ng, iNLM, FISH(ng)%name,                    &
     &                      'Fish_num_free',                            &
     &                      FISHES(ng) % num_free,                      &
     &                      (/1/), (/Nspecies(ng)/),                    &
     &                      ncid = FISH(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fish next free.
!
      CALL netcdf_put_ivar (ng, iNLM, FISH(ng)%name,                    &
     &                      'Fish_next_free',                           &
     &                      FISHES(ng) % next_free,                     &
     &                      (/1/), (/Nspecies(ng)/),                    &
     &                      ncid = FISH(ng)%ncid)
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out superindividuals per spawn.
!
      CALL netcdf_put_ivar (ng, iNLM, FISH(ng)%name,                    &
     &                      'Fish_num_super',                           &
     &                      FISHES(ng) % num_super,                     &
     &                      (/1/), (/Nspecies(ng)/),                    &
     &                      ncid = FISH(ng)%ncid)
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
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      'Xgrid', Tout,                              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idXgrd))
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
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      'Ygrid', Tout,                              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idYgrd))
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
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      'Zgrid', Tout,                              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idZgrd))
      IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Write out fish (lon,lat) or (x,y) locations.
!
      DO l=1,Nfish(ng)
        Tout(l)=FISHES(ng)%track(iflon,nf(ng),l)
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                  &
     &                        'lon', Tout,                              &
     &                        (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),   &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        varid = FISH(ng)%Vid(idglon))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                  &
     &                        'x', Tout,                                &
     &                        (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),   &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        varid = FISH(ng)%Vid(idglon))
      END IF
      IF (exit_flag.ne.NoError) RETURN
!
      DO l=1,Nfish(ng)
        Tout(l)=FISHES(ng)%track(iflat,nf(ng),l)
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                  &
     &                        'lat', Tout,                              &
     &                        (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),   &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        varid = FISH(ng)%Vid(idglat))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                  &
     &                        'y', Tout,                                &
     &                        (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),   &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        varid = FISH(ng)%Vid(idglat))
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
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      'depth', Tout,                              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(iddpth))
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
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idDano)), Tout,                &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idDano))
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
        CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                  &
     &                        TRIM(Vname(1,idTvar(itrc))), Tout,        &
     &                        (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),   &
     &                        ncid = FISH(ng)%ncid,                     &
     &                        varid = FISH(ng)%Tid(itrc))
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
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFwwt)), Tout,                &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFwwt))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish worth.
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifworth,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFworth)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFworth))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish length
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(iflngth,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFlngth)),                    &
     &                      Tout,                                       &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFlngth))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish P-value.
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifpval,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFpval)), Tout,               &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFpval))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish age
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifage,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFage)), Tout,                &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFage))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish birthday
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifbday,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFbday)), Tout,               &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFbday))
      IF (exit_flag.ne.NoError) RETURN
!  Write PS consumption
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmPS,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFcsmPS)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFcsmPS))
      IF (exit_flag.ne.NoError) RETURN
!  Write PL consumption
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmPL,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFcsmPL)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFcsmPL))
      IF (exit_flag.ne.NoError) RETURN
!  Write ZS consumption
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmZS,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFcsmZS)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFcsmZS))
      IF (exit_flag.ne.NoError) RETURN
!  Write ZL consumption
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmZL,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFcsmZL)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFcsmZL))
      IF (exit_flag.ne.NoError) RETURN
!  Write ZP consumption
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifcsmZP,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFcsmZP)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFcsmZP))
      IF (exit_flag.ne.NoError) RETURN
!  Write Fish consumption on species 1
!      DO l=1,Nfish(ng)
!        IF (FISHES(ng)%bounded(l)) THEN
!          Tout(l)=FISHES(ng)%bioenergy(ifcsmF1,l)
!        ELSE
!          Tout(l)=spval
!        END IF
!      END DO
!      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
!     &                      TRIM(Vname(1,idFcsmF1)), Tout,              &
!     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
!     &                      ncid = FISH(ng)%ncid,                       &
!     &                      varid = FISH(ng)%Vid(idFcsmF1))
!      IF (exit_flag.ne.NoError) RETURN
!  Write Fish consumption on species 2
!      DO l=1,Nfish(ng)
!        IF (FISHES(ng)%bounded(l)) THEN
!          Tout(l)=FISHES(ng)%bioenergy(ifcsmF2,l)
!        ELSE
!          Tout(l)=spval
!        END IF
!      END DO
!      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
!     &                      TRIM(Vname(1,idFcsmF2)), Tout,              &
!     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
!     &                      ncid = FISH(ng)%ncid,                       &
!     &                      varid = FISH(ng)%Vid(idFcsmF2))
!      IF (exit_flag.ne.NoError) RETURN
!  Write Fish consumption on species 3
!      DO l=1,Nfish(ng)
!        IF (FISHES(ng)%bounded(l)) THEN
!          Tout(l)=FISHES(ng)%bioenergy(ifcsmF3,l)
!        ELSE
!          Tout(l)=spval
!        END IF
!      END DO
!      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
!     &                      TRIM(Vname(1,idFcsmF3)), Tout,              &
!     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
!     &                      ncid = FISH(ng)%ncid,                       &
!     &                      varid = FISH(ng)%Vid(idFcsmF3))
!      IF (exit_flag.ne.NoError) RETURN
!  Write Fish consumption on species 4
!      DO l=1,Nfish(ng)
!        IF (FISHES(ng)%bounded(l)) THEN
!          Tout(l)=FISHES(ng)%bioenergy(ifcsmF4,l)
!        ELSE
!          Tout(l)=spval
!        END IF
!      END DO
!      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
!     &                      TRIM(Vname(1,idFcsmF4)), Tout,              &
!     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
!     &                      ncid = FISH(ng)%ncid,                       &
!     &                      varid = FISH(ng)%Vid(idFcsmF4))
!      IF (exit_flag.ne.NoError) RETURN
!  Write Fish consumption on species 5
!      DO l=1,Nfish(ng)
!        IF (FISHES(ng)%bounded(l)) THEN
!          Tout(l)=FISHES(ng)%bioenergy(ifcsmF5,l)
!        ELSE
!          Tout(l)=spval
!        END IF
!      END DO
!      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
!     &                      TRIM(Vname(1,idFcsmF5)), Tout,              &
!     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
!     &                      ncid = FISH(ng)%ncid,                       &
!     &                      varid = FISH(ng)%Vid(idFcsmF5))
!      IF (exit_flag.ne.NoError) RETURN
!  Write respiration
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifresp,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFresp)), Tout,               &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFresp))
      IF (exit_flag.ne.NoError) RETURN
!  Write egg duration
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%egg_dur(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFdegg)), Tout,               &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFdegg))
      IF (exit_flag.ne.NoError) RETURN
!  Write number of eggs
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%egg_num(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFnegg)), Tout,               &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFnegg))
      IF (exit_flag.ne.NoError) RETURN
!  Write yolk sac duration
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%ysac_dur(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFdysac)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFdysac))
      IF (exit_flag.ne.NoError) RETURN
!  Write number of yolk sacs
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%ysac_num(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFnysac)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFnysac))
      IF (exit_flag.ne.NoError) RETURN
!  Write larval duration
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%larv_dur(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFdlarv)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFdlarv))
      IF (exit_flag.ne.NoError) RETURN
!  Write number of larvae
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%larv_num(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFnlarv)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFnlarv))
      IF (exit_flag.ne.NoError) RETURN
!  Write juvenile duration
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%juv_dur(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFdjuv)), Tout,               &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFdjuv))
      IF (exit_flag.ne.NoError) RETURN
!  Write number of juveniles
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%juv_num(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFnjuv)), Tout,               &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFnjuv))
      IF (exit_flag.ne.NoError) RETURN
!  Write number of subadults
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%suba_num(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFnsuba)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFnsuba))
      IF (exit_flag.ne.NoError) RETURN
!  Write natural mortality
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%fmortN(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFmortN)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFmortN))
      IF (exit_flag.ne.NoError) RETURN
!  Write starvation mortality
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%fmortS(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFmortS)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFmortS))
      IF (exit_flag.ne.NoError) RETURN
!  Write fishing mortality
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%fmortF(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFmortF)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFmortF))
      IF (exit_flag.ne.NoError) RETURN
!  Write predation mortality
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%fmortP(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFmortP)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFmortP))
      IF (exit_flag.ne.NoError) RETURN
!  Write summed predation mortality
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%fmortPsum(l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFmrtPs)),                    &
     &                      Tout,                                       &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFmrtPs))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish batch.
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifbatch,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFbatch)), Tout,              &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFbatch))
      IF (exit_flag.ne.NoError) RETURN
!  Write total fish eggs.
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifteggs,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISH(ng)%name,                    &
     &                      TRIM(Vname(1,idFeggs)),                     &
     &                      Tout,                                       &
     &                      (/1,FISH(ng)%Rindex/), (/Nfish(ng),1/),     &
     &                      ncid = FISH(ng)%ncid,                       &
     &                      varid = FISH(ng)%Vid(idFeggs))
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize fish NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, FISH(ng)%name, FISH(ng)%ncid)
!
#else
      SUBROUTINE wrt_fish
#endif
      RETURN
      END SUBROUTINE wrt_fish
