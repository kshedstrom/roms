#include "cppdefs.h"
#ifdef NEMURO_SAN
      SUBROUTINE wrt_fish (ng)
!
!svn $Id: wrt_fish.F 975 2009-05-05 22:51:13Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
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
     &                      TRIM(Vname(1,idfwwt)), Tout,                &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idfwwt,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write fish batch number.
      DO l=1,Nfish(ng)
        IF (FISHES(ng)%bounded(l)) THEN
          Tout(l)=FISHES(ng)%bioenergy(ifbatch,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FISHname(ng),                     &
     &                      TRIM(Vname(1,idbatch)), Tout,               &
     &                      (/1,tFISHindx(ng)/), (/Nfish(ng),1/),       &
     &                      ncid = ncFISHid(ng),                        &
     &                      varid = fishVid(idbatch,ng))
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
