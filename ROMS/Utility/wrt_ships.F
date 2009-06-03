#include "cppdefs.h"
#if defined NEMURO_SAN && defined FLEET
      SUBROUTINE wrt_ships (ng)
!
!svn $Id: wrt_ships.F 975 2009-05-05 22:51:13Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes simulated drifter trajectories into ships    !
!  NetCDF file.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ships
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

      real(r8), dimension(Nships(ng)) :: Tout
!
      SourceFile='wrt_ships.F'
!
!-----------------------------------------------------------------------
!  Write out ship data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set time record index.
!
      tSHIPSindx(ng)=tSHIPSindx(ng)+1
      NrecSHIPS(ng)=NrecSHIPS(ng)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, SHIPname(ng),                     &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/tSHIPSindx(ng)/), (/1/),                  &
     &                      ncid = ncSHIPSid(ng),                       &
     &                      varid = shipsVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out ships X-grid locations.
!
      DO l=1,Nships(ng)
        IF (SHIPS(ng)%bounded(l)) THEN
          Tout(l)=SHIPS(ng)%track(ixgrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, SHIPname(ng),                     &
     &                      'Xgrid', Tout,                              &
     &                      (/1,tSHIPSindx(ng)/), (/Nships(ng),1/),     &
     &                      ncid = ncSHIPSid(ng),                       &
     &                      varid = shipsVid(idXgrd,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out ships Y-grid locations.
!
      DO l=1,Nships(ng)
        IF (SHIPS(ng)%bounded(l)) THEN
          Tout(l)=SHIPS(ng)%track(iygrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, SHIPname(ng),                     &
     &                      'Ygrid', Tout,                              &
     &                      (/1,tSHIPSindx(ng)/), (/Nships(ng),1/),     &
     &                      ncid = ncSHIPSid(ng),                       &
     &                      varid = shipsVid(idYgrd,ng))
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Write out ships Z-grid locations.
!
      DO l=1,Nships(ng)
        IF (SHIPS(ng)%bounded(l)) THEN
          Tout(l)=SHIPS(ng)%track(izgrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, SHIPname(ng),                     &
     &                      'Zgrid', Tout,                              &
     &                      (/1,tSHIPSindx(ng)/), (/Nships(ng),1/),     &
     &                      ncid = ncSHIPSid(ng),                       &
     &                      varid = shipsVid(idZgrd,ng))
      IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Write out ships (lon,lat) or (x,y) locations.
!
      DO l=1,Nships(ng)
        Tout(l)=SHIPS(ng)%track(iflon,nf(ng),l)
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, SHIPname(ng),                   &
     &                        'lon', Tout,                              &
     &                        (/1,tSHIPSindx(ng)/), (/Nships(ng),1/),   &
     &                        ncid = ncSHIPSid(ng),                     &
     &                        varid = shipsVid(idglon,ng))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, SHIPname(ng),                   &
     &                        'x', Tout,                                &
     &                        (/1,tSHIPSindx(ng)/), (/Nships(ng),1/),   &
     &                        ncid = ncSHIPSid(ng),                     &
     &                        varid = shipsVid(idglon,ng))
      END IF
      IF (exit_flag.ne.NoError) RETURN
!
      DO l=1,Nships(ng)
        Tout(l)=SHIPS(ng)%track(iflat,nf(ng),l)
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, SHIPname(ng),                   &
     &                        'lat', Tout,                              &
     &                        (/1,tSHIPSindx(ng)/), (/Nships(ng),1/),   &
     &                        ncid = ncSHIPSid(ng),                     &
     &                        varid = shipsVid(idglat,ng))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, SHIPname(ng),                   &
     &                        'y', Tout,                                &
     &                        (/1,tSHIPSindx(ng)/), (/Nships(ng),1/),   &
     &                        ncid = ncSHIPSid(ng),                     &
     &                        varid = shipsVid(idglat,ng))
      END IF
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Write out ships depths.
!
      DO l=1,Nships(ng)
        IF (SHIPS(ng)%bounded(l)) THEN
          Tout(l)=SHIPS(ng)%track(idpth,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, SHIPname(ng),                     &
     &                      'depth', Tout,                              &
     &                      (/1,tSHIPSindx(ng)/), (/Nships(ng),1/),     &
     &                      ncid = ncSHIPSid(ng),                       &
     &                      varid = shipsVid(iddpth,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out density anomaly.
!
      DO l=1,Nships(ng)
        IF (SHIPS(ng)%bounded(l)) THEN
          Tout(l)=SHIPS(ng)%track(ifden,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, SHIPname(ng),                     &
     &                      TRIM(Vname(1,idDano)), Tout,                &
     &                      (/1,tSHIPSindx(ng)/), (/Nships(ng),1/),     &
     &                      ncid = ncSHIPSid(ng),                       &
     &                      varid = shipsVid(idDano,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        DO l=1,Nships(ng)
          IF (SHIPS(ng)%bounded(l)) THEN
            Tout(l)=SHIPS(ng)%track(itrc+10,nf(ng),l)
          ELSE
            Tout(l)=spval
          END IF
        END DO
        CALL netcdf_put_fvar (ng, iNLM, SHIPname(ng),                   &
     &                        TRIM(Vname(1,idTvar(itrc))), Tout,        &
     &                        (/1,tSHIPSindx(ng)/), (/Nships(ng),1/),   &
     &                        ncid = ncSHIPSid(ng),                     &
     &                        varid = shipsTid(itrc,ng))
        IF (exit_flag.ne.NoError) RETURN
      END DO
#  ifdef NEMURO_SAN
!  Write ships wet weight.
!
!      DO l=1,Nships(ng)
!        IF (SHIPS(ng)%bounded(l)) THEN
!          Tout(l)=SHIPS(ng)%bioenergy(ifwwt,l)
!        ELSE
!          Tout(l)=spval
!        END IF
!      END DO
!      CALL netcdf_put_fvar (ng, iNLM, SHIPname(ng),                     &
!     &                      TRIM(Vname(1,idfwwt)), Tout,                &
!     &                      (/1,tSHIPSindx(ng)/), (/Nships(ng),1/),     &
!     &                      ncid = ncSHIPSid(ng),                       &
!     &                      varid = shipsVid(idfwwt,ng))
!      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize ships NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, SHIPname(ng), ncSHIPSid(ng))
!
#else
      SUBROUTINE wrt_ships
#endif
      RETURN
      END SUBROUTINE wrt_ships
