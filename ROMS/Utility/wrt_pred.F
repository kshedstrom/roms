#include "cppdefs.h"
#if defined NEMURO_SAN && defined PREDATOR
      SUBROUTINE wrt_pred (ng)
!
!svn $Id: wrt_pred.F 975 2009-05-05 22:51:13Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes simulated drifter trajectories into pred     !
!  NetCDF file.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_pred
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

      real(r8), dimension(Npred(ng)) :: Tout
!
      SourceFile='wrt_pred.F'
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set time record index.
!
      tPREDindx(ng)=tPREDindx(ng)+1
      NrecPRED(ng)=NrecPRED(ng)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/tPREDindx(ng)/), (/1/),                   &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out pred species.
!
      CALL netcdf_put_ivar (ng, iNLM, PREDname(ng),                     &
     &                      TRIM(Vname(1,idPspc)),                      &
     &                      PREDS(ng) % species,                        &
     &                      (/1/), (/Npred(ng)/),                       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idPspc,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out pred aliveness.
!
      CALL netcdf_put_lvar (ng, iNLM, PREDname(ng),                     &
     &                      TRIM(Vname(1,idPalv)),                      &
     &                      PREDS(ng) % alive,                          &
     &                      (/1/), (/Npred(ng)/),                       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idPalv,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out pred X-grid locations.
!
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%track(ixgrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      'Xgrid', Tout,                              &
     &                      (/1,tPREDindx(ng)/), (/Npred(ng),1/),       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idXgrd,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out pred Y-grid locations.
!
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%track(iygrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      'Ygrid', Tout,                              &
     &                      (/1,tPREDindx(ng)/), (/Npred(ng),1/),       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idYgrd,ng))
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Write out pred Z-grid locations.
!
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%track(izgrd,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      'Zgrid', Tout,                              &
     &                      (/1,tPREDindx(ng)/), (/Npred(ng),1/),       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idZgrd,ng))
      IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Write out pred (lon,lat) or (x,y) locations.
!
      DO l=1,Npred(ng)
        Tout(l)=PREDS(ng)%track(iflon,nf(ng),l)
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                   &
     &                        'lon', Tout,                              &
     &                        (/1,tPREDindx(ng)/), (/Npred(ng),1/),     &
     &                        ncid = ncPREDid(ng),                      &
     &                        varid = predVid(idglon,ng))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                   &
     &                        'x', Tout,                                &
     &                        (/1,tPREDindx(ng)/), (/Npred(ng),1/),     &
     &                        ncid = ncPREDid(ng),                      &
     &                        varid = predVid(idglon,ng))
      END IF
      IF (exit_flag.ne.NoError) RETURN
!
      DO l=1,Npred(ng)
        Tout(l)=PREDS(ng)%track(iflat,nf(ng),l)
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                   &
     &                        'lat', Tout,                              &
     &                        (/1,tPREDindx(ng)/), (/Npred(ng),1/),     &
     &                        ncid = ncPREDid(ng),                      &
     &                        varid = predVid(idglat,ng))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                   &
     &                        'y', Tout,                                &
     &                        (/1,tPREDindx(ng)/), (/Npred(ng),1/),     &
     &                        ncid = ncPREDid(ng),                      &
     &                        varid = predVid(idglat,ng))
      END IF
      IF (exit_flag.ne.NoError) RETURN

# ifdef SOLVE3D
!
!  Write out pred depths.
!
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%track(idpth,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      'depth', Tout,                              &
     &                      (/1,tPREDindx(ng)/), (/Npred(ng),1/),       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(iddpth,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out density anomaly.
!
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%track(ifden,nf(ng),l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      TRIM(Vname(1,idDano)), Tout,                &
     &                      (/1,tPREDindx(ng)/), (/Npred(ng),1/),       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idDano,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        DO l=1,Npred(ng)
          IF (PREDS(ng)%bounded(l)) THEN
            Tout(l)=PREDS(ng)%track(itrc+10,nf(ng),l)
          ELSE
            Tout(l)=spval
          END IF
        END DO
        CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                   &
     &                        TRIM(Vname(1,idTvar(itrc))), Tout,        &
     &                        (/1,tPREDindx(ng)/), (/Npred(ng),1/),     &
     &                        ncid = ncPREDid(ng),                      &
     &                        varid = predTid(itrc,ng))
        IF (exit_flag.ne.NoError) RETURN
      END DO
#  if defined NEMURO_SAN && defined PREDATOR
!  Write pred wet weight.
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ipwwt,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      TRIM(Vname(1,idPwwt)), Tout,                &
     &                      (/1,tPREDindx(ng)/), (/Npred(ng),1/),       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idPwwt,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write pred worth.
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ipworth,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      TRIM(Vname(1,idPworth)), Tout,              &
     &                      (/1,tPREDindx(ng)/), (/Npred(ng),1/),       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idPworth,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write pred P-value.
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ippval,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      TRIM(Vname(1,idPpval)), Tout,               &
     &                      (/1,tPREDindx(ng)/), (/Npred(ng),1/),       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idPpval,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write Pred consumption on species 1
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ipcsmF1,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      TRIM(Vname(1,idPcsmF1)), Tout,              &
     &                      (/1,tPREDindx(ng)/), (/Npred(ng),1/),       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idPcsmF1,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write Pred consumption on species 2
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ipcsmF2,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      TRIM(Vname(1,idPcsmF2)), Tout,              &
     &                      (/1,tPREDindx(ng)/), (/Npred(ng),1/),       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idPcsmF2,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write Pred consumption on species 3
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ipcsmF3,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      TRIM(Vname(1,idPcsmF3)), Tout,              &
     &                      (/1,tPREDindx(ng)/), (/Npred(ng),1/),       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idPcsmF3,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write Pred consumption on species 4
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ipcsmF4,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      TRIM(Vname(1,idPcsmF4)), Tout,              &
     &                      (/1,tPREDindx(ng)/), (/Npred(ng),1/),       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idPcsmF4,ng))
      IF (exit_flag.ne.NoError) RETURN
!  Write Pred consumption on species 5
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ipcsmF5,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PREDname(ng),                     &
     &                      TRIM(Vname(1,idPcsmF5)), Tout,              &
     &                      (/1,tPREDindx(ng)/), (/Npred(ng),1/),       &
     &                      ncid = ncPREDid(ng),                        &
     &                      varid = predVid(idPcsmF5,ng))
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize pred NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, PREDname(ng), ncPREDid(ng))
!
#else
      SUBROUTINE wrt_pred
#endif
      RETURN
      END SUBROUTINE wrt_pred
