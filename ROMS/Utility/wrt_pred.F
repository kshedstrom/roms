#include "cppdefs.h"
#if defined NEMURO_SAN && defined PREDATOR
      SUBROUTINE wrt_pred (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
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
      PRED(ng)%Rindex=PRED(ng)%Rindex+1
      PRED(ng)%Nrec=PRED(ng)%Nrec+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/PRED(ng)%Rindex/), (/1/),                 &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out pred species.
!
      CALL netcdf_put_ivar (ng, iNLM, PRED(ng)%name,                    &
     &                      TRIM(Vname(1,idPspc)),                      &
     &                      PREDS(ng) % species,                        &
     &                      (/1/), (/Npred(ng)/),                       &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idPspc))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out pred aliveness.
!
      CALL netcdf_put_lvar (ng, iNLM, PRED(ng)%name,                    &
     &                      TRIM(Vname(1,idPalv)),                      &
     &                      PREDS(ng) % alive,                          &
     &                      (/1/), (/Npred(ng)/),                       &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idPalv))
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
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      'Xgrid', Tout,                              &
     &                      (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),     &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idXgrd))
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
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      'Ygrid', Tout,                              &
     &                      (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),     &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idYgrd))
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
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      'Zgrid', Tout,                              &
     &                      (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),     &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idZgrd))
      IF (exit_flag.ne.NoError) RETURN
# endif
!
!  Write out pred (lon,lat) or (x,y) locations.
!
      DO l=1,Npred(ng)
        Tout(l)=PREDS(ng)%track(iflon,nf(ng),l)
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                  &
     &                        'lon', Tout,                              &
     &                        (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),   &
     &                        ncid = PRED(ng)%ncid,                     &
     &                        varid = PRED(ng)%Vid(idglon))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                  &
     &                        'x', Tout,                                &
     &                        (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),   &
     &                        ncid = PRED(ng)%ncid,                     &
     &                        varid = PRED(ng)%Vid(idglon))
      END IF
      IF (exit_flag.ne.NoError) RETURN
!
      DO l=1,Npred(ng)
        Tout(l)=PREDS(ng)%track(iflat,nf(ng),l)
      END DO
      IF (spherical) THEN
        CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                  &
     &                        'lat', Tout,                              &
     &                        (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),   &
     &                        ncid = PRED(ng)%ncid,                     &
     &                        varid = PRED(ng)%Vid(idglat))
      ELSE
        CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                  &
     &                        'y', Tout,                                &
     &                        (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),   &
     &                        ncid = PRED(ng)%ncid,                     &
     &                        varid = PRED(ng)%Vid(idglat))
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
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      'depth', Tout,                              &
     &                      (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),     &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(iddpth))
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
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      TRIM(Vname(1,idDano)), Tout,                &
     &                      (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),     &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idDano))
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
        CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                  &
     &                        TRIM(Vname(1,idTvar(itrc))), Tout,        &
     &                        (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),   &
     &                        ncid = PRED(ng)%ncid,                     &
     &                        varid = PRED(ng)%Tid(itrc))
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
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      TRIM(Vname(1,idPwwt)), Tout,                &
     &                      (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),     &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idPwwt))
      IF (exit_flag.ne.NoError) RETURN
!  Write pred worth.
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ipworth,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      TRIM(Vname(1,idPworth)), Tout,              &
     &                      (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),     &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idPworth))
      IF (exit_flag.ne.NoError) RETURN
!  Write pred P-value.
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ippval,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      TRIM(Vname(1,idPpval)), Tout,               &
     &                      (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),     &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idPpval))
      IF (exit_flag.ne.NoError) RETURN
!  Write Pred consumption on species 1
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ipcsmF1,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      TRIM(Vname(1,idPcsmF1)), Tout,              &
     &                      (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),     &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idPcsmF1))
      IF (exit_flag.ne.NoError) RETURN
!  Write Pred consumption on species 2
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ipcsmF2,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      TRIM(Vname(1,idPcsmF2)), Tout,              &
     &                      (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),     &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idPcsmF2))
      IF (exit_flag.ne.NoError) RETURN
!  Write Pred consumption on species 3
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ipcsmF3,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      TRIM(Vname(1,idPcsmF3)), Tout,              &
     &                      (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),     &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idPcsmF3))
      IF (exit_flag.ne.NoError) RETURN
!  Write Pred consumption on species 4
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ipcsmF4,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      TRIM(Vname(1,idPcsmF4)), Tout,              &
     &                      (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),     &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idPcsmF4))
      IF (exit_flag.ne.NoError) RETURN
!  Write Pred consumption on species 5
      DO l=1,Npred(ng)
        IF (PREDS(ng)%bounded(l)) THEN
          Tout(l)=PREDS(ng)%bioenergy(ipcsmF5,l)
        ELSE
          Tout(l)=spval
        END IF
      END DO
      CALL netcdf_put_fvar (ng, iNLM, PRED(ng)%name,                    &
     &                      TRIM(Vname(1,idPcsmF5)), Tout,              &
     &                      (/1,PRED(ng)%Rindex/), (/Npred(ng),1/),     &
     &                      ncid = PRED(ng)%ncid,                       &
     &                      varid = PRED(ng)%Vid(idPcsmF5))
      IF (exit_flag.ne.NoError) RETURN
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Synchronize pred NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, PRED(ng)%name, PRED(ng)%ncid)
!
#else
      SUBROUTINE wrt_pred
#endif
      RETURN
      END SUBROUTINE wrt_pred
