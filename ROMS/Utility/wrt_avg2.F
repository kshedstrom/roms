#include "cppdefs.h"
#ifdef AVERAGES2
#undef AVERAGES_DETIDE
      SUBROUTINE wrt_avg2 (ng)
!
!svn $Id: wrt_avg.F 702 2008-08-12 16:44:47Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes model time-averaged fields into averages     !
!  NetCDF file.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_average2
      USE mod_forces
      USE mod_grid
#ifdef ICE_MODEL
      USE mod_ice
#endif
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
# if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
# endif
# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
      USE mod_tides
# endif
# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
# endif
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: gfactor, gtype, i, itrc, status

      real(r8) :: scale
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out time-averaged fields when appropriate.
!-----------------------------------------------------------------------
!
      if (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
        gfactor=-1
# else
        gfactor=1
# endif
!
!  Set time record index.
!
      tAVG2indx(ng)=tAVG2indx(ng)+1
      NrecAVG2(ng)=NrecAVG2(ng)+1
!
!  Write out averaged time.
!
      CALL netcdf_put_fvar (ng, iNLM, AVG2name(ng),                      &
     &                      TRIM(Vname(idtime,ng)), AVG2time(ng:),       &
     &                      (/tAVG2indx(ng)/), (/1/),                    &
     &                      ncid = ncAVG2id(ng),                         &
     &                      varid = avg2Vid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out free-surface (m).
!
      IF (Hout2(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idFsur,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask,                   &
# endif
     &                     AVERAGE2(ng) % avgzeta)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided free-surface (m).
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idFsuD,ng),      &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#  endif
     &                   TIDES(ng) % zeta_detided)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idFsuD)), tAVG2indx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout2(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idUbar,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask,                   &
# endif
     &                     AVERAGE2(ng) % avgu2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idu2dD,ng),      &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask,                     &
#  endif
     &                   TIDES(ng) % ubar_detided)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idu2dD)), tAVG2indx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout2(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idVbar,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask,                   &
# endif
     &                     AVERAGE2(ng) % avgv2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idv2dD,ng),      &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                     &
#  endif
     &                   TIDES(ng) % vbar_detided)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idv2dD)), tAVG2indx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# if defined FORWARD_WRITE && defined SOLVE3D
!
!  Write out time-averaged mass fluxes for 3D momentum coupling.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idUfx1,ng),      &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask,                     &
#  endif
     &                   AVERAGE2(ng) % avgDU_avg1)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUfx1)), tAVG2indx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idUfx2,ng),      &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask,                     &
#  endif
     &                   AVERAGE2(ng) % avgDU_avg2)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUfx2)), tAVG2indx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idVfx1,ng),      &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % vmask,                     &
#  endif
     &                   AVERAGE2(ng) % avgDV_avg1)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVfx1)), tAVG2indx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idVfx2,ng),      &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#  ifdef MASKING
     &                   GRID(ng) % umask,                     &
#  endif
     &                   AVERAGE2(ng) % avgDV_avg2)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVfx2)), tAVG2indx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef AVERAGES_NEARSHORE
!
!  Write out 2D radiation stress, Sxx-component.
!
      IF (Hout2(idW2xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idW2xx,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#  endif
     &                     AVERAGE2(ng) % avgSxx2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xx)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Sxy-component.
!
      IF (Hout2(idW2xy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idW2xy,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#  endif
     &                     AVERAGE2(ng) % avgSxy2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xy)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Syy-component.
!
      IF (Hout2(idW2yy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idW2yy,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#  endif
     &                     AVERAGE2(ng) % avgSyy2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2yy)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress in the XI-direction.
!
      IF (Hout2(idU2rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idU2rs,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                   &
#  endif
     &                     AVERAGE2(ng) % avgu2RS)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2rs)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress in the ETA-direction
!
      IF (Hout2(idV2rs,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idV2rs,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                   &
#  endif
     &                     AVERAGE2(ng) % avgv2RS)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2rs)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D stokes momentum component (m/s) in the XI-direction.
!
      IF (Hout2(idU2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idU2Sd,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                   &
#  endif
     &                     AVERAGE2(ng) % avgu2Sd)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2Sd)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D stokes momentum component (m/s) in the ETA-direction.
!
      IF (Hout2(idV2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idV2Sd,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                   &
#  endif
     &                     AVERAGE2(ng) % avgv2Sd)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2Sd)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Hout2(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idUvel,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                   &
#  endif
     &                     AVERAGE2(ng) % avgu3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

#  if defined AVERAGES_DETIDE && (defined SSH_TIDES || defined UV_TIDES)
!
!  Write out detided 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idu3dD,ng),      &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % umask,                     &
#   endif
     &                   TIDES(ng) % u_detided)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idu3dD)), tAVG2indx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Hout2(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idVvel,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                   &
#  endif
     &                     AVERAGE2(ng) % avgv3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NT(ng)
        IF (Hout2(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Tid(itrc,ng),    &
     &                       tAVG2indx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                 &
#  endif
     &                       AVERAGE2(ng) % avgt(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                        tAVG2indx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# ifdef BIO_GOANPZ
!---------------------------------------------
!  Write out stationary tracer type variables.
!---------------------------------------------
      DO itrc=1,NTS(ng)
	!print *, 'XXX',ncAVG2id(ng),avg2TSid(itrc,ng),tAVG2indx(ng)
	!print *, 'xxx',Hout2(idTSvar(itrc),ng)
        IF (Hout2(idTSvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2TSid(itrc,ng),    &
     &                       tAVG2indx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                 &
#  endif
     &                       AVERAGE2(ng) % avgst(:,:,1,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTSvar(itrc))),            &
     &                        tAVG2indx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
!
!  Write out density anomaly.
!
      IF (Hout2(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idDano,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#  endif
     &                     AVERAGE2(ng) % avgrho)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      IF (Hout2(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idHsbl,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#   endif
     &                     AVERAGE2(ng) % avghsbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      IF (Hout2(idHbbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idHbbl,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#   endif
     &                     AVERAGE2(ng) % avghbbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHbbl)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef ICE_MODEL
!
!  Write out ice 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout2(idUice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idUice,ng),    &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % umask,                     &
#   endif
     &                   AVERAGE2(ng) % avguice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUice)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
!  Write out ice 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout2(idVice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idVice,ng),    &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % vmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgvice)
       IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVice)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice concentration.
!
      IF (Hout2(idAice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idAice,ng),    &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgaice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAice)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice thickness.
!
      IF (Hout2(idHice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idHice,ng),    &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avghice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHice)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice/snow surface temperature.
!
      IF (Hout2(idTice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idTice,ng),    &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgtice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTice)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice internal temperature.
!
      IF (Hout2(idTimid,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idTimid,ng),   &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgtimid)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTimid)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out snow thickness.
!
      IF (Hout2(idHsno,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idHsno,ng),    &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avghsno)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsno)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out surface water thickness (on ice).
!
      IF (Hout2(idSfwat,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idSfwat,ng),   &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgsfwat)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSfwat)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice age.
!
      IF (Hout2(idAgeice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idAgeice,ng),   &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgageice)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAgeice)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice-ocean mass flux
!
      IF (Hout2(idIomflx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idIomflx,ng),  &

     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgiomflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIomflx)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out internal ice stress component 11
!
      IF (Hout2(idSig11,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idSig11,ng),   &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgsig11)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig11)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out internal ice stress component 12
!
      IF (Hout2(idSig12,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idSig12,ng),   &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgsig12)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig12)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out internal ice stress component 22
!
      IF (Hout2(idSig22,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idSig22,ng),   &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgsig22)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig22)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out temperature of molecular sublayer under ice.
!
      IF (Hout2(idT0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idT0mk,ng),    &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgT0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idT0mk)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out salinity of molecular sublayer under ice.
!
      IF (Hout2(idS0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idS0mk,ng),    &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgS0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idS0mk)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice-water friction velocity
!
      IF (Hout2(idTauiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idTauiw,ng),    &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgutau_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTauiw)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
!
! Write out ice-water momentum transfer coefficient
!
      IF (Hout2(idChuiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idChuiw,ng),    &
     &                   tAVG2indx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                     &
#   endif
     &                   AVERAGE2(ng) % avgchu_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idChuiw)), tAVG2indx(ng)
          END IF
          exit_flag=3
          RETURN
        END IF
      END IF
#  endif
# endif
# ifdef SOLVE3D
#  ifdef AVERAGES_FLUXES
!
!  Write out surface net heat flux.
!
      IF (Hout2(idTsur(itemp),ng)) THEN
#   ifdef ADJOINT
#    if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
        defined OPT_OBSERVATIONS
        scale=1.0_r8/(rho0*Cp)
#    else
        scale=1.0_r8
#    endif
#   else
        scale=rho0*Cp
#   endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng),                       &
     &                     avg2Vid(idTsur(itemp),ng),                    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#   endif
     &                     AVERAGE2(ng) % avgstf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface salt flux  (PSU m/s = kg salt/m2/s).
!
      IF (Hout2(idTsur(isalt),ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng),                       &
     &                     avg2Vid(idTsur(isalt),ng),                    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#   endif
     &                     AVERAGE2(ng) % avgswf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(isalt))), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Hout2(idLhea,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_r8
#    else
        scale=rho0*Cp
#    endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idLhea,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#    endif
     &                     AVERAGE2(ng) % avglhf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout2(idShea,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_r8
#    else
        scale=rho0*Cp
#    endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idShea,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#    endif
     &                     AVERAGE2(ng) % avgshf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Hout2(idLrad,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_r8
#    else
        scale=rho0*Cp
#    endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idLrad,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#    endif
     &                     AVERAGE2(ng) % avglrf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface u-wind.
!
      IF (Hout2(idUair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idUair,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#  endif
     &                     AVERAGE2(ng) % avguwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUair)), tAVG2indx(ng)
          END IF                    
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-wind.
!
      IF (Hout2(idVair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idVair,ng),  &
     &                     tAVG2indx(ng), gtype,                        &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#  endif
     &                     AVERAGE2(ng) % avgvwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVair)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#    ifdef EMINUSP
!
!  Write out evaportaion rate (kg/m2/s).
!
      IF (Hout2(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idevap,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#     ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#     endif
     &                     AVERAGE2(ng) % avgevap)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout2(idrain,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idrain,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#     ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#     endif
     &                     AVERAGE2(ng) % avgrain)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#    endif
#   endif
#   ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Hout2(idSrad,ng)) THEN
#    ifdef ADJOINT
        scale=1.0_r8
#    else
        scale=rho0*Cp
#    endif
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idSrad,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#    ifdef MASKING
     &                     GRID(ng) % rmask,                   &
#    endif
     &                     AVERAGE2(ng) % avgsrf)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), tAVG2indx(ng)
          END IF          
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
!
!  Write out surface u-momentun stress.
!
      IF (Hout2(idUsms,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#   else
        scale=1.0_r8
#   endif
#  else
        scale=rho0
#  endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idUsms,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                   &
#  endif
     &                     AVERAGE2(ng) % avgsus)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), tAVG2indx(ng)
          END IF                    
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentun stress.
!
      IF (Hout2(idVsms,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#   else
        scale=1.0_r8
#   endif
#  else
        scale=rho0
#  endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idVsms,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                   &
#  endif
     &                     AVERAGE2(ng) % avgsvs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom u-momentun stress.
!
      IF (Hout2(idUbms,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#   else
        scale=1.0_r8
#   endif
#  else
        scale=rho0
#  endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idUbms,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask,                   &
#  endif
     &                     AVERAGE2(ng) % avgbus)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), tAVG2indx(ng)
          END IF                    
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom v-momentun stress.
!
      IF (Hout2(idVbms,ng)) THEN
#  ifdef ADJOINT
#   if defined AD_SENSITIVITY   || defined OBS_SENSITIVITY || \
       defined OPT_OBSERVATIONS
        scale=1.0_r8/rho0
#   else
        scale=1.0_r8
#   endif
#  else
        scale=rho0
#  endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng), avg2Vid(idVbms,ng),    &
     &                     tAVG2indx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask,                   &
#  endif
     &                     AVERAGE2(ng) % avgbvs)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), tAVG2indx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined SEDIMENT && defined BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO itrc=1,NST
        IF (Hout2(idUbld(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*u2dvar
          status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng),                     &
     &                       avg2Vid(idUbld(itrc),ng),                   &
     &                       tAVG2indx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % umask,                 &
#  endif
     &                       AVERAGE2(ng) % avgbedldu(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbld(itrc))),            &
     &                          tAVG2indx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out bed load transport in V-direction.
!
        IF (Hout2(idVbld(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*v2dvar
          status=nf_fwrite2d(ng, iNLM, ncAVG2id(ng),                     &
     &                       avg2Vid(idVbld(itrc),ng),                   &
     &                       tAVG2indx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % vmask,                 &
#  endif
     &                       AVERAGE2(ng) % avgbedldv(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbld(itrc))),            &
     &                          tAVG2indx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
!
!-----------------------------------------------------------------------
!  Synchronize time-average NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, AVG2name(ng), ncAVG2id(ng))
      IF (exit_flag.ne.NoError) RETURN

      IF (Master) WRITE (stdout,20) tAVG2indx(ng)
!
  10  FORMAT (/,' WRT_AVG2 - error while writing variable: ',a,/,11x,    &
     &        'into averages NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_AVG2   - wrote averaged fields into time ',        &
     &        'record =',t72,i7.7)
#else
      SUBROUTINE wrt_avg2
#endif
      RETURN
      END SUBROUTINE wrt_avg2
