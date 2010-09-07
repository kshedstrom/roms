#include "cppdefs.h"
#ifdef FILTERED_RST
      SUBROUTINE wrt_filt (ng, filterLevel)
!
!=======================================================================
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine writes fields into filter restart NetCDF files.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_filter
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, filterLevel
!
!  Local variable declarations.
!
      integer :: LBi, UBi, LBj, UBj
      integer :: gfactor, gtype, i, itrc, status, ifile, varid
      integer :: start(10), total(10)
      real(r8) :: scale

      character (len=256) :: ncname
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out restart fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
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
      tFILindx(ng)=filterLevel
      NrecFIL(ng)=NrecFIL(ng)+1
!
!  Write out model time (s).
!
      DO ifile=1,nfile
        ncname=FILname(ifile,ng)
        CALL netcdf_put_ivar (ng, iNLM, ncname, 'filterLevel',          &
     &                        tFILindx(ng), (/tFILindx(ng)/), (/1/),    &
     &                        ncid = ncFILid(ifile,ng))
        IF (exit_flag.ne.NoError) RETURN
!
!  Write out averaged time.
!
        CALL netcdf_put_fvar (ng, iNLM, ncname,                         &
     &                        TRIM(Vname(idtime,ng)), time(ng),         &
     &                        (/tFILindx(ng)/), (/1/),                  &
     &                        ncid = ncFILid(ifile,ng) )
        IF (exit_flag.ne.NoError) RETURN

        CALL netcdf_put_ivar (ng, iNLM, ncname, 'fcount',               &
     &                        fcount(tFILindx(ng)), (/tFILindx(ng)/),   &
     &                        (/1/), ncid = ncFILid(ifile,ng))
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Write out free-surface (m)
!
      IF (Aout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), filVid(idFsur,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
# ifdef MASKING
     &                 GRID(ng) % rmask,                                &
# endif
     &                 FILTER(ng) % filzeta(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), filVid(idUbar,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
# ifdef MASKING
     &                 GRID(ng) % umask,                                &
# endif
     &                 FILTER(ng) % filu2d(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), filVid(idVbar,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
# ifdef MASKING
     &                 GRID(ng) % vmask,                                &
# endif
     &                 FILTER(ng) % filv2d(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Aout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncFILid(2,ng), filVid(idUvel,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % umask,                                &
#  endif
     &                 FILTER(ng) % filu3d(:,:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Aout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncFILid(3,ng), filVid(idVvel,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % vmask,                                &
#  endif
     &                 FILTER(ng) % filv3d(:,:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out density anomaly.
!
      IF (Aout(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncFILid(7,ng), filVid(idDano,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filrho(:,:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NAT
        IF (Aout(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, ncFILid(8,ng), filTid(itrc,ng),  &
     &                   tFILindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#  ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#  endif
     &                   FILTER(ng) % filt(:,:,:,itrc,tFILindx(ng)))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), tFILindx(ng)
              PRINT *, trim(nf90_strerror(status))
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out vertical (omega) velocity
!
      IF (Aout(idOvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncFILid(9,ng), filVid(idOvel,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filw3d(:,:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
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
      IF (Aout(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), filVid(idHsbl,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filhsbl(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
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
      IF (Aout(idHbbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), filVid(idHbbl,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filhbbl(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHbbl)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Aout(idVvis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncFILid(6,ng), filVid(idVvis,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filAKv(:,:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Aout(idTdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncFILid(4,ng), filVid(idTdif,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filAKt(:,:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Aout(idSdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncFILid(5,ng), filVid(idSdif,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filAKs(:,:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface net heat flux.
!
      IF (Aout(idTsur(itemp),ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng),                     &
     &                 filVid(idTsur(itemp),ng),                        &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filstf(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface net salt flux.
!
      IF (Aout(idTsur(isalt),ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng),                     &
     &                 filVid(idTsur(isalt),ng),                        &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filswf(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTsur(isalt))), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Aout(idLhea,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), filVid(idLhea,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % fillhf(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Aout(idShea,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), filVid(idShea,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filshf(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Aout(idLrad,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), filVid(idLrad,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % fillrf(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
#  ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Aout(idSrad,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), filVid(idSrad,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filsrf(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
!
!  Write out surface u-momentum stress.
!
      IF (Aout(idUsms,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), filVid(idUsms,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
# ifdef MASKING
     &                 GRID(ng) % umask,                                &
# endif
     &                 FILTER(ng) % filsus(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface v-momentum stress.
!
      IF (Aout(idVsms,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), filVid(idVsms,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
# ifdef MASKING
     &                 GRID(ng) % vmask,                                &
# endif
     &                 FILTER(ng) % filsvs(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef ICE_MODEL
!
!  Write out ice u-velocity
!
      IF (Aout(idUice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idUice,ng), &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % umask,                                &
#  endif
     &                 FILTER(ng) % filui(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUice)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice v-velocity
!
      IF (Aout(idVice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idVice,ng), &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % vmask,                                &
#  endif
     &                 FILTER(ng) % filvi(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVice)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice concentration
!
      IF (Aout(idAice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idAice,ng), &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filai(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAice)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice thickness
!
      IF (Aout(idHice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idHice,ng), &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filhi(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHice)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out snow thickness
!
      IF (Aout(idHsno,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idHsno,ng), &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filhsn(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsno)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface water
!
      IF (Aout(idSfwat,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idSfwat,ng),&
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filsfwat(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSfwat)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface ice temperature
!
      IF (Aout(idTice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idTice,ng), &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filtis(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTice)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out interior ice temperature
!
      IF (Aout(idTimid,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idTimid,ng),&
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filti(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTimid)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice age
!
      IF (Aout(idAgeice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng),                    &
     &                 filVid(idAgeice,ng), &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filAgeice(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAgeice)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out under-ice salinity
!
      IF (Aout(idS0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idS0mk,ng), &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % fils0mk(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idS0mk)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out under-ice temperature
!
      IF (Aout(idT0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idT0mk,ng), &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filt0mk(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idT0mk)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out frazil ice growth
!
      IF (Aout(idWfr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idWfr,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filWfr(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWfr)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice growth/melt
!
      IF (Aout(idWai,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idWai,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filWai(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWai)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice growth/melt
!
      IF (Aout(idWao,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idWao,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filWao(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWao)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice growth/melt
!
      IF (Aout(idWio,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idWio,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filWio(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWio)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice melt runoff
!
      IF (Aout(idWro,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idWro,ng),  &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filWro(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWro)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sig11 of ice stress tensor
!
      IF (Aout(idSig11,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idSig11,ng),&
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filsig11(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig11)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sig12 of ice stress tensor
!
      IF (Aout(idSig12,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idSig12,ng),&
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filsig12(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig12)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sig22 of ice stress tensor
!
      IF (Aout(idSig22,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idSig22,ng),&
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filsig22(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig22)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice-ocean mass flux
!
      IF (Aout(idIomflx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng),                    &
     &                 filVid(idIomflx,ng),                             &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filio_mflux(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIomflx)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice-water tau
!
      IF (Aout(idTauiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idTauiw,ng),&
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filutau_iw(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTauiw)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice-water chu
!
      IF (Aout(idChuiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), filVid(idChuiw,ng),&
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#  endif
     &                 FILTER(ng) % filchu_iw(:,:,tFILindx(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idChuiw)), tFILindx(ng)
            PRINT *, trim(nf90_strerror(status))
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Synchronize history NetCDF file to disk.
!
      DO ifile=1,nfile
        CALL netcdf_sync(ng, iNLM, FILname(ifile,ng), ncFILid(ifile,ng))
        IF (exit_flag.ne.NoError) RETURN
      END DO
      IF (Master) WRITE (stdout,20) tFILindx(ng)
!
  10  FORMAT (/,' WRT_FILT  - error while writing variable: ',a,/,11x,  &
     &        'into filter NetCDF file for time record: ',i4)
  20  FORMAT (6x,'WRT_FILT  - wrote filter fields into time ',          &
     &        'record =',t72,i7.7)
      RETURN
      END SUBROUTINE wrt_filt
#else
      SUBROUTINE wrt_filt
      RETURN
      END SUBROUTINE wrt_filt
#endif
