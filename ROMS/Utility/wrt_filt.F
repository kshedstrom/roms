#include "cppdefs.h"
#ifdef FILTERED_RST
      SUBROUTINE wrt_filt (ng, filterLevel)
!
!=======================================================================
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
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
# ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_bcasti
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
        IF (OutThread) THEN
          status=nf90_inq_varid(ncFILid(ifile,ng), TRIM('filterLevel'), &
     &                        varid)
          start(1)=tFILindx(ng)
          total(1)=1
          status=nf90_put_var(ncFILid(ifile,ng), varid, tFILindx(ng),   &
     &                        start)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'filterLevel', tFILindx(ng)
            exit_flag=3
            ioerror=status
          END IF
          status=nf90_inq_varid(ncFILid(ifile,ng), TRIM('ocean_time'),  &
     &                        varid)
          status=nf90_put_var(ncFILid(ifile,ng), varid, time(ng))
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'ocean_time', tFILindx(ng)
            exit_flag=3
            ioerror=status
          END IF
          status=nf90_inq_varid(ncFILid(ifile,ng), TRIM('fcount'),      &
     &                        varid)
          start(1)=tFILindx(ng)
          total(1)=1
          status=nf90_put_var(ncFILid(ifile,ng), varid,                 &
     &                        fcount(tFILindx(ng)), start)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) 'fcount', tFILindx(ng)
            exit_flag=3
            ioerror=status
          END IF
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
        IF (exit_flag.ne.NoError) RETURN
      END DO
!
!  Write out free-surface (m)
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idFsur)), varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
# ifdef MASKING
     &                 GRID(ng) % rmask,                                &
# endif
     &                 FILTER(ng) % filzeta(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idFsur)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idUbar)), varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
# ifdef MASKING
     &                 GRID(ng) % umask,                                &
# endif
     &                 FILTER(ng) % filu2d(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUbar)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idVbar)), varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
# ifdef MASKING
     &                 GRID(ng) % vmask,                                &
# endif
     &                 FILTER(ng) % filv2d(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVbar)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# ifdef SOLVE3D
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      scale=1.0_r8
      gtype=gfactor*u3dvar
      status=nf90_inq_varid(ncFILid(2,ng), TRIM(Vname(1,idUvel)), varid)
      status=nf_fwrite3d(ng, iNLM, ncFILid(2,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % umask,                                &
#  endif
     &                 FILTER(ng) % filu3d(:,:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUvel)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      scale=1.0_r8
      gtype=gfactor*v3dvar
      status=nf90_inq_varid(ncFILid(3,ng), TRIM(Vname(1,idVvel)), varid)
      status=nf_fwrite3d(ng, iNLM, ncFILid(3,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#  ifdef MASKING
     &                 GRID(ng) % vmask,                                &
#  endif
     &                 FILTER(ng) % filv3d(:,:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVvel)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  ifndef FILTRIM
!
!  Write out density anomaly.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf90_inq_varid(ncFILid(7,ng), TRIM(Vname(1,idDano)), varid)
      status=nf_fwrite3d(ng, iNLM, ncFILid(7,ng), filVid(idDano,ng),    &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 1, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filrho(:,:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idDano)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO itrc=1,NAT
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncFILid(8,ng), filTid(itrc,ng),    &
     &                   tFILindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, 1, N(ng), scale,           &
#   ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#   endif
     &                   FILTER(ng) % filt(:,:,:,itrc,tFILindx(ng)))
        IF (OutThread.and.(status.ne.nf90_noerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))), tFILindx(ng)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END DO
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
!
!  Write out vertical (omega) velocity
!
      scale=1.0_r8
      gtype=gfactor*w3dvar
      status=nf90_inq_varid(ncFILid(9,ng), TRIM(Vname(1,idOvel)), varid)
      status=nf_fwrite3d(ng, iNLM, ncFILid(9,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filw3d(:,:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idOvel)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef LMD_SKPP
!
!  Write out depth of surface boundary layer.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idHsbl)), varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filhsbl(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Write out depth of bottom boundary layer.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idHbbl)), varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filhbbl(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idHbbl)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_AKV
!
!  Write out vertical viscosity coefficient.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf90_inq_varid(ncFILid(6,ng), TRIM(Vname(1,idVvis)), varid)
      status=nf_fwrite3d(ng, iNLM, ncFILid(6,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filAKv(:,:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVvis)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_AKT
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf90_inq_varid(ncFILid(4,ng), TRIM(Vname(1,idTdif)), varid)
      status=nf_fwrite3d(ng, iNLM, ncFILid(4,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filAKt(:,:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTdif)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_AKS
!
!  Write out vertical diffusion coefficient for salinity.
!
      scale=1.0_r8
      gtype=gfactor*r3dvar
      status=nf90_inq_varid(ncFILid(5,ng), TRIM(Vname(1,idSdif)), varid)
      status=nf_fwrite3d(ng, iNLM, ncFILid(5,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, 0, N(ng), scale,             &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filAKs(:,:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSdif)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_FLUXES
!
!  Write out surface net heat flux.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idTsur(itemp))),&
     &                    varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filstf(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTsur(itemp))), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface net salt flux.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idTsur(isalt))),&
     &                    varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#   ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#   endif
     &                 FILTER(ng) % filswf(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTsur(isalt))), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#   ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idLhea)), varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % fillhf(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idLhea)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out sensible heat flux.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idShea)), varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filshf(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idShea)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out longwave radiation flux.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idLrad)), varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % fillrf(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idLrad)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#   endif
#   ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idSrad)), varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filsrf(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSrad)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
!
!  Write out surface u-momentum stress.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idUsms)), varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % umask,                                &
#  endif
     &                 FILTER(ng) % filsus(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUsms)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface v-momentum stress.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idVsms)), varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(1,ng), varid,                &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % vmask,                                &
#  endif
     &                 FILTER(ng) % filsvs(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVsms)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
# endif
# ifdef ICE_MODEL
!
!  Write out ice u-velocity
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idUice)),      &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % umask,                                &
#  endif
     &                 FILTER(ng) % filui(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idUice)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice v-velocity
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idVice)),      &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#  ifdef MASKING
     &                 GRID(ng) % vmask,                                &
#  endif
     &                 FILTER(ng) % filvi(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idVice)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice concentration
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idAice)),      &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filai(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idAice)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice thickness
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idHice)),      &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filhi(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idHice)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out snow thickness
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idHsno)),      &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filhsn(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idHsno)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface water
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idSfwat)),     &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filsfwat(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSfwat)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out surface ice temperature
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idTice)),      &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filtis(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTice)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out interior ice temperature
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idTimid)),     &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filti(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTimid)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice age
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idAgeice)),    &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filAgeice(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idAgeice)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out under-ice salinity
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idS0mk)),      &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % fils0mk(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idS0mk)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out under-ice temperature
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idT0mk)),      &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filt0mk(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idT0mk)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out frazil ice growth
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idWfr)),       &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filWfr(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idWfr)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice growth/melt
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idWai)),       &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filWai(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idWai)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice growth/melt
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idWao)),       &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filWao(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idWao)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice growth/melt
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idWio)),       &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filWio(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idWio)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice melt runoff
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idWro)),       &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filWro(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idWro)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out sig11 of ice stress tensor
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idSig11)),     &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filsig11(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSig11)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out sig12 of ice stress tensor
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idSig12)),     &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filsig12(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSig12)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out sig22 of ice stress tensor
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idSig22)),     &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filsig22(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idSig22)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice-ocean mass flux
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idIomflx)),    &
     &             varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filio_mflux(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idIomflx)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice-water tau
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idTauiw)),     &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filutau_iw(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idTauiw)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out ice-water chu
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idChuiw)),     &
     &                 varid)
      status=nf_fwrite2d(ng, iNLM, ncFILid(10,ng), varid,               &
     &                 tFILindx(ng), gtype,                             &
     &                 LBi, UBi, LBj, UBj, scale,                       &
#    ifdef MASKING
     &                 GRID(ng) % rmask,                                &
#    endif
     &                 FILTER(ng) % filchu_iw(:,:,tFILindx(ng)))
      IF (OutThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,10) TRIM(Vname(1,idChuiw)), tFILindx(ng)
        exit_flag=3
        ioerror=status
        RETURN
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
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, iNLM, exit_flag)
# endif
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
