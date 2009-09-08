#include "cppdefs.h"
#ifdef FILTERED_RST
      SUBROUTINE get_filter (ng, filterLevel)
!
!=======================================================================
!  Copyright (c) 2002-2009 ROMS/TOMS Group                             !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine reads in primitive variables initial conditions     !
!  from an initial, restart, or history NetCDF file.                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_filter
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE nf_fread2d_mod, ONLY : nf_fread2d
      USE nf_fread3d_mod, ONLY : nf_fread3d

# if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
!
      USE exchange_2d_mod
#  ifdef SOLVE3D
      USE exchange_3d_mod
#  endif
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
      logical, dimension(NV) :: got_var

      integer :: LBi, UBi, LBj, UBj
      integer :: gtype, i, ifield, itrc, ifile, ifil
      integer :: nvatts, nrec, status, tindx, varid, vartype
      integer :: Vsize(4)

      integer :: start(4), total(4)

      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin, Tmax, time_scale, tmptime(Ngrids)

      character (len=15) attnam, tvarnam
      character (len=40) tunits
      character (len=80) fname, ncname
!
!-----------------------------------------------------------------------
!  Don't inquire about the contents of input NetCDF file,
!  the dimensions and variables, since multiple files
!-----------------------------------------------------------------------
!
      tindx = filterLevel
      ncname=FILname(1,ng)
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!  Set Vsize to zero to deativate interpolation of input data to model
!  grid in "nf_fread2d" and "nf_fread3d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
!
!-----------------------------------------------------------------------
!  Read initial conditions from input NetCDF file.
!-----------------------------------------------------------------------
!
!
!  Open input NetCDF file.
!
      DO ifile = 1,nfile
        ncname=FILname(ifile,ng)
        status=nf90_open(TRIM(ncname), nf90_nowrite, ncFILid(ifile,ng))
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,20) TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
      END DO
      IF (InpThread) THEN
        WRITE (stdout,90) 'Read filter restart values,', tindx
      END IF
!
!  Get the status of filter sums.
!
      status=nf90_inq_varid(ncFILid(1,ng), 'fcount', varid)
      start(1) = tindx
      total(1) = 1
      status=nf90_get_var(ncFILid(1,ng), varid, fcount(tindx),          &
     &              start)
!      print *, 'GET_FILTER fcount = ', fcount
!
!  Read in free-surface (m).
!
      ncname=FILname(1,ng)
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idFsur)), varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(1,ng),                &
     &                Vname(1,idFsur), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % rmask,                                 &
# endif
     &                FILTER(ng) % filzeta(:,:,tindx))
      IF (InpThread) THEN
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(Vname(1,idFsur)), tindx, TRIM(ncname)
          exit_flag=2
          RETURN
        ELSE
          WRITE (stdout,100) TRIM(Vname(2,idFsur)), Fmin, Fmax
        END IF
      END IF
!
!  Read in 2D momentum component (m/s) in the XI-direction.
!
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idUbar)), varid)
      gtype=var_flag(varid)*u2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(1,ng),                &
     &                Vname(1,idUbar), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % umask,                                 &
# endif
     &                FILTER(ng) % filu2d(:,:,tindx))
      IF (InpThread) THEN
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(Vname(1,idUbar)), tindx, TRIM(ncname)
          exit_flag=2
          RETURN
        ELSE
          WRITE (stdout,100) TRIM(Vname(2,idUbar)), Fmin, Fmax
        END IF
      END IF
!
!  Read in 2D momentum component (m/s) in the ETA-direction.
!
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idVbar)), varid)
      gtype=var_flag(varid)*v2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(1,ng),                &
     &                Vname(1,idVbar), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
# ifdef MASKING
     &                GRID(ng) % vmask,                                 &
# endif
     &                FILTER(ng) % filv2d(:,:,tindx))
      IF (InpThread) THEN
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(Vname(1,idVbar)), tindx, TRIM(ncname)
          exit_flag=2
          RETURN
        ELSE
          WRITE (stdout,100) TRIM(Vname(2,idVbar)), Fmin, Fmax
        END IF
      END IF
# ifdef SOLVE3D
!
!  Read in 3D momentum component (m/s) in the XI-direction.
!
      ncname=FILname(2,ng)
      status=nf90_inq_varid(ncFILid(2,ng), TRIM(Vname(1,idUvel)), varid)
      gtype=var_flag(varid)*u3dvar
      status=nf_fread3d(ng, iNLM, ncname, ncFILid(2,ng),                &
     &                Vname(1,idUvel), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % umask,                                 &
#  endif
     &                FILTER(ng) % filu3d(:,:,:,tindx))
      IF (InpThread) THEN
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(Vname(1,idUvel)), tindx, TRIM(ncname)
          exit_flag=2
          RETURN
        ELSE
          WRITE (stdout,100) TRIM(Vname(2,idUvel)), Fmin, Fmax
        END IF
      END IF
!
!  Read in 3D momentum component (m/s) in the ETA-direction.
!
      ncname=FILname(3,ng)
      status=nf90_inq_varid(ncFILid(3,ng), TRIM(Vname(1,idVvel)), varid)
      gtype=var_flag(varid)*v3dvar
      status=nf_fread3d(ng, iNLM, ncname, ncFILid(3,ng),                &
     &                Vname(1,idVvel), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#  ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#  endif
     &                FILTER(ng) % filv3d(:,:,:,tindx))
      IF (InpThread) THEN
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,30) TRIM(Vname(1,idVvel)), tindx, TRIM(ncname)
          exit_flag=2
          RETURN
        ELSE
          WRITE (stdout,100) TRIM(Vname(2,idVvel)), Fmin, Fmax
        END IF
      END IF
#  ifndef FILTRIM
!
!  Read in density anomaly.
!
      ncname=FILname(7,ng)
      status=nf90_inq_varid(ncFILid(7,ng), TRIM(Vname(1,idDano)), varid)
      gtype=var_flag(varid)*r3dvar
      status=nf_fread3d(ng, iNLM, ncname, ncFILid(7,ng),                &
     &                Vname(1,idDano), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filrho(:,:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idDano)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in tracer type variables.
!
      ncname=FILname(8,ng)
      DO itrc=1,NAT
        status=nf90_inq_varid(ncFILid(8,ng),                            &
     &                  TRIM(Vname(1,idTvar(itrc))), varid)
        gtype=var_flag(varid)*r3dvar
        status=nf_fread3d(ng, iNLM, ncname, ncFILid(8,ng),              &
     &                Vname(1,idTvar(itrc)), varid, tindx, gtype,       &
     &                  Vsize, LBi, UBi, LBj, UBj, 1, N(ng),            &
     &                  Fscl, Fmin, Fmax,                               &
#   ifdef MASKING
     &                  GRID(ng) % rmask,                               &
#   endif
     &                  FILTER(ng) % filt(:,:,:,itrc,tindx))
        IF (InpThread) THEN
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) TRIM(Vname(1,idTvar(itrc))), tindx,       &
     &                        TRIM(ncname)
            exit_flag=2
            RETURN
          ELSE
            WRITE (stdout,100) TRIM(Vname(2,idTvar(itrc))), Fmin, Fmax
          END IF
        END IF
      END DO
!
!  Read in vertical (omega) velocity.
!
      ncname=FILname(9,ng)
      status=nf90_inq_varid(ncFILid(9,ng), TRIM(Vname(1,idOvel)), varid)
      gtype=var_flag(varid)*w3dvar
      status=nf_fread3d(ng, iNLM, ncname, ncFILid(9,ng),                &
     &                Vname(1,idOvel), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filw3d(:,:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idOvel)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
#  ifdef LMD_SKPP
!
!  Read in depth of surface boundary layer.
!
      ncname=FILname(1,ng)
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idHsbl)), varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(1,ng),                &
     &                Vname(1,idHsbl), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filhsbl(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idHsbl)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
#  ifdef LMD_BKPP
!
!  Read in depth of bottom boundary layer.
!
      ncname=FILname(1,ng)
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idHbbl)), varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(1,ng),                &
     &                Vname(1,idHbbl), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filhbbl(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idHbbl)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_AKV
!
!  Read in vertical viscosity coefficient.
!
      ncname=FILname(6,ng)
      status=nf90_inq_varid(ncFILid(6,ng), TRIM(Vname(1,idVvis)), varid)
      gtype=var_flag(varid)*w3dvar
      status=nf_fread3d(ng, iNLM, ncname, ncFILid(6,ng),                &
     &                Vname(1,idVvis), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filAKv(:,:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idVvis)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_AKT
!
!  Read in vertical diffusion coefficient for potential temperature.
!
      ncname=FILname(4,ng)
      status=nf90_inq_varid(ncFILid(4,ng), TRIM(Vname(1,idTdif)), varid)
      gtype=var_flag(varid)*w3dvar
      status=nf_fread3d(ng, iNLM, ncname, ncFILid(4,ng),                &
     &                Vname(1,idTdif), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filAKt(:,:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idTdif)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_AKS
!
!  Read in vertical diffusion coefficient for salinity.
!
      ncname=FILname(5,ng)
      status=nf90_inq_varid(ncFILid(5,ng), TRIM(Vname(1,idSdif)), varid)
      gtype=var_flag(varid)*w3dvar
      status=nf_fread3d(ng, iNLM, ncname, ncFILid(5,ng),                &
     &                Vname(1,idSdif), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj, 0, N(ng),              &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filAKs(:,:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idSdif)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#  endif
#  ifdef AVERAGES_FLUXES
!
!  Read in surface net heat flux.
!
      ncname=FILname(1,ng)
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idTsur(itemp))),&
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(1,ng),                &
     &                Vname(1,idTsur(itemp)), varid, tindx, gtype,      &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filstf(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idTsur(itemp))), tindx,          &
     &                         TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in surface net salt flux.
!
      ncname=FILname(1,ng)
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idTsur(isalt))),&
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(1,ng),                &
     &                Vname(1,idTsur(isalt)), varid, tindx, gtype,      &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filswf(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idTsur(isalt))), tindx,          &
     &                         TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#   ifdef BULK_FLUXES
!
!  Read in latent heat flux.
!
      ncname=FILname(1,ng)
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idLhea)), varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(1,ng),                &
     &                Vname(1,idLhea), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % fillhf(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idLhea)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in sensible heat flux.
!
      ncname=FILname(1,ng)
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idShea)), varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(1,ng),                &
     &                Vname(1,idShea), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filshf(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idShea)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in longwave radiation flux.
!
      ncname=FILname(1,ng)
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idLrad)), varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(1,ng),                &
     &                Vname(1,idLrad), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % fillrf(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idLrad)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#   endif
#   ifdef SHORTWAVE
!
!  Read in shortwave radiation flux.
!
      ncname=FILname(1,ng)
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idSrad)), varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(1,ng),                &
     &                Vname(1,idSrad), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filsrf(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idSrad)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
#   endif
#  endif
# endif
# ifdef AVERAGES_FLUXES
!
!  Read in surface u-momentum stress.
!
      ncname=FILname(1,ng)
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idUsms)), varid)
      gtype=var_flag(varid)*u2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(1,ng),                &
     &                Vname(1,idUsms), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                FILTER(ng) % filsus(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idUsms)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in surface v-momentum stress.
!
      ncname=FILname(1,ng)
      status=nf90_inq_varid(ncFILid(1,ng), TRIM(Vname(1,idVsms)), varid)
      gtype=var_flag(varid)*v2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(1,ng),                &
     &                Vname(1,idVsms), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#   endif
     &                FILTER(ng) % filsvs(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idVsms)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif
# ifdef ICE_MODEL
!
!  Read in ice u-velocity
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idUice)),      &
     &                varid)
      gtype=var_flag(varid)*u2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idUice), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % umask,                                 &
#   endif
     &                FILTER(ng) % filui(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idUice)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in ice v-velocity
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idVice)),      &
     &                varid)
      gtype=var_flag(varid)*v2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idVice), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % vmask,                                 &
#   endif
     &                FILTER(ng) % filvi(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idVice)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in ice concentration
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idAice)),      &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idAice), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filai(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idAice)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in ice thickness
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idHice)),      &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idHice), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filhi(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idHice)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in snow thickness
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idHsno)),      &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idHsno), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filhsn(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idHsno)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in surface water
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idSfwat)),     &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idSfwat), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filsfwat(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idSfwat)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in surface ice temperature
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idTice)),      &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idTice), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filtis(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idTice)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in interior ice temperature
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idTimid)),     &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idTimid), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filti(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idTimid)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in ice age
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idAgeice)),    &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idAgeice), varid, tindx, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filAgeice(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idAgeice)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in ice age
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idAgeice)),    &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idAgeice), varid, tindx, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filAgeice(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idAgeice)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in under-ice salinity
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idS0mk)),      &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idS0mk), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % fils0mk(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idS0mk)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in under-ice temperature
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idT0mk)),      &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idT0mk), varid, tindx, gtype,             &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filt0mk(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idT0mk)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in frazil ice growth rate
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idWfr)),       &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idWfr), varid, tindx, gtype,              &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filWfr(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idWfr)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in ice growth/melt rate
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idWai)),       &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idWai), varid, tindx, gtype,              &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filWai(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idWai)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in ice growth/melt rate
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idWao)),       &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idWao), varid, tindx, gtype,              &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filWao(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idWao)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in ice growth/melt rate
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idWio)),       &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idWio), varid, tindx, gtype,              &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filWio(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idWio)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in ice melt runoff rate
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idWro)),       &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idWro), varid, tindx, gtype,              &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filWro(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idWro)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in sig11 of ice stress tensor
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idSig11)),     &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idSig11), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filsig11(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idSig11)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in sig12 of ice stress tensor
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idSig12)),     &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idSig12), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filsig12(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idSig12)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in sig22 of ice stress tensor
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idSig22)),     &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idSig22), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filsig22(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idSig22)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in ice-ocean mass flux
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idIomflx)),    &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idIomflx), varid, tindx, gtype,           &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filio_mflux(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idIomflx)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in ice-water tau
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idTauiw)),     &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idTauiw), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filutau_iw(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idTauiw)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
!
!  Read in ice-water chu
!
      ncname=FILname(10,ng)
      status=nf90_inq_varid(ncFILid(10,ng), TRIM(Vname(1,idChuiw)),     &
     &                    varid)
      gtype=var_flag(varid)*r2dvar
      status=nf_fread2d(ng, iNLM, ncname, ncFILid(10,ng),               &
     &                Vname(1,idChuiw), varid, tindx, gtype,            &
     &                Vsize, LBi, UBi, LBj, UBj,                        &
     &                Fscl, Fmin, Fmax,                                 &
#   ifdef MASKING
     &                GRID(ng) % rmask,                                 &
#   endif
     &                FILTER(ng) % filchu_iw(:,:,tindx))
      IF (InpThread.and.(status.ne.nf90_noerr)) THEN
        WRITE (stdout,30) TRIM(Vname(1,idChuiw)), tindx, TRIM(ncname)
        exit_flag=2
        RETURN
      END IF
# endif

# ifdef FOO
! Don't think I need this because these are just running sums, and aren't
! used dynamically
# if (defined EW_PERIODIC || defined NS_PERIODIC) && !defined DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set periodic boundary conditions.
!-----------------------------------------------------------------------
!
      CALL exchange_r2d_tile (ng, iNLM, 1, Lm(ng), 1, Mm(ng),           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        FILTER(ng) % filzeta(:,:,tindx))
      CALL exchange_u2d_tile (ng, iNLM, 1, Lm(ng), 1, Mm(ng),           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        FILTER(ng) % filu2d(:,:,tindx))
      CALL exchange_v2d_tile (ng, iNLM, 1, Lm(ng), 1, Mm(ng),           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        FILTER(ng) % filv2d(:,:,tindx))
#  ifdef SOLVE3D
      CALL exchange_u3d_tile (ng, iNLM, 1, Lm(ng), 1, Mm(ng),           &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        FILTER(ng) % filu3d(:,:,:,tindx))
      CALL exchange_v3d_tile (ng, iNLM, 1, Lm(ng), 1, Mm(ng),           &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        NghostPoints,                             &
     &                        FILTER(ng) % filv3d(:,:,:,tindx))

#   ifndef FILTRIM
      DO itrc=1,NAT
        CALL exchange_r3d_tile (ng, iNLM, 1, Lm(ng), 1, Mm(ng),         &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          FILTER(ng) % filt(:,:,:,itrc,tindx))
      END DO

      CALL exchange_w3d_tile (ng, iNLM, 1, Lm(ng), 1, Mm(ng),           &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        FILTER(ng) % filw3d(:,:,:,tindx))
#   endif
#   ifdef LMD_SKPP
      CALL exchange_r2d_tile (ng, iNLM, 1, Lm(ng), 1, Mm(ng),           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        FILTER(ng) % filhsbl(:,:,tindx))
#   endif
#   ifdef LMD_BKPP
      CALL exchange_r2d_tile (ng, iNLM, 1, Lm(ng), 1, Mm(ng),           &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        NghostPoints,                             &
     &                        FILTER(ng) % filhbbl(:,:,tindx))
#   endif
#   ifdef AVERAGES_AKV
      CALL exchange_w3d_tile (ng, iNLM, 1, Lm(ng), 1, Mm(ng),           &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        FILTER(ng) % filAKv(:,:,:,tindx))
#   endif
#   ifdef AVERAGES_AKT
      CALL exchange_w3d_tile (ng, iNLM, 1, Lm(ng), 1, Mm(ng),           &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        FILTER(ng) % filAKt(:,:,:,tindx))
#   endif
#   ifdef AVERAGES_AKS
      CALL exchange_w3d_tile (ng, iNLM, 1, Lm(ng), 1, Mm(ng),           &
     &                        LBi, UBi, LBj, UBj, 0, N(ng),             &
     &                        NghostPoints,                             &
     &                        FILTER(ng) % filAKs(:,:,:,tindx))
#   endif
#   ifdef AVERAGES_FLUXES
#    ifdef BULK_FLUXES
#    endif
#    ifdef SHORTWAVE
#    endif
#   endif
#  endif
#  ifdef AVERAGES_FLUXES
#  endif
# endif
# endif

!
!  Close input NetCDF file.
!
      DO ifile=1,nfile
        status=nf90_close(ncFILid(ifile,ng))
      END DO
      IF (Master) WRITE (stdout,80) tindx
!
  10  FORMAT (/,' GET_FILTER - unable to find model variable: ',a,      &
     &        /,15x,'in filter NetCDF file: ',a)
  20  FORMAT (/,' GET_FILTER - unable to open filter NetCDF file: ',a)
  30  FORMAT (/,' GET_FILTER - error while reading variable: ',a,2x,    &
     &        'at time record = ',i3,/,15x,'in filter NetCDF file: ',a)
  40  FORMAT (/,' GET_FILTER - requested filter time record = ',i3,/,   &
     &        15x,'not found in input NetCDF: ',a,/,                    &
     &        15x,'number of available records = ',i3)
  50  FORMAT (/,' GET_FILTER - error while reading attribute: ',a)
  60  FORMAT (/,' GET_FILTER - error while inquiring attributes for',   &
     &          ' variable: ',a)
  70  FORMAT (/,' GET_FILTER - cannot inquire about time variable in',  &
     &          ' filter NetCDF file: ',a)
  80  FORMAT (' GET_FILTER - Processed filter fields,',                 &
     &        t48,'filterLevel = ',i3)
  90  FORMAT (6x,'GET_FILTER - ',a,t64,'t = ',i12)
 100  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')
      RETURN
      END SUBROUTINE get_filter
#else
      SUBROUTINE get_filter
      RETURN
      END SUBROUTINE get_filter
#endif
