#include "cppdefs.h"
      SUBROUTINE wrt_his (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes requested model fields at requested levels      !
!  into history NetCDF file.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BBL_MODEL
      USE mod_bbl
#endif
#if defined BEST_NPZ || defined BIO_GOANPZ
      USE mod_biology
#endif
#ifdef ADJUST_BOUNDARY
      USE mod_boundary
#endif
#ifdef SOLVE3D
      USE mod_coupling
#endif
      USE mod_forces
      USE mod_grid
#ifdef ICE_MODEL
      USE mod_ice
#endif
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed
      USE mod_sediment
#endif
      USE mod_stepping
!
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
#ifdef ADJUST_BOUNDARY
      USE nf_fwrite2d_bry_mod, ONLY : nf_fwrite2d_bry
#endif
#ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# ifdef ADJUST_BOUNDARY
      USE nf_fwrite3d_bry_mod, ONLY : nf_fwrite3d_bry
# endif
      USE omega_mod, ONLY : scale_omega
#endif
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
#ifdef ADJUST_BOUNDARY
      integer :: LBij, UBij
#endif
      integer :: gfactor, gtype, status
#ifdef SOLVE3D
      integer :: i, itrc, j, k, tile
#endif
      real(r8) :: scale

#ifdef SOLVE3D
      real(r8), allocatable :: wrk(:,:,:)
#endif
!
      SourceFile='wrt_his.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
#ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
#endif
!
!-----------------------------------------------------------------------
!  Write out history fields.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#if defined WRITE_WATER && defined MASKING
      gfactor=-1
#else
      gfactor=1
#endif
!
!  Set time record index.
!
      tHISindx(ng)=tHISindx(ng)+1
      NrecHIS(ng)=NrecHIS(ng)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, HISname(ng),                      &
     &                      TRIM(Vname(idtime,ng)), time(ng:),          &
     &                      (/tHISindx(ng)/), (/1/),                    &
     &                      ncid = ncHISid(ng),                         &
     &                      varid = hisVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN

#if defined SEDIMENT && defined SED_MORPH
!
!  Write out time-dependent bathymetry (m)
!
      IF (Hout(idBath,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idbath,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     GRID(ng) % h,                                &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idbath)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef WET_DRY
!
!  Write out wet/dry mask at RHO-points.
!
      scale=1.0_r8
      gtype=gfactor*r2dvar
      status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idRwet,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % rmask_io,                           &
# endif
     &                   GRID(ng) % rmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idRwet)), tHISindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at U-points.
!
      scale=1.0_r8
      gtype=gfactor*u2dvar
      status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUwet,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % umask_io,                           &
# endif
     &                   GRID(ng) % umask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idUwet)), tHISindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out wet/dry mask at V-points.
!
      scale=1.0_r8
      gtype=gfactor*v2dvar
      status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVwet,ng),      &
     &                   tHISindx(ng), gtype,                           &
     &                   LBi, UBi, LBj, UBj, scale,                     &
# ifdef MASKING
     &                   GRID(ng) % vmask_io,                           &
# endif
     &                   GRID(ng) % vmask_wet,                          &
     &                   SetFillVal = .FALSE.)
      IF (status.ne.nf90_noerr) THEN
        IF (Master) THEN
          WRITE (stdout,10) TRIM(Vname(1,idVwet)), tHISindx(ng)
        END IF
        exit_flag=3
        ioerror=status
        RETURN
      END IF
#endif
!
!  Write out free-surface (m)
!
      IF (Hout(idFsur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idFsur,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#endif
#ifdef WET_DRY
     &                     OCEAN(ng) % zeta(:,:,KOUT),                  &
     &                     SetFillVal = .FALSE.)
#else
     &                     OCEAN(ng) % zeta(:,:,KOUT))
#endif
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idFsur)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idRzet,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     OCEAN(ng) % rzeta(:,:,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRzet)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#endif
      END IF
#ifdef ADJUST_BOUNDARY
!
!  Write out free-surface open boundaries.
!
      IF (ANY(Lobc(:,isFsur,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iNLM, HISname(ng), ncHISid(ng),     &
     &                          Vname(1,idSbry(isFsur)),                &
     &                          hisVid(idSbry(isFsur),ng),              &
     &                          tHISindx(ng), r2dvar,                   &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % zeta_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isFsur))),            &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
!
!  Write out 2D U-momentum component (m/s).
!
      IF (Hout(idUbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbar,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
#endif
     &                     OCEAN(ng) % ubar(:,:,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbar)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idRu2d,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
#  endif
     &                     OCEAN(ng) % rubar(:,:,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRu2d)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idRuct,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
#   endif
     &                     COUPLING(ng) % rufrc)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRuct)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUfx1,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
#  endif
     &                     COUPLING(ng) % DU_avg1)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx1)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUfx2,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
#  endif
     &                     COUPLING(ng) % DU_avg2)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUfx2)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
#endif
      END IF
#ifdef ADJUST_BOUNDARY
!
!  Write out 2D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUbar,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iNLM, HISname(ng), ncHISid(ng),     &
     &                          Vname(1,idSbry(isUbar)),                &
     &                          hisVid(idSbry(isUbar),ng),              &
     &                          tHISindx(ng), u2dvar,                   &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % ubar_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUbar))),            &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
!
!  Write out 2D V-momentum component (m/s).
!
      IF (Hout(idVbar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbar,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
#endif
     &                     OCEAN(ng) % vbar(:,:,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbar)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#ifdef FORWARD_WRITE
# ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idRv2d,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
#  endif
     &                     OCEAN(ng) % rvbar(:,:,KOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRv2d)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idRvct,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
#   endif
     &                     COUPLING(ng) % rvfrc)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRvct)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVfx1,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
#  endif
     &                     COUPLING(ng) % DV_avg1)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx1)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVfx2,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
#  endif
     &                     COUPLING(ng) % DV_avg2)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVfx2)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
#endif
      END IF
#ifdef ADJUST_BOUNDARY
!
!  Write out 2D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVbar,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite2d_bry (ng, iNLM, HISname(ng), ncHISid(ng),     &
     &                          Vname(1,idSbry(isVbar)),                &
     &                          hisVid(idSbry(isVbar),ng),              &
     &                          tHISindx(ng), v2dvar,                   &
     &                          LBij, UBij, Nbrec(ng), scale,           &
     &                          BOUNDARY(ng) % vbar_obc(LBij:,:,:,      &
     &                                                  Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVbar))),            &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef SOLVE3D
!
!  Write out 3D U-momentum component (m/s).
!
      IF (Hout(idUvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idUvel,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
# endif
     &                     OCEAN(ng) % u(:,:,:,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUvel)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idRu3d,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
#  endif
     &                     OCEAN(ng) % ru(:,:,:,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRu3d)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
      END IF
# ifdef ADJUST_BOUNDARY
!
!  Write out 3D U-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isUvel,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite3d_bry (ng, iNLM, HISname(ng), ncHISid(ng),     &
     &                          Vname(1,idSbry(isUvel)),                &
     &                          hisVid(idSbry(isUvel),ng),              &
     &                          tHISindx(ng), u3dvar,                   &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % u_obc(LBij:,:,:,:,       &
     &                                               Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isUvel))),            &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 3D V-momentum component (m/s).
!
      IF (Hout(idVvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idVvel,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
# endif
     &                     OCEAN(ng) % v(:,:,:,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvel)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# if defined FORWARD_WRITE && defined FORWARD_RHS
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idRv3d,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
#  endif
     &                     OCEAN(ng) % rv(:,:,:,NOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idRv3d)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
# endif
      END IF
# ifdef ADJUST_BOUNDARY
!
!  Write out 3D V-momentum component open boundaries.
!
      IF (ANY(Lobc(:,isVvel,ng))) THEN
        scale=1.0_r8
        status=nf_fwrite3d_bry (ng, iNLM, HISname(ng), ncHISid(ng),     &
     &                          Vname(1,idSbry(isVvel)),                &
     &                          hisVid(idSbry(isVvel),ng),              &
     &                          tHISindx(ng), v3dvar,                   &
     &                          LBij, UBij, 1, N(ng), Nbrec(ng), scale, &
     &                          BOUNDARY(ng) % v_obc(LBij:,:,:,:,       &
     &                                               Lbout(ng)))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSbry(isVvel))),            &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      IF (Hout(idOvel,ng)) THEN
        IF (.not.allocated(wrk)) THEN
          allocate (wrk(LBi:UBi,LBj:UBj,0:N(ng)))
          wrk(LBi:UBi,LBj:UBj,0:N(ng))=0.0_r8
        END IF
        scale=1.0_r8
        gtype=gfactor*w3dvar
        DO tile=0,NtileX(ng)*NtileE(ng)-1
          CALL scale_omega (ng, TILE, LBi, UBi, LBj, UBj, 0, N(ng),     &
     &                      GRID(ng) % pm,                              &
     &                      GRID(ng) % pn,                              &
     &                      OCEAN(ng) % W,                              &
     &                      wrk)
        END DO
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idOvel,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     wrk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOvel)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
        deallocate (wrk)
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Hout(idWvel,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idWvel,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     OCEAN(ng) % wvel)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWvel)), tHISindx(ng)
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
        IF (Hout(idTvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisTid(itrc,ng),    &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
# ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
# endif
     &                       OCEAN(ng) % t(:,:,:,NOUT,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTvar(itrc))),            &
     &                          tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

# ifdef BEST_NPZ
#  ifdef STATIONARY
!--------------------------------------
!  Write out 3D stationary variable array
!--------------------------------------
      DO itrc=1,NTS(ng)
        IF (Hout(idTSvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
!           print *,'st=',OCEAN(ng) % st(:,:,:,NOUT,itrc)
!          if (hisTSid(itrc,ng) .lt. 0) hisTSid(itrc,ng) = 0
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisTSid(itrc,ng),   &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
#   endif
     &                       OCEAN(ng) % st(:,:,:,NOUT,itrc))
          !---------------------
          !Zero the output array
          !---------------------

!           print *,'st=',OCEAN(ng) % st(:,:,:,1,1)
!           print *,'NOUT=',NOUT
!           print *,'itrc=',itrc

          OCEAN(ng) % st(:,:,:,NOUT,itrc) = 0.0_r8
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN

              WRITE (stdout,10) TRIM(Vname(1,idTSvar(itrc))),           &
     &                          tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif

#  ifdef STATIONARY2
!--------------------------------------
!  Write out 2D stationary variable array
!--------------------------------------
      DO itrc=1,NTS2(ng)
        IF (Hout(idTS2var(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          if (hisTS2id(itrc,ng) .lt. 0) hisTS2id(itrc,ng) = 0
          status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisTS2id(itrc,ng),  &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj,  scale,                &
#   ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
#   endif
     &                       OCEAN(ng) % st2(:,:,NOUT,itrc))
          !---------------------
          !Zero the output array
          !---------------------
          OCEAN(ng) % st2(:,:,NOUT,itrc) = 0.0_r8
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTS2var(itrc))),          &
     &                          tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif

!--------------------------------------
!  Write out 2D production array
!--------------------------------------

#  ifdef PROD2
      DO itrc=1,NPT2(ng)
        IF (Hout(idPT2var(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          if (hisPT2id(itrc,ng) .lt. 0) hisPT2id(itrc,ng) = 0
          status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisPT2id(itrc,ng),  &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
#   endif
     &                       OCEAN(ng) % pt2(:,:,NOUT,itrc))
          !---------------------
          !Zero the output array
          !---------------------
          OCEAN(ng) % pt2(:,:,NOUT,itrc) = 0.0_r8
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idPT2var(itrc))),          &
     &                          tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!--------------------------------------
!  Write out 3D production array
!--------------------------------------

#  ifdef PROD3
      DO itrc=1,NPT3(ng)
        IF (Hout(idPT3var(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          if (hisPT3id(itrc,ng) .lt. 0) hisPT3id(itrc,ng) = 0
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisPT3id(itrc,ng),  &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#   ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
#   endif
     &                       OCEAN(ng) % pt3(:,:,:,NOUT,itrc))
          !---------------------
          !Zero the output array
          !---------------------
          OCEAN(ng) % pt3(:,:,:,NOUT,itrc) = 0.0_r8
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idPT3var(itrc))),          &
     &                          tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!------------------------------
!Write out benthic variables
!-----------------------------

#  if defined BENTHIC
              DO itrc=1,NBeT(ng)
              IF (Hout(idBvar(itrc),ng)) THEN
                 scale=1.0_r8
                 gtype=gfactor*r2dvar
             status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisBid(itrc,ng), &
     &                          tHISindx(ng), gtype,                    &
     &                          LBi, UBi, LBj, UBj, scale,              &

!Will need to switch this if have more than one depth level for benthos
!      status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisBid(itrc,ng),       &
!     &                          tHISindx(ng), gtype,                   &
!     &                          LBi, UBi, LBj, UBj, 1, NBL(ng), scale, &

#   ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
#   endif
     &                       OCEAN(ng) % bt(:,:,1,NOUT,itrc))
!             print *,OCEAN(ng) % bt(LBi,LBj,1,NOUT,itrc)
!             print *, 'LBi=',LBi,'LBj=',LBj,'NOUT=',NOUT,'itrc=',itrc

           IF (status.ne.nf90_noerr) THEN
           IF (Master) THEN
            WRITE (stdout,10)TRIM(Vname(1,idBvar(itrc))),               &
     &                       tHISindx(ng)
           END IF
              exit_flag=3
              ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif

!------------------------------
!Write out ice bio variables
!-----------------------------

#  if defined ICE_BIO
#   ifdef CLIM_ICE_1D
      DO itrc=1,NIceT(ng)
        IF (Hout(idIceBvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisIceBid(itrc,ng), &
     &                          tHISindx(ng), gtype,                    &
     &                          LBi, UBi, LBj, UBj, scale,              &

!Will need to switch this if have more than one depth level for benthos
!      status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisIceBid(itrc,ng), &
!     &                          tHISindx(ng), gtype,                    &
!     &                          LBi, UBi, LBj, UBj, 1, NBL(ng), scale,  &

#   ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
#   endif
     &                       OCEAN(ng) % it(:,:,NOUT,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10)TRIM(Vname(1,idIceBvar(itrc))),          &
     &                       tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#   else
!
!  Write out Ice Algae
!
      IF (Hout(idIcePhL,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idIcePhL,ng),                         &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % IcePhL(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIcePhL)),                  &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

      !
!  Write out Ice Nitrate
!
      IF (Hout(idIceNO3,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idIceNO3,ng),                         &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % IceNO3(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIceNO3)),                  &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

!  Write out Ice Ammonium
!
      IF (Hout(idIceNH4,ng)) THEN
!      IF (Hout(idIceBvar,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idIceNH4,ng),                         &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % IceNH4(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIceNH4)),                  &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

!  Write out IceLog (doesn't work with IceLog being int)
!
#    ifdef FOO
      IF (Hout(idIceLog,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idIceLog,ng),                         &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % IceLog(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIceLog)),                  &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#    endif
#   endif
#  endif
# endif

# ifdef BIO_GOANPZ
!--------------------------------------
!  Write out stationary tracer variable
!--------------------------------------
      DO itrc=1,NTS(ng)
        IF (Hout(idTSvar(itrc),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r3dvar
          if (hisTSid(itrc,ng) .lt. 0) hisTSid(itrc,ng) = 0
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisTSid(itrc,ng),   &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
#  endif
     &                       OCEAN(ng) % st(:,:,1,NOUT,itrc))
          !---------------------
          !Zero the output array
          !---------------------
          OCEAN(ng) % st(:,:,:,NOUT,itrc) = 0.0_r8
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTSvar(itrc))),           &
     &                          tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
!----------------------------
# ifdef ADJUST_BOUNDARY
!
!  Write out 3D tracers open boundaries.
!
      DO itrc=1,NT(ng)
        IF (ANY(Lobc(:,isTvar(itrc),ng))) THEN
          scale=1.0_r8
          status=nf_fwrite3d_bry (ng, iNLM, HISname(ng), ncHISid(ng),   &
     &                            Vname(1,idSbry(isTvar(itrc))),        &
     &                            hisVid(idSbry(isTvar(itrc)),ng),      &
     &                            tHISindx(ng), r3dvar,                 &
     &                            LBij, UBij, 1, N(ng), Nbrec(ng),      &
     &                            scale,                                &
     &                            BOUNDARY(ng) % t_obc(LBij:,:,:,:,     &
     &                                                 Lbout(ng),itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbry(isTvar(itrc)))),    &
     &                          tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
!----------------------------
!  Write out density anomaly.
!----------------------------
      IF (Hout(idDano,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idDano,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     OCEAN(ng) % rho)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDano)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef NEMURO_SED1
!
!  Write out PON in sediment.
!
      IF (Hout(idPONsed,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idPONsed,ng),  &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     OCEAN(ng) % PONsed)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPONsed)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out OPAL in sediment.
!
      IF (Hout(idOPALsed,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idOPALsed,ng), &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     OCEAN(ng) % OPALsed)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOPALsed)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out DENIT in sediment.
!
      IF (Hout(idDENITsed,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idDENITsed,ng),&
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     OCEAN(ng) % DENITsed)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idDENITsed)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out buried PON in sediment.
!
      IF (Hout(idPONbur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idPONbur,ng),  &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     OCEAN(ng) % PON_burial)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPONbur)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out buried OPAL in sediment.
!
      IF (Hout(idOPALbur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idOPALbur,ng), &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     OCEAN(ng) % OPAL_burial)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idOPALbur)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef LMD_SKPP
!
!  Write out depth surface boundary layer.
!
      IF (Hout(idHsbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idHsbl,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     MIXING(ng) % hsbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsbl)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef LMD_BKPP
!
!  Write out depth surface boundary layer.
!
      IF (Hout(idHbbl,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idHbbl,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     MIXING(ng) % hbbl)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHbbl)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef SSSFLX
!
!  Write out sea surface salinity correction salt flux
!
      IF (Hout(idSSSf,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idSSSf,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     FORCES(ng) % sssflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSSSf)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif

# if defined FORWARD_WRITE && defined LMD_NONLOCAL
!
!  Write out KPP nonlocal transport.
!
      DO i=1,NAT
        IF (Hout(idGhat(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*w3dvar
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idGhat(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 0, N(ng), scale,       &
#  ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
#  endif
     &                       MIXING(ng) % ghats(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idGhat(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
!
!  Write out vertical viscosity coefficient.
!
      IF (Hout(idVvis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idVvis,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     MIXING(ng) % Akv,                            &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVvis)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Hout(idTdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idTdif,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     MIXING(ng) % Akt(:,:,:,itemp),               &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTdif)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      IF (Hout(idSdif,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idSdif,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     MIXING(ng) % Akt(:,:,:,isalt),               &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSdif)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
!
!  Write out turbulent kinetic energy.
!
      IF (Hout(idMtke,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idMtke,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     MIXING(ng) % tke(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtke)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  ifdef FORWARD_WRITE
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idVmKK,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     MIXING(ng) % Akk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmKK)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF
!
!  Write out turbulent length scale field.
!
      IF (Hout(idMtls,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idMtls,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     MIXING(ng) % gls(:,:,:,NOUT),                &
     &                     SetFillVal = .FALSE.)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idMtls)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  ifdef FORWARD_WRITE
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idVmLS,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     MIXING(ng) % Lscale)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmLS)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
#  if defined FORWARD_WRITE && defined GLS_MIXING
        scale=1.0_r8
        gtype=gfactor*w3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idVmKP,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 0, N(ng), scale,         &
#   ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#   endif
     &                     MIXING(ng) % Akp)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVmKP)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
#  endif
      END IF
# endif
# ifdef ICE_MODEL
!
!  Write out ice 2D momentum component (m/s) in the XI-direction.
!
      IF (Hout(idUice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idUice,ng),                           &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
#  endif
     &                     ICE(ng) % ui(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUice)),                    &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice 2D momentum component (m/s) in the ETA-direction.
!
      IF (Hout(idVice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idVice,ng),                           &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
#  endif
     &                     ICE(ng) % vi(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVice)),                    &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice concentration
!
      IF (Hout(idAice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idAice,ng),                           &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % ai(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAice)),                    &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice average thickness
!
      IF (Hout(idHice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idHice,ng),                           &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % hi(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHice)),                    &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out snow average thickness
!
      IF (Hout(idHsno,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idHsno,ng),                           &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % hsn(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idHsno)),                    &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface water thickness (on ice)
!
      IF (Hout(idSfwat,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idSfwat,ng),                          &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % sfwat(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSfwat)),                   &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice age.
!
      IF (Hout(idAgeice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idAgeice,ng),                         &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % ageice(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idAgeice)),                  &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice-ocean mass flux
!
      IF (Hout(idIomflx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idIomflx,ng),                         &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % io_mflux)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idIomflx)),                  &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice/snow surface temperature
!
      IF (Hout(idTice,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idTice,ng),                           &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % tis)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTice)),                    &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice interior temperature
!
      IF (Hout(idTimid,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idTimid,ng),                          &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % ti(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTimid)),                   &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out internal ice stress component 11
!
      IF (Hout(idSig11,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idSig11,ng),                          &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % sig11(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig11)),                   &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out internal ice stress component 12
!
      IF (Hout(idSig12,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idSig12,ng),                          &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % sig12(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig12)),                   &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out internal ice stress component 22
!
      IF (Hout(idSig22,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idSig22,ng),                          &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % sig22(:,:,IUOUT))
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSig22)),                   &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice-ocean friction velocity
!
      IF (Hout(idTauiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idTauiw,ng),                          &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % utau_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTauiw)),                   &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice-ocean momentum transfer coefficient
!
      IF (Hout(idChuiw,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idChuiw,ng),                          &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % chu_iw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idChuiw)),                   &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out temperature of molecular sublayer under ice
!
      IF (Hout(idT0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idT0mk,ng),                           &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % t0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idT0mk)),                    &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out salinity of molecular sublayer under ice
!
      IF (Hout(idS0mk,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idS0mk,ng),                           &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % s0mk)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idS0mk)),                    &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice freeze Wfr
!
      IF (Hout(idWfr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idWfr,ng),                            &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % wfr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWfr)),                     &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice melt/freeze wai
!
      IF (Hout(idWai,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idWai,ng),                            &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % wai)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWai)),                     &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice melt/freeze Wao
!
      IF (Hout(idWao,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idWao,ng),                            &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % wao)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWao)),                     &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice melt/freeze wio
!
      IF (Hout(idWio,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idWio,ng),                            &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % wio)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWio)),                     &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out ice melt/freeze wro
!
      IF (Hout(idWro,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                       &
     &                     hisVid(idWro,ng),                            &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     ICE(ng) % wro)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWro)),                     &
     &                        tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM || defined ATM_PRESS
!
!  Write out surface air pressure.
!
      IF (Hout(idPair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idPair,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     FORCES(ng) % Pair)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idPair)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# if defined BULK_FLUXES || defined ECOSIM
!
!  Write out surface winds.
!
      IF (Hout(idUair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUair,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     FORCES(ng) % Uwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUair)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

      IF (Hout(idVair,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVair,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     FORCES(ng) % Vwind)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVair)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out surface active traces fluxes.
!
      DO itrc=1,NAT
        IF (Hout(idTsur(itrc),ng)) THEN
          IF (itrc.eq.itemp) THEN
# ifdef SO_SEMI
            scale=1.0_r8
# else
            scale=rho0*Cp                   ! Celsius m/s to W/m2
# endif
          ELSE IF (itrc.eq.isalt) THEN
            scale=1.0_r8
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idTsur(itrc),ng),                   &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
# ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
# endif
     &                       FORCES(ng) % stflx(:,:,itrc))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idTsur(itrc))),            &
     &                          tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

# ifdef BULK_FLUXES
!
!  Write out latent heat flux.
!
      IF (Hout(idLhea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idLhea,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     FORCES(ng) % lhflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLhea)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out sensible heat flux.
!
      IF (Hout(idShea,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idShea,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     FORCES(ng) % shflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idShea)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out longwave radiation flux.
!
      IF (Hout(idLrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idLrad,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     FORCES(ng) % lrflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idLrad)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  ifdef EMINUSP
!
!  Write out E-P (m/s).
!
      IF (Hout(idEmPf,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idEmPf,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#   endif
     &                     FORCES(ng) % EminusP)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idEmPf)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out evaporation rate (kg/m2/s).
!
      IF (Hout(idevap,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idevap,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#   endif
     &                     FORCES(ng) % evap)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idevap)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out precipitation rate (kg/m2/s).
!
      IF (Hout(idrain,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idrain,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#   ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#   endif
     &                     FORCES(ng) % rain)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idrain)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#  endif
# endif
# ifdef SHORTWAVE
!
!  Write out shortwave radiation flux.
!
      IF (Hout(idSrad,ng)) THEN
        scale=rho0*Cp
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idSrad,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     FORCES(ng) % srflx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idSrad)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
#endif
!
!  Write out surface U-momentum stress.
!
      IF (Hout(idUsms,ng)) THEN
#ifdef SO_SEMI
        scale=1.0_r8
#else
        scale=rho0                          ! m2/s2 to Pa
#endif
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUsms,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
#endif
     &                     FORCES(ng) % sustr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUsms)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Hout(idVsms,ng)) THEN
#ifdef SO_SEMI
        scale=1.0_r8
#else
        scale=rho0
#endif
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVsms,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
#endif
     &                     FORCES(ng) % svstr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVsms)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Hout(idUbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbms,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
#endif
     &                     FORCES(ng) % bustr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbms)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Hout(idVbms,ng)) THEN
        scale=-rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbms,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
#endif
     &                     FORCES(ng) % bvstr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbms)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#ifdef SOLVE3D
# ifdef BBL_MODEL
!
!  Write out current-induced, bottom U-stress at RHO-points.
!
      IF (Hout(idUbrs,ng)) THEN
        scale=-rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbrs,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     BBL(ng) % bustrc)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbrs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out current-induced, bottom V-stress at RHO-points.
!
      IF (Hout(idVbrs,ng)) THEN
        scale=-rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbrs,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     BBL(ng) % bvstrc)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbrs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bottom U-stress at RHO-points.
!
      IF (Hout(idUbws,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbws,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     BBL(ng) % bustrw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbws)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bottom V-stress at RHO-points.
!
      IF (Hout(idVbws,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbws,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     BBL(ng) % bvstrw)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbws)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wind and current, bottom U-stress at RHO-points.
!
      IF (Hout(idUbcs,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbcs,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     BBL(ng) % bustrcwmax)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbcs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out maximum wind and current, bottom V-stress at RHO-points.
!
      IF (Hout(idVbcs,ng)) THEN
        scale=rho0
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbcs,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     BBL(ng) % bvstrcwmax)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbcs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity at RHO-points.
!
      IF (Hout(idUbot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbot,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     BBL(ng) % Ubot)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbot)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity at RHO-points
!
      IF (Hout(idVbot,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbot,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     BBL(ng) % Vbot)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbot)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom U-velocity above bed at RHO-points.
!
      IF (Hout(idUbur,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idUbur,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     BBL(ng) % Ur)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idUbur)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out bottom V-velocity above bed at RHO-points.
!
      IF (Hout(idVbvr,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idVbvr,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     BBL(ng) % Vr)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idVbvr)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
# ifdef SEDIMENT
#  ifdef BEDLOAD
!
!  Write out bed load transport in U-direction.
!
      DO i=1,NST
        IF (Hout(idUbld(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*u2dvar
          status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idUbld(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % umask_io,                       &
#   endif
     &                       SEDBED(ng) % bedldu(:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idUbld(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
!
!  Write out bed load transport in V-direction.
!
        IF (Hout(idVbld(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*v2dvar
          status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idVbld(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#   ifdef MASKING
     &                       GRID(ng) % vmask_io,                       &
#   endif
     &                       SEDBED(ng) % bedldv(:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idVbld(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
#  endif
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        IF (Hout(idfrac(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*b3dvar
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idfrac(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, Nbed, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
#  endif
     &                       SEDBED(ng) % bed_frac(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idfrac(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment mass of each size class in each bed layer.
!
      DO i=1,NST
        IF (Hout(idBmas(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*b3dvar
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idBmas(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, Nbed, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
#  endif
     &                       SEDBED(ng) % bed_mass(:,:,:,NOUT,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBmas(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (Hout(idSbed(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*b3dvar
          status=nf_fwrite3d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idSbed(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, Nbed, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
#  endif
     &                       SEDBED(ng) % bed(:,:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idSbed(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBOTP
        IF (Hout(idBott(i),ng)) THEN
          IF (i.eq.itauc) THEN
            scale=rho0
          ELSE
            scale=1.0_r8
          END IF
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, iNLM, ncHISid(ng),                     &
     &                       hisVid(idBott(i),ng),                      &
     &                       tHISindx(ng), gtype,                       &
     &                       LBi, UBi, LBj, UBj, scale,                 &
#  ifdef MASKING
     &                       GRID(ng) % rmask_io,                       &
#  endif
     &                       SEDBED(ng) % bottom(:,:,i))
          IF (status.ne.nf90_noerr) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idBott(i))), tHISindx(ng)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
# endif
#endif
#ifdef NEARSHORE_MELLOR
!
!  Write out 2D radiation stress, Sxx-component.
!
      IF (Hout(idW2xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idW2xx,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     MIXING(ng) % Sxx_bar)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xx)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Sxy-component.
!
      IF (Hout(idW2xy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idW2xy,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     MIXING(ng) % Sxy_bar)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2xy)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D radiation stress, Syy-component.
!
      IF (Hout(idW2yy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idW2yy,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     MIXING(ng) % Syy_bar)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW2yy)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out total 2D U-radiation stress.
!
      IF (Hout(idU2rs,ng)) THEN
        scale=rho0
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idU2rs,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
# endif
     &                     MIXING(ng) % rustr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2rs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out total 2D V-radiation stress.
!
      IF (Hout(idV2rs,ng)) THEN
        scale=rho0
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idV2rs,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
# endif
     &                     MIXING(ng) % rvstr2d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2rs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D U-momentum Stokes drift velocity.
!
      IF (Hout(idU2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idU2sd,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
# endif
     &                     OCEAN(ng) % ubar_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU2Sd)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 2D V-momentum Stokes drift velocity.
!
      IF (Hout(idV2Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idV2Sd,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
# endif
     &                     OCEAN(ng) % vbar_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV2Sd)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# ifdef SOLVE3D
!
!  Write out 3D radiation stress, Sxx-horizontal component.
!
      IF (Hout(idW3xx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idW3xx,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     MIXING(ng) % Sxx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3xx)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Sxy-horizontal component.
!
      IF (Hout(idW3xy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idW3xy,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     MIXING(ng) % Sxy)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3xy)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Syy-horizontal component.
!
      IF (Hout(idW3yy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idW3yy,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     MIXING(ng) % Syy)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3yy)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# ifdef NEARSHORE_MELLOR05
!
!  Write out 3D radiation stress, Szx-vertical component.
!
      IF (Hout(idW3zx,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idW3zx,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     MIXING(ng) % Szx)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3zx)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D radiation stress, Szy-vertical component.
!
      IF (Hout(idW3zy,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng), hisVid(idW3zy,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
#  endif
     &                     MIXING(ng) % Szy)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idW3zy)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out 3D total U-radiation stress.
!
      IF (Hout(idU3rs,ng)) THEN
        scale=rho0
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng),hisVid(idU3rs,ng),     &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
#  endif
     &                     MIXING(ng) % rustr3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU3rs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D total V-radiation stress.
!
      IF (Hout(idV3rs,ng)) THEN
        scale=rho0
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng),hisVid(idV3rs,ng),     &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
#  endif
     &                     MIXING(ng) % rvstr3d)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV3rs)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D U-momentum Stokes drift velocity.
!
      IF (Hout(idU3Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*u3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng),hisVid(idU3Sd,ng),     &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % umask_io,                         &
#  endif
     &                     OCEAN(ng) % u_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idU3Sd)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out 3D V-momentum stokes velocity.
!
      IF (Hout(idV3Sd,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*v3dvar
        status=nf_fwrite3d(ng, iNLM, ncHISid(ng),hisVid(idV3Sd,ng),     &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, 1, N(ng), scale,         &
#  ifdef MASKING
     &                     GRID(ng) % vmask_io,                         &
#  endif
     &                     OCEAN(ng) % v_stokes)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idV3Sd)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
#endif
#ifdef WAVES_HEIGHT
!
!  Write out wind-induced wave height.
!
      IF (Hout(idWamp,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWamp,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     FORCES(ng) % Hwave)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWamp)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef WAVES_LENGTH
!
!  Write out wind-induced wave length.
!
      IF (Hout(idWlen,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWlen,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     FORCES(ng) % Lwave)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWlen)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef WAVES_DIR
!
!  Write out wind-induced wave direction.
!
      IF (Hout(idWdir,ng)) THEN
        scale=rad2deg
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWdir,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     FORCES(ng) % Dwave)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdir)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef WAVES_TOP_PERIOD
!
!  Write out wind-induced surface wave period.
!
      IF (Hout(idWptp,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWptp,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     FORCES(ng) % Pwave_top)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWptp)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef WAVES_BOT_PERIOD
!
!  Write out wind-induced bottom wave period.
!
      IF (Hout(idWpbt,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWpbt,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     FORCES(ng) % Pwave_bot)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWpbt)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#ifdef WAVES_UB
!
!  Write out wind-induced wave bottom orbital velocity.
!
      IF (Hout(idWorb,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWorb,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     FORCES(ng) % Ub_swan)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWorb)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
#if defined TKE_WAVEDISS || defined WAVES_OCEAN
!
!  Write out wave dissipation.
!
      IF (Hout(idWdis,ng)) THEN
        scale=1.0_r8
        gtype=gfactor*r2dvar
        status=nf_fwrite2d(ng, iNLM, ncHISid(ng), hisVid(idWdis,ng),    &
     &                     tHISindx(ng), gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
# ifdef MASKING
     &                     GRID(ng) % rmask_io,                         &
# endif
     &                     FORCES(ng) % wave_dissip)
        IF (status.ne.nf90_noerr) THEN
          IF (Master) THEN
            WRITE (stdout,10) TRIM(Vname(1,idWdis)), tHISindx(ng)
          END IF
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Synchronize history NetCDF file to disk to allow other processes
!  to access data immediately after it is written.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, HISname(ng), ncHISid(ng))
      IF (exit_flag.ne.NoError) RETURN

#ifdef SOLVE3D
      IF (Master) WRITE (stdout,20) KOUT, NOUT, tHISindx(ng)
#ifdef OFFLINE_FLOATS
      IF (Master) print *, '  WRT_HIS model time (days), t=', tdays
#endif
#else
      IF (Master) WRITE (stdout,20) KOUT, tHISindx(ng)
#endif
!
  10  FORMAT (/,' WRT_HIS - error while writing variable: ',a,/,11x,    &
     &        'into history NetCDF file for time record: ',i4)
#ifdef SOLVE3D
  20  FORMAT (6x,'WRT_HIS   - wrote history  fields (Index=', i1,       &
     &        ',',i1,') into time record = ',i7.7)
#else
  20  FORMAT (6x,'WRT_HIS   - wrote history  fields (Index=', i1,       &
     &        ') into time record = ',i7.7)
#endif
      RETURN
      END SUBROUTINE wrt_his
