#include "cppdefs.h"
#if defined NEMURO_SAN && defined FLEET
      SUBROUTINE wrt_fleet (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine writes simulated drifter trajectories into fleet    !
!  NetCDF file.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fleet
      USE mod_fish
      USE mod_biology
      USE mod_grid
      USE mod_ocean
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_stepping
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
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
      integer :: LBi, UBi, LBj, UBj

      real(r8), dimension(Nboats(ng)) :: Tout
!
      SourceFile='wrt_fleet.F'
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out fleet data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set time record index.
!
      tFLEETindx(ng)=tFLEETindx(ng)+1
      NrecFLEET(ng)=NrecFLEET(ng)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, FLEETname(ng),                    &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/tFLEETindx(ng)/), (/1/),                  &
     &                      ncid = ncFLEETid(ng),                       &
     &                      varid = fleetVid(idtime,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet X-grid locations.
!
      DO l=1,Nboats(ng)
        Tout(l)=REAL(BOATS(ng)%boat(ibiloc,l),r8)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEETname(ng),                    &
     &                      TRIM(Vname(1,idBxprt)), Tout,               &
     &                      (/1,tFLEETindx(ng)/), (/Nboats(ng),1/),     &
     &                      ncid = ncFLEETid(ng),                       &
     &                      varid = fleetVid(idBxprt,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet Y-grid locations.
!
      DO l=1,Nboats(ng)
        Tout(l)=REAL(BOATS(ng)%boat(ibjloc,l),r8)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEETname(ng),                    &
     &                      TRIM(Vname(1,idByprt)), Tout,               &
     &                      (/1,tFLEETindx(ng)/), (/Nboats(ng),1/),     &
     &                      ncid = ncFLEETid(ng),                       &
     &                      varid = fleetVid(idByprt,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet port locations.
!
      DO l=1,Nboats(ng)
        Tout(l)=REAL(BOATS(ng)%boat(ibport,l),r8)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEETname(ng),                    &
     &                      TRIM(Vname(1,idBnprt)), Tout,               &
     &                      (/1,tFLEETindx(ng)/), (/Nboats(ng),1/),     &
     &                      ncid = ncFLEETid(ng),                       &
     &                      varid = fleetVid(idBnprt,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet number of fished locations.
!
      DO l=1,Nboats(ng)
        Tout(l)=REAL(BOATS(ng)%boat(ibnloc,l),r8)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEETname(ng),                    &
     &                      TRIM(Vname(1,idBnfloc)), Tout,              &
     &                      (/1,tFLEETindx(ng)/), (/Nboats(ng),1/),     &
     &                      ncid = ncFLEETid(ng),                       &
     &                      varid = fleetVid(idBnfloc,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet daily catch.
!
      DO l=1,Nboats(ng)
        Tout(l)=BOATS(ng)%catch(l)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEETname(ng),                    &
     &                      TRIM(Vname(1,idBcatch)), Tout,              &
     &                      (/1,tFLEETindx(ng)/), (/Nboats(ng),1/),     &
     &                      ncid = ncFLEETid(ng),                       &
     &                      varid = fleetVid(idBcatch,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet cummulative annual catch.
!
      DO l=1,Nboats(ng)
        Tout(l)=BOATS(ng)%cumcatch(l)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEETname(ng),                    &
     &                      TRIM(Vname(1,idBcmctch)), Tout,             &
     &                      (/1,tFLEETindx(ng)/), (/Nboats(ng),1/),     &
     &                      ncid = ncFLEETid(ng),                       &
     &                      varid = fleetVid(idBcmctch,ng))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out average catch per unit effort (kg)
!
      status=nf_fwrite2d(ng, iNLM, ncFLEETid(ng), fleetVid(idBacpue,ng),&
     &                   tFLEETindx(ng), r2dvar,                        &
     &                   LBi, UBi, LBj, UBj, 1.0_r8,                    &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   OCEAN(ng) % avgcpue(:,:))
      IF (status.ne.nf90_noerr) RETURN
!
!  Write out perceived catch per unit effort (kg)
!
      status=nf_fwrite2d(ng, iNLM, ncFLEETid(ng), fleetVid(idBpcpue,ng),&
     &                   tFLEETindx(ng), r2dvar,                        &
     &                   LBi, UBi, LBj, UBj, 1.0_r8,                    &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   OCEAN(ng) % percpue(:,:))
      IF (status.ne.nf90_noerr) RETURN
!
!  Write out distance to nearest port (km)
!
      status=nf_fwrite2d(ng, iNLM, ncFLEETid(ng), fleetVid(idBpdist,ng),&
     &                   tFLEETindx(ng), r2dvar,                        &
     &                   LBi, UBi, LBj, UBj, 1.0_r8,                    &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   OCEAN(ng) % distport(:,:))
      IF (status.ne.nf90_noerr) RETURN
!
!  Write out nearest port number
!
      status=nf_fwrite2d(ng, iNLM, ncFLEETid(ng), fleetVid(idBpnear,ng),&
     &                   tFLEETindx(ng), r2dvar,                        &
     &                   LBi, UBi, LBj, UBj, 1.0_r8,                    &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   REAL(OCEAN(ng) % nearport(:,:),r8))
      IF (status.ne.nf90_noerr) RETURN
!
!-----------------------------------------------------------------------
!  Synchronize fleet NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, FLEETname(ng), ncFLEETid(ng))
!
#else
      SUBROUTINE wrt_fleet
#endif
      RETURN
      END SUBROUTINE wrt_fleet
