#include "cppdefs.h"
#if defined NEMURO_SAN && defined FISHING_FLEET
      SUBROUTINE wrt_fleet (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
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
      FLEET(ng)%Rindex=FLEET(ng)%Rindex+1
      FLEET(ng)%Nrec=FLEET(ng)%Nrec+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, FLEET(ng)%name,                   &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/FLEET(ng)%Rindex/), (/1/),                &
     &                      ncid = FLEET(ng)%ncid,                      &
     &                      varid = FLEET(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet X-grid locations.
!
      DO l=1,Nboats(ng)
        Tout(l)=REAL(BOATS(ng)%boat(ibiloc,l),r8)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEET(ng)%name,                   &
     &                      TRIM(Vname(1,idBxprt)), Tout,               &
     &                      (/1,FLEET(ng)%Rindex/), (/Nboats(ng),1/),   &
     &                      ncid = FLEET(ng)%ncid,                      &
     &                      varid = FLEET(ng)%Vid(idBxprt))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet Y-grid locations.
!
      DO l=1,Nboats(ng)
        Tout(l)=REAL(BOATS(ng)%boat(ibjloc,l),r8)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEET(ng)%name,                   &
     &                      TRIM(Vname(1,idByprt)), Tout,               &
     &                      (/1,FLEET(ng)%Rindex/), (/Nboats(ng),1/),   &
     &                      ncid = FLEET(ng)%ncid,                      &
     &                      varid = FLEET(ng)%Vid(idByprt))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet port locations.
!
      DO l=1,Nboats(ng)
        Tout(l)=REAL(BOATS(ng)%boat(ibport,l),r8)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEET(ng)%name,                   &
     &                      TRIM(Vname(1,idBnprt)), Tout,               &
     &                      (/1,FLEET(ng)%Rindex/), (/Nboats(ng),1/),   &
     &                      ncid = FLEET(ng)%ncid,                      &
     &                      varid = FLEET(ng)%Vid(idBnprt))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet number of fished locations.
!
      DO l=1,Nboats(ng)
        Tout(l)=REAL(BOATS(ng)%boat(ibnloc,l),r8)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEET(ng)%name,                   &
     &                      TRIM(Vname(1,idBnfloc)), Tout,              &
     &                      (/1,FLEET(ng)%Rindex/), (/Nboats(ng),1/),   &
     &                      ncid = FLEET(ng)%ncid,                      &
     &                      varid = FLEET(ng)%Vid(idBnfloc))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet X-grid 1st fishing locations.
!
      DO l=1,Nboats(ng)
        Tout(l)=REAL(BOATS(ng)%boat(ibifsh1,l),r8)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEET(ng)%name,                   &
     &                      TRIM(Vname(1,idBxfsh1)), Tout,              &
     &                      (/1,FLEET(ng)%Rindex/), (/Nboats(ng),1/),   &
     &                      ncid = FLEET(ng)%ncid,                      &
     &                      varid = FLEET(ng)%Vid(idBxfsh1))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet Y-grid 1st fishing locations.
!
      DO l=1,Nboats(ng)
        Tout(l)=REAL(BOATS(ng)%boat(ibjfsh1,l),r8)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEET(ng)%name,                   &
     &                      TRIM(Vname(1,idByfsh1)), Tout,              &
     &                      (/1,FLEET(ng)%Rindex/), (/Nboats(ng),1/),   &
     &                      ncid = FLEET(ng)%ncid,                      &
     &                      varid = FLEET(ng)%Vid(idByfsh1))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet X-grid 2nd fishing locations.
!
      DO l=1,Nboats(ng)
        Tout(l)=REAL(BOATS(ng)%boat(ibifsh2,l),r8)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEET(ng)%name,                   &
     &                      TRIM(Vname(1,idBxfsh2)), Tout,              &
     &                      (/1,FLEET(ng)%Rindex/), (/Nboats(ng),1/),   &
     &                      ncid = FLEET(ng)%ncid,                      &
     &                      varid = FLEET(ng)%Vid(idBxfsh2))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet Y-grid 2nd fishing locations.
!
      DO l=1,Nboats(ng)
        Tout(l)=REAL(BOATS(ng)%boat(ibjfsh2,l),r8)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEET(ng)%name,                   &
     &                      TRIM(Vname(1,idByfsh2)), Tout,              &
     &                      (/1,FLEET(ng)%Rindex/), (/Nboats(ng),1/),   &
     &                      ncid = FLEET(ng)%ncid,                      &
     &                      varid = FLEET(ng)%Vid(idByfsh2))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet daily catch.
!
      DO l=1,Nboats(ng)
        Tout(l)=BOATS(ng)%catch(l)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEET(ng)%name,                   &
     &                      TRIM(Vname(1,idBcatch)), Tout,              &
     &                      (/1,FLEET(ng)%Rindex/), (/Nboats(ng),1/),   &
     &                      ncid = FLEET(ng)%ncid,                      &
     &                      varid = FLEET(ng)%Vid(idBcatch))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out fleet cummulative annual catch.
!
      DO l=1,Nboats(ng)
        Tout(l)=BOATS(ng)%cumcatch(l)
      END DO
      CALL netcdf_put_fvar (ng, iNLM, FLEET(ng)%name,                   &
     &                      TRIM(Vname(1,idBcmctch)), Tout,             &
     &                      (/1,FLEET(ng)%Rindex/), (/Nboats(ng),1/),   &
     &                      ncid = FLEET(ng)%ncid,                      &
     &                      varid = FLEET(ng)%Vid(idBcmctch))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out average catch per unit effort (kg)
!
      status=nf_fwrite2d(ng, iNLM, FLEET(ng)%ncid,                      &
     &                   FLEET(ng)%Vid(idBacpue),                       &
     &                   FLEET(ng)%Rindex, r2dvar,                      &
     &                   LBi, UBi, LBj, UBj, 1.0_r8,                    &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   OCEAN(ng) % avgcpue(:,:))
      IF (status.ne.nf90_noerr) RETURN
!
!  Write out perceived catch per unit effort (kg)
!
      status=nf_fwrite2d(ng, iNLM, FLEET(ng)%ncid,                      &
     &                   FLEET(ng)%Vid(idBpcpue),                       &
     &                   FLEET(ng)%Rindex, r2dvar,                      &
     &                   LBi, UBi, LBj, UBj, 1.0_r8,                    &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   OCEAN(ng) % percpue(:,:))
      IF (status.ne.nf90_noerr) RETURN
!
!  Write out distance to nearest port (km)
!
      status=nf_fwrite2d(ng, iNLM, FLEET(ng)%ncid,                      &
     &                   FLEET(ng)%Vid(idBpdist),                       &
     &                   FLEET(ng)%Rindex, r2dvar,                      &
     &                   LBi, UBi, LBj, UBj, 1.0_r8,                    &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   OCEAN(ng) % distport(:,:))
      IF (status.ne.nf90_noerr) RETURN
!
!  Write out nearest port number
!
      status=nf_fwrite2d(ng, iNLM, FLEET(ng)%ncid,                      &
     &                   FLEET(ng)%Vid(idBpnear),                       &
     &                   FLEET(ng)%Rindex, r2dvar,                      &
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
      CALL netcdf_sync (ng, iNLM, FLEET(ng)%name, FLEET(ng)%ncid)
!
#else
      SUBROUTINE wrt_fleet
#endif
      RETURN
      END SUBROUTINE wrt_fleet
