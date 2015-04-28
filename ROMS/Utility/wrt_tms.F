#include "cppdefs.h"
#ifdef TIMESERIES
      SUBROUTINE wrt_tms (ng)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes out data into stations NetCDF file.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: Fcount, i, ifield
      real(r8) :: scale

      SourceFile='wrt_tms.F'
!
!-----------------------------------------------------------------------
!  Write out timeserie data 
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Set time record index.
!
      TMS(ng)%Rindex=TMS(ng)%Rindex+1
      Fcount=TMS(ng)%Fcount
      TMS(ng)%Nrec(Fcount)=TMS(ng)%Nrec(Fcount)+1
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, TMS(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/TMS(ng)%Rindex/), (/1/),                  &
     &                      ncid = TMS(ng)%ncid,                        &
     &                      varid = TMS(ng)%Vid(idtime))
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out tracer type variables.
!
      DO i=1,NT(ng)
        ifield=idTvar(i)
          scale=1.0_r8
          CALL netcdf_put_fvar (ng, iNLM, TMS(ng)%name,                 &
     &                          TRIM(Vname(1,idTvar(i))),               &
     &                          OCEAN(ng)%tms(1,1,i) ,                  &
     &                          (/TMS(ng)%Rindex/),                     &
     &                          (/1/),                                  &
     &                          ncid = TMS(ng)%ncid,                    &
     &                          varid = TMS(ng)%Tid(i))
          IF (exit_flag.ne.NoError) RETURN
      END DO
!
!
!-----------------------------------------------------------------------
!  Synchronize stations NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, TMS(ng)%name, TMS(ng)%ncid)

#else
      SUBROUTINE wrt_tms
#endif
      RETURN
      END SUBROUTINE wrt_tms
