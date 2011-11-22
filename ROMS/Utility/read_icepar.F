#include "cppdefs.h"
#ifdef ICE_MODEL
      SUBROUTINE read_IcePar (inp, out, Lwrite)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads and reports ice model input parameters.          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: inp, out
!
!  Local variable declarations.
!
      integer :: Lstr, Lval, Npts, Nval, i, ng, itrc, status

      integer :: decode_line, lenstr, load_i, load_l, load_r

      real(r8), dimension(100) :: Rval

      character (len=40) :: KeyWord
      character (len=80) :: line
      character (len=80), dimension(100) :: Cval

!
!-----------------------------------------------------------------------
!  Read in ice model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.true.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lice') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lice)
          ELSE IF (TRIM(KeyWord).eq.'DTICE') THEN
            Npts=load_r(Nval, Rval, Ngrids, dtice)
          ELSE IF (TRIM(KeyWord).eq.'DTICE_EQ') THEN
            Npts=load_r(Nval, Rval, Ngrids, dtice_eq)
          ELSE IF (TRIM(KeyWord).eq.'nstrs') THEN
            Npts=load_i(Nval, Rval, Ngrids, nstrs)
          ELSE IF (TRIM(KeyWord).eq.'nevp') THEN
            Npts=load_i(Nval, Rval, Ngrids, nevp)
          ELSE IF (TRIM(KeyWord).eq.'rhoice') THEN
            Npts=load_r(Nval, Rval, Ngrids, rhoice)
          ELSE IF (TRIM(KeyWord).eq.'cdiw') THEN
            Npts=load_r(Nval, Rval, Ngrids, cdiw)
          ELSE IF (TRIM(KeyWord).eq.'cdai') THEN
            Npts=load_r(Nval, Rval, Ngrids, cdai)
          ELSE IF (TRIM(KeyWord).eq.'rho_air') THEN
            Npts=load_r(Nval, Rval, Ngrids, rho_air)
          ELSE IF (TRIM(KeyWord).eq.'rhosnow_dry') THEN
            Npts=load_r(Nval, Rval, Ngrids, rhosnow_dry)
          ELSE IF (TRIM(KeyWord).eq.'rhosnow_wet') THEN
            Npts=load_r(Nval, Rval, Ngrids, rhosnow_wet)
          ELSE IF (TRIM(KeyWord).eq.'pstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, pstar)
          ELSE IF (TRIM(KeyWord).eq.'astren') THEN
            Npts=load_r(Nval, Rval, Ngrids, astren)
          ELSE IF (TRIM(KeyWord).eq.'zetamax') THEN
            Npts=load_r(Nval, Rval, Ngrids, zetamax)
          ELSE IF (TRIM(KeyWord).eq.'zetamin') THEN
            Npts=load_r(Nval, Rval, Ngrids, zetamin)
          ELSE IF (TRIM(KeyWord).eq.'ellip_sq') THEN
            Npts=load_r(Nval, Rval, Ngrids, ellip_sq)
          ELSE IF (TRIM(KeyWord).eq.'alphai') THEN
            Npts=load_r(Nval, Rval, Ngrids, alphai)
            do ng=1,Ngrids
               alphai(ng) = alphai(ng)*deg2rad
            enddo
          ELSE IF (TRIM(KeyWord).eq.'tol') THEN
            Npts=load_r(Nval, Rval, 1, tol)
          ELSE IF (TRIM(KeyWord).eq.'min_h') THEN
            Npts=load_r(Nval, Rval, Ngrids, min_h)
          ELSE IF (TRIM(KeyWord).eq.'min_a') THEN
            Npts=load_r(Nval, Rval, Ngrids, min_a)
          ELSE IF (TRIM(KeyWord).eq.'max_a') THEN
            Npts=load_r(Nval, Rval, Ngrids, max_a)
          ELSE IF (TRIM(KeyWord).eq.'sfwat_max') THEN
            Npts=load_r(Nval, Rval, Ngrids, sfwat_max)
          ELSE IF (TRIM(KeyWord).eq.'stressang') THEN
            Npts=load_r(Nval, Rval, Ngrids, stressang)
            do ng=1,Ngrids
               stressang(ng) = stressang(ng)*deg2rad
            enddo
          ELSE IF (TRIM(KeyWord).eq.'ice_emiss') THEN
            Npts=load_r(Nval, Rval, 1, ice_emiss)
          ELSE IF (TRIM(KeyWord).eq.'spec_heat_air') THEN
            Npts=load_r(Nval, Rval, 1, spec_heat_air)
          ELSE IF (TRIM(KeyWord).eq.'trans_coeff') THEN
            Npts=load_r(Nval, Rval, 1, trans_coeff)
          ELSE IF (TRIM(KeyWord).eq.'sublim_latent_heat') THEN
            Npts=load_r(Nval, Rval, 1, sublim_latent_heat)
          ELSE IF (TRIM(KeyWord).eq.'t0deg') THEN
            Npts=load_r(Nval, Rval, 1, t0deg)
          END IF
        END IF
      END DO
  10  IF (Master) WRITE (out,30) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)

! Set ice time step to ocean time step

      DO ng = 1,Ngrids
        dtice(ng) = dt(ng)
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lice(ng)) THEN
            WRITE (out,40) ng
            WRITE(out,*) 'Ice time step = ocean time step'
            WRITE (out,100) dtice(ng), 'DTICE',                         &
     &            'Ice model time step (s).'
            WRITE (out,100) dtice_eq(ng), 'DTICE_EQ',                   &
     &            'Ice drift update (equilibrium) time step (s).'
            WRITE (out,50) nstrs(ng), 'nstrs',                          &
     &            'Number of iterations for nonlinear ice dynamics.'
            WRITE (out,50) nevp(ng), 'nevp',                            &
     &            'Number of elastic steps per plastic step in EVP.'
            WRITE (out,100) rhoice(ng), 'rhoice',                       &
     &            'Density of sea ice (kg/m3).'
            WRITE (out,100) cdiw(ng), 'cdiw',                           &
     &            'Ice-water drag coefficient (nondimensional).'
            WRITE (out,100) cdai(ng), 'cdai',                           &
     &            'Air-ice drag coefficient (nondimensional).'
            WRITE (out,100) rho_air(ng), 'rho_air',                     &
     &            'Air density (kg/m3).'
            WRITE (out,100) rhosnow_dry(ng), 'rhosnow_dry',             &
     &            'Dry snow density (kg/m3).'
            WRITE (out,100) rhosnow_wet(ng), 'rhosnow_wet',             &
     &            'Wet snow density (kg/m3).'
            WRITE (out,100) alphai(ng)*rad2deg, 'alphai',               &
     &            'Mohr-Coulomb stress angle (degrees).'
            WRITE (out,100) min_h(ng), 'min_h',                         &
     &            'Minimum average ice thickness (m).'
            WRITE (out,100) min_a(ng), 'min_a',                         &
     &            'Minimum ice concentration (nondimensional).'
            WRITE (out,100) max_a(ng), 'max_a',                         &
     &            'Maximum ice concentration (nondimensional).'
            WRITE (out,100) sfwat_max(ng)*rad2deg, 'sfwat_max',         &
     &            'Maximum surface fresh water (nondimensional).'
            WRITE (out,100) stressang(ng)*rad2deg, 'stressang',         &
     &            'Turning angle for ice-water drag (degrees).'
            WRITE (out,100) tol, 'tol',                                 &
     &            'Numerical tolerance in rheology calculations .'
            WRITE (out,100) ice_emiss, 'ice_emiss',                     &
     &            'Ice emissivity.'
            WRITE (out,100) spec_heat_air, 'spec_heat_air',             &
     &            'Specific heat of air.'
            WRITE (out,100) trans_coeff, 'trans_coeff',                 &
     &            'Transfer coefficient.'
            WRITE (out,100) sublim_latent_heat, 'sublim_latent_heat',   &
     &            'Latent_heat of sublimation.'
            WRITE (out,100) t0deg, 't0deg',                             &
     &            'Zero degrees Celsius in degrees Kelvin.'
          END IF
        END DO
      END IF

  30  FORMAT (/,' READ_IcePar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' Ice Parameters, Grid: ',i2.2,                   &
     &        /,' ============================',/)
  50  FORMAT (1x,i10,2x,a,t28,a)
  60  FORMAT (10x,l1,2x,a,t28,a,i2.2,':',1x,a)
  70  FORMAT (f11.3,2x,a,t28,a)
  80  FORMAT (f11.3,2x,a,t28,a,/,t30,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t28,a,/,t30,a,i2.2,':',1x,a)
 100  FORMAT (1p,e11.4,2x,a,t28,a)
 110  FORMAT (1p,e11.4,2x,a,t28,a,/,t30,a)
 120  FORMAT (10x,l1,2x,a,t28,a)
 280  FORMAT (/,' READ_IcePar - variable info not yet loaded, ', a)
      RETURN
      END SUBROUTINE read_IcePar
#endif
