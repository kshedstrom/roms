#include "cppdefs.h"

MODULE optic_manizza_mod

#if defined OPTIC_MANIZZA && defined SOLVE3D
!
!svn $Id$
!================================================== Raphael Dussin   ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  From GFDL? See below for more comments
!
      implicit none

      PRIVATE
      PUBLIC  :: optic_manizza

      CONTAINS

      SUBROUTINE optic_manizza (ng, tile)

      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping, only: nstp
#ifdef BIO_COBALT
     USE mod_biology, only: iochl
#elif defined BIO_UMAINE
     USE mod_biology, only: iS1CH, iS2CH, iS3CH
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.

# include "tile.h"

      CALL optic_manizza_tile (ng, tile,                                &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            GRID(ng) % z_w,                       &
#ifdef MASKING
     &                            GRID(ng) % rmask,                     &
#endif
#ifdef BIO_COBALT
     &                            OCEAN(ng) % obgc(:,:,:,nstp,iochl),   &
#else
     &                            FORCES(ng) % chl,                     &
#endif
     &                            OCEAN(ng) % decayW)

      RETURN

      END SUBROUTINE optic_manizza

      SUBROUTINE optic_manizza_tile (ng, tile,                          &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            z_w,                                  &
#ifdef MASKING
     &                            rmask,                                &
#endif
     &                            chl,                                  &
     &                            decayW)
!
!=======================================================================
!                                                                      !
!  This subroutine computes the  fraction  of  solar shortwave flux    !
!  penetrating to specified depth (times Zscale) due to exponential    !
!  decay using Manizza scheme (using chlorophyll data or model)        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Zscale   scale factor to apply to depth array.                   !
!     TZ       vertical depth (meters, positive) for                   !
!              desired solar short-wave fraction.                      !
!     z_w       depth array 3d                                         !
!     chl     chlorophyl concentration                                 !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     decayW    shortwave (radiation) fractional decay.                !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Paulson, C.A., and J.J. Simpson, 1977: Irradiance meassurements     !
!     in the upper ocean, J. Phys. Oceanogr., 7, 952-956.              !
!                                                                      !
!  Manizza, M et al, 2005 : Bio-optical feedbacks among phytoplankton, !
!  upper ocean physics and sea-ice in a global model, Geophys. Rev.    !
!  Letters                                                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_mixing
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj 
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

      real(r8), intent(in)  :: z_w(LBi:UBi,LBj:UBj,0:N(ng))
#ifdef MASKING
      real(r8), intent(in)  :: rmask(LBi:UBi,LBj:UBj)
#endif
      real(r8), intent(in)  :: chl(LBi:UBi,LBj:UBj,N(ng)) 
      real(r8), intent(out) :: decayW(LBi:UBi,LBj:UBj,0:N(ng),4)
!
!  Local variable declarations.
!
      integer :: Jindex, i, j, k, nb

      real(r8) :: visib_frac        ! fraction of incoming shortwave going to visible
      real(r8) :: ir_frac           ! fraction of incoming shortwave going to near infra-red (not thermal)
      real(r8) :: red_frac          ! fraction of incoming visible going to red band
      real(r8) :: blue_frac         ! fraction of incoming visible going to blue band

      real(r8) :: k_ir              ! total inverse e-folding scale for infrared
      real(r8) :: k_red, k_blue     ! total inverse e-folding scale (red and blue)
      real(r8) :: ksw_red, ksw_blue ! inverse e-folding scale for pure seawter (red and blue)
      real(r8) :: xsi_red, xsi_blue ! coeficient for chlorophyll contribution (red and blue)
      real(r8) :: eps_red, eps_blue ! exponent for chlorophyll contribution (red and blue)

      real(r8) :: dz0, dz1, dz2
      real(r8) :: expo1, expo2

      real(r8) :: swdk_blue(IminS:ImaxS,JminS:JmaxS,0:N(ng))
      real(r8) :: swdk_red( IminS:ImaxS,JminS:JmaxS,0:N(ng))

      real(r8) :: swdk_ir( IminS:ImaxS,JminS:JmaxS,0:N(ng))
      real(r8) :: chl_conc( IminS:ImaxS,JminS:JmaxS,0:N(ng))
      real(r8) :: depth_w( IminS:ImaxS,JminS:JmaxS,0:N(ng))

# include "set_bounds.h"

    ! According to Manizza paper the fraction R is defined as :
    ! I_{IR} = R x I_{0} = 0.58 x I_{0} and I_{vis} = (1-R) x I_{0} = 0.42 x I_{0}
    ! using R=0.58 which correspond to type I Jerlov Water type.
    ! here R will be set by the water type at each point (variable water type)

    ! The paper also states an equi-repartition of the radiation between red and blue
    ! bands

      red_frac=0.5
      blue_frac=0.5

    ! These are the values for different coeficients :

      ksw_red  = 0.225  ! m-1
      ksw_blue = 0.0232 ! m-1

      xsi_red  = 0.037  ! m-2 mgChl m-3
      xsi_blue = 0.074  ! m-2 mgChl m-3

      eps_red  = 0.629  ! no units
      eps_blue = 0.674  ! no units


#ifdef BIO_COBALT
   ! the rho0 / 1000 is a unit conversion from ug Chl.kg-1 of Cobalt to
   ! mg Chl.m-3 used in the light penetration ( manizza ) routine

   DO k=1,N(ng)
     DO j=Jstr,Jend
       DO i=Istr,Iend
          chl_conc(i,j,k) = chl(i,j,k) * rho0 / 1000.0d0
       END DO
     END DO
   ENDDO
#else
   DO k=1,N(ng)
     DO j=Jstr,Jend
       DO i=Istr,Iend
          chl_conc(i,j,k) = chl(i,j,k)
       END DO
     END DO
   ENDDO
#endif

   ! 4. Compute depth at w-point
   DO k=0,N(ng)
     DO j=Jstr,Jend
       DO i=Istr,Iend
         depth_w(i,j,k) = z_w(i,j,N(ng)) - z_w(i,j,k)
       ENDDO
     ENDDO
   ENDDO
!
!-----------------------------------------------------------------------
!  Use Paulson and Simpson (1977) two wavelength bands solar
!  absorption model.
!-----------------------------------------------------------------------
!

      ! Compute fraction at the surface k = N(ng)
      ! total fraction will be = 1.
      ! this is where we set the values for blue and red fractions at the
      ! surface that depends on the Jerlov water types (at least for the
      ! fraction that goes into visible)
      DO j=Jstr,Jend
        DO i=Istr,Iend
           ! define visible and IR fraction according to variable water type
           Jindex     = INT(MIXING(ng)%Jwtype(i,j))
           ir_frac    = lmd_r1(Jindex)
           ! define InfraRed decay scale according to variable water type
           k_ir = 1.0 / lmd_mu1(Jindex)
           visib_frac = 1.0 - ir_frac
           ! compute the values of red and blue bands according to the Jerlov
           ! water type
           swdk_blue(i,j,N(ng)) = visib_frac * blue_frac
           swdk_red(i,j,N(ng))  = visib_frac * red_frac

           ! 1. this part is about the InfraRed part that still use the Jerlov
           ! variable water type scheme
           ! define visible and IR fraction according to variable water type
           ! compute the IR decay that will be added at the end

           decayW(i,j,N(ng),1) = 1.
           decayW(i,j,N(ng),2) = ir_frac
           decayW(i,j,N(ng),3) = visib_frac * blue_frac
           decayW(i,j,N(ng),4) = visib_frac * red_frac
        END DO
      END DO

      ! From subsurface to bottom
      DO k=N(ng)-1,0,-1

       DO j=Jstr,Jend
        DO i=Istr,Iend
           ! define visible and IR fraction according to variable water type
           Jindex     = INT(MIXING(ng)%Jwtype(i,j))
           ir_frac    = lmd_r1(Jindex)
           ! define InfraRed decay scale according to variable water type
           k_ir = 1.0 / lmd_mu1(Jindex)
           ! define visible decay scale using chlorophyll
           ! we use the chlorophyll concentration between two interfaces
           k_red  = ksw_red  + xsi_red  * ( max( chl_conc(i,j,k+1), 0. ) ** eps_red  )
           k_blue = ksw_blue + xsi_blue * ( max( chl_conc(i,j,k+1), 0. ) ** eps_blue )

           ! compute the decay in the red and blue bands
           dz0 = z_w(i,j,k+1) - z_w(i,j,k)

           swdk_red(i,j,k)  = swdk_red(i,j,k+1)  * exp( - k_red  * dz0 )
           swdk_blue(i,j,k) = swdk_blue(i,j,k+1) * exp( - k_blue * dz0 )

           ! IR
           decayW(i,j,k,2) = ir_frac * exp( - k_ir * depth_w(i,j,k) )
           ! Blue band
           decayW(i,j,k,3) = swdk_blue(i,j,k) 
           ! Red band
           decayW(i,j,k,4) = swdk_red(i,j,k) 

           ! Total
           decayW(i,j,k,1) = decayW(i,j,k,2) + decayW(i,j,k,3) + decayW(i,j,k,4)
        END DO
       END DO
      END DO

#ifdef MASKING
      DO nb=1,4
       DO k=0,N(ng)
        DO j=Jstr,Jend
         DO i=Istr,Iend
           decayW(i,j,k,nb) = decayW(i,j,k,nb) * rmask(i,j)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif

      RETURN
      END SUBROUTINE optic_manizza_tile
#endif

END MODULE optic_manizza_mod
