#include "cppdefs.h"

#if defined OPTIC_MANIZZA && defined SOLVE3D

      SUBROUTINE optic_manizza_tile (ng, tile,                          &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Zscale, TZ, ZW, CCHL, SSWDK)
!
!svn $Id$
!================================================== Raphael Dussin   ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
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
!     ZW       depth array 3d                                          !
!     CCHL     chlorophyl concentration                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     SSWDK    shortwave (radiation) fractional decay.                 !
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

      real(r8), intent(in) :: Zscale

      real(r8), intent(in)  :: TZ(IminS:ImaxS,JminS:JmaxS)
      real(r8), intent(in)  :: ZW(IminS:ImaxS,JminS:JmaxS,0:N(ng))
      real(r8), intent(in)  :: CCHL(IminS:ImaxS,JminS:JmaxS,0:N(ng))
      real(r8), intent(out) :: SSWDK(IminS:ImaxS,JminS:JmaxS)
!
!  Local variable declarations.
!
      integer :: Jindex, i, j, k

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

      real(r8) :: swdk_blue(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: swdk_red( IminS:ImaxS,JminS:JmaxS)

      real(r8) :: swdk_ir( IminS:ImaxS,JminS:JmaxS)
      real(r8) :: zref_swf( IminS:ImaxS,JminS:JmaxS)

# include "set_bounds.h"

!    PRINT *, 'IJmin/max in tile ', tile, ' : (', ImaxS - IminS, ';', JmaxS - JminS, ')'
!    PRINT *, 'IJstr/end in tile ', tile, ' : (', Iend - Istr, ';', Jend - Jstr, ')'

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
           ! Compute the reference depth at which we want the shortwave decay fraction
           ! this is basically getting back the depth from thickness given in input
           ! this will be useful to find out on which layer we stop the recursive
           ! computation of blue and red decays
           zref_swf(i,j) = ZW(i,j,N(ng)) - TZ(i,j) 

           ! define visible and IR fraction according to variable water type
           Jindex     = INT(MIXING(ng)%Jwtype(i,j))
           ir_frac    = lmd_r1(Jindex)
           ! define InfraRed decay scale according to variable water type
           k_ir = 1.0 / lmd_mu1(Jindex)
           visib_frac = 1.0 - ir_frac
           ! compute the values of red and blue bands according to the Jerlov
           ! water type
           swdk_blue(i,j) = visib_frac * blue_frac
           swdk_red(i,j)  = visib_frac * red_frac

           ! 1. this part is about the InfraRed part that still use the Jerlov
           ! variable water type scheme
           ! define visible and IR fraction according to variable water type
           ! compute the IR decay that will be added at the end
           swdk_ir(i,j) = ir_frac * exp( Zscale * k_ir * TZ(i,j) )

        END DO
      END DO

      DO k=N(ng)-1,0,-1 
       DO j=Jstr,Jend
        DO i=Istr,Iend
            ! 2. Compute recursively the blue+red fractional decays to apply at
            ! depth = zref_swf. When we reach the desired value, loop continues
            ! but multiplies by one the obtained results.
            ! We loop from subsurface to bottom

            ! define visible decay scale at level k using chlorophyll at W-point
            k_red  = ksw_red  + xsi_red  * ( CCHL(i,j,k) ** eps_red  )
            k_blue = ksw_blue + xsi_blue * ( CCHL(i,j,k) ** eps_blue )

            ! dz will be the thickness to use in the exponential. If desired
            ! depth is below level k then dz = thickness of the layer k / k+1.
            ! If desired depth is within the layer then dz = fraction of the
            ! layer k / k+1. If it is above, then dz = 0 and we multiply the
            ! fractional decay by one down to the bottom.

            ! dz1 = layer thickness (negative)
            dz1 = ZW(i,j,k) - ZW(i,j,k+1)
            ! dz2 = distance between zref_swf and upper W-point (negative or 0 if
            ! zref is above)
            dz2 = MIN(zref_swf(i,j) - ZW(i,j,k+1), 0.)
            ! we take the smallest value (MAX because values are negative or zero)
            ! if dz = 0, we just multiply by one the final array
            dz0 = MAX(dz1,dz2)

            swdk_red(i,j)  = swdk_red(i,j)  * EXP( k_red  * dz0 )
            swdk_blue(i,j) = swdk_blue(i,j) * EXP( k_blue * dz0 )

        END DO
       END DO
      END DO

      ! compute the sum of visible blue/red and InfraRed fractions.
      DO j=Jstr,Jend
        DO i=Istr,Iend
          SSWDK(i,j) = swdk_red(i,j) + swdk_blue(i,j) + swdk_ir(i,j)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE optic_manizza_tile

      SUBROUTINE optic_manizza3_tile (ng, tile,                          &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            Zscale, ZW, CCHL, SSWDK)
!
!svn $Id$
!================================================== Raphael Dussin   ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
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
!     ZW       depth array 3d                                          !
!     CCHL     chlorophyl concentration                                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     SSWDK    shortwave (radiation) fractional decay.                 !
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
      integer, parameter  :: nbands=4 ! RD
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj 
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS

      real(r8), intent(in) :: Zscale

      real(r8), intent(in)  :: ZW(IminS:ImaxS,JminS:JmaxS,0:N(ng))
      real(r8), intent(in)  :: CCHL(IminS:ImaxS,JminS:JmaxS,0:N(ng))
      real(r8), intent(out) :: SSWDK(IminS:ImaxS,JminS:JmaxS,0:N(ng),nbands)
!
!  Local variable declarations.
!
      integer :: Jindex, i, j, k

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
!      real(r8) :: zref_swf( IminS:ImaxS,JminS:JmaxS)

# include "set_bounds.h"

    ! According to Manizza paper the fraction R is defined as :
    ! I_{IR} = R x I_{0} = 0.58 x I_{0} and I_{vis} = (1-R) x I_{0} = 0.42 x I_{0}
    ! using R=0.58 which correspond to type I Jerlov Water type.
    ! here R will be set by the water type at each point (variable water type)

    ! The paper also states an equi-repartition of the radiation between red and blue
    ! bands

    !PRINT *, 'in tile ', tile, ' : (', IminS, '-', ImaxS, ';', JminS, '-', JmaxS, ')'
    !PRINT *, 'IJmin/max in tile ', tile, ' : (', ImaxS - IminS, ';', JmaxS - JminS, ')'
    !PRINT *, 'IJstr/end in tile ', tile, ' : (', Iend - Istr, ';', Jend - Jstr, ')'

      red_frac=0.5
      blue_frac=0.5

    ! These are the values for different coeficients :

      ksw_red  = 0.225  ! m-1
      ksw_blue = 0.0232 ! m-1

      xsi_red  = 0.037  ! m-2 mgChl m-3
      xsi_blue = 0.074  ! m-2 mgChl m-3

      eps_red  = 0.629  ! no units
      eps_blue = 0.674  ! no units

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
           ! Compute the reference depth at which we want the shortwave decay fraction
           ! this is basically getting back the depth from thickness given in input
           ! this will be useful to find out on which layer we stop the recursive
           ! computation of blue and red decays
           !zref_swf(i,j) = ZW(i,j,N(ng)) - TZ(i,j) 

           ! define visible and IR fraction according to variable water type
           Jindex     = INT(MIXING(ng)%Jwtype(i,j))
           ir_frac    = lmd_r1(Jindex)
           ! define InfraRed decay scale according to variable water type
           k_ir = 1.0 / lmd_mu1(Jindex)
           visib_frac = 1.0 - ir_frac
           ! compute the values of red and blue bands according to the Jerlov
           ! water type
           !swdk_blue(i,j,0) = visib_frac * blue_frac
           !swdk_red(i,j,0)  = visib_frac * red_frac
           ! RD was bugged ?!
           swdk_blue(i,j,N(ng)) = visib_frac * blue_frac
           swdk_red(i,j,N(ng))  = visib_frac * red_frac

           ! 1. this part is about the InfraRed part that still use the Jerlov
           ! variable water type scheme
           ! define visible and IR fraction according to variable water type
           ! compute the IR decay that will be added at the end
           !swdk_ir(i,j) = ir_frac * exp( Zscale * k_ir * TZ(i,j) )

           SSWDK(i,j,N(ng),1) = 1.
           SSWDK(i,j,N(ng),2) = ir_frac
           SSWDK(i,j,N(ng),3) = visib_frac * blue_frac
           SSWDK(i,j,N(ng),4) = visib_frac * red_frac
        END DO
      END DO

      ! From subsurface to bottom
      DO k=N(ng)-1,0,-1

       DO j=Jstr,Jend
        DO i=Istr,Iend
           ! define visible and IR fraction according to variable water type
           Jindex     = INT(MIXING(ng)%Jwtype(i,j))
           ir_frac    = lmd_r1(Jindex)
!           visib_frac = 1.0 - ir_frac ! has to be used somewhere !!!!!!! OK !
           ! define InfraRed decay scale according to variable water type
           k_ir = 1.0 / lmd_mu1(Jindex)
           ! define visible decay scale using chlorophyll
           ! replace chl_test by real array !!!
           k_red  = ksw_red  + xsi_red  * ( max( CCHL(i,j,k), 0. ) ** eps_red  )
           k_blue = ksw_blue + xsi_blue * ( max( CCHL(i,j,k), 0. ) ** eps_blue )

           ! compute the decay in the red and blue bands
           !dz0 = ZW(i,j,k+1) - ZW(i,j,k)
           dz0 = ZW(i,j,k) - ZW(i,j,k+1)

           swdk_red(i,j,k)  = swdk_red(i,j,k+1)  * exp( - k_red  * dz0 )
           swdk_blue(i,j,k) = swdk_blue(i,j,k+1) * exp( - k_blue * dz0 )

           ! Total
           SSWDK(i,j,k,1) = swdk_red(i,j,k) + swdk_blue(i,j,k) + &
       &                 ir_frac * exp( - k_ir * ZW(i,j,k) )
           ! IR
           SSWDK(i,j,k,2) = ir_frac * exp( - k_ir * ZW(i,j,k) )
           ! Blue band
           SSWDK(i,j,k,3) = swdk_blue(i,j,k) 
           ! Red band
           SSWDK(i,j,k,4) = swdk_red(i,j,k) 

        END DO
       END DO
      END DO

!      PRINT *, 'RD in optic3 total', MINVAL(SSWDK(:,:,:,1)), MAXVAL(SSWDK(:,:,:,1))
!      PRINT *, 'RD in optic3 IR ', MINVAL(SSWDK(:,:,:,2)), MAXVAL(SSWDK(:,:,:,2))
!      PRINT *, 'RD in optic3 blue ', MINVAL(SSWDK(:,:,:,3)), MAXVAL(SSWDK(:,:,:,3))
!      PRINT *, 'RD in optic3 red ', MINVAL(SSWDK(:,:,:,4)), MAXVAL(SSWDK(:,:,:,4))


      RETURN
      END SUBROUTINE optic_manizza3_tile
#else
      SUBROUTINE optic_manizza
      RETURN
      END SUBROUTINE optic_manizza

      SUBROUTINE optic_manizza3
      RETURN
      END SUBROUTINE optic_manizza3
#endif
