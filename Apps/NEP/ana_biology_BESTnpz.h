!
! svn $Id$
!!================================================== Georgina Gibson ===
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!                                                                      !
!  This routine sets initial conditions for biological tracer fields   !
!  using analytical expressions for the BEST NPZ model .               !
!                                                                      !
!=======================================================================
!

# undef DEPAVG      /* TEST case: set all equal to depth averaged values */

# ifdef DEPAVG
!
! Set all points equal to the same number
!
      DO i=IstrT,IendT
        DO j=JstrT,JendT
          DO k=1,N(ng)
            t(i,j,k,1,iNO3) = 18.0_r8
            t(i,j,k,1,iNH4) = 0.20_r8
            t(i,j,k,1,iPhS) = 0.10_r8
            t(i,j,k,1,iPhL) = 0.10_r8
            t(i,j,k,1,iMZS) = 0.00_r8
            t(i,j,k,1,iMZL) = 0.10_r8
            t(i,j,k,1,iCop) = 2.041_r8
            t(i,j,k,1,iNCa) = 2.417_r8
            t(i,j,k,1,iEup) = 1.337_r8
            t(i,j,k,1,iDet) = eps
            t(i,j,k,1,iDetF) = eps
#  ifdef IRON_LIMIT
            t(i,j,k,1,iFe) = 0.6_r8
#  endif
          END DO
        END DO
      END DO
# else
!
! Make the curves.  First come parameters that are constant through
! entire water column.
!
      deepval(iNO3) = 30.0_r8
      deepval(iNH4) = eps
      deepval(iPhS) = eps
      deepval(iPhL) = eps
      deepval(iMZS) = 0.0_r8
      deepval(iMZL) = eps
      deepval(iCop) = 0.00863_r8
      deepval(iNCaS) = eps
      deepval(iEupS) = eps
      deepval(iNCaO) = 0.0_r8
      deepval(iEupO) = 0.0_r8
      deepval(iDet) = eps
      deepval(iDetF) = eps
      deepval(iJel) = eps
      DO k=1,N(ng)
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            t(i,j,k,1,iNO3) = 14.0_r8 !  18.0_r8
            t(i,j,k,1,iNH4) =  0.0_r8 ! 2.0_r8
            t(i,j,k,1,iDet) =  eps
            t(i,j,k,1,iDetF) =  eps
            biod(i,j,k) = -1 * ( z_r(i,j,k) + 2.5_r8 )
#ifdef JELLY
            t(i,j,k,1,iJel) =  0.1_r8
#endif
          END DO
        END DO
      END DO
!
! PS - a combination of 2 curves: 2nd order polynomial curve to create
! a subsurface maximum at ~20 m and an exponentially decreasing curve to
! assure values aren't negative at depth.
!
      var1 = 25.58_r8
      var2 = -0.250_r8 / (5.0_r8**2)
      var3 =  0.008_r8 / (5.0_r8**3)
      var4 = 3.82_r8 * 5.0_r8
      var5 = 75.0_r8
      var6 = var5 - var4
      var7 = var1 + var2*(var6**2) + var3*(var6**3)
      DO i=IstrT,IendT
        DO j=JstrT,JendT
          DO k=1,N(ng)
            IF ( biod(i,j,k) .le. var5 ) THEN
              var6 = biod(i,j,k) - var4
              t(i,j,k,1,iPhS) =                                         &
     &          var1 + var2*(var6**2) + var3*(var6**3)
            ELSE
              t(i,j,k,1,iPhS) =                                         &
     &          var7 * exp( ( -1.0_r8*biod(i,j,k) + var5 ) / 5.0_r8 )
            END IF
          END DO
        END DO
      END DO
!
! PL - exponentially decreasing with depth.
!
      var1 = 8.25_r8
      var2 = 0.322_r8 / 5.0_r8
      DO i=IstrT,IendT
        DO j=JstrT,JendT
          DO k=1,N(ng)
            t(i,j,k,1,iPhL) = var1 * exp( -1.0_r8 * var2 * biod(i,j,k) )
          END DO
        END DO
      END DO
!
! Microzooplankton, from Howell-Kubler, 1996
! approximated with a straight line and an exponentially decreasing
! curve to assure values aren't negative at depth.  Curves meet ~60m
!
      var1 = 3.1714_r8
      var2 = -0.1865_r8 / 5.0_r8
      var3 = 60.0_r8
      var4 = 0.5_r8
      var5 = var1 + var2 * var3
      DO i=IstrT,IendT
        DO j=JstrT,JendT
          DO k=1,N(ng)
            IF ( biod(i,j,k) .le. var3 ) THEN
              t(i,j,k,1,iMZL) = var1 + var2 * biod(i,j,k)
            ELSE
              t(i,j,k,1,iMZL) = var4 +                                  &
     &          ( var5 - var4 ) * exp( ( var3 - biod(i,j,k) ) / 5.0_r8 )
            END IF
!            t(i,j,k,1,iMZS) = t(i,j,k,1,iMZL)
            t(i,j,k,1,iMZS) = 0.0_r8
          END DO
        END DO
      END DO
!
! Pseudocalanus Copepods - from Shelikof data via Shelikof NPZ
! Step at ~ 50 m
!
      var1 = pi / 2.0_r8
      var2 = 2.876_r8 / pi
      var3 = 5.0_r8 / 5.0_r8
      var4 = 9.151_r8 * 5.0_r8
      DO i=IstrT,IendT
        DO j=JstrT,JendT
          DO k=1,N(ng)
            var5 = var3 * ( biod(i,j,k) - var4 )
            t(i,j,k,1,iCop) = var1 - var2 * atan( var5 )
          END DO
        END DO
      END DO
!
! Neocalanus, from Shelikof NPZ. Step at ~ 30m
!
      var1 = 4.0_r8
      var2 = 1.3_r8
      var3 = 5.0_r8 / 5.0_r8
      var4 = 5.2_r8 * 5.0_r8
      DO i=IstrT,IendT
        DO j=JstrT,JendT
          DO k=1,N(ng)
            var5 = var3 * ( biod(i,j,k) - var4 )
            IF ( -z_r(i,j,1) .gt. 500.0_r8 ) THEN
              t(i,j,k,1,iNCaO) =  0.1_r8
              t(i,j,k,1,iNCaS) = 0._r8
            ELSE
              t(i,j,k,1,iNCaS) =  0.1_r8
              t(i,j,k,1,iNCaO) = 0._r8
            END IF
          END DO
        END DO
      END DO
!
! Euphausiids, Wild guesses. Step at ~ 30m
!
      var1 = 1.78_r8
      var2 = 0.8_r8
      var3 = 5.0_r8 / 5.0_r8
      var4 = 5.2_r8 * 5.0_r8
      DO i=IstrT,IendT
        DO j=JstrT,JendT
          DO k=1,N(ng)
!           var5 = var3 * ( biod(i,j,k) - var4)
            IF (-z_r(i,j,1) .gt. 100.0_r8 ) THEN
              t(i,j,k,1,iEupO) =  0.1_r8
              t(i,j,k,1,iEupS) =  0_r8
            ELSE
              t(i,j,k,1,iEupS) =  0.1_r8
              t(i,j,k,1,iEupO) =  0.0_r8
            END IF
          END DO
        END DO
      END DO
!
#  ifdef IRON_LIMIT
! Iron - linear from surface value to value at 100m and increase onshore
      DO i=IstrT,IendT
        DO j=JstrT,JendT
          var1 = MAX(0._r8,MIN(1._r8,                                   &
     &             (h(i,j)-Feinh)/(Feoffh-Feinh)))
          FeSurf = Feinlo + var1*(Feofflo-Feinlo)
          FeDeep = Feinhi + var1*(Feoffhi-Feinhi)
          var1 = (FeDeep-FeSurf) / 100._r8
          DO k=1,N(ng)
            t(i,j,k,1,iFe) = MIN(FeDeep, FeSurf - z_r(i,j,k)*var1)
          END DO
        END DO
      END DO
#  endif
!
! Concentrations of everything below 100m - i.e. below
! depths where calculations are performed.  Have linear slope
! between values above and below.
! Iron deep values have already been determined.
!
      DO i=IstrT,IendT
        DO j=JstrT,JendT
          DO k=N(ng),1,-1
            IF ( biod(i,j,k) .gt. 300.0_r8 ) THEN       !was 120 for GOA
              t(i,j,k,1,iNO3) = deepval(iNO3)
              t(i,j,k,1,iNH4) = deepval(iNH4)
              t(i,j,k,1,iPhS) = deepval(iPhS)
              t(i,j,k,1,iPhL) = deepval(iPhL)
              t(i,j,k,1,iMZS) = deepval(iMZS)
              t(i,j,k,1,iMZL) = deepval(iMZL)
              t(i,j,k,1,iCop) = deepval(iCop)
              t(i,j,k,1,iNCaS) = deepval(iNCaS)
              t(i,j,k,1,iEupS) = deepval(iEupS)
              t(i,j,k,1,iNCaO) = deepval(iNCaO)
              t(i,j,k,1,iEupO) = deepval(iEupO)
              t(i,j,k,1,iDet) = deepval(iDet)
              t(i,j,k,1,iDetF) = deepval(iDetF)
#ifdef JELLY
              t(i,j,k,1,iJel) = deepval(iJel)
#endif
            ELSE IF ( biod(i,j,k) .gt. 100.0_r8 .and.                   &
     &                biod(i,j,k) .le. 300.0_r8) THEN
              var1 = ( 100.0_r8 - biod(i,j,k) ) / ( 100.0_r8-300.0_r8 )
              t(i,j,k,1,iNO3) = loval(iNO3) +                           &
     &                          ( deepval(iNO3) - loval(iNO3) ) * var1
              t(i,j,k,1,iNH4) = loval(iNH4) +                           &
     &                          ( deepval(iNH4) - loval(iNH4) ) * var1
              t(i,j,k,1,iPhS) = loval(iPhS) +                           &
     &                          ( deepval(iPhS) - loval(iPhS) ) * var1
              t(i,j,k,1,iPhL) = loval(iPhL) +                           &
     &                          ( deepval(iPhL) - loval(iPhL) ) * var1
              t(i,j,k,1,iMZS) = loval(iMZS) +                           &
     &                          ( deepval(iMZS) - loval(iMZS) ) * var1
              t(i,j,k,1,iMZL) = loval(iMZL) +                           &
     &                          ( deepval(iMZL) - loval(iMZL) ) * var1
              t(i,j,k,1,iCop) = loval(iCop) +                           &
     &                          ( deepval(iCop) - loval(iCop) ) * var1
              t(i,j,k,1,iNCaS) = loval(iNCaS) +                         &
     &                          ( deepval(iNCaS) - loval(iNCaS) ) * var1
              t(i,j,k,1,iEupS) = loval(iEupS) +                         &
     &                          ( deepval(iEupS) - loval(iEupS) ) * var1
              t(i,j,k,1,iNCaO) = loval(iNCaO) +                         &
     &                          ( deepval(iNCaO) - loval(iNCaO) ) * var1
              t(i,j,k,1,iEupO) = loval(iEupO) +                         &
     &                          ( deepval(iEupO) - loval(iEupO) ) * var1
              t(i,j,k,1,iDet) = loval(iDet) +                           &
     &                          ( deepval(iDet) - loval(iDet) ) * var1
              t(i,j,k,1,iDetF) = loval(iDetF) +                         &
     &                          ( deepval(iDetF) - loval(iDetF) ) * var1
#ifdef JELLY
              t(i,j,k,1,iJel) = loval(iJel) +                           &
                                ( deepval(iJel) - loval(iJel) ) * var1
#endif
            ELSE
              loval(iNO3) = t(i,j,k,1,iNO3)
              loval(iNH4) = t(i,j,k,1,iNH4)
              loval(iPhS) = t(i,j,k,1,iPhS)
              loval(iPhL) = t(i,j,k,1,iPhL)
              loval(iMZS) = t(i,j,k,1,iMZS)
              loval(iMZL) = t(i,j,k,1,iMZL)
              loval(iCop) = t(i,j,k,1,iCop)
              loval(iNCaS) = t(i,j,k,1,iNCaS)
              loval(iEupS) = t(i,j,k,1,iEupS)
              loval(iNCaO) = t(i,j,k,1,iNCaO)
              loval(iEupO) = t(i,j,k,1,iEupO)
              loval(iDet) = t(i,j,k,1,iDet)
              loval(iDetF) = t(i,j,k,1,iDetF)
              loval(iJel) = t(i,j,k,1,iJel)
            END IF
          END DO
        END DO
      END DO
# endif /* DEPAVG */
#ifdef GAK1D
!
!  This is a hack for sensitivity studies - to test warmer or cooler water
!
!      DO i=IstrT,IendT
!        DO j=JstrT,JendT
!          DO k=1,N(ng)
!            t(i,j,k,1,itemp) = t(i,j,k,1,itemp) + 2.0_r8
!          END DO
!        END DO
!      END DO
#endif

#ifdef ICE_BIO
      DO i=IstrT,IendT
        DO j=JstrT,JendT
          it(i,j,1,iIcePhL) =0.0_r8       !was eps
          it(i,j,1,iIceNO3) =0.0_r8       !was eps
          it(i,j,1,iIceNH4) =0.0_r8       !was eps
          itL(i,j,1,iIceLog) =-1.0_r8     !was eps
        END DO
      END DO
# if defined BERING_10K
      DO i=IstrT,IendT
        DO j=JstrT,JendT  
          IcePhL(i,j,1) =0.0_r8        !was eps
          IceNO3(i,j,1) =0.0_r8        !was eps
          IceNH4(i,j,1) =0.0_r8        !was eps
          IceLog(i,j,1) =-1            !was eps
          IcePhL(i,j,2) =0.0_r8        !was eps
          IceNO3(i,j,2) =0.0_r8        !was eps
          IceNH4(i,j,2) =0.0_r8        !was eps
          IceLog(i,j,2) =-1            !was eps
        END DO
      END DO
# endif
#endif

!reduce the initial concs by 2 orders of magnitude
!test for Bering Sea - initialize in Jan

      DO i=IstrT,IendT
        DO j=JstrT,JendT
          DO k=1,N(ng)
            t(i,j,k,1,iPhS) = 5._r8 !t(i,j,k,1,iPhS)*.01_r8
            t(i,j,k,1,iPhL) = 1._r8 !t(i,j,k,1,iPhL)*.01_r8
            t(i,j,k,1,iMZS) = 0.0_r8
            t(i,j,k,1,iMZL) = 1._r8
            t(i,j,k,1,iCop) = 0.1_r8
!            t(i,j,k,1,iNCaS) = t(i,j,k,1,iNCaS)*.0001_r8
!            t(i,j,k,1,iEupS) = t(i,j,k,1,iEupS)*.0001_r8
!            t(i,j,k,1,iNCaO) = t(i,j,k,1,iNCaO)*.0001_r8
!            t(i,j,k,1,iEupO) = t(i,j,k,1,iEupO)*.0001_r8
          END DO
        END DO
      END DO
!
! Check for size, set other time index, and periodic BC's
!
      DO i=IstrT,IendT
        DO j=JstrT,JendT
          DO k=1,N(ng)
            DO is=1,NBT
              itrc=idbio(is)
              t(i,j,k,1,itrc) = MAX(t(i,j,k,1,itrc),eps)
              t(i,j,k,2,itrc) = t(i,j,k,1,itrc)
            END DO
          END DO
        END DO
      END DO


#ifdef BENTHIC
      DO i=IstrT,IendT
        DO j=JstrT,JendT
          DO k=1,NBL(ng)
            bt(i,j,k,1,iBen) = eps
            bt(i,j,k,1,iBenDet) = eps
            bt(i,j,k,2,iBen) = bt(i,j,k,1,iBen)
            bt(i,j,k,2,iBenDet) = bt(i,j,k,1,iBenDet)
          END DO
        END DO
      END DO
#endif

#ifdef STATIONARY
      DO i=IstrT,IendT
        DO j=JstrT,JendT
          DO k=1,N(ng)
            st(i,j,k,1,1) = 0.0_r8
            st(i,j,k,1,2) = 0.0_r8
            st(i,j,k,1,3) = 0.0_r8
            st(i,j,k,1,4) = 0.0_r8
            st(i,j,k,1,5) = 0.0_r8
            st(i,j,k,1,6) = 0.0_r8
            st(i,j,k,1,7) = 0.0_r8
            st(i,j,k,1,8) = 0.0_r8
            st(i,j,k,2,1) = 0.0_r8
            st(i,j,k,2,2) = 0.0_r8
            st(i,j,k,2,3) = 0.0_r8
            st(i,j,k,2,4) = 0.0_r8
            st(i,j,k,2,5) = 0.0_r8
            st(i,j,k,2,6) = 0.0_r8
            st(i,j,k,2,7) = 0.0_r8
            st(i,j,k,2,8) = 0.0_r8

            st(i,j,k,3,1) = 0.0_r8
            st(i,j,k,3,2) = 0.0_r8
            st(i,j,k,3,3) = 0.0_r8
            st(i,j,k,3,4) = 0.0_r8
            st(i,j,k,3,5) = 0.0_r8
            st(i,j,k,3,6) = 0.0_r8
            st(i,j,k,3,7) = 0.0_r8
            st(i,j,k,3,8) = 0.0_r8
          END DO
        END DO
      END DO
#endif

