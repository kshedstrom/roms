!
! svn $Id$
!!================================================== Georgina Gibson ===
!! Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!                                                                      !
!  This routine sets initial conditions for biological tracer fields   !
!  using analytical expressions for line 8 in GOA for GOANPZ.          !
!                                                                      !
!=======================================================================
!

# undef DEPAVG      /* TEST case: set all equal to depth averaged values */

# ifdef DEPAVG
!
! Set all points equal to the same number
!
      do i=IstrT,IendT
        do j=JstrT,JendT
          do k=1,N(ng)
            t(i,j,k,1,iNO3) = 18.0_r8
            t(i,j,k,1,iNH4) = 0.20_r8
            t(i,j,k,1,iPhS) = 0.10_r8
            t(i,j,k,1,iPhL) = 0.10_r8
            t(i,j,k,1,iMZS) = 0.10_r8
            t(i,j,k,1,iMZL) = 0.10_r8
            t(i,j,k,1,iCop) = 2.041_r8
            t(i,j,k,1,iNCa) = 2.417_r8
            t(i,j,k,1,iEup) = 1.337_r8
            t(i,j,k,1,iDet) = eps
#  ifdef IRON_LIMIT
            t(i,j,k,1,iFe) = 0.6_r8
#  endif
          enddo
        enddo
      enddo
# else
!
! Make the curves.  First come parameters that are constant through 
! entire water column.
!
      deepval(iNO3) = 30.0_r8
      deepval(iNH4) = eps
      deepval(iPhS) = eps
      deepval(iPhL) = eps
      deepval(iMZS) = eps
      deepval(iMZL) = eps
      deepval(iCop) = 0.00863_r8
      deepval(iNCa) = eps
      deepval(iEup) = 0.08_r8
      deepval(iDet) = eps
      do k=1,N(ng)
         do j=JstrT,JendT
            do i=IstrT,IendT
               t(i,j,k,1,iNO3) = 18.0_r8
               t(i,j,k,1,iNH4) =  1.0_r8
               t(i,j,k,1,iDet) =  eps
               biod(i,j,k) = -1 * ( z_r(i,j,k) + 2.5_r8 )
            enddo
         enddo
      enddo
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
      do i=IstrT,IendT
        do j=JstrT,JendT
          do k=1,N(ng)
            if ( biod(i,j,k) .le. var5 ) then
              var6 = biod(i,j,k) - var4
              t(i,j,k,1,iPhS) =                                         &
     &          var1 + var2*(var6**2) + var3*(var6**3)
            else
              t(i,j,k,1,iPhS) =                                         &
     &          var7 * exp( ( -1.0_r8*biod(i,j,k) + var5 ) / 5.0_r8 )
            endif
          enddo
        enddo
      enddo
!
! PL - exponentially decreasing with depth.
!
      var1 = 8.25_r8
      var2 = 0.322_r8 / 5.0_r8
      do i=IstrT,IendT
        do j=JstrT,JendT
          do k=1,N(ng)
            t(i,j,k,1,iPhL) = var1 * exp( -1.0_r8 * var2 * biod(i,j,k) )
          enddo
        enddo
      enddo
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
      do i=IstrT,IendT
        do j=JstrT,JendT
          do k=1,N(ng)
            if ( biod(i,j,k) .le. var3 ) then
              t(i,j,k,1,iMZL) = var1 + var2 * biod(i,j,k)
            else
              t(i,j,k,1,iMZL) = var4 +                                  &
     &          ( var5 - var4 ) * exp( ( var3 - biod(i,j,k) ) / 5.0_r8 )
            endif
            t(i,j,k,1,iMZS) = t(i,j,k,1,iMZL)
          enddo
        enddo
      enddo
!
! Pseudocalanus Copepods - from Shelikof data via Shelikof NPZ
! Step at ~ 50 m
!
      var1 = pi / 2.0_r8
      var2 = 2.876_r8 / pi
      var3 = 5.0_r8 / 5.0_r8
      var4 = 9.151_r8 * 5.0_r8
      do i=IstrT,IendT
        do j=JstrT,JendT
          do k=1,N(ng)
            var5 = var3 * ( biod(i,j,k) - var4 )
            t(i,j,k,1,iCop) = var1 - var2 * atan( var5 )
          enddo
        enddo
      enddo
!
! Neocalanus, from Shelikof NPZ. Step at ~ 30m
!
      var1 = 4.0_r8
      var2 = 1.3_r8
      var3 = 5.0_r8 / 5.0_r8
      var4 = 5.2_r8 * 5.0_r8
      do i=IstrT,IendT
        do j=JstrT,JendT
          do k=1,N(ng)
            var5 = var3 * ( biod(i,j,k) - var4 )
            t(i,j,k,1,iNCa) = var1 - var2 * atan( var5 )
          enddo
        enddo
      enddo
!
! Euphausiids, Wild guesses. Step at ~ 30m
!
      var1 = 1.78_r8
      var2 = 0.8_r8
      var3 = 5.0_r8 / 5.0_r8
      var4 = 5.2_r8 * 5.0_r8
      do i=IstrT,IendT
        do j=JstrT,JendT
          do k=1,N(ng)
            var5 = var3 * ( biod(i,j,k) - var4)
            t(i,j,k,1,iEup) =  var1 - var2 * atan( var5 )
          enddo
        enddo
      enddo
!
#  ifdef IRON_LIMIT
! Iron - linear from surface value to value at 100m and increase onshore
      do i=IstrT,IendT
        do j=JstrT,JendT
          var1 = MAX(0._r8,MIN(1._r8,                                   &
     &             (GRID(ng)%h(i,j)-Feinh)/(Feoffh-Feinh)))
          FeSurf = Feinlo + var1*(Feofflo-Feinlo)
          FeDeep = Feinhi + var1*(Feoffhi-Feinhi)
          var1 = (FeDeep-FeSurf) / 100._r8
          do k=1,N(ng)
            t(i,j,k,1,iFe) = MIN(FeDeep, FeSurf - z_r(i,j,k)*var1)
          enddo
        enddo
      enddo
#  endif
!
! Concentrations of everything below 100m - i.e. below
! depths where calculations are performed.  Have linear slope
! between values above and below.
! Iron deep values have already been determined.
!
      do i=IstrT,IendT
        do j=JstrT,JendT
          do k=N(ng),1,-1
            if ( biod(i,j,k) .gt. 120.0_r8 ) then
              t(i,j,k,1,iNO3) = deepval(iNO3)
              t(i,j,k,1,iNH4) = deepval(iNH4)
              t(i,j,k,1,iPhS) = deepval(iPhS)
              t(i,j,k,1,iPhL) = deepval(iPhL)
              t(i,j,k,1,iMZS) = deepval(iMZS)
              t(i,j,k,1,iMZL) = deepval(iMZL)
              t(i,j,k,1,iCop) = deepval(iCop)
              t(i,j,k,1,iNCa) = deepval(iNCa)
              t(i,j,k,1,iEup) = deepval(iEup)
              t(i,j,k,1,iDet) = deepval(iDet)
            else if ( biod(i,j,k) .gt. 100.0_r8 .and.                   &
     &                biod(i,j,k) .le. 120.0_r8) then
              var1 = ( 100.0_r8 - biod(i,j,k) ) / ( 100.0_r8-120.0_r8 )
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
              t(i,j,k,1,iNCa) = loval(iNCa) +                           &
     &                          ( deepval(iNCa) - loval(iNCa) ) * var1
              t(i,j,k,1,iEup) = loval(iEup) +                           &
     &                          ( deepval(iEup) - loval(iEup) ) * var1
              t(i,j,k,1,iDet) = loval(iDet) +                           &
     &                          ( deepval(iDet) - loval(iDet) ) * var1
            else
              loval(iNO3) = t(i,j,k,1,iNO3)
              loval(iNH4) = t(i,j,k,1,iNH4)
              loval(iPhS) = t(i,j,k,1,iPhS)
              loval(iPhL) = t(i,j,k,1,iPhL)
              loval(iMZS) = t(i,j,k,1,iMZS)
              loval(iMZL) = t(i,j,k,1,iMZL)
              loval(iCop) = t(i,j,k,1,iCop)
              loval(iNCa) = t(i,j,k,1,iNCa)
              loval(iEup) = t(i,j,k,1,iEup)
              loval(iDet) = t(i,j,k,1,iDet)
            endif
          enddo
        enddo
      enddo
# endif /* DEPAVG */
#ifdef GAK1D
!
!  This is a hack for sensitivity studies - to test warmer or cooler water
!
!      do i=IstrT,IendT
!        do j=JstrT,JendT
!          do k=1,N(ng)
!            t(i,j,k,1,itemp) = t(i,j,k,1,itemp) + 2.0_r8
!          enddo
!        enddo
!      enddo
#endif
!
! Check for size, set other time index, and periodic BC's
!
      do i=IstrT,IendT
         do j=JstrT,JendT
            do k=1,N(ng)
               DO is=1,NBT
                  itrc=idbio(is)
                  t(i,j,k,1,itrc) = MAX(t(i,j,k,1,itrc),eps)
                  t(i,j,k,2,itrc) = t(i,j,k,1,itrc)
               enddo
            enddo
         enddo
      enddo

