#include "cppdefs.h"
#ifdef STOPERTURB
      MODULE set_stoperturb
!
!=======================================================================
!  Copyright (c) 2002-2014 ROMS/TOMS Group                             !
!================================================== Hernan G. Arango ===
!                                                                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
      USE mod_ncparam
      USE pseudo_rand2D


      CONTAINS

      SUBROUTINE get_pseudo2D_field(ng,ifield,               &
                                    LBi, UBi, LBj, UBj,      &
                                    ILB,IUB,JLB,JUB,         &
                                    Iout, Irec,              &
                                    rcorr,tcorr,             &
                                    dx,dy,RanField,          &
                                    var)
!
!=======================================================================
! Main routine that generate the random field for ROMS and
! select the tile corresponding to the processor
!=======================================================================
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng,ifield,Iout, Irec
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: ILB, IUB, JLB, JUB

      real(r8), intent(in) :: rcorr           ! Horizontal decorrelation length
      real(r8), intent(in) :: tcorr           ! Temporal decorrelation period
      real(r8), intent(in) :: dx,dy           ! horizontal size of input grid
      real(r8), intent(in) :: var             ! desired variance for the perturbation

      real(r8), intent(inout) :: RanField(LBi:UBi,LBj:UBj,Iout) ! random field on tile
!
!  Local variable declarations.
!
      integer              :: nx,ny            ! horizontal dimensions
      integer              :: Tindex           ! index of the current field (for 2 field interpolation)
      real(r8)             :: autocorr=0.75    ! autocorrelation for dt=pert_tcorr
      real(r8)             :: alpha            ! coefficient for time dependance
      real(r8)             :: dt               ! dt between perturbation fields
      real(r8),allocatable :: ran_tmp(:,:)     ! Temporary random field
      real(r8) :: RanFieldAll(ILB:IUB,JLB:JUB) ! random field on whole grid
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      Tindex=Iinfo(8,ifield,ng)
      nx=IUB-ILB+1
      ny=JUB-JLB+1
      ALLOCATE(ran_tmp(nx,ny))
!
!-----------------------------------------------------------------------
!  Get the random field
!-----------------------------------------------------------------------
!
      CALL pseudo2D(ran_tmp,nx,ny,rcorr,dx,dy)
      RanFieldAll(ILB:IUB,JLB:JUB) = ran_tmp(:,:)
!
!-----------------------------------------------------------------------
!  Select tile
!-----------------------------------------------------------------------
!
      RanField(LBi:UBi,LBj:UBj,Tindex) = RanFieldAll(LBi:UBi,LBj:UBj)
!
!-----------------------------------------------------------------------
!  Time dependence
!-----------------------------------------------------------------------
!
      dt = ABS(Vtime(Tindex,ifield,ng) - Vtime(3-Tindex,ifield,ng))
      WRITE(6,*) 'RE: dt is :', dt
      alpha = autocorr**(dt/tcorr)

      RanField(LBi:UBi,LBj:UBj,Tindex) = &
        alpha*RanField(LBi:UBi,LBj:UBj,3-Tindex) &
        + sqrt(1-alpha*alpha)*RanField(LBi:UBi,LBj:UBj,Tindex)
!
!-----------------------------------------------------------------------
!  Normalize by desired variance
!-----------------------------------------------------------------------
!
      RanField(LBi:UBi,LBj:UBj,Tindex) = &
        RanField(LBi:UBi,LBj:UBj,Tindex) * Fscale(ifield,ng)*SQRT(var)


      END SUBROUTINE get_pseudo2D_field

      END MODULE set_stoperturb

#endif
