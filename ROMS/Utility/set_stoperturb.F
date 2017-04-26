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
      USE mod_ncparam
      USE pseudo_rand2D
      USE mod_parallel
      USE mod_stoperturb


      CONTAINS

      SUBROUTINE get_pseudo2D_forc(ng,ifield,                           &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             rcorr,tcorr,                         &
     &                             dx,dy,RanField,                      &
                                   var)
!
!=======================================================================
! Routine to generate a pseudo random field for forcings
!=======================================================================
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng,ifield             ! grid id, field id
      integer, intent(in) :: LBi, UBi, LBj, UBj    ! Size of tile

      real(r8), intent(in) :: rcorr           ! Horizontal decorrelation length
      real(r8), intent(in) :: tcorr           ! Temporal decorrelation period
      real(r8), intent(in) :: dx,dy           ! horizontal size of input grid
      real(r8), intent(in) :: var             ! desired variance for the perturbation

      real(r8), intent(inout) :: RanField(LBi:UBi,LBj:UBj,2) ! random field on tile
!
!  Local variable declarations.
!
      integer  :: ILB, IUB, JLB, JUB  ! Size of whole ROMS grid
      integer  :: Tindex       ! index of the current field (for 2 field interpolation)
      integer  :: my_tile      ! tile index
      real(r8) :: scal_fac     ! scale factor for current field
      real(r8) :: dt           ! dt between perturbation fields (days)
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 57)
!
!-----------------------------------------------------------------------
!  Get ROMS grid size
!-----------------------------------------------------------------------
!
      my_tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(my_tile)
      IUB=BOUNDS(ng)%UBi(my_tile)
      JLB=BOUNDS(ng)%LBj(my_tile)
      JUB=BOUNDS(ng)%UBj(my_tile)
!
!-----------------------------------------------------------------------
!  Get infos on the forcing variable
!-----------------------------------------------------------------------
!
      scal_fac = Fscale(ifield,ng)
      Tindex   = Iinfo(8,ifield,ng)
      dt = ABS(Vtime(Tindex,ifield,ng) - Vtime(3-Tindex,ifield,ng))
!
!-----------------------------------------------------------------------
!  Get pseudo random 2D field
!-----------------------------------------------------------------------
!
      IF ((dt .EQ. 0.0) .AND. (iic(ng) .EQ. 0)) THEN
         ! Read the perturbation from the restart file on first time
         ! step (iic=0) except for the first job (nrrec=0). In this case
         ! we generate a random field without time dependance

         IF (nrrec(ng) .EQ. 0) THEN
            CALL get_pseudo2D_field_notimedep(LBi, UBi, LBj, UBj,       &
     &               ILB,IUB,JLB,JUB,                                   &
     &               rcorr,                                             &
     &               dx,dy,var,                                         &
     &               RanField(LBi:UBi,LBj:UBj,Tindex),                  &
     &               PERTURB(ng) % plan_fftw)
!
!-----------------------------------------------------------------------
!  Scale field by the variable scaling factor
!-----------------------------------------------------------------------
!
            RanField(LBi:UBi,LBj:UBj,Tindex) =                          &
     &          RanField(LBi:UBi,LBj:UBj,Tindex) * scal_fac 

         END IF
         
      ELSE
         CALL get_pseudo2D_field(LBi, UBi, LBj, UBj,                    &
     &               ILB,IUB,JLB,JUB,                                   &
     &               rcorr,tcorr,                                       &
     &               dx,dy,dt,var,                                      &
     &               RanField(LBi:UBi,LBj:UBj,3-Tindex)/scal_fac,       &
     &               RanField(LBi:UBi,LBj:UBj,Tindex),                  &
     &               PERTURB(ng) % plan_fftw)

!
!-----------------------------------------------------------------------
!  Scale field by the variable scaling factor
!-----------------------------------------------------------------------
!
         RanField(LBi:UBi,LBj:UBj,Tindex) =                             &
     &       RanField(LBi:UBi,LBj:UBj,Tindex) * scal_fac 
      END IF
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 57)

      IF (Master) WRITE(6,*) "RE: (get_pseudo) max( ",                  &
     &       TRIM(Vname(1,ifield)),",t) = ",                            &
     &       MAXVAL(ABS(RanField(LBi:UBi,LBj:UBj,Tindex))),             &
     &       "| dt,iic,ntstart,nrrec = ",dt,iic(ng),ntstart(ng),        &
     &       nrrec(ng)

      IF (Master) WRITE(6,*) "RE: (get_pseudo) max( ",                  &
     &       TRIM(Vname(1,ifield)),",t-1) = ",                          &
     &       MAXVAL(ABS(RanField(LBi:UBi,LBj:UBj,3-Tindex))),           &
     &       "| dt,iic,ntstart,nrrec = ",dt,iic(ng),ntstart(ng),        &
     &       nrrec(ng)

      END SUBROUTINE get_pseudo2D_forc





      SUBROUTINE get_pseudo2D_field(LBi, UBi, LBj, UBj,                 &
     &               ILB,IUB,JLB,JUB,                                   &
     &               rcorr,tcorr,                                       &
     &               dx,dy,dt,var,                                      &
     &               randfield_prev,                                    &
     &               randfield_new, plan)
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
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: ILB, IUB, JLB, JUB

      real(r8), intent(in) :: rcorr           ! Horizontal decorrelation length
      real(r8), intent(in) :: tcorr           ! Temporal decorrelation period
      real(r8), intent(in) :: dx,dy           ! horizontal size of input grid
      real(r8), intent(in) :: dt              ! dt between perturbation fields
      real(r8), intent(in) :: var             ! desired variance for the perturbation

      TYPE(C_PTR),intent(in) :: plan          ! plan for FFTW

      real(r8), intent(in)  :: randfield_prev(LBi:UBi,LBj:UBj)
      real(r8), intent(out) :: randfield_new(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer              :: nx,ny            ! horizontal dimensions
      real(r8)             :: autocorr=0.5     ! autocorrelation for dt=pert_tcorr
      real(r8)             :: alpha            ! coefficient for time dependance
      real(r8),allocatable :: ran_tmp(:,:)     ! Temporary random field
      real(r8) :: RanFieldAll(ILB:IUB,JLB:JUB) ! random field on whole grid
!
!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------
!
      ! Get size of ROMS grid
      nx=IUB-ILB+1
      ny=JUB-JLB+1
      ! Allocate temporary random field at ROMS size
      ALLOCATE(ran_tmp(nx,ny))
!
!-----------------------------------------------------------------------
!  Get the random field
!-----------------------------------------------------------------------
!
      CALL pseudo2D(ran_tmp,nx,ny,rcorr,dx,dy,plan)
      RanFieldAll(ILB:IUB,JLB:JUB) = ran_tmp(:,:)
!
!-----------------------------------------------------------------------
!  Select tile
!-----------------------------------------------------------------------
!
      randfield_new(LBi:UBi,LBj:UBj) = RanFieldAll(LBi:UBi,LBj:UBj)
!
!-----------------------------------------------------------------------
!  Time dependence
!-----------------------------------------------------------------------
!
      alpha = autocorr**(dt/tcorr)
      randfield_new(LBi:UBi,LBj:UBj) =                                  &
     &        alpha*randfield_prev(LBi:UBi,LBj:UBj) / SQRT(var) +       &
     &        sqrt(1-alpha*alpha)*randfield_new(LBi:UBi,LBj:UBj)
!
!-----------------------------------------------------------------------
!  Normalize by desired variance
!-----------------------------------------------------------------------
!
      randfield_new(LBi:UBi,LBj:UBj) =                                  &
     &        randfield_new(LBi:UBi,LBj:UBj) * SQRT(var)

      END SUBROUTINE get_pseudo2D_field




      SUBROUTINE get_pseudo2D_field_notimedep(LBi, UBi, LBj, UBj,       &
     &               ILB,IUB,JLB,JUB,                                   &
     &               rcorr,                                             &
     &               dx,dy,var,                                         &
     &               randfield_new, plan)
!
!=======================================================================
! Main routine that generate the random field for ROMS and
! select the tile corresponding to the processor witout time dependance
!=======================================================================
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: ILB, IUB, JLB, JUB

      real(r8), intent(in) :: rcorr           ! Horizontal decorrelation length
      real(r8), intent(in) :: dx,dy           ! horizontal size of input grid
      real(r8), intent(in) :: var             ! desired variance for the perturbation

      TYPE(C_PTR),intent(in) :: plan          ! plan for FFTW

      real(r8), intent(inout) :: randfield_new(LBi:UBi,LBj:UBj)
!
!  Local variable declarations.
!
      integer              :: nx,ny            ! horizontal dimensions
      real(r8)             :: autocorr=0.5     ! autocorrelation for dt=pert_tcorr
      real(r8)             :: alpha            ! coefficient for time dependance
      real(r8),allocatable :: ran_tmp(:,:)     ! Temporary random field
      real(r8) :: RanFieldAll(ILB:IUB,JLB:JUB) ! random field on whole grid
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      ! Get size of ROMS grid
      nx=IUB-ILB+1
      ny=JUB-JLB+1
      ! Allocate temporary random field at ROMS size
      ALLOCATE(ran_tmp(nx,ny))
!
!-----------------------------------------------------------------------
!  Get the random field
!-----------------------------------------------------------------------
!
      CALL pseudo2D(ran_tmp,nx,ny,rcorr,dx,dy,plan)
      RanFieldAll(ILB:IUB,JLB:JUB) = ran_tmp(:,:)
!
!-----------------------------------------------------------------------
!  Select tile
!-----------------------------------------------------------------------
!
      randfield_new(LBi:UBi,LBj:UBj) = RanFieldAll(LBi:UBi,LBj:UBj)
!
!-----------------------------------------------------------------------
!  Normalize by desired variance
!-----------------------------------------------------------------------
!
      randfield_new(LBi:UBi,LBj:UBj) =                                  &
     &        randfield_new(LBi:UBi,LBj:UBj) * SQRT(var)

      END SUBROUTINE get_pseudo2D_field_notimedep



      END MODULE set_stoperturb

#endif
