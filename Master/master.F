#include "cppdefs.h"
/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2017 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Master program to run ROMS/TOMS as single ocean model or coupled  **
**  to other models using the MCT or ESMF libraries.                  **
**                                                                    **
************************************************************************
*/
#if defined MODEL_COUPLING

      PROGRAM mct_coupler
!
!svn $Id: mct_coupler.h 553 2011-04-22 21:30:04Z arango $
!=======================================================================
!  Copyright (c) 2002-2011 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!==================================================== John C. Warner ===
!                                                                      !
!  Master program to couple ROMS/TOMS to other models using the Model  !
!  Coupling Toolkit (MCT) library.                                     !
!                                                                      !
!  The following models are coupled to ROMS/TOMS:                      !
!                                                                      !
#ifdef WRF_COUPLING
!  WRF, Weather Research and Forecasting model:                        !
!       http://www.wrf-model.org                                       !
!                                                                      !
#endif
#ifdef SWAN_COUPLING
!  SWAN, Simulating WAves Nearshore model:                             !
!        http://vlm089.citg.tudelft.nl/swan/index.htm                  !
!                                                                      !
#endif
#ifdef CICE_COUPLING
!  CICE, The Los Alamos sea ice model:                                 !
!        http://http://oceans11.lanl.gov/drupal/CICE                   !
!                                                                      !
#endif
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_iounits
      USE mod_scalars
!
      USE m_MCTWorld, ONLY : MCTWorld_clean => clean

      USE ocean_control_mod, ONLY : ROMS_initialize
      USE ocean_control_mod, ONLY : ROMS_run
      USE ocean_control_mod, ONLY : ROMS_finalize
#ifdef WRF_COUPLING
      USE ocean_coupler_mod, ONLY : finalize_ocn2atm_coupling
#endif
#if defined SWAN_COUPLING || defined REFDIF_COUPLING
      USE ocean_coupler_mod, ONLY : finalize_ocn2wav_coupling
#endif
#if defined CICE_COUPLING
      USE ocean_coupler_mod, ONLY : finalize_ocn2cice_coupling
      USE ice_communicate, ONLY: MPI_COMM_ICE
      USE CICE_MCT, ONLY: ice_Nmodels => Nmodels
      USE CICE_MCT, ONLY: ice_CICEid => CICEid
      USE CICE_MCT, ONLY: ice_OCNid => OCNid
      USE cice_initmod, ONLY: cice_initialize
      USE cice_runmod, ONLY: cice_run
      USE cice_finalmod, ONLY: cice_finalize
#endif
!
      implicit none
!
!  Local variable declarations.
!
      logical, save :: first

      integer :: MyColor, MyCOMM, MyError, MyKey, Nnodes
      integer :: ng

      real(r4) :: CouplingTime             ! single precision
!
!-----------------------------------------------------------------------
!  Initialize distributed-memory (MPI) configuration
!-----------------------------------------------------------------------
!
!  Initialize MPI execution environment.
!
      CALL mpi_init (MyError)
!
!  Get rank of the local process in the group associated with the
!  comminicator.
!
      CALL mpi_comm_size (MPI_COMM_WORLD, Nnodes, MyError)
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
!
!  Set temporarily the ocean communicator to current handle before
!  splitting so the input coupling script name can be broadcasted to
!  all the nodes.
!
      OCN_COMM_WORLD=MPI_COMM_WORLD
!
!  Read in coupled model parameters from standard input.
!
!
      if (master) WRITE (stdout,*) ' Calling couple_par'
      CALL read_CouplePar (iNLM)
!
!  Allocate several coupling variables.
!
      if (master) WRITE (stdout,*) ' Allocate coupler'
      CALL allocate_coupler (Nnodes)
!
!  Split the communicator into coupled models sub-groups based
!  on color and key.
!
!jd      WRITE (stdout,*) ' Assigning color'
      MyKey=0
      IF ((pets(Iocean)%val(1).le.MyRank).and.                          &
     &    (MyRank.le.pets(Iocean)%val(Nthreads(Iocean)))) THEN
        MyColor=OCNid
      END IF
#ifdef AIR_OCEAN
      IF ((pets(Iatmos)%val(1).le.MyRank).and.                          &
     &    (MyRank.le.pets(Iatmos)%val(Nthreads(Iatmos)))) THEN
        MyColor=ATMid
      END IF
#endif
#ifdef WAVES_OCEAN
      IF ((pets(Iwaves)%val(1).le.MyRank).and.                          &
     &    (MyRank.le.pets(Iwaves)%val(Nthreads(Iwaves)))) THEN
        MyColor=WAVid
      END IF
#endif
#if defined CICE_COUPLING
      IF ((pets(Icice)%val(1).le.MyRank).and.                          &
     &    (MyRank.le.pets(Icice)%val(Nthreads(Icice)))) THEN
        MyColor=CICEid
      END IF
#endif


      CALL mpi_comm_split (MPI_COMM_WORLD, MyColor, MyKey, MyCOMM,      &
     &                     MyError)
!
!-----------------------------------------------------------------------
!  Run coupled models according to the processor rank.
!-----------------------------------------------------------------------
!
#if defined SWAN_COUPLING
      IF (MyColor.eq.WAVid) THEN
        CouplingTime=REAL(TimeInterval(Iocean,Iwaves))
        CALL SWAN_INITIALIZE (MyCOMM, INPname(Iwaves))
        CALL SWAN_RUN (CouplingTime)
        CALL SWAN_FINALIZE
      END IF
#elif defined REFDIF_COUPLING
      IF (MyColor.eq.WAVid) THEN
        CouplingTime=REAL(TimeInterval(Iocean,Iwaves))
        CALL refdif_initialize (MyCOMM)
        CALL refdif_run (CouplingTime, INPname(Iwaves))
        CALL refdif_finalize
      END IF
#endif
#ifdef WRF_COUPLING
      IF (MyColor.eq.ATMid) THEN
        CouplingTime=REAL(TimeInterval(Iocean,Iatmos))
!!      CALL module_wrf_top_mp_wrf_init (MyCOMM)
!!      CALL module_wrf_top_mp_wrf_run (TimeInterval(Iocean,Iwaves))
!!      CALL module_wrf_top_mp_wrf_finalize
        CALL module_wrf_top_wrf_init (MyCOMM)
        CALL module_wrf_top_wrf_run (CouplingTime)
        CALL module_wrf_top_wrf_finalize
      END IF
#endif
#if defined CICE_COUPLING
      IF (MyColor.eq.CICEid) THEN
        WRITE (stdout,*) ' I am CICE node: calling CICE_INITIALIZE: ', MyRank
        MPI_COMM_ICE=MyCOMM
        ice_Nmodels=Nmodels
        ice_CICEid=CICEid
        ice_OCNid=OCNid
        CALL CICE_INITIALIZE
        CALL CICE_RUN(TimeInterval(Iocean,Icice))
        CALL CICE_FINALIZE
      END IF
#endif
      IF (MyColor.eq.OCNid) THEN
        WRITE (stdout,*) ' I am Ocean node: ', MyRank
        first=.TRUE.
        IF (exit_flag.eq.NoError) THEN
!jd        WRITE (stdout,*) ' Calling ROMS_INITIALIZE: ', MyRank
          CALL ROMS_initialize (first, mpiCOMM=MyCOMM)
          run_time=0.0_r8
          DO ng=1,Ngrids
            run_time=MAX(run_time, dt(ng)*ntimes(ng))
          END DO
        ELSE
          CALL mpi_abort(MPI_COMM_WORLD, MyError)
        END IF
        IF (exit_flag.eq.NoError) THEN
!jd          WRITE (stdout,*) ' Calling ROMS_run: ', MyRank
          CALL ROMS_run (run_time)
        ELSE
          CALL mpi_abort(MPI_COMM_WORLD, MyError)
        END IF
        IF (exit_flag.eq.NoError) THEN
           CALL ROMS_finalize
        ELSE
          CALL ROMS_finalize
          CALL mpi_abort(MPI_COMM_WORLD, MyError)
        END IF
#if defined SWAN_COUPLING || defined REFDIF_COUPLING
        CALL finalize_ocn2wav_coupling
#endif
#ifdef WRF_COUPLING
        CALL finalize_ocn2atm_coupling
#endif
#if defined CICE_COUPLING
        CALL finalize_ocn2cice_coupling
#endif
      END IF
!
!-----------------------------------------------------------------------
!  Terminates all the mpi-processing and coupling.
!-----------------------------------------------------------------------
!
      CALL mpi_barrier (MPI_COMM_WORLD, MyError)
      CALL MCTWorld_clean ()
      CALL mpi_finalize (MyError)

      STOP

      END PROGRAM mct_coupler
#else
# include "ocean.h"
#endif
