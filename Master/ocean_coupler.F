#include "cppdefs.h"

!#define REPORT_COUPLING_ALL

      MODULE ocean_coupler_mod

#if defined MODEL_COUPLING && defined MCT_LIB
!
!svn $Id$
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module is used to communicate and exchange data between        !
!  ROMS/TOMS and other coupled model(s)  via the Model Coupling        !
!  Toolkit (MCT), developed at the Argonne National Laboratory.        !
!                                                                      !
!=======================================================================
!
!  Component Model Registry.
!
      USE m_MCTWorld, ONLY : MCTWorld_init => init
      USE m_MCTWorld, ONLY : MCTWorld_clean => clean
!
!  Domain Decomposition Descriptor DataType and associated methods.
!
      USE m_GlobalSegMap, ONLY : GlobalSegMap
      USE m_GlobalSegMap, ONLY : GlobalSegMap_init => init
      USE m_GlobalSegMap, ONLY : GlobalSegMap_lsize => lsize
      USE m_GlobalSegMap, ONLY : GlobalSegMap_clean => clean
      USE m_GlobalSegMap, ONLY : GlobalSegMap_Ordpnts => OrderedPoints
!
!  Field Storage DataType and associated methods.
!
      USE m_AttrVect, ONLY : AttrVect
      USE m_AttrVect, ONLY : AttrVect_init => init
      USE m_AttrVect, ONLY : AttrVect_zero => zero
      USE m_AttrVect, ONLY : AttrVect_lsize => lsize
      USE m_AttrVect, ONLY : AttrVect_clean => clean
      USE m_AttrVect, ONLY : AttrVect_importRAttr => importRAttr
      USE m_AttrVect, ONLY : AttrVect_exportRAttr => exportRAttr
!
!  Intercomponent communications scheduler.
!
      USE m_Router, ONLY : Router
      USE m_Router, ONLY : Router_init => init
      USE m_Router, ONLY : Router_clean => clean
!
!  Intercomponent transfer.
!
      USE m_Transfer, ONLY: MCT_Send => send
      USE m_Transfer, ONLY: MCT_Recv => recv
!
!  Sparse Matrix DataType and associated methods.
!
      USE m_SparseMatrix, ONLY : SparseMatrix
      USE m_SparseMatrix, ONLY : SparseMatrix_init => init
      USE m_SparseMatrix, ONLY : SparseMatrix_importGRowInd =>          &
     &                           importGlobalRowIndices
      USE m_SparseMatrix, ONLY : SparseMatrix_importGColInd =>          &
     &                           importGlobalColumnIndices
      USE m_SparseMatrix, ONLY : SparseMatrix_importMatrixElts =>       &
     &                           importMatrixElements
      USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus
      USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus_init => init
      USE m_SparseMatrixPlus, ONLY : SparseMatrixPlus_clean => clean
!
!  Decompose matrix by row.
!
      USE m_SparseMatrixPlus, ONLY : Xonly
!
!  Matrix-Vector multiply methods.
!
      USE m_MatAttrVectMul, ONLY : MCT_MatVecMul => sMatAvMult
!
      implicit none
!
      PRIVATE

# ifdef SWAN_COUPLING
      PUBLIC :: initialize_ocn2wav_coupling
      PUBLIC :: ocn2wav_coupling
      PUBLIC :: finalize_ocn2wav_coupling
# endif
# ifdef WRF_COUPLING
      PUBLIC :: initialize_ocn2atm_coupling
      PUBLIC :: ocn2atm_coupling
      PUBLIC :: finalize_ocn2atm_coupling
# endif
# ifdef CICE_COUPLING
      PUBLIC :: initialize_ocn2cice_coupling
      PUBLIC :: ocn2cice_coupling
      PUBLIC :: finalize_ocn2cice_coupling
# endif
!
!  Declarations.
!
# ifdef SWAN_COUPLING
      TYPE(AttrVect) :: wav2ocn_AV            ! AttrVect variables
      TYPE(AttrVect) :: ocn2wav_AV
      TYPE(Router)   :: ROMStoSWAN            ! Router variables
# endif
# ifdef WRF_COUPLING
      TYPE(AttrVect) :: atm2ocn_AV            ! AttrVect variables
      TYPE(AttrVect) :: ocn2atm_AV
      TYPE(Router)   :: ROMStoWRF             ! Router variables
# endif
# ifdef CICE_COUPLING
      TYPE(AttrVect) :: cice2ocn_AV           ! AttrVect variables
      TYPE(AttrVect) :: ocn2cice_AV
      TYPE(Router)   :: ROMStoCICE            ! Router variables
# endif
      TYPE(GlobalSegMap) :: GSMapROMS         ! GloabalSegMap variables

      CONTAINS
!************************************************************************
!*  Include model specific communication routines.
!************************************************************************

# ifdef SWAN_COUPLING
#  include "mct_roms_swan.h"
# endif
# ifdef REFDIF_COUPLING
#  include "mct_roms_refdif.h"
# endif
# ifdef WRF_COUPLING
#  include "mct_roms_wrf.h"
# endif
# ifdef CICE_COUPLING

      SUBROUTINE initialize_ocn2cice_coupling (ng, tile)
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_kinds
      USE mod_scalars
      USE mod_iounits, only: stdout

!
!  Imported variable definitions.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Istr, Iend, Jstr, Jend
      integer :: Asize, Jsize, MyError
      integer :: j, jc, nprocs

      integer, allocatable :: length(:)
      integer, allocatable :: start(:)

#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 36)
#endif

!
!-----------------------------------------------------------------------
!  Begin initialization phase.
!-----------------------------------------------------------------------
!
!  Get communicator local rank and size.
!
#ifdef REPORT_COUPLING_ALL
      WRITE (stdout,*) ' ROMS Ocean Coupler: mpi_comm_rank'
#endif
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
#ifdef REPORT_COUPLING_ALL
      WRITE (stdout,*) ' ROMS Ocean Coupler: mpi_comm_size'
#endif
      CALL mpi_comm_size (OCN_COMM_WORLD, nprocs, MyError)
!
!  Initialize MCT coupled model registry.
!
#ifdef REPORT_COUPLING_ALL
      WRITE (stdout,*) ' ROMS Ocean Coupler: calling MCTWorld_init'
#endif
      CALL MCTWorld_init (Nmodels, MPI_COMM_WORLD, OCN_COMM_WORLD,      &
     &                    OCNid)
#ifdef REPORT_COUPLING_ALL
      WRITE(stdout,*) ' Done: ROMS Ocean Coupler: calling MCTWorld_init'
#else
      IF (Master)                                                       &
     &     WRITE(stdout,*)' ROMS Ocean Coupler: MCTWorld_init called'
#endif
!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables. Notice that "set_bounds.h" is
!  not used here because of implementation of periodicity in other
!  models.
!-----------------------------------------------------------------------
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
!  Determine start and lengths for domain decomposition.
!
      Jsize=Jend-Jstr+1
#ifdef REPORT_COUPLING_ALL
      WRITE (stdout,*) ' ROMS Ocean Coupler: Jsize=', Jsize
#endif
      IF (.not.allocated(start)) THEN
        allocate ( start(Jsize) )
      END IF
      IF (.not.allocated(length)) THEN
        allocate ( length(Jsize) )
      END IF
      jc=0
      DO j=Jstr,Jend
        jc=jc+1
        start (jc)=(j-1)*Lm(ng)+Istr
        length(jc)=(Iend-Istr+1)
      END DO
#ifdef REPORT_COUPLING_ALL
      WRITE (stdout,*) ' ROMS Ocean Coupler: calling GlobalSegMap_init'
#else
      if (Master) WRITE (stdout,*)                                      &
     & ' ROMS Ocean Coupler: calling GlobalSegMap_init'
#endif
      CALL GlobalSegMap_init (GSMapROMS, start, length, 0,              &
     &     OCN_COMM_WORLD, OCNid)
!
!  Initialize attribute vector holding the export data code strings of
!  the cice model. The Asize is the number of grid point on this
!  processor.
!
      Asize=GlobalSegMap_lsize(GSMapROMS, OCN_COMM_WORLD)

#ifdef REPORT_COUPLING_ALL
      WRITE (stdout,*) ' ROMS Ocean Coupler: calling AttrVect_init.',   &
     &     ' Asize=', Asize
#endif

      if (master) write(stdout,*)                                       &
     &     'cice2ocn_AV, Exportlist ',trim(ExportList(Icice))

      CALL AttrVect_init (cice2ocn_AV, rList=TRIM(ExportList(Icice)),   &
     &     lsize=Asize)
!
!  Initialize attribute vector holding the export data code string of
!  the ocean model.
!
      if (master) write(stdout,*)                                       &
     &     'ocn2cice_AV, Exportlist ',trim(ExportList(Iocean))
      CALL AttrVect_init (ocn2cice_AV, rList=TRIM(ExportList(Iocean)),  &
     &     lsize=Asize)
      CALL AttrVect_zero (ocn2cice_AV)

#ifdef REPORT_COUPLING_ALL
      WRITE (stdout,*) ' ROMS Ocean Coupler: calling Router_init'
#endif
      CALL Router_init (CICEid, GSMapROMS, OCN_COMM_WORLD, ROMStoCICE)

#ifdef REPORT_COUPLING_ALL
      WRITE (stdout,*) ' ROMS Ocean Coupler: calling Router_init. Done'
#else
      if (master) write(stdout,*)                                       &
     &     ' ROMS Ocean Coupler: calling Router_init. Done'
#endif
!
!  Deallocate working arrays.
!
      IF (allocated(start)) THEN
        deallocate (start)
      END IF
      IF (allocated(length)) THEN
        deallocate (length)
      END IF

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 36)
#endif

      RETURN
      END SUBROUTINE initialize_ocn2cice_coupling

      SUBROUTINE ocn2cice_coupling (ng, tile, ncouple)
      USE mod_param
      USE mod_parallel
      USE mod_coupler
      USE mod_forces
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_iounits, only: stdout
      USE mod_grid
!      USE distribute_mod, ONLY : mp_reduce
      USE ROMS_import_mod, ONLY : ROMS_import2d
      USE ROMS_export_mod, ONLY : ROMS_export2d
#  ifdef LMD_SKPP
      USE mod_mixing
#  endif
      USE bc_2d_mod
      USE exchange_2d_mod
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
      USE mod_ice
!jd
      USE frazil_ice_prod_mod, only : t_freeze
!jd

      implicit none

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, ncouple
!
!  Local variable declarations.
!
      integer :: Istr, Iend, Jstr, Jend
      integer :: Asize, Iimport, Iexport, MyError
      integer :: gtype, i, id, ifield, ij, j, status

      real(r8) :: add_offset, scale

      real(r8), pointer :: A(:)

      ! intermediates for ice stress interpolation
      real(r8) :: aice_psip, aice_psi, sap, sa

      character (len=3 ), dimension(2) :: op_handle
      character (len=40) :: code


      integer :: k
      integer, dimension(:,:), allocatable :: nbot
      real(r8), dimension(:,:), allocatable :: uw
      real(r8), dimension(:,:), allocatable :: vw
      real(r8), dimension(:,:), allocatable :: frzmlt

      real(r8) :: Sold, t_fr, dz, mlio, totml, dml

#include "tile.h"

#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 36)
#endif

      allocate(nbot(LBi:UBi,LBj:UBj))
      allocate(uw(LBi:UBi,LBj:UBj))
      allocate(vw(LBi:UBi,LBj:UBj))
      allocate(frzmlt(LBi:UBi,LBj:UBj))

      uw=0.0_r8
      vw=0.0_r8
      frzmlt=0.0_r8

      IF (Master) THEN
          write(stdout,*) ''
          write(stdout,*) '  ***************************************** '
          write(stdout,*)                                               &
     &     '    Ocean - CICE: coupling routine called from ROMS'
          write(stdout,*) &
          '    Time : ', time_code, time
          write(stdout,*) '  ****************************************  '

          write(stdout,*) ''
      END IF

!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables. Notice that "set_bounds.h" is
!  not used here because of implementation of periodicity in other
!  models.
!-----------------------------------------------------------------------
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
!-----------------------------------------------------------------------
!  Allocate communications array.
!-----------------------------------------------------------------------
!
      Asize=GlobalSegMap_lsize(GSMapROMS, OCN_COMM_WORLD)
      allocate(A(Asize))
      A=0.0_r8

!
!  Schedule receiving fields from cice model.
!
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)

      CALL MCT_Recv (cice2ocn_AV, ROMStoCICE)
!
!  Receive fields from cice model.
!
      Iimport=0
      DO ifield=1,Nimport(Iocean)
        id=ImportID(Iocean)%val(ifield)
        code=ADJUSTL(Fields(id)%code)
        gtype=Fields(id)%GridType
        scale=Fields(id)%scale
        add_offset=Fields(id)%AddOffset

        SELECT CASE (TRIM(code))



        CASE ('AICE')           ! Sea ice concentration

           CALL AttrVect_exportRAttr (cice2ocn_AV, TRIM(code),A,Asize)
           Iimport=Iimport+1
#ifdef REPORT_COUPLING_ALL
           write(stdout,*) 'ROMS rank ', MyRank,                        &
     &          ' received ',trim(code),' field (max/min): ',           &
     &           maxval(A),' ',minval(A),' gtype : ', gtype
#endif
           scale=1.0_r8
           add_offset=0.0_r8
           CALL ROMS_import2d (ng, tile,                                &
     &                         id, gtype, scale, add_offset,            &
     &                         Asize, A,                                &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Fields(id)%ImpMin, Fields(id)%ImpMax,    &
     &                         ICE(ng)%aice,                            &
     &                         status)

!jd Transfere to u- and v- points for weigthing with stress

           DO j=Jstr,Jend
              DO i=Istr,Iend
                 ICE(ng)%aice_u(i,j)=                                   &
     &                0.5_r8*(ICE(ng)%aice(i-1,j)+ICE(ng)%aice(i,j))
!# ifdef MASKING
!                 ICE(ng)%aice_u(i,j)=ICE(ng)%aice_u(i,j)*umask(i,j)
!# endif
              END DO
           END DO
           DO j=Jstr,Jend
              DO i=Istr,Iend
                 ICE(ng)%aice_v(i,j)=                                   &
     &                0.5_r8*(ICE(ng)%aice(i,j-1)+ICE(ng)%aice(i,j))
!# ifdef MASKING
!                 ICE(ng)%aice_v(i,j)=ICE(ng)%aice_v(i,j)*vmask(i,j)
!# endif
              END DO
           END DO

           IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &             LBi, UBi, LBj, UBj,                                  &
     &             ICE(ng)%aice_u)
              CALL exchange_v2d_tile (ng, tile,                         &
     &             LBi, UBi, LBj, UBj,                                  &
     &             ICE(ng)%aice_v)
           END IF

#ifdef DISTRIBUTE
           CALL mp_exchange2d (ng, tile, iNLM, 2,                       &
     &          LBi, UBi, LBj, UBj,                                     &
     &          NghostPoints,                                           &
     &          EWperiodic(ng), NSperiodic(ng),                         &
     &          ICE(ng)%aice_u,ICE(ng)%aice_v )
#endif


        CASE ('freshAI')        ! Fresh water flux

           CALL AttrVect_exportRAttr (cice2ocn_AV, TRIM(code), A, Asize)
           Iimport=Iimport+1
#ifdef REPORT_COUPLING_ALL
           write(stdout,*) 'ROMS rank ', MyRank,                        &
     &          ' received ',trim(code),' field (max/min): ',           &
     &          maxval(A),' ',minval(A),' gtype : ', gtype
#endif
           scale=1.0_r8
           add_offset=0.0_r8
           CALL ROMS_import2d (ng, tile,                                &
     &                         id, gtype, scale, add_offset,            &
     &                         Asize, A,                                &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Fields(id)%ImpMin, Fields(id)%ImpMax,    &
     &                         ICE(ng)%freshAI,                         &
     &                         status)

        CASE ('fsaltAI')        ! Salt flux

           CALL AttrVect_exportRAttr (cice2ocn_AV, TRIM(code), A, Asize)
           Iimport=Iimport+1
#ifdef REPORT_COUPLING_ALL
           write(stdout,*) 'ROMS rank ', MyRank,                        &
     &          ' received ',trim(code),' field (max/min): ',           &
     &          maxval(A),' ',minval(A),' gtype : ', gtype
#endif
           scale=1.0_r8
           add_offset=0.0_r8
           CALL ROMS_import2d (ng, tile,                                &
     &                         id, gtype, scale, add_offset,            &
     &                         Asize, A,                                &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Fields(id)%ImpMin, Fields(id)%ImpMax,    &
     &                         ICE(ng)%fsaltAI,                         &
     &                         status)

        CASE ('fhocnAI')        ! Nonradiative heat flux

           CALL AttrVect_exportRAttr(cice2ocn_AV, TRIM(code), A, Asize)
           Iimport=Iimport+1
#ifdef REPORT_COUPLING_ALL
           write(stdout,*) 'ROMS rank ', MyRank,                        &
     &          ' received ',trim(code),' field (max/min): ',           &
     &          maxval(A),' ',minval(A),' gtype : ', gtype
#endif
           scale=1.0_r8
           add_offset=0.0_r8
           CALL ROMS_import2d (ng, tile,                                &
     &                         id, gtype, scale, add_offset,            &
     &                         Asize, A,                                &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Fields(id)%ImpMin, Fields(id)%ImpMax,    &
     &                         ICE(ng)%fhocnAI,                         &
     &                         status)

        CASE ('fswthruAI')      ! Radiative heat flux through ice

           CALL AttrVect_exportRAttr (cice2ocn_AV, TRIM(code), A, Asize)
           Iimport=Iimport+1
#ifdef REPORT_COUPLING_ALL
           write(stdout,*) 'ROMS rank ', MyRank,                        &
     &          ' received ',trim(code),' field (max/min): ',           &
     &          maxval(A),' ',minval(A),' gtype : ', gtype
#endif
           scale=1.0_r8
           add_offset=0.0_r8
           CALL ROMS_import2d (ng, tile,                                &
     &                         id, gtype, scale, add_offset,            &
     &                         Asize, A,                                &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Fields(id)%ImpMin, Fields(id)%ImpMax,    &
     &                         ICE(ng)%fswthruAI,                       &
     &                         status)

        CASE ('strocnx')        ! ice-ocean stress x-dir

           CALL AttrVect_exportRAttr (cice2ocn_AV, TRIM(code),A,Asize)
           Iimport=Iimport+1
#ifdef REPORT_COUPLING_ALL
           write(stdout,*) 'ROMS rank ', MyRank,                        &
     &          ' received strocnx field (max/min): ',                  &
     &          maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
           scale=1.0_r8
           add_offset=0.0_r8

           CALL ROMS_import2d (ng, tile,                                &
     &                         id, gtype, scale, add_offset,            &
     &                         Asize, A,                                &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Fields(id)%ImpMin, Fields(id)%ImpMax,    &
     &                         ICE(ng)%strx,                            &
     &                         status)


           DO j=Jstr,Jend
              DO i=Istr,Iend
                 ICE(ng)%stru(i,j) = 0.5_r8*                            &
     &                (ICE(ng)%strx(i-1,j) + ICE(ng)%strx(i,j))

!# ifdef MASKING
!                 ICE(ng)%stru(i,j)=ICE(ng)%stru(i,j)*umask(i,j)
!# endif
              END DO
           END DO

           IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &             LBi, UBi, LBj, UBj,                                  &
     &             ICE(ng)%stru)
           END IF
#ifdef DISTRIBUTE
           CALL mp_exchange2d (ng, tile, iNLM, 1,                       &
     &          LBi, UBi, LBj, UBj,                                     &
     &          NghostPoints,                                           &
     &          EWperiodic(ng), NSperiodic(ng),                         &
     &          ICE(ng)%stru )
#endif


        CASE ('strocny')        ! ice-ocean stress y-dir

           CALL AttrVect_exportRAttr (cice2ocn_AV, TRIM(code), A, Asize)
           Iimport=Iimport+1
#ifdef REPORT_COUPLING_ALL
           write(stdout,*) 'ROMS rank ', MyRank,                        &
     &          ' received strocny field (max/min): ',                  &
     &          maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
           scale=1.0_r8
           add_offset=0.0_r8
           CALL ROMS_import2d (ng, tile,                                &
     &                         id, gtype, scale, add_offset,            &
     &                         Asize, A,                                &
     &                         Istr, Iend, Jstr, Jend,                  &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         Fields(id)%ImpMin, Fields(id)%ImpMax,    &
     &                         ICE(ng)%stry,                            &
     &                         status)

           DO j=Jstr,Jend
              DO i=Istr,Iend
                 ICE(ng)%strv(i,j)=0.5_r8*                              &
     &                (ICE(ng)%stry(i,j-1)+ICE(ng)%stry(i,j))
!# ifdef MASKING
!                 ICE(ng)%strv(i,j)=ICE(ng)%strv(i,j)*vmask(i,j)
!# endif
              END DO
           END DO
           IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &             LBi, UBi, LBj, UBj,                                  &
     &             ICE(ng)%strv)
           END IF
#ifdef DISTRIBUTE
           CALL mp_exchange2d (ng, tile, iNLM, 1,                       &
     &          LBi, UBi, LBj, UBj,                                     &
     &          NghostPoints,                                           &
     &          EWperiodic(ng), NSperiodic(ng),                         &
     &          ICE(ng)%strv )
#endif

        END SELECT
      END DO

!
!-----------------------------------------------------------------------
!  Export fields from ocean (ROMS) to sea ice (CICE) model.
!-----------------------------------------------------------------------
!

!  Prepare for depth-averaging.

      do j=Jstr,Jend
         do i=Istr,Iend
!           mlio = min(MIXING(ng)%hsbl(i,j),-5.0r8)
            mlio = -5.0_r8
            nbot(i,j) = 1
            do k=N(ng),1,-1
               if(GRID(ng)%z_r(i,j,k).lt.mlio) then
                  nbot(i,j) = min(k,N(ng))
                  nbot(i,j) = max(nbot(i,j),1)
                  goto 1111
               endif
            enddo
 1111       continue
         enddo
      enddo

!  Schedule sending fields to the cice model.
!
      Iexport=0
      DO ifield=1,Nexport(Iocean)
         id=ExportID(Iocean)%val(ifield)
         code=ADJUSTL(Fields(id)%code)
         gtype=Fields(id)%GridType
         scale=Fields(id)%scale
         add_offset=Fields(id)%AddOffset

         SELECT CASE (TRIM(code))

         CASE ('SST')

            CALL ROMS_export2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(ng)%t(:,:,N(ng),NOUT,itemp),      &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
#ifdef REPORT_COUPLING_ALL
            write(stdout,*)                                             &
     &           'ROMS rank ', MyRank, ' sending sst field (max/min): ',&
     &           maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
            CALL AttrVect_importRAttr(ocn2cice_AV, TRIM(code), A, Asize)
            Iexport=Iexport+1

         CASE ('SSS')

            CALL ROMS_export2d (ng, tile,                               &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(ng)%t(:,:,N(ng),NOUT,isalt),      &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
#ifdef REPORT_COUPLING_ALL
            write(stdout,*) &
     &           'ROMS rank ', MyRank, ' sending sss field (max/min): ',&
     &           maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
            CALL AttrVect_importRAttr(ocn2cice_AV, TRIM(code), A, Asize)
            Iexport=Iexport+1

         CASE ('FRZMLT')
           do j=Jstr,Jend
            do i=Istr,Iend
# ifdef ICESHELF
             ! Don't touch the ice shelf cavities
             if (GRID(ng)%zice(i,j) .eq. 0.0_r8) then
# endif
             if (ICE(ng)%qfraz_accum(i,j).gt.0.0_r8) then
              ! Ice was formed
              ! Convert from K*m to W/m^2
              frzmlt(i,j) = ICE(ng)%qfraz_accum(i,j)*rho0*Cp/           &
     &                      (ncouple*dt(ng))
             else
              ! No ice was formed; recompute melt potential in surface
              Sold = max(0.0_r8,OCEAN(ng)%t(i,j,N(ng),NOUT,isalt))
              t_fr = t_freeze(Sold,GRID(ng)%z_r(i,j,N(ng)))
              dz = GRID(ng)%Hz(i,j,N(ng))
              frzmlt(i,j) = min(0.0_r8,                                 &
     &                   t_fr - OCEAN(ng)%t(i,j,N(ng),NOUT,itemp))*     &
     &                   dz*rho0*Cp/(ncouple*dt(ng))
             end if
# ifdef ICESHELF
             else
               frzmlt(i,j) = 0.0_r8
             end if
# endif
            end do
           end do

           CALL ROMS_export2d (ng, tile,                                &
     &                         id, gtype, scale, add_offset,            &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         frzmlt,                                  &
     &                         Fields(id)%ExpMin, Fields(id)%ExpMax,    &
     &                         Asize, A,                                &
     &                         status)
#ifdef REPORT_COUPLING_ALL
           write(stdout,*) &
     &          'ROMS rank ',MyRank,' sending frzmlt field (max/min):', &
     &          maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
           CALL AttrVect_importRAttr(ocn2cice_AV, TRIM(code), A, Asize)
           Iexport=Iexport+1

!jd Reset accumulation array
           ICE(ng)%qfraz_accum(:,:) = 0

! KSH: this is a bit of a mess in that we average to rho points, then
! in CICE we average this to psi points.
        CASE ('u')
           do j=Jstr,Jend
              do i=Istr,Iend
                 uw(i,j) = 0.0_r8
                 totml = 0.0_r8
                 do k=N(ng),nbot(i,j),-1
                    dml = GRID(ng)%z_w(i,j,k)-GRID(ng)%z_w(i,j,k-1)
                    uw(i,j) = uw(i,j) + 0.5_r8*dml*                     &
     &                        (OCEAN(ng)%u(i+1,j,k,NOUT) +          &
     &                         OCEAN(ng)%u(i,j,k,NOUT))
                    totml = totml + dml
                 enddo
                 uw(i,j) = uw(i,j)/totml
              enddo
           enddo
!  Apply boundary conditions.
           CALL bc_r2d_tile (ng, tile,                                  &
     &          LBi, UBi, LBj, UBj,                                     &
     &          uw)
#ifdef DISTRIBUTE
           CALL mp_exchange2d (ng, tile, iNLM, 1,                       &
     &          LBi, UBi, LBj, UBj,                                     &
     &          NghostPoints, EWperiodic(ng), NSperiodic(ng),           &
     &          uw)
#endif
!            write(stdout,*) TRIM(code),id,gtype
           CALL ROMS_export2d (ng, tile,                                &
     &                         id, gtype, scale, add_offset,            &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         uw,                                      &
     &                         Fields(id)%ExpMin, Fields(id)%ExpMax,    &
     &                         Asize, A,                                &
     &                         status)

#ifdef REPORT_COUPLING_ALL
           write(stdout,*)                                              &
     &          'ROMS rank ',MyRank,' sending u field (max/min):',      &
     &          maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
           CALL AttrVect_importRAttr(ocn2cice_AV, TRIM(code), A, Asize)
           Iexport=Iexport+1

        CASE ('v')
           do j=Jstr,Jend
              do i=Istr,Iend
                 vw(i,j) = 0.0_r8
                 totml = 0.0_r8
                 do k=N(ng),nbot(i,j),-1
                    dml = GRID(ng)%z_w(i,j,k)-GRID(ng)%z_w(i,j,k-1)
                    vw(i,j) = vw(i,j) + 0.5_r8*dml*                     &
     &                        (OCEAN(ng)%v(i,j+1,k,NOUT) +          &
     &                         OCEAN(ng)%v(i,j,k,NOUT))
                    totml = totml + dml
                 enddo
                 vw(i,j) = vw(i,j)/totml
              enddo
           enddo
!  Apply boundary conditions.
             CALL bc_r2d_tile (ng, tile,                                &
     &            LBi, UBi, LBj, UBj,                                   &
     &            vw)
#ifdef DISTRIBUTE
             CALL mp_exchange2d (ng, tile, iNLM, 1,                     &
     &            LBi, UBi, LBj, UBj,                                   &
     &            NghostPoints, EWperiodic(ng), NSperiodic(ng),         &
     &            vw)
#endif
             CALL ROMS_export2d (ng, tile,                              &
     &                           id, gtype, scale, add_offset,          &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           vw,                                    &
     &                           Fields(id)%ExpMin, Fields(id)%ExpMax,  &
     &                           Asize, A,                              &
     &                           status)
#ifdef REPORT_COUPLING_ALL
             write(stdout,*)                                            &
     &            'ROMS rank ',MyRank,' sending v field (max/min):',    &
     &            maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
             CALL AttrVect_importRAttr(ocn2cice_AV, TRIM(code),A,Asize)
             Iexport=Iexport+1

          CASE ('SSH')

             CALL ROMS_export2d (ng, tile,                              &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          OCEAN(ng)%zeta(:,:,KOUT),               &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
#ifdef REPORT_COUPLING_ALL
             write(stdout,*)                                            &
     &           'ROMS rank ',MyRank,' sending SSH field (max/min):',   &
     &            maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
             CALL AttrVect_importRAttr(ocn2cice_AV,TRIM(code),A,Asize)
             Iexport=Iexport+1

          CASE ('TAIR')

           do j=Jstr,Jend
              do i=Istr,Iend
                 vw(i,j) = FORCES(ng)%Tair(i,j) + 273.16
              enddo
           enddo
!  Apply boundary conditions.
           CALL bc_r2d_tile (ng, tile,                                  &
     &          LBi, UBi, LBj, UBj,                                     &
     &          vw)
#ifdef DISTRIBUTE
           CALL mp_exchange2d (ng, tile, iNLM, 1,                       &
     &          LBi, UBi, LBj, UBj,                                     &
     &          NghostPoints, EWperiodic(ng), NSperiodic(ng),           &
     &          vw)
#endif
             CALL ROMS_export2d (ng, tile,                              &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vw,                                     &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
#ifdef REPORT_COUPLING_ALL
             write(stdout,*)                                            &
     &           'ROMS rank ',MyRank,' sending Tair field (max/min):',  &
     &            maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
             CALL AttrVect_importRAttr(ocn2cice_AV,TRIM(code),A,Asize)
             Iexport=Iexport+1

          CASE ('QAIR')

             CALL ROMS_export2d (ng, tile,                              &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%Hair,                        &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
#ifdef REPORT_COUPLING_ALL
             write(stdout,*)                                            &
     &           'ROMS rank ',MyRank,' sending Qair field (max/min):',  &
     &            maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
             CALL AttrVect_importRAttr(ocn2cice_AV,TRIM(code),A,Asize)
             Iexport=Iexport+1

          CASE ('UWIND')

             CALL ROMS_export2d (ng, tile,                              &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%Uwind,                       &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
#ifdef REPORT_COUPLING_ALL
             write(stdout,*)                                            &
     &           'ROMS rank ',MyRank,' sending Uwind field (max/min):', &
     &            maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
             CALL AttrVect_importRAttr(ocn2cice_AV,TRIM(code),A,Asize)
             Iexport=Iexport+1

          CASE ('VWIND')

             CALL ROMS_export2d (ng, tile,                              &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%Vwind,                       &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
#ifdef REPORT_COUPLING_ALL
             write(stdout,*)                                            &
     &           'ROMS rank ',MyRank,' sending Vwind field (max/min):', &
     &            maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
             CALL AttrVect_importRAttr(ocn2cice_AV,TRIM(code),A,Asize)
             Iexport=Iexport+1

          CASE ('RAIN')

             CALL ROMS_export2d (ng, tile,                              &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%rain,                        &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
#ifdef REPORT_COUPLING_ALL
             write(stdout,*)                                            &
     &           'ROMS rank ',MyRank,' sending rain field (max/min):',  &
     &            maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
             CALL AttrVect_importRAttr(ocn2cice_AV,TRIM(code),A,Asize)
             Iexport=Iexport+1

          CASE ('SNOW')

             CALL ROMS_export2d (ng, tile,                              &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%snow,                        &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
#ifdef REPORT_COUPLING_ALL
             write(stdout,*)                                            &
     &           'ROMS rank ',MyRank,' sending snow field (max/min):',  &
     &            maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
             CALL AttrVect_importRAttr(ocn2cice_AV,TRIM(code),A,Asize)
             Iexport=Iexport+1

          CASE ('SW_DN')

             CALL ROMS_export2d (ng, tile,                              &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%SW_down,                     &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
#ifdef REPORT_COUPLING_ALL
             write(stdout,*)                                            &
     &           'ROMS rank ',MyRank,' sending SW_DN field (max/min):', &
     &            maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
             CALL AttrVect_importRAttr(ocn2cice_AV,TRIM(code),A,Asize)
             Iexport=Iexport+1

          CASE ('LW_DN')

             CALL ROMS_export2d (ng, tile,                              &
     &                          id, gtype, scale, add_offset,           &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          FORCES(ng)%LW_down,                     &
     &                          Fields(id)%ExpMin, Fields(id)%ExpMax,   &
     &                          Asize, A,                               &
     &                          status)
#ifdef REPORT_COUPLING_ALL
             write(stdout,*)                                            &
     &           'ROMS rank ',MyRank,' sending LW_DN field (max/min):', &
     &            maxval(A), ' ', minval(A),' gtype : ', gtype
#endif
             CALL AttrVect_importRAttr(ocn2cice_AV,TRIM(code),A,Asize)
             Iexport=Iexport+1


        END SELECT
      END DO

!
!  Send ocean fields
!
      IF (Iexport.gt.0) THEN

#ifdef REPORT_COUPLING_ALL
         if (master)                                                    &
     &        write(stdout,*)'ROMS sends', iexport,' fields to CICE'
#endif
         CALL MCT_Send(ocn2cice_AV, ROMStoCICE)
      END IF

!
!  Deallocate communication arrays.
!
      deallocate (A,uw,vw)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 36)
#endif

      RETURN
      END SUBROUTINE ocn2cice_coupling

      SUBROUTINE finalize_ocn2cice_coupling
!
!========================================================================
!                                                                       !
!  This routine finalizes ocean and ice models coupling data streams.  !
!                                                                       !
!========================================================================
!
!  Local variable declarations.
!
      integer :: MyError

!
!-----------------------------------------------------------------------
!  Deallocate MCT environment.
!-----------------------------------------------------------------------
!
      CALL Router_clean (ROMStoCICE, MyError)
      CALL AttrVect_clean (ocn2cice_AV, MyError)
      CALL GlobalSegMap_clean (GSMapROMS, MyError)

      RETURN

      END SUBROUTINE finalize_ocn2cice_coupling

# endif

#endif
      END MODULE ocean_coupler_mod
