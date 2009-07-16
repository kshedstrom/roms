#include "cppdefs.h"
      MODULE mod_iounits
!
!svn $Id: mod_iounits.F 1020 2009-07-10 23:10:30Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  ADJbase     Output adjoint base file name.                          !
!  ADJname     Output adjoint file name.                               !
!  ADSname     Input adjoint sensitivity functional file.              !
!  AVGbase     Output averages base file name.                         !
!  AVGname     Output averages file name.                              !
!  AVG2base    Output secondary base file name.                        !
!  AVG2name    Output secondary file name.                             !
!  BRYname     Input boundary data file name.                          !
!  CLMname     Input climatology file name.                            !
!  DIAbase     Output diagnostics base file name.                      !
!  DIAname     Output diagnostics file name.                           !
!  ERRname     Output 4DVar posterior error file name.                 !
!  FILname     Output filter file name.                                !
!  FLTname     Output floats file name.                                !
!  FISHname    Output fish file name.                                  !
!  FRCids      NetCDF ID descriptors associated with FRCname.          ! 
!  FRCname     Input forcing fields file name(s).                      !
!  FWDbase     Input forward solution base file name.                  !
!  FWDname     Input forward solution file name.                       !
!  GRDname     Input grid file name.                                   !
!  GSTname     Input/output GST analysis check pointing NetCDF file.   !
!  HISbase     Output history base file name.                          !
!  HISname     Output history file name.                               !
!  HSSbase     Input/output Hessian eigenvectors base file name.       !
!  HSSname     Input/output Hessian eigenvectors file name.            !
!  IADname     Input adjoint initial conditions file name.             !
!  INIname     Input nonlinear initial conditions file name.           !
!  IPRname     Input representer initial conditions file name.         !
!  ITLname     Input tangent linear initial conditions file name.      !
!  Iname       Physical parameters standard input script file name.    !
!  LCZname     Input/output Lanczos vectors file name.                 !
!  MODname     Output 4DVAR processed fields file name.                !
!  MyAppCPP    Application C-preprocessing flag.                       !
!  NRMname     Input/output error covariance normalization file name:  !
!                NRMname(1)  initial conditions                        !
!                NRMname(2)  model                                     !
!                NRMname(3)  boundary conditions                       !
!                NRMname(4)  surface forcing                           !
!  OBSname     Input/output datum observations file name.              !
!  REPname     Input/output representer coefficients file name.        !
!  Rerror      Running error messages.                                 !
!  RSTname     Output restart file name.                               !
!  SSHname     Input SSH observations file name.                       !
!  SSTname     Input SST observations file name.                       !
!  TIDEname    Input tide forcing file name.                           !
!  TLFname     Input/output tangent linear impulse forcing file name.  !
!  TLMbase     Output tangent linear base file name.                   !
!  TLMname     Output tangent linear file name.                        !
!  TOBSname    Input tracer observations file name.                    !
!  USRname     USER input/output generic file name.                    !
!  VSURname    Input surface currents observations file name.          !
!  VOBSname    Input horizontal currents observations file name.       !
!  Wname       Wave model stadard input file name.                     !
!  aparnam     Input assimilation parameters file name.                !
!  bparnam     Input biology parameters file name.                     !
!  fposnam     Input initial floats positions file name.               !
!  iparnam     Input ice parameters file name.                         !
!  ioerror     IO error flag.                                          !
!  ncfile      Current NetCDF file name being processed.               !
!  nFfiles     Number of forcing files.                                !
!  sparnam     Input sediment transport parameters file name.          !
!  sposnam     Input station positions file name.                      !
!  SourceFile  Current executed file name. It is used for IO error     !
!                purposes.                                             !
!  STAname     Output station data file name.                          !
!  STDname     Input error covariance standard deviations file name:   !
!                STDname(1)  initial conditions                        !
!                STDname(2)  model                                     !
!                STDname(3)  boundary conditions                       !
!                STDname(4)  surface forcing                           !
!  stdinp      Unit number for standard input (often 5).               !
!  stdout      Unit number for standard output (often 6).              !
!  usrout      Unit number for generic USER output.                    !
!  varname     Input IO variables information file name.               !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none

        integer, parameter :: stdinp = 5
        integer, parameter :: stdout = 6
        integer, parameter :: usrout = 10
        integer :: ioerror = 0

        integer, dimension(Ngrids) :: nFfiles

        integer, allocatable :: FRCids(:,:)

        integer, parameter :: fname_len=128

        character (len=80) :: SourceFile

        character (len=50), dimension(8) :: Rerror =                    &
     &       (/ ' ROMS/TOMS - Blows up ................ exit_flag: ',   &
     &          ' ROMS/TOMS - Input error ............. exit_flag: ',   &
     &          ' ROMS/TOMS - Output error ............ exit_flag: ',   &
     &          ' ROMS/TOMS - I/O error ............... exit_flag: ',   &
     &          ' ROMS/TOMS - Configuration error ..... exit_flag: ',   &
     &          ' ROMS/TOMS - Partition error ......... exit_flag: ',   &
     &          ' ROMS/TOMS - Illegal input parameter . exit_flag: ',   &
     &          ' ROMS/TOMS - Fatal algorithm result .. exit_flag: ' /)

        character (len=fname_len), allocatable :: FRCname(:,:)

        character (len=fname_len), dimension(Ngrids) :: ADJbase
        character (len=fname_len), dimension(Ngrids) :: ADJname
        character (len=fname_len), dimension(Ngrids) :: ADSname
        character (len=fname_len), dimension(Ngrids) :: AVGbase
        character (len=fname_len), dimension(Ngrids) :: AVGname
        character (len=fname_len), dimension(Ngrids) :: AVG2base
        character (len=fname_len), dimension(Ngrids) :: AVG2name
        character (len=fname_len), dimension(Ngrids) :: BRYname
        character (len=fname_len), dimension(Ngrids) :: CLMname
        character (len=fname_len), dimension(Ngrids) :: DIAbase
        character (len=fname_len), dimension(Ngrids) :: DIAname
        character (len=fname_len), dimension(Ngrids) :: ERRname
        character (len=fname_len), dimension(Ngrids) :: FLTname
        character (len=fname_len), dimension(Ngrids) :: FISHname
#ifdef FILTERED_RST
        character (len=fname_len), dimension(nfile,Ngrids) :: FILname
#endif
        character (len=fname_len), dimension(Ngrids) :: FWDbase
        character (len=fname_len), dimension(Ngrids) :: FWDname
        character (len=fname_len), dimension(Ngrids) :: GRDname
        character (len=fname_len), dimension(Ngrids) :: GSTname
        character (len=fname_len), dimension(Ngrids) :: HISbase
        character (len=fname_len), dimension(Ngrids) :: HISname
        character (len=fname_len), dimension(Ngrids) :: HSSbase
        character (len=fname_len), dimension(Ngrids) :: HSSname
        character (len=fname_len), dimension(Ngrids) :: IADname
        character (len=fname_len), dimension(Ngrids) :: INIname
        character (len=fname_len), dimension(Ngrids) :: IRPname
        character (len=fname_len), dimension(Ngrids) :: ITLname
        character (len=fname_len), dimension(Ngrids) :: LCZname
        character (len=fname_len), dimension(Ngrids) :: MODname
        character (len=fname_len), dimension(Ngrids) :: OBSname
        character (len=fname_len), dimension(Ngrids) :: REPname
        character (len=fname_len), dimension(Ngrids) :: SSHname
        character (len=fname_len), dimension(Ngrids) :: SSTname
        character (len=fname_len), dimension(Ngrids) :: TIDEname
        character (len=fname_len), dimension(Ngrids) :: TLFname
        character (len=fname_len), dimension(Ngrids) :: TLMbase
        character (len=fname_len), dimension(Ngrids) :: TLMname
        character (len=fname_len), dimension(Ngrids) :: TOBSname
        character (len=fname_len), dimension(Ngrids) :: VSURname
        character (len=fname_len), dimension(Ngrids) :: VOBSname
        character (len=fname_len), dimension(Ngrids) :: RSTname
        character (len=fname_len), dimension(Ngrids) :: STAname

        character (len=fname_len), dimension(4,Ngrids) :: NRMname
        character (len=fname_len), dimension(4,Ngrids) :: STDname

        character (len=fname_len) :: Iname
        character (len=fname_len) :: Wname
        character (len=fname_len) :: MyAppCPP
        character (len=fname_len) :: USRname
        character (len=fname_len) :: aparnam
        character (len=fname_len) :: bparnam
        character (len=fname_len) :: fposnam
        character (len=fname_len) :: iparnam
        character (len=fname_len) :: ncfile
        character (len=fname_len) :: sparnam
        character (len=fname_len) :: sposnam
        character (len=fname_len) :: varname

      END MODULE mod_iounits
