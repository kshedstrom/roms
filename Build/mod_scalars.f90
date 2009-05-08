      MODULE mod_scalars
!
!svn $Id: mod_scalars.F 975 2009-05-05 22:51:13Z kate $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
        USE mod_param
!
        implicit none
!
!-----------------------------------------------------------------------
! Multiple grid structure.
!-----------------------------------------------------------------------
!
!    Lstate        Logical switches to control computations of the
!                    model state.
!    SOstate       Logical switches to control computations of the
!                    stochastic optimals.
!    Cs_r          Set of S-curves used to stretch the vertical grid
!                    that follows the bathymetry at vertical RHO-points.
!    Cs_w          Set of S-curves used to stretch the vertical grid
!                    that follows the bathymetry at vertical W-points.
!    sc_r          S-coordinate independent variable, [-1 < sc < 0] at
!                    vertical RHO-points.
!    sc_w          S-coordinate independent variable, [-1 < sc < 0] at
!                    vertical W-points.
!
        TYPE T_SCALARS
          logical, pointer :: Lstate(:)
          logical, pointer :: SOstate(:)
          real(r8), pointer :: Cs_r(:)
          real(r8), pointer :: Cs_w(:)
          real(r8), pointer :: sc_r(:)
          real(r8), pointer :: sc_w(:)
        END TYPE T_SCALARS
        TYPE (T_SCALARS), allocatable :: SCALARS(:)
!
!-----------------------------------------------------------------------
!  Tracer identification indices.
!-----------------------------------------------------------------------
!
        integer :: itemp              ! Potential temperature
        integer :: isalt              ! Salinity
        integer, pointer :: idbio(:)  ! Biological tracers
        integer :: iLphy              ! Large Phytoplankton biomass
        integer :: iSphy              ! Small Phytoplankton biomass
        integer :: iLzoo              ! Large Zooplankton biomass
        integer :: iSzoo              ! Small Zooplankton biomass
        integer :: iPzoo              ! Predator Zooplankton biomass
        integer :: iNO3_              ! Nitrate concentration
        integer :: iNH4_              ! Ammonium concentration
        integer :: iPON_              ! Particulate Organic Nitrogen
        integer :: iDON_              ! Dissolved Organic Nitrogen
        integer :: iSiOH              ! Silicate concentration
        integer :: iopal              ! Particulate organic silica
!
!-----------------------------------------------------------------------
!  Diagnostic fields identification indices.
!-----------------------------------------------------------------------
!
        integer :: iTrate            ! Tracer, time rate of change
        integer :: iThadv            ! Tracer, horizontal advection
        integer :: iThdif            ! Tracer, horizontal diffusion
        integer :: iTvadv            ! Tracer, vertical advection
        integer :: iTvdif            ! Tracer, vertical diffusion
        integer :: M2pgrd            ! 2D momentum, pressure gradient
        integer :: M2fcor            ! 2D momentum, Coriolis
        integer :: M2rate            ! 2D momentum, time rate of change
        integer :: M2hadv            ! 2D momentum, horizontal advection
        integer :: M2hvis            ! 2D momentum, horizontal viscosity
        integer :: M2sstr            ! 2D momentum, surface stress
        integer :: M2bstr            ! 2D momentum, bottom stress
        integer :: M3fcor            ! 3D momentum, Coriolis
        integer :: M3rate            ! 3D momentum, time rate of change
        integer :: M3hadv            ! 3D momentum, horizontal advection
        integer :: M3hvis            ! 3D momentum, horizontal viscosity
        integer :: M3pgrd            ! 3D momentum, pressure gradient
        integer :: M3vadv            ! 3D momentum, vertical advection
        integer :: M3vvis            ! 3D momentum, vertical viscosity
!
!-----------------------------------------------------------------------
!  Floats fields identification indices.
!-----------------------------------------------------------------------
!
        integer, parameter :: itstr = 0          ! release time
        integer, parameter :: ixgrd = 1          ! x-grid location
        integer, parameter :: iygrd = 2          ! y-grid location
        integer, parameter :: izgrd = 3          ! z-grid location
        integer, parameter :: iflon = 4          ! longitude location
        integer, parameter :: iflat = 5          ! latitude location
        integer, parameter :: idpth = 6          ! depth
        integer, parameter :: ixrhs = 7          ! x-slope
        integer, parameter :: iyrhs = 8          ! y-slope
        integer, parameter :: izrhs = 9          ! z-slope
        integer, parameter :: ifden = 10         ! density anomaly
        integer, parameter :: ifwwt = 1          ! fish wet weight
        integer, parameter :: ifworth = 2        ! fish worth
! Species IDs
	integer, parameter :: if_none = 0
	integer, parameter :: if_anchovy = 1
	integer, parameter :: if_sardine = 2
	integer, parameter :: if_herring = 3
	integer, parameter :: if_polluck = 4
	integer, parameter :: if_pinksalmon = 5
!
!  Set float tracjectory types:
!
!     flt_Lagran:  3D Lagrangian floats
!     flt_Isobar:  Isobaric floats, p=g*(z+zeta)=constant
!     flt_Geopot:  Geopotential floats, constant depth
!
      integer, parameter :: flt_Lagran = 1
      integer, parameter :: flt_Isobar = 2
      integer, parameter :: flt_Geopot = 3
!
!-----------------------------------------------------------------------
!  Time stepping indices, variables, and clocks.
!-----------------------------------------------------------------------
!
!    indx1         2D timestep rolling counter.
!    iic           Timestep counter for 3D primitive equations.
!    iif           Timestep counter for 2D primitive equations.
!    ndtfast       Number of barotropic timesteps between each
!                    baroclinic timestep.
!    nfast         Number of barotropic timesteps needed to compute
!                    time-averaged barotropic variables centered at
!                    time level n+1.
!    dt            Size baroclinic timestep (s).
!    dtfast        Size barotropic timestep (s).
!    tdays         Model time clock (days).
!    time          Model time clock (s).
!    time_code     Model time clock (string, Day HH:MM:SS)
!    AVGtime       Model time clock for averages output (s).
!    AVG2time      Model time clock for averages output (s).
!    DIAtime       Model time clock for diagnostics output (s).
!    IMPtime       Impulse forcing time (s) to process.
!    ObsTime       Observation time (s) to process.
!    FrcTime       Adjoint or tangent linear Impulse forcing time (s).
!    r_date        Model initialization reference date (vector):
!                    r_date(1) => reference date (YYYYMMDD.dd).
!                    r_date(2) => year.
!                    r_date(3) => year day.
!                    r_date(4) => month.
!                    r_date(5) => day.
!                    r_date(6) => hour.
!                    r_date(7) => minute.
!                    r_date(8) => second.
!    dstart        Time stamp assigned to model initialization (usually
!                    a Calendar day, like modified Julian Day).
!    tide_start    Reference time for tidal forcing (days).
!    time_ref      Reference time for "units" attribute (YYYYMMDD.dd).
!    r_text        Model initialization reference text (string).
!
        logical, dimension(Ngrids) :: PerfectRST
        logical, dimension(Ngrids) :: PREDICTOR_2D_STEP
        integer, dimension(Ngrids) :: indx1
        integer, dimension(Ngrids) :: iic
        integer, dimension(Ngrids) :: iif
        integer, dimension(Ngrids) :: ndtfast
        integer, dimension(Ngrids) :: nfast
        real(r8), dimension(Ngrids) :: dt                ! seconds
        real(r8), dimension(Ngrids) :: dtfast            ! seconds
        real(r8), dimension(Ngrids) :: tdays             ! days
        real(r8), dimension(Ngrids) :: time              ! seconds
        real(r8), dimension(Ngrids) :: AVGtime           ! seconds
        real(r8), dimension(Ngrids) :: AVG2time           ! seconds
        real(r8), dimension(Ngrids) :: DIAtime           ! seconds
        real(r8), dimension(Ngrids) :: IMPtime           ! seconds
        real(r8), dimension(Ngrids) :: ObsTime           ! seconds
        real(r8), dimension(Ngrids) :: FrcTime           ! seconds
        integer, dimension(Ngrids) :: idays
        integer, dimension(Ngrids) :: iyesterday
        real(r8), dimension(8) :: r_date
        real(r8) :: dstart = 0.0_r8                      ! days
        real(r8) :: tide_start = 0.0_r8                  ! days
        real(r8) :: time_ref = 0.0_r8                    ! YYYYMMDD.dd
        character (len=14) :: time_code(Ngrids)          ! DD HH:MM:SS
        character (len=19) :: r_text
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  Fields.  The power-law shape filters are given by:
!
!     F(xi)=xi^Falpha*(1-xi^Fbeta)-Fgamma*xi 
!
!  Possible settings of parameters to yield the second-order accuracy:
!
!     Falpha  Fbeta      Fgamma
!     ------------------------------
!      2.0     1.0    0.1181  0.169     The problem here is setting
!      2.0     2.0    0.1576  0.234     Fgamma. Its value here is
!      2.0     3.0    0.1772  0.266     understood as the MAXIMUM
!      2.0     4.0    0.1892  0.284     allowed. It is computed using
!      2.0     5.0    0.1976  0.296     a Newton iteration scheme.
!      2.0     6.0    0.2039  0.304
!      2.0     8.0    0.2129  0.314
!
!  NOTE: Theoretical values of Fgamma presented in the table above are
!  derived assuming "exact" barotropic mode stepping. Consequently, it
!  does not account for effects caused by Forward-Euler (FE) startup
!  of the barotropic mode at every 3D time step.  As the result, the
!  code may become unstable if the theoretical value of Fgamma is used
!  when mode splitting ratio "ndtfast" is small, thus yielding non-
!  negligible start up effects.  To compensate this, the accepted
!  value of Fgamma is reduced relatively to theoretical one, depending
!  on splitting ratio "ndtfast".  This measure is empirical. It is
!  shown to work with setting of "ndtfast" as low as 15, which is
!  more robust that the Hamming Window the squared cosine weights
!  options in "set_weights".
!
        real(r8) :: Falpha = 2.0_r8
        real(r8) :: Fbeta  = 4.0_r8
        real(r8) :: Fgamma = 0.284_r8
!
!  Total number timesteps in current run. In 3D configurations, "ntimes"
!  is the total of baroclinic timesteps. In 2D configuration, "ntimes"
!  is the total of barotropic timesteps.
!
        integer, dimension(Ngrids) :: ntimes
!
!  Number of time interval divisions for stochastic optimals.  It must
!  a multiple of "ntimes".
!
        integer :: Nintervals = 1
!
!  Starting, current, and ending ensemble run parameters.
!
        integer :: ERstr = 1                    ! Starting value
        integer :: ERend = 1                    ! Ending value
        integer :: Ninner = 1                   ! number of inner loops
        integer :: Nouter = 1                   ! number of outer loops
        integer :: Nrun = 1                     ! Current counter
        integer :: inner = 0                    ! inner loop counter
        integer :: outer = 0                    ! outer loop counter
!
!  First, starting, and ending timestepping parameters
!
        integer, dimension(Ngrids) :: ntfirst   ! Forward-Euler step
        integer, dimension(Ngrids) :: ntstart   ! Start step
        integer, dimension(Ngrids) :: ntend     ! End step
!
!  Adjoint model or tangent linear model impulse forcing time record
!  counter and number of records available.
!
        integer, dimension(Ngrids) :: FrcRec
        integer, dimension(Ngrids) :: NrecFrc
!
!-----------------------------------------------------------------------
!  Control switches.
!-----------------------------------------------------------------------
!
!  These switches are designed to control computational options within
!  nested and/or multiple connected grids.  They are .TRUE. by default.
!  They can turned off for a particular grind in input scripts.
!
        logical, dimension(Ngrids) :: Lassimilate
        logical, dimension(Ngrids) :: Lbiology
        logical, dimension(Ngrids) :: Lfloats
        logical, dimension(Ngrids) :: Lsediment
        logical, dimension(Ngrids) :: Lstations
!
!-----------------------------------------------------------------------
!  Physical constants.   
!-----------------------------------------------------------------------
!
!    Cp            Specific heat for seawater (Joules/Kg/degC).
!    Csolar        Solar irradiantion constant (W/m2).
!    Eradius       Earth equatorial radius (m).
!    StefBo        Stefan-Boltzmann constant (W/m2/K4).
!    emmiss        Infrared emmissivity.
!    g             Acceleration due to gravity (m/s2).
!    gorho0        gravity divided by mean density anomaly.
!    rhow          fresh water density (kg/m3).
!    vonKar        von Karman constant.
!
        real(r8) :: Cp = 3985.0_r8              ! Joules/kg/degC
        real(r8) :: Csolar = 1353.0_r8          ! 1360-1380 W/m2
        real(r8) :: Eradius = 6371315.0_r8      ! m
        real(r8) :: StefBo = 5.67E-8_r8         ! Watts/m2/K4
        real(r8) :: emmiss = 0.97_r8            ! non_dimensional
        real(r8) :: rhow = 1000.0_r8            ! kg/m3         
        real(r8) :: g = 9.81_r8                 ! m/s2
        real(r8) :: gorho0                      ! m4/s2/kg
        real(r8) :: vonKar = 0.41_r8            ! non-dimensional
!
!-----------------------------------------------------------------------
!  Various model parameters.  Some of these parameters are overwritten
!  with the values provided from model standard input script.
!-----------------------------------------------------------------------
!
!  Switch for spherical grid (lon,lat) configurations.
!
        logical :: spherical = .FALSE.
!
!  Switch to compute the grid stiffness.
!
        logical :: Lstiffness = .TRUE.
!
!  Execution termination flag.
!
!    exit_flag = 0   No error
!    exit_flag = 1   Blows up
!    exit_flag = 2   Input error
!    exit_flag = 3   Output error
!    exit_flag = 4   IO error
!    exit_flag = 5   Configuration error
!    exit_flag = 6   Partition error
!    exit_flag = 7   Illegal input parameter
!    exit_flag = 8   Fatal algorithm result
!     
        integer :: exit_flag = 0
        integer :: blowup = 0
        integer :: NoError = 0
!
!  Set threshold maximum speed (m/s) and density anomaly (kg/m3) to
!  test if the model is blowing-up.
!
        real(r8), dimension(Ngrids) :: maxspeed
        real(r8), dimension(Ngrids) :: maxrho
!
        real(r8) :: max_speed = 20.0_r8         ! m/s
        real(r8) :: max_rho = 200.0_r8          ! kg/m3
        real(r8), dimension(NBT, Ngrids) :: maxbio
        real(r8), dimension(NBT) :: max_bio = 10000.0_r8   ! bio units
!
!  Interpolation scheme.
!
        integer, parameter :: linear = 0        ! linear interpolation
        integer, parameter :: cubic  = 1        ! cubic  interpolation
!
        integer :: InterpFlag = cubic          ! interpolation flag
!
!  Shallowest and Deepest levels to apply bottom momemtum stresses as
!  a bodyforce
!
        integer, dimension(Ngrids) :: levsfrc
        integer, dimension(Ngrids) :: levbfrc
!
!  Vertical coordinates transform.  Currently, there are two vertical
!  transformation equations (see set_scoord.F for details):
!
!    Original transform (Vtransform=1):
!
!         z_r(x,y,s,t) = Zo_r + zeta(x,y,t) * [1.0 + Zo_r / h(x,y)]
!
!                 Zo_r = hc * [s(k) - C(k)] + C(k) * h(x,y)
!
!    New transform (Vtransform=2):
!
!         z_r(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_r
!
!                 Zo_r = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)] 
!
        integer :: Vtransform(Ngrids)
!
!  Vertical grid stretching function flag:
!
!    Vstretcing = 1   Original function (Song and Haidvogel, 1994)
!               = 2   A. Shchepetkin (ROMS-UCLA) function
!               = 3   R. Geyer BBL function
!
        integer :: Vstretching(Ngrids)
!
!  Vertical grid stretching parameters.
!
!    Tcline        Width (m) of surface or bottom boundary layer in
!                    which higher vertical resolution is required
!                    during stretching.
!    hc            S-coordinate critical depth, hc=MIN(hmin,Tcline).
!    theta_s       S-coordinate surface control parameter.
!    theta_b       S-coordinate bottom control parameter.
!
        real(r8), dimension(Ngrids) :: Tcline   ! m, positive
        real(r8), dimension(Ngrids) :: hc       ! m, positive
        real(r8), dimension(Ngrids) :: theta_s  ! 0 < theta_s < 20
        real(r8), dimension(Ngrids) :: theta_b  ! 0 < theta_b < 1
!
!  Bathymetry range values.
!
        real(r8), dimension(Ngrids) :: hmin     ! m, positive
        real(r8), dimension(Ngrids) :: hmax     ! m, positive
!
!  Length (m) of domain box in the XI- and ETA-directions.
!
        real(r8), dimension(Ngrids) :: xl       ! m
        real(r8), dimension(Ngrids) :: el       ! m
!
!  Minimum and Maximum longitude and latitude at RHO-points
!
        real(r8), dimension(Ngrids) :: LonMin   ! degrees east
        real(r8), dimension(Ngrids) :: LonMax   ! degrees east
        real(r8), dimension(Ngrids) :: LatMin   ! degrees north
        real(r8), dimension(Ngrids) :: LatMax   ! degrees north
!
!  Diagnostic volume averaged variables.
!
        integer, dimension(Ngrids) :: first_time
        real(r8) :: avgke = 0.0_r8              ! Kinetic energy
        real(r8) :: avgpe = 0.0_r8              ! Potential energy
        real(r8) :: avgkp = 0.0_r8              ! Total energy
        real(r8) :: volume = 0.0_r8             ! diagnostics volume
        real(r8) :: ad_volume = 0.0_r8          ! adjoint volume
        real(r8) :: TotVolume = 0.0_r8          ! Total volume
        real(r8) :: MaxVolume = 0.0_r8          ! Minimum cell volume
        real(r8) :: MinVolume = 0.0_r8          ! Maximum cell volume
!
!  Minimun and maximum grid spacing
!
        real(r8), dimension(Ngrids) :: DXmin
        real(r8), dimension(Ngrids) :: DXmax
        real(r8), dimension(Ngrids) :: DYmin
        real(r8), dimension(Ngrids) :: DYmax
        real(r8), dimension(Ngrids) :: DZmin
        real(r8), dimension(Ngrids) :: DZmax
!
!  Maximum size of a grid node (m) over the whole curvilinear grid
!  application. Used for scaling horizontal mixing by the grid size.
!
        real(r8), dimension(Ngrids) :: grdmax
!
!  Courant Numbers.
!
        real(r8) :: Cu_min = 0.0_r8             ! Minimun barotropic
        real(r8) :: Cu_max = 0.0_r8             ! Maximun barotropic
        real(r8) :: Cu_Cor = 0.0_r8             ! Maximun Coriolis
!
!  Linear equation of state parameters.
!
!    R0            Background constant density anomaly (kg/m3).
!    Tcoef         Thermal expansion coefficient (1/Celsius).
!    Scoef         Saline contraction coefficient (1/PSU).
!
        real(r8), dimension(Ngrids) :: R0
        real(r8), dimension(Ngrids) :: Tcoef
        real(r8), dimension(Ngrids) :: Scoef
!
!  Background potential temperature (Celsius) and salinity (PSU) values
!  used in analytical initializations.
!                              
        real(r8), dimension(Ngrids) :: T0
        real(r8), dimension(Ngrids) :: S0
!
!  Slipperiness variable, either 1.0 (free slip) or -1.0 (no slip).
!
        real(r8), dimension(Ngrids) :: gamma2
!
!  Weighting coefficient for the newest (implicit) time step derivatives
!  in the Crack-Nicolson implicit scheme (usually, lambda=0.5).
!
        real(r8) :: lambda = 1.0_r8             ! backward implicit
!
!  Jerlov water type to assign everywhere, range values: 1 - 5.
!
        integer, dimension(Ngrids) :: lmd_Jwt
!
!  Grid r-factor (non-dimensional).
!
        real(r8) :: rx0 = 0.0_r8                ! Beckmann and Haidvogel
        real(r8) :: rx1 = 0.0_r8                ! Haney
!
!  Linear (m/s) and quadratic (nondimensional) bottom drag coefficients.
!
        real(r8), dimension(Ngrids) :: rdrg
        real(r8), dimension(Ngrids) :: rdrg2
!
!  Minimum and maximum threshold for transfer coefficient of momentum.
!
        real(r8) :: Cdb_min = 0.000001_r8
        real(r8) :: Cdb_max = 0.5_r8
!
!  Surface and bottom roughness (m)
!
        real(r8), dimension(Ngrids) :: Zos
        real(r8), dimension(Ngrids) :: Zob
!
!  Minimum depth for wetting and drying (m).
!
        real(r8), dimension(Ngrids) :: Dcrit
!
!  Mean density (Kg/m3) used when the Boussinesq approximation is
!  inferred.
!
        real(r8) :: rho0 = 1025.0_r8
!
!  Background Brunt-Vaisala frequency (1/s2)
!
        real(r8) :: bvf_bak = 0.00001_r8
!
!  Vector containing USER generic parameters.
!
        integer :: Nuser
        real(r8), dimension(25) :: user(25)
!
!  Weights for the time average of 2D fields.
!
        real(r8), dimension(2,0:256,Ngrids) :: weight
!
!  Constants.
!
        real(r8), parameter :: pi = 3.14159265358979323846_r8
        real(r8), parameter :: deg2rad = pi / 180.0_r8
        real(r8), parameter :: rad2deg = 180.0_r8 / pi
        real(r8), parameter :: day2sec = 86400.0_r8
        real(r8), parameter :: sec2day = 1.0_r8 / 86400.0_r8
        real(r8), parameter :: spval = 1.0E+37_r8
        real(r8), parameter :: jul_off = 2440000.0_r8
!
!  Set special check value.  Notice that a smaller value is assigned
!  to account for both NetCDF fill value and roundoff. There are
!  many Matlab scripts out there that do not inquire correctly
!  the spval from the _FillValue attribute in single/double
!  precision.
!
        real(r8), parameter :: spval_check = 1.0E+35_r8
!
!-----------------------------------------------------------------------
!  Horizontal and vertical constant mixing coefficients.
!-----------------------------------------------------------------------
!
!    Akk_bak       Background vertical mixing coefficient (m2/s) for
!                    turbulent energy.
!    Akp_bak       Background vertical mixing coefficient (m2/s) for
!                    generic statistical field "psi".
!    Akt_bak       Background vertical mixing coefficient (m2/s) for
!                    tracers.
!    Akv_bak       Background vertical mixing coefficient (m2/s) for
!                    momentum.
!    Kdiff         Isopycnal mixing thickness diffusivity (m2/s) for
!                    tracers.
!    visc2         Lateral harmonic constant mixing coefficient
!                    (m2/s) for momentum.
!    visc4         Square root lateral biharmonic constant mixing
!                    coefficient (m2 s^-1/2) for momentum.
!    tkenu2        Lateral harmonic constant mixing coefficient
!                    (m2/s) for turbulent energy.
!    tkenu4        Square root lateral biharmonic constant mixing
!                    coefficient (m2 s^-1/2) for turbulent energy.
!    tnu2          Lateral harmonic constant mixing coefficient
!                    (m2/s) for tracer type variables.
!    tnu4          Square root lateral biharmonic constant mixing
!                    coefficient (m2 s^-1/2) for tracer type variables.
!
        real(r8), dimension(Ngrids) :: Akk_bak       ! m2/s
        real(r8), dimension(Ngrids) :: Akp_bak       ! m2/s
        real(r8), dimension(Ngrids) :: Akv_bak       ! m2/s
        real(r8), dimension(Ngrids) :: visc2         ! m2/s
        real(r8), dimension(Ngrids) :: visc4         ! m2 s-1/2
        real(r8), dimension(Ngrids) :: tkenu2        ! m2/s
        real(r8), dimension(Ngrids) :: tkenu4        ! m2 s-1/2
        real(r8), allocatable :: Akt_bak(:,:)        ! m2/s
        real(r8), allocatable :: Kdiff(:,:)          ! m2/s
        real(r8), allocatable :: tnu2(:,:)           ! m2/s
        real(r8), allocatable :: tnu4(:,:)           ! m2 s-1/2
!
!  Horizontal diffusive relaxation coefficients (m2/s) used to smooth
!  representer tangent linear solution during Picard iterations to
!  improve stability and convergence.
!
        real(r8), dimension(Ngrids) :: tl_M2diff     ! 2D momentum
        real(r8), dimension(Ngrids) :: tl_M3diff     ! 3D momentum
        real(r8), allocatable :: tl_Tdiff(:,:)       ! tracers
!
!-----------------------------------------------------------------------
!  IO parameters.
!-----------------------------------------------------------------------
!
!  Switches to activate creation and writing of output NetCDF files.
!
        logical, dimension(Ngrids) :: LdefADJ    ! Adjoint file
        logical, dimension(Ngrids) :: LdefAVG    ! Average file
        logical, dimension(Ngrids) :: LdefAVG2   ! Average file
        logical, dimension(Ngrids) :: LdefDIA    ! Diagnostics file
        logical, dimension(Ngrids) :: LdefFLT    ! Floats file
        logical, dimension(Ngrids) :: LdefFISH   ! Fish file
        logical, dimension(Ngrids) :: LdefHIS    ! History file
        logical, dimension(Ngrids) :: LdefHSS    ! Hessian file
        logical, dimension(Ngrids) :: LdefINI    ! Initial file
        logical, dimension(Ngrids) :: LdefIRP    ! Initial RPM file
        logical, dimension(Ngrids) :: LdefITL    ! Initial TLM file
        logical, dimension(Ngrids) :: LdefLCZ    ! Lanczos file
        logical, dimension(Ngrids) :: LdefMOD    ! 4DVAR file
        logical, dimension(Ngrids) :: LdefRST    ! Restart file
        logical, dimension(Ngrids) :: LdefSTA    ! Stations file
        logical, dimension(Ngrids) :: LdefTIDE   ! tide forcing file
        logical, dimension(Ngrids) :: LdefTLM    ! Tangent linear file
        logical, dimension(Ngrids) :: LdefTLF    ! TLM/RPM impulse file
        logical, dimension(Ngrids) :: LwrtADJ    ! Write adjoint file
        logical, dimension(Ngrids) :: LwrtAVG    ! Write average file
        logical, dimension(Ngrids) :: LwrtAVG2   ! Write average file
        logical, dimension(Ngrids) :: LwrtDIA    ! Write diagnostic file
        logical, dimension(Ngrids) :: LwrtHIS    ! Write history file
        logical, dimension(Ngrids) :: LwrtPER    ! Write during ensemble
        logical, dimension(Ngrids) :: LwrtRST    ! Write restart file
        logical, dimension(Ngrids) :: LwrtTLM    ! Write tangent file
        logical, dimension(Ngrids) :: LwrtTLF    ! Write impulse file
        logical, dimension(4,Ngrids) :: LdefNRM  ! Norm file
        logical, dimension(4,Ngrids) :: LwrtNRM  ! Write norm file
!
!  Switch to write out adjoint 2D state arrays instead of IO solution
!  arrays and adjoint ocean time. This is used in 4DVAR for IO
!  maniputations.
!
        logical, dimension(Ngrids) :: LwrtState2d
        logical, dimension(Ngrids) :: LwrtTime
!
!  Switch to write out adjoint surface forcing fields adjusted by the
!  4DVAR algorithms.
!
        logical, dimension(Ngrids) :: Ladjusted
!
!  Switch to write application set-up information to standard output.
!
        logical, dimension(Ngrids) :: LwrtInfo
!
!  Switch used to create new output NetCDF files. If TRUE, new output
!  files are created. If FALSE, data is appended to an existing output
!  files.  Used only for history, average and station files.
!
        logical, dimension(Ngrids) :: ldefout    ! New output files
!
!  Number of timesteps between creation of new output files.
!
        integer, dimension(Ngrids) :: ndefADJ    ! Adjoint file
        integer, dimension(Ngrids) :: ndefAVG    ! Average file
        integer, dimension(Ngrids) :: ndefAVG2   ! Average file
        integer, dimension(Ngrids) :: ndefDIA    ! Diagnostics file
        integer, dimension(Ngrids) :: ndefHIS    ! History file
        integer, dimension(Ngrids) :: ndefTLM    ! Tangent linear file
!
!  Starting timestep for accumulation of output.
!
        integer, dimension(Ngrids) :: ntsAVG     ! Average file
        integer, dimension(Ngrids) :: ntsAVG2    ! Average file
        integer, dimension(Ngrids) :: ntsDIA     ! Diagnostics file
!
!  Number of timesteps between writing of output data.
!
        integer, dimension(Ngrids) :: nADJ       ! Adjoint file
        integer, dimension(Ngrids) :: nAVG       ! Average file
        integer, dimension(Ngrids) :: nAVG2      ! Average file
        integer, dimension(Ngrids) :: nDIA       ! Diagnostics file
        integer, dimension(Ngrids) :: nFLT       ! Floats file
        integer, dimension(Ngrids) :: nFISHES    ! Fish file
        integer, dimension(Ngrids) :: nHIS       ! History file
        integer, dimension(Ngrids) :: nRST       ! Restart file
        integer, dimension(Ngrids) :: nSTA       ! Stations file
        integer, dimension(Ngrids) :: nTLM       ! Tangent linear file
!
!  Number of time records written into output files.
!
        integer, dimension(Ngrids) :: NrecADJ    ! Adjoint file
        integer, dimension(Ngrids) :: NrecAVG    ! Average file
        integer, dimension(Ngrids) :: NrecAVG2   ! Average file
        integer, dimension(Ngrids) :: NrecDIA    ! Diagnostics file
        integer, dimension(Ngrids) :: NrecFIL    ! Filters file
        integer, dimension(Ngrids) :: NrecFLT    ! Floats file
        integer, dimension(Ngrids) :: NrecFISH   ! Fish file
        integer, dimension(Ngrids) :: NrecHIS    ! History file
        integer, dimension(Ngrids) :: NrecHSS    ! Hessian file
        integer, dimension(Ngrids) :: NrecINI    ! NLM Initial file
        integer, dimension(Ngrids) :: NrecITL    ! TLM Initial file
        integer, dimension(Ngrids) :: NrecLCZ    ! Lanczos file
        integer, dimension(Ngrids) :: NrecRST    ! Restart file
        integer, dimension(Ngrids) :: NrecSTA    ! Station file
        integer, dimension(Ngrids) :: NrecTLM    ! Tangent linear file
        integer, dimension(4,Ngrids) :: NrecNRM  ! norm file
!
!  Number of timesteps between print of single line information to
!  standard output.
!
        integer, dimension(Ngrids) :: ninfo
!
!  Number of timesteps between 4DVAR adjustment of open boundaries.
!  In strong constraint 4DVAR, it is possible to open bounadies at
!  other intervals in addition to initial time. These parameters are
!  used to store the appropriate number of open boundary records in
!  output history NetCDF files.
!
!    Nbrec(:) = 1 + ntimes(:) / nOBC(:)
!
!  Here, it is assumed that nOBC is a multiple of NTIMES or greater
!  than NTIMES. If nOBC > NTIMES, only one record is stored in the
!  output history NetCDF files and the adjustment is for constant
!  open boundaries with constant correction.
!
        integer, dimension(Ngrids) :: nOBC       ! number of timesteps
        integer, dimension(Ngrids) :: Nbrec      ! number of records
        integer, dimension(Ngrids) :: OBCcount   ! record counter
!
!  Number of timesteps between adjustment of 4DVAR surface forcing
!  fields. In strong constraint 4DVAR, it is possible to adjust surface
!  forcing fields at other intervals in addition to initial time.
!  These parameters are used to store the appropriate number of
!  surface forcing records in output history NetCDF files.
!
!    Nfrec(:) = 1 + ntimes(:) / nSFF(:)
!
!  Here, it is assumed that nSFF is a multiple of NTIMES or greater
!  than NTIMES. If nSFF > NTIMES, only one record is stored in the
!  output history NetCDF files and the adjustment is for constant
!  forcing with constant correction.
!
        integer, dimension(Ngrids) :: nSFF       ! number of timesteps
        integer, dimension(Ngrids) :: Nfrec      ! number of records
        integer, dimension(Ngrids) :: SFcount    ! record counter
!
!  Restart time record to read from disk and use as the initial
!  conditions. Use nrrec=0 for new solutions. If nrrec is negative
!  (say, nrrec=-1), the model will restart from the most recent
!  time record. That is, the initialization record is assigned
!  internally.
!
        integer, dimension(Ngrids) :: nrrec
!
!  Fish restart switch.
!
        integer, dimension(Ngrids) :: fishrrec
!
!  Switch to activate processing of input data.  This switch becomes
!  very useful when reading input data serially in parallel
!  applications.
!
        logical, dimension(Ngrids) :: synchro_flag
!
!  Switch to inialize model with latest time record from initial
!  (restart/history) NetCDF file.
!
        logical, dimension(Ngrids) :: LastRec
!
!  Generalized Statbility Theory (GST) parameters.
!
        logical :: LrstGST                  ! restart switch
        integer :: MaxIterGST               ! Number of iterations
        integer :: nGST                     ! check pointing interval
!
!  Switches used to recycle time records in some output file. If TRUE,
!  only the latest two time records are maintained.  If FALSE, all
!  field records are saved.
!
        logical, dimension(Ngrids) :: LcycleADJ
        logical, dimension(Ngrids) :: LcycleRST
        logical, dimension(Ngrids) :: LcycleTLM
!
!-----------------------------------------------------------------------
!  Adjoint sensitivity parameters.
!-----------------------------------------------------------------------
!
!  Starting and ending vertical levels of the 3D adjoint state whose
!  sensitivity is required.
!
        integer, dimension(Ngrids) :: KstrS        ! starting level
        integer, dimension(Ngrids) :: KendS        ! ending level
!
!  Starting and ending day for adjoint sensitivity forcing.
!
        real(r8), dimension(Ngrids) :: DstrS       ! starting day
        real(r8), dimension(Ngrids) :: DendS       ! ending day
!
!-----------------------------------------------------------------------
!  Stochastic optimals parameters.
!-----------------------------------------------------------------------
!
!  Trace of stochastic optimals matrix.
!
        real(r8), dimension(Ngrids) :: TRnorm
!
!  Stochastic optimals time decorrelation scale (days) assumed for
!  red noise processes.
!
        real(r8), dimension(Ngrids) :: SO_decay
!
!  Stochastic optimals surface forcing standard deviation for
!  dimensionalization.
!
        real(r8), allocatable :: SO_sdev(:,:)
!
!-----------------------------------------------------------------------
!  Nudging variables for passive (outflow) and active (inflow) oepn
!  boundary conditions.
!-----------------------------------------------------------------------
!
!    iwest         West  identification index in boundary arrays.
!    isouth        South identification index in boundary arrays.
!    ieast         East  identification index in boundary arrays.
!    inorth        North identification index in boundary arrays.
!    obcfac        Factor between passive and active open boundary
!                    conditions (nondimensional and greater than one).  
!                    The nudging time scales for the active conditions
!                    are obtained by multiplying the passive values by
!                    factor.
!    FSobc_in      Active and strong time-scale (1/sec) coefficients
!                    for nudging towards free-surface data at  inflow.
!    FSobc_out     Passive and weak  time-scale (1/sec) coefficients
!                    for nudging towards free-surface data at outflow.
!    M2obc_in      Active and strong time-scale (1/sec) coefficients
!                    for nudging towards 2D momentum data at  inflow.
!    M2obc_out     Passive and weak  time-scale (1/sec) coefficients
!                    for nudging towards 2D momentum data at outflow.
!    M3obc_in      Active and strong time-scale (1/sec) coefficients
!                    for nudging towards 3D momentum data at  inflow.
!    M3obc_out     Passive and weak  time-scale (1/sec) coefficients
!                    for nudging towards 3D momentum data at outflow.
!    Tobc_in       Active and strong time-scale (1/sec) coefficients
!                    for nudging towards tracer data at  inflow.
!    Tobc_out      Passive and weak  time-scale (1/sec) coefficients
!                    for nudging towards tracer data at outflow.
!
        integer, parameter :: iwest = 1
        integer, parameter :: isouth = 2
        integer, parameter :: ieast = 3
        integer, parameter :: inorth = 4
        real(r8), dimension(Ngrids) :: obcfac
        real(r8), dimension(Ngrids,4) :: FSobc_in
        real(r8), dimension(Ngrids,4) :: FSobc_out
        real(r8), dimension(Ngrids,4) :: M2obc_in
        real(r8), dimension(Ngrids,4) :: M2obc_out
        real(r8), dimension(Ngrids,4) :: M3obc_in
        real(r8), dimension(Ngrids,4) :: M3obc_out
        real(r8), allocatable :: Tobc_in(:,:,:)
        real(r8), allocatable :: Tobc_out(:,:,:)
!
!  Inverse time-scales (1/s) for nudging at open boundaries and sponge
!  areas.
!
        real(r8), dimension(Ngrids) :: Znudg       ! Free-surface
        real(r8), dimension(Ngrids) :: M2nudg      ! 2D momentum
        real(r8), dimension(Ngrids) :: M3nudg      ! 3D momentum
        real(r8), allocatable :: Tnudg (:,:)       ! Tracers
!
!  Inverse time-scales (1/s) for assimilation via nudging.
!
        real(r8), dimension(Ngrids) :: Znudass     ! Free-surface
        real(r8), dimension(Ngrids) :: M2nudass    ! 2D momentum
        real(r8), dimension(Ngrids) :: M3nudass    ! 3D momentum
        real(r8), allocatable :: Tnudass(:,:)      ! Tracers
!
!  Variables used to impose mass flux conservation in open boundary
!  configurations.
!
        real(r8) :: bc_area = 0.0_r8
        real(r8) :: bc_flux = 0.0_r8
        real(r8) :: ubar_xs = 0.0_r8
!
!-----------------------------------------------------------------------
!  Generic Length Scale parameters.
!-----------------------------------------------------------------------
!
!    gls_Gh0
!    gls_Ghcri
!    gls_Ghmin
!    gls_Kmin      Minimum value of specific turbulent kinetic energy.
!    gls_Pmin      Minimum Value of dissipation.
!    gls_cmu0      Stability coefficient (non-dimensional).
!    gls_c1        Shear production coefficient (non-dimensional).
!    gls_c2        Dissipation coefficient (non-dimensional).
!    gls_c3m       Buoyancy production coefficient (minus).
!    gls_c3p       Buoyancy production coefficient (plus).
!    gls_E2
!    gls_m         Turbulent kinetic energy exponent (non-dimensional).
!    gls_n         Turbulent length scale exponent (non-dimensional).
!    gls_p         Stability exponent (non-dimensional).
!    gls_sigk      Constant Schmidt number (non-dimensional) for
!                    turbulent kinetic energy diffusivity.
!    gls_sigp      Constant Schmidt number (non-dimensional) for
!                    turbulent generic statistical field, "psi".
!
        real(r8), dimension(Ngrids) :: gls_m
        real(r8), dimension(Ngrids) :: gls_n
        real(r8), dimension(Ngrids) :: gls_p
        real(r8), dimension(Ngrids) :: gls_sigk
        real(r8), dimension(Ngrids) :: gls_sigp
        real(r8), dimension(Ngrids) :: gls_cmu0
        real(r8), dimension(Ngrids) :: gls_cmupr
        real(r8), dimension(Ngrids) :: gls_c1
        real(r8), dimension(Ngrids) :: gls_c2
        real(r8), dimension(Ngrids) :: gls_c3m
        real(r8), dimension(Ngrids) :: gls_c3p
        real(r8), dimension(Ngrids) :: gls_Kmin
        real(r8), dimension(Ngrids) :: gls_Pmin
!
! Constants used in the various formulation of surface flux boundary
! conditions for the GLS vertical turbulence closure in terms of
! Charnok surface roughness (CHARNOK_ALPHA), roughness from wave
! amplitude (zos_hsig_alpha), wave dissipation (SZ_ALPHA), and
! Craig and Banner wave breaking (CRGBAN_CW).
!
        real(r8), dimension(Ngrids) :: charnok_alpha
        real(r8), dimension(Ngrids) :: zos_hsig_alpha
        real(r8), dimension(Ngrids) :: sz_alpha
        real(r8), dimension(Ngrids) :: crgban_cw
!
!-----------------------------------------------------------------------
!  Tangent linear and adjoint model parameters.
!-----------------------------------------------------------------------
!
!  Tangent linear and adjoint model control switches.
!
        logical :: TLmodel = .FALSE.
        logical :: ADmodel = .FALSE.
      CONTAINS
      SUBROUTINE initialize_scalars
!
!=======================================================================
!                                                                      !
!  This routine initializes several variables in module "mod_scalars"  !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: i, ic, j, ng, itrc
      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate and initialize variables in module structure.
!-----------------------------------------------------------------------
!
      allocate ( SCALARS(Ngrids) )
      DO ng=1,Ngrids
        allocate ( SCALARS(ng) % Lstate(5+MT) )
        SCALARS(ng) % Lstate(1:5+MT) = .FALSE.
        allocate ( SCALARS(ng) % SOstate(2+MT) )
        SCALARS(ng) % SOstate(1:2+MT) = .FALSE.      
        allocate ( SCALARS(ng) % Cs_r(N(ng)) )
        SCALARS(ng) % Cs_r(1:N(ng)) = IniVal
        allocate ( SCALARS(ng) % Cs_w(0:N(ng)) )
        SCALARS(ng) % Cs_w(0:N(ng)) = IniVal
        allocate ( SCALARS(ng) % sc_r(N(ng)) )
        SCALARS(ng) % sc_r(1:N(ng)) = IniVal
        allocate ( SCALARS(ng) % sc_w(0:N(ng)) )
        SCALARS(ng) % sc_w(0:N(ng)) = IniVal
      END DO                
!
!-----------------------------------------------------------------------
!  Select the bigger value of number of tracer between all grids. Then,
!  allocate tracer type variables.
!-----------------------------------------------------------------------
!
      allocate ( Akt_bak(MT,Ngrids) )
      allocate ( Kdiff(MT,Ngrids) )
      allocate ( SO_sdev(2+MT,Ngrids) )
      allocate ( tnu2(MT,Ngrids) )
      allocate ( tnu4(MT,Ngrids) )
      allocate ( tl_Tdiff(MT,Ngrids) )
      allocate ( Tobc_in(MT,Ngrids,4) )
      Tobc_in = IniVal
      allocate ( Tobc_out(MT,Ngrids,4) )
      Tobc_out = IniVal
      allocate ( idbio(NBT) )
      allocate ( Tnudg(MT,Ngrids) )
      allocate ( Tnudass(MT,Ngrids) )
!
!---------------------------------------------------------------------
!  Set tracer identification indices.
!---------------------------------------------------------------------
!
      itemp=1
      isalt=2
      ic=NAT
!
!  Set Nemuro ecosystem tracer indices.
!
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      iSphy=ic+1
      iLphy=ic+2
      iSzoo=ic+3
      iLzoo=ic+4
      iPzoo=ic+5
      iNO3_=ic+6
      iNH4_=ic+7
      iPON_=ic+8
      iDON_=ic+9
      iSiOH=ic+10
      iopal=ic+11
      ic=ic+11
!
!---------------------------------------------------------------------
!  Set diagnostic fields identification indices.
!---------------------------------------------------------------------
!
      iThadv=1
      iTvadv=2
      ic=2
      iThdif=ic+1
      ic=ic+1
      iTvdif=ic+1
      iTrate=ic+2
!
!  Indices for 2D momentum diagnostic variables.
!  Order is important!
!
      M2pgrd=1
      ic=1
      M2fcor=ic+1
      ic=ic+1
      M2hadv=ic+1
      ic=ic+1
      M2sstr=ic+1
      ic=ic+1
      M2bstr=ic+1
      M2hvis=ic+2
      M2rate=NDM2d
!
!  Indices for 3D momentum diagnostic variables.
!  Order is important!
!
      M3pgrd=1
      ic=1
      M3fcor=ic+1
      ic=ic+1
      M3hadv=ic+1
      M3vadv=ic+2
      ic=ic+2
      M3vvis=ic+1
      M3hvis=ic+2
      M3rate=NDM3d
!
!-----------------------------------------------------------------------
!  Activate all computation control switches.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        LastRec(ng)=.FALSE.
        Lassimilate(ng)=.TRUE.
        Lbiology(ng)=.TRUE.
        LcycleADJ(ng)=.FALSE.
        LcycleRST(ng)=.FALSE.
        LcycleTLM(ng)=.FALSE.
        Lfloats(ng)=.TRUE.
        Lsediment(ng)=.TRUE.
        Lstations(ng)=.TRUE.
      END DO
!
!-----------------------------------------------------------------------
!  Initialize several scalar variables.
!-----------------------------------------------------------------------
!
      gorho0=g/rho0
      DO ng=1,Ngrids
        ntfirst(ng)=1
        ntstart(ng)=1
        ntend(ng)=0
        gamma2(ng)=-1.0_r8
        Vtransform(ng)=1
        Vstretching(ng)=1
        first_time(ng)=0
        DO itrc=1,NAT+NPT
          tnu2(itrc,ng)=IniVal
          tnu4(itrc,ng)=IniVal
        END DO
        visc2(ng)=IniVal
        visc4(ng)=IniVal
        DO i=1,4
          FSobc_in (ng,i)=IniVal
          FSobc_out(ng,i)=IniVal
          M2obc_in (ng,i)=IniVal
          M2obc_out(ng,i)=IniVal
          M3obc_in (ng,i)=IniVal
          M3obc_out(ng,i)=IniVal
        END DO
      END DO
!
!  Initialize several IO flags.
!
      DO ng=1,Ngrids
        PerfectRST(ng)=.FALSE.
        Ladjusted(ng)=.FALSE.
        LdefADJ(ng)=.FALSE.
        LdefAVG(ng)=.TRUE.
        LdefAVG2(ng)=.TRUE.
        LdefDIA(ng)=.TRUE.
        LdefFISH(ng)=.TRUE.
        LdefHIS(ng)=.TRUE.
        LdefINI(ng)=.FALSE.
        LdefIRP(ng)=.FALSE.
        LdefITL(ng)=.FALSE.
        LdefMOD(ng)=.FALSE.
        LdefRST(ng)=.TRUE.
        LdefSTA(ng)=.TRUE.
        LdefTLM(ng)=.FALSE.
        LwrtADJ(ng)=.FALSE.
        LwrtAVG(ng)=.FALSE.
        LwrtAVG2(ng)=.FALSE.
        LwrtDIA(ng)=.FALSE.
        LwrtHIS(ng)=.FALSE.
        LwrtPER(ng)=.FALSE.
        LwrtRST(ng)=.FALSE.
        LwrtTLM(ng)=.FALSE.
        LwrtInfo(ng)=.TRUE.
        LwrtState2d(ng)=.FALSE.
        LwrtTime(ng)=.TRUE.
        ldefout(ng)=.FALSE.
        synchro_flag(ng)=.FALSE.
        NrecADJ(ng)=0
        NrecAVG(ng)=0
        NrecAVG2(ng)=0
        NrecDIA(ng)=0
        NrecFIL(ng)=0
        NrecFLT(ng)=0
        NrecINI(ng)=0
        NrecITL(ng)=0
        NrecHIS(ng)=0
        NrecRST(ng)=0
        NrecSTA(ng)=0
        NrecTLM(ng)=0
        NrecNRM(1:4,ng)=0
      END DO
      RETURN
      END SUBROUTINE initialize_scalars
      END MODULE mod_scalars
