!                                                                      !
!  Parameters for UMaine CoSiNE model:                                 !
!                                                                      !
!   reg1     Microzooplankton excretion rate to ammonium [1/day].      !
!   reg2     Mesozooplankton excretion rate to ammonium [1/day].       !
!   gmaxs1   Maximum specific growth rate of small phytoplankton       !
!              [1/day]                                                 !
!   gmaxs2   Maximum specific growth rate of diatom [1/day]            !
!   beta1    Microzooplankton maximum grazing rate [1/day]             !
!   beta2    Mesozooplankton maximum grazing rate [1/day]              !
!   akz1     Half saturation constant for microzooplankton grazing     !
!              [mmol_N/m3]                                             !
!   akz2     Half saturation constant for mesozooplankton grazing      !
!              [mmol_N/m3]                                             !
!   PARfrac  Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!   amaxs2   Initial slope of P-I curve of small phytoplankton         !
!              [1/(Watts/m2)/day]                                      !
!   parsats1 PAR saturation onset parameter of small phytoplankton     !
!              [Watts/m2]                                              !
!   parsats2 PAR saturation onset parameter of diatom [Watts/m2]       !
!              [Watts/m2]                                              !
!   pis1     Ammonium inhibition parameter for small phytoplankton     !
!              [mmol_N/m3]                                             !
!   pis2     Ammonium inhibition parameter for diatom [mmol_N/m3]      !
!   akno3s1  Half saturation concentration for nitrate uptake by       !
!              small phytoplankton [mmol_N/m3].                        !
!   akno3s2  Half saturation concentration for nitrate uptake by       !
!              diatom [mmol_N/m3].                                     !
!   aknh4s1  Half saturation concentration for ammonium uptake by      !
!              small phytoplankton [mmol_N/m3].                        !
!   aknh4s2  Half saturation concentration for ammonium uptake by      !
!              diatom [mmol_N/m3].                                     !
!   akpo4s1  Half saturation concentration for phosphate uptake by     !
!              small phytoplankton [mmol_P/m3].                        !
!   akpo4s2  Half saturation concentration for phosphate uptake by     !
!              diatom [mmol_P/m3].                                     !
!   akco2s1  Half saturation concentration for co2 uptake by           !
!              small phytoplankton [mmol_C/m3].                        !
!   akco2s2  Half saturation concentration for co2 uptake by           !
!              diatom [mmol_C/m3].                                     !
!   aksio4s2 Half saturation constant for silicate uptake by           !
!              diatom [mmol_Si/m3].                                    !
!   ak1      Light attenuation coefficient of water [1/m]              !
!   ak2      Specific light attenuation coefficient for                !
!              phytoplankton [1/m/(mmol_N/m3)].                        !
!   bgamma0   Mesozooplankton specific mortality rate [1/day].         !
!   bgamma1   Grazing efficiency of microzooplankton [nondimensional]. !
!   bgamma2   Grazing efficiency of mesozooplankton [nondimensional].  !
!   bgamma3   Death rate of small phytoplankton [1/day].               !
!   bgamma4   Death rate of large phytoplankton [1/day].               !
!   bgamma5   Decay rate of detritus [1/day].                          !
!   bgamma6                                                            !
!   bgamma7   Nitrafication rate [1/day].                              !
!   wsd      Sinking velocity of detritus [m/day].                     !
!   wsdsi    Sinking velocity of detritus silicate [m/day].            !
!   wsp      Sinking velocity of large phytoplankton [m/day].          !
!   pco2a    Air pCO2 [ppmv].                                          !
!   si2n     Silicate to nitrogen ratio [mol_Si/mol_N].                !
!   p2n      Phosphorus to nitrogen ratio [mol_P/mol_N].               !
!   o2no     Oxygen to nitrate ratio [mol_O2/mol_NO3].                 !
!   o2nh     Oxygen to ammonium ratio [mol_O2/mol_NH4].                !
!   c2n      Carbon to nitrogen ratio [mol_C/mol_N].                   !
!   ro5      Grazing preference for diatom [nondimensional].           !
!   ro6      Grazing preference for mesozooplankton [nondimensional]   !
!   ro7      Grazing preference for detritus [nondimensional].         !
!=======================================================================
!
        USE mod_param

        implicit none
        
        integer, dimension(Ngrids) :: BioIter

        real(r8), dimension(Ngrids) :: reg1          ! 1/day
        real(r8), dimension(Ngrids) :: reg2          ! 1/day
        real(r8), dimension(Ngrids) :: gmaxs1        ! 1/day
        real(r8), dimension(Ngrids) :: gmaxs2        ! 1/day
        real(r8), dimension(Ngrids) :: beta1         ! 1/day
        real(r8), dimension(Ngrids) :: beta2         ! 1/day
        real(r8), dimension(Ngrids) :: akz1          ! mmol_N/m3
        real(r8), dimension(Ngrids) :: akz2          ! mmol_N/m3
        real(r8), dimension(Ngrids) :: PARfrac       ! nondimensional
        real(r8), dimension(Ngrids) :: amaxs2        ! 1/(Watts/m2)/day
        real(r8), dimension(Ngrids) :: pis1          ! m3/mmol_N
        real(r8), dimension(Ngrids) :: pis2          ! m3/mmol_N
        real(r8), dimension(Ngrids) :: akno3s1       ! mmol_N/m3
        real(r8), dimension(Ngrids) :: aknh4s1       ! mmol_N/m3
        real(r8), dimension(Ngrids) :: akpo4s1       ! mmol_P/m3
        real(r8), dimension(Ngrids) :: akco2s1       ! mmol_C/m3
        real(r8), dimension(Ngrids) :: akno3s2       ! mmol_N/m3
        real(r8), dimension(Ngrids) :: aknh4s2       ! mmol_N/m3
        real(r8), dimension(Ngrids) :: aksio4s2      ! mmol_Si/m3
        real(r8), dimension(Ngrids) :: akpo4s2       ! mmol_P/m3
        real(r8), dimension(Ngrids) :: akco2s2       ! mmol_C/m3
        real(r8), dimension(Ngrids) :: ak1           ! 1/m
        real(r8), dimension(Ngrids) :: ak2           ! 1/m/(mmol_N/m3)
        real(r8), dimension(Ngrids) :: parsats1      ! Watts/m2
        real(r8), dimension(Ngrids) :: parsats2      ! Watts/m2
        real(r8), dimension(Ngrids) :: bgamma0       ! 1/day
        real(r8), dimension(Ngrids) :: bgamma1       ! [nondimensional]
        real(r8), dimension(Ngrids) :: bgamma2       ! [nondimensional]
        real(r8), dimension(Ngrids) :: bgamma3       ! 1/day
        real(r8), dimension(Ngrids) :: bgamma4       ! 1/day
        real(r8), dimension(Ngrids) :: bgamma5       ! 1/day
        real(r8), dimension(Ngrids) :: bgamma6       ! 
        real(r8), dimension(Ngrids) :: bgamma7       ! 1/day
        real(r8), dimension(Ngrids) :: wsd           ! m/day
        real(r8), dimension(Ngrids) :: wsdsi         ! m/day
        real(r8), dimension(Ngrids) :: wsp           ! m/day
        real(r8), dimension(Ngrids) :: si2n          ! mol_Si/mol_N
        real(r8), dimension(Ngrids) :: pco2a         ! ppmv
        real(r8), dimension(Ngrids) :: p2n           ! mol_P/mol_N
        real(r8), dimension(Ngrids) :: o2no          ! mol_O2/mol_NO3
        real(r8), dimension(Ngrids) :: o2nh          ! mol_O2/mol_NH4
        real(r8), dimension(Ngrids) :: c2n           ! mol_C/mol_N
        real(r8), dimension(Ngrids) :: ro5           ! nondimensional
        real(r8), dimension(Ngrids) :: ro6           ! nondimensional
        real(r8), dimension(Ngrids) :: ro7           ! nondimensional
