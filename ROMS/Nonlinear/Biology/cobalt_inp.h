      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in biological model input parameters.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval, i, itrc, ng, status,ifield

      integer :: decode_line, load_i, load_l, load_r,load_lbc

      integer ::  igrid, itracer,iTrcStr, iTrcEnd,nline

      logical, dimension(NBT,Ngrids) :: Ltrc
#ifdef DIAGNOSTICS_BIO
      logical, dimension(Ngrids)         :: Lbio
      logical, dimension(NDbio2d,Ngrids) :: Lbio2d
      logical, dimension(NDbio3d,Ngrids) :: Lbio3d
#endif
      real(r8), dimension(NBT,Ngrids) :: Rbio

      real(r8), dimension(200) :: Rval

      character (len=40 ) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(200) :: Cval
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      igrid=1                            ! nested grid counter
      itracer=0                          ! LBC tracer counter
      iTrcStr=1                          ! first LBC tracer to process
      iTrcEnd=NBT                        ! last  LBC tracer to process
      nline=0                            ! LBC multi-line counter
!-----------------------------------------------------------------------
!  Read in Cobalt biological model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)

        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Lbiology')
              Npts=load_l(Nval, Cval, Ngrids, Lbiology)
            CASE ('BioIter')
              Npts=load_i(Nval, Rval, Ngrids, BioIter)
             CASE ('htotal_scale_lo')
               Npts=load_r(Nval, Rval, Ngrids, htotal_scale_lo)
             CASE ('htotal_scale_hi')
               Npts=load_r(Nval, Rval, Ngrids, htotal_scale_hi)
             CASE ('RHO_0')
               Npts=load_r(Nval, Rval, Ngrids, RHO_0)
             CASE ('NKML')
               Npts=load_i(Nval, Rval, Ngrids, NKML)
             CASE ('a_0')
               Npts=load_r(Nval, Rval, Ngrids, a_0)
             CASE ('a_1')
               Npts=load_r(Nval, Rval, Ngrids, a_1)
             CASE ('a_2')
               Npts=load_r(Nval, Rval, Ngrids, a_2)
             CASE ('a_3')
               Npts=load_r(Nval, Rval, Ngrids, a_3)
             CASE ('a_4')
               Npts=load_r(Nval, Rval, Ngrids, a_4)
             CASE ('a_5')
               Npts=load_r(Nval, Rval, Ngrids, a_5)
             CASE ('b_0')
               Npts=load_r(Nval, Rval, Ngrids, b_0)
             CASE ('b_1')
               Npts=load_r(Nval, Rval, Ngrids, b_1)
             CASE ('b_2')
               Npts=load_r(Nval, Rval, Ngrids, b_2)
             CASE ('b_3')
               Npts=load_r(Nval, Rval, Ngrids, b_3)
             CASE ('c_0')
               Npts=load_r(Nval, Rval, Ngrids, c_0)
             CASE ('a1_co2')
               Npts=load_r(Nval, Rval, Ngrids, a1_co2)
             CASE ('a2_co2')
               Npts=load_r(Nval, Rval, Ngrids, a2_co2)
             CASE ('a3_co2')
               Npts=load_r(Nval, Rval, Ngrids, a3_co2)
             CASE ('a4_co2')
               Npts=load_r(Nval, Rval, Ngrids, a4_co2)
             CASE ('a1_o2')
               Npts=load_r(Nval, Rval, Ngrids, a1_o2)
             CASE ('a2_o2')
               Npts=load_r(Nval, Rval, Ngrids, a2_o2)
             CASE ('a3_o2')
               Npts=load_r(Nval, Rval, Ngrids, a3_o2)
             CASE ('a4_o2')
               Npts=load_r(Nval, Rval, Ngrids, a4_o2)
             CASE ('mass_2_n')
               Npts=load_r(Nval, Rval, Ngrids, mass_2_n)
             CASE ('n_2_n_denit')
               Npts=load_r(Nval, Rval, Ngrids, n_2_n_denit)
             CASE ('o2_2_c')
               Npts=load_r(Nval, Rval, Ngrids, o2_2_c)
             CASE ('o2_2_nfix')
               Npts=load_r(Nval, Rval, Ngrids, o2_2_nfix)
             CASE ('o2_2_nh4')
               Npts=load_r(Nval, Rval, Ngrids, o2_2_nh4)
             CASE ('o2_2_nitrif')
               Npts=load_r(Nval, Rval, Ngrids, o2_2_nitrif)
             CASE ('o2_2_no3')
               Npts=load_r(Nval, Rval, Ngrids, o2_2_no3)
             CASE ('k_fed_Di')
               Npts=load_r(Nval, Rval, Ngrids, k_fed_Di)
             CASE ('k_fed_Lg')
               Npts=load_r(Nval, Rval, Ngrids, k_fed_Lg)
             CASE ('k_fed_Sm')
               Npts=load_r(Nval, Rval, Ngrids, k_fed_Sm)
             CASE ('k_nh4_Lg')
               Npts=load_r(Nval, Rval, Ngrids, k_nh4_Lg)
             CASE ('k_nh4_Sm')
               Npts=load_r(Nval, Rval, Ngrids, k_nh4_Sm)
             CASE ('k_nh4_Di')
               Npts=load_r(Nval, Rval, Ngrids, k_nh4_Di)
             CASE ('k_no3_Lg')
               Npts=load_r(Nval, Rval, Ngrids, k_no3_Lg)
             CASE ('k_no3_Sm')
               Npts=load_r(Nval, Rval, Ngrids, k_no3_Sm)
             CASE ('k_no3_Di')
               Npts=load_r(Nval, Rval, Ngrids, k_no3_Di)
             CASE ('k_po4_Di')
               Npts=load_r(Nval, Rval, Ngrids, k_po4_Di)
             CASE ('k_po4_Lg')
               Npts=load_r(Nval, Rval, Ngrids, k_po4_Lg)
             CASE ('k_po4_Sm')
               Npts=load_r(Nval, Rval, Ngrids, k_po4_Sm)
             CASE ('k_sio4_Lg')
               Npts=load_r(Nval, Rval, Ngrids, k_sio4_Lg)
             CASE ('k_fe_2_n_Di')
               Npts=load_r(Nval, Rval, Ngrids, k_fe_2_n_Di)
             CASE ('k_fe_2_n_Lg')
               Npts=load_r(Nval, Rval, Ngrids, k_fe_2_n_Lg)
             CASE ('k_fe_2_n_Sm')
               Npts=load_r(Nval, Rval, Ngrids, k_fe_2_n_Sm)
             CASE ('fe_2_n_max_Sm')
               Npts=load_r(Nval, Rval, Ngrids, fe_2_n_max_Sm)
             CASE ('fe_2_n_max_Lg')
               Npts=load_r(Nval, Rval, Ngrids, fe_2_n_max_Lg)
             CASE ('fe_2_n_max_Di')
               Npts=load_r(Nval, Rval, Ngrids, fe_2_n_max_Di)
             CASE ('fe_2_n_upt_fac')
               Npts=load_r(Nval, Rval, Ngrids, fe_2_n_upt_fac)
             CASE ('alpha_Di')
               Npts=load_r(Nval, Rval, Ngrids, alpha_Di)
             CASE ('alpha_Lg')
               Npts=load_r(Nval, Rval, Ngrids, alpha_Lg)
             CASE ('alpha_Sm')
               Npts=load_r(Nval, Rval, Ngrids, alpha_Sm)
             CASE ('kappa_eppley')
               Npts=load_r(Nval, Rval, Ngrids, kappa_eppley)
             CASE ('P_C_max_Di')
               Npts=load_r(Nval, Rval, Ngrids, P_C_max_Di)
             CASE ('P_C_max_Lg')
               Npts=load_r(Nval, Rval, Ngrids, P_C_max_Lg)
             CASE ('P_C_max_Sm')
               Npts=load_r(Nval, Rval, Ngrids, P_C_max_Sm)
             CASE ('thetamax_Di')
               Npts=load_r(Nval, Rval, Ngrids, thetamax_Di)
             CASE ('thetamax_Lg')
               Npts=load_r(Nval, Rval, Ngrids, thetamax_Lg)
             CASE ('thetamax_Sm')
               Npts=load_r(Nval, Rval, Ngrids, thetamax_Sm)
             CASE ('bresp_Di')
               Npts=load_r(Nval, Rval, Ngrids, bresp_Di)
             CASE ('bresp_Lg')
               Npts=load_r(Nval, Rval, Ngrids, bresp_Lg)
             CASE ('bresp_Sm')
               Npts=load_r(Nval, Rval, Ngrids, bresp_Sm)
             CASE ('thetamin')
               Npts=load_r(Nval, Rval, Ngrids, thetamin)
             CASE ('thetamin_nolim')
               Npts=load_r(Nval, Rval, Ngrids, thetamin_nolim)
             CASE ('zpllgr')
               Npts=load_r(Nval, Rval, Ngrids, zpllgr)
             CASE ('gamma_irr_mem')
               Npts=load_r(Nval, Rval, Ngrids, gamma_irr_mem)
             CASE ('gamma_mu_mem')
               Npts=load_r(Nval, Rval, Ngrids, gamma_mu_mem)
             CASE ('k_n_inhib_Di')
               Npts=load_r(Nval, Rval, Ngrids, k_n_inhib_Di)
             CASE ('o2_inhib_Di_pow')
               Npts=load_r(Nval, Rval, Ngrids, o2_inhib_Di_pow)
             CASE ('o2_inhib_Di_sat')
               Npts=load_r(Nval, Rval, Ngrids, o2_inhib_Di_sat)
             CASE ('p_2_n_static')
               Npts=load_i(Nval, Rval, Ngrids, p_2_n_static)
             CASE ('c_2_n')
               Npts=load_r(Nval, Rval, Ngrids, c_2_n)
             CASE ('alk_2_n_denit')
               Npts=load_r(Nval, Rval, Ngrids, alk_2_n_denit)
             CASE ('p_2_n_static_Di')
               Npts=load_r(Nval, Rval, Ngrids, p_2_n_static_Di)
             CASE ('p_2_n_static_Lg')
               Npts=load_r(Nval, Rval, Ngrids, p_2_n_static_Lg)
             CASE ('p_2_n_static_Sm')
               Npts=load_r(Nval, Rval, Ngrids, p_2_n_static_Sm)
             CASE ('si_2_n_static_Lg')
               Npts=load_r(Nval, Rval, Ngrids, si_2_n_static_Lg)
             CASE ('si_2_n_max_Lg')
               Npts=load_r(Nval, Rval, Ngrids, si_2_n_max_Lg)
             CASE ('ca_2_n_arag')
               Npts=load_r(Nval, Rval, Ngrids, ca_2_n_arag)
             CASE ('ca_2_n_calc')
               Npts=load_r(Nval, Rval, Ngrids, ca_2_n_calc)
             CASE ('caco3_sat_max')
               Npts=load_r(Nval, Rval, Ngrids, caco3_sat_max)
             CASE ('q_p_2_n_smz')
               Npts=load_r(Nval, Rval, Ngrids, q_p_2_n_smz)
             CASE ('q_p_2_n_mdz')
               Npts=load_r(Nval, Rval, Ngrids, q_p_2_n_mdz)
             CASE ('q_p_2_n_lgz')
               Npts=load_r(Nval, Rval, Ngrids, q_p_2_n_lgz)
             CASE ('q_p_2_n_bact')
               Npts=load_r(Nval, Rval, Ngrids, q_p_2_n_bact)
             CASE ('agg_Sm')
               Npts=load_r(Nval, Rval, Ngrids, agg_Sm)
             CASE ('agg_Di')
               Npts=load_r(Nval, Rval, Ngrids, agg_Di)
             CASE ('agg_Lg')
               Npts=load_r(Nval, Rval, Ngrids, agg_Lg)
             CASE ('vir_Sm')
               Npts=load_r(Nval, Rval, Ngrids, vir_Sm)
             CASE ('vir_Di')
               Npts=load_r(Nval, Rval, Ngrids, vir_Di)
             CASE ('vir_Lg')
               Npts=load_r(Nval, Rval, Ngrids, vir_Lg)
             CASE ('vir_Bact')
               Npts=load_r(Nval, Rval, Ngrids, vir_Bact)
             CASE ('ktemp_vir')
               Npts=load_r(Nval, Rval, Ngrids, ktemp_vir)
             CASE ('exu_Sm')
               Npts=load_r(Nval, Rval, Ngrids, exu_Sm)
             CASE ('exu_Di')
               Npts=load_r(Nval, Rval, Ngrids, exu_Di)
             CASE ('exu_Lg')
               Npts=load_r(Nval, Rval, Ngrids, exu_Lg)
             CASE ('imax_smz')
               Npts=load_r(Nval, Rval, Ngrids, imax_smz)
             CASE ('imax_mdz')
               Npts=load_r(Nval, Rval, Ngrids, imax_mdz)
             CASE ('imax_lgz')
               Npts=load_r(Nval, Rval, Ngrids, imax_lgz)
             CASE ('ki_smz')
               Npts=load_r(Nval, Rval, Ngrids, ki_smz)
             CASE ('ki_mdz')
               Npts=load_r(Nval, Rval, Ngrids, ki_mdz)
             CASE ('ki_lgz')
               Npts=load_r(Nval, Rval, Ngrids, ki_lgz)
             CASE ('ktemp_smz')
               Npts=load_r(Nval, Rval, Ngrids, ktemp_smz)
             CASE ('ktemp_mdz')
               Npts=load_r(Nval, Rval, Ngrids, ktemp_mdz)
             CASE ('ktemp_lgz')
               Npts=load_r(Nval, Rval, Ngrids, ktemp_lgz)
             CASE ('mu_max_bact')
               Npts=load_r(Nval, Rval, Ngrids, mu_max_bact)
             CASE ('k_ldon_bact')
               Npts=load_r(Nval, Rval, Ngrids, k_ldon_bact)
             CASE ('ktemp_bact')
               Npts=load_r(Nval, Rval, Ngrids, ktemp_bact)
             CASE ('nswitch_smz')
               Npts=load_r(Nval, Rval, Ngrids, nswitch_smz)
             CASE ('nswitch_mdz')
               Npts=load_r(Nval, Rval, Ngrids, nswitch_mdz)
             CASE ('nswitch_lgz')
               Npts=load_r(Nval, Rval, Ngrids, nswitch_lgz)
             CASE ('mswitch_smz')
               Npts=load_r(Nval, Rval, Ngrids, mswitch_smz)
             CASE ('mswitch_mdz')
               Npts=load_r(Nval, Rval, Ngrids, mswitch_mdz)
             CASE ('mswitch_lgz')
               Npts=load_r(Nval, Rval, Ngrids, mswitch_lgz)
             CASE ('smz_ipa_smp')
               Npts=load_r(Nval, Rval, Ngrids, smz_ipa_smp)
             CASE ('smz_ipa_lgp')
               Npts=load_r(Nval, Rval, Ngrids, smz_ipa_lgp)
             CASE ('smz_ipa_diaz')
               Npts=load_r(Nval, Rval, Ngrids, smz_ipa_diaz)
             CASE ('smz_ipa_smz')
               Npts=load_r(Nval, Rval, Ngrids, smz_ipa_smz)
             CASE ('smz_ipa_mdz')
               Npts=load_r(Nval, Rval, Ngrids, smz_ipa_mdz)
             CASE ('smz_ipa_lgz')
               Npts=load_r(Nval, Rval, Ngrids, smz_ipa_lgz)
             CASE ('smz_ipa_bact')
               Npts=load_r(Nval, Rval, Ngrids, smz_ipa_bact)
             CASE ('smz_ipa_det')
               Npts=load_r(Nval, Rval, Ngrids, smz_ipa_det)
             CASE ('mdz_ipa_smp')
               Npts=load_r(Nval, Rval, Ngrids, mdz_ipa_smp)
             CASE ('mdz_ipa_lgp')
               Npts=load_r(Nval, Rval, Ngrids, mdz_ipa_lgp)
             CASE ('mdz_ipa_diaz')
               Npts=load_r(Nval, Rval, Ngrids, mdz_ipa_diaz)
             CASE ('mdz_ipa_smz')
               Npts=load_r(Nval, Rval, Ngrids, mdz_ipa_smz)
             CASE ('mdz_ipa_mdz')
               Npts=load_r(Nval, Rval, Ngrids, mdz_ipa_mdz)
             CASE ('mdz_ipa_lgz')
               Npts=load_r(Nval, Rval, Ngrids, mdz_ipa_lgz)
             CASE ('mdz_ipa_bact')
               Npts=load_r(Nval, Rval, Ngrids, mdz_ipa_bact)
             CASE ('mdz_ipa_det')
               Npts=load_r(Nval, Rval, Ngrids, mdz_ipa_det)
             CASE ('lgz_ipa_smp')
               Npts=load_r(Nval, Rval, Ngrids, lgz_ipa_smp)
             CASE ('lgz_ipa_lgp')
               Npts=load_r(Nval, Rval, Ngrids, lgz_ipa_lgp)
             CASE ('lgz_ipa_diaz')
               Npts=load_r(Nval, Rval, Ngrids, lgz_ipa_diaz)
             CASE ('lgz_ipa_smz')
               Npts=load_r(Nval, Rval, Ngrids, lgz_ipa_smz)
             CASE ('lgz_ipa_mdz')
               Npts=load_r(Nval, Rval, Ngrids, lgz_ipa_mdz)
             CASE ('lgz_ipa_lgz')
               Npts=load_r(Nval, Rval, Ngrids, lgz_ipa_lgz)
             CASE ('lgz_ipa_bact')
               Npts=load_r(Nval, Rval, Ngrids, lgz_ipa_bact)
             CASE ('lgz_ipa_det')
               Npts=load_r(Nval, Rval, Ngrids, lgz_ipa_det)
             CASE ('gge_max_smz')
               Npts=load_r(Nval, Rval, Ngrids, gge_max_smz)
             CASE ('gge_max_mdz')
               Npts=load_r(Nval, Rval, Ngrids, gge_max_mdz)
             CASE ('gge_max_lgz')
               Npts=load_r(Nval, Rval, Ngrids, gge_max_lgz)
             CASE ('bresp_smz')
               Npts=load_r(Nval, Rval, Ngrids, bresp_smz)
             CASE ('bresp_mdz')
               Npts=load_r(Nval, Rval, Ngrids, bresp_mdz)
             CASE ('bresp_lgz')
               Npts=load_r(Nval, Rval, Ngrids, bresp_lgz)
             CASE ('gge_max_bact')
               Npts=load_r(Nval, Rval, Ngrids, gge_max_bact)
             CASE ('bresp_bact')
               Npts=load_r(Nval, Rval, Ngrids, bresp_bact)
             CASE ('phi_det_smz')
               Npts=load_r(Nval, Rval, Ngrids, phi_det_smz)
             CASE ('phi_det_mdz')
               Npts=load_r(Nval, Rval, Ngrids, phi_det_mdz)
             CASE ('phi_det_lgz')
               Npts=load_r(Nval, Rval, Ngrids, phi_det_lgz)
             CASE ('phi_ldon_smz')
               Npts=load_r(Nval, Rval, Ngrids, phi_ldon_smz)
             CASE ('phi_ldon_mdz')
               Npts=load_r(Nval, Rval, Ngrids, phi_ldon_mdz)
             CASE ('phi_ldon_lgz')
               Npts=load_r(Nval, Rval, Ngrids, phi_ldon_lgz)
             CASE ('phi_ldop_smz')
               Npts=load_r(Nval, Rval, Ngrids, phi_ldop_smz)
             CASE ('phi_ldop_mdz')
               Npts=load_r(Nval, Rval, Ngrids, phi_ldop_mdz)
             CASE ('phi_ldop_lgz')
               Npts=load_r(Nval, Rval, Ngrids, phi_ldop_lgz)
             CASE ('phi_srdon_smz')
               Npts=load_r(Nval, Rval, Ngrids, phi_srdon_smz)
             CASE ('phi_srdon_mdz')
               Npts=load_r(Nval, Rval, Ngrids, phi_srdon_mdz)
             CASE ('phi_srdon_lgz')
               Npts=load_r(Nval, Rval, Ngrids, phi_srdon_lgz)
             CASE ('phi_srdop_smz')
               Npts=load_r(Nval, Rval, Ngrids, phi_srdop_smz)
             CASE ('phi_srdop_mdz')
               Npts=load_r(Nval, Rval, Ngrids, phi_srdop_mdz)
             CASE ('phi_srdop_lgz')
               Npts=load_r(Nval, Rval, Ngrids, phi_srdop_lgz)
             CASE ('phi_sldon_smz')
               Npts=load_r(Nval, Rval, Ngrids, phi_sldon_smz)
             CASE ('phi_sldon_mdz')
               Npts=load_r(Nval, Rval, Ngrids, phi_sldon_mdz)
             CASE ('phi_sldon_lgz')
               Npts=load_r(Nval, Rval, Ngrids, phi_sldon_lgz)
             CASE ('phi_sldop_smz')
               Npts=load_r(Nval, Rval, Ngrids, phi_sldop_smz)
             CASE ('phi_sldop_mdz')
               Npts=load_r(Nval, Rval, Ngrids, phi_sldop_mdz)
             CASE ('phi_sldop_lgz')
               Npts=load_r(Nval, Rval, Ngrids, phi_sldop_lgz)
             CASE ('phi_nh4_smz')
               Npts=load_r(Nval, Rval, Ngrids, phi_nh4_smz)
             CASE ('phi_nh4_mdz')
               Npts=load_r(Nval, Rval, Ngrids, phi_nh4_mdz)
             CASE ('phi_nh4_lgz')
               Npts=load_r(Nval, Rval, Ngrids, phi_nh4_lgz)
             CASE ('phi_po4_smz')
               Npts=load_r(Nval, Rval, Ngrids, phi_po4_smz)
             CASE ('phi_po4_mdz')
               Npts=load_r(Nval, Rval, Ngrids, phi_po4_mdz)
             CASE ('phi_po4_lgz')
               Npts=load_r(Nval, Rval, Ngrids, phi_po4_lgz)
             CASE ('phi_ldon_vir')
               Npts=load_r(Nval, Rval, Ngrids, phi_ldon_vir)
             CASE ('phi_srdon_vir')
               Npts=load_r(Nval, Rval, Ngrids, phi_srdon_vir)
             CASE ('phi_sldon_vir')
               Npts=load_r(Nval, Rval, Ngrids, phi_sldon_vir)
             CASE ('phi_ldop_vir')
               Npts=load_r(Nval, Rval, Ngrids, phi_ldop_vir)
             CASE ('phi_srdop_vir')
               Npts=load_r(Nval, Rval, Ngrids, phi_srdop_vir)
             CASE ('phi_sldop_vir')
               Npts=load_r(Nval, Rval, Ngrids, phi_sldop_vir)
             CASE ('imax_hp')
               Npts=load_r(Nval, Rval, Ngrids, imax_hp)
             CASE ('ki_hp')
               Npts=load_r(Nval, Rval, Ngrids, ki_hp)
             CASE ('coef_hp')
               Npts=load_r(Nval, Rval, Ngrids, coef_hp)
             CASE ('ktemp_hp')
               Npts=load_r(Nval, Rval, Ngrids, ktemp_hp)
             CASE ('nswitch_hp')
               Npts=load_r(Nval, Rval, Ngrids, nswitch_hp)
             CASE ('mswitch_hp')
               Npts=load_r(Nval, Rval, Ngrids, mswitch_hp)
             CASE ('hp_ipa_smp')
               Npts=load_r(Nval, Rval, Ngrids, hp_ipa_smp)
             CASE ('hp_ipa_lgp')
               Npts=load_r(Nval, Rval, Ngrids, hp_ipa_lgp)
             CASE ('hp_ipa_diaz')
               Npts=load_r(Nval, Rval, Ngrids, hp_ipa_diaz)
             CASE ('hp_ipa_smz')
               Npts=load_r(Nval, Rval, Ngrids, hp_ipa_smz)
             CASE ('hp_ipa_mdz')
               Npts=load_r(Nval, Rval, Ngrids, hp_ipa_mdz)
             CASE ('hp_ipa_lgz')
               Npts=load_r(Nval, Rval, Ngrids, hp_ipa_lgz)
             CASE ('hp_ipa_bact')
               Npts=load_r(Nval, Rval, Ngrids, hp_ipa_bact)
             CASE ('hp_ipa_det')
               Npts=load_r(Nval, Rval, Ngrids, hp_ipa_det)
             CASE ('hp_phi_det')
               Npts=load_r(Nval, Rval, Ngrids, hp_phi_det)
             CASE ('hp_phi_ldon')
               Npts=load_r(Nval, Rval, Ngrids, hp_phi_ldon)
             CASE ('hp_phi_ldop')
               Npts=load_r(Nval, Rval, Ngrids, hp_phi_ldop)
             CASE ('hp_phi_srdon')
               Npts=load_r(Nval, Rval, Ngrids, hp_phi_srdon)
             CASE ('hp_phi_srdop')
               Npts=load_r(Nval, Rval, Ngrids, hp_phi_srdop)
             CASE ('hp_phi_sldon')
               Npts=load_r(Nval, Rval, Ngrids, hp_phi_sldon)
             CASE ('hp_phi_sldop')
               Npts=load_r(Nval, Rval, Ngrids, hp_phi_sldop)
             CASE ('hp_phi_nh4')
               Npts=load_r(Nval, Rval, Ngrids, hp_phi_nh4)
             CASE ('hp_phi_po4')
               Npts=load_r(Nval, Rval, Ngrids, hp_phi_po4)
             CASE ('felig_bkg')
               Npts=load_r(Nval, Rval, Ngrids, felig_bkg)
             CASE ('felig_2_don')
               Npts=load_r(Nval, Rval, Ngrids, felig_2_don)
             CASE ('fe_2_n_sed')
               Npts=load_r(Nval, Rval, Ngrids, fe_2_n_sed)
             CASE ('fe_coast')
               Npts=load_r(Nval, Rval, Ngrids, fe_coast)
             CASE ('alpha_fescav')
               Npts=load_r(Nval, Rval, Ngrids, alpha_fescav)
             CASE ('beta_fescav')
               Npts=load_r(Nval, Rval, Ngrids, beta_fescav)
             CASE ('remin_eff_fedet')
               Npts=load_r(Nval, Rval, Ngrids, remin_eff_fedet)
             CASE ('ki_fescav')
               Npts=load_r(Nval, Rval, Ngrids, ki_fescav)
             CASE ('io_fescav')
               Npts=load_r(Nval, Rval, Ngrids, io_fescav)
             CASE ('gamma_fescav')
               Npts=load_r(Nval, Rval, Ngrids, gamma_fescav)
             CASE ('kfe_eq_lig_ll')
               Npts=load_r(Nval, Rval, Ngrids, kfe_eq_lig_ll)
             CASE ('kfe_eq_lig_hl')
               Npts=load_r(Nval, Rval, Ngrids, kfe_eq_lig_hl)
             CASE ('k_o2')
               Npts=load_r(Nval, Rval, Ngrids, k_o2)
             CASE ('o2_min')
               Npts=load_r(Nval, Rval, Ngrids, o2_min)
             CASE ('rpcaco3')
               Npts=load_r(Nval, Rval, Ngrids, rpcaco3)
             CASE ('rplith')
               Npts=load_r(Nval, Rval, Ngrids, rplith)
             CASE ('rpsio2')
               Npts=load_r(Nval, Rval, Ngrids, rpsio2)
             CASE ('gamma_ndet')
               Npts=load_r(Nval, Rval, Ngrids, gamma_ndet)
             CASE ('gamma_cadet_arag')
               Npts=load_r(Nval, Rval, Ngrids, gamma_cadet_arag)
             CASE ('gamma_cadet_calc')
               Npts=load_r(Nval, Rval, Ngrids, gamma_cadet_calc)
             CASE ('gamma_sidet')
               Npts=load_r(Nval, Rval, Ngrids, gamma_sidet)
             CASE ('phi_lith')
               Npts=load_r(Nval, Rval, Ngrids, phi_lith)
             CASE ('k_lith')
               Npts=load_r(Nval, Rval, Ngrids, k_lith)
             CASE ('z_sed')
               Npts=load_r(Nval, Rval, Ngrids, z_sed)
             CASE ('k_no3_denit')
               Npts=load_r(Nval, Rval, Ngrids, k_no3_denit)
             CASE ('gamma_srdon')
               Npts=load_r(Nval, Rval, Ngrids, gamma_srdon)
             CASE ('gamma_srdop')
               Npts=load_r(Nval, Rval, Ngrids, gamma_srdop)
             CASE ('gamma_sldon')
               Npts=load_r(Nval, Rval, Ngrids, gamma_sldon)
             CASE ('gamma_sldop')
               Npts=load_r(Nval, Rval, Ngrids, gamma_sldop)
             CASE ('gamma_nitrif')
               Npts=load_r(Nval, Rval, Ngrids, gamma_nitrif)
             CASE ('irr_inhibit')
               Npts=load_r(Nval, Rval, Ngrids, irr_inhibit)
             CASE ('htotal_in')
               Npts=load_r(Nval, Rval, Ngrids, htotal_in)
             CASE ('wsink')
               Npts=load_r(Nval, Rval, Ngrids, wsink)
#ifdef COASTDIAT
             CASE ('k_fed_Md')
               Npts=load_r(Nval, Rval, Ngrids, k_fed_Md)
             CASE ('k_nh4_Md')
               Npts=load_r(Nval, Rval, Ngrids, k_nh4_Md)
             CASE ('k_no3_Md')
               Npts=load_r(Nval, Rval, Ngrids, k_no3_Md)
             CASE ('k_po4_Md')
               Npts=load_r(Nval, Rval, Ngrids, k_po4_Md)
             CASE ('k_sio4_Md')
               Npts=load_r(Nval, Rval, Ngrids, k_sio4_Md)
             CASE ('k_fe_2_n_Md')
               Npts=load_r(Nval, Rval, Ngrids, k_fe_2_n_Md)
             CASE ('fe_2_n_max_Md')
               Npts=load_r(Nval, Rval, Ngrids, fe_2_n_max_Md)
             CASE ('alpha_Md')
               Npts=load_r(Nval, Rval, Ngrids, alpha_Md)
             CASE ('P_C_max_Md')
               Npts=load_r(Nval, Rval, Ngrids, P_C_max_Md)
             CASE ('thetamax_Md')
               Npts=load_r(Nval, Rval, Ngrids, thetamax_Md)
             CASE ('bresp_Md')
               Npts=load_r(Nval, Rval, Ngrids, bresp_Md)
             CASE ('p_2_n_static_Md')
               Npts=load_r(Nval, Rval, Ngrids, p_2_n_static_Md)
             CASE ('si_2_n_static_Md')
               Npts=load_r(Nval, Rval, Ngrids, si_2_n_static_Md)
             CASE ('si_2_n_max_Md')
               Npts=load_r(Nval, Rval, Ngrids, si_2_n_max_Md)
             CASE ('agg_Md')
               Npts=load_r(Nval, Rval, Ngrids, agg_Md)
             CASE ('vir_Md')
               Npts=load_r(Nval, Rval, Ngrids, vir_Md)
             CASE ('exu_Md')
               Npts=load_r(Nval, Rval, Ngrids, exu_Md)
             CASE ('smz_ipa_mdp')
               Npts=load_r(Nval, Rval, Ngrids, smz_ipa_mdp)
             CASE ('mdz_ipa_mdp')
               Npts=load_r(Nval, Rval, Ngrids, mdz_ipa_mdp)
             CASE ('lgz_ipa_mdp')
               Npts=load_r(Nval, Rval, Ngrids, lgz_ipa_mdp)
#endif
            CASE ('TNU2')
              Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  tnu2(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('TNU4')
              Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  tnu4(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('AKT_BAK')
              Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  Akt_bak(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ('TNUDG')
              Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  Tnudg(i,ng)=Rbio(itrc,ng)
                END DO
              END DO
            CASE ( 'LBC(isTvar)')
              IF (itracer.lt.NBT) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(idbio(itracer))
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      iTrcStr, iTrcEnd,                           &
     &                      Vname(1,idTvar(idbio(itracer))), LBC)
            CASE ('LtracerCLM')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerCLM(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO

            CASE ('LnudgeTCLM')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LnudgeTCLM(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('LtracerSrc')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idbio(itrc)
                  LtracerSrc(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idTvar)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTvar(idbio(itrc))
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,120)                           &
       &                      'idTvar(idbio(', itrc, '))'
                    exit_flag=5
                    RETURN
                  END IF
                  Hout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idTsur)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTsur(idbio(itrc))
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,120)                           &
       &                      'idTsur(idbio(', itrc, '))'
                    exit_flag=5
                    RETURN
                  END IF
                  Hout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
#ifdef BENTHIC
          CASE ('Hout(idBeTvar)')
            Npts=load_l(Nval, Cval,NBEN*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBEN
                i=idBeTvar(idben(itrc))
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
#endif
#ifdef DIAGNOSTICS_BIO
            CASE ('Dout(iDbio2)')
              Npts=load_l(Nval, Cval, NDbio2d*Ngrids, Lbio2d)
              DO ng=1,Ngrids
                DO itrc=1,NDbio2d
                  i=iDbio2(itrc)
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,120)                           &
       &                      'iDbio2(', itrc, ')'
                    exit_flag=5
                    RETURN
                  END IF
                  Dout(i,ng)=Lbio2d(itrc,ng)
                END DO
              END DO
! obsolete but kept for backward compatibility with old namelist (for now)
            CASE ('Dout(iDbio3)')
              Npts=load_l(Nval, Cval, NDbio3d*Ngrids, Lbio3d)
              DO ng=1,Ngrids
                DO itrc=1,NDbio3d
                  i=iDbio3(itrc)
                  IF (i.eq.0) THEN
                    IF (Master) WRITE (out,120)                           &
       &                      'iDbio3(', itrc, ')'
                    exit_flag=5
                    RETURN
                  END IF
                  Dout(i,ng)=Lbio3d(itrc,ng)
                END DO
              END DO
! starts here
            CASE ('Dout(ichl)')
              IF (iDbio3(ichl).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ichl)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ichl)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iirr_mem)')
              IF (iDbio3(iirr_mem).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iirr_mem)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iirr_mem)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ico3_ion)')
              IF (iDbio3(ico3_ion).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ico3_ion)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ico3_ion)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ihtotal)')
              IF (iDbio3(ihtotal).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ihtotal)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ihtotal)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ico3_sol_arag)')
              IF (iDbio3(ico3_sol_arag).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ico3_sol_arag)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ico3_sol_arag)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ico3_sol_calc)')
              IF (iDbio3(ico3_sol_calc).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ico3_sol_calc)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ico3_sol_calc)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iomega_cadet_arag)')
              IF (iDbio3(iomega_cadet_arag).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iomega_cadet_arag)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iomega_cadet_arag)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iomega_cadet_calc)')
              IF (iDbio3(iomega_cadet_calc).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iomega_cadet_calc)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iomega_cadet_calc)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iswdk)')
              IF (iDbio3(iswdk).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iswdk)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iswdk)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iirr_mix)')
              IF (iDbio3(iirr_mix).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iirr_mix)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iirr_mix)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iirr_inst)')
              IF (iDbio3(iirr_inst).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iirr_inst)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iirr_inst)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ife_bulk_flx)')
              IF (iDbio3(ife_bulk_flx).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ife_bulk_flx)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ife_bulk_flx)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(imu_mem_sm)')
              IF (iDbio3(imu_mem_sm).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(imu_mem_sm)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(imu_mem_sm)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(imu_mem_di)')
              IF (iDbio3(imu_mem_di).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(imu_mem_di)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(imu_mem_di)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(imu_mem_lg)')
              IF (iDbio3(imu_mem_lg).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(imu_mem_lg)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(imu_mem_lg)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iagg_lim_sm)')
              IF (iDbio3(iagg_lim_sm).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iagg_lim_sm)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iagg_lim_sm)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iagg_lim_di)')
              IF (iDbio3(iagg_lim_di).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iagg_lim_di)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iagg_lim_di)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iagg_lim_lg)')
              IF (iDbio3(iagg_lim_lg).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iagg_lim_lg)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iagg_lim_lg)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iaggloss_sm)')
              IF (iDbio3(iaggloss_sm).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iaggloss_sm)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iaggloss_sm)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iaggloss_di)')
              IF (iDbio3(iaggloss_di).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iaggloss_di)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iaggloss_di)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iaggloss_lg)')
              IF (iDbio3(iaggloss_lg).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iaggloss_lg)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iaggloss_lg)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ivirloss_sm)')
              IF (iDbio3(ivirloss_sm).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ivirloss_sm)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ivirloss_sm)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ivirloss_di)')
              IF (iDbio3(ivirloss_di).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ivirloss_di)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ivirloss_di)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ivirloss_lg)')
              IF (iDbio3(ivirloss_lg).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ivirloss_lg)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ivirloss_lg)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(izloss_sm)')
              IF (iDbio3(izloss_sm).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(izloss_sm)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(izloss_sm)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(izloss_di)')
              IF (iDbio3(izloss_di).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(izloss_di)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(izloss_di)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(izloss_lg)')
              IF (iDbio3(izloss_lg).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(izloss_lg)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(izloss_lg)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(idef_fe_sm)')
              IF (iDbio3(idef_fe_sm).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(idef_fe_sm)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(idef_fe_sm)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(idef_fe_di)')
              IF (iDbio3(idef_fe_di).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(idef_fe_di)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(idef_fe_di)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(idef_fe_lg)')
              IF (iDbio3(idef_fe_lg).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(idef_fe_lg)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(idef_fe_lg)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ifelim_sm)')
              IF (iDbio3(ifelim_sm).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ifelim_sm)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ifelim_sm)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ifelim_di)')
              IF (iDbio3(ifelim_di).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ifelim_di)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ifelim_di)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ifelim_lg)')
              IF (iDbio3(ifelim_lg).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ifelim_lg)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ifelim_lg)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ino3lim_sm)')
              IF (iDbio3(ino3lim_sm).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ino3lim_sm)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ino3lim_sm)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ino3lim_di)')
              IF (iDbio3(ino3lim_di).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ino3lim_di)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ino3lim_di)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ino3lim_lg)')
              IF (iDbio3(ino3lim_lg).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ino3lim_lg)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ino3lim_lg)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(inh4lim_sm)')
              IF (iDbio3(inh4lim_sm).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(inh4lim_sm)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(inh4lim_sm)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(inh4lim_di)')
              IF (iDbio3(inh4lim_di).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(inh4lim_di)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(inh4lim_di)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(inh4lim_lg)')
              IF (iDbio3(inh4lim_lg).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(inh4lim_lg)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(inh4lim_lg)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ipo4lim_sm)')
              IF (iDbio3(ipo4lim_sm).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ipo4lim_sm)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ipo4lim_sm)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ipo4lim_di)')
              IF (iDbio3(ipo4lim_di).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ipo4lim_di)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ipo4lim_di)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ipo4lim_lg)')
              IF (iDbio3(ipo4lim_lg).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ipo4lim_lg)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ipo4lim_lg)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ichl_di)')
              IF (iDbio3(ichl_di).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ichl_di)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ichl_di)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iC_2_chl_di)')
              IF (iDbio3(iC_2_chl_di).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iC_2_chl_di)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iC_2_chl_di)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ichl_sm)')
              IF (iDbio3(ichl_sm).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ichl_sm)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ichl_sm)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iC_2_chl_sm)')
              IF (iDbio3(iC_2_chl_sm).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iC_2_chl_sm)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iC_2_chl_sm)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ichl_lg)')
              IF (iDbio3(ichl_lg).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ichl_lg)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ichl_lg)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iC_2_chl_lg)')
              IF (iDbio3(iC_2_chl_lg).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iC_2_chl_lg)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iC_2_chl_lg)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
# ifdef COASTDIAT
            CASE ('Dout(imu_mem_md)')
              IF (iDbio3(imu_mem_md).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(imu_mem_md)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(imu_mem_md)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iaggloss_md)')
              IF (iDbio3(iaggloss_md).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iaggloss_md)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iaggloss_md)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ivirloss_md)')
              IF (iDbio3(ivirloss_md).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ivirloss_md)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ivirloss_md)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(izloss_md)')
              IF (iDbio3(izloss_md).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(izloss_md)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(izloss_md)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iagg_lim_md)')
              IF (iDbio3(iagg_lim_md).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iagg_lim_md)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iagg_lim_md)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(idef_fe_md)')
              IF (iDbio3(idef_fe_md).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(idef_fe_md)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(idef_fe_md)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ifelim_md)')
              IF (iDbio3(ifelim_md).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ifelim_md)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ifelim_md)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ino3lim_md)')
              IF (iDbio3(ino3lim_md).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ino3lim_md)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ino3lim_md)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(inh4lim_md)')
              IF (iDbio3(inh4lim_md).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(inh4lim_md)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(inh4lim_md)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ipo4lim_md)')
              IF (iDbio3(ipo4lim_md).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ipo4lim_md)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ipo4lim_md)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(ichl_md)')
              IF (iDbio3(ichl_md).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(ichl_md)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(ichl_md)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
            CASE ('Dout(iC_2_chl_md)')
              IF (iDbio3(iC_2_chl_md).eq.0) THEN
                IF (Master) WRITE (out,120) 'iDbio3(iC_2_chl_md)'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Lbio)
              i=iDbio3(iC_2_chl_md)
              DO ng=1,Ngrids
                Dout(i,ng)=Lbio(ng)
              END DO
# endif
!            CASE ('Dout(inewdiag)')
!              IF (iDbio3(inewdiag).eq.0) THEN
!                IF (Master) WRITE (out,120) 'iDbio3(inewdiag)'
!                exit_flag=5
!                RETURN
!              END IF
!              Npts=load_l(Nval, Cval, Ngrids, Lbio)
!              i=iDbio3(inewdiag)
!              DO ng=1,Ngrids
!                Dout(i,ng)=Lbio(ng)
!              END DO
#endif
#ifdef PRIMARY_PROD
            CASE ('Hout(idNPP)')
              IF (idNPP.eq.0) THEN
                IF (Master) WRITE (out,120) 'idNPP'
                exit_flag=5
                RETURN
              END IF
              Npts=load_l(Nval, Cval, Ngrids, Hout(idNPP,:))
#endif
#ifdef AVERAGES
            CASE ('Aout(idTvar)')
              Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
              DO ng=1,Ngrids
                DO itrc=1,NBT
                  i=idTvar(idbio(itrc))
                  Aout(i,ng)=Ltrc(itrc,ng)
                END DO
              END DO
#endif
# ifdef PRIMARY_PROD
            CASE ('Aout(idNPP)')
              Npts=load_l(Nval, Cval, Ngrids, Aout(idNPP,:))
# endif
# ifdef BENTHIC
            CASE ('Aout(idBeTvar)')
            Npts=load_l(Nval, Cval,NBEN*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBEN
                i=idBeTvar(idben(itrc))
                Aout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
# endif

          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,30) line
      exit_flag=4
      RETURN
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,40) ng
            WRITE (out,50) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
!             WRITE (out,100) init(ng), 'init',                          &
!     &            'initialisation[N/A].'
             WRITE (out,100) htotal_scale_lo(ng), 'htotal_scale_lo',    &
     &            'htotal_scale_lo ?[N/A].'
             WRITE (out,100) htotal_scale_hi(ng), 'htotal_scale_hi',    &
     &            'htotal_scale_hi ?[N/A].'
             WRITE (out,100) RHO_0(ng), 'RHO_0',                        &
     &            'sea water density[N/A].'
! RD: not a real, change the format
!             WRITE (out,100) NKML(ng), 'NKML',                          &
!     &            'NKML ?[N/A].'
             WRITE (out,100) a_0(ng), 'a_0',                            &
     &            'a_0 : coefficients for O2 saturation[N/A].'
             WRITE (out,100) a_1(ng), 'a_1',                            &
     &            'a_1 : coefficients for O2 saturation[N/A].'
             WRITE (out,100) a_2(ng), 'a_2',                            &
     &            'a_2 : coefficients for O2 saturation[N/A].'
             WRITE (out,100) a_3(ng), 'a_3',                            &
     &            'a_3 : coefficients for O2 saturation[N/A].'
             WRITE (out,100) a_4(ng), 'a_4',                            &
     &            'a_4 : coefficients for O2 saturation[N/A].'
             WRITE (out,100) a_5(ng), 'a_5',                            &
     &            'a_5 : coefficients for O2 saturation[N/A].'
             WRITE (out,100) b_0(ng), 'b_0',                            &
     &            'b_0 : coefficients for O2 saturation[N/A].'
             WRITE (out,100) b_1(ng), 'b_1',                            &
     &            'b_1 : coefficients for O2 saturation[N/A].'
             WRITE (out,100) b_2(ng), 'b_2',                            &
     &            'b_2 : coefficients for O2 saturation[N/A].'
             WRITE (out,100) b_3(ng), 'b_3',                            &
     &            'b_3 : coefficients for O2 saturation[N/A].'
             WRITE (out,100) c_0(ng), 'c_0',                            &
     &            'c_0 : coefficients for O2 saturation[N/A].'
             WRITE (out,100) a1_co2(ng), 'a1_co2',                      &
     &            'a1_co2 : Compute the Schmidt number of CO2 in        &
     &            seawater[N/A].'
             WRITE (out,100) a2_co2(ng), 'a2_co2',                      &
     &            'a2_co2 : Compute the Schmidt number of CO2 in        &
     &            seawater[N/A].'
             WRITE (out,100) a3_co2(ng), 'a3_co2',                      &
     &            'a3_co2 : Compute the Schmidt number of CO2 in        &
     &            seawater[N/A].'
             WRITE (out,100) a4_co2(ng), 'a4_co2',                      &
     &            'a4_co2 : Compute the Schmidt number of CO2 in        &
     &            seawater[N/A].'
             WRITE (out,100) a1_o2(ng), 'a1_o2',                        &
     &            'a1_o2 : Compute the Schmidt number of O2 in          &
     &            seawater[N/A].'
             WRITE (out,100) a2_o2(ng), 'a2_o2',                        &
     &            'a2_o2 : Compute the Schmidt number of O2 in          &
     &            seawater[N/A].'
             WRITE (out,100) a3_o2(ng), 'a3_o2',                        &
     &            'a3_o2 : Compute the Schmidt number of O2 in          &
     &            seawater[N/A].'
             WRITE (out,100) a4_o2(ng), 'a4_o2',                        &
     &            'a4_o2 : Compute the Schmidt number of O2 in          &
     &            seawater[N/A].'
             WRITE (out,100) mass_2_n(ng), 'mass_2_n',                  &
     &            'mass_2_n: Stoichiometry[g mol N-1].'
             WRITE (out,100) n_2_n_denit(ng), 'n_2_n_denit',            &
     &            'n_2_n_denit: Stoichiometry[mol N NO3 mol N org-1].'
             WRITE (out,100) o2_2_c(ng), 'o2_2_c',                      &
     &            'o2_2_c: Stoichiometry[mol O2 mol C-1].'
             WRITE (out,100) o2_2_nfix(ng), 'o2_2_nfix',                &
     &            'o2_2_nfix: Stoichiometry[mol O2 mol N-1].'
             WRITE (out,100) o2_2_nh4(ng), 'o2_2_nh4',                  &
     &            'o2_2_nh4: Stoichiometry[mol O2 mol N-1].'
             WRITE (out,100) o2_2_nitrif(ng), 'o2_2_nitrif',            &
     &            'o2_2_nitrif: Stoichiometry[mol O2 mol N-1].'
             WRITE (out,100) o2_2_no3(ng), 'o2_2_no3',                  &
     &            'o2_2_no3: Stoichiometry[mol O2 mol N-1].'
             WRITE (out,100) k_fed_Di(ng), 'k_fed_Di',                  &
     &            'k_fed_Di: Nutrient Limitation Parameters             &
     &            (phytoplankton)[mol Fed kg-1].'
             WRITE (out,100) k_fed_Lg(ng), 'k_fed_Lg',                  &
     &            'k_fed_Lg: Nutrient Limitation Parameters             &
     &            (phytoplankton)[mol Fed kg-1].'
             WRITE (out,100) k_fed_Sm(ng), 'k_fed_Sm',                  &
     &            'k_fed_Sm: Nutrient Limitation Parameters             &
     &            (phytoplankton)[mol Fed kg-1].'
             WRITE (out,100) k_nh4_Lg(ng), 'k_nh4_Lg',                  &
     &            'k_nh4_Lg: Nutrient Limitation Parameters             &
     &            (phytoplankton)[mol NH4 kg-1].'
             WRITE (out,100) k_nh4_Sm(ng), 'k_nh4_Sm',                  &
     &            'k_nh4_Sm: Nutrient Limitation Parameters             &
     &            (phytoplankton)[mol NH4 kg-1].'
             WRITE (out,100) k_nh4_Di(ng), 'k_nh4_Di',                  &
     &            'k_nh4_Di: Nutrient Limitation Parameters             &
     &            (phytoplankton)[mol NH4 kg-1].'
             WRITE (out,100) k_no3_Lg(ng), 'k_no3_Lg',                  &
     &            'k_no3_Lg: Nutrient Limitation Parameters             &
     &            (phytoplankton)[mol NO3 kg-1].'
             WRITE (out,100) k_no3_Sm(ng), 'k_no3_Sm',                  &
     &            'k_no3_Sm: Nutrient Limitation Parameters             &
     &            (phytoplankton) [mol NO3 kg-1].'
             WRITE (out,100) k_no3_Di(ng), 'k_no3_Di',                  &
     &            'k_no3_Di: Nutrient Limitation Parameters             &
     &            (phytoplankton) [mol NO3 kg-1].'
             WRITE (out,100) k_po4_Di(ng), 'k_po4_Di',                  &
     &            'k_po4_Di: Nutrient Limitation Parameters             &
     &            (phytoplankton) [mol PO4 kg-1].'
             WRITE (out,100) k_po4_Lg(ng), 'k_po4_Lg',                  &
     &            'k_po4_Lg: Nutrient Limitation Parameters             &
     &            (phytoplankton) [mol PO4 kg-1].'
             WRITE (out,100) k_po4_Sm(ng), 'k_po4_Sm',                  &
     &            'k_po4_Sm: Nutrient Limitation Parameters             &
     &            (phytoplankton) [mol PO4 kg-1].'
             WRITE (out,100) k_sio4_Lg(ng), 'k_sio4_Lg',                &
     &            'k_sio4_Lg: Nutrient Limitation Parameters            &
     &            (phytoplankton) [mol SiO4 kg-1].'
             WRITE (out,100) k_fe_2_n_Di(ng), 'k_fe_2_n_Di',            &
     &            'k_fe_2_n_Di: Nutrient Limitation Parameters          &
     &            (phytoplankton) [mol Fe mol N-1].'
             WRITE (out,100) k_fe_2_n_Lg(ng), 'k_fe_2_n_Lg',                           &
     &            'k_fe_2_n_Lg: Nutrient Limitation Parameters (phytoplankton) [mol Fe mol N-1].'
             WRITE (out,100) k_fe_2_n_Sm(ng), 'k_fe_2_n_Sm',                           &
     &            'k_fe_2_n_Sm: Nutrient Limitation Parameters (phytoplankton) [mol Fe mol N-1].'
             WRITE (out,100) fe_2_n_max_Sm(ng), 'fe_2_n_max_Sm',                           &
     &            'fe_2_n_max_Sm: Nutrient Limitation Parameters (phytoplankton)[mol Fe mol N-1].'
             WRITE (out,100) fe_2_n_max_Lg(ng), 'fe_2_n_max_Lg',                           &
     &            'fe_2_n_max_Lg: Nutrient Limitation Parameters (phytoplankton)[mol Fe mol N-1].'
             WRITE (out,100) fe_2_n_max_Di(ng), 'fe_2_n_max_Di',                           &
     &            'fe_2_n_max_Di: Nutrient Limitation Parameters (phytoplankton)[mol Fe mol N-1].'
             WRITE (out,100) fe_2_n_upt_fac(ng), 'fe_2_n_upt_fac',                           &
     &            'fe_2_n_upt_fac: Nutrient Limitation Parameters (phytoplankton)[mol Fe mol N-1].'
             WRITE (out,100) alpha_Di(ng), 'alpha_Di',                           &
     &            'alpha_Di: Phytoplankton light limitation/growth rate[g C g Chl-1 m2 W-1 s-1].'
             WRITE (out,100) alpha_Lg(ng), 'alpha_Lg',                           &
     &            'alpha_Lg: Phytoplankton light limitation/growth rate[g C g Chl-1 m2 W-1 s-1].'
             WRITE (out,100) alpha_Sm(ng), 'alpha_Sm',                           &
     &            'alpha_Sm: Phytoplankton light limitation/growth rate[g C g Chl-1 m2 W-1 s-1].'
             WRITE (out,100) kappa_eppley(ng), 'kappa_eppley',                           &
     &            'kappa_eppley: Phytoplankton light limitation/growth rate[deg C-1].'
             WRITE (out,100) P_C_max_Di(ng), 'P_C_max_Di',                           &
     &            'P_C_max_Di: Phytoplankton light limitation/growth rate[s-1].'
             WRITE (out,100) P_C_max_Lg(ng), 'P_C_max_Lg',                           &
     &            'P_C_max_Lg: Phytoplankton light limitation/growth rate[s-1].'
             WRITE (out,100) P_C_max_Sm(ng), 'P_C_max_Sm',                           &
     &            'P_C_max_Sm: Phytoplankton light limitation/growth rate[s-1].'
             WRITE (out,100) thetamax_Di(ng), 'thetamax_Di',                           &
     &            'thetamax_Di: Phytoplankton light limitation/growth rate[g Chl g C-1].'
             WRITE (out,100) thetamax_Lg(ng), 'thetamax_Lg',                           &
     &            'thetamax_Lg: Phytoplankton light limitation/growth rate[g Chl g C-1].'
             WRITE (out,100) thetamax_Sm(ng), 'thetamax_Sm',                           &
     &            'thetamax_Sm: Phytoplankton light limitation/growth rate[g Chl g C-1].'
             WRITE (out,100) bresp_Di(ng), 'bresp_Di',                           &
     &            'bresp_Di: Phytoplankton light limitation/growth rate[sec-1].'
             WRITE (out,100) bresp_Lg(ng), 'bresp_Lg',                           &
     &            'bresp_Lg: Phytoplankton light limitation/growth rate[sec-1].'
             WRITE (out,100) bresp_Sm(ng), 'bresp_Sm',                           &
     &            'bresp_Sm: Phytoplankton light limitation/growth rate[sec-1].'
             WRITE (out,100) thetamin(ng), 'thetamin',                           &
     &            'thetamin: Phytoplankton light limitation/growth rate[g Chl g C-1].'
             WRITE (out,100) thetamin_nolim(ng), 'thetamin_nolim',                           &
     &            'thetamin_nolim: Phytoplankton light limitation/growth rate[g Chl g C-1].'
             WRITE (out,100) zpllgr(ng), 'zpllgr',                           &
     &            'zpllgr: Phytoplankton light limitation/growth rate[dimensionless].'
             WRITE (out,100) gamma_irr_mem(ng), 'gamma_irr_mem',                           &
     &            'gamma_irr_mem: Phytoplankton light limitation/growth rate[s-1].'
             WRITE (out,100) gamma_mu_mem(ng), 'gamma_mu_mem',                           &
     &            'gamma_mu_mem: Phytoplankton aggregation rate[s-1].'
             WRITE (out,100) k_n_inhib_Di(ng), 'k_n_inhib_Di',                           &
     &            'k_n_inhib_Di: Nitrogen fixation inhibition parameters[mol NO3 kg-1].'
             WRITE (out,100) o2_inhib_Di_pow(ng), 'o2_inhib_Di_pow',                           &
     &            'o2_inhib_Di_pow: Nitrogen fixation inhibition parameters[mol O2-1 m3].'
             WRITE (out,100) o2_inhib_Di_sat(ng), 'o2_inhib_Di_sat',                           &
     &            'o2_inhib_Di_sat: Nitrogen fixation inhibition parameters[mol O2 kg-1].'
!RD: same that NMKL
!             WRITE (out,100) p_2_n_static(ng), 'p_2_n_static',                           &
!     &            'p_2_n_static: Other stoichiometry[N/A].'
             WRITE (out,100) c_2_n(ng), 'c_2_n',                           &
     &            'c_2_n: Other stoichiometry[N/A].'
             WRITE (out,100) alk_2_n_denit(ng), 'alk_2_n_denit',                           &
     &            'alk_2_n_denit: Other stoichiometry[eq. alk mol NO3-1].'
             WRITE (out,100) p_2_n_static_Di(ng), 'p_2_n_static_Di',                           &
     &            'p_2_n_static_Di: Other stoichiometry[mol P mol N-1].'
             WRITE (out,100) p_2_n_static_Lg(ng), 'p_2_n_static_Lg',                           &
     &            'p_2_n_static_Lg: Other stoichiometry[mol P mol N-1].'
             WRITE (out,100) p_2_n_static_Sm(ng), 'p_2_n_static_Sm',                           &
     &            'p_2_n_static_Sm: Other stoichiometry[mol P mol N-1].'
             WRITE (out,100) si_2_n_static_Lg(ng), 'si_2_n_static_Lg',                           &
     &            'si_2_n_static_Lg: Other stoichiometry[mol Si mol N-1].'
             WRITE (out,100) si_2_n_max_Lg(ng), 'si_2_n_max_Lg',                           &
     &            'si_2_n_max_Lg: Other stoichiometry[mol Si mol N-1].'
             WRITE (out,100) ca_2_n_arag(ng), 'ca_2_n_arag',                           &
     &            'ca_2_n_arag: Other stoichiometry[mol Ca mol N-1].'
             WRITE (out,100) ca_2_n_calc(ng), 'ca_2_n_calc',                           &
     &            'ca_2_n_calc: Other stoichiometry[mol Ca mol N-1].'
             WRITE (out,100) caco3_sat_max(ng), 'caco3_sat_max',                           &
     &            'caco3_sat_max: Other stoichiometry[dimensionless].'
             WRITE (out,100) q_p_2_n_smz(ng), 'q_p_2_n_smz',                           &
     &            'q_p_2_n_smz: Zooplankton Stoichiometry - presently static[mol P mol N-1].'
             WRITE (out,100) q_p_2_n_mdz(ng), 'q_p_2_n_mdz',                           &
     &            'q_p_2_n_mdz: Zooplankton Stoichiometry - presently static[mol P mol N-1].'
             WRITE (out,100) q_p_2_n_lgz(ng), 'q_p_2_n_lgz',                           &
     &            'q_p_2_n_lgz: Zooplankton Stoichiometry - presently static[mol P mol N-1].'
             WRITE (out,100) q_p_2_n_bact(ng), 'q_p_2_n_bact',                           &
     &            'q_p_2_n_bact: Bacteria Stoichiometry - presently static[mol P mol N-1].'
             WRITE (out,100) agg_Sm(ng), 'agg_Sm',                           &
     &            'agg_Sm: Phytoplankton aggregation[s-1 (mole N kg)-1].'
             WRITE (out,100) agg_Di(ng), 'agg_Di',                           &
     &            'agg_Di: Phytoplankton aggregation[s-1 (mole N kg)-1].'
             WRITE (out,100) agg_Lg(ng), 'agg_Lg',                           &
     &            'agg_Lg: Phytoplankton aggregation[s-1 (mole N kg)-1].'
             WRITE (out,100) vir_Sm(ng), 'vir_Sm',                           &
     &            'vir_Sm: Phytoplankton and bacterial losses to viruses[s-1 (mole N kg)-1].'
             WRITE (out,100) vir_Di(ng), 'vir_Di',                           &
     &            'vir_Di: Phytoplankton and bacterial losses to viruses[s-1 (mole N kg)-1].'
             WRITE (out,100) vir_Lg(ng), 'vir_Lg',                           &
     &            'vir_Lg: Phytoplankton and bacterial losses to viruses[s-1 (mole N kg)-1].'
             WRITE (out,100) vir_Bact(ng), 'vir_Bact',                           &
     &            'vir_Bact: Phytoplankton and bacterial losses to viruses[s-1 (mole N kg)-1].'
             WRITE (out,100) ktemp_vir(ng), 'ktemp_vir',                           &
     &            'ktemp_vir: Phytoplankton and bacterial losses to viruses[C-1].'
             WRITE (out,100) exu_Sm(ng), 'exu_Sm',                           &
     &            'exu_Sm: Phytoplankton losses to exudation[dimensionless (fraction of NPP)].'
             WRITE (out,100) exu_Di(ng), 'exu_Di',                           &
     &            'exu_Di: Phytoplankton losses to exudation[dimensionless (fraction of NPP)].'
             WRITE (out,100) exu_Lg(ng), 'exu_Lg',                           &
     &            'exu_Lg: Phytoplankton losses to exudation[dimensionless (fraction of NPP)].'
             WRITE (out,100) imax_smz(ng), 'imax_smz',                           &
     &            'imax_smz: Zooplankton ingestion parameterization and temperature dependence[s-1].'
             WRITE (out,100) imax_mdz(ng), 'imax_mdz',                           &
     &            'imax_mdz: Zooplankton ingestion parameterization and temperature dependence[s-1].'
             WRITE (out,100) imax_lgz(ng), 'imax_lgz',                           &
     &            'imax_lgz: Zooplankton ingestion parameterization and temperature dependence[s-1].'
             WRITE (out,100) ki_smz(ng), 'ki_smz',                           &
     &            'ki_smz: Zooplankton ingestion parameterization and temperature dependence[moles N kg-1].'
             WRITE (out,100) ki_mdz(ng), 'ki_mdz',                           &
     &            'ki_mdz: Zooplankton ingestion parameterization and temperature dependence[moles N kg-1].'
             WRITE (out,100) ki_lgz(ng), 'ki_lgz',                           &
     &            'ki_lgz: Zooplankton ingestion parameterization and temperature dependence[moles N kg-1].'
             WRITE (out,100) ktemp_smz(ng), 'ktemp_smz',                           &
     &            'ktemp_smz: Zooplankton ingestion parameterization and temperature dependence[C-1].'
             WRITE (out,100) ktemp_mdz(ng), 'ktemp_mdz',                           &
     &            'ktemp_mdz: Zooplankton ingestion parameterization and temperature dependence[C-1].'
             WRITE (out,100) ktemp_lgz(ng), 'ktemp_lgz',                           &
     &            'ktemp_lgz: Zooplankton ingestion parameterization and temperature dependence[C-1].'
             WRITE (out,100) mu_max_bact(ng), 'mu_max_bact',                           &
     &            'mu_max_bact: Bacterial growth and uptake parameters[s-1].'
             WRITE (out,100) k_ldon_bact(ng), 'k_ldon_bact',                           &
     &            'k_ldon_bact: Bacterial growth and uptake parameters[mol ldon kg-1].'
             WRITE (out,100) ktemp_bact(ng), 'ktemp_bact',                           &
     &            'ktemp_bact: Bacterial growth and uptake parameters[C-1].'
             WRITE (out,100) nswitch_smz(ng), 'nswitch_smz',                           &
     &            'nswitch_smz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) nswitch_mdz(ng), 'nswitch_mdz',                           &
     &            'nswitch_mdz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) nswitch_lgz(ng), 'nswitch_lgz',                           &
     &            'nswitch_lgz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) mswitch_smz(ng), 'mswitch_smz',                           &
     &            'mswitch_smz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) mswitch_mdz(ng), 'mswitch_mdz',                           &
     &            'mswitch_mdz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) mswitch_lgz(ng), 'mswitch_lgz',                           &
     &            'mswitch_lgz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) smz_ipa_smp(ng), 'smz_ipa_smp',                           &
     &            'smz_ipa_smp: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) smz_ipa_lgp(ng), 'smz_ipa_lgp',                           &
     &            'smz_ipa_lgp: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) smz_ipa_diaz(ng), 'smz_ipa_diaz',                           &
     &            'smz_ipa_diaz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) smz_ipa_smz(ng), 'smz_ipa_smz',                           &
     &            'smz_ipa_smz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) smz_ipa_mdz(ng), 'smz_ipa_mdz',                           &
     &            'smz_ipa_mdz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) smz_ipa_lgz(ng), 'smz_ipa_lgz',                           &
     &            'smz_ipa_lgz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) smz_ipa_bact(ng), 'smz_ipa_bact',                           &
     &            'smz_ipa_bact: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) smz_ipa_det(ng), 'smz_ipa_det',                           &
     &            'smz_ipa_det: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) mdz_ipa_smp(ng), 'mdz_ipa_smp',                           &
     &            'mdz_ipa_smp: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) mdz_ipa_lgp(ng), 'mdz_ipa_lgp',                           &
     &            'mdz_ipa_lgp: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) mdz_ipa_diaz(ng), 'mdz_ipa_diaz',                           &
     &            'mdz_ipa_diaz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) mdz_ipa_smz(ng), 'mdz_ipa_smz',                           &
     &            'mdz_ipa_smz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) mdz_ipa_mdz(ng), 'mdz_ipa_mdz',                           &
     &            'mdz_ipa_mdz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) mdz_ipa_lgz(ng), 'mdz_ipa_lgz',                           &
     &            'mdz_ipa_lgz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) mdz_ipa_bact(ng), 'mdz_ipa_bact',                           &
     &            'mdz_ipa_bact: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) mdz_ipa_det(ng), 'mdz_ipa_det',                           &
     &            'mdz_ipa_det: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) lgz_ipa_smp(ng), 'lgz_ipa_smp',                           &
     &            'lgz_ipa_smp: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) lgz_ipa_lgp(ng), 'lgz_ipa_lgp',                           &
     &            'lgz_ipa_lgp: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) lgz_ipa_diaz(ng), 'lgz_ipa_diaz',                           &
     &            'lgz_ipa_diaz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) lgz_ipa_smz(ng), 'lgz_ipa_smz',                           &
     &            'lgz_ipa_smz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) lgz_ipa_mdz(ng), 'lgz_ipa_mdz',                           &
     &            'lgz_ipa_mdz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) lgz_ipa_lgz(ng), 'lgz_ipa_lgz',                           &
     &            'lgz_ipa_lgz: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) lgz_ipa_bact(ng), 'lgz_ipa_bact',                           &
     &            'lgz_ipa_bact: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) lgz_ipa_det(ng), 'lgz_ipa_det',                           &
     &            'lgz_ipa_det: Zooplankton switching and prey preference parameters[dimensionless].'
             WRITE (out,100) gge_max_smz(ng), 'gge_max_smz',                           &
     &            'gge_max_smz: Zooplankton bioenergetics[dimensionless].'
             WRITE (out,100) gge_max_mdz(ng), 'gge_max_mdz',                           &
     &            'gge_max_mdz: Zooplankton bioenergetics[dimensionless].'
             WRITE (out,100) gge_max_lgz(ng), 'gge_max_lgz',                           &
     &            'gge_max_lgz: Zooplankton bioenergetics[dimensionless].'
             WRITE (out,100) bresp_smz(ng), 'bresp_smz',                           &
     &            'bresp_smz: Zooplankton bioenergetics[s-1].'
             WRITE (out,100) bresp_mdz(ng), 'bresp_mdz',                           &
     &            'bresp_mdz: Zooplankton bioenergetics[s-1].'
             WRITE (out,100) bresp_lgz(ng), 'bresp_lgz',                           &
     &            'bresp_lgz: Zooplankton bioenergetics[s-1].'
             WRITE (out,100) gge_max_bact(ng), 'gge_max_bact',                           &
     &            'gge_max_bact: Bacterial bioenergetics[dimensionless].'
             WRITE (out,100) bresp_bact(ng), 'bresp_bact',                           &
     &            'bresp_bact: Bacterial bioenergetics[s-1].'
             WRITE (out,100) phi_det_smz(ng), 'phi_det_smz',                           &
     &            'phi_det_smz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_det_mdz(ng), 'phi_det_mdz',                           &
     &            'phi_det_mdz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_det_lgz(ng), 'phi_det_lgz',                           &
     &            'phi_det_lgz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_ldon_smz(ng), 'phi_ldon_smz',                           &
     &            'phi_ldon_smz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_ldon_mdz(ng), 'phi_ldon_mdz',                           &
     &            'phi_ldon_mdz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_ldon_lgz(ng), 'phi_ldon_lgz',                           &
     &            'phi_ldon_lgz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_ldop_smz(ng), 'phi_ldop_smz',                           &
     &            'phi_ldop_smz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_ldop_mdz(ng), 'phi_ldop_mdz',                           &
     &            'phi_ldop_mdz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_ldop_lgz(ng), 'phi_ldop_lgz',                           &
     &            'phi_ldop_lgz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_srdon_smz(ng), 'phi_srdon_smz',                           &
     &            'phi_srdon_smz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_srdon_mdz(ng), 'phi_srdon_mdz',                           &
     &            'phi_srdon_mdz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_srdon_lgz(ng), 'phi_srdon_lgz',                           &
     &            'phi_srdon_lgz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_srdop_smz(ng), 'phi_srdop_smz',                           &
     &            'phi_srdop_smz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_srdop_mdz(ng), 'phi_srdop_mdz',                           &
     &            'phi_srdop_mdz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_srdop_lgz(ng), 'phi_srdop_lgz',                           &
     &            'phi_srdop_lgz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_sldon_smz(ng), 'phi_sldon_smz',                           &
     &            'phi_sldon_smz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_sldon_mdz(ng), 'phi_sldon_mdz',                           &
     &            'phi_sldon_mdz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_sldon_lgz(ng), 'phi_sldon_lgz',                           &
     &            'phi_sldon_lgz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_sldop_smz(ng), 'phi_sldop_smz',                           &
     &            'phi_sldop_smz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_sldop_mdz(ng), 'phi_sldop_mdz',                           &
     &            'phi_sldop_mdz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_sldop_lgz(ng), 'phi_sldop_lgz',                           &
     &            'phi_sldop_lgz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_nh4_smz(ng), 'phi_nh4_smz',                           &
     &            'phi_nh4_smz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_nh4_mdz(ng), 'phi_nh4_mdz',                           &
     &            'phi_nh4_mdz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_nh4_lgz(ng), 'phi_nh4_lgz',                           &
     &            'phi_nh4_lgz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_po4_smz(ng), 'phi_po4_smz',                           &
     &            'phi_po4_smz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_po4_mdz(ng), 'phi_po4_mdz',                           &
     &            'phi_po4_mdz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_po4_lgz(ng), 'phi_po4_lgz',                           &
     &            'phi_po4_lgz: Partitioning of zooplankton ingestion to other compartments[dimensionless].'
             WRITE (out,100) phi_ldon_vir(ng), 'phi_ldon_vir',                           &
     &            'phi_ldon_vir: Partitioning of viral losses to various dissolved pools[dimensionless].'
             WRITE (out,100) phi_srdon_vir(ng), 'phi_srdon_vir',                           &
     &            'phi_srdon_vir: Partitioning of viral losses to various dissolved pools[dimensionless].'
             WRITE (out,100) phi_sldon_vir(ng), 'phi_sldon_vir',                           &
     &            'phi_sldon_vir: Partitioning of viral losses to various dissolved pools[dimensionless].'
             WRITE (out,100) phi_ldop_vir(ng), 'phi_ldop_vir',                           &
     &            'phi_ldop_vir: Partitioning of viral losses to various dissolved pools[dimensionless].'
             WRITE (out,100) phi_srdop_vir(ng), 'phi_srdop_vir',                           &
     &            'phi_srdop_vir: Partitioning of viral losses to various dissolved pools[dimensionless].'
             WRITE (out,100) phi_sldop_vir(ng), 'phi_sldop_vir',                           &
     &            'phi_sldop_vir: Partitioning of viral losses to various dissolved pools[dimensionless].'
             WRITE (out,100) imax_hp(ng), 'imax_hp',                           &
     &            'imax_hp: Parameters for unresolved higher predators[s-1].'
             WRITE (out,100) ki_hp(ng), 'ki_hp',                           &
     &            'ki_hp: Parameters for unresolved higher predators[mol N kg-1].'
             WRITE (out,100) coef_hp(ng), 'coef_hp',                           &
     &            'coef_hp: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) ktemp_hp(ng), 'ktemp_hp',                           &
     &            'ktemp_hp: Parameters for unresolved higher predators[C-1].'
             WRITE (out,100) nswitch_hp(ng), 'nswitch_hp',                           &
     &            'nswitch_hp: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) mswitch_hp(ng), 'mswitch_hp',                           &
     &            'mswitch_hp: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_ipa_smp(ng), 'hp_ipa_smp',                           &
     &            'hp_ipa_smp: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_ipa_lgp(ng), 'hp_ipa_lgp',                           &
     &            'hp_ipa_lgp: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_ipa_diaz(ng), 'hp_ipa_diaz',                           &
     &            'hp_ipa_diaz: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_ipa_smz(ng), 'hp_ipa_smz',                           &
     &            'hp_ipa_smz: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_ipa_mdz(ng), 'hp_ipa_mdz',                           &
     &            'hp_ipa_mdz: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_ipa_lgz(ng), 'hp_ipa_lgz',                           &
     &            'hp_ipa_lgz: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_ipa_bact(ng), 'hp_ipa_bact',                           &
     &            'hp_ipa_bact: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_ipa_det(ng), 'hp_ipa_det',                           &
     &            'hp_ipa_det: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_phi_det(ng), 'hp_phi_det',                           &
     &            'hp_phi_det: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_phi_ldon(ng), 'hp_phi_ldon',                           &
     &            'hp_phi_ldon: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_phi_ldop(ng), 'hp_phi_ldop',                           &
     &            'hp_phi_ldop: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_phi_srdon(ng), 'hp_phi_srdon',                           &
     &            'hp_phi_srdon: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_phi_srdop(ng), 'hp_phi_srdop',                           &
     &            'hp_phi_srdop: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_phi_sldon(ng), 'hp_phi_sldon',                           &
     &            'hp_phi_sldon: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_phi_sldop(ng), 'hp_phi_sldop',                           &
     &            'hp_phi_sldop: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_phi_nh4(ng), 'hp_phi_nh4',                           &
     &            'hp_phi_nh4: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) hp_phi_po4(ng), 'hp_phi_po4',                           &
     &            'hp_phi_po4: Parameters for unresolved higher predators[dimensionless].'
             WRITE (out,100) felig_bkg(ng), 'felig_bkg',                           &
     &            'felig_bkg: Iron chemistry[mol Fe kg-1].'
             WRITE (out,100) felig_2_don(ng), 'felig_2_don',                           &
     &            'felig_2_don: Iron chemistry[mol Fe mol N-1].'
             WRITE (out,100) fe_2_n_sed(ng), 'fe_2_n_sed',                           &
     &            'fe_2_n_sed: Iron chemistry[mol Fe mol N-1].'
             WRITE (out,100) fe_coast(ng), 'fe_coast',                           &
     &            'fe_coast: Iron chemistry[mol Fe m kg-1 s-1].'
             WRITE (out,100) alpha_fescav(ng), 'alpha_fescav',                           &
     &            'alpha_fescav: Iron chemistry[sec-1].'
             WRITE (out,100) beta_fescav(ng), 'beta_fescav',                           &
     &            'beta_fescav: Iron chemistry[mol N-1 sec-1].'
             WRITE (out,100) remin_eff_fedet(ng), 'remin_eff_fedet',                           &
     &            'remin_eff_fedet: Iron chemistry[unitless].'
             WRITE (out,100) ki_fescav(ng), 'ki_fescav',                           &
     &            'ki_fescav: Iron chemistry[watts m-2].'
             WRITE (out,100) io_fescav(ng), 'io_fescav',                           &
     &            'io_fescav: Iron chemistry[watts m-2].'
             WRITE (out,100) gamma_fescav(ng), 'gamma_fescav',                           &
     &            'gamma_fescav: Iron chemistry[watts m-2].'
             WRITE (out,100) kfe_eq_lig_ll(ng), 'kfe_eq_lig_ll',                           &
     &            'kfe_eq_lig_ll: Iron chemistry[mol lig-1 kg].'
             WRITE (out,100) kfe_eq_lig_hl(ng), 'kfe_eq_lig_hl',                           &
     &            'kfe_eq_lig_hl: Iron chemistry[mol lig-1 kg].'
             WRITE (out,100) k_o2(ng), 'k_o2',                           &
     &            'k_o2: Remineralization[mol O2 kg-1].'
             WRITE (out,100) o2_min(ng), 'o2_min',                           &
     &            'o2_min: Remineralization[mol O2 kg-1].'
             WRITE (out,100) rpcaco3(ng), 'rpcaco3',                           &
     &            'rpcaco3: Remineralization[mol N mol Ca-1].'
             WRITE (out,100) rplith(ng), 'rplith',                           &
     &            'rplith: Remineralization[mol N g lith-1].'
             WRITE (out,100) rpsio2(ng), 'rpsio2',                           &
     &            'rpsio2: Remineralization[mol N mol Si-1].'
             WRITE (out,100) gamma_ndet(ng), 'gamma_ndet',                           &
     &            'gamma_ndet: Remineralization[s-1].'
             WRITE (out,100) gamma_cadet_arag(ng), 'gamma_cadet_arag',                           &
     &            'gamma_cadet_arag: Remineralization[s-1].'
             WRITE (out,100) gamma_cadet_calc(ng), 'gamma_cadet_calc',                           &
     &            'gamma_cadet_calc: Remineralization[s-1].'
             WRITE (out,100) gamma_sidet(ng), 'gamma_sidet',                           &
     &            'gamma_sidet: Remineralization[s-1].'
             WRITE (out,100) phi_lith(ng), 'phi_lith',                           &
     &            'phi_lith: Remineralization[kg mol-1].'
             WRITE (out,100) k_lith(ng), 'k_lith',                           &
     &            'k_lith: Remineralization[s-1].'
             WRITE (out,100) z_sed(ng), 'z_sed',                           &
     &            'z_sed: Remineralization[m].'
             WRITE (out,100) k_no3_denit(ng), 'k_no3_denit',                           &
     &            'k_no3_denit: Remineralization[mol NO3 kg-1].'
             WRITE (out,100) gamma_srdon(ng), 'gamma_srdon',                           &
     &            'gamma_srdon: Dissolved Organic Material[s-1].'
             WRITE (out,100) gamma_srdop(ng), 'gamma_srdop',                           &
     &            'gamma_srdop: Dissolved Organic Material[s-1].'
             WRITE (out,100) gamma_sldon(ng), 'gamma_sldon',                           &
     &            'gamma_sldon: Dissolved Organic Material[s-1].'
             WRITE (out,100) gamma_sldop(ng), 'gamma_sldop',                           &
     &            'gamma_sldop: Dissolved Organic Material[s-1].'
             WRITE (out,100) gamma_nitrif(ng), 'gamma_nitrif',                           &
     &            'gamma_nitrif: Nitrification[s-1].'
             WRITE (out,100) irr_inhibit(ng), 'irr_inhibit',                           &
     &            'irr_inhibit: Nitrification[m2 W-1].'
!             WRITE (out,100) tracer_debug(ng), 'tracer_debug',                           &
!     &            'tracer_debug: MISC[N/A].'
             WRITE (out,100) htotal_in(ng), 'htotal_in',                           &
     &            'htotal_in: ?[N/A].'
             WRITE (out,100) wsink(ng), 'wsink',                           &
     &            'wsink: Sinking velocity of detritus[m s-1].'
#ifdef COASTDIAT
             WRITE (out,100) k_fed_Md(ng), 'k_fed_Md',                  &
     &            'k_fed_Md: Nutrient Limitation Parameters             &
     &            (phytoplankton)[mol Fed kg-1].'
             WRITE (out,100) k_nh4_Md(ng), 'k_nh4_Md',                  &
     &            'k_nh4_Md: Nutrient Limitation Parameters             &
     &            (phytoplankton)[mol NH4 kg-1].'
             WRITE (out,100) k_no3_Md(ng), 'k_no3_Md',                  &
     &            'k_no3_Md: Nutrient Limitation Parameters             &
     &            (phytoplankton) [mol NO3 kg-1].'
             WRITE (out,100) k_po4_Md(ng), 'k_po4_Md',                  &
     &            'k_po4_Md: Nutrient Limitation Parameters             &
     &            (phytoplankton) [mol PO4 kg-1].'
             WRITE (out,100) k_sio4_Md(ng), 'k_sio4_Md',                &
     &            'k_sio4_Md: Nutrient Limitation Parameters            &
     &            (phytoplankton) [mol SiO4 kg-1].'
             WRITE (out,100) k_fe_2_n_Md(ng), 'k_fe_2_n_Md',            &
     &            'k_fe_2_n_Md: Nutrient Limitation Parameters          &
     &            (phytoplankton) [mol Fe mol N-1].'
#endif
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Rescale biological tracer parameters
!-----------------------------------------------------------------------
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          tnu4(i,ng)=SQRT(ABS(tnu4(i,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' GFDL Cobalt Model Parameters, Grid: ',i2.2,         &
     &        /,  ' =================================',/)
  50  FORMAT (1x,i10,2x,a,t30,a)
  60  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
  70  FORMAT (f11.3,2x,a,t30,a)
  80  FORMAT (f11.3,2x,a,t30,a,/,t32,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 100  FORMAT (1p,e11.4,2x,a,t30,a)
 110  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 120  FORMAT (/,' read_BioPar - variable info not yet loaded, ',        &
     &        a,i2.2,a)
 130  FORMAT (/,' read_BioPar - variable info not yet loaded, ',a)
 140  FORMAT (10x,l1,2x,a,t30,a,1x,a)

      RETURN
      END SUBROUTINE read_BioPar
