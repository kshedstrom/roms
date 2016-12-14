!
!  Write out GFDL Cobalt ecosystem parameters.
!
!       CALL netcdf_put_ivar (ng, model, ncname, 'init',                 &
!                             init(ng), (/0/), (/0/),                    &
!      &                      ncid = ncid)
!       IF (exit_flag.ne.NoError) RETURN

      CALL netcdf_put_ivar (ng, model, ncname, 'BioIter',               &
     &                      BioIter(ng), (/0/), (/0/),                  &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'htotal_scale_lo',      &
      &                      htotal_scale_lo(ng), (/0/), (/0/),         &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'htotal_scale_hi',      &
      &                      htotal_scale_hi(ng), (/0/), (/0/),         &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'RHO_0',                &
      &                      RHO_0(ng), (/0/), (/0/),                   &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_ivar (ng, model, ncname, 'NKML',                 &
      &                      NKML(ng), (/0/), (/0/),                    &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a_0',                  &
      &                      a_0(ng), (/0/), (/0/),                     &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a_1',                  &
      &                      a_1(ng), (/0/), (/0/),                     &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a_2',                  &
      &                      a_2(ng), (/0/), (/0/),                     &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a_3',                  &
      &                      a_3(ng), (/0/), (/0/),                     &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a_4',                  &
      &                      a_4(ng), (/0/), (/0/),                     &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a_5',                  &
      &                      a_5(ng), (/0/), (/0/),                     &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'b_0',                  &
      &                      b_0(ng), (/0/), (/0/),                     &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'b_1',                  &
      &                      b_1(ng), (/0/), (/0/),                     &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'b_2',                  &
      &                      b_2(ng), (/0/), (/0/),                     &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'b_3',                  &
      &                      b_3(ng), (/0/), (/0/),                     &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'c_0',                  &
      &                      c_0(ng), (/0/), (/0/),                     &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a1_co2',               &
      &                      a1_co2(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a2_co2',               &
      &                      a2_co2(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a3_co2',               &
      &                      a3_co2(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a4_co2',               &
      &                      a4_co2(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a1_o2',                &
      &                      a1_o2(ng), (/0/), (/0/),                   &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a2_o2',                &
      &                      a2_o2(ng), (/0/), (/0/),                   &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a3_o2',                &
      &                      a3_o2(ng), (/0/), (/0/),                   &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'a4_o2',                &
      &                      a4_o2(ng), (/0/), (/0/),                   &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mass_2_n',             &
      &                      mass_2_n(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'n_2_n_denit',          &
      &                      n_2_n_denit(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'o2_2_c',               &
      &                      o2_2_c(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'o2_2_nfix',            &
      &                      o2_2_nfix(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'o2_2_nh4',             &
      &                      o2_2_nh4(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'o2_2_nitrif',          &
      &                      o2_2_nitrif(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'o2_2_no3',             &
      &                      o2_2_no3(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_fed_Di',             &
      &                      k_fed_Di(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_fed_Lg',             &
      &                      k_fed_Lg(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_fed_Sm',             &
      &                      k_fed_Sm(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_nh4_Lg',             &
      &                      k_nh4_Lg(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_nh4_Sm',             &
      &                      k_nh4_Sm(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_nh4_Di',             &
      &                      k_nh4_Di(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_no3_Lg',             &
      &                      k_no3_Lg(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_no3_Sm',             &
      &                      k_no3_Sm(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_no3_Di',             &
      &                      k_no3_Di(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_po4_Di',             &
      &                      k_po4_Di(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_po4_Lg',             &
      &                      k_po4_Lg(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_po4_Sm',             &
      &                      k_po4_Sm(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_sio4_Lg',            &
      &                      k_sio4_Lg(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_fe_2_n_Di',          &
      &                      k_fe_2_n_Di(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_fe_2_n_Lg',          &
      &                      k_fe_2_n_Lg(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_fe_2_n_Sm',          &
      &                      k_fe_2_n_Sm(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'fe_2_n_max_Sm',        &
      &                      fe_2_n_max_Sm(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'fe_2_n_max_Lg',        &
      &                      fe_2_n_max_Lg(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'fe_2_n_max_Di',        &
      &                      fe_2_n_max_Di(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'fe_2_n_upt_fac',       &
      &                      fe_2_n_upt_fac(ng), (/0/), (/0/),          &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'alpha_Di',             &
      &                      alpha_Di(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'alpha_Lg',             &
      &                      alpha_Lg(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'alpha_Sm',             &
      &                      alpha_Sm(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'kappa_eppley',         &
      &                      kappa_eppley(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'P_C_max_Di',           &
      &                      P_C_max_Di(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'P_C_max_Lg',           &
      &                      P_C_max_Lg(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'P_C_max_Sm',           &
      &                      P_C_max_Sm(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'thetamax_Di',          &
      &                      thetamax_Di(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'thetamax_Lg',          &
      &                      thetamax_Lg(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'thetamax_Sm',          &
      &                      thetamax_Sm(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'bresp_Di',             &
      &                      bresp_Di(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'bresp_Lg',             &
      &                      bresp_Lg(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'bresp_Sm',             &
      &                      bresp_Sm(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'thetamin',             &
      &                      thetamin(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'thetamin_nolim',       &
      &                      thetamin_nolim(ng), (/0/), (/0/),          &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'zpllgr',                 &
      &                      zpllgr(ng), (/0/), (/0/),                    &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gamma_irr_mem',        &
      &                      gamma_irr_mem(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gamma_mu_mem',        &
      &                      gamma_mu_mem(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_n_inhib_Di',         &
      &                      k_n_inhib_Di(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'o2_inhib_Di_pow',      &
      &                      o2_inhib_Di_pow(ng), (/0/), (/0/),         &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'o2_inhib_Di_sat',      &
      &                      o2_inhib_Di_sat(ng), (/0/), (/0/),         &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_ivar (ng, model, ncname, 'p_2_n_static',         &
      &                      p_2_n_static(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'c_2_n',                &
      &                      c_2_n(ng), (/0/), (/0/),                   &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'alk_2_n_denit',        &
      &                      alk_2_n_denit(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'p_2_n_static_Di',      &
      &                      p_2_n_static_Di(ng), (/0/), (/0/),         &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'p_2_n_static_Lg',      &
      &                      p_2_n_static_Lg(ng), (/0/), (/0/),         &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'p_2_n_static_Sm',      &
      &                      p_2_n_static_Sm(ng), (/0/), (/0/),         &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'si_2_n_static_Lg',     &
      &                      si_2_n_static_Lg(ng), (/0/), (/0/),        &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'si_2_n_max_Lg',        &
      &                      si_2_n_max_Lg(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'ca_2_n_arag',          &
      &                      ca_2_n_arag(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'ca_2_n_calc',          &
      &                      ca_2_n_calc(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'caco3_sat_max',        &
      &                      caco3_sat_max(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'q_p_2_n_smz',          &
      &                      q_p_2_n_smz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'q_p_2_n_mdz',          &
      &                      q_p_2_n_mdz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'q_p_2_n_lgz',          &
      &                      q_p_2_n_lgz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'q_p_2_n_bact',         &
      &                      q_p_2_n_bact(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'agg_Sm',               &
      &                      agg_Sm(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'agg_Di',               &
      &                      agg_Di(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'agg_Lg',               &
      &                      agg_Lg(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'vir_Sm',               &
      &                      vir_Sm(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'vir_Di',               &
      &                      vir_Di(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'vir_Lg',               &
      &                      vir_Lg(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'vir_Bact',             &
      &                      vir_Bact(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'ktemp_vir',            &
      &                      ktemp_vir(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'exu_Sm',               &
      &                      exu_Sm(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'exu_Di',               &
      &                      exu_Di(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'exu_Lg',               &
      &                      exu_Lg(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'imax_smz',             &
      &                      imax_smz(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'imax_mdz',             &
      &                      imax_mdz(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'imax_lgz',             &
      &                      imax_lgz(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'ki_smz',               &
      &                      ki_smz(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'ki_mdz',               &
      &                      ki_mdz(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'ki_lgz',               &
      &                      ki_lgz(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'ktemp_smz',            &
      &                      ktemp_smz(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'ktemp_mdz',            &
      &                      ktemp_mdz(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'ktemp_lgz',            &
      &                      ktemp_lgz(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mu_max_bact',          &
      &                      mu_max_bact(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_ldon_bact',          &
      &                      k_ldon_bact(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'ktemp_bact',           &
      &                      ktemp_bact(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'nswitch_smz',          &
      &                      nswitch_smz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'nswitch_mdz',          &
      &                      nswitch_mdz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'nswitch_lgz',          &
      &                      nswitch_lgz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mswitch_smz',          &
      &                      mswitch_smz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mswitch_mdz',          &
      &                      mswitch_mdz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mswitch_lgz',          &
      &                      mswitch_lgz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'smz_ipa_smp',          &
      &                      smz_ipa_smp(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'smz_ipa_lgp',          &
      &                      smz_ipa_lgp(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'smz_ipa_diaz',         &
      &                      smz_ipa_diaz(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'smz_ipa_smz',          &
      &                      smz_ipa_smz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'smz_ipa_mdz',          &
      &                      smz_ipa_mdz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'smz_ipa_lgz',          &
      &                      smz_ipa_lgz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'smz_ipa_bact',         &
      &                      smz_ipa_bact(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'smz_ipa_det',          &
      &                      smz_ipa_det(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mdz_ipa_smp',          &
      &                      mdz_ipa_smp(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mdz_ipa_lgp',          &
      &                      mdz_ipa_lgp(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mdz_ipa_diaz',         &
      &                      mdz_ipa_diaz(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mdz_ipa_smz',          &
      &                      mdz_ipa_smz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mdz_ipa_mdz',          &
      &                      mdz_ipa_mdz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mdz_ipa_lgz',          &
      &                      mdz_ipa_lgz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mdz_ipa_bact',         &
      &                      mdz_ipa_bact(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mdz_ipa_det',          &
      &                      mdz_ipa_det(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'lgz_ipa_smp',          &
      &                      lgz_ipa_smp(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'lgz_ipa_lgp',          &
      &                      lgz_ipa_lgp(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'lgz_ipa_diaz',         &
      &                      lgz_ipa_diaz(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'lgz_ipa_smz',          &
      &                      lgz_ipa_smz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'lgz_ipa_mdz',          &
      &                      lgz_ipa_mdz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'lgz_ipa_lgz',          &
      &                      lgz_ipa_lgz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'lgz_ipa_bact',         &
      &                      lgz_ipa_bact(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'lgz_ipa_det',          &
      &                      lgz_ipa_det(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gge_max_smz',          &
      &                      gge_max_smz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gge_max_mdz',          &
      &                      gge_max_mdz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gge_max_lgz',          &
      &                      gge_max_lgz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'bresp_smz',            &
      &                      bresp_smz(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'bresp_mdz',            &
      &                      bresp_mdz(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'bresp_lgz',            &
      &                      bresp_lgz(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gge_max_bact',         &
      &                      gge_max_bact(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'bresp_bact',           &
      &                      bresp_bact(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_det_smz',          &
      &                      phi_det_smz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_det_mdz',          &
      &                      phi_det_mdz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_det_lgz',          &
      &                      phi_det_lgz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_ldon_smz',         &
      &                      phi_ldon_smz(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_ldon_mdz',         &
      &                      phi_ldon_mdz(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_ldon_lgz',         &
      &                      phi_ldon_lgz(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_ldop_smz',         &
      &                      phi_ldop_smz(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_ldop_mdz',         &
      &                      phi_ldop_mdz(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_ldop_lgz',         &
      &                      phi_ldop_lgz(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_srdon_smz',        &
      &                      phi_srdon_smz(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_srdon_mdz',        &
      &                      phi_srdon_mdz(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_srdon_lgz',        &
      &                      phi_srdon_lgz(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_srdop_smz',        &
      &                      phi_srdop_smz(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_srdop_mdz',        &
      &                      phi_srdop_mdz(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_srdop_lgz',        &
      &                      phi_srdop_lgz(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_sldon_smz',        &
      &                      phi_sldon_smz(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_sldon_mdz',        &
      &                      phi_sldon_mdz(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_sldon_lgz',        &
      &                      phi_sldon_lgz(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_sldop_smz',        &
      &                      phi_sldop_smz(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_sldop_mdz',        &
      &                      phi_sldop_mdz(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_sldop_lgz',        &
      &                      phi_sldop_lgz(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_nh4_smz',          &
      &                      phi_nh4_smz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_nh4_mdz',          &
      &                      phi_nh4_mdz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_nh4_lgz',          &
      &                      phi_nh4_lgz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_po4_smz',          &
      &                      phi_po4_smz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_po4_mdz',          &
      &                      phi_po4_mdz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_po4_lgz',          &
      &                      phi_po4_lgz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_ldon_vir',         &
      &                      phi_ldon_vir(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_srdon_vir',        &
      &                      phi_srdon_vir(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_sldon_vir',        &
      &                      phi_sldon_vir(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_ldop_vir',         &
      &                      phi_ldop_vir(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_srdop_vir',        &
      &                      phi_srdop_vir(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_sldop_vir',        &
      &                      phi_sldop_vir(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'imax_hp',              &
      &                      imax_hp(ng), (/0/), (/0/),                 &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'ki_hp',                &
      &                      ki_hp(ng), (/0/), (/0/),                   &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'coef_hp',              &
      &                      coef_hp(ng), (/0/), (/0/),                 &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'ktemp_hp',             &
      &                      ktemp_hp(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'nswitch_hp',           &
      &                      nswitch_hp(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mswitch_hp',           &
      &                      mswitch_hp(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_ipa_smp',           &
      &                      hp_ipa_smp(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_ipa_lgp',           &
      &                      hp_ipa_lgp(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_ipa_diaz',          &
      &                      hp_ipa_diaz(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_ipa_smz',           &
      &                      hp_ipa_smz(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_ipa_mdz',           &
      &                      hp_ipa_mdz(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_ipa_lgz',           &
      &                      hp_ipa_lgz(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_ipa_bact',          &
      &                      hp_ipa_bact(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_ipa_det',           &
      &                      hp_ipa_det(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_phi_det',           &
      &                      hp_phi_det(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_phi_ldon',          &
      &                      hp_phi_ldon(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_phi_ldop',          &
      &                      hp_phi_ldop(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_phi_srdon',         &
      &                      hp_phi_srdon(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_phi_srdop',         &
      &                      hp_phi_srdop(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_phi_sldon',         &
      &                      hp_phi_sldon(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_phi_sldop',         &
      &                      hp_phi_sldop(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_phi_nh4',           &
      &                      hp_phi_nh4(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'hp_phi_po4',           &
      &                      hp_phi_po4(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'felig_bkg',            &
      &                      felig_bkg(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'felig_2_don',          &
      &                      felig_2_don(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'fe_2_n_sed',           &
      &                      fe_2_n_sed(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'fe_coast',             &
      &                      fe_coast(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'alpha_fescav',         &
      &                      alpha_fescav(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'beta_fescav',          &
      &                      beta_fescav(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'remin_eff_fedet',      &
      &                      remin_eff_fedet(ng), (/0/), (/0/),         &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'ki_fescav',            &
      &                      ki_fescav(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'io_fescav',            &
      &                      io_fescav(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gamma_fescav',         &
      &                      gamma_fescav(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'kfe_eq_lig_ll',        &
      &                      kfe_eq_lig_ll(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'kfe_eq_lig_hl',        &
      &                      kfe_eq_lig_hl(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_o2',                 &
      &                      k_o2(ng), (/0/), (/0/),                    &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'o2_min',               &
      &                      o2_min(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'rpcaco3',              &
      &                      rpcaco3(ng), (/0/), (/0/),                 &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'rplith',               &
      &                      rplith(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'rpsio2',               &
      &                      rpsio2(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gamma_ndet',           &
      &                      gamma_ndet(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gamma_cadet_arag',     &
      &                      gamma_cadet_arag(ng), (/0/), (/0/),        &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gamma_cadet_calc',     &
      &                      gamma_cadet_calc(ng), (/0/), (/0/),        &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gamma_sidet',          &
      &                      gamma_sidet(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'phi_lith',             &
      &                      phi_lith(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_lith',               &
      &                      k_lith(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'z_sed',                &
      &                      z_sed(ng), (/0/), (/0/),                   &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_no3_denit',          &
      &                      k_no3_denit(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gamma_srdon',          &
      &                      gamma_srdon(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gamma_srdop',          &
      &                      gamma_srdop(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gamma_sldon',          &
      &                      gamma_sldon(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gamma_sldop',          &
      &                      gamma_sldop(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'gamma_nitrif',         &
      &                      gamma_nitrif(ng), (/0/), (/0/),            &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'irr_inhibit',          &
      &                      irr_inhibit(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

     !  CALL netcdf_put_fvar (ng, model, ncname, 'tracer_debug',         &
     ! &                      tracer_debug(ng), (/0/), (/0/),            &
     ! &                      ncid = ncid)
     !  IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'htotal_in',            &
      &                      htotal_in(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'wsink',                &
      &                      wsink(ng), (/0/), (/0/),                   &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

#ifdef COASTDIAT
       CALL netcdf_put_fvar (ng, model, ncname, 'k_fed_Md',             &
      &                      k_fed_Md(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_nh4_Md',             &
      &                      k_nh4_Md(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_no3_Md',             &
      &                      k_no3_Md(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_po4_Md',             &
      &                      k_po4_Md(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_sio4_Md',            &
      &                      k_sio4_Md(ng), (/0/), (/0/),               &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'k_fe_2_n_Md',          &
      &                      k_fe_2_n_Md(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'fe_2_n_max_Md',        &
      &                      fe_2_n_max_Md(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'alpha_Md',             &
      &                      alpha_Md(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'P_C_max_Md',           &
      &                      P_C_max_Md(ng), (/0/), (/0/),              &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'thetamax_Md',          &
      &                      thetamax_Md(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'bresp_Md',             &
      &                      bresp_Md(ng), (/0/), (/0/),                &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'p_2_n_static_Md',      &
      &                      p_2_n_static_Md(ng), (/0/), (/0/),         &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'si_2_n_static_Md',     &
      &                      si_2_n_static_Md(ng), (/0/), (/0/),        &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'si_2_n_max_Md',        &
      &                      si_2_n_max_Md(ng), (/0/), (/0/),           &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'agg_Md',               &
      &                      agg_Md(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'vir_Md',               &
      &                      vir_Md(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'exu_Md',               &
      &                      exu_Md(ng), (/0/), (/0/),                  &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'smz_ipa_mdp',          &
      &                      smz_ipa_mdp(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'mdz_ipa_mdp',          &
      &                      mdz_ipa_mdp(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN

       CALL netcdf_put_fvar (ng, model, ncname, 'lgz_ipa_mdp',          &
      &                      lgz_ipa_mdp(ng), (/0/), (/0/),             &
      &                      ncid = ncid)
       IF (exit_flag.ne.NoError) RETURN
#endif
