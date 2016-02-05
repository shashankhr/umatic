#ifdef OLDSOL
      /* Dliq */
    } else if (strcasecmp (token, "DliqAlloy") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->Dliq_alloy = atof (token);
      /* Dsol */
    } else if (strcasecmp (token, "DsolAlloy") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->Dsol_alloy = atof (token);
      /* part_coef alloy */
    } else if (strcasecmp (token, "part_coef_Alloy") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->part_coef_alloy = atof (token);
      /* phase_diagram slope */
    } else if (strcasecmp (token, "m_Alloy") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->m_alloy = atof (token);
      /*superheat for melt back xly */
    } else if (strcasecmp (token, "Dt_Melt") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->dt_melt = atof (token);
    } else if (strcasecmp (token, "InitialConcAlloy") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->Cinit_alloy = atof (token);
      /* phase diagram ref. temperature (pure A) */
    } else if (strcasecmp (token, "t_pure") == 0) {
      if ((token = strtok (NULL, sep)) != NULL) {
        tp_flag = 1;
        mp->tp = atof (token);
      }
>>> edit marker <<< Take polyphase out of oldsol section
    } else if (strcasecmp (token, "t_purePoly") == 0) {

      tp_flag = 1;

      for (ieq = 0; ieq < ieq_tot; ieq++) {
        if ((token = strtok (NULL, sep)) != NULL) {
          mp->tp_poly[ieq] = atof (token);
        }
      }

      /* phase diagram eutectic temperature  */
    } else if (strcasecmp (token, "T_eut") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->T_eut = atof (token);
      /* phase diagram eutectic composition  */
    } else if (strcasecmp (token, "c_eut") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        ce_flag = 1;
      mp->c_eut = atof (token);
      /* Surf_Tens_Alloy *//*by Wei WANG 11-07-02 */
    } else if (strcasecmp (token, "Surf_Tens_Alloy") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->surf_tens_alloy = atof (token);

      /* Mould_Source_Alloy  type of source */
    } else if (strcasecmp (token, "MouldSourceAlloy") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->alloyprops.mould_src = (SrcFn_T) atoi (token);
      /* Mould_Source_Value_Alloy  strength of source */
    } else if (strcasecmp (token, "MouldSourceValueAlloy") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->alloyprops.mould_source_value = atof (token);

    } else if (strcasecmp (token, "DoLiqAlloy") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->alloyprops.DoLiq = atof (token);
    } else if (strcasecmp (token, "QaLiqAlloy") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->alloyprops.QaLiq = atof (token);
    } else if (strcasecmp (token, "DoSolAlloy") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->alloyprops.DoSol = atof (token);
    } else if (strcasecmp (token, "QaSolAlloy") == 0) {
      if ((token = strtok (NULL, sep)) != NULL)
        mp->alloyprops.QaSol = atof (token);
#endif /*OLDSOL*/
