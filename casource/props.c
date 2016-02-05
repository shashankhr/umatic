
/****************************************************************/
/*      Copyright (c) 1993 Peter D Lee                          */
/*      Copyright (c) 1998 Dept. of Materials, ICSTM            */
/*      All Rights Reserved                                     */
/*      The copyright notice above does not evidence any        */
/*      actual or intended publication of such source code,     */
/*      and is an unpublished work by Dept. of Materials, ICSTM.*/
/*      continuing D Phil work from University of Oxford        */
/****************************************************************/
/* This code is part of the umats routines developed at in the  */
/* Materials Processing Group, Dept. of Materials, ICSTM.       */
/*      email p.d.lee or r.atwood @imperial.ac.uk for details   */
/****************************************************************/

/********************************************************************************/
/*  This version is distributed under a BSD style public license, as follows:   */
/*                                                                              */
/*  Copyright (c) 2007, Dept. of Materials, Imperial College London             */
/*  All rights reserved.                                                        */
/*  Redistribution and use in source and binary forms, with or without          */
/*  modification, are permitted provided that the following conditions          */
/*  are met:                                                                    */
/*                                                                              */
/*  * Redistributions of source code must retain the above copyright            */
/*  notice, this list of conditions and the following disclaimer.               */
/*                                                                              */
/*  * Redistributions in binary form must reproduce the above                   */
/*  copyright notice, this list of conditions and the following                 */
/*  disclaimer in the documentation and/or other materials provided             */
/*  with the distribution.                                                      */
/*                                                                              */
/*  * Neither the name of the Dept. of Materials, Imperial College London, nor  */
/*  the names of its contributors may be used to endorse or promote products    */
/*  derived from this software without specific prior written permission.       */
/*                                                                              */
/*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS         */
/*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT           */
/*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR       */
/*  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT        */
/*  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,       */
/*  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    */
/*  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR      */
/*  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF      */
/*  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING        */
/*  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS          */
/*  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                */
/********************************************************************************/
/*END of LICENSE NOTICE*/

/* RCS ID:$Id: props.c 1373 2008-08-27 20:51:52Z  $*/

/*Props subroutines.*/
/*Required definitions.*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "machine.h"
#include "blocks.h"
/*Diffusion Coefficient Calculation.*/

/*************************************/
/*Calculate the diffusion coefficient*/
/*depending upon temperature.        */
/*-- uses degrees C as input.        */
/*-- returns in m^2/ s$.             */
/*************************************/
extern CA_FLOAT global_pressure;

CA_FLOAT get_d_arr (CA_FLOAT temp, CA_FLOAT Qa, CA_FLOAT D_o)
{
      /************************/
  /* apply arrhenius eqn  */
  /* d = D0 exp(-Q/RT)    */
  /* input Q in J/Mol     */
  /*       D_o in M^2/sec */
  /*       temp in deg C  */
  /* output in M^2 / s    */
      /************************/

  CA_FLOAT tempK, D;

  tempK = temp + STD_TEMP;
  D = (CA_FLOAT) (D_o * (exp ((double) (-Qa / (8.31441 * tempK)))));
  /*printf("%.5g,%.5g,%.5g,%.5g,%.5g\n",Qa,D_o,temp,tempK,D) */

  return (D);
}

/* hard coded encapsulations for specific diffusion models */

CA_FLOAT get_ds_ti (CA_FLOAT temp)
{
  CA_FLOAT Qa, Do;

#define W13
#ifdef R61
  Do = 8.3e-2;                  /* in cm^2/s */
  Qa = 130;                     /* in KJ/mol */
#elif defined W13
  Do = 1.6;                     /* in cm^2/s */
  Qa = 201.7;                   /* in KJ/mol */
#endif
  return get_d_arr (temp, Qa * 1e3, Do * 1e-4);
}

CA_FLOAT get_dg_ti (CA_FLOAT temp)
{
  CA_FLOAT Qa, Do;

#define W13
#ifdef R61
  Do = 8.3e-2;                  /* in cm^2/s */
  Qa = 130;                     /* in KJ/mol */
#elif defined W13
  Do = 1.6;                     /* in cm^2/s */
  Qa = 201.7;                   /* in KJ/mol */
#endif
  Qa = Qa / 10;
  return get_d_arr (temp, Qa * 1e3, Do * 1e-4);
}

CA_FLOAT get_dl_ti (CA_FLOAT temp)
{
  CA_FLOAT Qa, Do;

  Do = 100.6;                   /* in cm^2/s */
  Qa = 201.7;                   /* in KJ/mol */
  return get_d_arr (temp, Qa * 1e3, Do * 1e-4);
}

/************************************************************/
CA_FLOAT get_reac_rate (CA_FLOAT temp)
{
  CA_FLOAT tempK, A, Q;
  CA_FLOAT rate;

   /*************************/
  A = 1e-5;
  Q = 274;
  tempK = temp + STD_TEMP;
   /*************************/
  rate = A * (EXP (Q / (GAS_CONST_SI * tempK)));
  return (rate);
}

/***********************************************************/

/* see Sung, Poirer et al, J Cr. Gro 226:363 (2001) for some */
/* values that may be useful for IN718 */
#if defined CELL_DIFF_COEFF
CA_FLOAT get_ds (CA_FLOAT temp)
{

  CA_FLOAT tempK, ds_cm, ds;

  tempK = temp + STD_TEMP;
  ds_cm = 1.1e-1 * EXP (-4922 / tempK);
  ds = ds_cm * 1e-4;
  return (ds);
}                               /*end of |get_ds| funciton */

CA_FLOAT get_dl (CA_FLOAT temp)
{

  CA_FLOAT tempK, dl_cm, dl;

  tempK = temp + STD_TEMP;
  dl_cm = 3.8e-2 * EXP (-2315 / tempK);
  dl = dl_cm * 1e-4;
  return (dl);
}                               /*end of |get_dl| function */
#elif defined CELL_DIFF_ARR
CA_FLOAT get_ds (CA_FLOAT temp, Solute_props * sp)
{

  CA_FLOAT ds;

  ds = get_d_arr (temp, sp->QaSol[0], sp->DoSol[0]);
  return (ds);
}                               /*end of |get_ds| funciton */

CA_FLOAT get_dl (CA_FLOAT temp, Solute_props * sp)
{

  CA_FLOAT dl;

  dl = get_d_arr (temp, sp->QaLiq, sp->DoLiq);
  return (dl);
}                               /*end of |get_dl| function */

#else

/* Clones of the function for */
/* testing with constand D coeff.*/

CA_FLOAT get_ds (CA_FLOAT temp)
{

  return (TEST_DS);

}                               /*end of TEST |get_ds| function */

CA_FLOAT get_dl (CA_FLOAT temp)
{

  return (TEST_DL);

}                               /*end of TEST |get_dl| funciton */

#endif /*|CELL_DIFF_COEFF| */

/*Average Diff Coeff. */

/*Uses a mixture rule to calculate averaged diffusion coefficient*/
/*for the partly solid region.*/

/*  Rewritten to reduce number of division operations*/
/*  this line is the biggest single bottleneck*/
/*  it would seem!*/
CA_FLOAT getav_d (CA_FLOAT dl, CA_FLOAT ds, CA_FLOAT fs)
{
  CA_FLOAT av_d;

  /*|av_d = 1/(fs/ds + ((1-fs)/dl));| */
  av_d = (ds * dl) / (fs * dl + (1 - fs) * ds);

  return (av_d);
}                               /*end of |getav_d| function */

/*Schiel Fraction solid calculation.*/
/*The Schiel equation for fraction solid as function of temperature:*/
/*$$ f_l = \biggl({T-T_p \over mC_o}\biggr)^{1 \over {k-1}}  $$*/
/*is used here.*/

/*Terminates at eutectic $f_s$, by switching to linear over a 2 degree range.*/

/*Eutectic $f_s$ should be calculated by calling this subroutine with the*/
/*input eutectic temperature as the |Tcell|, and 1.0 as the |eut_fs|, this*/
/*should be done outside the loop and stored in the block*/
/*as it does not change!.*/

CA_FLOAT schiel (CA_FLOAT Tcell)
{
  fprintf (stderr, "Sorry Schiel mode not available");
  exit (0);
#ifdef JUNK
  CA_FLOAT schiel_fs, schiel_fl;

  if (Tcell > T_LIQ) {
    return (0.0);
  }
#ifdef NLM_FORMULA
  else if (Tcell < T_SOL) {
    schiel_fs = 1.0;
  } else {
    schiel_fs = POW ((T_LIQ - Tcell) / (T_LIQ - T_SOL), ALLOY_EXPONENT);
  }
#else /*nomral scheil equation */
  if (Tcell > T_EUT) {
    schiel_fl = POW (((Tcell - TP) / (PD_SLOPE * CINITB)), (1 / (KB - 1)));
    schiel_fs = 1.0 - schiel_fl;
  } else if (Tcell > T_EUT - 2) {
    schiel_fs = FS_EUT + (T_EUT - Tcell) * ((1 - FS_EUT) / 2);
  } else {
    schiel_fs = 1.0;
  }
#endif /*NLM_FORMULA */
  return (schiel_fs);
#endif
}                               /*end of schiel function */

/*Schiel concentration eqn.*/
/* Find the concentration of solute B (alloy) according to the Schiel eqn. */
/*assuming the fl is known. This may be calculated from subroutine |schiel| but */
/*more generally may come out of some other growth model. */
CA_FLOAT find_sch_conc (CA_FLOAT tempK, CA_FLOAT fl)
{
  fprintf (stderr, "Sorry Schiel mode not available");
  exit (0);
#ifdef JUNK
  CA_FLOAT sch_c;
  CA_FLOAT Tcell;

  Tcell = tempK - STD_TEMP;
#ifdef NLM_FORMULA
  sch_c = CINITB - (T_LIQ - Tcell) / PD_SLOPE;
#else
  sch_c = CINITB * POW ((fl), (KB - 1));
#endif
  if (sch_c < CINITB)
    sch_c = CINITB;
  if (sch_c > MAX_B_CONC)
    sch_c = MAX_B_CONC;
  return (sch_c);
#endif
}

/* Growth subroutine.*/

/*This holds the equation that is used to calculate the growth. It may*/
/*be modified..$$ \delta{}x=GT_u^2\delta{}t$$.*/
/****************************************************************/
/* Subroutine to calc. growth as fn of undercooling             */
/* and return in mm                                             */
/****************************************************************/
/*dt is the time increment $\delta{}t$ and Tunder is the undercooling $T_u$ */
/* |gg_const| is the grain growth constant $G$ supplied by the user. */
/* From Jackson and Hunt */
CA_FLOAT growth (CA_FLOAT gg_const, CA_FLOAT gg_cub, CA_FLOAT dt, CA_FLOAT Tunder)
{

  if (gg_const == 0) {
    return (dt * FLAT_GROWTH);
  } else {
    return (dt * (gg_const * Tunder * Tunder + gg_cub * Tunder * Tunder * Tunder));
  }

}                               /* end of growth function */

/*  multi-component growth */
CA_FLOAT growth_primary (CA_FLOAT gg_const, CA_FLOAT gg_cub, CA_FLOAT dt, CA_FLOAT Tunder, int eut_flag, CA_FLOAT cell_fs_n_eut,
                         CA_FLOAT cell_fs_b_eut)
{

  CA_FLOAT ins_growth_p;

  if ((gg_const == 0) && (gg_cub == 0)) {
    ins_growth_p = dt * FLAT_GROWTH;
  } else {
    ins_growth_p = dt * (gg_const * Tunder * Tunder);
  }

  return (ins_growth_p);

}                               /* end of primary growth function */

/************************************************************/

CA_FLOAT growth_eutectic (CA_FLOAT gg_const, CA_FLOAT gg_cub, CA_FLOAT dt, CA_FLOAT Tunder, int eut_flag, CA_FLOAT cell_fs_n_eut,
                          CA_FLOAT cell_fs_b_eut)
{

  CA_FLOAT ins_growth_e;

  if ((gg_const == 0) && (gg_cub == 0)) {
    ins_growth_e = dt * FLAT_GROWTH;
  } else {
    ins_growth_e = dt * ((gg_const * Tunder * Tunder) + (gg_cub * Tunder * Tunder * Tunder));
  }

  return (ins_growth_e);

}                               /* end of eutectic growth function */

/************************************************************/

/**********************/
/*Find Local Liquidus.*/
/**********************/
/*Use appropriately simplified phase diagram to find the liquidus for the cell.*/
/*Used if |phase_diag_on| option set*/
CA_FLOAT cell_liq_calc (CA_FLOAT cell_conc, Mat_str * mp)
{
  CA_FLOAT TliqCell;
 
  /* assuming the gas has no effect */
  if (mp->gasprops.m_solute[0] == 0) {
    TliqCell = mp->tp + (mp->alloyprops[0].m_solute[0] * cell_conc);
    if (TliqCell <= mp->alloyprops[0].T_eut)
    TliqCell = mp->alloyprops[0].T_eut;
  } else {
    fprintf (stderr, "SORRY:cell_liq_calc: gas effect on liquidus not implemented.\n");
    fprintf (stderr, "Please rewrite cell_liq_calc if that is what you need.\n");
    fprintf (stderr, "Otherwise set m_solute (liquidus slope) to zero for the gas.\n");
    exit (0);
  }
  return (TliqCell);
}

/* THUINET 0106 */
CA_FLOAT cell_liq_calc_poly (CA_FLOAT cell_conc, Mat_str * mp, int ele_1)
{
  CA_FLOAT TliqCell;
  int i;

  /* assuming the gas has no effect */
  if (mp->gasprops.m_solute[0] == 0) {
    TliqCell = mp->tp_poly[0];
    for (i = 0; i < ele_1; i++) {
      TliqCell += (mp->alloyprops[i].m_solute[0] * cell_conc);
    }
  } else {
    fprintf (stderr, "SORRY:cell_liq_calc: gas effect on liquidus not implemented.\n");
    fprintf (stderr, "Please rewrite cell_liq_calc if that is what you need.\n");
    fprintf (stderr, "Otherwise set m_solute (liquidus slope) to zero for the gas.\n");
    exit (0);
  }
  return (TliqCell);
}




/*Find local Liquidus using the approximated multi component phase diagram*/
/* a hyper surface is used */


/*Find H saturation.*/

/*This uses the formula */
/*$$ \log S_l = - {2760 \over T} + 2.796 $$ which */
/*is on page 3 (equation 2.2) of PDL thesis,*/
/*referenced to Ransley and Neufeld [1948].*/
/*and the correction factor for si:*/
/*$$ f = 10^{(-0.0119 * \%Si)}$$*/
/*eqn 2.5 PDL thesis, ref. to Doutre(1991).*/

/*Currently assuming constant partition coefficient for H hence ignoring */
/*the equation for soulubility in the solid. */
CA_FLOAT find_sat (Mat_str * mp, CA_FLOAT cell_tempK, CA_FLOAT cell_si, CA_FLOAT cell_fs)
{
  CA_FLOAT logss, sl, ss, logsl, f, cell_sat, cell_satMCUB;
  const CA_FLOAT mpm_cub = MPMETERCUB;

  f = POW (10, (DOUTRE_F * (MIN (cell_si, mp->alloyprops[0].c_eut))));
#ifdef RANSLEY
  logsl = -(2760 / cell_tempK) + 2.796;
  sl = POW (10, logsl);         /* ml/100g STP */
#else /*IMABAYASHI*/
    logsl = -(2392 / cell_tempK) + 3.256;
  sl = POW (10, logsl);         /* cc/1000g STP */
  sl *= 0.1;
#endif /*Ransley/Imabayashi */
  cell_sat = f * sl;
  cell_satMCUB = cell_sat * MPMETERCUB; /*conv to SI $mol\over m^3$ */

  cell_satMCUB *= SQRT (global_pressure);

  return (cell_satMCUB);
}

/************************the modified correction factor for********/
/************************multi component diffusion*****************/
#ifdef SENS
CA_FLOAT find_sat_multi (BB_struct * bp, CA_FLOAT cell_tempK, CA_FLOAT ** conc_multi, CA_FLOAT cell_fs, CA_FLOAT sens_P, int index)
#else
CA_FLOAT find_sat_multi (BB_struct * bp, CA_FLOAT cell_tempK, CA_FLOAT ** conc_multi, CA_FLOAT cell_fs, int index)
#endif
{
  CA_FLOAT logss, sl, ss, logsl, f, cell_sat, cell_satMCUB;
  MultiS_struct *ms;
  Ctrl_str *cp;
  int ele_num, ele_1;
  int i, j;
  CA_FLOAT sum;
  CA_FLOAT doutre[3], local_max[3];

#ifdef MULTICOMP
  local_max[0] = MAX_B_CONC_1;
  local_max[1] = MAX_B_CONC_2;
  local_max[2] = MAX_B_CONC_2;
#endif

  doutre[0] = -0.0119; /* silicon */
  doutre[1] = -0.0269; /* copper  */
  doutre[2] = 0.017;   /* magnesium*/ /* from Alcan alscan documentation */

  cp = bp->ctrl;
  ms = &(bp->MultiSvals);

  ele_num = cp->NUM_COMP;
  ele_1 = ele_num - 1;

  sum = 0.0;

/************************************************************/
  for (i = 0; i < ele_1; i++) {
    sum += doutre[i] * (MIN (conc_multi[i][index], local_max[i]));
  }
  f = POW (10, sum);
#ifdef RANSLEY
  logsl = -(2760 / cell_tempK) + 2.796;
  sl = POW (10, logsl);         /* ml/100g STP */
#else /*Imbayashi */
  logsl = -(2392 / cell_tempK) + 3.256;
  sl = POW (10, logsl);         /* cc/1000g STP */
  sl *= 0.1;
#endif /*Ransley/Imabayashi */
  cell_sat = 1.0 * f * sl;
  cell_satMCUB = cell_sat * MPMETERCUB; /*conv to SI $mol\over m^3$ */

  cell_satMCUB *= SQRT (global_pressure);
#ifdef SENS
  cell_satMCUB *= SQRT (sens_P);
#endif

  return (cell_satMCUB);
}


CA_FLOAT find_sat_poly (BB_struct * bp, int sbnum, CA_FLOAT cell_tempK)
{
  CA_FLOAT logss, sl, ss, logsl, f, cell_sat, cell_satMCUB;
  Ctrl_str *cp;
  SB_struct *sp;
  int ele_num, ele_1;
  int isol;
  CA_FLOAT sum;
  CA_FLOAT *ncl_poly[NSOLMAX];

  sp = bp->sb[sbnum];           /*subblock pointer */

  cp = bp->ctrl;
  ele_num = cp->NUM_COMP;       /* number of elements in the alloy */
  ele_1 = ele_num - 1;

  for (isol = 0; isol < ele_1; isol++) {        /*loop on solutes by Ludovic THUINET */
    ncl_poly[isol] = sp->c_sol_poly[isol];
  }  /* end loop on solutes */

  sum = 0.0;

  for (isol= 0; isol < ele_1; isol++) {
    sum +=bp->mprops.alloyprops[isol].coef_doutre *(*ncl_poly[isol]);
  }
  f = POW (10, sum);
#ifdef RANSLEY
  logsl = -(2760 / cell_tempK) + 2.796;
  sl = POW (10, logsl);         /* ml/100g STP */
#else /*Imbayashi */
  logsl = -(2392 / cell_tempK) + 3.256;
  sl = POW (10, logsl);         /* cc/1000g STP */
  sl *= 0.1;
#endif /*Ransley/Imabayashi */
  cell_sat = 1.0 * f * sl;
  cell_satMCUB = cell_sat * MPMETERCUB; /*conv to SI $mol\over m^3$ */

  cell_satMCUB *= SQRT (global_pressure);
#ifdef SENS
  cell_satMCUB *= SQRT (sens_P);
#endif

  return (cell_satMCUB);
}

/*Find cubic root.*/

/*Radius finding routine by Peter Lee using Mathematica.*/
/* DISUSED : replaced by iterative (Newton method) in findroot.c */

/* Solve for the one real root of $ar^3+br^2+c=0$   		*/
/* as given analytically from Mathematica.\par{}*/
/*{\bf Inputs:} temperature (${\rm^\circ K}$), amount of gas*/
/*(mol),surface tension $(\gamma)$ (${N\over m}$ or ${pa . m}$), old radius (m).\par{}*/
/*{\bf Output:} new radius (m).*/

/* if only it did a better job of reducing...\par */
/* {\obeylines\parindent = 1cm \codefont*/
/* x1 = -b/(3*a) + cbrt(2.0)*pow(b,2)/ (3*a*cbrt(-2*pow(b,3) 	 */
/* 	- 27*pow(a,2)*c + pow(3,1.5)*a*sqrt(4*pow(b,3)*c + 	     */
/*	27*pow(a,2)*pow(c,2)))) + cbrt(-2*pow(b,3) -		              */
/*	27*pow(a,2)*c + pow(3,1.5)*a*sqrt(4*pow(b,3)*c +	            */
/*	27*pow(a,2)*pow(c,2)))/(3*cbrt(2.0)*a);			               */
/*}*/
/* Local definitions: (Numerical constants moved to |machine.h|) */
#ifdef SENS
CA_FLOAT r_from_ideal_gas (CA_FLOAT temperature, CA_FLOAT n_gmol, CA_FLOAT gammaSI, CA_FLOAT old_rSI, CA_FLOAT sens_PRESS)
#else
CA_FLOAT r_from_ideal_gas (CA_FLOAT temperature, CA_FLOAT n_gmol, CA_FLOAT gammaSI, CA_FLOAT old_rSI)
#endif
{
#define P_ATM		1.0     /* ambient pressure in atm */
#define TWO_GAMMA_ATM	0.00177646      /* surf.ten in Atm x cm =|G_PA/101325| */
#define LMT_SQRT	1.0e-18 /* exact soln only if |tmp>LMT_SQRT| */

  double a, b, c, tmp, tcbrt, ansCM, old_r;     /* must stay as double */
  extern double sqrt (double);
  extern double cbrt (double);

/*NOTE: if |P_ATM| and |GAMMA| remain constant, this can be speed	*/
/*up quite a bit, most becomes const...			*/

/* {\bf NOTE:} If |tmp| = $(4b^3c + 27a^2c^2)\to 0$, only an approximate	*/
/* 	root is used. Also, the value may be inexact near the	*/
/*	transition!		 				*/
  old_r = (double) (old_rSI * 100);

#ifdef JUNK
  a = P_ATM;
#endif
  a = global_pressure;
#ifdef SENS
  a *= sens_PRESS;
#endif
  b = 2 * gammaSI * CONV_ATM_CM;
  c = -THREE_BY_4PI * GAS_CONST_ATM * n_gmol * temperature;
  tmp = 4 * b * b * b * c + 27 * a * a * c * c;
  if (tmp > LMT_SQRT) {
    tcbrt = cbrt (-2 * b * b * b - 27 * a * a * c + POW_3_TO_1PT5 * a * sqrt (tmp));
    ansCM = (CA_FLOAT) (-b / (3.0 * a) + CBRT2 * b * b / (3.0 * a * tcbrt) + tcbrt / (3.0 * CBRT2 * a));
  } else
    ansCM = (CA_FLOAT) sqrt ((-c - a * old_r * old_r * old_r) / b);
  return (CA_FLOAT) (ansCM / 100);
}

/*Find maturity function.*/
/*not currentl used.*/

CA_FLOAT get_mat (CA_FLOAT thresh, CA_FLOAT delta_t, CA_FLOAT ss, CA_FLOAT oldmat)
{
  CA_FLOAT mat;

  mat = oldmat + delta_t * (ss - thresh);
  return (mat);
}

#ifdef TEST_PROPS
void main (void)
{
  CA_FLOAT t;

  for (t = 300; t < 1500; t += 50) {
    printf (t, get_ds_ti (t));
  }
}
#endif /*TEST_PROPS */
/* rcs id subroutine.*/
/*Little subroutine to include the |rcs Id| in the object code.*/
/*This can also be called to print or include the Id in a file.*/

char const *rcs_id_props_c ()
{
  static char const rcsid[] = "$Id: props.c 1373 2008-08-27 20:51:52Z  $";

  return (rcsid);
}

/* RCS ID:$Id: props.c 1373 2008-08-27 20:51:52Z  $*/
/*
*/
