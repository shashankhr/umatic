
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

/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/****************************************************************/
/*      sb_decentred_step_poly.c:                                           */
/*                                                              */
/* The main subroutines of decentred square/octahedron algorithm */
/* adapted for poly-component multiphase system                 */
/* These routines initially written by L. THUINET               */
/*                                                              */
/*                                                              */
/*                                                              */
/****************************************************************/
/*RCS ID: $Id: sb_decentred_step_poly.c 1356 2008-08-18 13:41:15Z  $*/
/*****************************************************************/
/*Extension of stochastique cell nucleation to multiphase system */
/*by THUINET 03/05                                               */
/*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "machine.h"
#include "blocks.h"
#include "umat_matrix.h"
#include "props.h"
/* prototypes for the nucleation function choices */
#include "nucfuncs.h"

/* Prototype defs */
/* sb_nuc.c*/
extern int init_fixed_angle_grain (BB_struct * bp, int igr, int sbnum, int icell, int nc, CA_FLOAT a0, CA_FLOAT a1, CA_FLOAT a2);

/*??*/
extern void add_to_grain (Ind_grain * gp, int i, int j, int k);

/*from trans_interp_calc.c*/
extern CA_FLOAT trans_interp_calc (FGrid_str * fg, BB_struct * bp, int sbnum, int x, int y);

/* from sb_temp_calc.c */
extern CA_FLOAT fg_temp_calc (BB_struct * bp, int sbnum, int x, int y);

/*from sb_nuc.c*/
extern int init_new_grain (BB_struct * bp, int igr, int sbnum, int xcell, int ycell, int zcell, int nc);

/*from props*/
extern CA_FLOAT growth (CA_FLOAT gg_const, CA_FLOAT gg_cub, CA_FLOAT dt, CA_FLOAT Tunder);

/*for seed melt_back */
extern int cut_seed (BB_struct * bp, int sbnum);

int cell_nucleation_poly (BB_struct * bp, int sbnum)
{
  int errors = 0;
  int global_undercooling = 0;

/*THUINET 02/05*/
  int isol, ele_num, ele_1;
  int iphs, iphs_tot;
  CA_FLOAT *ncl_poly[NSOLMAX];
  CA_FLOAT C_ini_poly[NSOLMAX];

/*FIN THUINET 02/05*/
#ifdef OLD_TUNDER
  CA_FLOAT *otu, *otup;
#endif /*OLD_TUNDER */
  int i, j, k;
  int nx, ny, nz, skip;
  int fixed_nuc, i_nuc, icell, index;
  SB_struct *sp;
  Ctrl_str *cp = bp->ctrl;

/*THUINET 04/05*/
  /*CA_FLOAT Tliq,  Tunder,TliqFixed, min_temp; */
  CA_FLOAT Tliq_poly[NPHAMAX], Tunder, Tunder0, Tunder1;
  CA_FLOAT min_temp;
  int *nat_site;
  int *nat_cell_p, *onat_cell_p;
  int *nat_grain_p, *onat_grain_p;

/*END THUINET 04/05*/
  CA_FLOAT a0, a1, a2, liq_slope, c_ini;
  CA_FLOAT *thrp, *c_temp_p, *ncl, *nce;

/*THUINET 04/05*/
  int *ngr, *ogr;

/*END THUINET 04/05*/
  int phase_diag_on, cell_nuc_on, cell_temp_on, n_neigh;
  int **gv_p = bp->gr_array;    /* grain block value array pointer pointer :-/ */
  static int nucmsg = 0;
  int (*cell_nuc_func) (BB_struct * bp, CA_FLOAT Tunder, CA_FLOAT argthree);

/*added by THUINET 10/05*/
  int indic, n, i_nn, j_nn, k_nn, ijk_nn, gr_nn, nat_grain_nn;

/* nearest neighbours structure */
  struct nns {
    int nx;
    int ny;
    int nz;
  } nn[] = {
  1, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 1, 1, 0, -1, 1, 0, -1, -1, 0, 1, -1, 0};
  int NNB = 6;
  int nbskip[6];

/*End THUINET 10/05*/

  /*define variables for decentred algorithm */
  /*by Wei Wang on 09-01-01 */
  sp = bp->sb[sbnum];
  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  skip = 2 * (nx + 2);

/*THUINET 10/05*/

  for (n = 0; n < NNB; n++) {
    i_nn = nn[n].nx;
    j_nn = nn[n].ny;
    k_nn = nn[n].nz;
    nbskip[n] = i_nn + j_nn * (nx + 2) + k_nn * (nx + 2) * (ny + 2);
  }
/*END THUINET 10/05*/

  c_temp_p = sp->c_temp;

  ele_num = cp->NUM_COMP;       /* number of elements in the alloy */
  ele_1 = ele_num - 1;

/*THUINET 04/05*/
  iphs_tot = cp->NUM_PHS;
  nat_site = sp->nat_sol_site;

  onat_cell_p = bp->itmp_nat_cell;
  nat_cell_p = sp->nat_cell;
  icopy_matrix (PAD, onat_cell_p, nat_cell_p, bp, gv_p, sbnum);
  onat_cell_p += bp->cubeptr.flist[0][START];

  onat_grain_p = bp->itmp_nat_grain;
  nat_grain_p = sp->nat_grain;
  icopy_matrix (PAD, onat_grain_p, nat_grain_p, bp, gv_p, sbnum);
  onat_grain_p += bp->cubeptr.flist[0][START];


/*END THUINET*/

  ogr = bp->itmp_one;           /* array of old cell grain, size (nc[i]+2)^3 */
  ngr = sp->gr;                 /* array of new cell grain, size (nc[i])^3   */

  /* copy grain number to temporary buffer */
  icopy_matrix (PAD, ogr, ngr, bp, gv_p, sbnum);

  /*move the ptr to outside corner to inside corner */
  ogr += bp->cubeptr.flist[0][START];

  for (isol = 0; isol < ele_1; isol++) {
    C_ini_poly[isol] = bp->mprops.alloyprops[isol].Cinit;
    ncl_poly[isol] = sp->c_sol_poly[isol];
  }

  for (iphs = 0; iphs < iphs_tot; iphs++) {
    Tliq_poly[iphs] = bp->mprops.Tliq_poly[iphs];
  }

  /* set up phase diagram mode */
  phase_diag_on = bp->ctrl->phase_diag_on;
  global_undercooling = (phase_diag_on && bp->ctrl->global_undercooling);

  /* Set up local pointers */
#ifdef OLD_TUNDER
  otup = otu = bp->old_Tunder;
#endif /*OLD_TUNDER */

  /* Set up option flags */
  cell_nuc_on = cp->use_cell_nuc;
  cell_temp_on = cp->use_cell_temp;

  if (cell_nuc_on == TRUE) {    /*Select the cell nuc function -- should move to higher level! */
    switch (bp->nprops.nmodel) {
    case N_RAPPAZ:             /* use the PDL coded Rappaz based method, 'rolling dice' per cell... */
      cell_nuc_func = cell_nuc;
      break;
    case N_DIST:               /* use the DISTRIBUTION lookup-table method */
      cell_nuc_func = dist_cell_nuc;
      break;
    case N_HETERO:
      fprintf (stderr, "ERROR: sb_umat_step: cannot find cell_nuc function for N_HETERO model\n");
      exit (0);
      break;
    case N_RAP_DEL:
#           ifndef OLD_TUNDER
      fprintf (stderr, "ERROR:sb_umat_step: At present you have to define OLD_TUNDER to use N_RAP_DEL nucleation model.\n");
      exit (0);
#           endif /*OLD_TUNDER */
      cell_nuc_func = del_cell_nuc_norate;
      break;
    case N_OLDF_DEL:
      fprintf (stderr, "ERROR: sb_umat_step: cannot find cell_nuc function for N_OLDF_DEL model\n");
      exit (0);
      break;
    case N_BLOCK:              /* use method of setting a threshold for each cell */
      /* similar to the method used for pores */
      cell_nuc_func = block_cell_nuc;
      thrp = sp->c_nuc_thresh;
      break;
    default:
      fprintf (stderr, "ERROR: sb_umat_step: cannot find cell_nuc function; Nucleation model undefined.(%i)\n", (bp->nprops.nmodel));
      exit (0);
      break;
    }
  }

  /*end of set up cell nuc function */
 /***********************************************/
  /* Beginning of main loop(s)                   */
 /***********************************************/
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {

        if (*ngr == LIQUID) {   /* cell grain number =0 */

          if (*nat_site == -1) {
            Tunder = -1.0;
          }

          if (*nat_site == 0) { /* nucleate the primary */
            Tunder = Tliq_poly[(*nat_site)] - *c_temp_p;
            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes by Ludovic THUINET */
              Tunder += bp->mprops.alloyprops[isol].m_solute[(*nat_site)] * (*(ncl_poly[isol]) - C_ini_poly[isol]);
            }

          }

          if (*nat_site == 1) { /* nucleate the eutectic */

            Tunder0 = Tliq_poly[0] - *c_temp_p;
            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes by Ludovic THUINET */
              Tunder0 += bp->mprops.alloyprops[isol].m_solute[0] * (*(ncl_poly[isol]) - C_ini_poly[isol]);
            }

            Tunder1 = Tliq_poly[1] - *c_temp_p;
            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes by Ludovic THUINET */
              Tunder1 += bp->mprops.alloyprops[isol].m_solute[1] * (*(ncl_poly[isol]) - C_ini_poly[isol]);
            }

            if (Tunder0 > Tunder1) {
              Tunder = Tunder1;
            } else {
              Tunder = Tunder0;
            }

          }

          if (*nat_site == 2) { /* nucleation from the primary phase */

            indic = 0;

            for (n = 0; n < NNB; n++) { /* neighbour loop */
              i_nn = nn[n].nx;
              j_nn = nn[n].ny;
              k_nn = nn[n].nz;
              ijk_nn = nbskip[n];
              nat_grain_nn = *(onat_grain_p + ijk_nn);

              /* type 2 eutectic nucleation only activates if a primary is neighbouring */
              if ((nat_grain_nn == 1) && (indic == 0)) { /*one of the neighbour is the mushy primary phase */
                indic = 1;
                gr_nn = *(ogr + ijk_nn);

                /*Calculation of Tunder */

                Tunder0 = Tliq_poly[0] - *c_temp_p;
                for (isol = 0; isol < ele_1; isol++) {  /*loop on solutes by Ludovic THUINET */
                  Tunder0 += bp->mprops.alloyprops[isol].m_solute[0] * (*(ncl_poly[isol]) - C_ini_poly[isol]);
                }

                Tunder1 = Tliq_poly[1] - *c_temp_p;
                for (isol = 0; isol < ele_1; isol++) {  /*loop on solutes by Ludovic THUINET */
                  Tunder1 += bp->mprops.alloyprops[isol].m_solute[1] * (*(ncl_poly[isol]) - C_ini_poly[isol]);
                }

                if (Tunder0 > Tunder1) {
                  Tunder = Tunder1;
                } else {
                  Tunder = Tunder0;
                }
              }
            }                   /*End of loop on neighbour */

            if (indic == 0) {   /*No neighbour in the mushy primary phase */
              Tunder = -1;
            }

          }
          /*end if nat_site=2 */
          if (cell_nuc_on) {

            if (Tunder > 0 && (*cell_nuc_func) (bp, Tunder, *thrp) == 1) {
              /* threshold is exceeded so nucleate a grain */
              /* error checks */

              if (nucmsg < MAX_NUC_MSG) {
                fprintf (stderr, "cell temp[%g] cell Tliq[%g]\n", *c_temp_p, Tliq_poly[(*nat_site)]);
                nucmsg++;
              }

              if (bp->nprops.ngr == bp->nprops.gd_max_total) {
                static int emess=0;
                if (emess <= MAX_ERR_MSG){
                   fprintf (stderr,
                         "ERROR:sb_decentred_step_poly.c-cell_nucleation(): Max # grains set by user [%d] exceeded. Increase option MaxTotGrains.\n",
                         bp->nprops.gd_max_total);
                         emess ++;
                 }else{ /* maximum error messages exceeded */
                    fprintf(stderr,"ERROR: maximum error messages %i exceeded\nQuitting...\n",emess);
                    exit(0);
                }

              } else {
                /* no errors */
                /* create new grain and set up */
                /*record nucleation values in structure */
                bp->nprops.ngr++;
                bp->sb[sbnum]->ngr++;
                *ngr = bp->nprops.ngr;

                if (*nat_site == -1) {
                  *(nat_cell_p) = 0;
                  *(nat_grain_p) = 0;
                }
                if (*nat_site == 0) {
                  *(nat_cell_p) = 1;
                  *(nat_grain_p) = 1;
                }
                if (*nat_site == 1) {
                  *(nat_cell_p) = 2;
                  *(nat_grain_p) = 2; /* eutectic grain */
                }
                if (*nat_site == 2) {
                  *(nat_cell_p) = 2;
                  *(nat_grain_p) = 2; /* eutectic grain (nucleated from primary) */
                }

                icell = i + j * nx + k * nx * ny;
                fixed_nuc = 0;

                for (i_nuc = 0; i_nuc < cp->nfnuc; i_nuc++) {
                  index = (int) (cp->nsite[i_nuc][0] + cp->nsite[i_nuc][1] * nx + cp->nsite[i_nuc][2] * nx * ny);
                  if (icell == index) {
                    a0 = PI_BY_FOUR / 45.0 * cp->nsite[i_nuc][3];
                    a1 = PI_BY_FOUR / 45.0 * cp->nsite[i_nuc][4];
                    a2 = PI_BY_FOUR / 45.0 * cp->nsite[i_nuc][5];
                    init_fixed_angle_grain (bp, bp->nprops.ngr, sbnum, index, 1, a0, a1, a2);
                    fprintf (stderr, "sb_decentred_step.c-cell_nucleation(): nucleating fixed gr[%d] at cell %d.\n", bp->nprops.ngr,
                             index);
                    fixed_nuc = 1;
                  }
                }

/*                 if(fixed_nuc == 0){ 
                   init_new_grain(bp, bp->nprops.ngr, sbnum, i ,j,k,1);        
		 }*/

/*modified by thuinet 10/05*/

                if ((fixed_nuc == 0) && (*nat_site != 2)) {
                  init_new_grain (bp, bp->nprops.ngr, sbnum, i, j, k, 1);
                }

                if ((fixed_nuc == 0) && (*nat_site == 2)) {
                  a0 = bp->gr[gr_nn]->ang[0];
                  a1 = bp->gr[gr_nn]->ang[1];
                  a2 = bp->gr[gr_nn]->ang[2];
                  init_fixed_angle_grain (bp, bp->nprops.ngr, sbnum, index, 1, a0, a1, a2);
                }

/*End modification*/

                bp->gr[bp->nprops.ngr]->TNuc = *c_temp_p;
                bp->gr[bp->nprops.ngr]->TunderNuc = Tunder;
                /* bp->gr[bp->nprops.ngr]->CellConcNuc = phase_diag_on ? *cell_conc : 0; */

              }
            }
          }
        } /* end if liquid */
#ifdef OLD_TUNDER
        otup++;
#endif /*OLD_TUNDER */
        c_temp_p++;
        thrp++;

/*THUINET 05/05*/
        nat_site++;
        onat_cell_p++;
        nat_cell_p++;
        onat_grain_p++;
        nat_grain_p++;

/*END THUINET 05/05*/

        ogr++;                  /* ptr to outside block cell grain number */
        ngr++;                  /* ptr to inside  block cell grain number */

        for (isol = 0; isol < ele_1; isol++) {
          ncl_poly[isol]++;
        }

      }                         /*end of I loop */

      ogr += 2;
      onat_cell_p += 2;
      onat_grain_p += 2;

    }                           /* end of J loop */

    ogr += skip;
    onat_cell_p += skip;
    onat_grain_p += skip;

  }                             /* end of K loop */

  if (sp->ncsolid >= bp->ncsb) {
    fprintf (stderr, "sb_decentred_step.c-cell_nucleation(): SB#%d completely solid.\n", sbnum);
    sp->done = TRUE;
  }
  return (errors);
}


/******************************************************************/
/* fraction solid change with diffusion of multicomponent alloy   */
/* by Ludovic THUINET on 07-11-01                                 */
/******************************************************************/

int fs_change_diffuse_poly (BB_struct * bp, int sbnum) {
  int i, j, k, nx, ny, nz, skip;
  int *ngr, *ogr;

/* THUINET 02/05 and 04/05 */
  int i1, i2;
  int isol;
  int ele_num, ele_1, noninv;
  int ieq, ipheq_tot;
  int dim_tab;
  int iphs, iphs_tot;
  CA_FLOAT *oce_poly[NSOLMAX], *nce_poly[NSOLMAX], *ncl_poly[NSOLMAX];
  CA_FLOAT dce_poly[NSOLMAX];
  CA_FLOAT *nfs_poly[NPHAMAX];

/* End THUINET 02/05 and 04/05*/
  SB_struct *sp;
  Ctrl_str *cp;
  CA_FLOAT *ofs, *nfs, *ncv;

/*multiphase THUINET 04/05*/
  CA_FLOAT *ofs_poly[NPHAMAX];
  CA_FLOAT delta_T;
  int *nat_cell_p, *onat_cell_p;
  int *nat_grain_p, *onat_grain_p;
  CA_FLOAT oTliq0,oTliq1;

/*END THUINET 04/05*/
  CA_FLOAT *ncl, *nce;
  CA_FLOAT *oce;
  CA_FLOAT *nct;
  CA_FLOAT part_coef, liq_slope;
  CA_FLOAT Tcell, km, k_inv, km_inv, m_inv;
  CA_FLOAT temp_pure;
  CA_FLOAT coef_curv;
  CA_FLOAT D_x, D_xt, fs, c_int;
  CA_FLOAT dce, dcl, dfs;
  CA_FLOAT sum_del_fs;

/*variable declarations added by THUINET 03-02-05*/

  CA_FLOAT oTliq;
  void invsys (BB_struct * bp, int sbnum);

/*end variable declarations added by THUINET 03-02-05*/

  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  skip = 2 * (nx + 2);

  cp = bp->ctrl;
  ele_num = cp->NUM_COMP;       /* number of elements in the alloy */
  ele_1 = ele_num - 1;
  iphs_tot = cp->NUM_PHS;

  temp_pure = bp->mprops.tp;
  coef_curv = bp->mprops.gibbs_thomson / bp->size_c[0];

  sp = bp->sb[sbnum];

  ngr = sp->gr;

  onat_cell_p = bp->itmp_nat_cell;
  nat_cell_p = sp->nat_cell;

  onat_cell_p += bp->cubeptr.flist[0][START];

  onat_grain_p = bp->itmp_nat_grain;
  nat_grain_p = sp->nat_grain;

  onat_grain_p += bp->cubeptr.flist[0][START];

  nfs = sp->c_fs;
  /*   ncl = sp -> c_sol_alloy;*//*Ludovic THUINET 02-02-05 */
  /*   nce = sp -> c_eqv_alloy;*//*Ludovic THUINET 02-02-05 */
  ncv = sp->curv;

  nct = sp->c_temp;

  ofs = bp->ftmp_one;
  ofs += bp->cubeptr.flist[0][START];

  for (iphs = 0; iphs < iphs_tot; iphs++) {
    nfs_poly[iphs] = sp->c_fs_poly[iphs];
    ofs_poly[iphs] = bp->ftmp_one_poly[iphs];
    ofs_poly[iphs] += bp->cubeptr.flist[0][START];
  }

  /*   oce = bp->ftmp_five;*//*Ludovic THUINET 02-02-05 */
  ogr = bp->itmp_one;

  /*   oce += bp->cubeptr.flist[0][START];*//*Ludovic THUINET 02-02-05 */
  ogr += bp->cubeptr.flist[0][START];

  for (isol = 0; isol < ele_1; isol++) {        /*loop on solutes by Ludovic THUINET */
    ncl_poly[isol] = sp->c_sol_poly[isol];
    nce_poly[isol] = sp->c_eqv_poly[isol];
    oce_poly[isol] = bp->ftmp_ce_poly[isol];
    oce_poly[isol] += bp->cubeptr.flist[0][START];

  }                             /* end loop on solutes */

  sum_del_fs = 0.0;

  /*BEGIN loop through all cells */
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
        Tcell = *nct;
        fs = *ofs;
        if (*ngr == 0) {

          /* CASE 0 : the cell is LIQUID or NOT-CASTING */
          /** \todo SHOULD NOT calculate this if it is NOT-CASTING! */

          for (isol = 0; isol < ele_1; isol++) {        /*loop on solutes */
            *ncl_poly[isol] = *nce_poly[isol];  /* C_E = C_S*fs + C_L*(1-fs), fs=0 */
          }
        }

        if (*nat_cell_p == 1) {
          /* CASE 1 : liquid + primary phase */

          if (*(ofs_poly[0]) < 1) {
            /* initialize the matrix */
            bp->dim_tab = ele_1;
            for (i1 = 0; i1 < NPHAMAX; i1++) {
              for (i2 = 0; i2 < NPHAMAX; i2++) {
                bp->a[i1][i2] = 0.0;
              }
            }

            /*calculations of the cofficients */
            coef_curv = bp->mprops.gibbs_thomson / bp->size_c[0];
            /*mass conservation equation for each solute */

            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes */
              dce_poly[isol] = *nce_poly[isol] - *oce_poly[isol];
              bp->a[isol][isol] = 1.0;
              bp->a[isol][isol] -= (1.0 - bp->mprops.alloyprops[isol].part_coef[0]) * (*ofs_poly[0]);
              bp->a[isol][ele_1] = -(1.0 - bp->mprops.alloyprops[isol].part_coef[0]) * (*ncl_poly[isol]);
              bp->a[isol][ele_1 + 1] = dce_poly[isol];
            }

            /* Thermodynamic equilibria */

            /* first case : constant thermodynamic data */
            oTliq = bp->mprops.tp_poly[0];
            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes */
              oTliq += bp->mprops.alloyprops[isol].m_solute[0] * (*ncl_poly[isol]);
            }
            /* second case : coupling with Thermo-Calc */
            /* end liquidus temperature calculation */
            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes */
              bp->a[ele_1][isol] = bp->mprops.alloyprops[isol].m_solute[0];
            }
            bp->a[ele_1][ele_1 + 1] = Tcell + coef_curv * (*ncv) - oTliq;

            invsys (bp, sbnum);

            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes */
              *ncl_poly[isol] += bp->a[isol][ele_1 + 1];
            }

            *ofs_poly[0] += bp->a[ele_1][ele_1 + 1];
            *ofs_poly[0] = MIN (1.0, *ofs_poly[0]);

            /*test to know if the equivalent composition lies in the tie triangle  */
            oTliq1=bp->mprops.tp_poly[1];
            for (isol = 0; isol < ele_1; isol++) { /*loop on solutes */
              oTliq1 +=bp->mprops.alloyprops[isol].m_solute[1]*(*ncl_poly[isol]);
            }

            /* change the state to eutectic */
            if ((Tcell < oTliq1)&&(*nat_grain_p == 2)) {
              *nat_cell_p=2;
              *onat_cell_p=2;
            }
            /*end of test*/

            /* remove a cell due to remelting */
            /**todo: is solute conserved here? */
            if (*(ofs_poly[0]) < 0.) {
              *ofs_poly[0] = 0.;
              *ngr = 0;
              *ogr = 0;
              *nat_cell_p = 0;
              *onat_cell_p = 0;
              *nat_grain_p = 0;
              *onat_grain_p = 0;

              for (isol = 0; isol < ele_1; isol++) {    /*loop on solutes */
                *ncl_poly[isol] = *nce_poly[isol];
              }
              sp->ncsolid--;    /* number of solid cells decreases by one */
            }

            *ofs_poly[1] = 0.;

            *ofs = *ofs_poly[0];
            sum_del_fs += *ofs - fs;    /* increase of the fraction solid */

          } else {              /*ofs_poly[0]=1 */
          /* cell is already solid */
            for (isol = 0; isol < ele_1; isol++) {
              *ncl_poly[isol] = *nce_poly[isol] / bp->mprops.alloyprops[isol].part_coef[0];
              /* C_E = C_S, when fs = 1 */
            }
          }
        }/* end *nat_cell_p == 1 case */

        if (*nat_cell_p == 3) {
          /* CASE 3 : liquid + the secondary eutectic phase */
          if (*(ofs_poly[1]) < 1) {
            bp->dim_tab = ele_1;
            for (i1 = 0; i1 < NPHAMAX; i1++) {
              for (i2 = 0; i2 < NPHAMAX; i2++) {
                bp->a[i1][i2] = 0.0;
              }
            }
            /*calculations of the cofficients */
            coef_curv = bp->mprops.gibbs_thomson / bp->size_c[0];
            /*mass conservation equation for each solute */
            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes */
              dce_poly[isol] = *nce_poly[isol] - *oce_poly[isol];
              bp->a[isol][isol] = 1.0;

              if (bp->mprops.stoechio == FALSE){
                bp->a[isol][isol] -= (1.0 - bp->mprops.alloyprops[isol].part_coef[1]) * (*ofs_poly[1]);
                bp->a[isol][ele_1] = -(1.0 - bp->mprops.alloyprops[isol].part_coef[1]) * (*ncl_poly[isol]);
              
              }else{
                 bp->a[isol][isol]-=*ofs_poly[1];
                 bp->a[isol][ele_1]=-(*ncl_poly[isol]-bp->mprops.alloyprops[isol].cs_stoechio);
              }

              bp->a[isol][ele_1+1]=dce_poly[isol];

            }
            /* Thermodynamic equilibria */
            /* first case : constant thermodynamic data */
            oTliq = bp->mprops.tp_poly[1];
            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes */
              oTliq += bp->mprops.alloyprops[isol].m_solute[1] * (*ncl_poly[isol]);
            }
            /* second case : coupling with Thermo-Calc */
            /* end liquidus temperature calculation */
            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes */
              bp->a[ele_1][isol] = bp->mprops.alloyprops[isol].m_solute[1];
            }
            bp->a[ele_1][ele_1 + 1] = Tcell + coef_curv * (*ncv) - oTliq;
            invsys (bp, sbnum);
            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes */
              *ncl_poly[isol] += bp->a[isol][ele_1 + 1];
            }
            *ofs_poly[1] += bp->a[ele_1][ele_1 + 1];
            *ofs_poly[1] = MIN (1.0, *ofs_poly[1]);

            /*test to know if the equivalent composition lies in the tie triangle  */
            oTliq0=bp->mprops.tp_poly[0];
            for (isol = 0; isol < ele_1; isol++) { /*loop on solutes */
              oTliq0 +=bp->mprops.alloyprops[isol].m_solute[0]*(*ncl_poly[isol]);
            }

            if (Tcell < oTliq0) {
              *nat_cell_p=2;
              *onat_cell_p=2;
            }

            /*end test*/


            if (*(ofs_poly[1]) < 0.) {
              *ofs_poly[1] = 0.;
              *ngr = 0;
              *ogr = 0;
              *nat_cell_p = 0;
              *onat_cell_p = 0;
              *nat_grain_p = 0;
              *onat_grain_p = 0;

              for (isol = 0; isol < ele_1; isol++) {    /*loop on solutes */
                *ncl_poly[isol] = *nce_poly[isol];
              }

              sp->ncsolid--;    /* number of solid cells decreases by one */
            }
            *ofs_poly[0] = 0.;
            *ofs = *ofs_poly[1];
            sum_del_fs += *ofs - fs;    /* increase of the fraction solid */
          } else {              /*ofs_poly[1]=1 */

            for (isol = 0; isol < ele_1; isol++) {
                 /**ncl_poly[isol] = *nce_poly[isol]/bp->mprops.part_coef_poly[1][isol];*/
              *ncl_poly[isol] = *nce_poly[isol];
              /* C_E = C_S, when fs = 1 */
            }
          }
        }
        /****************************/
        /* end *nat_cell_p ==3 case */
        /****************************/


        /* CASE 2 : liquid + eutectique phase */
        if (*nat_cell_p == 2) {
          if (*ofs_poly[0] + *ofs_poly[1] < 1.) {

/* calculation of the coefficients */

            /* initialize the matrix of coefficients */
            for (i1 = 0; i1 < NPHAMAX; i1++) {

              for (i2 = 0; i2 < NPHAMAX; i2++) {

                bp->a[i1][i2] = 0.0;

              }
            }

            bp->iphs_eff = 2;
            bp->ipheq_eff = 2;
            bp->dim_tab = ele_1 + 1;

            coef_curv = bp->mprops.gibbs_thomson / bp->size_c[0];

            /*mass conservation equation for each solute */

            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes */
              dce_poly[isol] = *nce_poly[isol] - *oce_poly[isol];
              bp->a[isol][isol] = 1.0;
              bp->a[isol][isol] -= (1.0 - bp->mprops.alloyprops[isol].part_coef[0]) * (*ofs_poly[0]);

              if (bp->mprops.stoechio == FALSE){
                bp->a[isol][isol] -= (1.0 - bp->mprops.alloyprops[isol].part_coef[1]) * (*ofs_poly[1]);
                bp->a[isol][ele_1 + 1] = -(1.0 - bp->mprops.alloyprops[isol].part_coef[1]) * (*ncl_poly[isol]);

              }else{
                bp->a[isol][isol]-=*ofs_poly[1];
                bp->a[isol][ele_1+1]=-(*ncl_poly[isol] - bp->mprops.alloyprops[isol].cs_stoechio);
              }


              bp->a[isol][ele_1] = -(1.0 - bp->mprops.alloyprops[isol].part_coef[0]) * (*ncl_poly[isol]);

              bp->a[isol][ele_1 + 2] = dce_poly[isol];
            }

            /* Thermodynamic equilibria */

            oTliq = bp->mprops.tp_poly[1];
            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes */
              oTliq += bp->mprops.alloyprops[isol].m_solute[1] * (*ncl_poly[isol]);
            }

            delta_T = Tcell + coef_curv * (*ncv) - oTliq;

            for (ieq = 0; ieq < bp->ipheq_eff; ieq++) {
              /* second case : coupling with Thermo-Calc */
              /* end liquidus temperature calculation */
              for (isol = 0; isol < ele_1; isol++) {    /*loop on solutes */
                bp->a[ele_1 + ieq][isol] = bp->mprops.alloyprops[isol].m_solute[ieq];
              }
              bp->a[ele_1 + ieq][ele_1 + bp->iphs_eff] = delta_T;

            }

/*End calculation of the coefficients */

            invsys (bp, sbnum);

            for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes */
              *ncl_poly[isol] += bp->a[isol][ele_1 + bp->iphs_eff];
            }

            for (iphs = 0; iphs < bp->iphs_eff; iphs++) {
              *ofs_poly[iphs] += bp->a[ele_1 + iphs][ele_1 + bp->iphs_eff];
            }

            if (*ofs_poly[0] > 1.) {
              *ofs_poly[0] = 1.;
              *ofs_poly[1] = 0.;
              *nat_cell_p = 1;
              *onat_cell_p = 1;

              for (isol = 0; isol < ele_1; isol++) {
                *ncl_poly[isol] = *nce_poly[isol] / bp->mprops.alloyprops[isol].part_coef[0];
                /* C_E = C_S, when fs = 1 */
              }

            }

            if (*ofs_poly[1] > 1.) {
              *ofs_poly[1] = 1.;
              *ofs_poly[0] = 0.;
              /**todo: why does state change to 3 here ? */
              *nat_cell_p = 3;
              *onat_cell_p = 3;

              for (isol = 0; isol < ele_1; isol++) {
                /* *ncl_poly[isol] = *nce_poly[isol] / bp->mprops.alloyprops[isol].part_coef[1];*/
                *ncl_poly[isol]=*nce_poly[isol];

                /* C_E = C_S, when fs = 1 */
              }

              /*fprintf(stderr,"transition non geree\n"); */

            }

            if ((*ofs_poly[0] < 0.) && (*ofs_poly[1] < 1.) && (*ofs_poly[1] > 0.)) {
              *ofs_poly[0] = 0.;
              *nat_cell_p = 3;
              *onat_cell_p = 3;

              for (isol = 0; isol < ele_1; isol++) {
                /* *ncl_poly[isol] = *nce_poly[isol] / (1.0 - (1.0 - bp->mprops.alloyprops[isol].part_coef[1]) * (*ofs_poly[1]));*/

                *ncl_poly[isol]=(*nce_poly[isol]-(*ofs_poly[1])*bp->mprops.alloyprops[isol].cs_stoechio)/(1.-(*ofs_poly[1]));

                /* C_E = C_S, when fs = 1 */
              }

              /*fprintf(stderr,"transition non geree\n"); */

            }

            if ((*ofs_poly[1] < 0.) && (*ofs_poly[0] < 1.) && (*ofs_poly[0] > 0.)) {
              *ofs_poly[1] = 0.;
              *nat_cell_p = 1;
              *onat_cell_p = 1;
              for (isol = 0; isol < ele_1; isol++) {
                *ncl_poly[isol] = *nce_poly[isol] / (1.0 - (1.0 - bp->mprops.alloyprops[isol].part_coef[0]) * (*ofs_poly[1]));
                /* C_E = C_S, when fs = 1 */
              }
            }

            if ((*ofs_poly[0] < 0.) && (*ofs_poly[1] < 0.)) {   /*the cell is fully liquid */
              *ofs_poly[0] = 0.;
              *ofs_poly[1] = 0.;
              *ngr = 0;
              *ogr = 0;
              *nat_cell_p = 0;
              *onat_cell_p = 0;
              *nat_grain_p = 0;
              *onat_grain_p = 0;


              for (isol = 0; isol < ele_1; isol++) {
                *ncl_poly[isol] = *nce_poly[isol];
                /* C_E = C_L, when fs = 0 */
              }

              sp->ncsolid--;    /* number of solid cells decreases by one */

            }

            if ((*ofs_poly[0] + *ofs_poly[1] > 1.) && (*ofs_poly[0] < 1.) && (*ofs_poly[0] > 0.) && (*ofs_poly[1] < 1.)
                && (*ofs_poly[1] > 0.)) {

              *ofs_poly[0] = 1.0 - *ofs_poly[1];

              for (isol = 0; isol < ele_1; isol++) {
                /* *ncl_poly[isol] =
                  *nce_poly[isol] / (bp->mprops.alloyprops[isol].part_coef[0] * (*ofs_poly[0]) +
                                     bp->mprops.alloyprops[isol].part_coef[1] * (*ofs_poly[1]));*/

                *ncl_poly[isol] = *nce_poly[isol];

                /* C_E = C_S, when fs = 1 */
              }
            }

            *ofs = *ofs_poly[0] + *ofs_poly[1];
            sum_del_fs += *ofs - fs;

          } else {              /*(*ofs_poly[0]+*ofs_poly[1])=1 */

            for (isol = 0; isol < ele_1; isol++) {
              /* *ncl_poly[isol] =
                *nce_poly[isol] / (bp->mprops.alloyprops[isol].part_coef[0] * (*ofs_poly[0]) +
                                   bp->mprops.alloyprops[isol].part_coef[1] * (*ofs_poly[1]));*/

              *ncl_poly[isol] = *nce_poly[isol];

              /* C_E = C_S, when fs = 1 */
            }
          }

        } /* end *nat_cell_p == 2 case */

        *ofs = *ofs_poly[0] + *ofs_poly[1];

        /* check for some error conditions */
        if ((*ofs == 0.0) && (*ngr > 0)) {
          fprintf (stderr, "probleme:%s: ((*ofs == 0.0) && (*ngr > 0))\n",__func__);
        }

        if (*nat_cell_p > 3) {
          fprintf (stderr, "probleme:%s: *nat_cell_p > 3\n",__func__);
        }

        /* end of the real work, now */
        /* increment the pointers */
        nct++;
        ngr++;
        nat_cell_p++;
        onat_cell_p++;
        nat_grain_p++;
        onat_grain_p++;

        nfs++;
        for (iphs = 0; iphs < iphs_tot; iphs++) {
          nfs_poly[iphs]++;
        }

        for (isol = 0; isol < ele_1; isol++) {
          ncl_poly[isol]++;
          nce_poly[isol]++;
          oce_poly[isol]++;
        }

        ncv++;
        ofs++;
        for (iphs = 0; iphs < iphs_tot; iphs++) {
          ofs_poly[iphs]++;
        }
        ogr++;

      }                         /*end of the loop on x */

      ofs += 2;
      for (iphs = 0; iphs < iphs_tot; iphs++) {
        ofs_poly[iphs] += 2;
      }

      for (isol = 0; isol < ele_1; isol++) {
        oce_poly[isol] += 2;
      }
      ogr += 2;
      onat_cell_p += 2;
      onat_grain_p += 2;

    }                           /*end of the loop on y */

    ofs += skip;
    for (iphs = 0; iphs < iphs_tot; iphs++) {
      ofs_poly[iphs] += skip;
    }

    for (isol = 0; isol < ele_1; isol++) {
      oce_poly[isol] += skip;
    }

    ogr += skip;
    onat_cell_p += skip;
    onat_grain_p += skip;

  }                             /*end of the loop on z */

  /* added in by Wei WANG on 06-11-03 */
  sp->Tvals.del_fs = sum_del_fs / bp->ncsb;     /* increase of average fs in SubBlock */
  sp->Tvals.fsavg += sp->Tvals.del_fs;  /* average fs in SubBlock */

  return (1);

}/* END of fs_change_diffuse_poly */

/*************************************************/
/* Inversion of the system L. THUINET 07-02-05   */
/* Solution is written in the column a(i,n+1)    */
/* if the system is not reversible, noninv=1     */
/*************************************************/

void invsys (BB_struct * bp, int sbnum)
{
  int noninv, i_sys, pospiv, l_sys, k_sys, j_sys;
  CA_FLOAT pivot, temp, coef, somme;

  noninv = 0;
  i_sys = 0;

  while ((noninv == 0) && (i_sys < bp->dim_tab)) {
    pivot = fabs (bp->a[i_sys][i_sys]);
    pospiv = i_sys;

    for (l_sys = i_sys + 1; l_sys <= bp->dim_tab; l_sys++) {

      if (fabs (bp->a[l_sys][i_sys]) > pivot) {
        pivot = fabs (bp->a[l_sys][i_sys]);
        pospiv = l_sys;
      }
    }

    if (pivot > 0.) {
      if (pospiv != i_sys) {
        for (k_sys = i_sys; k_sys <= bp->dim_tab + 1; k_sys++) {
          temp = bp->a[i_sys][k_sys];
          bp->a[i_sys][k_sys] = bp->a[pospiv][k_sys];
          bp->a[pospiv][k_sys] = temp;
        }
      }

      for (j_sys = i_sys + 1; j_sys <= bp->dim_tab; j_sys++) {
        coef = bp->a[j_sys][i_sys] / bp->a[i_sys][i_sys];

        for (k_sys = i_sys + 1; k_sys <= bp->dim_tab + 1; k_sys++) {
          bp->a[j_sys][k_sys] = bp->a[j_sys][k_sys] - bp->a[i_sys][k_sys] * coef;
        }
      }
    } else {
      noninv = 1;
    }
    i_sys = i_sys + 1;
  }

  if (bp->a[bp->dim_tab][bp->dim_tab] != 0) {
    bp->a[bp->dim_tab][bp->dim_tab + 1] = bp->a[bp->dim_tab][bp->dim_tab + 1] / bp->a[bp->dim_tab][bp->dim_tab];
  } else {
    noninv = 1;
  }

  if (noninv == 0) {

    for (i_sys = bp->dim_tab - 1; i_sys >= 0; i_sys--) {
      somme = 0.;

      for (j_sys = i_sys + 1; j_sys <= bp->dim_tab; j_sys++) {
        somme = somme + bp->a[i_sys][j_sys] * (bp->a[j_sys][bp->dim_tab + 1]);
      }

      bp->a[i_sys][bp->dim_tab + 1] = (bp->a[i_sys][bp->dim_tab + 1] - somme) / bp->a[i_sys][i_sys];
    }

  }

}


/************************************************************/
/* Growth of square/octahedron for polycomponent systems    */
/************************************************************/
int grow_octahedron_poly (BB_struct * bp, int sbnum)
{
  SB_struct *sp;
  int i, j, k;                  /* tmp counter */

  int nx, ny, nz, skip;
  int *ngr;
  int *onat_grain_p;
  CA_FLOAT *npd, *npx, *npy, *npz;
  CA_FLOAT *ofs;
  CA_FLOAT g00, g10, g20, g01, g11, g21, g02, g12, g22;
  CA_FLOAT rx0, ry0, rz0, gx0, gy0, gz0, ds1;
  CA_FLOAT fs;

  sp = bp->sb[sbnum];

  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  skip = 2 * (nx + 2);

  /* Set up local pointers */
  ngr = sp->gr;                 /* array of old cell grain, size (nc[i]+2)^3 */

  npd = sp->dc_d;               /* array of decentred octahedron half diagonal */
  npx = sp->dc_x;               /* array of decentred octahedron centre */
  npy = sp->dc_y;
  npz = sp->dc_z;

  ofs = bp->ftmp_one;
  ofs += bp->cubeptr.flist[0][START];

  onat_grain_p = bp->itmp_nat_grain;
  onat_grain_p += bp->cubeptr.flist[0][START];

   /***********************************************/
  /* Beginning of main loop(s)                   */
   /***********************************************/

  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {

        if (*ofs <= 0) {

          /* cell grain number ==0, the cell is LIQUID */

        } else if ((*ofs > 0) && (*onat_grain_p == 1)) { /* the cell is captured by a grain */
          fs = *ofs;

#ifdef JUNK
          ca0 = bp->gr[*ngr]->cang[0];
          sa0 = bp->gr[*ngr]->sang[0];
          ca1 = bp->gr[*ngr]->cang[1];
          sa1 = bp->gr[*ngr]->sang[1];
          ca2 = bp->gr[*ngr]->cang[2];
          sa2 = bp->gr[*ngr]->sang[2];
#endif

          g00 = bp->gr[*ngr]->g[0][0];
          g10 = bp->gr[*ngr]->g[1][0];
          g20 = bp->gr[*ngr]->g[2][0];

          g01 = bp->gr[*ngr]->g[0][1];
          g11 = bp->gr[*ngr]->g[1][1];
          g21 = bp->gr[*ngr]->g[2][1];

          g02 = bp->gr[*ngr]->g[0][2];
          g12 = bp->gr[*ngr]->g[1][2];
          g22 = bp->gr[*ngr]->g[2][2];

          rx0 = 0. - *npx;
          ry0 = 0. - *npy;
          rz0 = 0. - *npz;

          gx0 = g00 * rx0 + g10 * ry0 + g20 * rz0;
          gy0 = g01 * rx0 + g11 * ry0 + g21 * rz0;
          gz0 = g02 * rx0 + g12 * ry0 + g22 * rz0;

          /* ds1 = ABS(gx0) + ABS(gy0); 2D */

          ds1 = ABS (gx0) + ABS (gy0) + ABS (gz0);      /*3D */

          /*
             if(fs > 1.- ds2) {
             *npd = ds1 + fs * (ca0+ABS(sa0));
             *npd = ds1 + fs * 1.414214 ;
             } else {
             *npd = ds1 + 0.707107 * (SQRT(ds2*ds2 + 4.*fs) - ds2);
             }
           */
          *npd = ds1 + fs * 1.414214;
          *ofs = MIN (*ofs, 1.0);
          /* *npd = ds1 + fs * (ca0 + ABS(sa0));  */
          /* *npd = ds1 + fs * (0.195262*(ca0+ABS(sa0)) + 1.138072);        */
          /*         *npd = ds1 + fs * (0.4*(ca0+ABS(sa0)) + 0.848528);     */
          /*          *npd = ds1 + fs*1.414214*1.1 ; */
        }

        ngr++;                  /* ptr to outside block cell grain number */
        npd++;                  /* ptr to decentred octahedron half diagonal */
        npx++;                  /* ptr to decentred octahedron centre */
        npy++;
        npz++;
        ofs++;
        onat_grain_p++;
      }                         /*end of I loop */
      ofs += 2;
      onat_grain_p += 2;
    }                           /* end of J loop */
    ofs += skip;
    onat_grain_p += skip;
  }                             /* end of K loop */

  return (1);
}


/************************************************************************************/
/* extension of capture new cell with diffusion in multicomponent by THUINET        */
/* on 17-02-05                                                                      */
/************************************************************************************/
int capture_octahedron_diffuse_poly (BB_struct * bp, int sbnum)
{
  SB_struct *sp;
  int i, j, k, n;               /* tmp counter */
  int i_nn, j_nn, k_nn, ijk_nn, gr_nn, capt;
  int nx, ny, nz, skip;
  int *ngr, *ogr;
  int *nid;
  CA_FLOAT *opd, *opx, *opy, *opz, *npd, *npx, *npy, *npz;
  CA_FLOAT *ofs, *nfs;
  CA_FLOAT *ncl, *nce;

/* THUINET 07-02-05 */
  int isol;
  int ele_num, ele_1;
  int iphs, iphs_tot;
  int noninv;
  CA_FLOAT *nce_poly[NSOLMAX], *ncl_poly[NSOLMAX];
  Ctrl_str *cp;
  int *nat_cell_p, *onat_cell_p;
  int nat_cell_nn;
  int *nat_grain_p, *onat_grain_p;
  int nat_grain_nn;

/* End THUINET 07-02-05*/
  CA_FLOAT *nct;
  CA_FLOAT g00, g10, g20, g01, g11, g21, g02, g12, g22;
  CA_FLOAT rx, ry, rz, dx, dy, dz, ds1;
  CA_FLOAT ncd, dd, ncx, ncy, ncz;
  CA_FLOAT part_coef, kInv, km;
  CA_FLOAT fs;

  /* nearest neighbours structure */
  struct nns {
    int nx;
    int ny;
    int nz;
  } nn[] = { 1, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 1, 1, 0, -1, 1, 0, -1, -1, 0, 1, -1, 0};
  int NNB = 6;
  int nbskip[6];

  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  skip = 2 * (nx + 2);
  for (n = 0; n < NNB; n++) {
    i_nn = nn[n].nx;
    j_nn = nn[n].ny;
    k_nn = nn[n].nz;
    nbskip[n] = i_nn + j_nn * (nx + 2) + k_nn * (nx + 2) * (ny + 2);
  }
  sp = bp->sb[sbnum];
  cp = bp->ctrl;
  ele_num = cp->NUM_COMP;       /* number of elements in the alloy */
  ele_1 = ele_num - 1;

  iphs_tot = cp->NUM_PHS;

  nat_cell_p = sp->nat_cell;
  onat_cell_p = bp->itmp_nat_cell;
  onat_cell_p += bp->cubeptr.flist[0][START];

  nat_grain_p = sp->nat_grain;
  onat_grain_p = bp->itmp_nat_grain;
  onat_grain_p += bp->cubeptr.flist[0][START];

/*   part_coef = bp->mprops.part_coef_alloy;*/
/*   km = 1-part_coef;*/
  /* Set up local pointers */
  nid = sp->index;
  ogr = bp->itmp_one;           /* array of old cell grain, size (nc[i]+2)^3 */
  ngr = sp->gr;                 /* array of new cell grain, size (nc[i])^3   */

  ofs = bp->ftmp_one;
  opd = bp->ftmp_dc_d;          /* array of decentred octahedron half diagonal */
  opx = bp->ftmp_dc_x;          /* array of decentred octahedron centre */
  opy = bp->ftmp_dc_y;
  opz = bp->ftmp_dc_z;
  npd = sp->dc_d;
  npx = sp->dc_x;
  npy = sp->dc_y;
  npz = sp->dc_z;
  nfs = sp->c_fs;
/*   ncl = sp->c_sol_alloy;*/
/*   nce = sp->c_eqv_alloy;*/
  for (isol = 0; isol < ele_1; isol++) {
    ncl_poly[isol] = sp->c_sol_poly[isol];
    nce_poly[isol] = sp->c_eqv_poly[isol];
  }
  nct = sp->c_temp;
  /*copy the half diagonal and coodinates of octahedron to temporary buffers */
  fcopy_matrix (PAD, opd, npd, bp, NULL, sbnum);
  fcopy_matrix (PAD, opx, npx, bp, NULL, sbnum);
  fcopy_matrix (PAD, opy, npy, bp, NULL, sbnum);
  fcopy_matrix (PAD, opz, npz, bp, NULL, sbnum);
  /*move the ptr to outside corner to inside corner */
  ogr += bp->cubeptr.flist[0][START];
  ofs += bp->cubeptr.flist[0][START];
  opd += bp->cubeptr.flist[0][START];
  opx += bp->cubeptr.flist[0][START];
  opy += bp->cubeptr.flist[0][START];
  opz += bp->cubeptr.flist[0][START];
   /***********************************************/
  /* Beginning of main loop(s)                   */
   /***********************************************/
  /* WEI_DC, a cell is captured when any part of it is trapped by a growing square */
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
        if (*nid == 1) {        /* the cell is LIQUID,and adjacent to s/l interface */
          capt = FALSE;

          for (n = 0; n < NNB && capt == FALSE; n++) {  /* neighbour loop */
            i_nn = nn[n].nx;
            j_nn = nn[n].ny;
            k_nn = nn[n].nz;
            ijk_nn = nbskip[n];

            gr_nn = *(ogr + ijk_nn);
            if (gr_nn > 0) {    /* the neighbour has a s/l interface */
              nat_grain_nn = *(onat_grain_p + ijk_nn);
              nat_cell_nn = *(onat_cell_p + ijk_nn);

              /* nature of the neighbour grain == 1 then it is a primary phase grain */
              if (nat_grain_nn == 1) {   /*application of the decentred_square algorithm */

                g00 = bp->gr[gr_nn]->g[0][0];
                g10 = bp->gr[gr_nn]->g[1][0];
                g20 = bp->gr[gr_nn]->g[2][0];
                g01 = bp->gr[gr_nn]->g[0][1];
                g11 = bp->gr[gr_nn]->g[1][1];
                g21 = bp->gr[gr_nn]->g[2][1];
                g02 = bp->gr[gr_nn]->g[0][2];
                g12 = bp->gr[gr_nn]->g[1][2];
                g22 = bp->gr[gr_nn]->g[2][2];
                rx = -1. * i_nn - *(opx + ijk_nn);
                ry = -1. * j_nn - *(opy + ijk_nn);
                rz = -1. * k_nn - *(opz + ijk_nn);
                dx = g00 * rx + g10 * ry + g20 * rz;
                dy = g01 * rx + g11 * ry + g21 * rz;
                dz = g02 * rx + g12 * ry + g22 * rz;
                /* ds1 =  ABS(dx) + ABS(dy); *//*2D */
                ds1 = ABS (dx) + ABS (dy) + ABS (dz);   /*3D */
                /*ds2 =  2.0 * MIN(ABS(dx), ABS(dy)); */
                /*ds2 = SQRT(dx*dx + dy*dy) + 0. * (ca0 + ABS(sa0) - 1.0);  *//* round the tip */
                /* fs = *(opd+ijk_nn) - MAX(ds1, ds2); */
                fs = *(opd + ijk_nn) - ds1;
                fs += 0.05 * drand48 ();        /* add some numerical noise */
                if (fs > 0.) {  /* probalistic growth */
                  /* if(ABS(dx) + ABS(dy) <= *(opd+ijk_nn)) *//* only for 2D */
                  /* if (*(ofs + ijk_nn) == 1.0) */
                  capt = TRUE;
                  ncd = MIN (*(opd + ijk_nn), 1.);
                  dd = *(opd + ijk_nn) - ncd;
                  /* determine which corner the newly captured cell is closest to */
                  if (dx + dy >= 0 && dx - dy > 0 && dx + dz > 0 && dx - dz >= 0) {
                    ncx = dd;
                    ncy = 0.;
                    ncz = 0.;
                  } else if (dx + dy <= 0 && dx - dy < 0 && dx + dz < 0 && dx - dz <= 0) {
                    ncx = -1. * dd;
                    ncy = 0.;
                    ncz = 0.;
                  } else if (dy + dx > 0 && dy - dx >= 0 && dy + dz >= 0 && dy - dz > 0) {
                    ncx = 0.;
                    ncy = dd;
                    ncz = 0.;
                  } else if (dy + dx < 0 && dy - dx <= 0 && dy + dz <= 0 && dy - dz < 0) {
                    ncx = 0.;
                    ncy = -1. * dd;
                    ncz = 0.;
                  } else if (dz + dx >= 0 && dz - dx > 0 && dz + dy > 0 && dz - dy >= 0) {
                    ncx = 0.;
                    ncy = 0.;
                    ncz = dd;
                  } else if (dz + dx <= 0 && dz - dx < 0 && dz + dy < 0 && dz - dy <= 0) {
                    ncx = 0.;
                    ncy = 0.;
                    ncz = -1. * dd;
                  } else {
                    ncx = 0.;
                    ncy = 0.;
                    ncz = 0.;
                  }
                  sp->ncsolid++;        /* number of solid cells in sb ++ */
                  *ngr = *(ogr + ijk_nn);
                  *nat_grain_p = nat_grain_nn;
                  *nat_cell_p = nat_cell_nn;

                  /* Robert added -- keep track of grain size and location */
                  add_to_grain (bp->gr[*ngr], i, j, k);
                    /*********************************************************/
                  *npd = ncd;
                  *npx = *(opx + ijk_nn) + i_nn + g00 * ncx + g01 * ncy + g02 * ncz;
                  *npy = *(opy + ijk_nn) + j_nn + g10 * ncx + g11 * ncy + g12 * ncz;
                  *npz = *(opz + ijk_nn) + k_nn + g20 * ncx + g21 * ncy + g22 * ncz;
                  fs = 0.;
                  /* fs = MIN(0.01, MAX(0.,fs)); */
                  *ofs = fs;
                  /* *ncl = (*nce - part_coef * c_int * fs) / (1-fs); */
                  for (isol = 0; isol < ele_1; isol++) {
                    *ncl_poly[isol] = *nce_poly[isol];
                  }
                  /* *ncl = *nce / (1-km*fs); */

                }               /*end of if fs is positive */
              } /*end of if nat_grain_nn==1 */

              /* for eutectic phases , do not use the decentered square */
              else {            /*application of a CA_rule for the growth of both eutectic phases */
                if (*(ofs + ijk_nn) >= 0.85) {  /*arbitrary value */
                  capt = TRUE;
                  sp->ncsolid++;        /* number of solid cells in sb ++ */
                  *ngr = *(ogr + ijk_nn);
                  *nat_grain_p = nat_grain_nn;
                  *nat_cell_p = nat_cell_nn;
                  add_to_grain (bp->gr[*ngr], i, j, k);
                  *ofs = 0.;
                  for (isol = 0; isol < ele_1; isol++) {
                    *ncl_poly[isol] = *nce_poly[isol];
                  }
                }               /*end of if fs_nn>=0.75 */
              }                 /*end of if nat_grain_nn==1 */
            }                   /*end of gr_nn pos */
          }                     /*end of loop on the neighbour */
        }
        /* end of nid equal to 1 */
        nid++;
        nct++;
        ogr++;                  /* ptr to outside block cell grain number */
        ngr++;                  /* ptr to inside  block cell grain number */

        onat_cell_p++;
        nat_cell_p++;

        onat_grain_p++;
        nat_grain_p++;

        npd++;
        npx++;
        npy++;
        npz++;
        opd++;                  /* ptr to decentred octahedron half diagonal */
        opx++;                  /* ptr to decentred octahedron centre */
        opy++;
        opz++;
        ofs++;
        nfs++;
        /*ncl++; */
        /*nce++; */
        for (isol = 0; isol < ele_1; isol++) {
          ncl_poly[isol]++;
          nce_poly[isol]++;
        }

      }                         /*end of I loop */

      onat_cell_p += 2;
      onat_grain_p += 2;
      ogr += 2;
      opd += 2;
      opx += 2;
      opy += 2;
      opz += 2;
      ofs += 2;
    }                           /* end of J loop */

    onat_cell_p += skip;
    onat_grain_p += skip;
    ogr += skip;

    opd += skip;
    opx += skip;
    opy += skip;
    opz += skip;
    ofs += skip;
  }                             /* end of K loop */

  if (sp->ncsolid >= bp->ncsb) {
    fprintf (stderr, "WARNING! sb_decentred_step.c-capture_octahedron_diffuse(): SB#%d completely solid.\n", sbnum);
    /* sp->done = TRUE; *//* only give a warning when all cells are fully/partially solid. Wei WANG */
  }
  return (1);
}

/*
*/

/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *sb_decentred_step_poly_c ()
{
  static char const rcsid[] = "$Id: sb_decentred_step_poly.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}
