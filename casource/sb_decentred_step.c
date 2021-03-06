
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
/* This code is part of the umats routines developed at in the  */
/* Materials Processing Group, Dept. of Materials, ICSTM.       */
/*      email p.d.lee or r.atwood @ic.ac.uk for details         */
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/****************************************************************/
/*      sb_decentred_step.c:                                           */
/*                                                              */
/* The main subroutines of decentred square/octahedron algorithm */
/* This file is added by Wei WANG on 11-07-02                   */
/*                                                              */
/*                                                              */
/*                                                              */
/****************************************************************/
/*RCS ID: $Id: sb_decentred_step.c 1356 2008-08-18 13:41:15Z  $*/
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

/***************************************************/
/* stochastic cell nucleation                      */
/* by Wei WANG 03-10-02                            */
/***************************************************/
int cell_nucleation (BB_struct * bp, int sbnum)
{
  int errors = 0;
  int global_undercooling = 0;

#ifdef OLD_TUNDER
  CA_FLOAT *otu, *otup;
#endif /*OLD_TUNDER */
  int i, j, k;
  int nx, ny, nz, skip;
  int fixed_nuc, i_nuc, icell, index;
  SB_struct *sp;
  Ctrl_str *cp = bp->ctrl;
  CA_FLOAT Tliq, Tunder, TliqFixed, min_temp;
  CA_FLOAT a0, a1, a2, liq_slope, c_ini;
  CA_FLOAT *thrp, *c_temp_p, *ncl, *nce;
  int *ngr, *ogr;
  int phase_diag_on, cell_nuc_on, cell_temp_on, n_neigh;
  int **gv_p = bp->gr_array;    /* grain block value array pointer pointer :-/ */
  static int nucmsg = 0;
  int (*cell_nuc_func) (BB_struct * bp, CA_FLOAT Tunder, CA_FLOAT argthree);

  /*define variables for decentred algorithm */
  /*by Wei Wang on 09-01-01 */
  sp = bp->sb[sbnum];
  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  skip = 2 * (nx + 2);

  c_temp_p = sp->c_temp;

  c_ini = bp->mprops.alloyprops[0].Cinit;
  ncl = sp->c_sol_alloy;
  nce = sp->c_eqv_alloy;
  liq_slope = bp->mprops.alloyprops[0].m_solute[0];
  /* todo needs to work for binary and poly */

  TliqFixed = bp->mprops.Tliq;
  Tliq = TliqFixed;

  /* set up phase diagram mode */
  phase_diag_on = bp->ctrl->phase_diag_on;
  global_undercooling = (phase_diag_on && bp->ctrl->global_undercooling);

  /* Set up local pointers */
#ifdef OLD_TUNDER
  otup = otu = bp->old_Tunder;
#endif /*OLD_TUNDER */

  ogr = bp->itmp_one;           /* array of old cell grain, size (nc[i]+2)^3 */
  ngr = sp->gr;                 /* array of new cell grain, size (nc[i])^3   */

  /* copy grain number to temporary buffer */
  icopy_matrix (PAD, ogr, ngr, bp, gv_p, sbnum);

  /*move the ptr to outside corner to inside corner */
  ogr += bp->cubeptr.flist[0][START];

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
/* TODO: Get liquidus for each cell seperately , not combined with undercooling calculation */
if (bp->ctrl->diffuse_alloy){
        Tunder = Tliq - *c_temp_p + liq_slope * (*ncl - c_ini); /* modified by Dong */
}else{
        Tunder = Tliq - *c_temp_p;   /*undercooling = Tliq-Tcell */
}

#         ifdef OLD_TUNDER
        if (bp->step <= 1)
          *otup = Tunder;
#         endif /*OLD_TUNDER */

        if (*ngr == LIQUID) {   /* cell grain number =0 */
          if (cell_nuc_on) {
#             ifdef GLOBAL_UND_NUC
            /*for global undercooling nucleation model */
            if (global_undercooling)
              Tunder = TliqFixed - *c_temp_p;
#             endif

#             ifdef OLD_TUNDER  /* only used when trying to calculate threshold by integrating */
            /* the undercooling from the previous step ..  */
            if (bp->ctrl->block_nuc) {
              fprintf (stderr, "You can't do that! You have to debug the code.\n");
              fprintf (stderr, "You can't use OLD_TUNDER and block_nuc at the same time\n.");
              exit (1);
            };
            /* then cell_nuc_func needs to have the old Tunder pointer otup defined */
            if (Tunder > 0 && (*cell_nuc_func) (bp, Tunder, *otup) == 1)
#             else /*OLD_TUNDER */
            if (Tunder > 0 && (*cell_nuc_func) (bp, Tunder, *thrp) == 1)
#             endif /*OLD_TUNDER */
            {                   /* trheshold is exceeded so nucleate a grain */
              /* error checks */
              if (nucmsg < MAX_NUC_MSG) {
                fprintf (stderr, "cell temp[%g] cell Tliq[%g]\n", *c_temp_p, Tliq);
                nucmsg++;
              }

              if (bp->nprops.ngr == bp->nprops.gd_max_total) {
                static int emess=0;
                if (emess <= MAX_ERR_MSG){
                   fprintf (stderr,
                         "ERROR:sb_decentred_step.c-cell_nucleation(): Max # grains set by user [%d] exceeded. Increase option MaxTotGrains.\n",
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
                (bp->nprops.ngr)++;
                (bp->sb[sbnum]->ngr)++;
                *ngr = bp->nprops.ngr;

                icell = i + j * nx + k * nx * ny;
                fixed_nuc = 0;
                for (i_nuc = 0; i_nuc < cp->nfnuc; i_nuc++) {
                  index = (int) (cp->nsite[i_nuc][0] + cp->nsite[i_nuc][1] * nx + cp->nsite[i_nuc][2] * nx * ny);
                  if (icell == index) {
                    a0 = PI_BY_FOUR / 45.0 * cp->nsite[i_nuc][3];
                    a1 = PI_BY_FOUR / 45.0 * cp->nsite[i_nuc][4];
                    a2 = PI_BY_FOUR / 45.0 * cp->nsite[i_nuc][5];
                    init_fixed_angle_grain (bp, bp->nprops.ngr, sbnum, index, 1, a0, a1, a2);   /*by Wei WANG 30-07-02 */
                    fprintf (stderr, "sb_decentred_step.c-cell_nucleation(): nucleating fixed gr[%d] at cell %d.\n", bp->nprops.ngr,
                             index);
                    fixed_nuc = 1;
                  }
                }
                if (fixed_nuc == 0) {
                  init_new_grain (bp, bp->nprops.ngr, sbnum, i, j, k, 1);
                }

                bp->gr[bp->nprops.ngr]->TNuc = *c_temp_p;
                bp->gr[bp->nprops.ngr]->TunderNuc = Tunder;
                /* bp->gr[bp->nprops.ngr]->CellConcNuc = phase_diag_on ? *cell_conc : 0; */
              }
            }
          }
        }
        /* end if liquid */
#ifdef OLD_TUNDER
        otup++;
#endif /*OLD_TUNDER */
        c_temp_p++;
        thrp++;
        ogr++;                  /* ptr to outside block cell grain number */
        ngr++;                  /* ptr to inside  block cell grain number */
        ncl++;
      }                         /*end of I loop */
      ogr += 2;

    }                           /* end of J loop */
    ogr += skip;

  }                             /* end of K loop */

  if (sp->ncsolid >= bp->ncsb) {
    fprintf (stderr, "sb_decentred_step.c-cell_nucleation(): SB#%d completely solid.\n", sbnum);
    sp->done = TRUE;
  }
  return (errors);
}

/***************************************************/
/* fraction solid change without diffusion:        */
/* by Wei WANG 16-09-02                            */
/***************************************************/
int fs_change_nodiffuse (BB_struct * bp, int sbnum)
{
  int errors = 0;
  int i, j, k;
  int nx, ny, nz, skip;
  int *ngr, *ogr;
  int cell_temp_on;
  int global_undercooling = 0;
  int **gv_p = bp->gr_array;    /* grain block value array pointer pointer :-/ */
  SB_struct *sp;
  CA_FLOAT *ofs, *c_temp_p;
  CA_FLOAT Tliq, Tunder, Tcell, TliqFixed, min_temp;
  CA_FLOAT del_fs, sum_del_fs, fs_tmp, sum_fs;

  sp = bp->sb[sbnum];
  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  skip = 2 * (nx + 2);

  TliqFixed = bp->mprops.Tliq;
  Tliq = TliqFixed;

  global_undercooling = bp->ctrl->global_undercooling;

  ogr = bp->itmp_one;           /* array of old cell grain, size (nc[i]+2)^3 */
  ofs = bp->ftmp_one;

  ngr = sp->gr;                 /* array of new cell grain, size (nc[i])^3   */
  c_temp_p = sp->c_temp;

  /* copy grain number to temporary buffer */
  icopy_matrix (PAD, ogr, ngr, bp, gv_p, sbnum);

  /*move the ptr to outside corner to inside corner */
  ofs += bp->cubeptr.flist[0][START];

  /* Set up option flags */
  cell_temp_on = bp->ctrl->use_cell_temp;

  sum_del_fs = 0.;
  sum_fs = 0.;
   /***********************************************/
  /* Beginning of main loop(s)                   */
   /***********************************************/
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {

        Tunder = Tliq - *c_temp_p;      /* undercooling = Tliq-Tcell */

        if (*ngr > 0) {         /* the cell is growing or solid */
          if (*ofs <= 1 && Tunder > 0) {
            fs_tmp = *ofs;
            del_fs = growth (bp->mprops.gg_const, bp->mprops.gg_cub, bp->delt, Tunder) / bp->size_c[0];
            *ofs += del_fs;     /* fs in growing cell could exceed 1, and then reduced to 1 in next subroutine */
            /* sum_del_fs += MIN(1.0,*ofs) - fs_tmp; */
          }
          sum_fs += MIN (1.0, *ofs);
        }
        c_temp_p++;
        ngr++;
        ofs++;
      }                         /*end of I loop */
      ofs += 2;
    }                           /* end of J loop */
    ofs += skip;
  }                             /* end of K loop */

  /*  sp->Tvals.del_fs = sum_del_fs / bp->ncsb;
     sp->Tvals.fsavg += sp->Tvals.del_fs; */
  sp->Tvals.fsavg = sum_fs / bp->ncsb;  /* need to modify later */

  /* if (sp->ncsolid >= bp->ncsb) {
     fprintf(stderr, "sb_decentred_step.c-fs_change_nodiffuse(): SB#%d completely solid.\n", sbnum);
     sp->done = TRUE;
     } */
  return (errors);
}

/***************************************************/
/* capture new cell without diffusion of alloy     */
/* by Wei WANG on 16-09-02                         */
/***************************************************/
int capture_octahedron (BB_struct * bp, int sbnum)
{
  SB_struct *sp;
  int i, j, k, n;               /* tmp counter */
  int i_nn, j_nn, k_nn, ijk_nn, capt;
  int nx, ny, nz, skip;
  int *ngr, *ogr;
  CA_FLOAT *opd, *opx, *opy, *opz, *npd, *npx, *npy, *npz;
  CA_FLOAT *ofs;

  CA_FLOAT g00, g10, g20, g01, g11, g21, g02, g12, g22;
  CA_FLOAT rx1, ry1, rz1;
  CA_FLOAT ds, dx, dy, dz, dd, dfs, sum_del_fs;

  CA_FLOAT ncd, ncx, ncy, ncz;

  /* nearest neighbours structure */
  struct nns
  {
    int nx;
    int ny;
    int nz;
  } nn[] = {
  1, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1};
  int NNB = 6;

  sp = bp->sb[sbnum];

  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  skip = 2 * (nx + 2);

  /* Set up local pointers */
  ogr = bp->itmp_one;           /* array of old cell grain, size (nc[i]+2)^3 */
  ngr = sp->gr;                 /* array of new cell grain, size (nc[i])^3   */

  opd = bp->ftmp_dc_d;          /* array of decentred octahedron half diagonal */
  opx = bp->ftmp_dc_x;          /* array of decentred octahedron centre */
  opy = bp->ftmp_dc_y;
  opz = bp->ftmp_dc_z;

  ofs = bp->ftmp_one;

  npd = sp->dc_d;
  npx = sp->dc_x;
  npy = sp->dc_y;
  npz = sp->dc_z;

  /*copy the half diagonal and coodinates of octahedron to temporary buffers */
  fcopy_matrix (PAD, opd, npd, bp, NULL, sbnum);
  fcopy_matrix (PAD, opx, npx, bp, NULL, sbnum);
  fcopy_matrix (PAD, opy, npy, bp, NULL, sbnum);
  fcopy_matrix (PAD, opz, npz, bp, NULL, sbnum);

  /*move the ptr to outside corner to inside corner */
  ogr += bp->cubeptr.flist[0][START];
  opd += bp->cubeptr.flist[0][START];
  opx += bp->cubeptr.flist[0][START];
  opy += bp->cubeptr.flist[0][START];
  opz += bp->cubeptr.flist[0][START];
  ofs += bp->cubeptr.flist[0][START];

  sum_del_fs = 0.;

   /***********************************************/
  /* Beginning of main loop(s)                   */
   /***********************************************/
  /* 2D: decentred square algorithm 3D: decentred octahedron algorithm */
  /* GANDIN97, a cell is captured when its centre touched by a growing square */
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {

        if (*ogr == 0) {        /* cell grain number ==0, the cell is LIQUID */

          capt = FALSE;

          for (n = 0; n < NNB && capt == FALSE; n++) {  /* neighbour loop */
            i_nn = nn[n].nx;
            j_nn = nn[n].ny;
            k_nn = nn[n].nz;
            ijk_nn = i_nn + j_nn * (nx + 2) + k_nn * (nx + 2) * (ny + 2);

            if (*(ogr + ijk_nn) > 0) {  /* the neighbour is solid */

              g00 = bp->gr[*(ogr + ijk_nn)]->g[0][0];
              g10 = bp->gr[*(ogr + ijk_nn)]->g[1][0];
              g20 = bp->gr[*(ogr + ijk_nn)]->g[2][0];

              g01 = bp->gr[*(ogr + ijk_nn)]->g[0][1];
              g11 = bp->gr[*(ogr + ijk_nn)]->g[1][1];
              g21 = bp->gr[*(ogr + ijk_nn)]->g[2][1];

              g02 = bp->gr[*(ogr + ijk_nn)]->g[0][2];
              g12 = bp->gr[*(ogr + ijk_nn)]->g[1][2];
              g22 = bp->gr[*(ogr + ijk_nn)]->g[2][2];

              rx1 = -1. * i_nn - *(opx + ijk_nn);
              ry1 = -1. * j_nn - *(opy + ijk_nn);
              rz1 = -1. * k_nn - *(opz + ijk_nn);

              dx = g00 * rx1 + g10 * ry1 + g20 * rz1;
              dy = g01 * rx1 + g11 * ry1 + g21 * rz1;
              dz = g02 * rx1 + g12 * ry1 + g22 * rz1;

              ds = ABS (dx) + ABS (dy);

              dfs = *(opd + ijk_nn) - ds;

              if (dfs >= 0) {
                capt = TRUE;

                ncd = MIN (*(opd + ijk_nn), 1.);

                dd = *(opd + ijk_nn) - ncd;

                /* determine which corner the new capture cell is close to */
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

                sp->ncsolid++;  /* number of solid cells in sb ++ */
                *ngr = *(ogr + ijk_nn);
                *npd = ncd;
                *npx = *(opx + ijk_nn) + i_nn + g00 * ncx + g01 * ncy + g02 * ncz;
                *npy = *(opy + ijk_nn) + j_nn + g10 * ncx + g11 * ncy + g12 * ncz;
                *npz = *(opz + ijk_nn) + k_nn + g20 * ncx + g21 * ncy + g22 * ncz;
                *ofs = MIN (0.2, dfs);
                /* Robert added -- keep track of grain size and location */
                add_to_grain (bp->gr[*ngr], i, j, k);
                /* sum_del_fs += dfs; */
              }

            }                   /* end of if the neighbour is solid */
          }                     /* end of neighbour loop */
        }
        ogr++;                  /* ptr to outside block cell grain number */
        ngr++;                  /* ptr to inside  block cell grain number */
        npd++;
        npx++;
        npy++;
        npz++;
        opd++;                  /* ptr to decentred octahedron half diagonal */
        opx++;                  /* ptr to decentred octahedron centre */
        opy++;
        opz++;
        ofs++;

      }                         /*end of I loop */
      ogr += 2;
      opd += 2;
      opx += 2;
      opy += 2;
      opz += 2;
      ofs += 2;

    }                           /* end of J loop */
    ogr += skip;
    opd += skip;
    opx += skip;
    opy += skip;
    opz += skip;
    ofs += skip;
  }                             /* end of K loop */

  /* sp->Tvals.del_fs = sum_del_fs / bp->ncsb;
     sp->Tvals.fsavg += sp->Tvals.del_fs; */

  if (sp->ncsolid >= bp->ncsb) {
    fprintf (stderr, "sb_decentred_step.c-capture_octahedron(): SB#%d completely solid.\n", sbnum);
    sp->done = TRUE;
  }

  return (1);
}

/***************************************************/
/* cell index for solidification:                  */
/* by Wei WANG 11-07-02                            */
/***************************************************/
int cell_index (BB_struct * bp, int sbnum)
{
  SB_struct *sp;
  int i, j, k;                  /* tmp counter */
  int nx, ny, nz, skip;
  int *ogr, *ogr_e, *ogr_w, *ogr_n, *ogr_s, *ogr_u, *ogr_d;
  int *nid;
  CA_FLOAT *ofs;  /* I think this is needed for melt back and added by xly */

  sp = bp->sb[sbnum];

  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  skip = 2 * (nx + 2);

  /*set local pointer */
  nid = sp->index;
  ofs = bp->ftmp_one;
  ogr = bp->itmp_one;           /* array of old cell grain, size (nc[i]+2)^3 */
  /*move the ptr to outside corner to inside corner */
  ofs += bp->cubeptr.flist[0][START];
  ogr += bp->cubeptr.flist[0][START];
  ogr_e = ogr + 1;              /* east,west,north,south,up,down: 6 neighbours */
  ogr_w = ogr - 1;
  ogr_n = ogr + nx + 2;
  ogr_s = ogr - nx - 2;
  ogr_u = ogr + (nx + 2) * (ny + 2);
  ogr_d = ogr - (nx + 2) * (ny + 2);

   /***********************************************/
  /* Beginning of main loop(s)                   */
   /***********************************************/
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
        if (*ogr == 0) {        /* this cell is LIQUID */
          *nid = (*ogr_e + *ogr_n + *ogr_w + *ogr_s + *ogr_u + *ogr_d) ? 1 : 0;
          /* 0: all neighbours are LIQUID 1: at least one is SOLID */
        } else if (*ofs == 1) {
          *nid = 3;
        } else {                /* this cell is SOLID or Growing */
          *nid = 2;             /* id=2, the cell is paritally or fully SOLID */
        }
        ofs++;
        nid++;
        ogr++;                  /* ptr to outside block cell grain number */
        ogr_e++;
        ogr_w++;
        ogr_n++;
        ogr_s++;
        ogr_u++;
        ogr_d++;

      }                         /*end of I loop */
      ofs += 2;
      ogr += 2;
      ogr_e += 2;
      ogr_w += 2;
      ogr_n += 2;
      ogr_s += 2;
      ogr_u += 2;
      ogr_d += 2;

    }                           /* end of J loop */
    ofs += skip;
    ogr += skip;
    ogr_e += skip;
    ogr_w += skip;
    ogr_n += skip;
    ogr_s += skip;
    ogr_u += skip;
    ogr_d += skip;
  }                             /* end of K loop */

  return (1);
}

/***************************************************/
/* Cell temperaure               */
/* by Wei Wang on 06-11-01                         */
/***************************************************/
int cell_temp (BB_struct * bp, int sbnum)
{
#ifdef JUNK
  /* This is not necessary - r.atwood 27 feb 2003 */
  /**  \todo  change the calling routines so that cell_temp is not used -- general - easy */
  /* and remove the funciton */
  SB_struct *sp;
  int cell_temp_on;
  int i, j, k;                  /* tmp counter */
  int nx, ny, nz;
  CA_FLOAT Tcell;
  CA_FLOAT *nct;

  sp = bp->sb[sbnum];

  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];

  nct = sp->c_temp;

  if (bp->ctrl->use_cell_temp == TRUE) {
    cell_temp_on = TRUE;
  } else {
    cell_temp_on = FALSE;
  }

  if (cell_temp_on == TRUE) {   /* non-isothermal condition */
    if (bp->ctrl->fgrid_input) {
      for (k = 0; k < nz; k++) {
        for (j = 0; j < ny; j++) {
          for (i = 0; i < nx; i++) {
            *nct = fg_temp_calc (bp, sbnum, i, j);
            nct++;
          }                     /* I */
        }                       /* J */
      }                         /* K */
    } else {
      for (k = 0; k < nz; k++) {
        for (j = 0; j < ny; j++) {
          for (i = 0; i < nx; i++) {
            *nct = cell_temp_calc_cc (bp, sbnum, i, j);
            nct++;
          }                     /* I */
        }                       /* J */
      }                         /* K */
    }

  } else {                      /* isothermal conditon  */
    Tcell = bp->sb[sbnum]->Tvals.Tavg;
    for (k = 0; k < nz; k++) {
      for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
          *nct = Tcell;
          nct++;
        }                       /* I */
      }                         /* J */
    }                           /* K */
  }

#endif /*JUNK*/
    return (1);
}


/***************************************************/
/* fraction solid change with diffusion of alloy   */
/* by Wei Wang on 07-11-01                         */
/***************************************************/

int fs_change_diffuse (BB_struct * bp, int sbnum)
{
  int i, j, k, nx, ny, nz, skip;
  int *ngr, *ogr;
  int use_curv;
  SB_struct *sp;
  CA_FLOAT *ofs, *nfs, *ncv = NULL;
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

  /**todo: Get rid of use_curv by simplifying the curvature options */
  use_curv=((bp->ctrl->curvature_3D != 0) || (bp->ctrl->curvature_2D != 0));

  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  skip = 2 * (nx + 2);

  part_coef = bp->mprops.alloyprops[0].part_coef[0];
  liq_slope = bp->mprops.alloyprops[0].m_solute[0];

  k_inv = 1. / part_coef;
  km = (1 - part_coef);
  km_inv = 1. / km;
  m_inv = 1 / liq_slope;        /* Inverse of liquidus slope M */

  temp_pure = bp->mprops.tp;


  D_x = bp->mprops.alloyprops[0].Dliq / bp->size_c[0] / bp->size_c[0];     /* D_L / (dx)^2  */
  D_xt = D_x * bp->delt;        /* D_L * dt / (dx)^2 */

  sp = bp->sb[sbnum];
  ngr = sp->gr;
  nfs = sp->c_fs;
  ncl = sp->c_sol_alloy;
  nce = sp->c_eqv_alloy;

  if(use_curv){
    ncv = sp->curv;
    coef_curv = bp->mprops.gibbs_thomson / bp->size_c[0];
  }

  nct = sp->c_temp;
  ofs = bp->ftmp_one;
  oce = bp->ftmp_five;
  ogr = bp->itmp_one;

  ofs += bp->cubeptr.flist[0][START];
  oce += bp->cubeptr.flist[0][START];
  ogr += bp->cubeptr.flist[0][START];

  sum_del_fs = 0.0;

  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
        Tcell = *nct;
        if (*ngr <= 0) {        /* the cell is LIQUID or NOT-CASTING */
          *ncl = *nce;          /* C_E = C_S*fs + C_L*(1-fs), fs=0 */
        } else {                /* the cell is paritially or fully SOLID */
          fs = *ofs;
          if(use_curv){
             c_int = m_inv * (Tcell - temp_pure + coef_curv * (*ncv));
          }else{
             c_int = m_inv * (Tcell - temp_pure);
          }

          if (fs < 1) {         /* 0 <= fs < 1 */
            dce = *nce - *oce;
            /* dcl = (c_int - *ncl) * D_xt / MAX(0.01, 1-fs) *0.1; */
            dcl = (c_int - *ncl) * D_xt * 2 + dce / (1 - km * fs);
            /*  dfs = -1.*dce * km_inv / *ncl; *//*fist guess */
            /* dcl = (c_int - *ncl) * MAX(0.005,dfs) / MAX(0.005, 1-fs); *//*derived from 1-D steady growth */
            dfs = (dcl * (1 - km * fs) - dce) * km_inv / *ncl;

            /*  dcl = (c_int + dce * x_Dt*0.125*(1-fs) * (ca0+ABS(sa0))- *ncl) * D_xt / MAX(0.001, 1-fs) / (ca0+ABS(sa0)); */
            *ncl += dcl;
            *ofs += dfs;
            *ofs = MIN (1.0, *ofs);
            if (*ofs < 0.) {
              *ofs = 0.;
              *ngr = 0;
              *ogr = 0;
              *ncl = *nce;
              sp->ncsolid--;    /* number of solid cells decreases by one */
            }
            sum_del_fs += *ofs - fs;    /* increase of the fraction solid */
          } else {              /* fs == 1 */
            *ncl = *nce * k_inv;        /* C_E = C_S, when fs = 1 */
          }

        }
        nct++;
        ngr++;
        nfs++;
        ncl++;
        nce++;
        if(use_curv) ncv++;
        ofs++;
        oce++;
        ogr++;
      }
      ofs += 2;
      oce += 2;
      ogr += 2;
    }
    ofs += skip;
    oce += skip;
    ogr += skip;
  }

  /* added in by Wei WANG on 06-11-03 */
  sp->Tvals.del_fs = sum_del_fs / bp->ncsb;     /* increase of average fs in SubBlock */
  sp->Tvals.fsavg += sp->Tvals.del_fs;  /* average fs in SubBlock */

  return (1);
}

/***************************************************/
/* Growth of square/octahedron                     */
/*   with/without diffuse of alloy                 */
/* by Wei WANG on 07-10-02                         */
/***************************************************/
int grow_octahedron (BB_struct * bp, int sbnum)
{
  SB_struct *sp;
  int i, j, k;                  /* tmp counter */

  int nx, ny, nz, skip;
  int *ngr;

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

   /***********************************************/
  /* Beginning of main loop(s)                   */
   /***********************************************/

  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {

        if (*ofs <= 0) {

          /* cell grain number ==0, the cell is LIQUID */

        } else if (*ngr > 0) {  /* the cell is captured by a grain */
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
      }                         /*end of I loop */
      ofs += 2;
    }                           /* end of J loop */
    ofs += skip;
  }                             /* end of K loop */

  return (1);
}

/***************************************************/
/* capture new cell with diffusion of alloy        */
/* by Wei WANG on 15-07-02                         */
/***************************************************/
int capture_octahedron_diffuse (BB_struct * bp, int sbnum)
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

  CA_FLOAT *nct;
  CA_FLOAT g00, g10, g20, g01, g11, g21, g02, g12, g22;
  CA_FLOAT rx, ry, rz, dx, dy, dz, ds1;

  CA_FLOAT ncd, dd, ncx, ncy, ncz;
  CA_FLOAT part_coef, kInv, km;
  CA_FLOAT fs;

/* stuff for meltback ? */
  int melt;
  CA_FLOAT melt_dt;
  CA_FLOAT liq_slope, liq_slope_1;
  CA_FLOAT Tliq, Tsol, Tpure, Tsol_ce;
  CA_FLOAT nco;

  /*CA_FLOAT t/*, t1, t2, v, v_coef1 */ ;

  /*t = bp->sim_time; */
  /*t1 = bp->time_velo1;
     v = bp->velocity;
     v_coef1 = bp->velo_coef1;
     t2 = t1 + v/v_coef1; */

  /* nearest neighbours structure */
  struct nns
  {
    int nx;
    int ny;
    int nz;
  } nn[] = {
  1, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 1, 1, 0, -1, 1, 0, -1, -1, 0, 1, -1, 0};
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

  part_coef = bp->mprops.alloyprops[0].part_coef[0];
  km = 1 - part_coef;
  melt_dt = bp->mprops.dt_melt;

  liq_slope = bp->mprops.alloyprops[0].m_solute[0];
  liq_slope_1 = 1. / liq_slope;
  Tliq = bp->mprops.Tliq;
  Tsol = bp->mprops.Tsol;
  nco = bp->mprops.alloyprops[0].Cinit;

  Tpure = liq_slope * nco * (-1) + Tliq;

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
  ncl = sp->c_sol_alloy;
  nce = sp->c_eqv_alloy;

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

            if (gr_nn > 0) {    /* the neighbour is solid */
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
              ds1 = ABS (dx) + ABS (dy) + ABS (dz);     /*3D */
              /*ds2 =  2.0 * MIN(ABS(dx), ABS(dy)); */
              /*ds2 = SQRT(dx*dx + dy*dy) + 0. * (ca0 + ABS(sa0) - 1.0);  *//* round the tip */

              /* fs = *(opd+ijk_nn) - MAX(ds1, ds2); */
              fs = *(opd + ijk_nn) - ds1;

	      /**todo: Convert to reproducible random sequence */
	      /* to enable restart with reproducible results */
	      /* it seems also that this always adds fs, it's biased */
	      /* doesn't this increase the growth rate compared to what's */
	      /* calculated ? */
              fs += 0.05 * drand48 ();  /* add some numerical noise */

              if (fs > 0.) {    /* probalistic growth */
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

                sp->ncsolid++;  /* number of solid cells in sb ++ */
                *ngr = *(ogr + ijk_nn);
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
                *ncl = *nce / (1 - km * fs);

              }
            }                   /* end of if the neighbour is solid */
          }                     /* end of neighbour loop */

        }

        else if (*nid == 3) {
          melt = FALSE;
          Tsol_ce = Tpure - *nce / part_coef * liq_slope * (-1);
          if ((*nct - Tsol_ce) >= melt_dt) {
            melt = TRUE;
            *ncl = (*nct - Tpure - melt_dt) * liq_slope_1;
            *ofs = (*ncl - *nce) / (*ncl * km);
            sp->ncsolid--;      /*number of solid cells decrease */
          }
        }
        nid++;
        nct++;
        ogr++;                  /* ptr to outside block cell grain number */
        ngr++;                  /* ptr to inside  block cell grain number */
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
        ncl++;
        nce++;

      }                         /*end of I loop */
      ogr += 2;
      opd += 2;
      opx += 2;
      opy += 2;
      opz += 2;
      ofs += 2;
    }                           /* end of J loop */
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
char const *sb_decentred_step_c ()
{
  static char const rcsid[] = "$Id: sb_decentred_step.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}
