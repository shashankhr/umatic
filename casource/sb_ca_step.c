
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
/*      sb_umat_step.c:                                           */
/*                                                              */
/* The main subroutine to calculate growth and extend into new  */
/* cells.:                                                      */
/*                                                              */
/*                                                              */
/*                                                              */
/****************************************************************/
/*RCS ID: $Id: sb_ca_step.c 1373 2008-08-27 20:51:52Z  $*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "machine.h"
#include "blocks.h"
#include "umat_matrix.h"
#include "props.h"
/* prototypes for the nucleation function choices */
#include "nucfuncs.h"

/* Prototype defs */
extern void add_to_grain (Ind_grain * gp, int i, int j, int k);

/* from p_growth.c */
extern CA_FLOAT growth_primary (CA_FLOAT gg_const, CA_FLOAT gg_cub, CA_FLOAT dt, CA_FLOAT Tunder, int eut_flag, CA_FLOAT cell_fs_n_eut,
                                CA_FLOAT cell_fs_b_eut);
extern CA_FLOAT growth_eutectic (CA_FLOAT gg_const, CA_FLOAT gg_cub, CA_FLOAT dt, CA_FLOAT Tunder, int eut_flag,
                                 CA_FLOAT cell_fs_n_eut, CA_FLOAT cell_fs_b_eut);
CA_FLOAT particle_growth (CA_FLOAT pg_const_a, CA_FLOAT pg_const_b, CA_FLOAT dt, CA_FLOAT Tunner, CA_FLOAT p_conc);

/*from trans_interp_calc.c*/
extern CA_FLOAT trans_interp_calc (FGrid_str * fg, BB_struct * bp, int sbnum, int x, int y);

/* from sb_temp_calc.c */


extern int init_new_grain (BB_struct * bp, int igr, int sbnum, int xcell, int ycell, int zcell, int nc);

/*beginning of sb_umat_step*/

/** 
 * Do a normal ca step
 *
 * @callgraph
 * @callergraph
*/
int sb_umat_step (BB_struct * bp, int sbnum)
{
  int errors = 0, ngrowing = 0;

#ifdef OLD_TUNDER
  CA_FLOAT *otu, *otup;
#endif /*OLD_TUNDER */
  int i, j, k, n;
  int *oni, *nni;
  int *tmp_oni, *tmp_nni, tmp_n_neigh;
  int nx, ny, nz, skip;
  int inter_count;
  int inter_flag;
  SB_struct *sp = bp->sb[sbnum];
  Ind_grain *this_grain;

/******added for the multi component****/
  CA_FLOAT **conc_multi;
  MultiS_struct *ms;

/****************************************/
  CA_FLOAT Tliq, Tunder, prob, TliqFixed;
  CA_FLOAT *cell_conc;
  CA_FLOAT *c_temp_p, *ofs, *nfs, *op, *np, *sol, *solp, *sfs, *sfsp, *thr, *thrp;
  CA_FLOAT *local_dfs_primary, *local_dfs_primary_start;
  CA_FLOAT *local_dfs_eutectic, *local_dfs_eutectic_start;
  int *bin_flag, *bin_flag_start;
  CA_FLOAT sum_del_fs, dr, dfs, dfs_tmp, extra_fs;
  CA_FLOAT dr_primary, dr_eutectic, dfs_primary, dfs_eutectic;
  CA_FLOAT oldfs = 0;
  int *ngr, *ogr;
  int oriented_on, cell_nuc_on, cell_temp_on, n_neigh;

  int particle_on = bp->ctrl->particle;
  int phase_diag_on = bp->ctrl->phase_diag_on;
  int global_undercooling = bp->ctrl->use_global_undercooling;

  int eut_flag = 0;
  static int nucmsg = 0;
  CA_FLOAT sch_sum;
  int (*cell_nuc_func) (BB_struct * bp, CA_FLOAT Tunder, CA_FLOAT argthree);
  Ctrl_str *cp = bp->ctrl;

  int index_ca;

  dfs_tmp = 0;
  ms = &(bp->MultiSvals);
  sfs = sfsp = bp->ftmp_three;

  /* set up phase diagram mode */
  phase_diag_on = bp->ctrl->phase_diag_on;
  global_undercooling = cp->use_global_undercooling;
  particle_on = bp->ctrl->particle;

  if (cp->use_csol_alloy) {
    cell_conc = sp->c_sol_alloy;
  }

  /* Set up local pointers */
#ifdef OLD_TUNDER
  otup = otu = bp->old_Tunder;
#endif /*OLD_TUNDER */
  local_dfs_primary = local_dfs_primary_start = sp->cell_dfs_primary;
  local_dfs_eutectic = local_dfs_eutectic_start = sp->cell_dfs_eutectic;
  op = ofs = bp->ftmp_one;
  bin_flag = bin_flag_start = bp->bin_flag;
  /*pop = pofs = sp->fs_n_eut; */
  /* eop = eofs = sp->fs_b_eut; */
  np = nfs = sp->c_fs;
  solp = sol = sp->c_sol;
  ogr = bp->itmp_one;
  ngr = sp->gr;
  c_temp_p = sp->c_temp;        /* cell temperature array pointer */

  /* dummy ref not used for cell_nuc_func,            */
  /* should be old Tunder when needed or threshold,   */
  /* but has to be something if neither mode is on!   */
  thr = thrp = sp->c_sol;

  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];

  /*set up local values */
  skip = 2 * (nx + 2);
  oriented_on = bp->nprops.oriented;
  TliqFixed = bp->mprops.Tliq;
  Tliq = TliqFixed;
  sum_del_fs = 0.0;

  /* Set up option flags */

  cell_nuc_on = cp->use_cell_nuc;

  if (cell_nuc_on) {            /*Select the cell nuc function -- should move to higher level! */
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
      thr = thrp = sp->c_nuc_thresh;
      break;
    default:
      fprintf (stderr, "ERROR: sb_umat_step: cannot find cell_nuc function; Nucleation model undefined.(%i)\n", (bp->nprops.nmodel));
      exit (0);
      break;
    }
  }

  /*end of set up cell nuc function */
  /*   if (bp->nprops.nmodel == N_DIST) cell_nuc_func = dist_cell_nuc; */
  /*   else cell_nuc_func = cell_nuc; */
  /* set up temporary neighbourhood */
  n_neigh = tmp_n_neigh = bp->ctrl->n_neigh;
  switch (n_neigh) {
  case NEIGH_6:
    oni = bp->nbhd.onq;
    nni = bp->nbhd.nnq;
    break;
  case NEIGH_8:
    oni = bp->nbhd.onhl;
    nni = bp->nbhd.nnhl;
    break;
  case NEIGH_10:
    oni = bp->nbhd.ono;
    nni = bp->nbhd.nno;
    break;
  case NEIGH_26:
    oni = bp->nbhd.onc;
    nni = bp->nbhd.nnc;
    break;
  }
  tmp_oni = oni;
  tmp_nni = nni;
  /* end set up temporary neighbourhood */

  /*sciel fraction solid */
  sfsp = sfs + bp->cubeptr.flist[0][START];

  /*outside frac solid pointer */
  op = ofs + bp->cubeptr.flist[0][START];
  bin_flag = bin_flag_start + bp->cubeptr.flist[0][START];
  local_dfs_primary = local_dfs_primary_start;
  local_dfs_eutectic = local_dfs_eutectic_start;
  /*pop = pofs; */
  /* eop = eofs; */

  /*outside grain */
  ogr = ogr + bp->cubeptr.flist[0][START];

  /*nside pointer */
  np = nfs;

  /*keep track of scheil frac solid seperately */
  if (bp->ctrl->scheil == TRUE) {
    sch_sum = 0;
  }

   /***********************************************/
  /* Beginning of main loop(s)                   */
   /***********************************************/

  sp->ncsolid = 0;

  for (k = 0, index_ca = 0; (k < nz && index_ca < bp->tnc[0] * bp->tnc[1] * bp->tnc[2]); k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {

        /*possibly unroll this section? */
        if (phase_diag_on) {
          eut_flag = FALSE;
          if ((cp->diffuse_alloy_multi == FALSE) && (cp->diffuse_alloy == TRUE)) {
            Tliq = cell_liq_calc (*cell_conc, &(bp->mprops));
          }
          /* Multi component */
          if (cp->diffuse_alloy_multi == TRUE) {

            /** \todo integrate Ludovic Thuinets' poly component into non-decentered version */
          }
          /* end Multi component */
          if (bp->ctrl->diffuse_alloy_multi == FALSE && cp->diffuse_alloy_poly == FALSE) {
            /** \todo  if the gas has an effect on the solution then this needs rethinking -- multicomponent */

            /* test if cell should solidify as eutectic */ 
            if (*c_temp_p <= bp->mprops.alloyprops[0].T_eut && *cell_conc >= bp->mprops.alloyprops[0].c_eut) {
              eut_flag = TRUE;
            }
          }
        }

        /*undercooling is always Tliq-Tcell */

        Tunder = Tliq - *c_temp_p;
#         ifdef OLD_TUNDER
        if (bp->step <= 1)
          *otup = Tunder;
#         endif /*OLD_TUNDER */

          /************************************/
        /*If the cell is GREATER THAN SOLID */
        /* this shouldnt be able to happen. */
          /************************************/
        if (*np > SOLID) {
          fprintf (stderr, "ERROR:sb_umat_step: cell > solid! sb %i, i %i, j %i, k %i, step %i\n", sbnum, i, j, k, bp->step);
#ifdef ERROR_EXIT
          fprintf (stderr, "EXITING!\n");
          exit (3);
#endif
        }
          /*********************************************/
        /* check the state of the cell and apply the */
        /* appropriate algorithm                     */
          /*********************************************/
        if (*np == NOT_CASTING) {
#ifdef EXTRA_TRAPS
          int dumb;

          /*******************************************/
          /*                                         */
          /*     if the cell is NOT IN THE CASTING   */
          /*                                         */
          /*******************************************/

                    /********************/
          /* do nothing       */
                    /********************/
          dumb = 0;             /* need to be able to set a trap here */
#endif /* EXTRA_TRAPS */
        } else if (*np == SOLID) {

          /*******************************************/
          /*                                         */
          /*     if the cell is SOLID                */
          /*                                         */
          /*******************************************/

          /*calculate schiel frac solid */
          if (bp->ctrl->scheil == TRUE) {
            *sfsp = schiel (*c_temp_p);
            sch_sum += *sfsp;
            if (*sfsp >= 1.0) {
              sp->ncsolid++;
            }                   /* end if schiel fully solid */
          } /*end if schiel */
          else {
            sp->ncsolid++;
          }

          /*******************************************/
          /*                                         */
          /*     if the cell is LIQUID               */
          /*                                         */
          /*******************************************/
        } else if (*np == LIQUID) {
          if (cell_nuc_on) {
#ifdef GLOBAL_UND_NUC
            /*for global undercooling nucleation model */
            if ((global_undercooling))
              Tunder = TliqFixed - *c_temp_p;
#endif

#              ifdef OLD_TUNDER /* only used when trying to calculate threshold by integrating */
            /* the undercooling from the previous step ..  */
            if (bp->ctrl->block_nuc) {
              fprintf (stderr, "You can't do that! You have to debug the code.\n");
              fprintf (stderr, "You can't use OLD_TUNDER and block_nuc at the same time\n.");
              exit (1);
            }
            /* then cell_nuc_func needs to have the old Tunder pointer otup defined */
            if (Tunder > 0 && (*cell_nuc_func) (bp, Tunder, *otup) == 1)
#              else /*OLD_TUNDER */
            {                   /*DEBUGGING section */
              static int nnucmsg = 0;

              if (nnucmsg < MAX_NUC_MSG && Tunder > 0) {
                if (bp->ctrl->diffuse_alloy_multi == FALSE) {
                  fprintf (stderr, "cell temp[%g] cell conc[%g] cell Tliq[%g]\n", *c_temp_p, phase_diag_on ? *cell_conc : 0, Tliq);
                }
                nnucmsg++;
              }
            }                   /*end DEBUGGING section */
            if (Tunder > 0 && (*cell_nuc_func) (bp, Tunder, *thrp) == 1)
#              endif /*OLD_TUNDER */
            {                   /* trheshold is exceeded so ... */

               /*********************/
              /* Nucleate a Grain  */
               /*********************/

              /* error checks */
              if (nucmsg < MAX_NUC_MSG) {
                if (bp->ctrl->diffuse_alloy_multi == FALSE) {
                  fprintf (stderr, "cell temp[%g] cell conc[%g] cell Tliq[%g]\n", *c_temp_p, phase_diag_on ? *cell_conc : 0, Tliq);
                }
                nucmsg++;
              }
              if (bp->nprops.ngr == bp->nprops.gd_max_total) {
                fprintf (stderr, "ERROR: sb_umat_step - Max # grains set by user [%d] exceeded. Increase option MaxTotGrains.\n",
                         bp->nprops.gd_max_total);
              } else {
                /* no errors */
                /* create new grain and set up */
                /*record nucleation values in structure */
                *op = EMBRYO;
                (bp->nprops.ngr)++;
                (bp->sb[sbnum]->ngr)++;
                *ogr = bp->nprops.ngr;
                /*
                   init_new_grain(bp, *ogr, sbnum, i + j * nx + k * nx * ny, 1);
                 */
                init_new_grain (bp, *ogr, sbnum, i, j, k, 1);
                this_grain = bp->gr[*ogr];
#ifdef DBM
                /* diagnostic output from dbmalloc (debug malloc) library */
                {
                  static int dbm_trapgrain = 0;

                  fprintf (stderr, "grain created,trap,num,cell,mem:,%i,%i,%i,%x\n", dbm_trapgrain, this_grain->num,
                           i + j * nx + k * nx * ny, this_grain);
                  dbm_trapgrain++;
                }
#endif           /*DBM*/
                  this_grain->TNuc = *c_temp_p;
                this_grain->TunderNuc = Tunder;
                this_grain->NucTh = *thrp;
                if (bp->ctrl->diffuse_alloy_multi == FALSE) {
                  this_grain->CellConcNuc = phase_diag_on ? *cell_conc : 0;
                }

              }                 /* end create new grain and set up */
            }                   /* end nucleation threshold test */
          }                     /* end if cell nuc on */
          if (bp->ctrl->scheil == TRUE) {
            *sfsp = 0;
          }
#ifdef OLD_TUNDER
          *otup = Tunder;       /*save the old Tunder */
#endif /*OLD_TUNDER */

        } else {                /* end of LIQUID section */
          /*******************************************/
          /*                                         */
          /*     if the cell is GROWING              */
          /*                                         */
          /*******************************************/
          /*The cell must be growing */
          ngrowing++;
          if (phase_diag_on) {
             /*******************************************/
            /* Phase Diagram -- Must consider growth   */
            /* fixing later by diffuse subroutine.     */
            /* Also use growth pct.-- grow when not    */
            /* fully solid.                            */
             /*******************************************/
#ifdef GLOBAL_UND_GRO
            /*for global undercooling growth model */
            if ((global_undercooling))
              Tunder = TliqFixed - *c_temp_p;
#endif

            /*growth will be fixed in diffuse subroutine */
            if (Tunder > 0) {   /*calculate growth rate per cell - repeated below */
              if (bp->ctrl->particle == TRUE) { /*replace with function pointer later */
                dr = particle_growth (bp->mprops.gg_const, bp->mprops.gg_cub, bp->delt, Tunder, *cell_conc);
                /* use fixed initial concentration of particles in this case */
              } else {
                if (bp->ctrl->diffuse_alloy_multi == 0) {
                  dr = growth (bp->mprops.gg_const, bp->mprops.gg_cub, bp->delt, Tunder);
                  dfs = dr / bp->size_c[0];
                  dfs_tmp = oriented_on ? (dfs * (bp->gr[*ngr]->gro_fact)) : dfs;
                  sum_del_fs += dfs_tmp;
                  *op += dfs_tmp;
                } else {
                  if (*bin_flag == 0) {
                    dr_primary =
                      growth_primary (bp->mprops.gg_const, bp->mprops.gg_cub, bp->delt, Tunder, *bin_flag, sp->fs_n_eut[index_ca],
                                      sp->fs_b_eut[index_ca]);
                    dfs_tmp = dfs_primary = dr_primary / bp->size_c[0];
                    *local_dfs_primary = dfs_tmp;
                    sum_del_fs += dfs_tmp;
                    /* *pop += dfs_primary; */
                    *op += dfs_tmp;
                    if (*op > SOLID) {
                      extra_fs = *op - SOLID;
                      sum_del_fs -= extra_fs;
                      /* *pop = SOLID; */
                      *op = SOLID;
                      sp->ncsolid++;
                      bp->gr[*ogr]->ngrow--;
                    }
                  } else {
                    dr_primary =
                      growth_primary (bp->mprops.gg_const, bp->mprops.gg_cub, bp->delt, Tunder, *bin_flag, sp->fs_n_eut[index_ca],
                                      sp->fs_b_eut[index_ca]);
                    dr_eutectic =
                      growth_eutectic (bp->mprops.gg_const, bp->mprops.gg_cub, bp->delt, Tunder, *bin_flag, sp->fs_n_eut[index_ca],
                                       sp->fs_b_eut[index_ca]);
                    dfs_primary = (dr_primary / (2.0 * bp->size_c[0]));
                    *local_dfs_primary = dfs_primary;
                    dfs_eutectic = (dr_eutectic / (2.0 * bp->size_c[0]));
                    *local_dfs_eutectic = dfs_eutectic;
                    dfs_tmp = dfs = dfs_primary + dfs_eutectic;
                    sum_del_fs += dfs_tmp;
                    /*      *pop += dfs_primary; */
                    /*      *eop += dfs_eutectic; */
                    *op += dfs_tmp;
                    if (*op > SOLID) {
                      extra_fs = *op - SOLID;
                      sum_del_fs -= extra_fs;
                      /*       *pop -=extra_fs/2.0; */
                      /*       *eop -=extra_fs; */
                      *op = SOLID;
                      sp->ncsolid++;
                      bp->gr[*ogr]->ngrow--;
                    }
                  }
                }
              }

              if (dfs > DFS_WARNING) {
                bp->dfs_warn++;
              }
              if (dfs > DFS_CAP) {
                dfs = DFS_CAP;
                bp->dfs_cap++;
              }
            } else {
              dr = 0.0;
              dfs = 0.0;
            }                   /*end caluclate growth rate per cell */

            if (*np >= bp->fsgrow) {
#include "oriented.inc"
              /*use EMBRYO if there is no extra fs */
              extra_fs = EMBRYO;
              /*special case if solid exceeded in this step */
              if (bp->ctrl->diffuse_alloy_multi == 0) {
                if (*op > SOLID) {
                  extra_fs = *op - SOLID;
                  sum_del_fs -= extra_fs;

                  if (extra_fs > CAP_GROWTH) {
                    extra_fs = CAP_GROWTH;
                    bp->extrafs_cap++;
                  }

                  *op = SOLID;
                  sp->ncsolid++;
                  bp->gr[*ogr]->ngrow--;
                }               /*end solid exceeded */
              }

              for (n = 0; n < tmp_n_neigh; n++) {       /* neighbour loop */
                if (*ogr == 0) {
                  continue;
                } else {
                  if (*(ogr + tmp_oni[n]) == LIQUID) {  /*activate neighbour */
                    *(op + tmp_oni[n]) = extra_fs;
                    *(ogr + tmp_oni[n]) = *ogr;
                    /* add the cell to the grain structure */
                    /* keeping track of its size and location */
                    add_to_grain (bp->gr[*ogr], i, j, k);
                    sum_del_fs += extra_fs;
                  }             /* end activate neigbour */
                }
              }                 /*end check neighbour loop */
            } else {            /* end if frac >= fsgrow */
              if (*op >= SOLID) {       /*annoying exception */
                /*solid has increased beyond 1.0 from below growth_pct */
                extra_fs = *op - 0.9999;        /*allow diffuse to correct */
                sum_del_fs -= extra_fs;
                *op = 0.9999;
                bp->gr[*ogr]->ngrow++;
              }                 /* end annoying exception */
            }                   /* end if frac ! >= fsgrow */
          } else {              /*phase diag not on */

             /*******************************************/
            /* Non Phase Diagram                       */
            /* Does not use growth pct. just now...    */
            /*                                         */
             /*******************************************/

            /*growth is finalised, eutectic doesn't exist, etc */
            if (cell_temp_on || bp->ctrl->external == TRUE) {
              if (Tunder > 0) { /*calculate growth rate per cell - repeated above */
                if (bp->ctrl->particle == TRUE) {       /*replace with function pointer later */
                  dr = particle_growth (bp->mprops.gg_const, bp->mprops.gg_cub, bp->delt, Tunder, *cell_conc);
                  /* use fixed initial concentration of particles in this case */
                } else {
                  dr = growth (bp->mprops.gg_const, bp->mprops.gg_cub, bp->delt, Tunder);
                }
                dfs = dr / bp->size_c[0];
                if (dfs > DFS_WARNING) {
                  bp->dfs_warn++;
                }
                if (dfs > DFS_CAP) {
                  dfs = DFS_CAP;
                  bp->dfs_cap++;
                }
              } else {
                dr = 0.0;
                dfs = 0.0;
              }                 /*end caluclate growth rate per cell */
            }                   /*end if cell temp on */
            dfs_tmp = oriented_on ? (dfs * (bp->gr[*ngr]->gro_fact)) : dfs;
            *op += dfs_tmp;
            sum_del_fs += dfs_tmp;
            if (*op >= SOLID) {
#include "oriented.inc"
              extra_fs = *op - 1.0;
              sum_del_fs -= extra_fs;
              *op = 1.0;
              bp->gr[*ogr]->ngrow--;
              if (bp->ctrl->scheil == FALSE) {
                sp->ncsolid++;
              }
              for (n = 0; n < tmp_n_neigh; n++) {       /* neighbour loop */
                if (*ogr == 0) {
                  continue;
                } else {
                  if (*(ogr + tmp_oni[n]) == LIQUID) {  /*activate neighbour */
                    *(op + tmp_oni[n]) = extra_fs;
                    *(ogr + tmp_oni[n]) = *ogr;
                    bp->gr[*ogr]->ncells++;
                    bp->gr[*ogr]->ngrow++;
                    sum_del_fs += extra_fs;
                    /*set min and max cells for grain */
                    bp->gr[*ogr]->max[0] = MAX (bp->gr[*ogr]->max[0], i);
                    bp->gr[*ogr]->max[1] = MAX (bp->gr[*ogr]->max[1], j);
                    bp->gr[*ogr]->max[2] = MAX (bp->gr[*ogr]->max[2], k);
                    bp->gr[*ogr]->min[0] = MIN (bp->gr[*ogr]->min[0], i);
                    bp->gr[*ogr]->min[1] = MIN (bp->gr[*ogr]->min[1], j);
                    bp->gr[*ogr]->min[2] = MIN (bp->gr[*ogr]->min[2], k);
                  }             /* end activate neigbour */
                }
              }                 /*end check neighbour loop */
            }
            /* end if frac >= SOLID */
            if (bp->ctrl->scheil == TRUE) {
              *sfsp = schiel (*c_temp_p);
              sch_sum += *sfsp;
              if (*sfsp >= 1.0) {
                sp->ncsolid++;
              }
            }                   /*end schiel fs calculation */
          }                     /*end phase diag not on */
        }                       /* end of GROWING section */

        if (*op > 1.0) {
          fprintf (stderr,
                   "WARNING:SB_CA_STEP: Something Weird is Happening!\n *op > 1.0\nsb %i, i %i, j %i, k %i, step %i, *op %.5g\n",
                   sbnum, i, j, k, bp->step, *op);
        }

        /* end of somethign wierd message */
        /* increment the pointers depending on   */
        /* whether they are padded or non-padded */
        /* arrays                                */
#ifdef OLD_TUNDER
        otup++;
#endif /*OLD_TUNDER */
        thrp++;                 /*nuc threshold */
        sfsp++;                 /*shciel fraction solid */
        np++;
        op++;
        bin_flag++;
        local_dfs_primary++;
        local_dfs_eutectic++;
        /* pop++; */
        /* eop++; */
        ngr++;
        index_ca++;
        solp++;                 /* gas concentration */
        ogr++;
        c_temp_p++;
        if (phase_diag_on || particle_on)
          cell_conc++;          /* alloy conc */
      }                         /*end of i loop -- x direction */
      sfsp += 2;
      op += 2;
      ogr += 2;
      bin_flag += 2;
    }                           /* end of j loop - y direction */
    sfsp += skip;
    op += skip;
    ogr += skip;
    bin_flag += skip;
  }                             /* end of k loop - z direction */

   /************************************************/
  /* end of horrible nested loops                 */
   /************************************************/

#ifdef PROBE_CELL
  {
    FILE *cell_file;
    char fname[64] = "";

    sprintf (fname, "cell%i.txt", PROBE_CELL);
    cell_file = fopen (fname, "a+");
    fprintf (cell_file, "%.5g,%.10g\n", bp->sim_time, sp->c_temp[PROBE_CELL]);
    fclose (cell_file);
  }
#endif /* PROBE_CELL */

  oldfs = sp->Tvals.fsavg;
  if (bp->ctrl->scheil == TRUE) {
    sp->Tvals.fsavg = sch_sum / bp->ncsb;
  } else {
    sp->Tvals.del_fs = sum_del_fs / (bp->ncsb - sp->nmould);
    sp->Tvals.fsavg += sp->Tvals.del_fs;
  }
  if (oldfs > sp->Tvals.fsavg) {
    static int fsmsg = 0;

    if (fsmsg < 5) {
      fprintf (stderr, "ERROR:sb_umat_step: fs decreasing!\n");
      fsmsg++;
    }
  }

  if (sp->ncsolid >= bp->ncsb) {
    fprintf (stderr, "SB_umat_step(): SB#%d completely solid.\n", sbnum);
    sp->done = TRUE;
  }
  return (errors);
}                               /* end of sb_umat_step */

/* subroutine to return rcs-id string */
char const *rcs_id_sb_umat_step_c ()
{
  static char const rcsid[] = "$Id: sb_ca_step.c 1373 2008-08-27 20:51:52Z  $";

  return (rcsid);
}

/*
*/
