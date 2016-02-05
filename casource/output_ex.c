
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
/*  output_excel.c:						*/
/*  All subroutines related to output in EXCEL format.          */
/*  including:							*/
/*   sb_dump2screen      dump a matrix to the screen            */
/*   write_excel_header  write a header describing run to file  */
/*   write_bb_excel      write a single step's info to file     */
/*   write_bb_excel_sum: write a final summary info to file     */
/*   init_stat_val:      init. an array for calc. stats         */
/*   add_stat_val:       add a single value to stats array      */
/*   calc_stat_val:      calc. stats from stats array           */
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Wed Jul  1 18:38:31 bst 1998                 		*/
/****************************************************************/
/* 	MODIFIED by:						*/
/*  PDL: Aug 16, 1998						*/
/*  PDL: Sept 5, 1998						*/
/****************************************************************/
/****** To Do List **********************************************/
/*General:							*/
/* 1)                 						*/
/*   print_bb:           print out the bigblock in cap form	*/
/****************************************************************/
/*RCS Id:$Id: output_ex.c 1382 2008-09-24 14:59:34Z  $*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "machine.h"
#include "blocks.h"
#include "find_max.h"
#include "props.h"
#include "interp.h"
#include "castats.h"
#include "checks.h"
#include "pore.h"
#include "write_values.h"

extern void write_grain_histo (BB_struct * bp, int stat_flag);
extern int sb_line_int (CA_FLOAT * line_res, BB_struct * bp, int sbnum);
extern void find_minmax (p_c_list * list, int cellminmax[]);

extern void write_header_ctrl (FILE * ffp, BB_struct * bp);

/* functions defined later in file  */
void write_conc_prof (int stat_flag, BB_struct * bp, Value_struct * vp);
void write_conc_prof_multi (int stat_flag, BB_struct * bp, Value_struct * vp, int i);
void write_temp_prof (int stat_flag, BB_struct * bp);
int write_excel_header (FILE * ffp, BB_struct * bp);

/* from sb_temp_calc.c */

/********************************************************/
/********************************************************/
/* sb_dump2screen                               	*/
/*      Subroutine to write a subblock to screen.       */
/********************************************************/
/********************************************************/
int sb_dump2screen (CA_FLOAT * a, int nx, int ny, int nz)
{
  int i, j, k, index;           /* tmp counters */

/************************************************/
/* Dump it out to screen for now...             */
/************************************************/
  index = 0;
  fprintf (stderr, "\nDump of CA_FLOAT[%d,%d]\n", nx, ny);
  for (k = 0; k < nz; k++) {    /* loop sb's in z direction */
    fprintf (stderr, "\nLayer[%d]\n", k);
    for (j = 0; j < ny; j++) {  /* loop sb's in z direction */
      for (i = 0; i < nx; i++) {        /* loop cells's in z direction */
        fprintf (stderr, "%2d ", (int) (a[index++] * 10));
      }
      fprintf (stderr, "\n");
    }
  }

  return (0);
}                               /* end of sb_dump2screen subroutine */

/********************************************************/
/********************************************************/
/* write_excel_header                                  	*/
/*   Subroutine to write HEADER SUMMARY in excel format.*/
/********************************************************/
/********************************************************/
int write_excel_header (FILE * ffp, BB_struct * bp)
{
  int i;                        /* tmp counters */
  int errors = 0;               /* error flag on return */
  int ele_num, ele_1;
  MultiS_struct *ms;
  Mat_str *mp;
  Ctrl_str *cp;
  Nuc_str *np;
  CA_FLOAT part, slop, liqdiff, soldiff, cinit;

  cp = (bp->ctrl);
  ele_num = cp->NUM_COMP;
  ele_1 = ele_num - 1;
  ms = &(bp->MultiSvals);
  mp = &(bp->mprops);
  np = &(bp->nprops);

  /*fprintf(stderr,"Hello world\n"); */
  fprintf (ffp, "Hello world\n");
  /* write out base filename used and date */
  fprintf (ffp, "Excel Data for run,%s\n", bp->ctrl->fn_base);
  fprintf (ffp, "Mat file name,%s\n", bp->ctrl->fn_mat);
  fprintf (ffp, "Geo file name,%s\n", bp->ctrl->fn_geo);
  fprintf (ffp, "Porosity,%i\nSchiel:,%i\n", bp->ctrl->pore, bp->ctrl->scheil);
  fprintf (ffp, "Number of components: %d \n", bp->ctrl->NUM_COMP);
  fprintf (ffp, "\n");

  fprintf (ffp, "alloy exponent,%f", ALLOY_EXPONENT);
  /* the diffusion info */
  if (bp->ctrl->diffuse == TRUE) {
#ifdef CONST_DIFF_COEFF
    fprintf (ffp, " CONST_DIFF_COEFF used, ");
#endif
#ifdef COND_DIFF_COEFF
    fprintf (ffp, " COND_DIFF_COEFF used, ");
#endif
#ifdef AV_DIFF_COEFF
    fprintf (ffp, " AV_DIFF_COEFF used, ");
#endif
#ifdef CELL_DIFF_COEFF
    fprintf (ffp, " CELL_DIFF_COEFF used, ");
#endif
#ifdef C_LIQ_OUTPUT
    fprintf (ffp, " C_LIQ_OUTPUT used, ");
#endif
    fprintf (ffp, "\n");
  }
  if (bp->ctrl->diffuse_alloy_multi == TRUE) {
    fprintf (ffp, "material property for multi component \n");
    for (i = 0; i < ele_1; i++) {
      liqdiff = ms->LDiff_multi[i];
      soldiff = ms->SDiff_multi[i];
      cinit = ms->Cinit_multi[i];
      part = ms->part_coef_multi[i];
      slop = ms->slope_multi[i];
      fprintf (ffp, "DliqAlloy for Component %d: %e \n", i + 1, liqdiff);
      fprintf (ffp, "DsolAlloy for component %d: %e \n", i + 1, soldiff);
      fprintf (ffp, "Initial Conc.Alloy for component %d: %f \n", i + 1, cinit);
      fprintf (ffp, "partcoefAlloy for component %d: %f \n", i + 1, part);
      fprintf (ffp, "slopeAlloy for component %d: %f \n", i + 1, slop);

    }
  }
#ifdef MULTICOMP
  fprintf (ffp, " T_EUT_1,%f, T_EUT_2, %f \n", T_EUT_1, T_EUT_2);
  fprintf (ffp, " TP_1,%f, TP_2, %f, \n", TP_1, TP_2);
  fprintf (ffp, " MAX_B_CONC_1,%f , MAX_B_CONC_2, %f \n ", MAX_B_CONC_1, MAX_B_CONC_2);
#endif

  fprintf (ffp, " DAS_PRE,%f\n", DAS_PRE);
  fprintf (ffp, " DAS_EXP,%f\n", DAS_EXP);
  fprintf (ffp, " DAS_COS_THETA,%f\n", DAS_COS_THETA);
  fprintf (ffp, " ALLOY_EXPONENT,%f\n", ALLOY_EXPONENT);
  fprintf (ffp, "\n");
  fprintf (ffp, " FLAT_GROWTH,%f\n", FLAT_GROWTH);
  fprintf (ffp, " P_NUC_MIN_SS,%f\n", P_NUC_MIN_SS);
  fprintf (ffp, " P_NUC_SIG_MULT,%f\n", P_NUC_SIG_MULT);
  fprintf (ffp, "\n");
  fprintf (ffp, "Nucleation Model,%i\n", np->nmodel);

  switch (np->nmodel) {
  case N_HETERO:               /* N_HETERO model */
    fprintf (ffp, "\nHeterogenious On [Species Nmax A rad B]: ");
    for (i = 0; i < np->nhet; i++)
      fprintf (ffp, "%d, %g, %g, %g, %g; ", i, np->nuc[i].Nmax, np->nuc[i].CA_A, np->nuc[i].rad, np->nuc[i].B);
    break;
  case N_DIST:
    fprintf (ffp, "nuc_dist_func_type,%i\n", np->nucdist.nuc_dist_func_type);
    fprintf (ffp, "Function Id string,%s\n", np->nucdist.Nuc_func_id[np->nucdist.nuc_dist_func_type]);
    fprintf (ffp, "n_U_bins,%i,Unuc_incr,%.5g,n_Tund_bins,%i,Tund_incr,%.5f\n",
             np->nucdist.n_U_bins, np->nucdist.Unuc_incr, np->nucdist.n_Tund_bins, np->nucdist.Tund_incr);
    fprintf (ffp, "\nNucleation Parameters");
    fprintf (ffp, "\ngd_max, %g; ", np->gd_max);
    fprintf (ffp, "Tn, %g; ", np->NucParams[0]);
    fprintf (ffp, "Tsigma, %g; ", np->NucParams[1]);
    break;
  case N_RAPPAZ:               /*fallthru */
  case N_RAP_DEL:              /*fallthru */
  case N_OLDF_DEL:

    fprintf (ffp, "\nNucleation Parameters");
    fprintf (ffp, "\ngd_max, %g; ", np->gd_max);
    fprintf (ffp, "Tn, %g; ", np->NucParams[0]);
    fprintf (ffp, "Tsigma, %g; ", np->NucParams[1]);
    break;
  case N_BLOCK:
    fprintf (ffp, "random function type,%i", np->nucdist.nuc_dist_func_type);
    fprintf (ffp, "\nNucleation Parameters");
    fprintf (ffp, "\ngd_max, %g; ", np->gd_max);
    fprintf (ffp, "Tn, %g; ", np->NucParams[0]);
    fprintf (ffp, "Tsigma, %g; ", np->NucParams[1]);
    fprintf (ffp, "GNParam2, %g; ", np->NucParams[2]);
    fprintf (ffp, "GNParam3, %g; ", np->NucParams[3]);
    break;

  default:
    fprintf (ffp, "\nERROR:write_excel_header: ********* Nucleation Model Undefined!!!!!!!!!!!![ %i ]\n", np->nmodel);
  }

/*******************************************/
/* Control structure values                */
/*******************************************/
  write_header_ctrl (ffp, bp);
  if (bp->ctrl->pore == TRUE) {
    fprintf (ffp, "\nPorosity ON\n");
    write_pprop_values (ffp, &(bp->pprops));
  }
  fprintf (ffp, "\n\n");

  fprintf (ffp, "\n*******Header Finished*******\n");

  return (errors);
}                               /* end of write_excel_header subroutine */

/********************************************************/
/********************************************************/
/* bb_excel_step                                	*/
/*      Subroutine to write data to file to be read by  */
/* excel.                                               */
/********************************************************/
/********************************************************/
int write_bb_excel (int stat_flag, BB_struct * bp)
{
  int i, j, k, index, err;      /* tmp counters */
  int poreminmax[6];
  Ctrl_str *cp;
  MultiS_struct *ms;
  int rtn_flag = 0, NTpores = 0;        /* error flag on return, total pores */
  int sbnum;                    /*subblock number */
  int ele_1, ele_num;
  Value_struct *local_v_struct;
  Value_struct **multi_conc_ptr;
  int t_index;

  char command[MAX_STRING_LEN];
  CA_FLOAT Nstat[MAX_NSP][4];
  CA_FLOAT *t_list;

  CA_FLOAT Tstat[4], fs_stat[4];
  CA_FLOAT statlist[N_T_LISTS][PORE_NTRAD][4], val;
  CA_FLOAT pctporlist[PORE_NTRAD];
  int count[N_T_LISTS][PORE_NTRAD];
  FILE *fd_ex, *fd_ex_pore, *fd_ex_pore_extent, *fd_ex_pore_stats;
  PORE_str *pl;
  T_struct *Tp;
  time_t realtime;
  char open_action[5];

  ms = &(bp->MultiSvals);
  cp = bp->ctrl;
  ele_num = cp->NUM_COMP;
  ele_1 = ele_num - 1;

  strcpy (open_action, "w");
/********************************************************/
/* if Flag = FIRST_CALL open file and write headers     */
/********************************************************/
  switch (stat_flag) {

  case RESTART_CALL:
    strcpy (open_action, "a");
  case FIRST_CALL:             /* first call to write_bb_excel, do all initilisation */
    /* open output file */
    bp->scr_dump_num = 0;
    sprintf (command, "%s.csv", bp->ctrl->fn_base);
    if ((bp->ctrl->fd_ex = fopen (command, open_action)) == NULL) {
      fprintf (stderr, "Error: can't open output file [%s]\n", command);
      rtn_flag = 1;
    } else {
      fd_ex = bp->ctrl->fd_ex;
      write_excel_header (fd_ex, bp);
      fprintf (fd_ex, "Realtime,");
      fprintf (fd_ex, "Step,Time");

      /* header for the nucleation model values */
      fprintf (fd_ex, ",GrainNcell,GrainNgrow,BB ngrains");
      if (bp->nprops.nmodel == N_HETERO) {
        for (i = 0; i < bp->nprops.nhet; i++)
          fprintf (fd_ex, ",Ntot[%d]", i);
      }

      /* header for the values per sb [T, ngr, fs] */
      for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {      /* write header for each sb */
        if (bp->sb[sbnum]->open) {
          fprintf (fd_ex, ",Temp[%d],ngrains[%d],fs[%d],fscheck,ngrowchk", sbnum, sbnum, sbnum);
          if (bp->ctrl->diffuse == TRUE) {
            fprintf (fd_ex, ",Ctot[%d],Cmax[%d],Cmin[%d]", sbnum, sbnum, sbnum);
            fprintf (fd_ex, ",SATmax[%d],SSmax[%d]", sbnum, sbnum);
            fprintf (fd_ex, ",CLmax[%d],CLmin[%d]", sbnum, sbnum);
          }
          if ((bp->ctrl->diffuse_alloy == TRUE) && (bp->ctrl->diffuse_alloy_multi == FALSE)) {
            fprintf (fd_ex, ",Atot[%d],Amax[%d],Amin[%d]", sbnum, sbnum, sbnum);
            fprintf (fd_ex, ",ALmax[%d],ALmin[%d]", sbnum, sbnum);
          }
          if (bp->ctrl->mould_src != 0) {
            fprintf (fd_ex, ",alloyaddsol[%d],alloyt_addsol[%d]", sbnum, sbnum);
          }
          if (bp->ctrl->mould_src != 0) {
            fprintf (fd_ex, ",gasaddsol[%d],gast_addsol[%d]", sbnum, sbnum);
          }
          if (bp->ctrl->pore == TRUE)
            fprintf (fd_ex, ",Npores[%d]/[%d],pore_nmols[%d]", sbnum, bp->sb[sbnum]->Npores, sbnum);
        }
      }
      fprintf (fd_ex, "\n");
    }

    /* output the frac solid info */
    write_conc_prof (stat_flag, bp, bp->c_fs_values);
    write_temp_prof (stat_flag, bp);

    /* output the gas info */
    if (bp->ctrl->diffuse == TRUE) {
      write_conc_prof (stat_flag, bp, bp->c_sol_values);
    }

    /* output the alloy info */
    if (bp->ctrl->diffuse_alloy == TRUE && bp->ctrl->diffuse_alloy_multi == FALSE) {
      write_conc_prof (stat_flag, bp, bp->c_sol_alloy_values);
    }


    /* output the pore info to seperate file */
    if (bp->ctrl->pore == TRUE) {
      /*open the main pore info file */
      sprintf (command, "P%s.csv", bp->ctrl->fn_base);
      if ((bp->ctrl->fd_ex_pore = fopen (command, open_action)) == NULL) {
        fprintf (stderr, "Error: can't open pore output file [%s]\n", command);
        rtn_flag = 1;
      } else {
        fd_ex_pore = bp->ctrl->fd_ex_pore;
        write_excel_header (fd_ex_pore, bp);
        fprintf (stderr, "Writing pore header\n");
        fprintf (fd_ex_pore, "Subblock:");

        for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {    /* write header for each sb */
          if (bp->sb[sbnum]->open) {

            NTpores += bp->sb[sbnum]->Npores;
            fprintf (fd_ex_pore, ",%i", sbnum);
          }
        }
        fprintf (fd_ex_pore, ",TotalPores,%i\n", NTpores);
        fprintf (fd_ex_pore, "MaxNpores: %i", sbnum);
        for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {    /* write header for each sb */
          if (bp->sb[sbnum]->open) {
            NTpores += bp->sb[sbnum]->Npores;
            fprintf (fd_ex_pore, ",%i", bp->sb[sbnum]->Npores);
          }
        }
        fprintf (fd_ex_pore, "\n");
        fprintf (fd_ex_pore, "\n");
      }
      /*open the pore min/max extent file */
      sprintf (command, "PE%s.csv", bp->ctrl->fn_base);
      if ((bp->ctrl->fd_ex_pore_extent = fopen (command, open_action)) == NULL) {
        fprintf (stderr, "Error: can't open pore E output file [%s]\n", command);
        rtn_flag = 1;
      } else {
        fd_ex_pore_extent = bp->ctrl->fd_ex_pore_extent;
        write_excel_header (fd_ex_pore_extent, bp);
        fprintf (stderr, "Writing pore header\n");
        fprintf (fd_ex_pore_extent, "Subblock:");

        for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {    /* write header for each sb */
          if (bp->sb[sbnum]->open) {
            NTpores += bp->sb[sbnum]->Npores;
            fprintf (fd_ex_pore, ",%i", sbnum);
          }
        }
        fprintf (fd_ex_pore_extent, ",TotalPores,%i\n", NTpores);
        fprintf (fd_ex_pore_extent, "MaxNpores: %i", sbnum);
        for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {    /* write header for each sb */
          if (bp->sb[sbnum]->open) {
            NTpores += bp->sb[sbnum]->Npores;
            fprintf (fd_ex_pore_extent, ",%i", bp->sb[sbnum]->Npores);
          }
        }
        fprintf (fd_ex_pore_extent, "\n");
        fprintf (fd_ex_pore_extent, "\n");
        fprintf (fd_ex_pore_extent, "Pore Extent:\n");
        fprintf (fd_ex_pore_extent, "PoreNum,MinX,MinY,MinZ,MaxX,MaxY,MaxZ\n");
      }
      /*open the pore stat only file */
      sprintf (command, "PS%s.csv", bp->ctrl->fn_base);
      if ((bp->ctrl->fd_ex_pore_stats = fopen (command, open_action)) == NULL) {
        fprintf (stderr, "Error: can't open pore S output file [%s]\n", command);
        rtn_flag = 1;
      } else {
        fd_ex_pore_stats = bp->ctrl->fd_ex_pore_stats;
        write_excel_header (fd_ex_pore_stats, bp);
        fprintf (stderr, "Writing pore header\n");
        fprintf (fd_ex_pore_stats, "Pore statistics\n");
      }

    }                           /* end of write pore info */
    write_grain_histo (bp, stat_flag);
    break;                      /* end of FIRST_CALL */

/********************************************************/
/* if Flag = STEP_CALL write out one timestep's info    */
/********************************************************/
  case STEP_CALL:
    fd_ex = bp->ctrl->fd_ex;
    realtime = time (0);
    /* loop through all sb's generating average values... */
    if (bp->ntsb == 1) {        /* if only one, use sb[0] values */
      Tstat[0] = Tstat[1] = Tstat[2] = bp->sb[0]->Tvals.Tavg;
      Tstat[3] = 0.0;
      fs_stat[0] = fs_stat[1] = fs_stat[2] = bp->sb[0]->Tvals.fsavg;
      fs_stat[3] = 0.0;
    } else {                    /* more than one, so loop it!    */
      init_stat_val (Tstat);
      init_stat_val (fs_stat);
      for (i = 0; i < bp->ntsb; i++) {  /* calc stats for each sb      */
        if (bp->sb[i]->open) {
          add_stat_val (Tstat, bp->sb[i]->Tvals.Tavg);
          add_stat_val (fs_stat, bp->sb[i]->Tvals.fsavg);
        }
      }
      calc_stat_val (Tstat, (CA_FLOAT) bp->ntsb);
      calc_stat_val (fs_stat, (CA_FLOAT) bp->ntsb);
    }                           /* end of stat generating loop */

      /******************************************************/
    /* Write out all values on one line for this timestep */
      /******************************************************/
    /* write out the BB stat values */
    fprintf (fd_ex, "%.10g,", (double) realtime);
    fprintf (fd_ex, "%d,%.10g,%.10g,%.10g", bp->step, bp->sim_time, grain_ncells (bp), grain_ngrow (bp));

    /* write out the nucleation model values */
    fprintf (fd_ex, ",%d", bp->nprops.ngr);
    if (bp->nprops.nmodel == N_HETERO) {
      for (i = 0; i < bp->nprops.nhet; i++)
        fprintf (fd_ex, ",%g", bp->nprops.nuc[i].Ntot);
    }

    /*write out frac_solid profile to seperate excel-type file */
    write_conc_prof (stat_flag, bp, bp->c_fs_values);

    /*write out concentration profile to seperate excel-type file */
    if (bp->ctrl->diffuse == TRUE) {
      write_conc_prof (stat_flag, bp, bp->c_sol_values);
    }

    /*write out concentration profile to seperate excel-type file */
    if ((bp->ctrl->diffuse_alloy == TRUE) && (bp->ctrl->diffuse_alloy_multi == FALSE)) {
      write_conc_prof (stat_flag, bp, bp->c_sol_alloy_values);
    }


    if (bp->ctrl->pore == TRUE) {
      /* collate the pore radius values */
      sbnum = bp->ctrl->pore_dump_sb /*=0*/ ;
      k = bp->scr_dump_num;
      fprintf (stderr, "  screendump # %i\n", k);
      if (bp->sb[sbnum]->open) {
        pl = bp->sb[sbnum]->porelist;
        if (k > bp->ctrl->nscrdumps) {
          fprintf (stderr, "ERROR: output_ex: Too many pore dumps! [%i]\n", k);
        } else {
          for (i = 0; i < bp->sb[sbnum]->Npores; i++) {
            if (pl[i].State != PORE_NONE) {
#ifdef MORE_PORE
              pl[i].radlist[k] = pl[i].Radius;
              pl[i].nmollist[k] = pl[i].Nmols;
              pl[i].templist[k] = pl[i].Temp;
#endif
            }                   /* end of if pore > 0 */
          }                     /* end of looping through pores */
        }
      }
    }

    /* end of if pore = true */
    /* append the histogram of grain sizes */
    write_grain_histo (bp, stat_flag);

    /* write out the values per sb [T, ngr, fs,npores] */
    for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {        /* write values for each sb */
      if (bp->sb[sbnum]->open) {

        /* print T, ngr, fs */
        Tp = &(bp->sb[sbnum]->Tvals);
        fprintf (fd_ex, ",%.10g,%d,%.10g,%.10g,%i", Tp->Tavg, bp->sb[sbnum]->ngr, Tp->fsavg, find_fs_av (bp, sbnum),
                 find_ngrow (bp, sbnum));
        if (bp->ctrl->diffuse == TRUE) {
          /*find the maximum and the total in the conc. array */
          err = find_max_conc (bp, bp->c_sol_values, sbnum);
          if (err != 0)
            fprintf (stderr, "ERROR: output_ex: step %i: %d negative gas concentration values!", bp->step, err);
          fprintf (fd_ex, ",%g,%g,%g", ((bp->vol_c) * (bp->c_sol_values->Ctot)), bp->c_sol_values->Cmax, bp->c_sol_values->Cmin);
          fprintf (fd_ex, ",%g,%g", bp->c_sol_values->SATmax, bp->c_sol_values->SSmax);
          fprintf (fd_ex, ",%g,%g", bp->c_sol_values->CLmax, bp->c_sol_values->CLmin);
        }                       /*end diffuse (gas) */
        if ((bp->ctrl->diffuse_alloy == TRUE) && (bp->ctrl->diffuse_alloy_multi == FALSE)) {
          /*find the maximum and the total in the conc. array */
          err = find_max_conc (bp, bp->c_sol_alloy_values, sbnum);
          if (err != 0)
            fprintf (stderr, "ERROR: output_ex: step %i:  %d negative alloy concentration values!", bp->step, err);
          fprintf (fd_ex, ",%g,%g,%g", bp->c_sol_alloy_values->Ctot, bp->c_sol_alloy_values->Cmax, bp->c_sol_alloy_values->Cmin);
          fprintf (fd_ex, ",%g,%g", bp->c_sol_alloy_values->CLmax, bp->c_sol_alloy_values->CLmin);
        }                       /*end diffuse alloy */
        if (bp->ctrl->mould_src != 0) {
          fprintf (fd_ex, ",%g,%g", bp->sb[sbnum]->Svals[ALLOY].addsol, bp->sb[sbnum]->Svals[ALLOY].t_addsol);
          fprintf (fd_ex, ",%g,%g", bp->sb[sbnum]->Svals[GAS], bp->sb[sbnum]->Svals[GAS].t_addsol);
        }

        /* and the amount of H2 in all pores */
        if (bp->ctrl->pore == TRUE) {
          CA_FLOAT pore_nmols = 0;

          pl = bp->sb[sbnum]->porelist;
          fprintf (fd_ex, ",%d", bp->sb[sbnum]->sb_npores);
          /* add up the n-mols in the pores */
          for (i = 0; i < bp->sb[sbnum]->Npores; i++) {
            switch (pl[i].State) {
            case PORE_OFF:     /*fallthrough */
            case PORE_NONE:    /*fallthrough */
            case PORE_LATENT:
            case NOT_CASTING:
              /* do nothing */
              break;

            case PORE_SPHERE:  /*fallthrough */
            case PORE_TUBE:
            case PORE_FROZEN:
            case PORE_MULTI:
              pore_nmols += pl[i].NmolsH;
              break;

            default:
              fprintf (stderr, "Output_ex: Pore state in hyperspace %i\n", pl[i].State);
              break;
            }                   /* end of pore switch */
          }                     /* end of looping through pores */
          fprintf (fd_ex, ",%.10g", pore_nmols);
#ifdef ONEPORE
          fprintf (fd_ex, ",%.10g", pl[0].EqRad);
#endif
        }
      }                         /*end check if block is open */
    }                           /*end loop trhough subblocksl */

    fprintf (fd_ex, "\n");
    fflush (NULL);
    bp->scr_dump_num++;
    break;

/********************************************************/
/* if Flag = LAST_CALL write out one timestep's info    */
/********************************************************/
  case LAST_CALL:

    /*write out frac solid profile to seperate excel-type file */
    write_conc_prof (stat_flag, bp, bp->c_fs_values);

    /*write out concentration profile to seperate excel-type file */
    if (bp->ctrl->diffuse == TRUE) {
      write_conc_prof (stat_flag, bp, bp->c_sol_values);
    }

    /*write out concentration profile to seperate excel-type file */
    if ((bp->ctrl->diffuse_alloy == TRUE) && (bp->ctrl->diffuse_alloy_multi == FALSE)) {
      write_conc_prof (stat_flag, bp, bp->c_sol_alloy_values);
    }


    /* write out the pore list */
    if (bp->ctrl->pore == TRUE) {
      fd_ex_pore_extent = bp->ctrl->fd_ex_pore_extent;
      sbnum = bp->ctrl->pore_dump_sb;
      if (bp->sb[sbnum]->open) {
        k = bp->scr_dump_num;
        pl = bp->sb[sbnum]->porelist;

       /*******************************************/
        /* the min,max extent of the pore to a seperate file */
       /*******************************************/
        for (i = 0; i < bp->sb[sbnum]->Npores; i++) {
          int ii;

          switch (pl[i].State) {
          case PORE_NONE:      /*fallthrough */
          case PORE_LATENT:
          case NOT_CASTING:
            /* do nothing */
            break;

          case PORE_OFF:       /*fallthrough */
          case PORE_SPHERE:    /*fallthrough */
          case PORE_TUBE:
          case PORE_MULTI:
          case PORE_FROZEN:
            fprintf (fd_ex_pore_extent, "%i", i);
            find_minmax (pl[i].boundary, poreminmax);
            for (ii = 0; ii < 6; ii++) {
              fprintf (fd_ex_pore_extent, ",%i", poreminmax[ii]);
            }
            fprintf (fd_ex_pore_extent, "\n");
            break;

          default:
            fprintf (stderr, "Output_ex: Pore state in hyperspace %i\n", pl[i].State);
            break;

          }                     /* end existing pore test */
        }                       /* end loop to find extent of each pore */

        /* produce trad list */
        fprintf (bp->ctrl->fd_ex_pore, "\n");
        fprintf (bp->ctrl->fd_ex_pore, "Pore Data lists:\n");
        /*******************************************/
        /* loop to print out all defined t_lists   */
        /*******************************************/
        for (t_index = 0; t_index < N_T_LISTS; t_index++) {     /*see pore.h */

          /*write the temperature line */
          fprintf (bp->ctrl->fd_ex_pore, "Temp:,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x");
          for (j = 0; j < PORE_NTRAD; j++) {
            fprintf (bp->ctrl->fd_ex_pore, ",%.10e", bp->pprops.P_trad_max - j * bp->pprops.P_trad_step);
          }
          fprintf (bp->ctrl->fd_ex_pore, "\n");
          /* initialize the stat values list for this t_list */
          for (j = 0; j < PORE_NTRAD; j++) {
            init_stat_val (statlist[t_index][j]);
            count[t_index][j] = 0;
          }
          /* process this data for all pores and output */
          for (i = 0; i < bp->sb[sbnum]->Npores; i++) {
            /* process existing pores */
            switch (pl[i].State) {
            case PORE_NONE:    /*fallthrough */
            case PORE_LATENT:
            case NOT_CASTING:
              /* do nothing */
              break;

            case PORE_OFF:     /*fallthrough */
            case PORE_SPHERE:  /*fallthrough */
            case PORE_TUBE:
            case PORE_MULTI:
            case PORE_FROZEN:
              t_list = pl[i].t_lists[t_index];  /* local version of current data t_list */

              /* row headers */
#                     ifndef STAT_ONLY  /* flag to limit amount of data to stat values only - for big runs */
              fprintf (bp->ctrl->fd_ex_pore, "%s:,%i", t_listnames[t_index], i);
              fprintf (bp->ctrl->fd_ex_pore, ",State:,%i", pl[i].State);
              fprintf (bp->ctrl->fd_ex_pore, ",Thr:,%.10g", pl[i].Thresh);
              fprintf (bp->ctrl->fd_ex_pore, ",NucSS,%.10g", pl[i].NucSS);
              fprintf (bp->ctrl->fd_ex_pore, ",NucTemp,%.10g", pl[i].NucTemp);
              fprintf (bp->ctrl->fd_ex_pore, ",PoreNcells:,%i", pl[i].ncells);
              fprintf (bp->ctrl->fd_ex_pore, ",PoreXYZ:,%i,%i,%i", pl[i].Cell[0], pl[i].Cell[1], pl[i].Cell[2]);
#                     endif /*STAT_ONLY */

              /* use the stat_values routines and print values */
              for (j = 0; j < PORE_NTRAD; j++) {
                val = t_list[j];
                if (val != 0.0) {
                  add_stat_val (statlist[t_index][j], val);
                  count[t_index][j]++;
                }
#                        ifndef STAT_ONLY       /* flag to limit amount of data to stat values only - for big runs */
                fprintf (bp->ctrl->fd_ex_pore, ",%.10e", val);
#                        endif /*STAT_ONLY */
              }                 /* end process current list (index j loop) */
#                     ifndef STAT_ONLY  /* flag to limit amount of data to stat values only - for big runs */
              fprintf (bp->ctrl->fd_ex_pore, "\n");
#                     endif /*STAT_ONLY */
              break;

            default:
              fprintf (stderr, "Output_ex: Pore state in hyperspace %i\n", pl[i].State);
              break;

            }                   /* end existing pore test */
          }                     /*end list out list for current pore (index i loop) */
#          ifndef STAT_ONLY     /* flag to limit amount of data to stat values only - for big runs */
          fprintf (bp->ctrl->fd_ex_pore, "\n");
#          endif /*STAT_ONLY */
        }                       /* end of t_list output loop (t_index loop) */

        fprintf (bp->ctrl->fd_ex_pore, "\n");

        /* caculate percent porosity from vol sums */
        /* Have to do this before calc-stat-val    */
        /* while statlist has the sum instead of the mean */
        for (j = 0; j < PORE_NTRAD; j++) {
          pctporlist[j] = statlist[T_VOL_LIST][j][0] / bp->vol_sb;
        }
        /* list outthe statistical values */
        for (t_index = 0; t_index < N_T_LISTS; t_index++) {

          fprintf (bp->ctrl->fd_ex_pore, "\n");
          fprintf (bp->ctrl->fd_ex_pore, "%s stat values:", t_listnames[t_index]);
          fprintf (bp->ctrl->fd_ex_pore_stats, "\n");
          fprintf (bp->ctrl->fd_ex_pore_stats, "%s stat values:", t_listnames[t_index]);
          for (j = 0; j < PORE_NTRAD; j++) {
            calc_stat_val (statlist[t_index][j], (CA_FLOAT) count[t_index][j]);
          }

          fprintf (bp->ctrl->fd_ex_pore, "\n");
          fprintf (bp->ctrl->fd_ex_pore, "Temp:,x,x,x,x,x,x,x,x,");
          fprintf (bp->ctrl->fd_ex_pore_stats, "\n");
          fprintf (bp->ctrl->fd_ex_pore_stats, "Temp:,x,x,x,x,x,x,x,x,");
          for (j = 0; j < PORE_NTRAD; j++) {
            fprintf (bp->ctrl->fd_ex_pore, ",%.10e", bp->pprops.P_trad_max - j * bp->pprops.P_trad_step);
            fprintf (bp->ctrl->fd_ex_pore_stats, ",%.10e", bp->pprops.P_trad_max - j * bp->pprops.P_trad_step);
          }
          /* loop through stat value array for each value */
          for (i = 0; i < 4; i++) {
            fprintf (bp->ctrl->fd_ex_pore, "\n");
            fprintf (bp->ctrl->fd_ex_pore, "statval %i:,x,x,x,x,x,x,x,x,", i);
            fprintf (bp->ctrl->fd_ex_pore_stats, "\n");
            fprintf (bp->ctrl->fd_ex_pore_stats, "statval %i:,x,x,x,x,x,x,x,x,", i);
            for (j = 0; j < PORE_NTRAD; j++) {
              fprintf (bp->ctrl->fd_ex_pore, ",%.10e", statlist[t_index][j][i]);
              fprintf (bp->ctrl->fd_ex_pore_stats, ",%.10e", statlist[t_index][j][i]);
            }
          }
          fprintf (bp->ctrl->fd_ex_pore, "\n");
          fprintf (bp->ctrl->fd_ex_pore, "count:,x,x,x,x,x,x,x,x,");
          fprintf (bp->ctrl->fd_ex_pore_stats, "\n");
          fprintf (bp->ctrl->fd_ex_pore_stats, "count:,x,x,x,x,x,x,x,x,");
          for (j = 0; j < PORE_NTRAD; j++) {
            fprintf (bp->ctrl->fd_ex_pore, ",%i", count[t_index][j]);
            fprintf (bp->ctrl->fd_ex_pore_stats, ",%i", count[t_index][j]);
          }
          fprintf (bp->ctrl->fd_ex_pore, "\n");
          fprintf (bp->ctrl->fd_ex_pore_stats, "\n");
        }                       /* end list out all statistical values */
        /* list out the percent porostiy */
        fprintf (bp->ctrl->fd_ex_pore, "pct por:,x,x,x,x,x,x,x,x,");
        fprintf (bp->ctrl->fd_ex_pore_stats, "pct por:,x,x,x,x,x,x,x,x,");
        for (j = 0; j < PORE_NTRAD; j++) {
          fprintf (bp->ctrl->fd_ex_pore, ",%.10e", pctporlist[j]);
          fprintf (bp->ctrl->fd_ex_pore_stats, ",%.10e", pctporlist[j]);
        }
        fprintf (bp->ctrl->fd_ex_pore, "\n");
        fprintf (bp->ctrl->fd_ex_pore_stats, "\n");

        /*add a recognizable end string to aid with post=processing */
        fprintf (bp->ctrl->fd_ex_pore, "EnDeNd EnDeNd\n");
        fprintf (bp->ctrl->fd_ex_pore_stats, "EnDeNd EnDeNd\n");
        fprintf (bp->ctrl->fd_ex_pore_extent, "EnDeNd EnDeNd\n");

        fclose (bp->ctrl->fd_ex_pore);
        fclose (bp->ctrl->fd_ex_pore_stats);
        fclose (bp->ctrl->fd_ex_pore_extent);
      }                         /* end check if sb is open */
    }
    /* end if  pore mode on */
    fclose (bp->ctrl->fd_ex);
    break;

  default:
    fprintf (stderr, "ERROR:output_ex: stat_flag value [%d] undefined.\n", stat_flag);
    break;
  }                             /* end of Flag switch */

  return (rtn_flag);
}                               /* end of write_bb_excel subroutine */

/********************************************************/
/********************************************************/
/* write_bb_excel_sum                                  	*/
/*    Subroutine to write SUMMARY table in excel format.*/
/********************************************************/
/********************************************************/
int write_bb_excel_sum (int stat_flag, BB_struct * bp)
{
  int i, j, k, index;           /* tmp counters */
  int rtn_flag = 0;             /* error flag on return */
  char command[MAX_STRING_LEN];
  int *c_per_gr, *gr, *zerograin;
  CA_FLOAT tmp, gr_stat[4], gr_Rstat[4];
  FILE *fd;
  Ind_grain *c_gr;

/********************************************************/
/* open the output file                                 */
/********************************************************/
  /* open output file */
  sprintf (command, "%s_sum.csv", bp->ctrl->fn_base);
  if ((fd = fopen (command, "w")) == NULL) {
    fprintf (stderr, "Error: can't open output file [%s]\n", command);
    rtn_flag = 1;
  } else {
    fprintf (fd, "test\n");

/********************************************************/
/* malloc array to hold grain # of cells and set = 0    */
/********************************************************/
    if (!(c_per_gr = (int *) calloc (1 + bp->nprops.ngr, sizeof (int)))) {
      fprintf (stderr, "ERROR: gr size counting array malloc failed\n");
      return (1);
    }
    if (!(zerograin = (int *) calloc (bp->ncsb, sizeof (int)))) {
      fprintf (stderr, "ERROR: zerogr array malloc failed\n");
      return (1);
    }

/********************************************************/
/* for each sb loop through and calc average gs and     */
/* write it out to the file ...                         */
/********************************************************/
    for (i = 0; i < bp->ntsb; i++) {    /*loop through all sb */
      if (bp->sb[i]->open == TRUE) {
        gr = bp->sb[i]->gr;
      } else {
        gr = zerograin;
      }
      for (j = 0; j < bp->ncsb; j++) {  /*loop thru all cells in sb */
        if (*gr >= 0) {
          c_per_gr[*(gr++)]++;
        }
      }
    }

/********************************************************/
/* calculate stats on grain sizes                       */
/********************************************************/
    init_stat_val (gr_stat);
    for (i = 0; i < bp->nprops.ngr; i++) {      /*write header for each sb */
      add_stat_val (gr_stat, (CA_FLOAT) c_per_gr[i]);
    }
    calc_stat_val (gr_stat, (CA_FLOAT) bp->nprops.ngr);
    for (i = 0; i < 4; i++) {   /* convert to mm^3 and microns */
      gr_stat[i] = bp->vol_c * gr_stat[i];
      tmp = gr_stat[i] * THREE_BY_4PI;
      gr_Rstat[i] = 1000.0 * POW (tmp, (CA_FLOAT) (1.0 / 3.0));
    }

/********************************************************/
/* write out the file                                   */
/********************************************************/
    fprintf (fd, "test\n");
    write_excel_header (fd, bp);
    if (bp->ctrl->diffuse_alloy_multi == TRUE) {
      fprintf (fd, "number of cells with eutectic phase is: %d  \n", bp->nbeut);
      fprintf (fd, "number of cells with ternary phase is: %d  \n", bp->nteut);
    }
    fprintf (fd, "Volume in m^3, Radius in Millimeters\n");
    fprintf (fd, "Ngr,Avg V,Min V,Max V,Sigma V");
    fprintf (fd, ",Avg R,Min R,Max R,Sigma R");
    fprintf (fd, ",LineIntGrSz\n");
    fprintf (fd, "%d", bp->nprops.ngr);
    for (i = 0; i < 4; i++)
      fprintf (fd, ",%g", gr_stat[i]);
    for (i = 0; i < 4; i++)
      fprintf (fd, ",%g", gr_Rstat[i]);
    for (i = 0; i < bp->ntsb; i++) {
      CA_FLOAT lineres;

      if (sb_line_int (&lineres, bp, i) == 0) {
        fprintf (fd, ",%.5g", lineres);
      } else {
        fprintf (fd, ",ERROR");
      }
    }
    fprintf (fd, "\n");
    fprintf (fd, "\nGrain,Ncells,NucTh,TNuc,TunderNuc,CellConcNuc,");
    fprintf (fd, "Xmax,Ymax,Zmax,Xmin,Ymin,Zmin,");
    fprintf (fd, "angle0,angle1,angle2\n");
    for (i = 1; i < bp->nprops.ngr; i++) {      /*write out the # cells/grain */
      c_gr = *(bp->gr + i);
      fprintf (fd, "%d,%d,%.10g,%.10g,%.10g,%.10g,", i, c_per_gr[i], c_gr->NucTh, c_gr->TNuc, c_gr->TunderNuc, c_gr->CellConcNuc);
      for (j = 0; j < 3; j++) {
        fprintf (fd, "%d,", c_gr->max[j]);
      }
      for (j = 0; j < 3; j++) {
        fprintf (fd, "%d,", c_gr->min[j]);
      }
      for (j = 0; j < 3; j++) {
	 fprintf(fd,"%.5g,",c_gr->ang[j]);
      }
      fprintf(fd,"\n");
    }


    if (bp->ctrl->pore == TRUE) {
      for (j = 0; j < bp->ntsb; j++) {
        fprintf (fd, "Pores:Subblock %i\n", j);
        fprintf (fd, "Cellnum,Radius,Dwelltime\n");
        for (i = 0; i < bp->sb[j]->Npores; i++) {
          if (!(bp->sb[j]->porelist[i].State == PORE_NONE || bp->sb[j]->porelist[i].State == PORE_LATENT)) {
            fprintf (fd, "%d,%g,%g\n", bp->sb[j]->porelist[i].Cellnum, bp->sb[j]->porelist[i].Radius, bp->sb[j]->porelist[i].Time);
          }
        }
      }
    }
    fclose (fd);
    free (c_per_gr);
    free (zerograin);
  }
  return (rtn_flag);
}                               /* end of write_bb_excel_sum subroutine */

void write_temp_prof (int stat_flag, BB_struct * bp)
{
}

/********************************************************/
/********************************************************/
/* write_conc_prof                                     	*/
/*    write a concentration profile to excel file        */
/********************************************************/
/********************************************************/
void write_conc_prof (int stat_flag, BB_struct * bp, Value_struct * output_v)
{
  if (!(bp->ctrl->do_conc_prof)) {
    return;                     /* return immediately if the flag for concentration profiles is not set */
  }

  {                             /* otherwise do the subroutine */
    int start, i, *conc_prof, sbnum;
    char command[MAX_STRING_LEN];

#ifdef PRINT_AV_PROFILE
    char avprofname[MAX_STRING_LEN];
    FILE *avfp;                 /*av profile file pointer */
#endif /*PRINT_AV_PROFILE */
    static int msg[3] = { 0, 0, 0 };
    CA_FLOAT *fs, *fsp, *sol, *sp, partcoef;
    FILE *fp;                   /*file pointer */
    char open_action[5];

    strcpy (open_action, "w");
/********************************************************/
/* if Flag = FIRST_CALL open file and write headers     */
/********************************************************/
    switch (stat_flag) {

    case RESTART_CALL:
      strcpy (open_action, "a");
      /* fall through */
    case FIRST_CALL:
      sprintf (command, "%s%s.txt", output_v->id_string, bp->ctrl->fn_base);
      if ((output_v->fd_exel = fopen (command, open_action)) == NULL) {
        fprintf (stderr, "Error: can't open %s. output file [%s]\n", output_v->id_string, command);
        return;
      } else {
        fp = output_v->fd_exel;
        fprintf (fp, "%s CONCENTRATION OUTPUT FOR MULTICOMP \n", output_v->id_string);
        write_excel_header (fp, bp);
        fprintf (fp, "Step,Time,");
        for (i = 0; i < bp->nc[0]; i++) {
          fprintf (fp, "%f,", i * bp->size_c[0]);
        }
        fprintf (fp, "\n");
      }
#ifdef PRINT_AV_PROFILE
      sprintf (avprofname, "%sAV%s.txt", output_v->id_string, bp->ctrl->fn_base);
      if ((avfp = fopen (avprofname, open_action)) == NULL) {
        fprintf (stderr, "Error: can't open %s. output file [%s]\n", output_v->id_string, command);
        return;
      } else {
        fprintf (avfp, "%s CONCENTRATION OUTPUT\n", output_v->id_string);
        write_excel_header (avfp, bp);
        fprintf (avfp, "Step,Time,");
        for (i = 0; i < bp->nc[0]; i++) {
          fprintf (avfp, "%f,", i * bp->size_c[0]);
        }
        fprintf (avfp, "\n");
      }
      fclose (avfp);
      if (strcmp ("G_", output_v->id_string) == 0) {
        /* do a supersat. profile also */
        sprintf (avprofname, "SS_%sAV%s.txt", output_v->id_string, bp->ctrl->fn_base);
        if ((avfp = fopen (avprofname, open_action)) == NULL) {
          fprintf (stderr, "Error: can't open %s. output file [%s]\n", output_v->id_string, command);
          return;
        } else {
          fprintf (avfp, "%s SUPERSATURATION OUTPUT\n", output_v->id_string);
          write_excel_header (avfp, bp);
          fprintf (avfp, "Step,Time,");
          for (i = 0; i < bp->nc[0]; i++) {
            fprintf (avfp, "%f,", i * bp->size_c[0]);
          }
          fprintf (avfp, "\n");
        }
        fclose (avfp);

      }                         /* end supersat profile bit */
#endif /*PRINT_AV_PROFILE */
      break;
    case STEP_CALL:
      conc_prof = bp->ctrl->conc_prof;
      sbnum = conc_prof[0];
      if (bp->sb[sbnum]->open == FALSE) /*sb not active yet */
        return;

      partcoef = output_v->part_coef;
      fp = output_v->fd_exel;
      /*fprintf(fp,"step,%i,time,%.5g",bp->step,bp->sim_time); */
      if (sbnum > bp->ntsb && !msg[0]) {
        fprintf (stderr, "ERROR: ConcProfile sb %i out of range\n", sbnum);
        fprintf (fp, "\nERROR: ConcProfile sb %i out of range\n", sbnum);
        msg[0] = TRUE;
        return;
      }
      if (conc_prof[1] >= bp->nc[2] && !msg[2]) {
        fprintf (stderr, "ERROR: ConcProfile slice %i out of range\n", conc_prof[1]);
        fprintf (fp, "\nERROR: ConcProfile slice %i out of range\n", conc_prof[1]);
        msg[2] = TRUE;
        return;
      }
      if (conc_prof[2] > bp->nc[1] && !msg[1]) {
        fprintf (stderr, "ERROR: ConcProfile row %i out of range\n", conc_prof[2]);
        fprintf (fp, "\nERROR: ConcProfile row %i out of range\n", conc_prof[2]);
        msg[1] = TRUE;
        return;
      }
      if (msg[0] || msg[1] || msg[2])
        return;
      sol = output_v->block_array[sbnum];
      if (bp->ctrl->scheil == TRUE) {
        fs = bp->sb[sbnum]->sch_fs;
      } else {
        fs = bp->sb[sbnum]->c_fs;
      }
      start = bp->nc[0] * bp->nc[1] * conc_prof[1] + bp->nc[0] * conc_prof[2];
      sp = sol + start;
      fsp = fs + start;
      fprintf (fp, "%i,%f,", bp->step, bp->sim_time);
      for (i = 0; i < bp->nc[0]; i++) {
        if (*fsp == NOT_CASTING) {
          fprintf (fp, "-1 ,");
        } else {
          fprintf (fp, "%.10f,", *sp);
        }
        sp++;
        fsp++;
      }
      fprintf (fp, "\n");
#ifdef PRINT_AV_PROFILE
      {
        int jj, kk;
        CA_FLOAT av = 0, ss = 0, gas = 0, alloy = 0, fs = 0, cell_temp = 0;
        int num = 0;
        int index_ca;

        sprintf (avprofname, "%sAV%s.txt", output_v->id_string, bp->ctrl->fn_base);
        if ((avfp = fopen (avprofname, "a")) == NULL) {
          fprintf (stderr, "Error: can't open %s. output file [%s]\n", output_v->id_string, command);
          return;
        }

        fprintf (avfp, "%i,%f,", bp->step, bp->sim_time);
        for (i = 0; i < bp->nc[0]; i++) {
          av = 0;
          num = 0;
          for (jj = 0; jj < bp->nc[1]; jj++) {
            for (kk = 0; kk < bp->nc[2]; kk++) {
              av += *(sol + i + jj * bp->nc[0] + kk * bp->nc[0] * bp->nc[1]);
              num++;
            }
          }
          av /= (CA_FLOAT) num;
          fprintf (avfp, "%.10f,", av);

        }
        fprintf (avfp, "\n");

        fclose (avfp);
        if (strcmp ("G_", output_v->id_string) == 0) {
          /* do a supersat. profile also */
          sprintf (avprofname, "SS_%sAV%s.txt", output_v->id_string, bp->ctrl->fn_base);
          if ((avfp = fopen (avprofname, "a")) == NULL) {
            fprintf (stderr, "Error: can't open %s. output file [%s]\n", output_v->id_string, command);
            return;
          }
          fprintf (avfp, "%i,%f,", bp->step, bp->sim_time);

          /* average over each i=const slice */
          for (i = 0; i < bp->nc[0]; i++) {
            av = 0;
            num = 0;
            for (jj = 0; jj < bp->nc[1]; jj++) {
              for (kk = 0; kk < bp->nc[2]; kk++) {
                /* the array index of the current cell */
                index_ca = i + jj * bp->nc[0] + kk * bp->nc[0] * bp->nc[1];
                cell_temp = bp->sb[sbnum]->c_temp[index_ca];
                gas = *(sol + index_ca);
                /*if (bp->ctrl->diffuse_alloy_multi == FALSE) {*/
                if((bp->ctrl->diffuse_alloy_multi == FALSE)&&(bp->ctrl->diffuse_alloy_poly == FALSE)){
                  alloy = *(bp->sb[sbnum]->c_sol_alloy + index_ca);
                  fs = *(bp->sb[sbnum]->c_fs + index_ca);
                  ss = gas / (find_sat (&(bp->mprops), cell_temp + STD_TEMP, alloy, fs));
                  av += ss;
                }
                num++;
              }
            }
            av /= (CA_FLOAT) num;
            fprintf (avfp, "%.10f,", av);

          }
          fprintf (avfp, "\n");
          fclose (avfp);

        }
        /* end supersat profile bit */
      }
#endif /*PRINT_AV_PROFILE */

      break;
    case LAST_CALL:
#ifdef PRINT_AV_PROFILE
#endif /*PRINT_AV_PROFILE */
      fclose (output_v->fd_exel);
      break;
    default:
      fprintf (stderr, "ERROR: write_conc_profile: it should be impossible to reach this line.\n");
      fprintf (stderr, "Something is seriously wrong.\n");
      exit (2);
      break;
    }                           /*end of switch stat flag */
  }
}                               /*end of write execl conc. */

/********************************************************/
/* write_conc_prof_multi                                      */
/*    write a concentration profile to excel file        */
/********************************************************/
/********************************************************/
void write_conc_prof_multi (int stat_flag, BB_struct * bp, Value_struct * output_v, int j)
{
  int start, i, *conc_prof, sbnum;
  char command[MAX_STRING_LEN];
  static int msg[3] = { 0, 0, 0 };
  CA_FLOAT *fs, *fsp, *sol, *sp, partcoef;
  FILE *fp;                     /*file pointer */
  MultiS_struct *ms;

  ms = &(bp->MultiSvals);

  switch (stat_flag) {

  case FIRST_CALL:
    sprintf (command, "%s%s_multi.txt", output_v->id_string, bp->ctrl->fn_base);
    if ((output_v->fd_exel = fopen (command, "w")) == NULL) {
      fprintf (stderr, "Error: can't open %s. output file [%s]\n", output_v->id_string, command);
      return;
    } else {
      fp = output_v->fd_exel;
      fprintf (fp, "%s CONCENTRATION OUTPUT for MULTICOMP\n", output_v->id_string);
      write_excel_header (fp, bp);
      fprintf (fp, "Step,Time,");
      for (i = 0; i < bp->nc[0]; i++) {
        fprintf (fp, "%f,", i * bp->size_c[0]);
      }
      fprintf (fp, "\n");
    }
    break;
  case STEP_CALL:
    conc_prof = bp->ctrl->conc_prof;
    sbnum = conc_prof[0];
    if (bp->sb[sbnum]->open == FALSE)   /*sb not active yet */
      return;

    fp = output_v->fd_exel;
    fprintf (fp, "Profile started at step,%i,time,%.5g \n", bp->step, bp->sim_time);
    if (sbnum > bp->ntsb && !msg[0]) {
      fprintf (stderr, "ERROR: ConcProfile sb %i out of range\n", sbnum);
      msg[0] = TRUE;
      return;
    }
    if (conc_prof[1] >= bp->nc[2] && !msg[2]) {
      fprintf (stderr, "ERROR: ConcProfile slice %i out of range\n", conc_prof[1]);
      fprintf (fp, "ERROR: ConcProfile slice %i out of range\n", conc_prof[1]);
      msg[2] = TRUE;
      return;
    }
    if (conc_prof[2] > bp->nc[1] && !msg[1]) {
      fprintf (stderr, "ERROR: ConcProfile row %i out of range\n", conc_prof[2]);
      fprintf (fp, "ERROR: ConcProfile row %i out of range\n", conc_prof[2]);
      msg[1] = TRUE;
      return;
    }
    if (msg[0] || msg[1] || msg[2])
      return;
    sol = output_v->block_array[sbnum];
    if (bp->ctrl->scheil == TRUE) {
      fs = bp->sb[sbnum]->sch_fs;
    } else {
      fs = bp->sb[sbnum]->c_fs;
    }
    start = bp->nc[0] * bp->nc[1] * conc_prof[1] + bp->nc[0] * conc_prof[2];
    sp = sol + start;
    fsp = fs + start;
    fprintf (fp, "%i,%f,", bp->step, bp->sim_time);
    for (i = 0; i < bp->nc[0]; i++) {
      fprintf (fp, "%.10f,", *sp);
      sp++;
    }
    fprintf (fp, "\n");
    break;
  case LAST_CALL:
    fclose (output_v->fd_exel);
    break;
  default:
    fprintf (stderr, "ERROR: write_conc_profile: it should be impossible to reach this line.\n");
    fprintf (stderr, "Something is seriously wrong.\n");
    exit (2);
    break;
  }                             /*end of switch stat flag */
}                               /*end of write execl conc multi. */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_output_ex_c ()
{
  static char const rcsid[] = "$Id: output_ex.c 1382 2008-09-24 14:59:34Z  $";

  return (rcsid);
}

/* end of rcs_id_output_ex_c subroutine */
/*
*/
