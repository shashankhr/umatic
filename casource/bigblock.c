
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
/*      bigblock.c:                                    */
/*  All subroutine related to the bigblock (sperblock),            */
/*  including:                                          */
/*      init_bb:      initialise the bigblock structure      */
/*      read_bb:      read in ca grid and details            */
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Jul 1, 1998                                                  */
/****************************************************************/
/*       MODIFIED by:                                           */
/****************************************************************/
/****** To Do List **********************************************/
/****************************************************************/

/*RCS Id:$Id: bigblock.c 1376 2008-09-08 13:41:09Z  $*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "machine.h"
#include "fem.h"
#include "read_ctrl.h"
#include "blocks.h"
#include "nuc_lookup.h"
#include "gaussdev.h"

extern CA_FLOAT global_pressure;

#ifdef EXTERNAL_CA
extern float t_final;
#endif

void init_pr_lookup (BB_struct * bp, Ctrl_str * cp);
void init_temp_lookup (BB_struct * bp, Ctrl_str * cp);
extern void alloc_ms (MultiS_struct * ms, int ele_1, int nc);

/* functions used from read_cap_ca.c */
extern int read_cap_ca (Ctrl_str *, BB_struct *);

/* functions used from subblock.c */
extern int read_geoplus (Ctrl_str *, BB_struct *);

/* functions used from readmat.c */
extern int read_matprop (Ctrl_str * cp, BB_struct * bp);

/* functions used from readphase.c */
extern int readphase (Ctrl_str * cp, BB_struct * bp);


/* functions used from thermo_trace_init.c */
extern int thermo_trace_init (TC_str * tc);

/* functions used from fidap_interp_init.c */
extern int init_output_img (RGB_struct *);

/* functions used from initface.c */
extern int init_facecode (Frame * cubeptr, int *ins, int dim);

/* functions used from initcube.c */
extern int init_cube (BB_struct * bp);

/*function from nbhd_def.c*/
extern void nbhd_def (BB_struct * bp);

/* pass in a blank BB structure, return inited        */
extern int alloc_bb (BB_struct * bp);

/**************************************************/
/*      Subroutine to initialize bigblock        */
/**************************************************/
/* pass in a blank BB structure, return inited        */
int init_bb (BB_struct * bp, Ctrl_str * cp)
{
  char command[MAX_STRING_LEN];
  int i, j, k, l, ele_num, ele_1;       /* tmp counters */
  int nc, local_index;
  int n_nuc_beut, n_nuc_teut;
  extern void srand48 (long);
  int errors = 0;
  MultiS_struct *ms;
  Nuc_str *np;
  FILE *fptr1, *fptr2;
  CA_FLOAT *cell_index_beut, *cell_index_teut;

  ms = &(bp->MultiSvals);
  np = &(bp->nprops);

   /********************set the number of components**/
  ele_num = cp->NUM_COMP;
  ele_1 = ele_num - 1;

   /************************************************/
  /* Start off the random number generator        */
   /************************************************/
  srand48 ((long) (cp->seed));
   /************************************************/
  /* If running SOLO, read geo from commands file, */
  /* or read it from CAP_CA file                  */
   /************************************************/
  if (cp->solo || cp->external) {
    read_geoplus (cp, bp);
      /************************************************/
    /* If running CAP, read geo from commands file, */
    /* or read it from CAP_CA file                  */
      /************************************************/
  } else if (cp->cap) {
    read_cap_ca (cp, bp);
    bp->dim = THREE_D;
  } else {
    /* do nothing */
  }

  /* Open up output file                        */

  fprintf (stderr, "calc. all subsiduary values...\n");
   /************************************************/
  /* finished reading all of the header, now      */
  /* calc. all subsiduary values...            */
   /************************************************/
  bp->ntsb = bp->nsb[0] * bp->nsb[1] * bp->nsb[2];
  bp->ncsb = (bp->tnc[0] * bp->tnc[1] * bp->tnc[2]) / bp->ntsb;
  if (bp->ctrl->external == 0) {
    for (i = 0; i < 3; i++)
      bp->size_bb[i] = (CA_FLOAT) bp->tnc[i] * bp->size_c[i];
  }
  if (cp->external == 0) {
    bp->nsteps = (int) (bp->finish_time / bp->delt);
  }
  if (cp->external == 0) {
    cp->nscrdumps = (int) bp->nsteps / cp->scr_dmp_freq;
  } else {
#ifdef EXTERNAL_CA
    cp->nscrdumps = (int) (t_final / (bp->delt * cp->scr_dmp_freq));
#endif
  }
  bp->scr_dump_num = 0;
  bp->blk_dump_num = 0;
  bp->slice_dump_num = 0;

  bp->pore_err = 0;


   /************************************************/
  /* Read in all of the material properties.      */
   /************************************************/
  fprintf (stderr, "calling read_matprop().\n");
  read_matprop (cp, bp);

      /*****************/
      /*****************/
  /* allocate all  */
  /* the memory    */
  /* for the big   */
  /* block struct  */
      /*****************/
      /*****************/

  alloc_bb (bp);

/*******************************************/
/* if there is a PRESSURE LOOKUP TABLE     */
/*******************************************/
  if (cp->pr_lookup == TRUE) {
    init_pr_lookup (bp, cp);
  }

/*******************************************/
/* if there is a TEMPERATURE LOOKUP TABLE     */
/*******************************************/
  if (cp->temp_lookup == TRUE) {
    init_temp_lookup (bp, cp);
  }

   /************************************************/
  /* Read in all of the phase diagram.      */
   /************************************************/
  if (bp->ctrl->thermocalc == TRUE) {
    fprintf (stderr, "calling readphase().\n");
    readphase (cp, bp);
  }
   /************************************************/
  /* precalculate tilt cosines      */
   /************************************************/
if(cp->gradtilt == 1 ){
/* 2-d rotation only */
   fprintf(stderr,"Calculating gradient angle cosines\n");
     bp->grad_ct=cos(((PI)/180.0) * cp->grad_angle[0]);
     bp->grad_st=sin(((PI)/180.0) * cp->grad_angle[0]);
   
  }else if (cp->gradtilt == 2 ) {
  /* 3d rotation */
  fprintf(stderr,"3d gradient not implemented yet! \n");
  exit(0);
  }
   

   /*******************************************/
  /*                                         */
  /* calculate DAS from cooling rate         */
   /*******************************************/
  if (cp->coolrate == TRUE) {
    CA_FLOAT das_pre, das_exp;

    das_pre = DAS_PRE;
    das_exp = DAS_EXP;
    bp->mprops.das = bp->mprops.das_factor * das_pre * POW (((-cp->delT)), das_exp);
  }

  bp->c_sol_values->part_coef = bp->mprops.gasprops.part_coef[0]; /*  gas only  */

  /** \todo This needs to be redone for poly-component , create multiblock arrays*/
  if (cp->diffuse_alloy_poly == FALSE){
     bp->c_sol_alloy_values->part_coef = bp->mprops.alloyprops[0].part_coef[0];
  }else{
     printf("WARNING:%s: -- initialization of polycomponent values is not complete\n",__func__);
     /* ?????????????? */
  }
  bp->mprops.rhocp = bp->mprops.rho * bp->mprops.cp;
  bp->vol_bb = bp->size_bb[0] * bp->size_bb[1] * bp->size_bb[2];
  bp->area_bb = 2.0 * ((bp->size_bb[0] * bp->size_bb[1])
                       + (bp->size_bb[0] * bp->size_bb[2]) + (bp->size_bb[1] * bp->size_bb[2]));
  bp->vol_sb = bp->ncsb * bp->vol_c;
  fprintf (stderr, "vol_sb:,%.10e\n", bp->vol_sb);
  /* QSV = Q * S / V */
#ifdef HEAT_FLUX_SURF
  bp->mprops.QSV = bp->mprops.CA_Q * bp->area_bb / bp->vol_bb;
#else
  bp->mprops.QSV = bp->mprops.CA_Q;
#endif /*HEAT_FLUX_SURF */
  bp->CR = bp->mprops.QSV / bp->mprops.rhocp;

   /********************************************************/
  /* If the nuc model is DIST (distribution from lookup)  */
   /********************************************************/
  if (bp->nprops.nmodel == N_DIST) {
    errors += init_nuc_func_table (&(bp->nprops));
    errors += init_nuc_values (&(bp->nprops));
    errors += init_nuc_dist (&(bp->nprops));
    errors += init_nuc_lookup (&(bp->nprops));
    errors += output_nuc_lookup (&(bp->nprops));
#ifdef MAKE_NUC_TABLE
    {
      extern int dist_cell_nuc (BB_struct * bp, CA_FLOAT Tunder);
      CA_FLOAT testund = 0;
      int teststep;

      for (teststep = 0; teststep < 1000; teststep++) {
        dist_cell_nuc (bp, testund);
        testund += 0.01;
      }
    }

    exit (0);
#endif /*MAKE_NUC_TABLE */
  }

   /********************************************************************/

   /************************************************/
  /* Init tctrace geom and results if required...   */
   /************************************************/
  if (bp->ctrl->t_input) {
    bp->tc.fn = strdup (bp->ctrl->fn_t_input);
    fprintf (stderr, "Calling thermo_trace_init for initialisation.\n");
    thermo_trace_init (&(bp->tc));
    bp->Tinit = bp->tc.Temp[0]; /* override Tinit with first temp */
    fprintf (stderr, "Finished thermo_trace_init for initialisation.\n");
    {
         /*******************************************/
      /*                                         */
      /* calculate DAS based on tc trace file    */
         /*******************************************/
      CA_FLOAT das_pre, das_exp;
      CA_FLOAT ftemp, ltemp, coolrate, ftime, ltime;

      ftemp = bp->tc.Temp[0];
      ltemp = bp->tc.Temp[bp->tc.Nlines - 1];
      ftime = bp->tc.Time[0];
      ltime = bp->tc.Time[bp->tc.Nlines - 1];
      coolrate = (ltemp - ftemp) / (ltime - ftime);
      das_pre = DAS_PRE;
      das_exp = DAS_EXP;
      bp->mprops.das = bp->mprops.das_factor * das_pre * POW (((-coolrate)), das_exp);
    }

  }
  /* init implemented neighbourhoods */
  nbhd_def (bp);

   /****************************************************/
  /* If solo_mode and auto_finish time, calcualte it! */
   /****************************************************/
  if (cp->solo && bp->auto_fin) {
    bp->autofin_time = -1.4 * (bp->mprops.rhocp * (bp->Tinit - bp->mprops.Tliq)
                               + bp->mprops.latentH) / bp->mprops.QSV;
    fprintf (stderr, "Autofin_time set = %g\n", bp->autofin_time);
  }

   /************************************************/
  /* Init output variables ....                   */
   /************************************************/
  fprintf (stderr, "Calling init_output_img for initialisation.\n");
  bp->ctrl->rgbp->grey = bp->ctrl->rgbgrey;
  init_output_img (bp->ctrl->rgbp);
  fprintf (stderr, "Finished init_output_img for initialisation.\n");

   /*************************************/
  /* Print out checks on input data... */
   /*************************************/
  fprintf (stderr, "Exiting init_bb().\n");
  fprintf (stderr, "sb_mask nsb: %d, %d, %d\n", bp->nsb[0], bp->nsb[1], bp->nsb[2]);

  return (0);
}                               /* end of init_bb subroutine */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_bigblock_c ()
{
  static char const rcsid[] = "$Id: bigblock.c 1376 2008-09-08 13:41:09Z  $";

  return (rcsid);
}

/* end of rcs_id_bigblock_c subroutine */
/*
*/
