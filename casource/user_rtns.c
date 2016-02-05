
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

/*******************************************************************/
/* this subroutine is to be used as the new wrapper for the ca code*/
/* for a coupled ca external mode************************************/
/*******************************************************************/

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <signal.h>
#include <setjmp.h>

#include "common.h"
#ifdef MYALLOC
float * float_alloc(size_t nn){
   float * retptr;
   retptr = (float*)calloc(nn,sizeof(float));
   return(retptr);
}
int * int_alloc(size_t nn){
   int * retptr;
   retptr = (int*)calloc(nn,sizeof(int));
   return(retptr);
}
double * double_alloc(size_t nn){
   double * retptr;
   retptr = (double*)calloc(nn,sizeof(double));
   return(retptr);
}
#endif
#include "blocks.h"             /* included for def. of INIT_BB, CALC_BB, FINISH_BB */
#include "read_ctrl.h"          /* included for def. of READ_CTRL subroutine */

/* header for enabling fpe traps in gnu glibc  library */
#ifdef GNU_TRAP_FPE
#  include <fenv.h>
#endif /* GNU_TRAP_FPE */

#include "handlers.h"
#ifdef DMALLOC
#include <dmalloc.h>
#endif

#define STD_TMP 273.16
#define DYN_ATM 1.0e-06
#define D_CA_EXTERNAL_FILE "umat_extern.in"
#define D_CTRL_FILE "umat_ctrl.in"
#define D_
#define TOTAL 1
#define LOCAL 2

/*************************************************************/
/* needed to handle signals defined in CA code by Robert Atwood */
jmp_buf env;                    /*jump environment buffer for signal handler */
int jflg = JFLG_END;            /* flag for behaviour of signal handler */
int the_signo;
static int sig_flag = 0;
int signal_change_freq = 0;

/*************************************************************/

extern int init_nuc_mould (BB_struct * bp, int sbnum, FILE * listfile);
extern int init_bb (BB_struct * bp, Ctrl_str * cp);

/*function from rcs_id.c*/
extern void print_rcs_id (FILE * outfile);

/* functions used from umat_solid.c */
extern int umat_solid (int stat_flag, CA_FLOAT time, CA_FLOAT delt, Ctrl_str * cp, BB_struct * bp);

/* functions used from read_ctrl.c */
extern int read_ctrl (char *filename, Ctrl_str * cp);

/* functions used from read_umat_extern.c */
extern int read_umat_extern (Ctrl_str * cp, BB_struct * bp);

extern void src_assem (); /* external function ? */
extern double time_mod ();/* external function ? */
extern double temp_mod ();/* external function ? */
extern double quad_mod ();/* external function ? */
extern double cool_mod ();/* external function ? */
CA_FLOAT global_pressure;

void external_sb_set_cells (BB_struct * bp, int sbnum)
{
  int i, j, k, ii;
  int bb_idx;                   /* the index of the cell within the whole bigblock */
  int bb_ijk[3];                /* the bigblock index components of the cell */
  CA_FLOAT *gasp, *alloyp;
  CA_FLOAT *fsp;
  int *gp;

  gp = bp->sb[sbnum]->gr;
  gasp = bp->sb[sbnum]->c_sol;
  alloyp = bp->sb[sbnum]->c_sol_alloy;
  fsp = bp->sb[sbnum]->c_fs;

/* run through the array setting up the mould cells (and multi-component pointer if necessary */
  for (k = 0; k < bp->nc[2]; k++) {
    for (j = 0; j < bp->nc[1]; j++) {
      for (i = 0; i < bp->nc[0]; i++) {
        /*transform the subblock index to a bigblock index for finding the Procast element */
        bb_ijk[2] = bp->sb[sbnum]->orig_bb_idx[2] + k;
        bb_ijk[1] = bp->sb[sbnum]->orig_bb_idx[1] + j;
        bb_ijk[0] = bp->sb[sbnum]->orig_bb_idx[0] + i;

        bb_idx = bb_ijk[2] * bp->tnc[0] * bp->tnc[1] + bb_ijk[1] * bp->tnc[0] + bb_ijk[0];

        if (bp->cell_element_array[bb_idx] == NOT_CASTING) {
          /* count the number of mould results in the subblock */
          bp->sb[sbnum]->nmould++;
          /* set the array elements to flag them */
          *gp = NOT_CASTING;
          *fsp = NOT_CASTING;
          if (bp->ctrl->diffuse == 1) {
            *gasp = NOT_CASTING;
          }
          if (bp->ctrl->diffuse_alloy == 1) {
            *alloyp = NOT_CASTING;
          }
        }                       /* end -- mark the not-casting cells */
        gp++;
        fsp++;
        if (bp->ctrl->diffuse == 1) {
          gasp++;
        }                       /* diffuse gas increment */
        if (bp->ctrl->diffuse_alloy == 1) {
          alloyp++;
        }                       /*diffuse alloy increment */
      }                         /*end -- i loop */
    }                           /*end -- j loop */
  }                             /*end -- k loop */
/* end -- run through cells looking for mould */
}
void cell_element_output (BB_struct * bp, int slice)
{
  FILE *fd;
  char filename[MAX_STRING_LEN];
  int *external_elem;
  int nxny;
  int i;

  nxny = bp->nc[0] * bp->nc[1]; /* number of cells in a slice */
  sprintf (filename, "I_ELEM_%s.dat", bp->ctrl->fn_base);
  if ((fd = fopen (filename, "w")) == NULL) {
    fprintf (stderr, "ERROR: cell_element_output: can't open  file %s\n", filename);
    return;
  }
  external_elem = bp->cell_element_array + slice * nxny;
  fprintf (fd, "VARIABLES = \"ProCast element\"\n");
  fprintf (fd, "ZONE i=%i,j=%i \n", bp->nc[0], bp->nc[1]);
  for (i = 0; i < nxny; i++) {
    fprintf (fd, "%i\n", *external_elem % 100);
    /*fprintf(fd,"%i\n",*external_elem ); */
    external_elem++;
  }
  fclose (fd);
}

/*************function for bilinear interpolation of temperature and *****/
/*************pressure values within a three-node or four-node 2D mesh****/
/*************using the mesh nodal values*********************************/

float interpolate_p (int element, float x, float y, float interpolated)
{
#define TRIANGLE 7
#define RECTANGLE 6

  int i;
  int current_node;
  int current_npe;

  float delta, test_s_1, test_s_2;
  float *xi, *yi, *value;
  float h[4], a[3], b[3], c[3];
  float aa, bb, cc, dd, aa2, bb2, cc2, dd2;
  float r, s;
  float c1, c2, c3, c4;
  float surface;

  current_npe = npe[el_type[element]];
  xi = float_alloc (current_npe);
  yi = float_alloc (current_npe);
  value = float_alloc (current_npe);

  for (i = 0; i < current_npe; i++) {
    current_node = ncon[element][i];
    xi[i] = x_cord[current_node];
    yi[i] = y_cord[current_node];
    value[i] = t1[current_node];
  }

  interpolated = 0.0;

  switch (el_type[element]) {
  case TRIANGLE:

    a[0] = xi[1] * yi[2] - xi[2] * yi[1];
    a[1] = xi[2] * yi[0] - xi[0] * yi[2];
    a[2] = xi[0] * yi[1] - xi[1] * yi[0];

    b[0] = yi[1] - yi[2];
    b[1] = yi[2] - yi[0];
    b[2] = yi[0] - yi[1];

    c[0] = xi[2] - xi[1];
    c[1] = xi[0] - xi[2];
    c[2] = xi[1] - xi[0];

    surface = (xi[0] * yi[1]) + (xi[1] * yi[2]) + (xi[2] * yi[0]) - (yi[0] * xi[1]) - (yi[1] * xi[2]) - (yi[2] * xi[0]);
    for (i = 0; i < current_npe; i++) {
      h[i] = (1.0 / surface) * (a[i] + x * b[i] + y * c[i]);
    }
    break;
  case RECTANGLE:

    aa = 0.0;
    aa2 = 0.0;

    for (i = 0; i < current_npe; i++) {
      aa += yi[i];
      aa2 += xi[i];
    }
    bb = yi[0] - yi[1] - yi[2] + yi[3];
    bb2 = xi[0] + xi[1] - xi[2] - xi[3];

    cc = yi[0] + yi[1] - yi[2] - yi[3];
    cc2 = xi[0] - xi[1] - xi[2] + xi[3];

    dd = yi[0] - yi[1] + yi[2] - yi[3];
    dd2 = xi[0] - xi[1] + xi[2] - xi[3];

    c1 = ((4 * y - aa) * cc2) - (4 * bb * x - aa2 * bb);
    c2 = (4 * y * dd2 - aa * dd2 + bb * bb2);
    c3 = cc * cc2 + 4 * dd * x - dd * aa2;
    c4 = cc * dd2 - dd * bb2;

    delta = ((c2 - c3) * (c2 - c3)) + 4 * c1 * c4;

    if (delta >= 0.0) {
      test_s_1 = ((c2 - c3) + powf (delta, 0.5)) / (2.0 * c4);
      test_s_2 = ((c2 - c3) - powf (delta, 0.5)) / (2.0 * c4);
    } else {
      fprintf (stderr, "error in the interpolation routine, delta is negative \n");
      return (-1);
    }

    if (test_s_1 >= -1.1 && test_s_1 <= 1.1) {
      s = test_s_1;
    } else {
      if (test_s_2 >= -1.1 && test_s_2 <= 1.1) {
        s = test_s_2;
      } else {
        fprintf (stderr, "error in the interpolation routine, value of s out of boundary \n");
        return (1);
      }
    }

    r = (4 * x - (aa2 + bb2 * s)) / (cc2 + dd2 * s);

    h[0] = 0.25 * (1 + r) * (1 + s);
    h[1] = 0.25 * (1 - r) * (1 + s);
    h[2] = 0.25 * (1 - r) * (1 - s);
    h[3] = 0.25 * (1 + r) * (1 - s);

    break;
  }

  for (i = 0; i < current_npe; i++) {
    interpolated += h[i] * value[i];
  }

  free (xi);
  free (yi);
  free (value);

  return (interpolated);

}                               /*end of interpolate_p */

/***********function for sorting an array of n elements using the shell****/
/***********sort method****************************************************/

void shell (int n, float *a)
{
  int i, j, inc;
  float v;

  inc = 1;

  do {
    inc *= 3;
    inc++;
  } while (inc <= n);
  do {
    inc /= 3;
    for (i = inc + 1; i <= n; i++) {
      v = a[i];
      j = i;
      while (a[j - inc] > v) {
        a[j] = a[j - inc];
        j -= inc;
        if (j <= inc)
          break;
      }
      a[j] = v;
    }

  } while (inc > 1);
  for (j = 1; j <= n; j++) {
  }
}

/*****************function for the cellular automata cell generation*****/
/*****************super imposed on the external mesh**********************/

void cell_generation (BB_struct * bp, Ctrl_str * cp, int *cell_number_total)
{
 /***********local variables****************/
  int j, k, h;
  float *shell_cord_x, *shell_cord_y, *shell_cord_z;
  float *cell_cord_x, *cell_cord_y, *cell_cord_z;
  float sizex, sizey, sizez, maxsize;

#ifdef EXTERNAL_CELL
  FILE *test1;
#endif
  int cell_counter;

  *cell_number_total = 1;
  cp->fn_umat_extern = strdup (D_CA_EXTERNAL_FILE);
  fprintf(stderr,"Running %s\n",__func__);

/*********************************************/
#ifdef EXTERNALCELL
  test1 = fopen ("external.cell", "w+");
#endif
/**********************************************/

/**********for 2D case allocate and set all z coordinates to zero**/
  if (TWO_D) {
    z_cord = double_alloc (nnod);
    for (j = 0; j < nnod; j++) {
      z_cord[j] = 0.0;
    }
  }

  /* call the subroutine to read in the umat_extern input file */
  /* get the control options */
  read_ctrl (D_CTRL_FILE, cp);
  bp->ctrl = cp;
  /* read the big block geometry */
  if (init_bb (bp, cp) != 0) {
    fprintf (stderr, "exiting due to init_bb failure");
    exit (4);
  }
  /* just read the cell option and bigblock origin now (robert) */
  read_umat_extern (cp, bp);
  if (cp->mould_source_freq != 0) {
    cp->mould_source_freq = (2 * PI * bp->size_c[0]) / (cp->mould_source_freq);
  }

/********total option for cell generation*********/

/* this option apparently resizes the cell, to cover the whole macromodel with the specified number of cells */
/* it would bbe more useful if the time step were also adjusted! */

  if (bp->cell_option == TOTAL) {
      /****generate cells for the whole cast*******************/

      /********allocation of local arrays and variables**********/
    shell_cord_x = float_alloc (nnod + 1);
    shell_cord_y = float_alloc (nnod + 1);
    shell_cord_z = float_alloc (nnod + 1);
      /*****************************************************************/

      /***************transform the length unit from cm to m and sort all coordinates****/
    for (j = 0; j < nnod; j++) {
      shell_cord_x[j + 1] = x_cord[j] * 0.01;
      shell_cord_y[j + 1] = y_cord[j] * 0.01;
      shell_cord_z[j + 1] = z_cord[j] * 0.01;
    }
    shell (nnod, shell_cord_x);
    shell (nnod, shell_cord_y);
    shell (nnod, shell_cord_z);
      /*********************************************************************************/

#ifdef EXTERNALCELL
    fprintf (test1, "min for x: \t %f \t max for x: \t %f \n", shell_cord_x[1], shell_cord_x[nnod]);
    fprintf (test1, "min for y: \t %f \t max for y: \t %f \n", shell_cord_y[1], shell_cord_y[nnod]);
    fprintf (test1, "min for z: \t %f \t max for z: \t %f \n", shell_cord_z[1], shell_cord_z[nnod]);
#endif

      /******************change the cell size***********/
    sizex = (fabs (shell_cord_x[nnod] - shell_cord_x[1])) / bp->tnc[0];
    sizey = (fabs (shell_cord_y[nnod] - shell_cord_y[1])) / bp->tnc[1];
    sizez = (fabs (shell_cord_z[nnod] - shell_cord_z[1])) / bp->tnc[2];
    /* but the cell size needs to be a cube */
    maxsize = MAX (MAX (sizex, sizey), sizez);

    fprintf (stderr, "CELL_GENERATION: cell_option 2 (TOTAL): changing cell size from %.5em to %.5em\n", bp->size_c[0], maxsize);
    if (maxsize <= 0) {
      fprintf (stderr, "ERROR: cell_generation: incorrect size generated %.10g\nExiting...", maxsize);
      exit (189);
    }

    for (j = 0; j < 2; j++) {
      bp->size_c[0] = maxsize;
    }

#ifdef JUNK
      /******************find number of cells in each direction***********/
    bp->tnc[0] = (int) (fabs (shell_cord_x[nnod] - shell_cord_x[1]) / bp->size_c[0]);
    bp->tnc[1] = (int) (fabs (shell_cord_y[nnod] - shell_cord_y[1]) / bp->size_c[1]);
    if (TWO_D) {
      bp->tnc[2] = 1;
    } else {
      bp->tnc[2] = (int) (fabs (shell_cord_z[nnod] - shell_cord_z[1]) / bp->size_c[2]);
    }
#endif /*JUNK*/
      /****************************/
    bp->orig_bb[0] = shell_cord_x[1];
    bp->orig_bb[1] = shell_cord_y[1];
    bp->orig_bb[2] = shell_cord_z[1];
      /****************************/
  } /*********end of total option for cell generation*******/
  else if (bp->cell_option > 2) {
    /*neither LOCAL nor TOTAL was selected, error condition */
    printf ("ERROR:%s: you entered a wrong option for cell_option, please try again. \n",__func__);
    exit (1);
  }

         /************end of user entry for cell generation*************/

/* If the TOTAL option is not selected, the given values are used (LOCAL option) */
 /********************************************************************************/
  for (j = 0; j < 3; j++) {
    *cell_number_total = (*cell_number_total) * (bp->tnc[j]);
#ifdef EXTERNALCELL
    fprintf (test1, "cell number for %d direction: %d \n", j, bp->tnc[j]);
#endif
  }
#ifdef EXTERNALCELL
  fprintf (test1, "total cell number: %d \n", *cell_number_total);
#endif
   /************************************************************************/

   /**************allocate cell coordinate arrays for ca code***************/
  bp->cell_cord_x = cell_cord_x = float_alloc ((*cell_number_total) + 1);
  bp->cell_cord_y = cell_cord_y = float_alloc ((*cell_number_total) + 1);
  bp->cell_cord_z = cell_cord_z = float_alloc ((*cell_number_total) + 1);
   /************************************************************************/

#ifdef EXTERNALCELL
   /****************output node coordinates for output control************/
   /* the arrays x_cord, y_cord, z_cord are provided by Procast */
  for (j = 0; j < nnod; j++) {
    fprintf (test1, "node_num: %d \t node_cord_x: %f \t node_cord_y: %f \t node_cord_z: %f \n", j, x_cord[j], y_cord[j], z_cord[j]);
  }
#endif
   /**************************generate cell coordinates for ca*****************/
  for (j = 0, cell_counter = 0; j < bp->tnc[2]; j++) {
    for (k = 0; k < bp->tnc[1]; k++) {
      for (h = 0; h < bp->tnc[0]; h++) {

        bp->cell_cord_x[cell_counter] = bp->orig_bb[0] + h * bp->size_c[0];
        bp->cell_cord_y[cell_counter] = bp->orig_bb[1] + k * bp->size_c[1];
        bp->cell_cord_z[cell_counter] = bp->orig_bb[2] + j * bp->size_c[2];

#ifdef EXTERNALCELL
        fprintf (test1, "cell number: %d \t x_cord: %f \t y_cord: %f \t z_cord: %f \n", cell_counter, bp->cell_cord_x[cell_counter],
                 bp->cell_cord_y[cell_counter], bp->cell_cord_z[cell_counter]);
#endif

        cell_counter++;
      }
    }
  }                             /* end of three-d loop through cells */
   /****************************************************************************/
#ifdef EXTERNALCELL
  fprintf (test1, "end of cell generation and before the allocation \n");
#endif

   /************assigne bigblock variables for block and cells size and number****/
  for (j = 0; j < 3; j++) {
    bp->size_bb[j] = (bp->tnc[j]) * bp->size_c[j];
    bp->nc[j] = (bp->tnc[j]) / (bp->nsb[j]);
  }
   /*****************************************************************************/

#ifdef EXTERNALCELL
  fclose (test1);
#endif

  return;
}                               /*end of cell_generation */

/********************the main routine for the ca external link calling the umat_solid from ca*****/
int umat_extern_wrapper (int step_num, BB_struct * bp, Ctrl_str * cp)
{

  extern int *cell_elem;

   /*************defined from user_rtn_cell.c*******************/
  extern void user_rtn_cell (int number_total,
                             int *nc,
                             float *origine,
                             float *cell_size_temp,
                             float *umat_t0_1,
                             float *umat_p0_1, float *umat_t1_1, float *umat_p1_1, float *xxx, float *yyy, float *zzz, int number);

  CA_FLOAT *gasp, *alloyp;
  CA_FLOAT *fsp;
  CA_FLOAT adj_delt = 0.0, old_delt = 0.0;
  int adj_flag = 0;
  int *gp;
  int error_umat_solid;
  char *ctrl_fname;
  int i, j;
  int ele_num, ele_1;
  int cell_number_total;
  int stat_flag;
  int startstep;

#ifdef VERBOSE_EXTERNAL
  FILE *test4, *test9;
#endif /* VERBOSE_EXTERNAL */
  int current_element, current_npe, current_node, node_number;
  float node_distance, nearest_node, nearest_x, nearest_y, nearest_z;
  CA_FLOAT del_ptime = 1.0;
  CA_FLOAT time = 0.0;
  CA_FLOAT delt = 1.0;
  FILE *rcs_id_file;
  FILE *listfile;
  char filename[MAX_STRING_LEN];

  ctrl_fname = strdup (D_CTRL_FILE);

/************************************************/

#ifdef VERBOSE_EXTERNAL
  test4 = fopen ("external.inter", "w+");
  test9 = fopen ("external.comp", "a+");
#endif /* VERBOSE_EXTERNAL */

/**************************************************************/

/**********for first macro time step call the cell generation routine******/
 #ifndef START_STEP
  startstep = inilev + 1;
 #else
  startstep = inilev + START_STEP;
  if (step_num < startstep){
     fprintf (stderr,"Step num %i is less than startstep %i\n",step_num,startstep);
     return(1);
  }
 #endif

  printf("step num %i, startstep %i, sim_time %f, t_final %f\n",step_num,startstep,bp->sim_time,t_final);

  if ((step_num == startstep) && bp->sim_time < t_final) {
  printf("step num %i, startstep %i, sim_time %f, t_final %f\n",step_num,startstep,bp->sim_time,t_final);

/*********if first time called print rcsid******/
#ifdef PRINT_RCS_ID
    print_rcs_id (stderr);
#endif /*PRINT_RCS_ID */
    rcs_id_file = fopen ("rcs_id.txt", "w");
    if (rcs_id_file == NULL) {
      fprintf (stderr, "ERROR: could not open file rcs_id.txt\n");
      exit (0);
    }
    print_rcs_id (rcs_id_file);
    fclose (rcs_id_file);

    /* read in the control files and generate the cells */
    /* this happens inside cell_generation */
    cell_generation (bp, cp, &cell_number_total);

    /* set the finish time to the time given by Procast */
    /**  \todo  straighten out so this only need to be done in one place  -- coupled*/
    fprintf (stderr, "step number is : %d \n", step_num);
    bp->finish_time = (double) (current_time);
                                          /**set the finish time for ca**/
    fprintf (stderr, "bp->finish_time %.10g, Procast current_time %.10g, Procast timestep %.10g\n", bp->finish_time, current_time, dt);
    bp->sim_time = (CA_FLOAT) (current_time - dt);
                                            /**initialise the simulation time for ca**/
#ifdef VERBOSE_EXTERNAL
    fprintf (test4, "the current finish time for CA is : %lf \n", bp->finish_time);
    fprintf (test4, "bp->nsteps=%d \n", bp->nsteps);
    fprintf (test4, "Calling umat_solid for one macro timestep for the first time.\n");
#endif /* VERBOSE_EXTERNAL */
   /**************************************************/

   /*********calcualte the number of microscale itterations*******/
    /* adjust time step if there are too few */
    bp->nsteps = (int) (FLOOR ((bp->finish_time - bp->sim_time) / bp->delt));
    if (bp->nsteps < EXTERNAL_MIN_STEPS) {
      fprintf (stderr, "%s: Small macro time step. %i < %i\n",__func__, bp->nsteps, EXTERNAL_MIN_STEPS);
      fprintf (stderr, "Adjusting micro time step ...\n");
      adj_delt = ((CA_FLOAT) dt) / ((CA_FLOAT) (EXTERNAL_MIN_STEPS));
      old_delt = bp->delt;
      bp->delt = adj_delt;
      if (adj_delt > old_delt) {
        fprintf (stderr, "ERROR: umat_extern_wrapper: bad adjustment %.5g %.5g\n", adj_delt, old_delt);
        exit (8);
      }
      adj_flag = 1;
      bp->nsteps = (int) (FLOOR ((bp->finish_time - bp->sim_time) / bp->delt));
      fprintf (stderr, "Original timestep old_delt is %.5g\n", old_delt);
      fprintf (stderr, "Adjusted timestep adj_delt is %.5g\n", bp->delt);
    }
#ifdef VERBOSE_EXTERNAL
    fprintf (test4, "number of micro time steps is : %d \n", bp->nsteps);
#endif /* VERBOSE_EXTERNAL */

    /* check some options */
    /* if fgrid or con_cast are set then the temperature function selection will be wrong */

    bp->ctrl->external = 1;
    if (bp->ctrl->fgrid_input) {
      fprintf (stderr, "WARNING: umat_extern_wrapper: fgrid_input should not be combined with external!\n");
      fprintf (stderr, "turning OFF the option!\n");
      bp->ctrl->fgrid_input = 0;
    }
    if (bp->ctrl->con_cast) {
      fprintf (stderr, "WARNING: umat_extern_wrapper: con_cast should not be combined with external!\n");
      fprintf (stderr, "turning OFF the option!\n");
      bp->ctrl->con_cast = 0;
    }

   /**********array location for cell's temperature and pressure and their gradients****/
    bp->current_cell_temp = float_alloc (cell_number_total);
    bp->current_cell_pres = float_alloc (cell_number_total);
    bp->cell_temp_extern = float_alloc (cell_number_total);
    bp->cell_temp_change_extern = float_alloc (cell_number_total);
    bp->cell_pres_extern = float_alloc (cell_number_total);
    bp->cell_pres_grad_extern = float_alloc (cell_number_total);
   /************************************************************************************/

    bp->cell_element_array = int_alloc (cell_number_total);
                                                        /****cell mesh correspondancy array***/
    bp->cell_node_array = int_alloc (cell_number_total);
                                                     /****cell nearest node correspondancy array****/

    bp->total_cell_number = cell_number_total;

  /****calling the interpolation routine using the cell mesh correspondancy****/

    /*
       extern void user_rtn_cell(int number_total,
       int *nc,
       float *origine,
       float *cell_size_temp,
       float *umat_t0_1, 
       float *umat_p0_1, 
       float *umat_t1_1, 
       float *umat_p1_1,
       float *xxx, 
       float *yyy, 
       float *zzz,
       int number);
     */

/* find the element/cell correspondance */
/* (put into global cell_elem array defined by ProCast) */

    user_rtn_cell (cell_number_total,
                   bp->tnc,
                   bp->orig_bb,
                   bp->size_c,
                   bp->cell_temp_extern,
                   bp->cell_pres_extern,
                   bp->cell_temp_change_extern,
                   bp->cell_pres_grad_extern, bp->cell_cord_x, bp->cell_cord_y, bp->cell_cord_z, step_num);

/* transfer the cell_elem data into the bigblock structure */

    for (j = 0; j < cell_number_total; j++) {
      bp->cell_element_array[j] = cell_elem[j];
      bp->cell_temp_change_extern[j] = ((bp->cell_temp_change_extern[j] - bp->cell_temp_extern[j]) / dt);
      bp->cell_pres_grad_extern[j] = ((bp->cell_pres_grad_extern[j] - bp->cell_pres_extern[j]) / dt);
#ifdef VERBOSE_EXTERNAL
      fprintf (test4, "%f \t %f \t %f \t %f \t %d \n",
               bp->cell_temp_extern[j], bp->cell_pres_extern[j], bp->cell_temp_change_extern[j], bp->cell_pres_grad_extern[j], j);
#endif /* VERBOSE_EXTERNAL */
    }
    /* write the slice of the element number information */
    for (j = 0; j < bp->ctrl->nsbslice; j++) {
      cell_element_output (bp, bp->ctrl->slice[j][1]);
    }

 /****************put in the routine for nearest neighbour node array******/
    for (j = 0; j < cell_number_total; j++) {
      if (bp->cell_element_array[j] >= 0) {
        current_element = bp->cell_element_array[j];
        current_npe = npe[el_type[current_element]];
        for (i = 0; i < current_npe; i++) {
          current_node = ncon[current_element][i];
          nearest_x = (bp->cell_cord_x[j] - (x_cord[current_node] * 0.01));
          nearest_y = (bp->cell_cord_y[j] - (y_cord[current_node] * 0.01));
          if (TWO_D) {
            nearest_z = 0.0;
          } else {
            nearest_z = (bp->cell_cord_z[j] - (z_cord[current_node] * 0.01));
          }
          nearest_node = sqrt (nearest_x * nearest_x + nearest_y * nearest_y + nearest_z * nearest_z);
          if (i == 0) {
            node_distance = nearest_node;
            node_number = current_node;
          } else {
            if (nearest_node < node_distance) {
              node_distance = nearest_node;
              node_number = current_node;
            }
          }
        }
       /*** end of i loop***/
        bp->cell_node_array[j] = node_number;
      } /****end of if****/
      else {
        bp->cell_node_array[j] = -1;
      }
     /***end of else*****/
    }
    /*****end of j loop****/
 /******************end of the nearest node finder routine*************/

#ifdef VERBOSE_EXTERNAL
    fprintf (test4, "Calling umat_solid for initialisation \n");
#endif /* VERBOSE_EXTERNAL */
  /*********at each macro time step reset all the temperatures and ***/
  /***********pressure arrays to zero before reassignment********/
    for (j = 0; j < bp->total_cell_number; j++) {
      bp->current_cell_temp[j] = 0.0;
      bp->current_cell_pres[j] = 0.0;
    }
  /**************************************************************/
    stat_flag = EXTERNAL_INIT_BB;
    umat_solid (stat_flag, time, delt, cp, bp);
                                           /*****call the umat_solid for initialisation***/
#ifdef VERBOSE_EXTERNAL
    fprintf (test4, "Finished umat_solid for initialisation.\n");
#endif /* VERBOSE_EXTERNAL */

  /*******************set the non casting elements***********/

    fsp = bp->sb[0]->c_fs;
    gp = bp->sb[0]->gr;
    bp->sb[0]->nmould = 0;
    if (bp->ctrl->diffuse == 1) {
      gasp = bp->sb[0]->c_sol;
    }
    if (bp->ctrl->diffuse_alloy == 1) {
      alloyp = bp->sb[0]->c_sol_alloy;
    }
#ifdef JUNK
/* run through the array setting up the mould cells (and multi-component pointer if necessary */
    for (i = 0; i < cell_number_total; i++) {
      if (bp->cell_element_array[i] == NOT_CASTING) {
        /* count the number of moudl reslts in th subblock */
        bp->sb[0]->nmould++;
        /* set the array elements to flag them */
        *gp = NOT_CASTING;
        *fsp = NOT_CASTING;
        if (bp->ctrl->diffuse == 1) {
          *gasp = NOT_CASTING;
        }
        if (bp->ctrl->diffuse_alloy == 1) {
          *alloyp = NOT_CASTING;
        }
        if (bp->ctrl->diffuse_alloy_multi == 1) {
          ele_num = cp->NUM_COMP;
          ele_1 = ele_num - 1;
          for (j = 0; j < ele_1; j++) {
            bp->sb[0]->c_sol_alloy_multi[j][i] = NOT_CASTING;
          }
        }
      }                         /* end -- mark the not-casting cells */
      gp++;
      fsp++;
      if (bp->ctrl->diffuse == 1) {
        gasp++;
      }                         /* diffuse gas increment */
      if (bp->ctrl->diffuse_alloy == 1) {
        alloyp++;
      }                         /*diffuse alloy increment */
    }                           /*end -- i loop */
/* end -- run through cells looking for mould */

    /* nucleate at mould bdy */
    /**  \todo  move transferring the NOT_CASTING information to the subblock level  - coupled */
#ifdef LIST_ALL_NUC
    /* note: sbnum is zero because subblocks are not implemented for external version */
    sprintf (filename, "G_L_sb%i_%s.csv", 0, bp->ctrl->fn_base);
    listfile = fopen (filename, "w+");
    fprintf (listfile, "Surface Mould Nuclei \n");
#endif /* LIST_ALL_NUC */
    if (bp->ctrl->mould_nuc) {
      init_nuc_mould (bp, 0, listfile);
    }
#ifdef LIST_ALL_NUC
    fclose (listfile);
#endif /* LIST_ALL_NUC */
#endif /*JUNK*/
      /* call umat_solid for first macro timestep */
      stat_flag = CALC_BB;
  /*********at each macro time step reset all the temperatures and ***/
  /***********pressure arrays to zero before reassignment********/
    for (j = 0; j < bp->total_cell_number; j++) {
      bp->current_cell_temp[j] = 0.0;
      bp->current_cell_pres[j] = 0.0;
    }
  /**************************************************************/
    error_umat_solid = umat_solid (stat_flag, time, delt, cp, bp);
    if (adj_flag) {
      bp->delt = old_delt;
      adj_delt = old_delt = 0.0;
      fprintf (stderr, "Restoring time step %.5g\n", bp->delt);
    }
#ifdef VERBOSE_EXTERNAL
    if (error_umat_solid == 0) {
      fprintf (test4, "Finished first macro-timestep in umat_solid.\n");
    } else {
      fprintf (test4, "no active sub blocks in CA returning to the wrapper \n");
    }
#endif /* VERBOSE_EXTERNAL */
  }

  /* end first timestep only routines */
  /* set the finish time to the time given by Procast */
  /**  \todo  straighten out so this only need to be done in one place  -- coupled*/
  fprintf (stderr, "step number is : %d \n", step_num);
  bp->finish_time = (double) (current_time);
                                          /**set the finish time for ca**/
  fprintf (stderr, "bp->finish_time %.10g, Procast current_time %.10g, Procast timestep %.10g\n", bp->finish_time, current_time, dt);
  bp->sim_time = (CA_FLOAT) (current_time - dt);
                                            /**initialise the simulation time for ca**/
#ifdef VERBOSE_EXTERNAL
  fprintf (test4, "the current finish time for CA is : %lf \n", bp->finish_time);
#endif /* VERBOSE_EXTERNAL */

/**************************************end of first macro time step in external***********/

  fprintf(stderr,"\nProcast parameters:\n");
  fprintf(stderr,"%s: step_num %i  inilev %i  nstep %i  bp->sim_time %g t_final %g\n",
                  __func__,step_num,inilev, nstep,  bp->sim_time, t_final ); 

  fprintf(stderr,"%s:dt %g dt_old %g dtmax %g t_final %g\n",
                  __func__,dt,dt_old,dtmax,t_final);
/***  normal time step, ongoing calculation *****/
  if (step_num > inilev + 1 && step_num < inilev + nstep && bp->sim_time < t_final) {

  /**********call the interpolation function*****/
    user_rtn_cell (cell_number_total,
                   bp->tnc,
                   bp->orig_bb,
                   bp->size_c,
                   bp->cell_temp_extern,
                   bp->cell_pres_extern,
                   bp->cell_temp_change_extern,
                   bp->cell_pres_grad_extern, bp->cell_cord_x, bp->cell_cord_y, bp->cell_cord_z, step_num);

    /* calculate the change per unit time from teh two end points */
    /* may be better doen inside user_rtn_cell ??? */
    for (j = 0; j < (bp->tnc[0] * bp->tnc[1] * bp->tnc[2]); j++) {
      bp->cell_temp_change_extern[j] = ((bp->cell_temp_change_extern[j] - bp->cell_temp_extern[j]) / dt);
      bp->cell_pres_grad_extern[j] = ((bp->cell_pres_grad_extern[j] - bp->cell_pres_extern[j]) / dt);
#ifdef VERBOSE_EXTERNAL
      fprintf (test4, "%f \t %f \t %f \t %f \t %d \n",
               bp->cell_temp_extern[j], bp->cell_pres_extern[j], bp->cell_temp_change_extern[j], bp->cell_pres_grad_extern[j], j);
#endif /* VERBOSE_EXTERNAL */
    }

 /**************put temp output for one cell at each macro step*********/
#ifdef VERBOSE_EXTERNAL
    for (j = 0; j < bp->total_cell_number; j++) {
      if (bp->cell_element_array[j] >= 0) {
        fprintf (test9, "%f \t %d \t %d \t %d \t %f \t %f \t %f \t %f \n",
                 bp->sim_time,
                 j,
                 bp->cell_node_array[j],
                 bp->cell_element_array[j],
                 bp->cell_temp_extern[j],
                 bp->cell_temp_change_extern[j],
                 t0[bp->cell_node_array[j]] - STD_TMP, t0[ncon[bp->cell_element_array[j]][0]] - STD_TMP);
        break;
      }
    }
#endif /* VERBOSE_EXTERNAL */
  /*********************************************************************/

    /* call umat_solid for macro timestep */
#ifdef VERBOSE_EXTERNAL
    fprintf (test4, "Calling umat_solid for one macro timestep.\n");
#endif /* VERBOSE_EXTERNAL */
    stat_flag = CALC_BB;
  /*********at each macro time step reset all the temperatures and ***/
  /***********pressure arrays to zero before reassignment********/
    for (j = 0; j < bp->total_cell_number; j++) {
      bp->current_cell_temp[j] = 0.0;
      bp->current_cell_pres[j] = 0.0;
    }
   /***********calculate number of microscale iterations******/
    bp->nsteps = (int) (FLOOR ((bp->finish_time - bp->sim_time) / bp->delt));
#ifdef VERBOSE_EXTERNAL
    fprintf (test4, "number of micro time steps is : %d \n", bp->nsteps);
#endif /* VERBOSE_EXTERNAL */
    if (bp->nsteps < EXTERNAL_MIN_STEPS) {
      fprintf (stderr, "%s: Small macro time step %i < %i\n",__func__, bp->nsteps, EXTERNAL_MIN_STEPS);
      fprintf (stderr, "Adjusting micro time step ...\n");
      adj_delt = ((CA_FLOAT) dt) / ((CA_FLOAT) (EXTERNAL_MIN_STEPS));
      old_delt = bp->delt;
      bp->delt = adj_delt;
      if (adj_delt > old_delt) {
        fprintf (stderr, "ERROR: umat_extern_wrapper: bad adjustment %.5g %.5g\n", adj_delt, old_delt);
        exit (8);
      }
      adj_flag = 1;
      bp->nsteps = (int) (FLOOR ((bp->finish_time - bp->sim_time) / bp->delt));
      fprintf (stderr, "Original timestep old_delt is %.5g\n", old_delt);
      fprintf (stderr, "Adjusted delt adj_delt is %.5g\n", bp->delt);
    }
   /********************************************/
   /**                                         */
   /**                                         */
   /**      call CA_SOLID                      */
   /**                                         */
   /**                                         */
   /********************************************/
/* handle signals defined in handlers.h and handlers.c */
/* if a signal defined to finish occurs, go to the end of this section */
    if (setjmp (env) != 0) {
      fprintf (stderr, "CA_EXTERNAL_WRAPPER: Handling singal %i. Finishing.", the_signo);
      sig_flag = 1;
    } else {
      error_umat_solid = umat_solid (stat_flag, time, delt, cp, bp);
    }

#ifdef VERBOSE_EXTERNAL
    if (error_umat_solid == 0) {
      fprintf (test4, "Finished one macro-timestep in umat_solid.\n");
    } else {
      fprintf (test4, "no active sub blocks in CA returning to the wrapper \n");
    }
#endif /* VERBOSE_EXTERNAL */
    if (bp->last_micro == 1 || sig_flag == 1) {
      /* micromodel has finished, tell EXTERNAL to stop soon */
      nstep = step_num + 2;
      bp->last_micro = 2;
      if (sig_flag == 1)
        return (1);

    }
    if (adj_flag) {
      bp->delt = old_delt;
      adj_delt = old_delt = 0.0;
      fprintf (stderr, "Restoring time step %.5g\n", bp->delt);
    }

  }
/***********************************************************************************/

/*******************************/
/*****                      ****/
/***** last macro time step ****/
/*****                      ****/
/*******************************/
  if (step_num == inilev + nstep || bp->finish_time >= t_final) {
     if (step_num == inilev + nstep ) fprintf(stderr,"%s: finishing due to (step_num == inilev + nstep )\n",__func__); 
     if (bp->finish_time >= t_final) fprintf(stderr,"%s: finishing due to (bp->finish_time >= t_final) \n",__func__); 

    user_rtn_cell (cell_number_total,
                   bp->tnc,
                   bp->orig_bb,
                   bp->size_c,
                   bp->cell_temp_extern,
                   bp->cell_pres_extern,
                   bp->cell_temp_change_extern,
                   bp->cell_pres_grad_extern, bp->cell_cord_x, bp->cell_cord_y, bp->cell_cord_z, step_num);

    for (j = 0; j < (bp->tnc[0] * bp->tnc[1] * bp->tnc[2]); j++) {

      bp->cell_temp_change_extern[j] = ((bp->cell_temp_change_extern[j] - bp->cell_temp_extern[j]) / dt);
      bp->cell_pres_grad_extern[j] = ((bp->cell_pres_grad_extern[j] - bp->cell_pres_extern[j]) / dt);
#ifdef VERBOSE_EXTERNAL
      fprintf (test4, "%lf \t %lf \t %lf \t %lf \t %d \n",
               bp->cell_temp_extern[j], bp->cell_pres_extern[j], bp->cell_temp_change_extern[j], bp->cell_pres_grad_extern[j], j);
#endif /* VERBOSE_EXTERNAL */
    }

    /* call umat_solid for macro timestep */
#ifdef VERBOSE_EXTERNAL
    fprintf (test4, "Calling umat_solid for one macro timestep.\n");
#endif /* VERBOSE_EXTERNAL */
    stat_flag = CALC_BB;
  /*********at each macro time step reset all the temperatures and ***/
  /***********pressure arrays to zero before reassignment********/
    for (j = 0; j < bp->total_cell_number; j++) {
      bp->current_cell_temp[j] = 0.0;
      bp->current_cell_pres[j] = 0.0;
    }
   /***************calculate the number of microscales itteration******/
    bp->nsteps = ((bp->finish_time - bp->sim_time) / bp->delt);
#ifdef VERBOSE_EXTERNAL
    fprintf (test4, "number of micro time steps is : %d \n", bp->nsteps);
#endif /* VERBOSE_EXTERNAL */
    error_umat_solid = umat_solid (stat_flag, time, delt, cp, bp);
#ifdef VERBOSE_EXTERNAL
    if (error_umat_solid == 0) {
      fprintf (test4, "Finished last macro-timestep in umat_solid.\n");
    } else {
      fprintf (test4, "no active sub blocks in CA returning to the wrapper \n");
    }
#endif /* VERBOSE_EXTERNAL */

    /* tell ca prgramme it has finished and to print out */
    stat_flag = FINISH_BB;
    error_umat_solid = umat_solid (stat_flag, time, delt, cp, bp);
#ifdef VERBOSE_EXTERNAL
    if (error_umat_solid == 0) {
      fprintf (test4, "Finished final call to umat_solid.\n");
    } else {
      fprintf (test4, "no active sub blocks in CA returning to the wrapper \n");
    }

    fprintf (test4, "\nfinished all calculations successfully.\n");
    fprintf (test4, "be seeing you...\n");
#endif /* VERBOSE_EXTERNAL */

  }
/**************end of last macro step***************/

#ifdef VERBOSE_EXTERNAL
  fclose (test9);
  fclose (test4);
#endif /* VERBOSE_EXTERNAL */
  return (1);

}

/************************************main function calling the umat_extern interface at each macro step********************/

void umat_extern_step ()
{

  static BB_struct bb;
  BB_struct *bp;
  static Ctrl_str *cp;

#ifdef VERBOSE_EXTERNAL
  FILE *test10;
#endif /* VERBOSE_EXTERNAL */

  static int user_option;
  static double xx_in_max, yy_in_max, zz_in_max;
  static double xx_in_min, yy_in_min, zz_in_min;

  FILE *tout;
  FILE *pout;
  FILE *Tempout;
  FILE *eut_pout;

  int errcode = 0;
  static int allocflag = 0;

  int j;
  int user_counter;
  float local_temp;
  FILE *my_temp_out, *my_node_out;

   /***********from umat_feedback.c***********/

  extern int umat_feedback (BB_struct *, int);

/* if micro was signalled to stop, return control immediately to Macro code */
  if (sig_flag == 1)
    return;

  /* Set up the signal handler to use SIGTERM */
  signal (SIGTERM, catchit);
  /* set up a signal handler to use SIGUSR1 (writeblock and exit) */
  signal (SIGUSR1, catchit);
  /* set up a signal handler to use SIGUSR2 (writeblock and continue) */
  signal (SIGUSR2, catchit);

#ifdef GNU_TRAP_FPE
  enable_fpe_traps ();
  signal (SIGFPE, float_error);
#endif

  if (allocflag == 0) {
    cp = (Ctrl_str *) calloc (1, sizeof (Ctrl_str));
    allocflag = 1;
  }
  bp = &bb;
  bp->ctrl = cp;
  if (t_final == 0) {
    t_final = 1000000;
    fprintf (stderr, "CA_EXTERNAL:%s:WARNING: t_final is zero, resetting to LARGE %.5g\n",__func__, t_final);
  }
#ifdef USER_1
  user_option = 1;
#endif

#ifdef USER_2
  user_option = 2;
#endif

#ifdef USER_3
  user_option = 3;
#endif

  user_option = 2;
  /* small routine to print out temperatures for a given region of casting */
  if (user_option == 1) {
    my_temp_out = fopen ("temperature.out", "a+");
    if (istep == inilev + 1) {
      my_node_out = fopen ("node_numbers.out", "w");
      printf ("give the range for x coordinates \n in CentiMeters Please!! \n");
      scanf ("%lf %lf", &xx_in_min, &xx_in_max);
      printf ("%lf \t %lf \n", xx_in_min, xx_in_max);
      printf ("give the range for y coordinates \n in CentiMeters Please!! \n");
      scanf ("%lf %lf", &yy_in_min, &yy_in_max);
      printf ("%lf \t %lf \n", yy_in_min, yy_in_max);
      if (TWO_D) {
        z_cord = double_alloc (nnod);
        for (j = 0; j < nnod; j++) {
          z_cord[j] = 0.0;
        }
        zz_in_min = zz_in_max = 0.0;
      } else {
        printf ("give the range for z coordinates \n in CentiMeters Please!! \n");
        scanf ("%lf %lf", &zz_in_min, &zz_in_max);
      }
      printf ("%lf \t %lf \n", zz_in_min, zz_in_max);
    }
    /* end of istep control */
    for (j = 0, user_counter = 0; j < nnod; j++) {
      if ((x_cord[j] >= xx_in_min && x_cord[j] <= xx_in_max) && (y_cord[j] >= yy_in_min && y_cord[j] <= yy_in_max)
          && (z_cord[j] >= zz_in_min && z_cord[j] <= zz_in_max)) {
        local_temp = t0[j];
        if (user_counter == 0) {
          fprintf (my_temp_out, " %f \t %d \t %f \t", current_time, j, local_temp);
          user_counter++;
        }
        if (user_counter > 0) {
          fprintf (my_temp_out, " %d \t %f \t", j, local_temp);
        }

        if (istep == inilev + 1) {
          fprintf (my_node_out, "%d \t %lf \t %lf \t %lf \n", j, x_cord[j], y_cord[j], z_cord[j]
            );
        }
      }
    }
    fprintf (my_temp_out, " \n ");

    fclose (my_temp_out);
    if (istep == inilev + 1) {
      fclose (my_node_out);
    }

  }
  /* end of option one */
  /****************************************/
  /*USER OPTION 2 -- call the ca routines */
  /****************************************/

  if (user_option == 2) {

#ifdef VERBOSE_EXTERNAL
    test10 = fopen ("node_cord.out", "w+");
      /****************output node coordinates for output control************/
    for (j = 0; j < nnod; j++) {
      fprintf (test10, "node_num: %d \t node_cord_x: %lf \t node_cord_y: %lf \n", j, x_cord[j], y_cord[j]);
    }
#endif /* VERBOSE_EXTERNAL */
      /************************************************************************/
    fprintf (stderr, "istep= %d \n", istep);

      /***********************************************************/
    /*                                                         */
    /*       Call the CA routine                               */
    /*                                                         */
      /***********************************************************/
    fprintf(stderr,"%s: Calling umat_extern_wrapper\n",__func__);
    errcode = umat_extern_wrapper (istep, bp, cp);
    if (errcode != 1) {
      fprintf (stderr, "The umat_extern_wrapper didn't execute");
      fprintf (stderr, "Return code: %i\n", errcode);
      exit (0);
    } else {
      if (istep < inilev + nstep && total_time < t_final) {
        fprintf (stderr, "the link was successful \n \n");
      }
      if (istep >= inilev + nstep || total_time >= t_final) {
        fprintf (stderr, "the last link was successful \n \n");
      }
    }
      /**************feedback*************************/

    if (bp->last_micro != 2) {
      if (cp->umat_feedback == TRUE) {
        /* call umat_feedback */
        if (umat_feedback (bp, 0) != 1) {
          fprintf (stderr, "the umat_feedback didn't execute");
          exit (0);
        } else {
          if (istep < inilev + nstep) {
            fprintf (stderr, "the umat_feedback was successful \n \n");
          } else {
            fprintf (stderr, "the last umat_feedback was successful \n \n");
          }
        }
      }                 /**********end of feedback test*********/
    }

    /* end last_micro test */
#ifdef VERBOSE_EXTERNAL
    fclose (test10);
#endif

  }
  /* end of option two */
  if (user_option == 3) {
      /******time step output for umat_pro_uc_filter ************/
    tout = fopen ("time_output.txt", "a+");
    pout = fopen ("pressure_output.unf", "a+");
    eut_pout = fopen ("eut_pressure_output.unf", "a+");
    Tempout = fopen ("temperature_output.unf", "a+");

    fprintf (tout, "%f\n", current_time);
    fwrite (p1, sizeof (float), nnod, pout);
    fwrite (t1, sizeof (float), nnod, Tempout);

    fclose (tout);
    fclose (pout);
    fclose (eut_pout);
    fclose (Tempout);
      /*******************************************/
  }

   /*****************end of option three *************/
  return;
}                               /* end_of_step */

/*************end of main link routine************/

double howmet (t, h)
     double t, h;

{

  printf ("Howmet Selected!\n");
  printf ("Howmet Selected!\n");
  printf ("Howmet Selected!\n");
  printf ("Howmet Selected!\n");
  printf ("Howmet Selected!\n");
  printf ("Howmet Selected!\n");
  exit (0);
  return (293.);
}

void src_mic_usr ()
{
#ifdef JUNK
  FILE *file_1;
  char output[60];

  register int i;
  int flt;
  extern int mfreq;

  static int FIRST_TIME = 1;
  float dummy[2];
  static float *ras, *dras;

  double q_value[10];

  flt = sizeof (float);

  if (FIRST_TIME) {
    /* dynamic allocation of arrays here */

    ras = float_alloc (n_mic_nod);
    dras = float_alloc (n_mic_nod);

    /* initialization of arrays */

    if (inilev == 0) {
      for (i = 0; i < n_mic_nod; i++)
        ras[i] = 0.8e-3;
    } else {
      /* for restart, read stored data */

      strcpy (output, prefix);
      strcat (output, "m1.unf");

      if ((file_1 = fopen (output, "r")) == NULL) {
        printf ("\n\nUnable to open output file %s\n\n", output);
        exit (1);
      }

      offset = ((inilev - mfreq) * ((n_mic_nod * 6 + 8) * flt) / mfreq);
      fseek (file_1, (long) offset, 0);

      fread (dummy, flt, 2, file_1);
      fread (ras, flt, n_mic_nod, file_1);

      for (i = 0; i < n_mic_nod; i++) {
        dras[i] = 0.;
      }

      fclose (file_1);
    }

    FIRST_TIME = 0;
  }

  for (i = 0; i < nnod; i++) {
    /* insert micro calculations on a nodal basis here */
    /* accumulate fraction solid in fs1 array */
  }

/*   for ( i = 0 ; i < nel ; i++ )
   { */
  /* for each element, load the latent heat released into q_value.
     src_assem assembles this information into the global rhs vector */
/*
      src_assem( i, q_value );
   } */

#endif /*JUNK*/
  return;
}

/*******************************************************************************/

double alcoa (elem, t_ave)

     int elem;
     double t_ave;

{
  register int i;
  double q_vol;

  /*

     The volumetric heat source term in ProCAST has the following form:

     q''' = q_const * time_func * temp_func * space_func

     Since ProCAST already has the time and temp functions built into it,
     this subroutine is used to provide the space function only.

     Setting USER to 4 in the prefixp.dat file will activate the option.

   */

/* The coding between the ****** is for illustration only.
   In this example, the space function modifier is equal to the
   average x coordinate within the element                          */

/* ******************************************************************** */
  q_vol = 0.;

  for (i = 0; i < npe[el_type[elem]]; i++)
    q_vol += x_cord[ncon[elem][i]];

  q_vol /= npe[el_type[elem]];

/* ******************************************************************** */

  return (q_vol);
}

int rolls_royce ()
{
  register int i;
  int node;
  float vtmp, target, tol;
  extern int *sg_list, n_sg;
  extern float *sg;
  extern float vel_mod;

  vel_mod = 1.;

  for (i = 0; i < n_sg; i++) {
    node = sg_list[i];
    if (sg[node] < target + tol) {
      vtmp = (sg[node] - target) / tol;
      if (vtmp < vel_mod)
        vel_mod = vtmp;
    }
  }

  if (vel_mod == 1.)
    vel_mod = 1.1;

  return (0);
}

/********************************************************end of user routine********************/
/************************************************/
/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
/************************************************/
char const *rcs_id_umat_wrapper_c ()
{
  static char const rcsid[] = "EXTERNAL VERSION MAIN ENTRY: $Id: user_rtns.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}

/* end of rcs_id_umat_wrapper_c subroutine */
