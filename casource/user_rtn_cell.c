
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

/*$Id: user_rtn_cell.c 1356 2008-08-18 13:41:15Z  $*/
#include <string.h>
#include <stdio.h>
#include "common.h"
#define STD_TMP 273.16
#define STD_ATM 1.0e-06

#ifdef MYALLOC
extern float * float_alloc(size_t nn);
extern int * int_alloc(size_t nn);
extern double * double_alloc(size_t nn);
#endif

#ifdef EXTERNAL_STUBS
 void cell_findelem (){};   /* external function ?*/
 float cell_findvalue (){return (1.0);}; /* external function ?*/
#else
extern void cell_findelem ();   /* external function ?*/
extern float cell_findvalue (); /* external function ?*/
#endif

/* external variables ??? not in common.h */
extern int Imax, Jmax, Kmax, ncell;
extern float deltaX, deltaY, deltaZ;
extern int *cell_elem; 
extern double Xmin, Ymin, Zmin;
/* end of external variables ??? */

void user_rtn_cell (int total_cell,     /* total number of cells */
                    int *nc,    /* pt to array of max in X, Y, Z */
                    float *origine,     /* array of origin in X, Y ,Z */
                    float *dd_cell,     /* sizeo f a cell */
                    float *umat_t0, float *umat_p0, float *umat_t1, float *umat_p1, float *x, float *y, float *z, int step_num)
{

  register int i, j, k, in, ic;
  char fname[64] = "";
  int mid;
  float xc, yc, zc;
  float value;
  int current_npe, node_number;

#ifdef VERBOSE_EXTERNAL
  FILE *fp_1, *fp_2, *fp_3;
#endif /* VERBOSE_EXTERNAL */

   /*********if the first time step assigne all global variables needed************/
  if (step_num == inilev + 1) {

#ifdef VERBOSE_EXTERNAL
    sprintf (fname, "user_cell_%i.out", step_num);
    fp_1 = fopen (fname, "w+");
    sprintf (fname, "pro_con_%i.out", step_num);
    fp_2 = fopen (fname, "w+");
#endif /* VERBOSE_EXTERNAL */

    ncell = total_cell;
    Imax = nc[0];
    Jmax = nc[1];
    Kmax = nc[2];

    Xmin = origine[0] * 100.0;
    Ymin = origine[1] * 100.0;
    if (TWO_D) {
      Zmin = 0.0;
    } else {
      Zmin = origine[2] * 100.0;
    }

    deltaX = dd_cell[0] * 100.0;
    deltaY = dd_cell[1] * 100.0;
    if (TWO_D) {
      deltaZ = 0.0;
    } else {
      deltaZ = dd_cell[2] * 100.0;
    }

#ifdef VERBOSE_EXTERNAL
    fprintf (fp_1, "\t Imax=%d \t Jmax=%d \t Kmax=%d \t \n 1=%lf \t 2=%lf \t 3=%lf \t \n d1=%lf \t d2=%lf \t d3=%lf \t \n", Imax, Jmax,
             Kmax, origine[0], origine[1], origine[2], dd_cell[0], dd_cell[1], dd_cell[2]);

    fprintf (fp_1,
             "\t Imax=%d \t Jmax=%d \t Kmax=%d \t \n Xmin=%lf \t Ymin=%lf \t Zmin=%lf \t \n deltaX=%lf \t deltaY=%lf \t deltaZ=%lf \t \n",
             Imax, Jmax, Kmax, Xmin, Ymin, Zmin, deltaX, deltaY, deltaZ);
#endif /* VERBOSE_EXTERNAL */

    cell_elem = int_alloc (ncell);
    /* locate the cells in the elements using coordinates stored in global array */
    cell_findelem ();

#ifdef VERBOSE_EXTERNAL
    ic = 0;
    for (k = 0; k < nc[2]; k++) {
      for (j = 0; j < nc[1]; j++) {
        for (i = 0; i < nc[0]; i++) {
          fprintf (fp_1, "%i %i %i %i %i \n", ic, i, j, k, cell_elem[ic]);
          ic++;
        }
      }
    }
#endif /* VERBOSE_EXTERNAL */
    for (ic = 0; ic < ncell; ic++) {
      if (cell_elem[ic] >= 0) {
        mid = mat_id[cell_elem[ic]];

#ifdef VERBOSE_EXTERNAL
        current_npe = npe[el_type[cell_elem[ic]]];
                  /*****print out the node conectivity data**********/
        fprintf (fp_2, "cell_num: %d \t elem_num: %d \n", ic, cell_elem[ic]);
        fprintf (fp_2, "mat_id: %d \t elem_type: %d \t npe: %d \n", mid, el_type[cell_elem[ic]], current_npe);
        for (i = 0; i < current_npe; i++) {
          node_number = ncon[cell_elem[ic]][i];
          fprintf (fp_2, "node_number %d is %d oldtemp %.5g newtemp $.5g \t", i, node_number, t0[node_number], t1[node_number]);
        }
        fprintf (fp_2, "\n \n");
#endif /* VERBOSE_EXTERNAL */
            /*************************************************/

        /* this is where we tell if the element is not in the casting */
        /* */
        if (fluid_state[mid] != 1)
          cell_elem[ic] = -1; /* set cell element to NOT_CASTING */
      }
    }
    if (umat_t0 == NULL || umat_t1 == NULL || umat_p0 == NULL || umat_p1 == NULL) {
      fprintf (stderr, "ERROR:user_rtn_cell: NULL pointers detected!\nExiting\n");
      exit (205);
    } else {
      /* initialize the data exchange arrays to zero */
      memset (umat_t0, 0, sizeof (float) * ncell);
      memset (umat_p0, 0, sizeof (float) * ncell);
      memset (umat_t1, 0, sizeof (float) * ncell);
      memset (umat_p1, 0, sizeof (float) * ncell);
    }
#ifdef VERBOSE_EXTERNAL
    fclose (fp_1);
    fclose (fp_2);
#endif /* VERBOSE_EXTERNAL */
  }
/***************end of first call*************/

#ifdef VERBOSE_EXTERNAL
  sprintf (fname, "external_temp_old%08i.dat", step_num);
  fp_1 = fopen (fname, "w");
  sprintf (fname, "external_temp_new%08i.dat", step_num);
  fp_2 = fopen (fname, "w");
  sprintf (fname, "node_vals_%i.out", step_num);
  fp_3 = fopen (fname, "w+");

  fprintf (fp_1, "VARIABLES = \"Oldtemp\"\n");
  fprintf (fp_1, "ZONE i=%i,j=%i,k=%i\n", nc[0], nc[1], nc[2]);
  fprintf (fp_2, "VARIABLES = \"Newtemp\"\n");
  fprintf (fp_2, "ZONE i=%i,j=%i,k=%i\n", nc[0], nc[1], nc[2]);
  for (ic = 0; ic < ncell; ic++) {
    if (cell_elem[ic] >= 0) {
      current_npe = npe[el_type[cell_elem[ic]]];
                  /*****print out the node conectivity data**********/
      for (i = 0; i < current_npe; i++) {
        node_number = ncon[cell_elem[ic]][i];
        fprintf (fp_3, "node_number %d is %d oldtemp %.5g newtemp %.5g \n", i, node_number, t0[node_number], t1[node_number]);
      }
      fprintf (fp_3, "\n \n");
            /*************************************************/

    }
  }
  fclose (fp_3);
#endif /* VERBOSE_EXTERNAL */

  for (ic = 0; ic < ncell; ic++) {
    if (cell_elem[ic] >= 0) {
      xc = x[ic] * 100.0;
      yc = y[ic] * 100.0;  /***transform the length unit from m to cm*****/
      if (TWO_D) {
        zc = 0.0;
      } else {
        zc = z[ic] * 100.0;
      }

         /***for each cell interpolate the T and P******/
      umat_t0[ic] = cell_findvalue (xc, yc, zc, cell_elem[ic], T, 0, 0);
      umat_t0[ic] -= STD_TMP;
      umat_t1[ic] = cell_findvalue (xc, yc, zc, cell_elem[ic], T, 1, 1);
      umat_t1[ic] -= STD_TMP;

      if (FLOW) {
        umat_p0[ic] = cell_findvalue (xc, yc, zc, cell_elem[ic], P, 1, 0);
        umat_p0[ic] *= STD_ATM;
        umat_p1[ic] = cell_findvalue (xc, yc, zc, cell_elem[ic], P, 1, 1);
        umat_p1[ic] *= STD_ATM;
      } else {
        umat_p0[ic] = umat_p1[ic] = 0.0;
      }
    }
#ifdef VERBOSE_EXTERNAL
    fprintf (fp_1, "%.5g\n", umat_t0[ic]);
    fprintf (fp_2, "%.5g\n", umat_t1[ic]);
#endif /*VERBOSE_EXTERNAL */
  }
#ifdef VERBOSE_EXTERNAL
  fclose (fp_1);
  fclose (fp_2);
#endif /*VERBOSE_EXTERNAL */

  return;
}

/***************************************************************************/
/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *user_rtn_cell_c ()
{
  static char const rcsid[] = "$Id: user_rtn_cell.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}
