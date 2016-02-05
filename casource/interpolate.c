
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
/* READPHASE.C:   (Part of CA code)                             */
/* Subroutine to interpolate values from phase diagram          */
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee, Robert C. Atwood and A. Chirazi     */
/*                                            Imperial College */
/* Mar 01, 2001                                                 */
/****************************************************************/

/****************************************************************/
/* Versions maintained with RCS                                 */
/* Version 1.0: Aug. 13, pdl                                    */
/****************************************************************/
/*RCS id $Id: interpolate.c 1339 2008-07-23 13:58:29Z  $*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "machine.h"
#include "read_ctrl.h"
#include "readmat.h"
#include "blocks.h"
#include "multi_diff_props.h"

/********************a degree N polynominal interpolation used in bilinear interpolate function************/

void polint (CA_FLOAT * xa, CA_FLOAT * ya, int n, CA_FLOAT x, CA_FLOAT * y, CA_FLOAT * dy)
{
  int i, m, ns = 1;
  CA_FLOAT den, dif, dift, ho, hp, w;
  CA_FLOAT *c, *d;
  FILE *temp;

/*  temp=fopen("liquidus1.out","a+"); */

  dif = fabs (x - xa[1]);
 /*************allocate local array***********/

  c = (CA_FLOAT *) calloc (n + 1, sizeof (CA_FLOAT));
  d = (CA_FLOAT *) calloc (n + 1, sizeof (CA_FLOAT));

 /********************************************/

  for (i = 1; i <= n; i++) {
    if (xa[i] != 0.0) {
      if ((dift = fabs (x - xa[i])) < dif) {
        ns = i;
        dif = dift;
      }
      c[i] = ya[i];
      d[i] = ya[i];
    } else {
      continue;
    }
  }
  for (m = 1; m < n; m++) {
    for (i = 1; i <= n - m; i++) {
      ho = xa[i] - x;
      hp = xa[i + m] - x;
      w = c[i + 1] - d[i];
      if ((den = ho - hp) == 0.0) {
        fprintf (stderr, "error in routine polint");
      }
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    *y += (*dy = (2 * ns < (n - m) ? c[ns + 1] : d[ns--]));
  }
/*   fprintf(temp,"%f \t %f\n",x, *y); */
/* fclose(temp); */

  free (c);
  free (d);
} /************end of polint subroutine*************/

/**********************bilinear interpolation function*****************/

void polint2 (CA_FLOAT * xa1, CA_FLOAT * xa2, CA_FLOAT ** ya, int m, int n, CA_FLOAT x1, CA_FLOAT x2, CA_FLOAT * y, CA_FLOAT * dy)
{
  int j;
  CA_FLOAT *ymtmp;

 /*****************allocate local array**********************/

  ymtmp = (CA_FLOAT *) calloc (m + 1, sizeof (CA_FLOAT));

/**********************************************************/

  for (j = 1; j <= m; j++) {
    polint (xa2, &ya[j][0], n, x2, &ymtmp[j], dy);
  }
  polint (xa1, ymtmp, m, x1, y, dy);

  free (ymtmp);

}

/******************end of polint2 subroutine**********************/

int interpolate (Ctrl_str * cp, BB_struct * bp, int cell_index, CA_FLOAT * cell_temp, int sbnum, int inter_flag, CA_FLOAT Tcell)
{
  char line[MAX_STRING_LEN];
  char *token;
  int i, j, k, l, h, itmp, n, npaded;
  int table_index, counter, initial;
  int numtietri;
  int nx, ny, nz, nc;
  int rflag = 0;
  int index = 0;
  int error = FALSE;
  Mat_str *mp;
  Nuc_str *np;
  P_str *pp;
  MultiS_struct *ms;
  SB_struct *sp;
  int ele_num, ele_1;
  CA_FLOAT AA, BB, CC;
  CA_FLOAT *Cinit_multi_ptr;
  CA_FLOAT **Clim_multi_ptr;
  FILE *fileout;
  CA_FLOAT yy, dyy;
  CA_FLOAT conc[5];
  CA_FLOAT xx1, xx2, xx3;
  CA_FLOAT lxx1, lxx2, lxx3;
  CA_FLOAT local_max[10];

  CA_FLOAT y1, y2, y3, y4, t, u;

  FILE *temp_new;

/*   temp_new=fopen("liquidus.out","w"); */

  ms = &(bp->MultiSvals);       /* local variable for the multi diff values */
  pp = &(bp->pprops);
  mp = &(bp->mprops);
  np = &(bp->nprops);
  sp = bp->sb[sbnum];

  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  nc = bp->ncsb;
  npaded = (nx + 2) * (ny + 2) * (nz + 2);
  bp->ntsb = bp->nsb[0] * bp->nsb[1] * bp->nsb[2];

  /*from read_ctrl.h */
  ele_num = D_NUM_COMP;
   /*******/
  ele_num = cp->NUM_COMP;       /*number of components in the alloy */
  ele_1 = ele_num - 1;

   /**************the bilinear interpolation scheme is used********/
   /**************three dimensional parameter space is considered**/
   /**************liquidus temperature is interpolated based on the solute***/
   /**************concentration and solute concentration for each component**/
   /**************is interpolated based on other solute values and the liquidus temperature***/

/*******************************begin interpolation***********************/

  numtietri = ms->numtietri;
                           /***number of lines in the table*****/

/***decide which interpolation is needed, liquidus or conc****/

  switch (inter_flag) {
  case LIQUIDUS:

    if (cp->interpolate == LINEAR2) {
      for (i = 0; i < ele_1; i++) {
        conc[i] = sp->c_sol_alloy_multi[i][cell_index];
      }

 /***************call the bilinear interpolation routine*******************/
      polint2 (&(ms->xxa[3][0]), &(ms->xxa[4][0]), ms->ya1, ms->counter[3], ms->counter[4], conc[0], conc[1], cell_temp, &dyy);
    }
    /* end of linear2 */
    if (cp->interpolate == LINEAR1) {

      for (i = 0; i < ele_1; i++) {
        conc[i] = sp->c_sol_alloy_multi[i][cell_index];
      }

      for (i = 1; i < ms->counter[3]; i++) {
        for (j = 1; j < ms->counter[4]; j++) {
          if ((conc[0] >= ms->xxa[3][i]) && (conc[0] <= ms->xxa[3][i + 1]) && (conc[1] >= ms->xxa[4][j])
              && (conc[1] <= ms->xxa[4][j + 1])) {
            y1 = ms->ya1[i][j];
            y2 = ms->ya1[i + 1][j];
            y3 = ms->ya1[i + 1][j + 1];
            y4 = ms->ya1[i][j + 1];

            t = (conc[0] - ms->xxa[3][i]) / (ms->xxa[3][i + 1] - ms->xxa[3][i]);
            u = (conc[1] - ms->xxa[4][j]) / (ms->xxa[4][j + 1] - ms->xxa[4][j]);

            *cell_temp = (1 - t) * (1 - u) * y1 + t * (1 - u) * y2 + t * u * y3 + (1 - t) * u * y4;
          }
        }
      }

    }
    /* end of linear1 */
    if (cp->interpolate == REGRESSION) {
      xx1 = (sp->c_sol_alloy_multi[0][cell_index]);
      xx2 = (sp->c_sol_alloy_multi[1][cell_index]);
      lxx1 = log (xx1);
      lxx2 = log (xx2);

/* if (Tcell>ms->bin_eut_temp[0]){ */
      *cell_temp = 672.37 - (7.648 * xx1) - (3.814 * xx2);
/*  }else{
  *cell_temp=exp(6.153+(0.349*lxx1)-(0.114*lxx1*lxx1)-(0.02320*lxx2*lxx2)+(0.01774*lxx1*lxx2));
 } */

    }

    /* end of regression */
    /*  fprintf(temp_new,"%f \t %f \t %f \n",conc[0],conc[1],*cell_temp);  */
    break;

  case CONC_MULTI:

    if (cp->interpolate == LINEAR2) {
      polint2 (&(ms->xxa[1][0]), &(ms->xxa[3][0]), ms->ya3, ms->counter[1], ms->counter[3], *cell_temp,
               sp->c_sol_alloy_multi[0][cell_index], &yy, &dyy);
      ms->Clim_multi[1][cell_index] = yy;
      polint2 (&(ms->xxa[1][0]), &(ms->xxa[4][0]), ms->ya2, ms->counter[1], ms->counter[4], *cell_temp,
               sp->c_sol_alloy_multi[1][cell_index], &yy, &dyy);
      ms->Clim_multi[0][cell_index] = yy;
    }
    /* end of linear2 */
    if (cp->interpolate == LINEAR1) {

      conc[0] = *cell_temp;
      conc[1] = sp->c_sol_alloy_multi[1][cell_index];

      for (i = 1; i < ms->counter[1]; i++) {
        for (j = 1; j < ms->counter[4]; j++) {
          if ((conc[0] >= ms->xxa[1][i]) && (conc[0] <= ms->xxa[1][i + 1]) && (conc[1] >= ms->xxa[4][j])
              && (conc[1] <= ms->xxa[4][j + 1])) {
            y1 = ms->ya2[i][j];
            y2 = ms->ya2[i + 1][j];
            y3 = ms->ya2[i + 1][j + 1];
            y4 = ms->ya2[i][j + 1];

            t = (conc[0] - ms->xxa[1][i]) / (ms->xxa[1][i + 1] - ms->xxa[1][i]);
            u = (conc[1] - ms->xxa[4][j]) / (ms->xxa[4][j + 1] - ms->xxa[4][j]);

            ms->Clim_multi[0][cell_index] = (1 - t) * (1 - u) * y1 + t * (1 - u) * y2 + t * u * y3 + (1 - t) * u * y4;
          }
        }
      }

      conc[0] = *cell_temp;
      conc[1] = sp->c_sol_alloy_multi[0][cell_index];

      for (i = 1; i < ms->counter[1]; i++) {
        for (j = 1; j < ms->counter[3]; j++) {
          if ((conc[0] >= ms->xxa[1][i]) && (conc[0] <= ms->xxa[1][i + 1]) && (conc[1] >= ms->xxa[3][j])
              && (conc[1] <= ms->xxa[3][j + 1])) {
            y1 = ms->ya3[i][j];
            y2 = ms->ya3[i + 1][j];
            y3 = ms->ya3[i + 1][j + 1];
            y4 = ms->ya3[i][j + 1];

            t = (conc[0] - ms->xxa[1][i]) / (ms->xxa[1][i + 1] - ms->xxa[1][i]);
            u = (conc[1] - ms->xxa[3][j]) / (ms->xxa[3][j + 1] - ms->xxa[3][j]);

            ms->Clim_multi[1][cell_index] = (1 - t) * (1 - u) * y1 + t * (1 - u) * y2 + t * u * y3 + (1 - t) * u * y4;
          }
        }
      }

    }
    /* end of linear1 */
    if (cp->interpolate == REGRESSION) {
      xx1 = (*cell_temp);
      lxx1 = log (xx1);
/*  ms->Clim_multi[0][cell_index]=exp(2727.3-(840.79*lxx1)-(25.04*lxx3)+(64.84*lxx1*lxx1)+(3.898*lxx1*lxx3)-(0.01413*lxx3*lxx3)); */
      ms->Clim_multi[0][cell_index] = xx2 = exp (1631.0 - (504.64 * lxx1) + (39.08 * lxx1 * lxx1));
      ms->Clim_multi[1][cell_index] =
        2896.4 - (8.833 * xx1) - (52.59 * xx2) + (0.00674 * xx1 * xx1) + (0.08078 * xx1 * xx2) + (0.202 * xx2 * xx2);
    }
    /* end of regression */
    break;

  case CONC_MULTI_MONO:

    if (cp->interpolate == REGRESSION) {
      xx1 = (*cell_temp);
      lxx1 = log (xx1);

/*  ms->Clim_multi[0][cell_index]=exp(-2380.2+(758.81*lxx1)-(0.216*lxx3*lxx3)-(60.42*lxx1*lxx1)+(0.08243*lxx1*lxx3));
  ms->Clim_multi[1][cell_index]=exp(14923.2-(4706.2*lxx1)+(370.09*lxx1*lxx1)-(8.353*lxx2*lxx2)+(36.5*lxx2)); */

      ms->Clim_multi[0][cell_index] = exp (-53.53 + (8.814 * lxx1));
      ms->Clim_multi[1][cell_index] = exp (-9228.9 + (2947.2 * lxx1) - (235.21 * lxx1 * lxx1));

    } else {
      fprintf (stderr, "option 3 must be used for the interpolation scheme");
    }                           /* end of regression */

    break;

  }
  /*********end of switch routine************/
/********check for the correct limiting value********/
  for (i = 0; i < ele_1; i++) {
    if (ms->Clim_multi[i][cell_index] < ms->Cinit_multi[i] && *cell_temp > ms->bin_eut_temp[1]) {
      ms->Clim_multi[i][cell_index] = ms->Cinit_multi[i];
    }
  }

  for (i = 0; i < ele_1; i++) {
/* fprintf(temp_new,"%f \t %f \t %f \n", *cell_temp,ms->Clim_multi[0][cell_index],ms->Clim_multi[1][cell_index]); */
  }

/* fclose(temp_new); */
  return (1);

}

/***end of interpolate subroutine***/
/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_interpolate_c ()
{
  static char const rcsid[] = "$Id: interpolate.c 1339 2008-07-23 13:58:29Z  $";

  return (rcsid);
}

/* end of rcs_id_interpolate_c subroutine */
