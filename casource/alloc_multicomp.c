
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

/*$Id: alloc_multicomp.c 1339 2008-07-23 13:58:29Z  $*/
#include <stdio.h>
#include "machine.h"
#include "blocks.h"
#include "readmat.h"

void free_ms (MultiS_struct * ms, int ele_1)
{
  int i;

  for (i = 0; i < ele_1; i++) {
    free (ms->Diff_matrix_liq[i]);
    free (ms->Diff_matrix_sol[i]);
    free (ms->Clim_multi[i]);
    free (ms->part_coef_matrix[i]);
  }
  free (ms->Cinit_multi);
  free (ms->LDiff_multi);
  free (ms->SDiff_multi);
  free (ms->part_coef_multi);
  free (ms->part_coef_matrix);
  free (ms->Clim_multi);
  free (ms->slope_multi);
  free (ms->Diff_matrix_liq);
  free (ms->Diff_matrix_sol);
  return;
}

void alloc_ms (MultiS_struct * ms, int ele_1, int nc)
{
  int errs = 0;
  int i, j;

   /******************allocating memory to array pointers******************/

  ms->Cinit_multi = (CA_FLOAT *) calloc (ele_1, sizeof (CA_FLOAT));
  if (ms->Cinit_multi == NULL)
    errs++;
  ms->LDiff_multi = (CA_FLOAT *) calloc (ele_1, sizeof (CA_FLOAT));
  if (ms->LDiff_multi == NULL)
    errs++;
  ms->SDiff_multi = (CA_FLOAT *) calloc (ele_1, sizeof (CA_FLOAT));
  if (ms->SDiff_multi == NULL)
    errs++;
  ms->part_coef_multi = (CA_FLOAT *) calloc (ele_1, sizeof (CA_FLOAT));
  if (ms->part_coef_multi == NULL)
    errs++;
  ms->part_coef_matrix = (CA_FLOAT **) calloc (ele_1, sizeof (CA_FLOAT *));
  if (ms->part_coef_matrix == NULL)
    errs++;
  ms->Clim_multi = (CA_FLOAT **) calloc (ele_1, sizeof (CA_FLOAT *));
  if (ms->Clim_multi == NULL)
    errs++;
  ms->slope_multi = (CA_FLOAT *) calloc (ele_1, sizeof (CA_FLOAT));
  if (ms->slope_multi == NULL)
    errs++;

  ms->Diff_matrix_liq = (CA_FLOAT **) calloc (ele_1, sizeof (CA_FLOAT *));
  if (ms->Diff_matrix_liq == NULL)
    errs++;

  ms->Diff_matrix_sol = (CA_FLOAT **) calloc (ele_1, sizeof (CA_FLOAT *));
  if (ms->Diff_matrix_sol == NULL)
    errs++;

  for (i = 0; i < ele_1; i++) {
    ms->Diff_matrix_liq[i] = (CA_FLOAT *) calloc (ele_1, sizeof (CA_FLOAT));
    if (ms->Diff_matrix_liq[i] == NULL)
      errs++;
    ms->Diff_matrix_sol[i] = (CA_FLOAT *) calloc (ele_1, sizeof (CA_FLOAT));
    if (ms->Diff_matrix_sol[i] == NULL)
      errs++;
    ms->Clim_multi[i] = (CA_FLOAT *) calloc (nc, sizeof (CA_FLOAT));
    if (ms->Clim_multi[i] == NULL)
      errs++;
    ms->part_coef_matrix[i] = (CA_FLOAT *) calloc (nc, sizeof (CA_FLOAT));
    if (ms->part_coef_matrix[i] == NULL)
      errs++;
  }

      /*********************************************************/
  /* Set the default values...                         */
      /*********************************************************/
  for (i = 0; i < ele_1; i++) { /*definning the default values for diff coeff *//*initial concentration and partitionning for */
    ms->Cinit_multi[i] = 6.0;   /*solute elements */
    ms->LDiff_multi[i] = 1.0E-08;
    ms->SDiff_multi[i] = 1.0E-13;
    ms->part_coef_multi[i] = 0.13;
    ms->slope_multi[i] = -7.0;
    ms->bin_eut_max[i] = 12;
    ms->bin_eut_temp[i] = 570.0;
    ms->ter_eut_max[i] = 27.0;
    ms->ter_eut_temp[i] = 520.0;
    for (j = 0; j < ele_1; j++) {
      if (j != i) {
        ms->Diff_matrix_liq[i][j] = (ms->LDiff_multi[i] + ms->LDiff_multi[j]) / 20.0;
        ms->Diff_matrix_sol[i][j] = (ms->SDiff_multi[i] + ms->SDiff_multi[j]) / 20.0;
      } else {
        ms->Diff_matrix_liq[i][j] = (ms->LDiff_multi[i] + ms->LDiff_multi[j]) / 2.0;
        ms->Diff_matrix_sol[i][j] = (ms->SDiff_multi[i] + ms->SDiff_multi[j]) / 2.0;
      }
    }
  }
   /*************************************************************************/
  if (errs != 0) {
    fprintf (stderr, "ERROR: alloc_multicomp: Allocation failed %i \n", errs);
  }
  return;
}

/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *alloc_multicomp_c ()
{
  static char const rcsid[] = "$Id: alloc_multicomp.c 1339 2008-07-23 13:58:29Z  $";

  return (rcsid);
}
