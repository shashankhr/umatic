
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
/*      pr_lookup.c:                                            */
/*  All subroutine related to the pressure lookup table         */
/****************************************************************/
/* Written by Robert C. Atwood, Imperial College                */
/* Jul 1, 2000                                                  */
/****************************************************************/
/*       MODIFIED by:                                           */
/****************************************************************/
/*RCS Id:$Id: pr_lookup.c 1341 2008-07-23 15:23:30Z  $*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "machine.h"
#include "fem.h"
#include "read_ctrl.h"
#include "blocks.h"
#include "nuc_lookup.h"

extern CA_FLOAT global_pressure;

/* pressure lookup table */
void init_pr_lookup (BB_struct * bp, Ctrl_str * cp)
{
  FILE *fp, *file_pr;
  static PR_str prlookup;
  PR_str *prl;
  int i, j;
  int *noutput;
  int noutput_new;
  float value1, value2;

  prl = &(prlookup);
  noutput = &(prl->npoints);

  if ((file_pr = fopen ("PrLookup.txt", "w")) == NULL) {
    fprintf (stderr, "ERROR init_prlookup: can't open file\n");
  }
#ifdef READ_PRESSURE
  if ((fp = fopen ("PrLookup.bin", "r")) == NULL) {
    fprintf (stderr, "ERROR init_prlookup: can't open file\n");
  } else {

    fread (noutput, sizeof (int), 1, fp);
    prl->npoints = *noutput;

    /* allocating memory for pressure double array */
    (prl->pr_points) = (float *) calloc (2 * (prl->npoints), sizeof (float));
    /* end of allocation */

    /* reading the time and pressure data from the pressure input file */
    fread (prl->pr_points, sizeof (float), 2 * (prl->npoints), fp);
    /* end of reading */

  }                             /* end of if */
  fclose (fp);
#endif

#ifdef SET_PRESSURE
  if ((fp = fopen ("PrLookup.dat", "r")) == NULL) {
    fprintf (stderr, "ERROR init_prlookup: can't open file\n");
  } else {

    fscanf (fp, "%d \n", &noutput_new);
    prl->npoints = noutput_new;

    /* allocating memory for pressure double array */
    (prl->pr_points) = (float *) calloc (2 * (prl->npoints), sizeof (float));
    /* end of allocation */

    /* reading the time and pressure data from the pressure input file */
    for (i = 0; i <= prl->npoints - 1; i++) {
      fscanf (fp, "%f %f \n", &value1, &value2);
      prl->pr_points[i] = value1;
      prl->pr_points[i + 1] = value2;
      fprintf (file_pr, "%f %f \n", prl->pr_points[i], prl->pr_points[i + 1]);
    }
    /* end of reading */

  }                             /* end of if */
  fclose (fp);
#endif

#ifdef READ_PRESSURE
  for (i = 0; i < prl->npoints - 1; i++) {
    fprintf (file_pr, "%f %f \n", prl->pr_points[i], prl->pr_points[i + 1]);
  }
#endif

  bp->prlookup = prl;

  fclose (file_pr);

  return;

}                               /*end of init_pr_lookup */

void set_global_pressure (BB_struct * bp)
{
  int i, j;
  PR_str *prl;
  CA_FLOAT new_time;

  prl = bp->prlookup;

  new_time = (bp->step) * (bp->delt);

  if (new_time <= prl->pr_points[0]) {
    bp->this_pr = prl->pr_points[1];
  } else if (new_time >= prl->pr_points[2 * (prl->npoints) - 2]) {
    bp->this_pr = prl->pr_points[2 * (prl->npoints) - 1];
  } else {
    for (i = 0; i < 2 * (prl->npoints) - 3; i += 2) {
      if (new_time > prl->pr_points[i] && new_time < prl->pr_points[i + 2]) {
        bp->this_pr =
          prl->pr_points[i + 1] +
          ((prl->pr_points[i + 3] - prl->pr_points[i + 1]) / (prl->pr_points[i + 2] - prl->pr_points[i])) * (new_time -
                                                                                                             prl->pr_points[i]);
      }
    }                           /* end of i loop */
  }

  return;
}

/* end of set_global_pressure */

/* temperature lookup table */
void init_temp_lookup (BB_struct * bp, Ctrl_str * cp)
{
  FILE *fp;
  static TEMP_str templookup;
  TEMP_str *templ;
  int i, j;
  int *noutput;

  templ = &(templookup);
  noutput = &(templ->npoints);

  if ((fp = fopen ("TempLookup.bin", "r")) == NULL) {
    fprintf (stderr, "ERROR init_templookup: can't open file\n");
  } else {

    fread (noutput, sizeof (int), 1, fp);
    templ->npoints = *noutput;

    /* allocating memory for temperature double array */
    (templ->temp_points) = (float *) calloc (2 * (templ->npoints), sizeof (float));
    /* end of allocation */

    /* reading the time and temperature data from the temperature input file */
    fread (templ->temp_points, sizeof (float), 2 * (templ->npoints), fp);
    /* end of reading */

  }                             /* end of if */

  bp->templookup = templ;

  fclose (fp);

  return;

}                               /*end of init_temp_lookup */

void set_global_temperature (BB_struct * bp)
{
  int i, j;
  TEMP_str *templ;
  CA_FLOAT new_time;

  templ = bp->templookup;

  new_time = (bp->step) * (bp->delt);

  if (new_time <= templ->temp_points[0]) {
    bp->this_temp = templ->temp_points[1];
  } else if (new_time >= templ->temp_points[2 * (templ->npoints) - 2]) {
    bp->this_temp = templ->temp_points[2 * (templ->npoints) - 1];
  } else {
    for (i = 0; i < 2 * (templ->npoints) - 3; i += 2) {
      if (new_time > templ->temp_points[i] && new_time < templ->temp_points[i + 2]) {
        bp->this_temp =
          templ->temp_points[i + 1] +
          ((templ->temp_points[i + 3] - templ->temp_points[i + 1]) / (templ->temp_points[i + 2] - templ->temp_points[i])) * (new_time -
                                                                                                                             templ->
                                                                                                                             temp_points
                                                                                                                             [i]);
      }
    }                           /* end of i loop */
  }

  return;
}

/* end od set_global_temperature */

#ifdef OLDJUNK
/* pressure lookup table */
void init_pr_lookup (BB_struct * bp, Ctrl_str * cp)
{
  FILE *fp;
  PR_str prlookup;
  PR_str *prl;
  int i, j;

  prl = &(prlookup);

  if ((fp = fopen ("PrLookup.txt", "r")) == NULL) {
    fprintf (stderr, "ERROR init_prlookup: can't open file\n");
  } else {
    prlookup.pr_init = global_pressure;
    fscanf (fp, "%i\n", &(prlookup.npoints));
    (prl->pr_points) = (CA_FLOAT *) calloc ((1 + prl->npoints), 2 * sizeof (CA_FLOAT));
    *(prl->pr_points) = 0;
    *(prl->pr_points + 1) = global_pressure;

    for (i = 1; i < (prl->npoints + 1); i++) {
      fscanf (fp, "%g,%g\n", (prl->pr_points + 2 * i), (prl->pr_points + 2 * i + 1));
    }
#ifdef PRINT_PR_LOOKUP
    for (i = 0; i < prl->npoints + 1; i++) {
      fprintf (stderr, "%g,%g\n", *(prl->pr_points + 2 * i), *(prl->pr_points + 2 * i + 1));
    }
#endif /*PRINT_PR_LOOKUP */
  }
  bp->prlookup = prl;
  fclose (fp);
}                               /*end of init_pr_lookup */

void set_global_pressure (BB_struct * bp, CA_FLOAT fs)
{
  int i;
  PR_str *prl;

  prl = bp->prlookup;

  i = 0;
  do {
    i += 2;
  } while (fs < *(prl->pr_points + i));

  global_pressure =
    *(prl->pr_points + i) +
    ((*(prl->pr_points + i + 3) - *(prl->pr_points + i + 1)) / (*(prl->pr_points + i + 2) - *(prl->pr_points + i)));
}
#endif /*OLDJUNK*/
/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_pr_lookup_c ()
{
  static char const rcsid[] = "$Id: pr_lookup.c 1341 2008-07-23 15:23:30Z  $";

  return (rcsid);
}

/* end of rcs_id_bigblock_c subroutine */
/*
*/
