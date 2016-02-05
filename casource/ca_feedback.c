
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

#include <stdio.h>
#include <math.h>
#include "common.h"
#include "blocks.h"

#define STD_TMP 273.16
#define STD_ATM 1.0e-06
/*$Id:*/

/* external variables ??? */
extern int Imax, Jmax, Kmax, ncell;
extern float deltaX, deltaY, deltaZ;
extern int *cell_elem;
extern double Xmin, Ymin, Zmin;
/*end external variables ??? */

/* external functions ??? */
extern double cool_mod ();
extern double conc_mod ();
extern double quad_mod ();
extern void src_assem ();
extern void src_lhs ();
extern double temp_mod ();
/*end external functions ??? */

/***************************************************************************/
/******this routine  feedbacks the values of pressure and             ******/
/******fraction solid calculated at microscale into the Procast*************/
/***************************************************************************/

extern float min_temp_calc_extern (BB_struct * bp);
int umat_feedback (BB_struct * bp, int sbnum)
{

  int current_element;
  int *numcell;
  int *elemflag;
  int *nodepor;
  int *numcellpor;
  float *nodefs;
  PORE_str *c_p;
  p_c_node *p_node;

  int offset, node, mat, curve, l;
  register int i, j, k, in, ic;
  float xc, yc, zc;
  int current_npe, node_number;
  float new_fs, old_fs, dfs;
  double q_value[10];
  float value, loc_rho, latent_heat;

#ifdef VERBOSE_EXTERNAL
  FILE *fp_1;
#endif

  numcell = (int *) calloc (nnod, sizeof (int));        /* number of cells the node */
  nodepor = (int *) calloc (nnod, sizeof (int));        /* number pore cells for the node */
  elemflag = (int *) calloc (nel, sizeof (int));        /* flag -- has the element been assebled? */
  nodefs = (float *) calloc (nnod, sizeof (int));       /* averaged fraction solid for the node */

  if (numcell == NULL || nodefs == NULL) {
    fprintf (stderr, "ERROR: umat_feedback: calloc failed! \n");
    exit (938);
  }
#ifdef VERBOSE_EXTERNAL
  fp_1 = fopen ("umat_feedback.out", "w+");
#endif

  /* find the average fs for the cells nearest each node */
  /* accumulate the sum */
  for (i = 0; i < bp->total_cell_number; i++) {
    if (bp->cell_node_array[i] >= 0) {
      node = bp->cell_node_array[i];
      nodefs[node] += (float) (bp->sb[sbnum]->c_fs[i]);
      numcell[node]++;

    }
  }

  if (bp->ctrl->pore) {
    c_p = bp->sb[sbnum]->porelist;
    for (i = 0; i < bp->sb[sbnum]->Npores; i++, c_p++) {
      switch (c_p->State) {
      case PORE_FROZEN:
      case PORE_MULTI:
      case PORE_TUBE:
      case PORE_SPHERE:
        for (p_node = c_p->boundary->first; p_node != NULL; p_node = p_node->next) {
          node = bp->cell_node_array[p_node->cellnum];
          nodepor[node]++;
        }
        break;
      case PORE_NONE:
      case PORE_OFF:
      case PORE_LATENT:
        break;
      default:
        break;
      }
    }
  }

  /* For each node, now ratio by the number of contributing cells */
  for (node = 0; node < nnod; node++) {
    if (numcell[node] > 0) {
#ifdef VERBOSE_EXTERNAL
      fprintf (fp_1, "(before) i-node-fs1-t1:\t %d \t %d \t %f \t %f \n", i, node, fs1[node], t1[node]);
#endif
      nodefs[node] /= (float) (numcell[node]);
      old_fs = fs0[node];
      new_fs = fs1[node];
      /* fs1 is the new fraction solid (global Procast array) */
#ifdef NOFS_EXTERNAL
      fs1[node] = 1;
#else
      fs1[node] = nodefs[node];
#ifdef VERBOSE_EXTERNAL
      if (nodefs[node] > 0) {
        static unsigned int counter = 0;

        if (counter < 50) {
          fprintf (stderr, "Feeding back fraction solid %.5g\n", nodefs[node]);
          counter++;
        }
      }
#endif /*VERBOSE_EXTERNAL */
#endif /*NOFS_EXTERNAL */
      if (bp->ctrl->pore) {
        vf_s1[node] = (float) (nodepor[node]) / (float) (numcell[node]);
      }
      /*fs1[node]=100.0;
         vf_s1[node]=0.1; */
#ifdef VERBOSE_EXTERNAL
      fprintf (fp_1, "         i-node-fs1-t1:\t %d \t %d \t %f \t %f \n", i, node, fs1[node], t1[node]);
#endif
    }
  }
  /* end of fraction solid feedback */

  /* porosity volume feedback */

#ifdef Q_FEEDBACK_ROUTINE
  /* for each element, load the latent heat released into q_value. */
  /* src_assem assembles this information into the global rhs vector */
  for (i = 0; i < ncell; i++) {
    if (bp->cell_element_array[i] >= 0 && elemflag[bp->cell_element_array[i]] == 0) {
      current_element = bp->cell_element_array[i];
      elemflag[current_element] = 1;

      mat = mat_num[mat_id[current_element]];
      latent_heat = lat_heat[mat];

      for (j = 0; j < 10; j++) {
        q_value[j] = (double) 0.0;
      }

      for (j = 0; j < npe[el_type[current_element]]; j++) {
        node = ncon[current_element][j];

        /* evaluate the variation in fraction solid at each node */
        old_fs = fs0[node];
        new_fs = fs1[node];
        if (new_fs > 1.0) {
          dfs = 1.0 - old_fs;
        } else {
          dfs = new_fs - old_fs;
        }
         /**********************************/

        curve = i_dens[mat];
        value = 1.0;

        /* Find the density as a function of temperature */
        /* using the information provided to Procast */
        if (curve > 0) {
          l = curve - 1;
          value = temp_mod (t1[node], l);
        } else if (curve < 0) {
          l = -(curve + 1);
          value = quad_mod (t1[node], l);
        }

        loc_rho = value * density[mat];

        q_value[j] = (double) (loc_rho * latent_heat * dfs / dt);

        /*q_value[j] = 10000000000000.00; */
      }
      src_assem (current_element, q_value);
    }                           /* end of main if loop */
  }                             /* end of cell loop */
#ifdef VERBOSE_EXTERNAL
  for (node = 0; node < nnod; node++) {
    if (numcell[node] > 0) {
      fprintf (fp_1, "(after)  i-node-fs1-t1:\t %d \t %d \t %f \t %f \n", i, node, fs1[node], t1[node]);
    }
  }
#endif /* VERBOSE_EXTERNAL */
#endif /*Q_FEEDBACK_ROUTINE */

#ifdef VERBOSE_EXTERNAL
  fclose (fp_1);
#endif

  free (numcell);
  free (elemflag);
  free (nodefs);
  free (nodepor);
  return (1);
}

/***************************************************************************/
/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *umat_feedback_c ()
{
  static char const rcsid[] = "$Id: ca_feedback.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}
