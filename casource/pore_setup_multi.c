
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
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include "machine.h"
#include "blocks.h"
#include "pore.h"
#include "props.h"
#include "rand_step.h"
#include "rand_square.h"
#include "gaussdev.h"
#include "umat_histo.h"
#include "constants.h"
#include "pore_routines.h"

/*******************************************/
/* Pore_setup: create the arrays and       */
/* allocate the space for pores, create    */
/* the pore nucleus distribution           */
/*******************************************/
/*$Id: pore_setup_multi.c 1356 2008-08-18 13:41:15Z  $*/
extern int getxyz (int cellnum, int *nc, int *Cell);
extern CA_FLOAT global_pressure;
int pore_histo (BB_struct * bp, int sbnum);

int pore_setup (BB_struct * bp, int sbnum)
{
  CA_FLOAT p_nmax, p_tn, p_tsig, binsize;
  CA_FLOAT inc, v, range, fract_pore, tau, tauold;
  CA_FLOAT tmp1, tmp2, tmp3;
  CA_FLOAT vtot, min_bin;
  SB_struct *sp;
  PORE_str *pl;
  CA_FLOAT gamma;
  CA_FLOAT *fs;
  int *templist, i, numnew, j;
  int errors = 0;
  int newpore, npores, nbins, last, nucmeth;
  FILE *listfile;
  char filename[MAX_STRING_LEN];
  Ctrl_str *cp;

  min_bin = P_NUC_MIN_SS;
  sp = bp->sb[sbnum];
  cp = bp->ctrl;
  p_nmax = bp->pprops.P_nmax;
  gamma = bp->mprops.surf_tens;
  p_tn = bp->pprops.P_tn;
  p_tsig = bp->pprops.P_tsig;
  binsize = bp->pprops.Binsize;
  nucmeth = bp->pprops.P_nmeth;
  fs = bp->sb[sbnum]->c_fs;
  /* for making the list of radius at temprature */

  bp->pprops.P_ntrad = PORE_NTRAD;
#ifdef BUG_NOV_00

  bp->pprops.P_trad_max = TRAD_START;
#else
  bp->pprops.P_trad_max = bp->mprops.Tliq + 10;
#endif
  bp->pprops.P_trad_min = TRAD_END;
  bp->pprops.P_trad_step = (bp->pprops.P_trad_max - bp->pprops.P_trad_min) / (CA_FLOAT) (bp->pprops.P_ntrad);

  fract_pore = tmp1 = tmp2 = 0;
  fprintf (stderr, "Creating pore nucleation threshold list...\n");

  vtot = bp->ncsb * bp->vol_c;
  npores = (int) (p_nmax * vtot);
  fprintf (stderr, "pore_setup: MAX PORES will be: %i\n", npores);
   /*******************************************/
  /* insert the fixed pores   */
   /*******************************************/
  if (cp->fixed_Pore == TRUE) {
    fprintf (stderr, "pore_setup: %i fixed pores\n", cp->nfPore);
    if (npores < cp->nfPore)
      npores = cp->nfPore;
  }

   /*******************************************/
  /* handle bad situations    */
   /*******************************************/
  if (npores >= bp->ncsb && cp->fixed_Pore != TRUE) {
    fprintf (stderr, "ERROR: pore_setup: EXITING: too many pores. Sorry..\n");
    exit (5);
  }

  if (npores == 0) {
    fprintf (stderr, "WARNING: pore_setup: Not enough pores! Creating one pore\n \
so it won't be a total loss!\n");
    npores = 1;
  }

   /*******************************************/
  /* allocate the list for temporary values  */
   /*******************************************/
  if (!(templist = (int *) calloc (bp->ncsb, sizeof (int)))) {
    fprintf (stderr, "ERROR: PL: templist calloc failed\n");
    return (1);
  }

  if (!(pl = (PORE_str *) calloc (npores, sizeof (PORE_str)))) {
    fprintf (stderr, "ERROR: PL_struct malloc failed\n");
    return (1);
  }
   /*******************************************/
  /* handle the various pore nuc distribution */
   /*******************************************/

  /* STEP distribution as in PETER LEE stochastic model */
  /* generalised to use a random function given in options. */
  if (nucmeth == PNUC_STEP || nucmeth == PNUC_GAUSSDEV || nucmeth == PNUC_SQUARE) {

    /* set which random generator to use */
    /* TO DO: select this in initialization routine not here. */
    switch (nucmeth) {
    case PNUC_STEP:
      bp->pprops.pore_rand_function = rand_step;
      break;
    case PNUC_SQUARE:
      bp->pprops.pore_rand_function = rand_square;
      break;

    case PNUC_GAUSSDEV:
      bp->pprops.pore_rand_function = gaussdev;
      break;

    default:
      fprintf (stderr, "ERROR:pore_setup: AAARRGGHH! Non existant option for nuc method: %i\n", nucmeth);
      break;
    }                           /* end set random generator */

    /* first generate and allocate the list of pores in the block */
    for (i = 0; i < npores; i++) {
      CA_FLOAT params[4];       /* store the params to pass to the nuc function */

      /* find a location that doesn't have a pore nuc already */
      do {
        newpore = (int) (drand48 () * bp->ncsb);
      } while (*(templist + newpore) == 1);
      /* define the location of the new pore */
      /* and the beginning of the member cells list */

      if (*(fs + newpore) == NOT_CASTING) {
        pl[i].State = NOT_CASTING;
      } else {
        pl[i].Cellnum = newpore;
        pl[i].boundary = calloc (1, sizeof (p_c_list));
        pl[i].body = calloc (1, sizeof (p_c_list));

        pl[i].boundary->first = new_node (newpore);
        pl[i].boundary->last = pl[i].boundary->first;

        /* put the cell x y z in the pore structure */
        /* pl[i].Cell is an array of [x,y,z] location */
        errors += getxyz (newpore, bp->nc, pl[i].Cell);
        errors += getxyz (newpore, bp->nc, pl[i].boundary->first->Cell);

        /* abort if pore is outside system */
        if (errors) {
          fprintf (stderr, "ERROR:pore_setup: tried to make a pore in Outer Space! %i\n", newpore);
          return (errors);
        }
        /* find the threshold */
        do {
          params[0] = p_tn;
          params[1] = p_tsig;
          params[2] = 0.3333333;
          /* using the pointer to the selected random function .. c syntax ugh */
          pl[i].Thresh = (*(bp->pprops.pore_rand_function)) (params);
        } while ((pl[i].Thresh < P_NUC_MIN_SS) || ((ABS (pl[i].Thresh - p_tn) > (p_tsig * P_NUC_SIG_MULT))));

        /*calculate starting radius and base limiting radius */
        pl[i].Startrad = 2 * gamma / (global_pressure * P_AP_PA * (pl[i].Thresh * pl[i].Thresh - P_AP_ATM));
        tmp1 = (1.0 - bp->pprops.P_limrad_perturb);
        tmp2 = (1.0 + bp->pprops.P_limrad_perturb);
        tmp3 = (tmp2 - tmp1);
        pl[i].base_Limrad = pl[i].Limrad = bp->mprops.das * (tmp1 + drand48 () * tmp3);

        /* allocate the array of t_lists -- to store the state of the pore at temperature benchmarks */
        if (!(pl[i].t_lists = (CA_FLOAT **) calloc (N_T_LISTS, sizeof (CA_FLOAT *)))) {
          fprintf (stderr, "ERROR: PL: t_list %i calloc failed\n", i);
          return (1);
        }

        /* put a flag in the temporary list to indicate a pore is here */
        *(templist + newpore) = 1;
      }                         /* end of NOT_CASTING test */
    }                           /* end --  generate list of pores */

  } else {                      /*other nuc models */
    v = min_bin;
    if (p_tn < min_bin) {
      p_tn = min_bin + 1.0;
      fprintf (stderr, "ERROR:pore_setup: p_tn too low, Setting to %g \n", p_tn);
    }
    range = (p_tn - min_bin) + P_NUC_SIG_MULT * p_tsig;
    nbins = range / binsize;
    last = 0;
    inc = v + binsize;
    fract_pore = 0;

    for (i = 0; i < nbins; i++) {
      /*Find the number of pores in this bin, somehow */
      switch (nucmeth) {

      case PNUC_GAUSS:         /* similar to RAPPAZ for grains */
        tau = v - p_tn;
        tauold = inc - p_tn;

        tmp1 = 1.0 / (SQRT2 * p_tsig);

        tmp1 = 0.5 * p_nmax * (erf ((double) tauold * tmp1) - erf ((double) tau * tmp1));
        tmp2 = tmp1 * bp->vol_c * bp->ncsb;
        tmp2 += (double) fract_pore;
        numnew = floor (tmp2);
        fract_pore = (CA_FLOAT) (tmp2 - (CA_FLOAT) numnew);
        break;

      case PNUC_TRUESTEP:      /*True step distribution */
        if (v > p_tn && inc < p_tn + p_tsig) {
          numnew = (int) ((CA_FLOAT) (npores) * binsize / p_tsig);
        } else {
          numnew = 0;
        }
        break;
      case PNUC_FUNCTION:
        {
          extern int pnuc_function (CA_FLOAT v, CA_FLOAT binsize, int npores, CA_FLOAT p_tsig, CA_FLOAT p_tn, CA_FLOAT par1,
                                    CA_FLOAT par2);

          numnew = pnuc_function (v, binsize, npores, p_tsig, p_tn, bp->pprops.P_par1, bp->pprops.P_par2);

        }
        break;

      case PNUC_INPUT:
        break;
      default:
        break;
      }

      /* avoid seg-fault by excessive pores */
      if ((last + numnew) >= npores) {
        fprintf (stderr, "WARNING: Pore-setup: Too many pores in distribution.Cutting off...\n");
        numnew = 0;
      }
      for (j = 0; j < numnew; j++) {
        do {
          newpore = (int) (drand48 () * bp->ncsb);
        } while (*(templist + newpore) == 1);

        if (*(fs + newpore) == NOT_CASTING) {
          pl[last + j].State = NOT_CASTING;
        } else {
          pl[last + j].Cellnum = newpore;

          /* put the cell x y z in the pore structure */
          errors = getxyz (newpore, bp->nc, pl[last + j].Cell);
          /* abort if pore is outside system */
          if (errors) {
            fprintf (stderr, "ERROR:pore_setup: tried to make a pore in Outer Space! %i\n", newpore);
            return (errors);
          }
          pl[last + j].Thresh = inc;
          pl[last + j].Startrad = 2 * gamma / (global_pressure * P_AP_PA * (pl[last + j].Thresh * pl[last + j].Thresh - P_AP_ATM));
          tmp1 = (1.0 - bp->pprops.P_limrad_perturb);
          tmp2 = (1.0 + bp->pprops.P_limrad_perturb);
          tmp3 = (tmp2 - tmp1);
          pl[i].base_Limrad = pl[i].Limrad = bp->mprops.das * (tmp1 + drand48 () * tmp3);
          /* allocate the array of t_lists */
          if (!(pl[last + j].t_lists = (CA_FLOAT **) calloc (N_T_LISTS, sizeof (CA_FLOAT *)))) {
            fprintf (stderr, "ERROR: PL: t_list %i calloc failed\n", last + j);
            return (1);
          }
          *(templist + newpore) = 1;
        }                       /* end of NOT_CASTING test */
      }                         /*end of loop through potential new pore nuclei (j loop ) */
      v = inc;
      inc += binsize;
      last += numnew;
    }
    if (last < PNUC_MIN) {
      fprintf (stderr, "WARNING: Pore-setup: Not enough pores in distribution. Using random step.\n");
#ifdef ERROR_EXIT
      fprintf (stderr, "ERROR_EXIT: Pore-setup: Actually, I decided to quit instead since ERROR_EXIT was defined.\n");
      exit (0);
#endif /*ERROR_EXIT */

      for (i = 0; i < npores; i++) {
        do {
          newpore = (int) (drand48 () * bp->ncsb);
        } while (*(templist + newpore) == 1);
        pl[i].Cellnum = newpore;
        /* put the cell x y z in the pore structure */
        errors = getxyz (newpore, bp->nc, pl[i].Cell);
        /* abort if pore is outside system */
        if (errors) {
          fprintf (stderr, "ERROR:pore_setup: tried to make a pore in Outer Space! %i\n", newpore);
          return (errors);
        }
        pl[i].Thresh = p_tn + ((drand48 () * p_tsig));
        pl[i].Startrad = 2 * gamma / (global_pressure * P_AP_PA * (pl[i].Thresh * pl[i].Thresh - P_AP_ATM));
        tmp1 = (1.0 - bp->pprops.P_limrad_perturb);
        tmp2 = (1.0 + bp->pprops.P_limrad_perturb);
        tmp3 = (tmp2 - tmp1);
        pl[i].base_Limrad = pl[i].Limrad = bp->mprops.das * (tmp1 + drand48 () * tmp3);
        /* allocate the array of t_lists */
        if (!(pl[i].t_lists = (CA_FLOAT **) calloc (N_T_LISTS, sizeof (CA_FLOAT *)))) {
          fprintf (stderr, "ERROR: PL: t_list %i calloc failed\n", i);
          return (1);
        }
        *(templist + newpore) = 1;
      }
    } else {
      npores = last;
    }

  }                             /*end loop through pore bins(i) */

  /*generate the trheshold for the fixed pores */
  if (bp->ctrl->fixed_Pore == TRUE) {
    for (i = 0; i < bp->ctrl->nfPore; i++) {
      newpore = (int) (bp->ctrl->nPsite[i][0] + bp->ctrl->nPsite[i][1] * bp->nc[0] + bp->ctrl->nPsite[i][2] * bp->nc[0] * bp->nc[1]);
      if (*(fs + newpore) == NOT_CASTING) {
        fprintf (stderr, "WARNING:  pore_setup_multi: Oops! Fixed pore not in casting. %i\n", i);
        pl[i].State = NOT_CASTING;
      } else {
        pl[i].Cellnum = newpore;
        if (pl[i].State == NOT_CASTING) {
          /* pore structure has not been allocated in this case */
          pl[i].State = PORE_NONE;
          pl[i].boundary = calloc (1, sizeof (p_c_list));
          pl[i].body = calloc (1, sizeof (p_c_list));
          pl[i].boundary->first = new_node (newpore);
          pl[i].boundary->last = pl[i].boundary->first;
          /* allocate the array of t_lists -- to store the state of the pore at temperature benchmarks */
          if (!(pl[i].t_lists = (CA_FLOAT **) calloc (N_T_LISTS, sizeof (CA_FLOAT *)))) {
            fprintf (stderr, "ERROR: PL: t_list %i calloc failed\n", i);
            return (1);
          }
        }
        pl[i].boundary->first->cellnum = newpore;
        /* put the cell x y z in the pore structure */
        errors = getxyz (newpore, bp->nc, pl[i].Cell);
        errors += getxyz (newpore, bp->nc, pl[i].boundary->first->Cell);
        /* abort if pore is outside system */
        if (errors) {
          fprintf (stderr, "ERROR:pore_setup: tried to make a pore in Outer Space! %i\n", newpore);
          return (errors);
        }
        pl[i].Thresh = bp->ctrl->nPsite[i][3];
        pl[i].Startrad = 2 * gamma / (global_pressure * P_AP_PA * (pl[i].Thresh * pl[i].Thresh - P_AP_ATM));
        tmp1 = (1.0 - bp->pprops.P_limrad_perturb);
        tmp2 = (1.0 + bp->pprops.P_limrad_perturb);
        tmp3 = (tmp2 - tmp1);
        pl[i].base_Limrad = pl[i].Limrad = bp->mprops.das * (tmp1 + drand48 () * tmp3);
      }                         /* end of NOT_CASTING test */
    }
  }

   /*******************************************/
  /* output the pore threshold list    */
   /*******************************************/
  sp->porelist = pl;
  sp->Npores = npores;

/*Don't create pore list P_L file if STAT_ONLY mode selected*/
#ifndef STAT_ONLY
  sprintf (filename, "P_L_%s.txt", bp->ctrl->fn_base);
  listfile = fopen (filename, "w");
  fprintf (listfile, "Pore List\n");
  fprintf (listfile, "Nuc Method: %i: (0=Gauss, 1=Step)\n", nucmeth);
  fprintf (listfile, "p_nmax\tp_tn\tp_tsig\n");
  fprintf (listfile, "%.5g\t%.5g\t%.5g\n", p_nmax, p_tn, p_tsig);
  fprintf (listfile, "N\tCellnum\tThresh\tStartrad\tbase_Limrad\n");
  for (i = 0; i < npores; i++) {
    fprintf (listfile, "%i\t%i\t%.10g\t%.10g\t%.10g\n", i, pl[i].Cellnum, pl[i].Thresh, pl[i].Startrad, pl[i].base_Limrad);
  }
  fclose (listfile);
#endif /*STAT_ONLY */

  free (templist);
  errors += pore_histo (bp, sbnum);
  return (errors);
}

/*******************************************/
/* pore_histo: set up parameters then call  */
/* umat_histo to make a histogram of the pore */
/* ss thresholds (input)                    */
/*******************************************/
int pore_histo (BB_struct * bp, int sbnum)
{
  int errors = 0;
  int i;
  int npores;
  CA_FLOAT *thr, *thrp;
  FILE *porehistofile;
  PORE_str *pl;
  char filename[MAX_STRING_LEN];
  Histo_struct histoparams;

  pl = bp->sb[sbnum]->porelist;
  npores = bp->sb[sbnum]->Npores;

  sprintf (filename, "P_L_H_sb%i_%s.csv", sbnum, bp->ctrl->fn_base);
  porehistofile = fopen (filename, "w");

  /* set up the histogram parameters */
  histoparams.binsize = P_BINSIZE;
  histoparams.minbin = P_MINBIN;
  histoparams.nbins = P_NBINS;
  histoparams.ndata = npores;

  thr = calloc (npores, sizeof (CA_FLOAT));

  /* collate the threshold values into an array */
  for (thrp = thr, i = 0; thrp < thr + npores; thrp++, i++) {
    *thrp = pl[i].Thresh;
  }

  errors += umat_histo (porehistofile, &histoparams, thr, FIRST_CALL);

  fclose (porehistofile);
  free (thr);
  return (errors);
}

char const *rcs_id_pore_setup_c ()
{
  static char const rcsid[] = "$Id: pore_setup_multi.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}

/*
*/
