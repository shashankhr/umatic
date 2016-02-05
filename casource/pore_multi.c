
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
#include "pore_routines.h"
#include "props.h"
#include "rand_step.h"
#include "gaussdev.h"
#include "umat_histo.h"

extern CA_FLOAT find_new_radius (BB_struct * bp, CA_FLOAT n, CA_FLOAT p, CA_FLOAT tempK, CA_FLOAT r);
extern CA_FLOAT global_pressure;
extern double cbrt (double);
extern int getxyz (int cellnum, int *nc, int *Cell);
extern int cell_find_nmols (BB_struct * bp, int sbnum, PORE_str * c_p, CA_FLOAT * fs, CA_FLOAT * sol, CA_FLOAT * alloy, int Cell[3],
                            int cellnum);

void trad_pore (BB_struct * bp, PORE_str * c_p, int trad_index);

/* sb_pore routine : nucleate or grow a pore and deduct the hydrogen used from the */
/* gas solute array*/
int sb_pore (BB_struct * bp, int sbnum)
{

  int errors = 0, ngrowing = 0;
  CA_FLOAT cell_temp, cell_tempK, cell_h, cell_si, cell_satSI, pore_rad, p_vol, pressure;
  CA_FLOAT ss, mat, logsl, p_a, gamma, nmols, new_nmols, new_rad;
  CA_FLOAT cell_fs, cell_fl, cell_sol, min_fs;
  CA_FLOAT pore_pressure, pore_sat;
  CA_FLOAT trad_max, trad_step;
  CA_FLOAT cap_vol, tube_length, extra_vol, eq_rad;
  CA_FLOAT test_vol, new_volume;
  int *nni, *nnip, *nniend;
  int trad_index, trad_last;
  int cellnum, min_cell;
  int dumb;
  int i, j, k, n;
  int nx, ny, nz;
  SB_struct *sp;
  CA_FLOAT Tliq, Tavg, Tunder, prob, extra_h, extra_nmols, f_por, phi_cell;
  CA_FLOAT pore_pressure_ideal, pore_pressure_gamma;
  CA_FLOAT *fs, *sol, *alloy, *c_temp_p;
  CA_FLOAT **c_t_lists;         /* current t_lists at benchmarks */
  int cell_temp_on, npores;
  PORE_str *pl, *c_p;
  p_c_node *node;

  sp = bp->sb[sbnum];           /*subblock pointer */
  nni = nnip = bp->nbhd.nnq;    /* non padded neighbour lookup */
  nniend = nni + 6;
#  ifdef NOSOLID
  /* stop processing the pores once sb is fully solid */
  if (sp->done == TRUE) {
    return (errors);
  }
#  endif /*NOSOLID*/
    /* some physical constants */
    gamma = bp->mprops.surf_tens;
  p_a = P_AP_PA;

  /*set up correct options according to control file */
  if (bp->ctrl->scheil == TRUE) {
    fs = sp->sch_fs;            /* use schiel calculated fraction solid */
  } else {
    fs = sp->c_fs;              /* use CA simulated fraction solid */
  }

  /* dimensions of subblok */
  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];

  /* how many benchmarks for storing the results */
  trad_max = bp->pprops.P_trad_max;
  trad_step = bp->pprops.P_trad_step;
  /* some local pointers into bigblock structure */

  sol = sp->c_sol;              /*cell hydrogen solute array */
  c_temp_p = sp->c_temp;
  alloy = sp->c_sol_alloy;      /*cell alloy solute arrray */
  pl = sp->porelist;            /* array of pore structures */
  npores = sp->Npores;          /*number of pores in subblock */

/**************************************************************/
   /*******************************/
  /* start looping through pores */
   /*******************************/
  for (i = 0; i < npores; i++) {
    /* set up some local variables for index i */
    c_p = pl + i;               /* the pore structure for the Current Pore */
      /**************************************************************************/
    cellnum = c_p->Cellnum;     /* the cell number of the home cell for pore */
    if (c_p->State == NOT_CASTING) {
      /* do nothing */
      int dumb;

      dumb = 0;
    } else {

     #ifdef TRAP_PORE
       {
          int dumb;
          if (i == TRAP_PORE ){
              if (c_p->State != 0){
              dumb = i;
              }
          }
       }
     #endif
      c_t_lists = c_p->t_lists; /* the current array of data lists for output */

      /* temperature -- needed here for assemblin the trad data */
      cell_temp = *(c_temp_p + cellnum);
      cell_tempK = cell_temp + STD_TEMP;        /* in Kelvin */

      /* find the index of the temperature benchmark */
      trad_index = (int) (FLOOR ((trad_max - cell_temp) / trad_step));

      /*warn if data overruns the arrays */
      if (trad_index >= bp->pprops.P_ntrad) {
        static int tradmsg = 0;

#           ifdef PORE_REALLOC
        int k;

        /*expand the arrays to store data for the pores */
        fprintf (stderr, "PORE_MULTI: ifdef PORE_REALLOC: Expanding pore temp. radius data (TRAD) lists from %i ", bp->pprops.P_ntrad);
        bp->pprops.P_ntrad += PORE_EXTRA;
        fprintf (stderr, "to %i\n", bp->pprops.P_ntrad);
        for (k = 0; k < npores; k++) {  /* doit for all pores */
          for (j = 0; j < N_T_LISTS; j++) {
            if (!((pl + k)->t_lists[j] = (CA_FLOAT *) realloc ((pl + k)->t_lists[j], bp->pprops.P_ntrad * sizeof (CA_FLOAT)))) {
              fprintf (stderr, "ERROR: PORE %i: t_list %i realloc failed.\n", k, j);
#                          ifdef PORE_EXIT
              fprintf (stderr, "pore_multi: ifdef PORE_EXIT: therefore finishing with signal 15\n");
              /* try to finish off nicely by raising signal  */
              /* to allow final output before crash and burn */
              tradmsg = raise (SIGTERM);
              if (tradmsg < 0) {        /* did not work, so crash! 8-< */
                fprintf (stderr, "pore_multi: ifdef PORE_EXIT: Signal failed. Exiting....\n");
                exit (tradmsg);
              }
#                          endif /*PORE_EXIT */
              return (1);
            }
          }
        }
#               else /* not PORE_REALLOC */

        /*and just put the results in zero */
        if (tradmsg < MAX_ERR_MSG) {
          fprintf (stderr, "ERROR: sb_pore: trad_index too big! %i\n", trad_index);
          tradmsg++;
#              ifdef PORE_EXIT
          fprintf (stderr, "pore_multi: ifdef PORE_EXIT: Finishing with signal 15\n");
          /* try to finish off nicely by raising signal  */
          /* to allow final output before crash and burn */
          tradmsg = raise (SIGTERM);
          if (tradmsg < 0)      /* did not work, so crash! 8-< */
            exit (tradmsg);
#              endif /*PORE_EXIT */
        }
        trad_index = 0;
#           endif /*PORE_REALLOC */
      }

      cell_fs = (*(fs + cellnum));      /* the fraction solid of the pore's home cell */

      /* set the initial Limrad */
      c_p->Limrad = c_p->base_Limrad * ((1 - cell_fs) / (2 * DAS_COS_THETA));   /* the limiting radius for the current pore */

         /************************************************/
      /*Nucleate or grow the pore depending upon state */
         /************************************************/
      switch (c_p->State) {
      case PORE_NONE:
        {       /***************************************/
          /* pore has not nucleated yet, check it */
                /***************************************/

          /* Nucleate a pore if ss is greatr than the threshold */
          /* but not if fs is greater than eutectic */

          /* find the supersaturation */
          /** \todo Needs updating for polyphase -- how to decide a minimum liquid fraction? */
          c_p->FracSol = 0;
          errors += cell_find_nmols (bp, sbnum, c_p, fs, sol, alloy, c_p->Cell, c_p->Cellnum);
          ss = c_p->supersat;

          if ((ss > c_p->Thresh) && (cell_fs < bp->mprops.alloyprops[0].Fs_eut) && (cell_fs < PORE_MAX_FS) && (cell_fs != NOT_CASTING)) {

                /******************/
            /* nucleate a pore */
                /******************/
            {                   /* protect from excessive messages */
              static int porenuc_msg = 0;

              if (porenuc_msg < MAX_NUC_MSG) {
                fprintf (stderr, "Nucleating pore %i at cell %i\n", i, c_p->Cellnum);
                fprintf (stderr, "Porenuc:cellnum,%i,cell_fs,%.10g,cell_temp,%.10g,,ss,%.10g\n", cellnum, cell_fs, cell_temp, ss);
                porenuc_msg++;
              }
            }                   /* end message protection */
            /* set all the initial values */
#ifdef USE_TUBE
            if (c_p->Startrad <= c_p->Limrad) {
              /* Pore can grow as Sphere initially */
              c_p->Radius = c_p->Oldrad = c_p->Startrad;
              c_p->State = PORE_SPHERE;
            } else {
              /* Pore is already a Tube */
              c_p->Radius = c_p->Oldrad = c_p->Limrad;
              c_p->State = PORE_TUBE;
            }
#else
            /* Pore can grow as Sphere initially */
            c_p->Radius = c_p->Oldrad = c_p->Startrad;
            c_p->State = PORE_SPHERE;
            c_p->EqRad = c_p->Radius;
#endif

            p_vol = FOURPI_BY_THREE * c_p->Radius * c_p->Radius * c_p->Radius;
            pressure = (global_pressure * p_a) + (2 * gamma) / (c_p->Radius);

            c_p->NucSS = ss;
            c_p->NucTemp = cell_temp;
            c_p->Temp = cell_temp;
            c_p->Volume = p_vol;
            c_p->Pressure = pressure;
            c_p->Nmols = (pressure * p_vol) / (GAS_CONST_SI * cell_temp + STD_TEMP);
            c_p->NmolsH = 0;
            c_p->extra_nmols = 0;
            c_p->ncells = 1;

            /*create the arrays to store data for the pores */
            /* and set all to zero by default.             */
            for (j = 0; j < N_T_LISTS; j++) {

              if (!(c_p->t_lists[j] = (CA_FLOAT *) calloc (bp->pprops.P_ntrad, sizeof (CA_FLOAT)))) {
                fprintf (stderr, "ERROR: PORE %i: t_list %i calloc failed.\n", i, j);
                return (1);
              }
            }

            /* increment the number of active pores */
            sp->sb_npores++;
            bp->bb_npores++;

            /* initialise the temp. rad. benchmark */
            /* last benchmark */
            c_p->trad_last = trad_index;
          }                     /* end of nucleate a pore section */
        }                       /*end of PORE_NONE case */
        break;

     /********************************************************/
        /* grow as TUBE with limiting radius, constant pressure */
     /********************************************************/
      case PORE_TUBE:
      case PORE_SPHERE:        /* fall through as cases are now identical */
      case PORE_MULTI:
        /* the pore is going to occupy more than one cell */

        /* paranoia check */
        if (c_p->extra_nmols != 0) {
          errors++;
          fprintf (stderr, "ERROR:pore_multi: Extra nmols still not zeroed. Debug.\n");
          c_p->extra_nmols = 0;
        }

        test_vol = bp->vol_c * c_p->ncells;

        if (c_p->ncells < 1) {
          /* this should not happen */
          fprintf (stderr, "ERROR:pore_multi: No cells in pore ???, %i \n", c_p->Cellnum);
          c_p->State = PORE_OFF;

#     ifdef ERROR_EXIT
          exit (5);
#     endif /*ERROR_EXIT */

        } else {                /* end of paranoia checks */
         /****************************************************/
          /* then the boundary list needs to be traversed to: */
          /*      get the extra h                             */
          /*      decide where to expand                      */
          /*      check if all boundary cells are surrounded by solid */
          /*      possibly check for an upper limit on the size */
          /*      or become spherical                          */
         /*****************************************************/
          /* cell_find_nmols also sums up the fraction solid in the listed cells */
          c_p->FracSol = 0;
          for (node = c_p->boundary->first; node != NULL; node = node->next) {
            errors += cell_find_nmols (bp, sbnum, c_p, fs, sol, alloy, node->Cell, node->cellnum);
            if (bp->ctrl->flow_on == TRUE) {
              global_pressure = bp->ctrl->ref_pres + bp->current_cell_pres[node->cellnum];
            }
            *(sol + node->cellnum) = c_p->sat;
          }
          /* get the average Fs for the pore boundary */
          c_p->FracSol /= c_p->ncells;

          /* converted to mols of H2 molecules */
          new_nmols = c_p->Nmols + (0.5 * c_p->extra_nmols);

          /* update the Limrad limiting radius for the pore */
          c_p->Limrad = MAX (((1 - c_p->FracSol) * c_p->base_Limrad) / (2 * DAS_COS_THETA), MIN_LIMRAD);

          /* find the new pressure and volume, needed for testing for expansion */
          if ((c_p->Radius < c_p->Limrad) || (c_p->FracSol <= PORE_MIN_TUBE_FS)) {
            /*find the new radius by Newton's method */
            c_p->State = PORE_SPHERE;
            new_rad = find_new_radius (bp, new_nmols, global_pressure * P_AP_PA, cell_tempK, c_p->Radius);
            c_p->Pressure = (global_pressure * P_AP_PA) + (2 * gamma / new_rad);
            c_p->Volume = new_nmols * GAS_CONST_SI * cell_tempK / c_p->Pressure;
            c_p->EqRad = cbrt (c_p->Volume * THREE_BY_4PI);
          } else {

            c_p->State = PORE_TUBE;
            new_rad = c_p->Limrad;
            /* assume tubular growth */
            pore_pressure_gamma = (global_pressure * p_a) + (2 * gamma) / (c_p->Limrad);
            /* find new volume using new pressure */
            pore_pressure_ideal = new_nmols * GAS_CONST_SI * cell_tempK / c_p->Volume;
            /* PETER LEE's test for growth */
            if (pore_pressure_ideal >= pore_pressure_gamma) {
              /* pore is growing */
              c_p->Pressure = pore_pressure_gamma;
              c_p->Volume = new_nmols * GAS_CONST_SI * cell_tempK / c_p->Pressure;
              c_p->EqRad = cbrt (c_p->Volume * THREE_BY_4PI);
            } else {
              /* hydrogen is building up in the pore */
              c_p->Pressure = pore_pressure_ideal;
              /* don't need to change the volume */
            }
          }

          /* update the values in the pore structure */
          c_p->Oldrad = c_p->Radius;
          c_p->Radius = c_p->Newrad = new_rad;
          c_p->Nmols = new_nmols;
          c_p->NmolsH += (0.5 * c_p->extra_nmols);
          c_p->Temp = cell_temp;

         /*****************************************/
          /* test if pore can expand into new cell */
         /*****************************************/
          if (c_p->Volume > test_vol) {
            /* expand into a new cell */
            /* find the cell */
            min_cell = find_new_cell (bp, fs, sol, c_p);
            if (min_cell == -1) {
              /* cannot add a cell - no liquid cells outside the pore */

            } else {
              /* make a node and add it to the boundary */
              node = new_node (min_cell);
              getxyz (min_cell, bp->nc, node->Cell);
              add_to_list (c_p->boundary, node);
              c_p->ncells++;
            }
          }

          /* zero out the extra nmols */
          c_p->extra_nmols = 0;

        }                       /* end of normal situation ( 1 or more cells in pore! ) */

        /*switch off if all H is gone */
        if (c_p->NmolsH <= 0.0) {
          {                     /* protect from excessive messages */
            static int switchoff_msg = 0;

            if (switchoff_msg < MAX_NUC_MSG) {
              fprintf (stderr, "Pore: switched off %i cell %i radius %.5g\n", i, c_p->Cellnum, c_p->Radius);
              switchoff_msg++;
            }
          }                     /* end protect messages */
          c_p->State = PORE_OFF;
          c_p->Radius = 0.0;
          c_p->Nmols = 0.0;
          c_p->NmolsH = 0.0;
          c_p->Volume = 0.0;
          c_p->Pressure = 0.0;
          c_p->EqRad = 0.0;
          c_p->FracSol = 0.0;

          /* decrement the number of active pores */
          sp->sb_npores--;
          bp->bb_npores--;
        }

        /*end turn off pore */
        /* if all cells are solid (or possibly nearly solid?) */
        if (c_p->FracSol >= PORE_MAX_FS) {
          c_p->State = PORE_FROZEN;
        }
        break;

      case PORE_FROZEN:
        c_p->FracSol = 0;
        c_p->extra_nmols = 0;
        for (node = c_p->boundary->first; node != NULL; node = node->next) {
          errors += cell_find_nmols (bp, sbnum, c_p, fs, sol, alloy, node->Cell, node->cellnum);
          *(sol + node->cellnum) = c_p->sat;
        }
        c_p->FracSol /= c_p->ncells;
        /* converted to mols of H2 molecules */
        new_nmols = c_p->Nmols + (0.5 * c_p->extra_nmols);
        c_p->Pressure = new_nmols * GAS_CONST_SI * cell_tempK / c_p->Volume;
        c_p->Nmols = new_nmols;
        c_p->NmolsH += (0.5 * c_p->extra_nmols);
        c_p->Temp = cell_temp;
        /*switch off if all H is gone */
        if (c_p->NmolsH <= 0.0) {
          {                     /* protect from excessive messages */
            static int switchoff_msg = 0;

            if (switchoff_msg < MAX_NUC_MSG) {
              fprintf (stderr, "Pore: switched off %i cell %i radius %.5g\n", i, c_p->Cellnum, c_p->Radius);
              switchoff_msg++;
            }
          }                     /* end protect messages */
          c_p->State = PORE_OFF;
          c_p->Radius = 0.0;
          c_p->Nmols = 0.0;
          c_p->NmolsH = 0.0;
          c_p->Volume = 0.0;
          c_p->Pressure = 0.0;
          c_p->EqRad = 0.0;
          c_p->FracSol = 0.0;

          /* decrement the number of active pores */
          sp->sb_npores--;
          bp->bb_npores--;
        }                       /*end turn off pore */
        break;

      case PORE_OFF:
        /* do nothing */
        break;

      default:
        fprintf (stderr, "ERROR:pore_multi: default state case activated (%i)\n", c_p->State);
        errors = raise (SIGTERM);
        exit (2);
        break;
      }                         /*end of pore state switch */

      /* call routine to put the info in the output list */
      trad_pore (bp, c_p, trad_index);

#ifdef  T_PORE
      {                         /* trace one pore */
        FILE *dumbfile;
        int dumb;

        if (i == T_PORE && c_p->State != PORE_NONE) {
          dumb = 0;
          dumbfile = fopen ("T_Pore.csv", "a");
          fprintf (dumbfile, "%.10e,%i, %.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n", cell_temp,
                   c_p->State, c_p->Radius, c_p->Volume, c_p->EqRad, c_p->Pressure, c_p->Nmols, c_p->NmolsH, c_p->Limrad, cell_fs,
                   cell_sol, cell_si, cell_satSI, ss);
          fclose (dumbfile);
          dumb = 0;
          if (eq_rad > 3e-5) {
            dumb = 0;
          }

        }
      }
#endif

    }                           /* end of NOT_CASTING else */
  }                             /*end of loop thru pores */

  /* debug by continuously listing one particular pore */
  return (errors);

}                               /* end of sb_pore_multi */

/***********************************/
/* check if trad_index has changed */
/* and insert the data if needed   */
/***********************************/
void trad_pore (BB_struct * bp, PORE_str * c_p, int trad_index)
{
  CA_FLOAT **c_t_lists;

  /* the pore doesn't exist yet */
  if (c_p->State == PORE_NONE)
    return;
  if (trad_index < 0)
    return;

  c_t_lists = c_p->t_lists;
  if ((trad_index != c_p->trad_last) && (trad_index >= 0) && (trad_index < bp->pprops.P_ntrad)) {
    *(c_t_lists[T_RAD_LIST] + trad_index) = c_p->Radius;
    *(c_t_lists[T_VOL_LIST] + trad_index) = c_p->Volume;
    *(c_t_lists[T_EQR_LIST] + trad_index) = c_p->EqRad;
    *(c_t_lists[T_PRES_LIST] + trad_index) = c_p->Pressure;
    *(c_t_lists[T_NMOL_LIST] + trad_index) = c_p->Nmols;
    *(c_t_lists[T_LIM_LIST] + trad_index) = c_p->Limrad;
    *(c_t_lists[T_NCELL_LIST] + trad_index) = c_p->ncells;
    *(c_t_lists[T_FRACSOL_LIST] + trad_index) = c_p->FracSol;
    c_p->trad_last = trad_index;
  }
}

/* insert the RCS id string in the program */
char const *rcs_id_pore_c ()
{
  static char const rcsid[] = "$Id: pore_multi.c 1373 2008-08-27 20:51:52Z  $";

  return (rcsid);
}

/*
*/
