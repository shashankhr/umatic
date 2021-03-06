
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
/* sb_temp_calc.c:						*/
/* Calculate the new temperature for a sub-block. This routine  */
/* is only to be used in solo mode, as the Temperatrue is       */
/* supplied when running with a FEM model, or in postprocessing */
/* mode.                                                        */
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Aug 22, 1998                                                  */
/****************************************************************/
/* 	MODIFIED by:						*/
/*  PDL: Aug 22, 1998						*/
/****************************************************************/
/****** To Do List **********************************************/
/*General:							*/
/* ONLY DOES CONT. Drop in Temp!!!                              */
/* 1) Convert to true enthalpy calc...                          */
/****************************************************************/
/*RCS Id:$Id: sb_temp_calc.c 1378 2008-09-10 20:21:14Z  $*/
#include <stdio.h>
#include <signal.h>
#include "machine.h"
#include "blocks.h"
#include "interp.h"

/* functions used from fidap_interp_calc.c */
extern CA_FLOAT fidap_interp_calc (FGrid_str * fg, CA_FLOAT x, CA_FLOAT y, CA_FLOAT time);

/* functions used from thermo_trace_calc.c */
extern CA_FLOAT thermo_trace_calc (TC_str * tc, CA_FLOAT tt);

/* functiosn used from bigblock.c */
extern void set_global_temperature (BB_struct * bp);
void sb_update_cell_temp (BB_struct * bp, int sbnum);

/***************************************************************/
/* Calculate the temperature if the finite grid method is used */
/***************************************************************/

CA_FLOAT fg_temp_calc (BB_struct * bp, int sbnum, int x, int y, int z)
{
  /* z is ignored but all temp calc functions must have the same args */
  int index_umat_2d;
  CA_FLOAT T_result, diff;

  index_umat_2d = x + bp->nc[0] * y;
  diff = bp->sb[sbnum]->c_fg_temp_next[index_umat_2d] - bp->sb[sbnum]->c_fg_temp[index_umat_2d];
  T_result = bp->sb[sbnum]->c_fg_temp[index_umat_2d] + bp->fg->wfact * diff;
  return (T_result);
}

CA_FLOAT const_temp_calc (BB_struct * bp, int sbnum, int x, int y, int z)
{
  return (bp->sb[sbnum]->Tvals.Tavg);
}

CA_FLOAT external_temp_calc (BB_struct * bp, int sbnum, int x, int y, int z)
{
/**  \todo  improve external temp array --  external - coupled */
/* put the array at the subblock level or */
/* quick fix -- link the single subblock array to point to */
/* the same array as bp->current_cell_temp */

  int index_ca;
  CA_FLOAT this_temp;

  index_ca = x + y * bp->nc[0] + z * bp->nc[0] * bp->nc[1];
  return (bp->current_cell_temp[index_ca]);

}                               /* end external_temp_calc */

CA_FLOAT dir_temp_calc (BB_struct * bp, int sbnum, int x, int y, int z) {
  CA_FLOAT T, rad, rn, axis;    /* tmp CA_FLOAT var. */
  CA_FLOAT v, vt, tgrad, slope, iso_coef1, iso_coef2;   /* tmp CA_FLOAT var. */
  CA_FLOAT velo_coef, grad_coef;        /*varying G. V */

  switch(bp->ctrl->swap_axes){
     case (0) :
    /* vertical growth */
    rad = bp->sb[sbnum]->orig_sb[0] + x * bp->size_c[0];
    axis = bp->sb[sbnum]->orig_sb[1] + y * bp->size_c[1];
     break;
     case(1):
/* horizontal growth */
    /* Swap x and y. by Wei WANG 07-08-02 */
    axis = bp->sb[sbnum]->orig_sb[0] + x * bp->size_c[0];
    rad = bp->sb[sbnum]->orig_sb[1] + y * bp->size_c[1];
     break;
     case(2):
        axis = bp->sb[sbnum]->orig_sb[0] + z * bp->size_c[0];
        rad = bp->sb[sbnum]->orig_sb[1] + x * bp->size_c[1];
     break;
     default:
    /* vertical growth */
    rad = bp->sb[sbnum]->orig_sb[0] + x * bp->size_c[0];
    axis = bp->sb[sbnum]->orig_sb[1] + y * bp->size_c[1];
     break;
     }
     
     #ifdef JUNK
  if (bp->ctrl->swap_axes  == ) {
/* horizontal growth */
    /* Swap x and y. by Wei WANG 07-08-02 */
    axis = bp->sb[sbnum]->orig_sb[0] + x * bp->size_c[0];
    rad = bp->sb[sbnum]->orig_sb[1] + z * bp->size_c[1];
  } else {
    /* vertical growth */
    rad = bp->sb[sbnum]->orig_sb[0] + x * bp->size_c[0];
    axis = bp->sb[sbnum]->orig_sb[1] + y * bp->size_c[1];
  }
  #endif /*JUNK*/

  rn = (CA_FLOAT) y / (CA_FLOAT) (bp->nc[1]);
/*****************************************************/
/* Calculate the temperature using analytic approx.  */
/*  uses a set gradient and velocity                 */
/* which is probably most useful for comparing with  */
/* directional solidification experiments.           */
/*****************************************************/
  velo_coef = bp->velo_coef;
  grad_coef = bp->grad_coef;
  tgrad = bp->gradient + grad_coef * bp->sim_time;
  v = bp->velocity + velo_coef * bp->sim_time;
  vt = v * bp->sim_time;
  iso_coef1 = bp->iso_coef1;
  iso_coef2 = bp->iso_coef2 / bp->size_c[0] / bp->nc[0];
  if (bp->grad_slope != 0) {
    slope = bp->grad_slope;     /* thermal gradient deviation slope */
    T = bp->Tinit + (tgrad) * (axis - vt) - ((rn < 0.5) ? 0 : (slope * (rn - 0.5)));
  } else {
    /* use a tilted gradient */
    if (bp->ctrl->gradtilt == 1 ){
	    double aprime,rprime;
	    double ct,st;
            rprime = (bp->grad_ct) * rad +bp->grad_st * axis;
            aprime = -(bp->grad_st) * rad +bp->grad_ct * axis;
	    axis=aprime;
	    rad=rprime;
	    }
    T = bp->Tinit + (tgrad) * (axis - vt);
    T -= tgrad * (iso_coef1 * rad + iso_coef2 * rad * rad);
  }
  T = MIN (T, bp->Tinit + CAP_TEMP_OFFSET);
  return (T);
}                               /* end of dir_temp_calc */

#ifdef JUNK
CA_FLOAT olddir_temp_calc (BB_struct * bp, int sbnum, int x, int y, int z)
{
  CA_FLOAT T, rad, rn, axis;    /* tmp CA_FLOAT var. */
  CA_FLOAT v, vt, tgrad, slope, iso_coef1, iso_coef2;   /* tmp CA_FLOAT var. */
  CA_FLOAT velo_coef, grad_coef;        /*varying G. V */

  if (bp->ctrl->swap_xy == 1) {
/* horizontal growth */
    axis = bp->sb[sbnum]->orig_sb[0] + x * bp->size_c[0];
    rad = bp->sb[sbnum]->orig_sb[1] + y * bp->size_c[1];
  } else {
    /* vertical growth */
    /* Swap x and y. by Wei WANG 07-08-02 */
    rad = bp->sb[sbnum]->orig_sb[0] + x * bp->size_c[0];
    axis = bp->sb[sbnum]->orig_sb[1] + y * bp->size_c[1];
  }

  rn = (CA_FLOAT) y / (CA_FLOAT) (bp->nc[1]);
/*****************************************************/
/* Calculate the temperature using analytic approx.  */
/*  uses a set gradient and velocity                 */
/* which is probably most useful for comparing with  */
/* directional solidification experiments.           */
/*****************************************************/
  velo_coef = bp->velo_coef;
  grad_coef = bp->grad_coef;
  tgrad = bp->gradient + grad_coef * bp->sim_time;
  v = bp->velocity + velo_coef * bp->sim_time;
  vt = v * bp->sim_time;
  iso_coef1 = bp->iso_coef1;
  iso_coef2 = bp->iso_coef2 / bp->size_c[0] / bp->nc[0];
  if (bp->grad_slope != 0) {
    slope = bp->grad_slope;     /* thermal gradient deviation slope */
    T = bp->Tinit + (tgrad) * (axis - vt) - ((rn < 0.5) ? 0 : (slope * (rn - 0.5)));
  } else {
    T = bp->Tinit + (tgrad) * (axis - vt);
    T -= tgrad * (iso_coef1 * rad + iso_coef2 * rad * rad);
                                                      /**dong*/
  }
  T = MIN (T, bp->Tinit + CAP_TEMP_OFFSET);
  return (T);
}                               /* end of dir_temp_calc */
#endif /*JUNK*/

/****************************************************************/
/****************************************************************/
/* Update the cell temperature for the whole array              */
/****************************************************************/
/****************************************************************/
void sb_update_cell_temp (BB_struct * bp, int sbnum)
{
  /* Update the entire cell temperature array */
  int i, j, k;
  CA_FLOAT x, y, z;
  CA_FLOAT *c_temp_p;           /* temperature pointer */
  int *gr_p;
  int n_casting = 0;
  CA_FLOAT t_sum = 0, tmin = 10000.0, tmax = 0.0;

  if (bp->sb[sbnum]->open != TRUE)
    return;
  /**  \todo  Figure out a way to test not yet open subblocks  -- multiblock */
  /* without needing the T array defined already */

  c_temp_p = bp->sb[sbnum]->c_temp;
  gr_p = bp->sb[sbnum]->gr;

  for (k = 0; k < bp->nc[2]; k++) {
    for (j = 0; j < bp->nc[1]; j++) {
      for (i = 0; i < bp->nc[0]; i++) {
        if (*gr_p != NOT_CASTING) {
          *c_temp_p = (bp->cell_temp_func) (bp, sbnum, i, j, k);
          /* update the block temperature stats */
          t_sum += *c_temp_p;
          tmin = MIN (tmin, *c_temp_p);
          tmax = MAX (tmax, *c_temp_p);
          n_casting++;
        }

        gr_p++;
        c_temp_p++;
      }                         /* end of i loop */
    }                           /* end of j loop */
  }                             /* end of k loop */

  /* copy the stat vals into the subblock */
  bp->sb[sbnum]->Tvals.Tavg = t_sum / (n_casting);
  bp->sb[sbnum]->Tvals.Tmin = tmin;
  bp->sb[sbnum]->Tvals.Tmax = tmax;
  bp->sb[sbnum]->Tvals.TminReached = MIN (tmin, bp->sb[sbnum]->Tvals.TminReached);
}                               /* end of update cell temp routine */

/****************************************************************/
/****************************************************************/
/* sb_temp_calc: Calculate the new temperature for a sub-block  */
/****************************************************************/
/****************************************************************/
int sb_temp_calc (BB_struct * bp, int sbnum)
{
  CA_FLOAT delT;                /* tmp CA_FLOAT var. */
  T_struct *Tp;
  Mat_str *mp;

  Tp = &(bp->sb[sbnum]->Tvals);
  mp = &(bp->mprops);

/* check for temperature lookup table */
/* used by Ali's version -- not sure what it is for!*/
/** \todo  ask Ali about this temperature lookup control - maybe obsolete - external */
  if (bp->ctrl->temp_lookup == FALSE) {
    if (!(bp->ctrl->use_cell_temp)) {
      /* Constant Cooling only */
      Tp->TminReached = MIN (Tp->TminReached, Tp->Tavg);        /* update minimum T reached first... */
      Tp->Tavg_old = Tp->Tavg;
      if (bp->ctrl->t_input) {
        /* Using thermocouple trace input */
        Tp->Tavg = thermo_trace_calc (&(bp->tc), bp->sim_time);

      } else if (bp->ctrl->coolrate == TRUE) {
        /* constant cooling rate */
        Tp->Tavg += (bp->ctrl->delT) * (bp->delt);
        Tp->del_fs = 0;

      } else {

        /* Using heatflow calculation */
        /* calc the change in temperature */
        /* delT = (QSV*delt + latentH*del_fs)/rhocp; */
        /*delT = (mp->QSV * bp->delt + (mp->latentH * Tp->del_fs*bp->vol_sb))/mp->rhocp; */
        delT = (mp->QSV * bp->delt + (mp->latentH * Tp->del_fs)) / mp->rhocp;
        Tp->Tavg += delT;
        Tp->del_fs = 0;
      }

    }

    /* update the cell temperature array, for both isothermal and cell temperature */
    sb_update_cell_temp (bp, sbnum);
  } else {
    /* Ali Chirazi global temperature option */
    set_global_temperature (bp);
    Tp->Tavg = bp->this_temp;
  }

  /*Failsafe finish if temperature goes below zero or FROZEN */
  if (Tp->Tavg < FROZEN) {
    fprintf (stderr, "ERROR: sb_temp_calc: Temperature has frozen! %.5g  Finishing ...\n", Tp->Tavg);
    return (raise (SIGTERM));
  }

  return (0);
}                               /* end of sb_temp_calc subroutine */

/*EXTERNAL routines here!!! */
/*EXTERNAL routines here!!! */
/*EXTERNAL routines here!!! */
/*EXTERNAL routines here!!! */
/*EXTERNAL routines here!!! */
/******************************************************************/
/* special case of Procast mode on calculation for temperature*****/
/*****************************************************************/
void cell_temp_calc_extern (BB_struct * bp, int index_ca)
{
  float t_elapsed, temp_change;

  if (bp->cell_element_array[index_ca] >= 0) {

    if (bp->first_micro == 0) {
      bp->current_cell_temp[index_ca] = bp->cell_temp_extern[index_ca];
    } else {
      t_elapsed = bp->delt * bp->micro_step;
      temp_change = bp->cell_temp_change_extern[index_ca] * t_elapsed;
      bp->current_cell_temp[index_ca] = bp->cell_temp_extern[index_ca] + temp_change;
      /*
         bp->current_cell_temp[index_ca]+=(bp->cell_temp_change_extern[index_ca]*bp->delt);
       */
    }
  } else {
    bp->current_cell_temp[index_ca] = 0.0;
  }

  return;
}

#ifdef AAA
void cell_temp_calc_extern (BB_struct * bp, int index_ca)
{
  double t_elapsed, temp_change;

  if (bp->cell_element_array[index_ca] >= 0) {
    if (bp->first_micro == 0) {
      bp->current_cell_temp[index_ca] = bp->cell_temp_extern[index_ca];
    } else {
      t_elapsed = (double) ((double) bp->delt * (double) bp->micro_step);
      temp_change = ((double) (bp->cell_temp_change_extern[index_ca]) * t_elapsed);
      bp->current_cell_temp[index_ca] = (float) ((double) (bp->cell_temp_extern[index_ca])
                                                 + temp_change);

      /*
         bp->current_cell_temp[index_ca]+=(bp->cell_temp_change_extern[index_ca]*bp->delt);
       */
    }
  } else {
    bp->current_cell_temp[index_ca] = 0.0;
  }

  return;
}
#endif
/******************************************************************/
/* special case of Procast mode on for min temp calculation */
/******************************************************************/
float min_temp_calc_extern (BB_struct * bp)
{
/*********from user_rtns.c**************/
#ifdef EXTERNAL_CA
  extern void shell (int n1, float *a);
#endif

 /****local variables*********/
  int k;
  float min_temp;
  float *local_temp;
  int local_counter;

    /*******allocate local temp array**************/
  local_temp = (float *) calloc (bp->total_cell_number + 1, sizeof (float));

    /***copy the temperature array into a local array for sort purposes***/
  for (k = 0, local_counter = 0; k < bp->total_cell_number; k++) {
    if (bp->cell_element_array[k] >= 0) {
      local_counter++;
      local_temp[local_counter] = bp->current_cell_temp[k];
    }
  }
    /*************************************************/

    /***********call shell function to sort the array*********/
#ifdef EXTERNAL_CA
  shell (local_counter, local_temp);
#endif

  min_temp = local_temp[1];
    /**********free the local array once the min temp is assigned******/
  free (local_temp);

  return (min_temp);
}

/****************************************************************/
/******************************************************************/
/* special case of Procast mode on calculation for pressure*****/
/*****************************************************************/
void find_cell_pressure (BB_struct * bp, int index_ca)
{
  if (bp->cell_element_array[index_ca] >= 0) {
    if (bp->first_micro == 0) {
      bp->current_cell_pres[index_ca] = bp->cell_pres_extern[index_ca];
    } else {
      bp->current_cell_pres[index_ca] += (bp->cell_pres_grad_extern[index_ca] * bp->delt);
    }
  } else {
    bp->current_cell_pres[index_ca] = 0.0;
  }

  return;
}

/*finite element (fgrid) routines here ! */
/*finite element (fgrid) routines here ! */
/*finite element (fgrid) routines here ! */
/*finite element (fgrid) routines here ! */
/*finite element (fgrid) routines here ! */

/*******************************************/
/* Set up the first temperature interp.    */
/*******************************************/
/* using the first two files */
/* TWO DIMENSIONAL temperature distribution only */
void sb_temp_setup (BB_struct * bp, int sbnum)
{
  int i, j, k, index_umat_2d;

  for (k = 0, index_umat_2d = 0; (k < bp->nc[2] && index_umat_2d < bp->nc[0] * bp->nc[1] * bp->nc[2]); k++) {
    for (j = 0; j < bp->nc[1]; j++) {
      for (i = 0; i < bp->nc[0]; i++) {
        /* get the temperature at this (x,y) location and store */
        bp->sb[sbnum]->c_fg_temp[index_umat_2d] = trans_interp_calc (bp->fg, bp->sb[sbnum]->nnd, bp, sbnum, i, j);
        bp->sb[sbnum]->c_fg_temp_next[index_umat_2d] = trans_interp_calc (bp->fg_next, bp->sb[sbnum]->nnd_next, bp, sbnum, i, j);

        index_umat_2d++;
      }
    }
  }

}                               /* end of sb_temp_setup */

/* update the start and finish temperature arrays */
/* after a new finite element (fgrid) file is read in */
void fg_temp_upd (BB_struct * bp, int sbnum)
{
  int i, j, k, index_ca;
  CA_FLOAT *holder;
  NODENB_str *n_holder;

  /* exchange the pointers for the temperatrue arrays */
  holder = bp->sb[sbnum]->c_fg_temp;
  bp->sb[sbnum]->c_fg_temp = bp->sb[sbnum]->c_fg_temp_next;
  bp->sb[sbnum]->c_fg_temp_next = holder;

  /* exchange the pointers for the interploation arrays */
  n_holder = bp->sb[sbnum]->nnd;
  bp->sb[sbnum]->nnd = bp->sb[sbnum]->nnd_next;
  bp->sb[sbnum]->nnd_next = n_holder;

   /*******************************************/
  /* Calculate the weighting factors in z    */
  /* only when a new fg has been read in      */
   /*******************************************/
  wfact_z_calc (bp->fg_next, bp->sb[sbnum]->nnd_next, bp, sbnum);
  /* update the temperature for the next step */
  for (k = 0, index_ca = 0; (k < bp->nc[2] && index_ca < bp->nc[0] * bp->nc[1] * bp->nc[2]); k++) {
    for (j = 0; j < bp->nc[1]; j++) {
      for (i = 0; i < bp->nc[0]; i++) {
        bp->sb[sbnum]->c_fg_temp_next[index_ca] = trans_interp_calc (bp->fg_next, bp->sb[sbnum]->nnd_next, bp, sbnum, i, j);

        index_ca++;
      }                         /*i */
    }                           /*j */
  }                             /*k */
}                               /* end of fg_temp_upd */

/********************************************************/
/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
/********************************************************/
char const *rcs_id_sb_temp_calc_c ()
{
  static char const rcsid[] = "$Id: sb_temp_calc.c 1378 2008-09-10 20:21:14Z  $";

  return (rcsid);
}

/* end of rcs_id_sb_temp_calc_c subroutine */
/*
*/
