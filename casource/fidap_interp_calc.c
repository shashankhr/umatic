
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
/* fidap_interp_calc.c:                                         */
/* Subroutine to interpolate the temperature at a point         */
/* in the CA code using values calculated in FIDAP in a         */
/* quasi-steady state solution of the heat, mass and            */
/* momentum transfer in VAR ingots.                             */
/****************************************************************/
/****************************************************************/
/* Written by X. Xu, P.D. Lee & R.C. Atwood, Imperial College   */
/* Nov 10, 1998                                                 */
/****************************************************************/
/*      MODIFIED by:                                            */
/*  PDL: Nov 10, 1998                                           */
/****************************************************************/
/****** To Do List **********************************************/
/*General:                                                      */
/* 1)                                                           */
/****************************************************************/
/*RCS Id:$Id: fidap_interp_calc.c 1341 2008-07-23 15:23:30Z  $*/
#include <stdio.h>
#include <math.h>

/* include header files requred by subroutines */
#include "machine.h"
#include "fidap.h"              /* included for def. of FGrid_struct */

/****************************************************************/
/****************************************************************/
/* fidap_interp_calc:                                           */
/* Subroutine to interpolate temp. of the cell                  */
/****************************************************************/
/* Input Variables:                                             */
/*      		*Note: currently FIDAP geom. is trans.  */
/*      		to [mm] from [m] since [mm] used in CA  */
/*   *fg:		ptr to the FGrid_str structure that     */
/*       		holds all FIDAP geometry and T's data.  */
/*   xi:		x location in the CA model (equiv. to   */
/*       		radial dist. in FIDAP.           [mm]   */
/*   yi:		y location in the CA model (equiv. to   */
/*       		height or z in FIDAP.            [mm]   */
/*   time:		time in the CA model             [s]    */
/*      		                                        */
/* Output Variables:    NONE                                    */
/*      		                                        */
/* Returned Value:                                              */
/*   T: 		Temperature for given x, y, time in     */
/*       		CA model space interpolated from FIDAP  */
/*       		results.                                */
/****************************************************************/
/****************************************************************/
CA_FLOAT fidap_interp_calc (FGrid_str * fg, CA_FLOAT xi, CA_FLOAT yi, CA_FLOAT time)
{
  CA_FLOAT r, z, h, rmax, T;
  CA_FLOAT tp1, tp2, tp3, tp4, tmpx1, tmpx2;
  int i, ib, ie, jb, je;
  CA_FLOAT mag, t_s[5], dt[5], dtstep;
  int j, nr1, nr2;

  h = fg->z[fg->nz - 1];
  rmax = fg->r[fg->nr - 1];
  /* fprintf(stderr," ***** %f %f ***** \n",h,rmax); */

/************************************************/
/* Check that the x value is in range           */
/************************************************/
  if ((xi > fg->r[fg->nr - 1]) || (xi < fg->r[0])) {
    return (fg->Tmin);
  }

  r = xi;

/* check radius range */
  for (i = 0; i < (fg->nr - 1); i++) {
    if (r > fg->r[i] && r <= fg->r[i + 1]) {
      jb = i;
      je = i + 1;
    }
    if (r == fg->r[i]) {
      jb = i;
      je = i;
    }
  }
  /* changed on nov.22,99 to see the growth direction with T */
  /* z = -yi + fg->v*time*(1+(fg->r[fg->nr-1]- xi)/fg->r[fg->nr-1])/10000.0; */
  z = -yi + fg->v * time;       /* PREVIOUS ONE */
/* first give a general situation, when treeing_on, change it */

/* Adjust z if Treerings on */
  if (fg->tring_on) {
    /*     if (!(tsta = (float *) calloc(fg->ntring, sizeof(float)))) {
       fprintf(stderr, "Error: t_s allay malloc failed \n"); }
       if (!(dtpe = (float *) calloc(fg->ntring, sizeof(float)))) {
       fprintf(stderr, "Error: dt_p allay malloc failed \n");}  */
    dtstep = time - fg->oldtime;
    /* fprintf(stderr," ** %d ***\n ",fg->ntring);  */
/************************************************/
/* Check that the x value is in range           */
/************************************************/
    for (i = 0; i < fg->ntring; i++) {
      t_s[i] = fg->tring[i].v_f[1];
      dt[i] = fg->tring[i].v_f[2];
    }                           /* end of read into time and duration for a pertubation */

#ifdef PDL_TREE

    if ((time >= t_s) && (time <= (dt + t_s))) {
      mag = fg->tring[i].v_f[0];
      nr1 = fg->tring[i].v_int[0];
      nr2 = fg->tring[i].v_int[1];
      /*  fprintf(stderr," ** %d %d ***\n ",fg->tring[i].v_int[0],fg->tring[i].v_int[1]); */

      for (j = nr1; j < nr2; j++) {
        /*            fprintf(stderr," **** %d %d ***\n ",nr1, nr2);  */
        fg->zoffset[j] += fg->v * mag * dtstep;
        /*      fprintf(stderr," **** zoff in each step %f ***\n ",fg->zoffset[j]); */
      }
    }                           /* end of pertubation time */
#endif /* end of def PDL */

    for (i = 0; i < fg->ntring; i++) {
      if (time >= t_s[i] && time < t_s[i] + dt[i]) {
        mag = fg->tring[i].v_f[0];
        z = -yi + fg->v * time * (1.0 + mag);
      }
    }                           /* end of ntring */
    fg->oldtime = time;
  }

  /* end of tring_on */
  /* If the z value is above the pool, set = Tmax */
  if (z > h) {
    T = fg->Tmin;
    return (T);

/* If the z value is below the pool, set = Tmax */
  } else if (z < 0.0) {
    T = fg->Tmax;
    return (T);

/* z in range of FIDAP calculation, so interpolate! */
  } else {

    /* check z range */
    for (i = 0; i < (fg->nz - 1); i++) {
      if (z > fg->z[i] && z <= fg->z[i + 1]) {
        ib = i;
        ie = i + 1;
      }
      if (z == fg->z[i]) {
        ib = i;
        ie = i;
      }
    }

    /* CA location identical to node location, return value at node */
    if (ib == ie && jb == je) {
      T = fg->Fidap_T[fg->nr * ib + jb];

      /* CA location identical to node location in i dir, interp j */
    } else if (ib == ie && jb != je) {
      tp1 = fg->Fidap_T[fg->nr * ib + jb];
      tp2 = fg->Fidap_T[fg->nr * ib + je];
      T = tp1 * (fg->r[je] - r) / (fg->r[je] - fg->r[jb])
        + tp2 * (r - fg->r[jb]) / (fg->r[je] - fg->r[jb]);

      /* CA location identical to node location in j dir, interp i */
    } else if (ib != ie && jb == je) {
      tp1 = fg->Fidap_T[fg->nr * ib + jb];
      tp2 = fg->Fidap_T[fg->nr * ie + jb];
      T = tp1 * (fg->z[ie] - z) / (fg->z[ie] - fg->z[ib])
        + tp2 * (z - fg->z[ib]) / (fg->z[ie] - fg->z[ib]);

      /* interp in both directions */
    } else {                    /*  if(ib!=ie && jb!=je) */
      tp1 = fg->Fidap_T[fg->nr * ib + jb];
      tp2 = fg->Fidap_T[fg->nr * ib + je];
      tp3 = fg->Fidap_T[fg->nr * ie + jb];
      tp4 = fg->Fidap_T[fg->nr * ie + je];
      tmpx1 = tp1 * (fg->r[je] - r) / (fg->r[je] - fg->r[jb])
        + tp2 * (r - fg->r[jb]) / (fg->r[je] - fg->r[jb]);
      tmpx2 = tp3 * (fg->r[je] - r) / (fg->r[je] - fg->r[jb])
        + tp4 * (r - fg->r[jb]) / (fg->r[je] - fg->r[jb]);
      T = tmpx1 * (fg->z[ie] - z) / (fg->z[ie] - fg->z[ib])
        + tmpx2 * (z - fg->z[ib]) / (fg->z[ie] - fg->z[ib]);
    }
  }

  return (T);

}                               /* End of subroutine: fidap_interp_calc */

/****************************************************************/
/****************************************************************/
/* fidap_calc_tr_offset:                                        */
/* Subroutine to calculate the z offset to simulate tree rings  */
/****************************************************************/
/* Input Variables:                                             */
/*   *fg:		ptr to the FGrid_str structure that     */
/*       		holds all FIDAP geometry and T's data.  */
/*   time:		time in the CA model             [s]    */
/*      		                                        */
/* Output Variables:    NONE                                    */
/*      		                                        */
/* Returned Value:                                              */
/*   zoffset (implicit):offset in z direction, rtn through str  */
/****************************************************************/
/****************************************************************/
void fidap_calc_tr_offset (FGrid_str * fg, CA_FLOAT time)
{
  CA_FLOAT mag, t_s, dt, dtstep;
  int i, j, nr1, nr2;

  dtstep = time - fg->oldtime;

/************************************************/
/* Check that the x value is in range           */
/************************************************/
  for (i = 0; i < fg->ntring; i++) {

    t_s = fg->tring[i].v_f[1];
    dt = fg->tring[i].v_f[2];
    if ((time >= t_s) && (time < (dt + t_s))) {
      mag = fg->tring[i].v_f[0];
      nr1 = fg->tring[i].v_int[0];
      nr2 = fg->tring[i].v_int[1];
      for (j = nr1; j < nr2; j++) {
        fg->zoffset[j] += fg->v * mag * dtstep;
      }

    }
  }

  fg->oldtime = time;

}                               /* End of subroutine: fidap_calc_tr_offset */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_fidap_interp_calc_c ()
{
  static char const rcsid[] = "$Id: fidap_interp_calc.c 1341 2008-07-23 15:23:30Z  $";

  return (rcsid);
}

/* end of rcs_id_fidap_interp_calc_c subroutine */

/*
*/
