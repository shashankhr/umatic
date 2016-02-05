
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
#include <string.h>
#include "blocks.h"
#include "machine.h"
void sb_mould_src (BB_struct * bp, Solute_props * sp, int sbnum, CA_FLOAT * osol)
{
  int c;
  CA_FLOAT *concp, *fsp, *tempp, add;
  int i, j, k;
  static int flag = 0;
  static int orow = 0, row = 0;
  static int olayer = 0, layer = 0;
  int offset;

  bp->sb[sbnum]->addsol = 0;

  if (flag == 0) {
    /* precalculate the offset shift for padded array */
    /**  \todo  store in block or cube structure!  general*/
    orow = (bp->nc[0] + 2);
    row = (bp->nc[0]);
    olayer = ((bp->nc[0] + 2) * (bp->nc[1] + 2));
    layer = ((bp->nc[0]) * (bp->nc[1]));
    flag = 1;
  }

  /* locate the surface cell from the surface array */
  for (c = 0; c < bp->sb[sbnum]->surface.ns_cell; c++) {
    i = (*(bp->sb[sbnum]->surface.surf_xyz + c))[0];
    j = (*(bp->sb[sbnum]->surface.surf_xyz + c))[1];
    k = (*(bp->sb[sbnum]->surface.surf_xyz + c))[2];
    /* get the index */
    offset = *(bp->sb[sbnum]->surface.c_surfp + c);

    fsp = bp->sb[sbnum]->c_fs + offset;
    tempp = bp->sb[sbnum]->c_temp + offset;

    /* osol is the PADDED array     */
    /* so offset indices are needed */
    concp = osol + (i + 1)
      + (j + 1) * orow + (k + 1) * olayer;

    /**  \todo  improve mould_src_func to account for fraction solid  -- mould */
    /* for now , it stops when (nearly) fully solid */
    if (*fsp >= 0.999) {
      add = 0;
    } else {
      /* the actual function is selected by an option */
      add = (*(sp->mould_src_func)) (bp, sp, *tempp, *concp, i, j, k);
    }
    *concp += add;
    bp->sb[sbnum]->Svals[sp->my_type].addsol += add;
  }

  bp->sb[sbnum]->Svals[sp->my_type].t_addsol += bp->sb[sbnum]->addsol;
  return;
}                               /* end of sb_mould_src */

CA_FLOAT mould_src (BB_struct * bp, Solute_props * sp, CA_FLOAT cell_temp, CA_FLOAT conc)
{
  /* mould flux */
  CA_FLOAT deltaC;

  deltaC = sp->mould_source_value * bp->size_c[0];
  return (conc + deltaC);
}

/********************************/
/* mould source functions *******/
/* to simplify the programming  */
/* all functions must have the  */
/* same arguments -- even though*/
/* they may get ignored for many*/
/* cases.                       */
/********************************/

/* constant mould concentration */
CA_FLOAT const_mould_src (BB_struct * bp, Solute_props * sp, CA_FLOAT cell_temp, CA_FLOAT conc, int i, int j, int k)
{
  /* constant mould concentration */
  CA_FLOAT sourceval;

  sourceval = (sp->mould_source_value - conc);
  return (sourceval);
}                               /* end of constant mould concentration */

CA_FLOAT flux_mould_src (BB_struct * bp, Solute_props * sp, CA_FLOAT cell_temp, CA_FLOAT conc, int i, int j, int k)
{
  /* constant mould concentration */
  CA_FLOAT sourceval;

  sourceval = (sp->mould_source_value * bp->delt / bp->size_c[0]);
  return (sourceval);
}                               /* end of constant mould concentration */

CA_FLOAT none_mould_src (BB_struct * bp, Solute_props * sp, CA_FLOAT cell_temp, CA_FLOAT conc, int i, int j, int k)
{
  /* no mould concentration */
  CA_FLOAT sourceval;

  sourceval = 0;
  return (sourceval);
}                               /* end of zero mould concentration */

/* mould concentration based on constant interface */
CA_FLOAT diff_mould_src (BB_struct * bp, Solute_props * sp, CA_FLOAT cell_temp, CA_FLOAT conc, int i, int j, int k)
{
  CA_FLOAT sourceval;
  CA_FLOAT rs;

  rs = sp->Dliq * (bp->delt / (bp->size_c[0] * bp->size_c[0]));
  sourceval = (sp->mould_source_value - conc) * rs;
  return (sourceval);
}                               /* end of constant mould concentration */

/* mould concentration based on piecewise linear approximation */

CA_FLOAT plin_mould_src (BB_struct * bp, Solute_props * sp, CA_FLOAT cell_temp, CA_FLOAT conc, int i, int j, int k)
{
  CA_FLOAT sourceval, econ;
  CA_FLOAT rs;

  /* hard coded TRIAL version */
  /**  \todo  create an input section for the  source/temperature relationshop -- mould */
  static const CA_FLOAT tpoints[5] = { 0, 800, 900, 1658, 3000 };
  static const CA_FLOAT cpoints[5] = { 0, .01, .05, 1, 5 };
  static const int npoints = 5;

  int nn;

  for (nn = 0; nn < npoints - 1; nn++) {
    if (cell_temp > tpoints[nn] && cell_temp <= tpoints[nn + 1]) {
      econ = cpoints[nn] + ((cpoints[nn + 1] - cpoints[nn]) / (tpoints[nn + 1] - tpoints[nn])) * (cell_temp - tpoints[nn]);
    }
  }
  if (cell_temp <= tpoints[0]) {
    econ = cpoints[0];
  }
  if (cell_temp > tpoints[npoints - 1]) {
    econ = cpoints[npoints - 1];
  }

  /* scale by the chosen factor from ctrl file */
  econ *= sp->mould_source_value;

  /**  \todo  get average diff.coeff. for source cells -- general - mould */
  rs = sp->Dliq * (bp->delt / (bp->size_c[0] * bp->size_c[0]));
  sourceval = (econ - conc) * rs;
/*
   if (bp->step%bp->ctrl->scr_dmp_freq ==0 ){
      fprintf(stderr,"plin: %i %i %.5g %.5g %.5g\n",i,j,cell_temp,econ,sourceval);
   }
*/

  return (sourceval);
}                               /* end of piecewise linear mould concentration */

/* perturbed mould concentration */
CA_FLOAT perturb_mould_src (BB_struct * bp, Solute_props * sp, CA_FLOAT cell_temp, CA_FLOAT conc, int i, int j, int k)
{
  /* perturbed mould concentration */
  CA_FLOAT add;
  CA_FLOAT sourceval;

  add = SIN (bp->ctrl->mould_source_freq * j) * bp->ctrl->mould_src_pert;
  /* the source cannot become a sink in this option! */
  sourceval = MAX (0, sp->mould_source_value + add);
  sourceval *= bp->delt;
  return (sourceval);
}                               /* end of perturbed mould concentration */
