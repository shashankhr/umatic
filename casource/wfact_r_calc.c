
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
/* T_interp_calc.c:                                             */
/* Subroutine to interpolate the temperature at a point         */
/* in the CA code using values calculated in a                  */
/* transient state solution of the heat, mass and               */
/* momentum transfer in VAR ingots.                             */
/****************************************************************/
/****************************************************************/
/* Written by X. Xu, P.D. Lee & R.C. Atwood, Imperial College   */
/* Feb. 18, 2000                                                */
/****************************************************************/
/*      MODIFIED by:                                            */
/****************************************************************/
/****** To Do List **********************************************/
/*General:                                                      */
/* 1)                                                           */
/****************************************************************/
#include <stdio.h>
#include <math.h>
/* include header files requred by subroutines */
#include "machine.h"
#include "blocks.h"
#include "fidap.h"              /* included for def. of FGrid_struct */
#include "nearnode.h"
int search_r (FGrid_str * fg, CA_FLOAT r);

/****************************************************************/
/****************************************************************/
/* wfact_r_calc.c:                                             */
/* Subroutine to calculation the weighting factor in r direction */
/****************************************************************/
/* Input Variables:                                             */
/*   *fg:		ptr to the FGrid_str structure that     */
/*       		holds all FIDAP geometry and T's data.  */
/*   *bp:+		hold location of CA cell                */
/*      		                                        */
/* Output Variables:    weight factor in r direction            */
/*      		                                        */
/* Returned Value:      none                                    */
/****************************************************************/
/****************************************************************/

void wfact_r_calc (FGrid_str * fg, NODENB_str * node_ptr, BB_struct * bp, int sbnum)
{
  CA_FLOAT rc, rm;
  int j, jj, jc, jmax;
  CA_FLOAT sub_off, cell_offset, size_cell, *fg_r_pt;

  /* set pointers and variables for local usage */

  fg_r_pt = fg->r;
  sub_off = bp->sb[sbnum]->orig_sb[1];  /* lower left corner of the subblock */
  size_cell = bp->size_c[1];    /* cell size in z dir, it equals to that in r dir. */
  cell_offset = sub_off - 0.5 * size_cell;      /* center of the previous cell. */

  jmax = bp->nc[1];
  rc = cell_offset;             /* Initial value  */

  /* calculate positions of r arrays */
  for (j = 0; j < jmax; j++) {
    rc += size_cell;
    /* calculate r arrays of cells in a CFD domain */
    rm = rc;

         /*******************************************************/
    /*  check if the rc of cell is out of the CFD domain,  */
    /*  and given an error message if it is.               */
         /*******************************************************/
    if (rm < 0.0 || rm > fg_r_pt[fg->nr - 1]) {
      fprintf (stderr, "WARNING r value %i  out of domain\n",j);
      fprintf (stderr, "WARNING value %f  does not lie between 0 and %f\n", rm , fg_r_pt[fg->nr - 1]);
    } else {

      jc = search_r (fg, rm);

      node_ptr->nl[j] = jc;     /* nearest node on the left of the cell. */
       /***************************************************************/
      /* factor that the cell is away from the left/right nodes,     */
      /* when it lies on the left node, wl[i] = 1.                   */
       /***************************************************************/

      node_ptr->wl[j] = (fg_r_pt[jc + 1] - rm) / (fg_r_pt[jc + 1] - fg_r_pt[jc]);
      node_ptr->wr[j] = 1 - node_ptr->wl[j];
    }
  }
  /* end of wfact_r_calc     */
}

/* subroutine of bisection search */

int search_r (FGrid_str * fg, CA_FLOAT r)
{
  int i, im, iu, il, JP;
  CA_FLOAT *fg_r_pt;

  fg_r_pt = fg->r;
  iu = fg->nr - 1;
  il = 0;
  while (iu - il > 1) {
    im = (iu + il) / 2;
    if ((r > fg_r_pt[il]) && (r > fg_r_pt[im])) {
      il = im;
    } else {
      iu = im;
    }
  }
  JP = il;
  return (JP);
}

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_wfact_r_calc_c ()
{
  static char const rcsid[] = "$Id: wfact_r_calc.c 1402 2008-11-20 15:36:41Z  $";

  return (rcsid);
}

/* end of rcs_id_fidap_interp_calc_c subroutine */

/*RCS Id:$Id: wfact_r_calc.c 1402 2008-11-20 15:36:41Z  $*/
/*
*/
