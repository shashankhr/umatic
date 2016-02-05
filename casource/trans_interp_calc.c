
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
/* trans_interp_calc.c:                                         */
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

/****************************************************************/
/****************************************************************/
/* trans_interp_calc.c                                          */
/* Subroutine to interpolate the temperature of cell            */
/****************************************************************/
/* Input Variables:                                             */
/*   *fg:		ptr to the FGrid_str structure that     */
/*       		holds all FIDAP geometry and T's data.  */
/*   r:	       	        radial coordinate of a cell centre      */
/*   z:                 Height coordination of a cell centre    */
/*      		                                        */
/* Output Variables:    NONE                                    */
/*      		                                        */
/* Returned Value:                                              */
/*   temperature in the centre of cellz                         */
/****************************************************************/

CA_FLOAT trans_interp_calc (FGrid_str * fg, NODENB_str * node_ptr, BB_struct * bp, int sbnum, int x, int y)
{
  int nr1, nr2, nz1, nz2, nnr;
  CA_FLOAT T, *fg_T_pt;

  /* pointer of val. for local usage. */
  fg_T_pt = fg->Fidap_T;
  nnr = fg->nr;

  /* In Z direction, check the node down the cell is over/below the CFD domain */

  nz1 = node_ptr->nd[x];

  if (nz1 == -1) {
    T = fg->Tmax;
  } else if (nz1 == fg->nz - 1) {
    T = fg->Tmin;
  } else {

    nr1 = node_ptr->nl[y];      /* near left node. */
    nr2 = nr1 + 1;
    nz2 = nz1 + 1;

    /* using pre-calculated interpolation weight functions (linear) */
    /*  T = wd*wl*tp1 + wd*wr*tp2 + wu*wl*tp3 + wu*wr*tp4; */
    T = node_ptr->wd[x] * (node_ptr->wl[y] * fg_T_pt[nnr * nz1 + nr1] +
                           node_ptr->wr[y] * fg_T_pt[nnr * nz1 + nr2]) +
      node_ptr->wu[x] * (node_ptr->wl[y] * fg_T_pt[nnr * nz2 + nr1] + node_ptr->wr[y] * fg_T_pt[nnr * nz2 + nr2]);
  }
  return (T);
}

/* end of trans_interp_calc(FGrid_str *fg, BB_struct *bp, int sbnum, int x, int y). */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_trans_interp_calc_c ()
{
  static char const rcsid[] = "$Id: trans_interp_calc.c 1341 2008-07-23 15:23:30Z  $";

  return (rcsid);
}

/* end of rcs_id_fidap_interp_calc_c subroutine */

/*RCS Id:$Id: trans_interp_calc.c 1341 2008-07-23 15:23:30Z  $*/
/*

*/
