
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
#include <string.h>
#include "machine.h"
#include "blocks.h"

int fcopy_surf (CA_FLOAT * new, CA_FLOAT * old, Frame * cubeptr, int fcode, int type, int dir);
int fcopy_neigh (BB_struct * bp, CA_FLOAT * new, int ci, CA_FLOAT * block_values[], int type, int dir);
int this_sb = 0;
int neigh_sb=0;


int fcopy_matrix (int type, CA_FLOAT * big, CA_FLOAT * lit, BB_struct * bp, CA_FLOAT * block_values[], int ci)
{
   this_sb=ci;

  int j, k, errors = 0;
  CA_FLOAT *plit, *pbig;
  int *ins, *outs;
  CA_FLOAT **bv_p;
  Frame *cubeptr;

  bv_p = block_values;
  cubeptr = &(bp->cubeptr);
  plit = lit;
  pbig = big;
  ins = cubeptr->ins;
  outs = cubeptr->outs;

  pbig += cubeptr->flist[0][START];

  for (k = 0; k < ins[ZAXIS]; k++) {
    for (j = 0; j < ins[YAXIS]; j++) {
      memcpy (pbig, plit, ins[XAXIS] * sizeof (CA_FLOAT));
      plit += ins[XAXIS];
      pbig += outs[XAXIS];
    }
    pbig += (outs[XAXIS] * 2);
  }

  if (bp->ntsb == 1) {
    for (j = 0; j < 6; j++) {
      type = cubeptr->facectrl[j];
      errors += fcopy_surf (big, lit, cubeptr, j, type, FORWARD);
    }
  }

  else
    errors += fcopy_neigh (bp, big, ci, bv_p, type, FORWARD);
  if (errors >= 1001)
    exit (1001);
  return (errors);
}


int fcopy_neigh (BB_struct * bp, CA_FLOAT * big, int ci, CA_FLOAT * block_values[], int type, int dir) {

  int code, nouts, errors = 0;
  int i, j;
  CA_FLOAT *lit;
  int ni;
  Frame *cubeptr;

  cubeptr = &(bp->cubeptr);
  code = bp->sb[ci]->code;
  nouts = bp->sb[ci]->nouts;

  switch (nouts) {
  case CM_INS:

    for (j = 0; j < 6; j++) {

      ni = ci + cubeptr->neigh[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:fcopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      lit = (block_values[ni]);
      if (lit == NULL) {
        fcopy_surf (big, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        fcopy_surf (big, lit, cubeptr, j, INS, dir);
      }
    }
    break;

  case CM_FACE:
    j = code;

    ni = ci + cubeptr->face[j];
    if (ni > bp->ntsb) {
      fprintf (stderr, "ERROR:fcopy_neigh: outs.neighbour index too high! ni %i ci %i\n", ni, ci);
      exit (100);
    }
    if (type == WRAP)
      lit = (block_values[ni]);
    else
      lit = (block_values[ci]);
    if (lit == NULL) {
      fcopy_surf (big, block_values[ci], cubeptr, j, PAD, dir);
    } else {
      fcopy_surf (big, lit, cubeptr, j, type, dir);
    }

    for (i = 1; i < 6; i++) {
      j = (code + i) % 6;
      ni = ci + cubeptr->neigh[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:fcopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      lit = (block_values[ni]);
      if (lit == NULL) {
        fcopy_surf (big, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        fcopy_surf (big, lit, cubeptr, j, INS, dir);
      }
    }
    break;

  case CM_EDGE:
    for (i = 0; i < nouts; i++) {
      j = (cubeptr->elist[code][i]);

      ni = ci + cubeptr->face[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:fcopy_neigh: outs.neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      if (type == WRAP)
        lit = (block_values[ni]);
      else
        lit = (block_values[ci]);
      if (lit == NULL) {
        fcopy_surf (big, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        fcopy_surf (big, lit, cubeptr, j, type, dir);
      }

    }

    for (i = nouts; i < 6; i++) {
      j = (cubeptr->elist[code][i]);

      ni = ci + cubeptr->neigh[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:fcopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      lit = (block_values[ni]);
      if (lit == NULL) {
        fcopy_surf (big, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        fcopy_surf (big, lit, cubeptr, j, INS, dir);
      }
    }
    break;

  case CM_CORN:
    for (i = 0; i < nouts; i++) {
      j = (cubeptr->clist[code][i]);
      ni = ci + cubeptr->face[j];

      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:fcopy_neigh: outs.neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }

      if (type == WRAP){
        lit = (block_values[ni]);
        neigh_sb=ni;
      }else{
        lit = (block_values[ci]);
        neigh_sb=ci;
      }

      if (lit == NULL) {
        neigh_sb=ci;
        fcopy_surf (big, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        fcopy_surf (big, lit, cubeptr, j, type, dir);
      }

    }

    for (i = nouts; i < 6; i++) {
      j = (cubeptr->clist[code][i]);
      ni = ci + cubeptr->neigh[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:fcopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      neigh_sb=ni;
      lit = (block_values[ni]);
      if (lit == NULL) {
        neigh_sb=ci;
        fcopy_surf (big, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        fcopy_surf (big, lit, cubeptr, j, INS, dir);
      }
    }
    break;

  case XY_CORN:
    for (i = 0; i < nouts; i++) {
      j = (cubeptr->dlist[code][i]);

      ni = ci + cubeptr->face[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:fcopy_neigh: outs.neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      if (type == WRAP)
        lit = (block_values[ni]);
      else
        lit = (block_values[ci]);
      if (lit == NULL) {
        fcopy_surf (big, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        fcopy_surf (big, lit, cubeptr, j, type, dir);
      }

    }

    for (i = nouts; i < 6; i++) {
      j = (cubeptr->dlist[code][i]);

      ni = ci + cubeptr->neigh[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:fcopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      lit = (block_values[ni]);
      if (lit == NULL) {
        fcopy_surf (big, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        fcopy_surf (big, lit, cubeptr, j, INS, dir);
      }
    }
    break;

  default:
    fprintf (stderr, "ERROR: fcopy_neigh: Confusion! Confusion!\n");
    errors++;
    break;
  }
  return (errors);
}


int fcopy_surf (CA_FLOAT * big, CA_FLOAT * lit, Frame * cubeptr, int fcode, int type, int dir)
{

  int i, j, m, n, errors = 0;
  int pindex, index, axis, nins;
  int bigskip, litskip, start, skip, jump, offset, ioffset;
  int istart, iskip, ijump;

#ifdef DB_COPY
  int l_index, b_index, xindex, yindex, zindex;
#endif
  int ndim;
  int jcount = 0, ijcount = 0;
  int *ins;
  CA_FLOAT *pbig, *plit, *bigst, *litst, bdy_value;

  pbig = big;
  plit = lit;
  start = cubeptr->flist[fcode][START];
  skip = cubeptr->flist[fcode][SKIP];
  jump = cubeptr->flist[fcode][JUMP];
  istart = cubeptr->flist[fcode][ISTART];
  iskip = cubeptr->flist[fcode][ISKIP];
  ijump = cubeptr->flist[fcode][IJUMP];

  ins = cubeptr->ins;
  bdy_value = cubeptr->ivalue;

  axis = (fcode / 2) % 3;
  m = ins[(axis + 1) % 3];
  n = ins[(axis + 2) % 3];
  nins = ins[0] * ins[1] * ins[2];

  switch (type) {
  case WRAP:
    offset = cubeptr->flist[fcode][WOFF];
    ioffset = 0;
    break;

  case PAD:
  case FLUX_BDY:
    bdy_value *= cubeptr->dtbydx;
  case FIX_BDY:
    offset = cubeptr->flist[fcode][POFF];
    ioffset = 0;
    break;

  case INIT:
    offset = cubeptr->flist[fcode][POFF];
    ioffset = 0;
    break;

  case INS:
    offset = cubeptr->flist[fcode][POFF];
    ioffset = cubeptr->flist[fcode][NBOFF];
    break;

  default:
    fprintf (stderr, "fcopy_surf: unknown type value %i\n", type);
    exit (10);
    break;
  }

  if (start - offset <= 0) {
    fprintf (stderr, "ERROR:fcopy_surf :-( index blew up. )-:  %i\n", offset);
    errors++;
  }
  pbig += start - offset;
  bigst = pbig;
  plit += (istart + ioffset);
  litst = plit;
  litskip = iskip;
  bigskip = skip;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < m; ++i) {
      if (dir == FORWARD) {
        if (type == FIX_BDY) {
          *pbig = bdy_value;
        } else if (type == FLUX_BDY) {
          *pbig += bdy_value;
        } else {
          *pbig = *plit;
        }
      } else {

        if (type == FIX_BDY) {
          *plit = bdy_value;
        } else {
          *plit = MAX (*plit, *pbig);
        }
      }
      plit += litskip;
      pbig += bigskip;
    }
    jcount += jump;
    ijcount += ijump;
    plit = litst + ijcount;
    pbig = bigst + jcount;
  }
  return (errors);
}


int fcopy_mat_back (int type, CA_FLOAT * lit, CA_FLOAT * big, BB_struct * bp, CA_FLOAT * block_values[], int ci)
{

  int j, k, errors = 0;
  CA_FLOAT *pbig, *plit;
  int *ins, *outs;
  Frame *cubeptr;

  cubeptr = &(bp->cubeptr);
  pbig = big;
  plit = lit;
  ins = cubeptr->ins;
  outs = cubeptr->outs;
  pbig += cubeptr->flist[0][START];

  for (k = 0; k < ins[ZAXIS]; k++) {
    for (j = 0; j < ins[YAXIS]; j++) {
      memcpy (plit, pbig, ins[XAXIS] * sizeof (CA_FLOAT));
      plit += ins[XAXIS];
      pbig += outs[XAXIS];
    }
    pbig += (outs[XAXIS] * 2);
  }

  if (bp->ntsb == 1) {
    for (j = 0; j < 6; j++) {
      type = cubeptr->facectrl[j];
      if (type == WRAP) {
        errors += fcopy_surf (big, lit, cubeptr, j, type, BACK);
      } else if (type == PAD || type == FIX_BDY || type == FLUX_BDY) {

      } else {
        fprintf (stderr, "ERROR:fcopy_mat_back: Unknown Flag %i\n", type);
        errors++;
      }
    }
  }

  else {
    errors += fcopy_neigh (bp, big, ci, block_values, type, BACK);
  }
  if (errors > 1000)
    exit (1000);
  return (errors);
}


char const *rcs_id_fcopymat_c ()
{
  static char const rcsid[] = "$Id: fcopymat.c 1386 2008-09-25 11:47:24Z  $";

  return (rcsid);
}

