
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

/*2:*/

#include <stdio.h>
#include <string.h>
#include "machine.h"
#include "blocks.h"
                   /*:2 *//*3: */

int icopy_surf (BB_struct * bp, int *new, int *old, Frame * cubeptr, int fcode, int type, int dir);
int icopy_neigh (BB_struct * bp, int *new, int ci, int *block_values[], int type, int dir);

                   /*:3 *//*4: */

int icopy_matrix (int type, int *new, int *old, BB_struct * bp, int *block_values[], int ci)
{
  int j, k, errors = 0;
  int *pold, *pnew;
  int *ins, *outs;
  int **bv_p;
  Frame *cubeptr;

  bv_p = block_values;
  cubeptr = &(bp->cubeptr);
  pold = old;
  pnew = new;
  ins = cubeptr->ins;
  outs = cubeptr->outs;
  pnew += cubeptr->flist[0][START];

  for (k = 0; k < ins[ZAXIS]; k++) {
    for (j = 0; j < ins[YAXIS]; j++) {
      memcpy (pnew, pold, ins[XAXIS] * sizeof (int));
      pold += ins[XAXIS];
      pnew += outs[XAXIS];
    }
    pnew += (outs[XAXIS] * 2);
  }

  if (bp->ntsb == 1) {
    for (j = 0; j < 6; j++) {
      type = cubeptr->facectrl[j];
      errors += icopy_surf (bp, new, old, cubeptr, j, type, FORWARD);
    }
  }

  else
    errors += icopy_neigh (bp, new, ci, bv_p, type, FORWARD);
  return (errors);
}

                   /*:4 *//*5: */

                   /*:5 *//*8: */

int icopy_neigh (BB_struct * bp, int *new, int ci, int *block_values[], int type, int dir)
{
  int code, nouts, errors = 0;
  int i, j;
  int *old;
  int ni;
  Frame *cubeptr;

  cubeptr = &(bp->cubeptr);
  code = bp->sb[ci]->code;
  nouts = bp->sb[ci]->nouts;

  switch (nouts) {
  case CM_INS:

    for (j = 0; j < 6; j++) {
/*6:*/

      ni = ci + cubeptr->neigh[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:icopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      old = (block_values[ni]);
      if (old == NULL) {
        icopy_surf (bp, new, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        icopy_surf (bp, new, old, cubeptr, j, INS, dir);
      }
/*:6*/
      ;
    }
    break;

  case CM_FACE:
    j = code;
/*7:*/

    ni = ci + cubeptr->face[j];
    if (ni > bp->ntsb) {
      fprintf (stderr, "ERROR:icopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
      exit (100);
    }
    if (type == WRAP)
      old = (block_values[ni]);
    else
      old = (block_values[ci]);
    if (old == NULL) {
      icopy_surf (bp, new, block_values[ci], cubeptr, j, PAD, dir);
    } else {
      icopy_surf (bp, new, old, cubeptr, j, type, dir);
    }

/*:7*/
    ;
    for (i = 1; i < 6; i++) {
      j = (code + i) % 6;
/*6:*/

      ni = ci + cubeptr->neigh[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:icopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      old = (block_values[ni]);
      if (old == NULL) {
        icopy_surf (bp, new, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        icopy_surf (bp, new, old, cubeptr, j, INS, dir);
      }
/*:6*/
      ;
    }
    break;

  case CM_EDGE:
    for (i = 0; i < nouts; i++) {
      j = (cubeptr->elist[code][i]);
/*7:*/

      ni = ci + cubeptr->face[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:icopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      if (type == WRAP)
        old = (block_values[ni]);
      else
        old = (block_values[ci]);
      if (old == NULL) {
        icopy_surf (bp, new, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        icopy_surf (bp, new, old, cubeptr, j, type, dir);
      }

/*:7*/
      ;
    }

    for (i = nouts; i < 6; i++) {
      j = (cubeptr->elist[code][i]);
/*6:*/

      ni = ci + cubeptr->neigh[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:icopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      old = (block_values[ni]);
      if (old == NULL) {
        icopy_surf (bp, new, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        icopy_surf (bp, new, old, cubeptr, j, INS, dir);
      }
/*:6*/
      ;
    }
    break;

  case CM_CORN:
    for (i = 0; i < nouts; i++) {
      j = (cubeptr->clist[code][i]);
/*7:*/

      ni = ci + cubeptr->face[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:icopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      if (type == WRAP)
        old = (block_values[ni]);
      else
        old = (block_values[ci]);
      if (old == NULL) {
        icopy_surf (bp, new, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        icopy_surf (bp, new, old, cubeptr, j, type, dir);
      }

/*:7*/
      ;
    }

    for (i = nouts; i < 6; i++) {
      j = (cubeptr->clist[code][i]);
/*6:*/

      ni = ci + cubeptr->neigh[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:icopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      old = (block_values[ni]);
      if (old == NULL) {
        icopy_surf (bp, new, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        icopy_surf (bp, new, old, cubeptr, j, INS, dir);
      }
/*:6*/
      ;
    }
    break;

  case XY_CORN:
    for (i = 0; i < nouts; i++) {
      j = (cubeptr->dlist[code][i]);
/*7:*/

      ni = ci + cubeptr->face[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:icopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      if (type == WRAP)
        old = (block_values[ni]);
      else
        old = (block_values[ci]);
      if (old == NULL) {
        icopy_surf (bp, new, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        icopy_surf (bp, new, old, cubeptr, j, type, dir);
      }

/*:7*/
      ;
    }

    for (i = nouts; i < 6; i++) {
      j = (cubeptr->dlist[code][i]);
/*6:*/

      ni = ci + cubeptr->neigh[j];
      if (ni > bp->ntsb) {
        fprintf (stderr, "ERROR:icopy_neigh: ins. neighbour index too high! ni %i ci %i\n", ni, ci);
        exit (100);
      }
      old = (block_values[ni]);
      if (old == NULL) {
        icopy_surf (bp, new, block_values[ci], cubeptr, j, PAD, dir);
      } else {
        icopy_surf (bp, new, old, cubeptr, j, INS, dir);
      }
/*:6*/
      ;
    }
    break;

  default:
    fprintf (stderr, "Confusion! Confusion! in icopy_neigh!\n");
    errors++;
    break;
  }
  return (errors);
}

                   /*:8 *//*10: */

int icopy_surf (BB_struct * bp, int *big, int *lit, Frame * cubeptr, int fcode, int type, int dir)
{

  int i, j, m, n, errors = 0;
  int pindex, index, axis, nins;
  int start, skip, jump, offset, ioffset;
  int istart, iskip, ijump;
  int bigskip, litskip;
  int ndim;
  int jcount = 0, ijcount = 0;
  int *ins;
  int *pbig, *plit, *bigst, *litst;

  pbig = big;
  plit = lit;
  start = cubeptr->flist[fcode][START];
  skip = cubeptr->flist[fcode][SKIP];
  jump = cubeptr->flist[fcode][JUMP];
  istart = cubeptr->flist[fcode][ISTART];
  iskip = cubeptr->flist[fcode][ISKIP];
  ijump = cubeptr->flist[fcode][IJUMP];

  ins = cubeptr->ins;

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
  case FIX_BDY:
  case FLUX_BDY:
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
    fprintf (stderr, "icopy_surf: unknown type value %i\n", type);
    exit (10);
    break;
  }

  if (start - offset <= 0) {
    fprintf (stderr, "Error :-( index blew up. )-:  %i\n", offset);
    errors++;
    return (errors);

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
        *pbig = *plit;
      } else {
        *plit = MAX (*plit, *pbig);
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

                      /*:10 *//*11: */

int icopy_mat_back (int type, int *lit, int *big, BB_struct * bp, int *block_values[], int ci)
{

  int j, k, errors = 0;
  int *pbig, *plit;
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
      memcpy (plit, pbig, ins[XAXIS] * sizeof (int));
      plit += ins[XAXIS];
      pbig += outs[XAXIS];
    }
    pbig += (outs[XAXIS] * 2);
  }

  if (bp->ntsb == 1) {

    for (j = 0; j < 6; j++) {
      type = cubeptr->facectrl[j];
      if (type == WRAP) {
        errors += icopy_surf (bp, big, lit, cubeptr, j, type, BACK);
      } else {
      }
    }
  }

  else {
    errors += icopy_neigh (bp, big, ci, block_values, type, BACK);
  }
  return (errors);
}

                      /*:11 *//*12: */

char const *rcs_id_icopymat_c ()
{
  static char const rcsid[] = "$Id: icopymat.c 1339 2008-07-23 13:58:29Z  $";

  return (rcsid);
}

/*:12*/
