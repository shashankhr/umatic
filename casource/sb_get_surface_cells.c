
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

/*$Id: sb_get_surface_cells.c 1386 2008-09-25 11:47:24Z  $*/

#include <stdio.h>
#include <stdlib.h>
#include "machine.h"
#include "blocks.h"
#include "umat_matrix.h"
#include "SurCell.h"
#include "SurCellRoutines.h"

void sb_get_surface_cells (BB_struct * bp, int sbnum)
{
  CA_FLOAT *padfs, *padfs_start, *fs, *fs_start;
  SB_struct *sp;
  int i, j, k;
  int umat_index = 0;
  int num_mould_nuc = 0, num_not_casting = 0;
  int errflg = 0;
  int *oni, *onip, *onend;

  /* stored in the subblock structure */
  SurCell *surface;

  fprintf (stderr, "sb_get_surface_cells: finding the cells at the surface between\n");
  fprintf (stderr, "         the casting and the mould for subblock %i.\n", sbnum);

  sp = bp->sb[sbnum];
  surface = &(sp->surface);

  /* set up local neighbourhood */
  /* use 6cell only for now     */
  oni = bp->nbhd.onq;           /*padded */
  onip = oni;
  onend = oni + 6;
  /* set up local values and pointers */
  padfs = padfs_start = bp->ftmp_one;
  fs = fs_start = sp->c_fs;

  /*****************************************/
  /* make a copy of the fs array           */
  /* this is needed to locate cells        */
  /* neighbouring on mould in the          */
  /* next subblock or in a periodic domain */
  /*****************************************/
  errflg += fcopy_matrix (PAD, padfs, fs, bp, bp->c_fs_values, sbnum);  /* (flag, :to, from, bp) */
  padfs = padfs_start + bp->cubeptr.flist[0][START];
  fs = fs_start;

  for (k = 0; k < bp->nc[2]; k++) {     /* loop cells in z direction */
    for (j = 0; j < bp->nc[1]; j++) {   /* loop cells in y direction */
      for (i = 0; i < bp->nc[0]; i++) { /* loop cells in x direction */
        /* find neighbour cells that are not in the casting */
        if (*fs == NOT_CASTING) {
          num_not_casting++;
        } else {
          for (onip = oni; onip < onend; onip++) {
            if ((*(padfs + *onip)) == NOT_CASTING) {
              /* Neigbour is not in casting so */
              /* set a threshold here */
              /* if all the allocated spaces are used up */
              if (surface->ns_cell >= surface->n_alloc) {
                /* expand the array */
                expand_SurCell (surface);
              }

              /* record the umat_index of the cell */
              surface->c_surfp[surface->ns_cell] = umat_index;
              /* record the grid location of the cell too! */
              (*(surface->surf_xyz + surface->ns_cell))[0] = i;
              (*(surface->surf_xyz + surface->ns_cell))[1] = j;
              (*(surface->surf_xyz + surface->ns_cell))[2] = k;

              surface->ns_cell++;
              /* then skip to the next cell */
              /* break out of the NEIGHBOUR loop */
              break;
            }                   /* end of Neighbour Not_casting test */
          }                     /* end of NEIGHBOUR loop */
        }                       /* end NOT_CASTING test */
        fs++;
        padfs++;
        umat_index++;
      }                         /*x */
      padfs += 2;
    }                           /*y */
    padfs += 2 * (bp->nc[0] + 2);
  }                             /*z */

  sp->nmould = num_not_casting;

  fprintf (stderr, "sb_get_surface_cells: %i has %i mould cells \n", sbnum, num_not_casting);
  fprintf (stderr, "sb_get_surface_cells: %i has %i interface cells \n", sbnum, surface->ns_cell);
  fprintf (stderr, "sb_get_surface_cells: finished subblock %i\n", sbnum);
  return;
}
