
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
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/****************************************************************/
/*      sb_line_int.c:                                          */
/*                                                              */
/* The subroutine to calculate grain size by line intercept     */
/*                                                              */
/*                                                              */
/*                                                              */
/*                                                              */
/****************************************************************/
/*RCS ID: $Id: sb_line_int.c 1342 2008-07-23 15:45:00Z  $*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "machine.h"
#include "blocks.h"

int sb_line_int (CA_FLOAT * line_res, BB_struct * bp, int sbnum)
{

  int errors = 0;
  int i, j;
  SB_struct *sp;
  int *gr, *gr_slice;
  int *grp, old_gr;
  int diag;
  CA_FLOAT length = 0;
  int nxny, nx, ny, nz, slicenum;
  int num_gr = 0;
  int nlines = 4;
  int *start[4];
  int stop[4];
  int step[4];
  CA_FLOAT linelen[4];

  if (bp->sb[sbnum]->open != SB_OPEN)
    return;

  sp = bp->sb[sbnum];
  gr = sp->gr;
  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  nxny = nx * ny;
  diag = MIN (nx, ny);
  slicenum = bp->ctrl->grainslice;      /*from ctrl file */
  if ((slicenum > nz) || (slicenum < 0)) {
    fprintf (stderr, "ERROR: sb_line_int: slicenum %i out of range %i\n", slicenum, nz);
    errors++;
    return (errors);
  }
  gr_slice = (gr + nxny * slicenum);

  /* 0 = x line */
  /* 1 = y line */
  /* 2 = diagonal ll - ur */
  /* 3 = diagonal lr - ul */

  start[0] = gr_slice + (int) FLOOR (0.5 * (CA_FLOAT) ny) * nx;
  start[1] = gr_slice + (int) FLOOR (0.5 * (CA_FLOAT) nx);
  start[2] = gr_slice;
  start[3] = gr_slice + nx;

  stop[0] = nx;
  stop[1] = ny;
  stop[2] = MIN (nx, ny);
  stop[3] = MIN (nx, ny);

  step[0] = 1;
  step[1] = nx;
  step[2] = nx + 1;
  step[3] = nx - 1;

  linelen[0] = (CA_FLOAT) nx;
  linelen[1] = (CA_FLOAT) ny;
  linelen[2] = (CA_FLOAT) (MIN (nx, ny)) * SQRT2;
  linelen[3] = (CA_FLOAT) (MIN (nx, ny)) * SQRT2;

  /* this loop could be used in subroutine which gets a structure     */
  /* defining an arbitrary number of various lines to do the intercept */
  /* method along.                                                    */

  for (j = 0; j < nlines; j++) {        /* loop through defined lines */
    grp = start[j];
    old_gr = *grp;
    for (i = 0; i < stop[j]; i++) {     /* traverse the line */
      if (*grp != old_gr) {     /* test for different grain */
        num_gr++;
        old_gr = *grp;
      }                         /* end test for different grain */
      grp += step[j];
    }                           /*end of i loop traverse line */
    length += linelen[j];
  }                             /* end of j loop through lines */
  if (num_gr != 0) {
    *line_res = (bp->size_c[0] * length) / (CA_FLOAT) num_gr;
  } else {
    *line_res = 0;
  }

  return (errors);
}                               /* end of sb_line_int */

/************************************************/
/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
/************************************************/
char const *rcs_id_sb_line_int_c ()
{
  static char const rcsid[] = "$Id: sb_line_int.c 1342 2008-07-23 15:45:00Z  $";

  return (rcsid);
}

/* end of rcs_id_sb_line_int_c subroutine */
/*
$Log$
Revision 11.1  2006/03/01 18:20:40  rcatwood
Merging polycomponent and gas with meltback

Revision 10.3  2005/12/01 14:38:02  rcatwood
Merged xly_05 changes into the main trunk
Primarily involving melt-back

Revision 10.1.2.2  2005/11/23 18:18:53  rcatwood
Result of merging mould_source and xly meltback+curvature 2d versions

Revision 10.1  2005/11/03 11:56:47  rcatwood
New version number -- using mould_src as base

Revision 8.1.12.2  2005/11/02 11:55:06  rcatwood
Fixing up the revision nubmer after loss of repository

Revision 9.1.4.1  2004/04/07 11:18:31  rcatwood
Fixed several division errors, added section to activate f.p.e. trapping

Revision 9.1  2003/08/14 14:38:40  rcatwood
Working merge with decentered/porosity/external, also including
Ali Chirazi's multicomponent (not tested in this version)

Revision 8.1.6.1  2003/01/22 16:53:47  rcatwood
Almost working read_fg version

Revision 8.1  2002/10/17 17:01:03  rcatwood
New version number! for decentered/porosity merge! Alpha Version!

Revision 7.5  2002/10/17 16:52:38  rcatwood
Merge from branch: combined Robert (porosity) and Wei (decentered octahedron) versions

Revision 7.4.10.1  2002/09/03 13:31:59  rcatwood
Merged with reorganized allocation routines, and adjusted nucleation to compartmentalize
the grain information assignment.

Revision 7.4  2002/02/14 13:15:02  rcatwood
Added write_block option instead of cpp definintion.

Revision 7.3  2001/03/13 11:48:33  rcatwood
fixed some comments and added sb_line_int id subroutine

Revision 7.2  2001/02/22 13:00:24  rcatwood
Included pore reallocation and stop temp (as cpp macro)
fixed x-y bug in line intercept routine

Revision 7.1  2001/02/19 19:28:47  rcatwood
fixed histo
for grains

and also make TcTrace mode override const. cooling rate

*/
