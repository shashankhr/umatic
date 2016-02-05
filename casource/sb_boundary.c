
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

/*RCS Id:$Id: sb_boundary.c 1356 2008-08-18 13:41:15Z  $*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "blocks.h"
#include "machine.h"
#include "sb_diffuse.h"
#include "umat_matrix.h"

int sb_boundary (BB_struct * bp, int sbnum)
{
  int errflg = 0, fileflag = 0, errors = 0;
  int nx, ny, nz, tsteps;
  int *oni, *onip, *onend;
  SB_struct *sp;
  int t, i, j, k, l;
  CA_FLOAT partcoef;
  char *np, *ngr;
  int *op, *ogr;
  CA_FLOAT rs, rl;
  FILE *fp;

/* set up local neighbourhood */
/* use 6cell only for now     */
  oni = bp->nbhd.onq;           /*padded */
  onip = oni;
  onend = oni + 6;
/* set up local values and pointers */
  sp = bp->sb[sbnum];
  /*don't do it if sb doesn't exist */
  if (sp->open == FALSE)
    return (0);

  bp->cubeptr.curr = sbnum;
  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];

  op = ogr = bp->itmp_one;
  np = ngr = (char *) calloc (bp->ncsb, sizeof (char));

  /* make a copy of the grain array */

  errflg += icopy_matrix (PAD, ogr, sp->gr, bp, bp->gr_array, sbnum);   /* (flag, to, from, bp) */

  op = ogr + bp->cubeptr.flist[0][START];       /*rewind */
/*  LOOP                               */
/* Run through all cells updating as needed.    */
/************************************************/
  for (k = 0; k < nz; k++) {    /* loop cells in z direction */
    for (j = 0; j < ny; j++) {  /* loop cells in y direction */
      for (i = 0; i < nx; i++) {        /* loop cells in x direction */

        /*loop neighbours */
        for (onip = oni; onip < onend; onip++) {
          if (*(op + *onip) != *op) {
            *np = (char) 255;
            break;
          }
        }
        np++;
        op++;
      }                         /*x */
      op += 2;
    }                           /*y */
    op += 2 * (nx + 2);
  }                             /*z */

  fp = fopen ("sb_bdy.bin", "w");
  fwrite (ngr, sizeof (char), bp->ncsb, fp);
  fclose (fp);
  free (ngr);

  return (errflg);
}                               /* end of */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_sb_boundary_c ()
{
  static char const rcsid[] = "$Id: sb_boundary.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}

/* end of rcs_id_sb_boundary_c subroutine */
/*
*/
