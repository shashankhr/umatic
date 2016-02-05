
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
/*      initface.c:                                             */
/* This subroutine calculates a list of offsets to step         */
/* a cell counter through the cells on a face of a cube         */
/* padded by one thickness of cells , and the offset necessary  */
/* to copy the  value to a cell within the padding, either      */
/* adjacent (pad offset) or opposite (wrap offset)              */
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Aug 25, 1998                                                 */
/****************************************************************/
/*      MODIFIED by:                                            */
/*  PDL: Aug 16, 1998                                           */
/*  RCA: Nov 9, 1998  -- fixed bugs in lookup table             */
/****************************************************************/
/****** To Do List **********************************************/
/****************************************************************/
/* 1) add stuff for copyback, and general size...               */
/****************************************************************/
/*RCS Id:$Id: initface.c 1341 2008-07-23 15:23:30Z  $*/
/*
*/
#include <stdio.h>
#include "machine.h"
#include "cube.h"

/******************************************************/
/*  quick subroutine to dump out the lookup table     */
/******************************************************/
int dump_facecode (Frame * cubeptr)
{
  int i, j;
  FILE *facecodefile;
  char facename[MAX_STRING_LEN];

  sprintf (facename, "facecode.txt");
  facecodefile = fopen (facename, "w");

  fprintf (facecodefile, "\n");

  for (j = 0; j < 9; ++j) {
    fprintf (facecodefile, "\t%5i ", j);
  }
  fprintf (facecodefile, "\n");

  for (i = 0; i < 6; ++i) {
    fprintf (facecodefile, "%3i: ", i);
    for (j = 0; j < 10; ++j) {
      fprintf (facecodefile, "\t%5i ", cubeptr->flist[i][j]);
    }
    fprintf (facecodefile, "\n");
  }
  fclose (facecodefile);
  return (0);
}                               /* end of dump_facecode */

/***************************************************************/
/* This subroutine calculates a list of offsets to step        */
/* a cell counter through the cells on a face of a cube        */
/* padded by one thickness of cells , and the offset necessary */
/* to copy the  value to a cell within the padding, either     */
/* adjacent (pad offset) or opposite (wrap offset)             */
/***************************************************************/
int init_facecode (Frame * cubeptr, int *ins, int dim)
{
/******************************************************/
/* local variables                                    */
/******************************************************/

                           /**********************************/
  int fcode;                    /* the face code 0=negx etc       */
  int i, j, ii, ndim, loop1;    /* loop counters                  */
  int nums[5], inums[5];        /* temporarily storecommon values */
  int row, level, irow, ilevel; /* precalculate internal row,level */
  int facelist[6][10];          /* temporarily store the facelist */
  int outs[3] = { 0, 0, 0 };    /* calculate dim of inside cube   */
                           /**********************************/

  fprintf (stderr, "Entering init_facecode...\n");
  fprintf (stderr, "ins: %d, %d, %d; dim: %d\n", ins[0], ins[1], ins[2], dim);
/******************************************************/
/* Set stuff to zero                                  */
/******************************************************/
  for (i = 0; i < 6; ++i) {
    for (j = 0; j < 9; ++j) {
      facelist[i][j] = 0;
    }
  }

/******************************************************/
/* Precalculate the dimensions of inside cube         */
/******************************************************/
  ndim = dim;

  for (i = 0; i < ndim; ++i) {
    if (ins[i] < 1) {
      printf ("ERROR: block < one cell thick!! - %i %i\n", ins[i], i);
      return (1);
    }
    outs[i] = ins[i] + 2;
  }

/******************************************************/
/* Precalculate values that get used over and over    */
/******************************************************/
  row = outs[XAXIS];
  irow = ins[XAXIS];

  switch (ndim) {
  case THREE_D:
    ilevel = ins[XAXIS] * ins[YAXIS];
    level = outs[XAXIS] * outs[YAXIS];
    break;
  case TWO_D_CA:
    level = 0;
    break;
  default:
    fprintf (stderr, "ERROR: Hyperspatial solid'n zone... %iD!\n", ndim);
    break;
  }

  nums[0] = 1;
  nums[1] = row;
  nums[2] = level;
  nums[3] = level + row + 1;    /* lowest inside cell index */
  nums[4] = ins[ZAXIS] * level + ins[YAXIS] * row + ins[XAXIS]; /* highest index inner block */

  inums[0] = 1;
  inums[1] = irow;
  inums[2] = ilevel;
  inums[3] = 0;
  inums[4] = ins[ZAXIS] * ins[YAXIS] * ins[XAXIS] - 1;

/*************************************************************/
/* set up the start values for all faces, even=neg odd = pos */
/*************************************************************/
  for (i = 0; i < 2 * ndim; i += 2) {
    facelist[i][START] = nums[3];       /* all neg faces start at the same place */
    facelist[i][ISTART] = inums[3];     /* all neg faces start at the same place */
    facelist[i + 1][START] = nums[4];
    facelist[i + 1][ISTART] = inums[4];
  }

/******************************************************/
/* set up the offsets for looping and copying         */
/* 1=skip,2=jump,3=pad offset, 4=wrap offset          */
/*                                                    */
/* The axis code is given by   (i/2)%3  where i is the*/
/* facecode.                                          */
/******************************************************/
  loop1 = 2 * ndim;
  for (i = 0; i < loop1; i += 2) {
    for (j = 1; j < 4; ++j) {
      ii = (i / 2 + j) % 3;
      facelist[i][j] = nums[ii];
      facelist[i + 1][j] = nums[ii] * (-1);

      facelist[i][j + 5] = inums[ii];
      facelist[i + 1][j + 5] = inums[ii] * (-1);
    }
    facelist[i][4] = facelist[i][3] * ins[(i / 2) % 3] * -1;
    facelist[i + 1][4] = facelist[i][3] * ins[(i / 2) % 3];
  }
  facelist[0][NBOFF] = irow - 1;
  facelist[1][NBOFF] = -(irow - 1);
  facelist[2][NBOFF] = (ins[XAXIS] * (ins[YAXIS] - 1));
  facelist[3][NBOFF] = -(ins[XAXIS] * (ins[YAXIS] - 1));
  facelist[4][NBOFF] = (ins[XAXIS] * ins[YAXIS] * (ins[ZAXIS] - 1));
  facelist[5][NBOFF] = -(ins[XAXIS] * ins[YAXIS] * (ins[ZAXIS] - 1));
/******************************************************/
/* Fix up a little problem in 2D case ;-)             */
/******************************************************/
  if (ndim == TWO_D_CA) {
    facelist[2][1] = facelist[2][2];
    facelist[2][2] = 0;
    facelist[3][1] = facelist[3][2];
    facelist[3][2] = 0;
  }

/******************************************************/
/*Now store the facecodes  in the cube structure      */
/******************************************************/
  for (i = 0; i < 6; ++i) {
    for (j = 0; j < 10; ++j) {
      cubeptr->flist[i][j] = facelist[i][j];
    }
  }

/**********************************************************/
/* and equate the cube structure members with the arrays. */
/**********************************************************/
  for (i = 0; i < 3; ++i) {
    cubeptr->outs[i] = outs[i];
    cubeptr->ins[i] = ins[i];
  }
  cubeptr->ndim = ndim;
  cubeptr->ivalue = 0.0;
  fprintf (stderr, "Exiting init_facecode...\n");
  dump_facecode (cubeptr);
  return (0);
}                               /* end of init_facelist */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_initface_c ()
{
  static char const rcsid[] = "$Id: initface.c 1341 2008-07-23 15:23:30Z  $";

  return (rcsid);
}

/* end of rcs_id_initface_c subroutine */
