
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
/*      read_cap_ca.c:						*/
/*  Read in the geo info from a CAP CA file			*/
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Jul 1, 1998                                                  */
/****************************************************************/
/* 	MODIFIED by:						*/
/*  PDL: Aug 19, 1998						*/
/****************************************************************/
/*RCS Id:$Id: read_cap_ca.c 1356 2008-08-18 13:41:15Z  $*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "machine.h"
#include "fem.h"
#include "read_ctrl.h"
#include "blocks.h"

#ifdef DEBUG
FILE *efd;
#endif /* DEBUG */

/****************************************************************/
/*      Subroutine to read cap file                             */
/* This subroutine only reads the header and mask               */
/* array, the cell info is read in per sub-block                */
/****************************************************************/
/* pass in a blank BB structure, return with header filled    */

int read_cap_ca (Ctrl_str * cp, BB_struct * ptr_bb)
{
  /* just a stub */
#ifdef JUNK
  char command[MAX_STRING_LEN];
  int echo = TRUE;
  int i, j;                     /* tmp counters */
  CA_FLOAT tmp;                 /* tmp CA_FLOAT var. */
  FILE *fp;                     /* tmp filehandle */

/************************************************/
/* Read in geo from *****************************/
/************************************************/
/************************************************/
/* open BB input file and create output file	*/
/************************************************/
/* Open up input file				*/
  if ((fp = fopen (cp->fn_cap, "r")) == NULL) {
    fprintf (stderr, "Error: can't open input file [%s]\n", cp->fn_cap);
    /* ifd = stdin; echo = FALSE; */
  }
  ptr_bb->ibb_fd = fp;
  cp->fd_cap = fp;

/************************************************/
/* Read all the header data from the data file	*/
/************************************************/
  /* header 1 */
  fread (ptr_bb->header, sizeof (char), 256, fp);

  /* num. sub-blocks x, y, z */
  for (i = 0; i < 3; i++)
    fread (&(ptr_bb->nsb[i]), sizeof (int), 1, fp);
  fprintf (stderr, "nsb: %d, %d, %d\n", ptr_bb->nsb[0], ptr_bb->nsb[1], ptr_bb->nsb[2]);

  /* num. cells x, y, z */
  for (i = 0; i < 3; i++)
    fread (&(ptr_bb->tnc[i]), sizeof (int), 1, fp);
  fprintf (stderr, "ncells: %d, %d, %d\n", ptr_bb->tnc[0], ptr_bb->tnc[1], ptr_bb->tnc[2]);

  /* big block origin */
  for (i = 0; i < 3; i++)
    fread (&(ptr_bb->orig_bb[i]), sizeof (int), 1, fp);
  fprintf (stderr, "orig_bb: %f, %f, %f\n", ptr_bb->orig_bb[0], ptr_bb->orig_bb[1], ptr_bb->orig_bb[2]);

  /* cell dimensions */
  for (i = 0; i < 3; i++)
    fread (&(ptr_bb->size_c[i]), sizeof (int), 1, fp);
  fprintf (stderr, "cell dim: %f, %f, %f\n", ptr_bb->size_c[0], ptr_bb->size_c[1], ptr_bb->size_c[2]);

  /* num. casting nodes, elements */
  fread (&(ptr_bb->fem.nnodes), sizeof (int), 1, fp);
  fread (&(ptr_bb->fem.nelm), sizeof (int), 1, fp);
  fprintf (stderr, "nnodes: %d, nelm: %d\n", ptr_bb->fem.nnodes, ptr_bb->fem.nelm);

  /* extra space */
  fread (ptr_bb->tailer, sizeof (char), 200, fp);

/************************************************/
/* Finished reading header files, now           */
/* calc. all subsiduary values...		*/
/************************************************/
  ptr_bb->ntsb = ptr_bb->nsb[0] * ptr_bb->nsb[1] * ptr_bb->nsb[2];
  ptr_bb->ncsb = (ptr_bb->tnc[0] * ptr_bb->tnc[1] * ptr_bb->tnc[2]) / ptr_bb->ntsb;
  fprintf (stderr, "total #sb: %d, # cells/sb: %d\n", ptr_bb->ntsb, ptr_bb->ncsb);
  for (i = 0; i < 3; i++)
    ptr_bb->nc[i] = ptr_bb->tnc[i] / ptr_bb->nsb[i];
  fprintf (stderr, "#cells in x,y,z: %d, %d, %d\n", ptr_bb->nc[0], ptr_bb->nc[1], ptr_bb->nc[2]);
  for (i = 0; i < 3; i++)
    ptr_bb->size_bb[i] = (CA_FLOAT) ptr_bb->tnc[i] * ptr_bb->size_c[i];
  fprintf (stderr, "bigblock size: %f, %f, %f\n", ptr_bb->size_bb[0], ptr_bb->size_bb[1], ptr_bb->size_bb[2]);

/************************************************/
/* Malloc all SB mask for indexing CAP_CA file  */
/************************************************/
/* malloc an array to hold submask */
  if (!(ptr_bb->sb_mask = (long *) malloc (ptr_bb->ntsb * sizeof (long *)))) {
    fprintf (stderr, "ERROR: sb_mask array malloc failed\n");
    return (1);
  }

  fprintf (stderr, "before sb_mask nsb: %d, %d, %d\n", ptr_bb->nsb[0], ptr_bb->nsb[1], ptr_bb->nsb[2]);

/* read sub-block mask */
  fread (ptr_bb->sb_mask, sizeof (long), ptr_bb->ntsb, fp);
  for (i = 0; i < ptr_bb->ntsb; i++)
    fprintf (stderr, "sb_mask[%d]: %ld\n", i, ptr_bb->sb_mask[i]);

/*************************************/
/* Print out checks on input data... */
/*************************************/
  fprintf (stderr, "sb_mask nsb: %d, %d, %d\n", ptr_bb->nsb[0], ptr_bb->nsb[1], ptr_bb->nsb[2]);

  fprintf (stderr, "Exiting read_cap_ca().\n");
#endif /*JUNK*/
    return (0);
}                               /* end of read_cap_ca subroutine */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_read_cap_umat_c ()
{
  static char const rcsid[] = "$Id: read_cap_ca.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}

/* end of rcs_id_read_cap_umat_c subroutine */
/*
*/
