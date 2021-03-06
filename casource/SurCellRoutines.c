
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

/*$Id: SurCellRoutines.c 1339 2008-07-23 13:58:29Z  $*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "machine.h"
#include "blocks.h"
#include "SurCell.h"

void alloc_SurCell (SurCell * surcp)
{
  fprintf (stderr, "alloc_SurCell: Initializing %i spaces ...\n", surcp->n_alloc);
  surcp->c_surfp = (int *) calloc (surcp->n_alloc, sizeof (int));
  if (surcp->c_surfp == NULL) {
    fprintf (stderr, "ERROR: init_SurCell: Allocation of %i memory units failed!\n", (int) (surcp->n_alloc * sizeof (int)));
    exit (95);
  }
  /* dynamic array of 3-element arrays ... */
  /* to store the xyz cell location */
  surcp->surf_xyz = (int (*)[3]) calloc (surcp->n_alloc, sizeof (int[3]));
  if (surcp->c_surfp == NULL) {
    fprintf (stderr, "ERROR: init_SurCell: Allocation of %i memory units failed!\n", (int) (surcp->n_alloc * sizeof (int)));
    exit (95);
  }
  fprintf (stderr, "alloc_SurCell: ... Done\n");
  return;
}

void init_SurCell (SurCell * surcp)
{
  surcp->n_alloc = SURF_NUC_ALLOC_START;
  surcp->ns_cell = 0;
}

void free_SurCell (SurCell * surcp)
{
  /* this routine is NOT responsible for freeing   */
  /* the threshold location storage, that is part  */
  /* of the subblock array collection and should be */
  /* freed there.                                  */
  if (surcp->n_alloc == 0) {
    fprintf (stderr, "ERROR:free_SurCell: tried to free an empty array!\n");
    if (surcp->c_surfp != NULL) {
      fprintf (stderr, "ERROR:free_SurCell: But the array is not NULL!\n");
      fprintf (stderr, "ERROR:free_SurCell: This is really, really bad!\n");
      exit (98);
    } else {
      fprintf (stderr, "Continuing anyways ...\nReturning ....\n");
      return;
    }
  }

  fprintf (stderr, "free_SurCell: Freeing the arrrays and zeroing the values\n");
  free (surcp->c_surfp);
  free (surcp->surf_xyz);
  surcp->c_surfp = NULL;
  surcp->surf_xyz = NULL;

  surcp->n_alloc = 0;
  surcp->ns_cell = 0;

  fprintf (stderr, "free_SurCell: ... Done\n");
  return;
}

void expand_SurCell (SurCell * surcp)
{
  static const int alloc_step = SURF_NUC_ALLOC_STEP;

  surcp->n_alloc += alloc_step;
  fprintf (stderr, "expand_SurCell: expanding to %i members\n", surcp->n_alloc);

  /* expand the array by the fixed step */
  surcp->c_surfp = (int *) realloc (surcp->c_surfp, surcp->n_alloc * sizeof (int));
  /* expand the array by the fixed step */
  surcp->surf_xyz = (int (*)[3]) realloc (surcp->surf_xyz, surcp->n_alloc * sizeof (int[3]));
  /* check if the realloc failed */
  if (surcp->c_surfp == NULL || surcp->surf_xyz == NULL) {
    fprintf (stderr, "ERROR: expand_SurCell: expansion to %i memory units failed!\n", (int) (surcp->n_alloc * sizeof (int)));
    exit (96);
  }

  /* zero the new elements */
  memset (surcp->c_surfp + (surcp->n_alloc - alloc_step), 0, alloc_step * sizeof (int));
  memset (surcp->surf_xyz + (surcp->n_alloc - alloc_step), 0, alloc_step * sizeof (int[3]));

  fprintf (stderr, "expand_SurCell: ... Done\n");
  return;
}

/************************************************/
/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
/************************************************/
char const *rcs_id_SurCellRoutines_c ()
{
  static char const rcsid[] = "$Id: SurCellRoutines.c 1339 2008-07-23 13:58:29Z  $";

  return (rcsid);
}

/* end of rcs_id subroutine */
