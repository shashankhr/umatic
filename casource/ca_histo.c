
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
#include <math.h>
#include "machine.h"
#include "umat_histo.h"

/* if stat flag is FIRST_CALL, print the bin values otherwise just the histo.*/
/* if there are more than 255 bins then you cannot easily read it into excel!*/
int umat_histo (FILE * listfile, Histo_struct * hp, CA_FLOAT * list, int stepnum)
{

  int errors = 0;
  int bin, i;
  CA_FLOAT *listp;
  CA_FLOAT maxbin;
  CA_FLOAT binsize;
  CA_FLOAT thisbin;
  int *g_nuc_histo;

  /* no grains */
  if (list == NULL)
    return (0);

  /*there are grains */
  maxbin = hp->minbin + hp->nbins * hp->binsize;
  g_nuc_histo = calloc (hp->nbins + 1, sizeof (int));

  for (listp = list; listp < list + hp->ndata; listp++) {
    bin = (int) FLOOR ((*listp - hp->minbin) / hp->binsize);
    bin = bin < 0 ? 0 : bin;
    bin = bin > hp->nbins ? hp->nbins : bin;
    (*(g_nuc_histo + bin))++;
  }

  fprintf (listfile, "%i", stepnum);
  /* print the statistical values */
  for (i = 0; i < 4; i++) {
    fprintf (listfile, ",%.5g", hp->stat[i]);
  }
  /* print the histogram values */
  for (bin = 0; bin <= hp->nbins; bin++) {
    fprintf (listfile, ",%i", (*(g_nuc_histo + bin)));
  }
  fprintf (listfile, "\n");

  free (g_nuc_histo);
  return (errors);
}

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file */
char const *rcs_id_umat_histo_c ()
{
  static char const rcsid[] = "$Id: ca_histo.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}

/* end of rcs_id_subroutine */
