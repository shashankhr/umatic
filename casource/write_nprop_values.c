
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

/*$Id: write_nprop_values.c 1339 2008-07-23 13:58:29Z  $*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "machine.h"
#include "blocks.h"

void write_nprop_values (FILE * fp, Nuc_str * np)
{
  fprintf (fp, "\n\nNUC PROP STRUCTURE VALUES\n");
  fprintf (fp, "nmodel, %i\n", np->nmodel);     /* type of nucleation model:            */
  fprintf (fp, "nhet, %i\n", np->nhet); /* number of different het. nuclei      */
  fprintf (fp, "nareanuc, %i\n", np->nareanuc); /* number of different area nucl. spec. */
  fprintf (fp, "ngr, %i\n", np->ngr);   /* number of grains                     */
  fprintf (fp, "gd_max, %.10g\n", np->gd_max);  /* max. grain density per [cm^3]        */
  fprintf (fp, "gd_max_total, %i\n", np->gd_max_total); /* total max. # grains for array size   */
  fprintf (fp, "oriented, %i\n", np->oriented); /* Orientation Calc: True/Flase         */
  fprintf (fp, "NucParams, %.10g, %.10g, %.10g, %.10g\n", np->NucParams[0], np->NucParams[1], np->NucParams[2], np->NucParams[3]);      /* the nucleation parameters array      */
  fprintf (fp, "perturb_on, %i\n", np->perturb_on);
  fprintf (fp, "n_perturb, %i\n", np->n_perturb);
}

/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *write_nprop_values_c ()
{
  static char const rcsid[] = "$Id: write_nprop_values.c 1339 2008-07-23 13:58:29Z  $";

  return (rcsid);
}
