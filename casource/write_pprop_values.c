
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

/*$Id: write_pprop_values.c 1339 2008-07-23 13:58:29Z  $*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "machine.h"
#include "blocks.h"
#include "pore.h"

void write_pprop_values (FILE * fp, P_str * pp)
{

  fprintf (fp, "\n\nPORE PROP STRUCTURE VALUES\n");
  fprintf (fp, "P_nmeth, %i\n", pp->P_nmeth);   /*GAUSS or STEP */
  fprintf (fp, "Binsize, %.10g\n", pp->Binsize);        /* size of temp. bins to use */
  fprintf (fp, "P_nmax, %.10g\n", pp->P_nmax);  /* nmax for the pore distrib. */
  fprintf (fp, "P_tn, %.10g\n", pp->P_tn);      /* center of gauss. dis. for pores */
  fprintf (fp, "P_tsig, %.10g\n", pp->P_tsig);  /*spread of pore distrib */
  fprintf (fp, "P_par1, %.10g\n", pp->P_par1);  /*user defined parameter */
  fprintf (fp, "P_par2, %.10g\n", pp->P_par2);  /*user defined parameter */
  /* for making the list of radius at steps of temperature */
  fprintf (fp, "P_ntrad, %i\n", pp->P_ntrad);
  fprintf (fp, "P_trad_max, %.10g\n", pp->P_trad_max);
  fprintf (fp, "P_trad_min, %.10g\n", pp->P_trad_min);
  fprintf (fp, "P_trad_step, %.10g\n", pp->P_trad_step);
  fprintf (fp, "P_limrad_perturb, %.10g\n", pp->P_limrad_perturb);
}

 /**/
/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *write_pprop_values_c ()
{
  static char const rcsid[] = "$Id: write_pprop_values.c 1339 2008-07-23 13:58:29Z  $";

  return (rcsid);
}
