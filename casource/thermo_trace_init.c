
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
/* thermo_trace_init.c:                                         */
/* Subroutine to initialise all variables and to                */
/* read in the values for a file.                               */
/*                                                              */
/* This subroutine is pard of the code to interpolate the       */
/* temperature at a point in the CA code using values           */
/* from a thermocouple trace                                    */
/****************************************************************/
/****************************************************************/
/* Written by X. Xu, P.D. Lee & R.C. Atwood, Imperial College   */
/* Dec 7 , 1998                                                  */
/****************************************************************/
/*      MODIFIED by:                                            */
/*                                                              */
/****************************************************************/
/****** To Do List **********************************************/
/*General:                                                      */
/* 1)                                                           */
/****************************************************************/
/*RCS Id:$Id: thermo_trace_init.c 1341 2008-07-23 15:23:30Z  $*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* include header files requred by subroutines */
#include "machine.h"
#include "tcpl.h"

/****************************************************************/
/****************************************************************/
/* thermo_trace_init:                                           */
/* Subroutine to init arrays and read in values.                */
/****************************************************************/
/* NOTE: the current version uses constant time step            */
/****************************************************************/
/** \todo  generalise to arbitrary (within reason) time steps -- input - coupled  */
/*                                                              */
/*                                                              */
/****************************************************************/
int thermo_trace_init (TC_str * tc)
{
  int errors = 0;
  FILE *fp;
  int i, j;                     /* temporary counters */

   /************************************************/
  /* Open the ascii trace infile to read          */
  /* the temperatures from.                       */
   /************************************************/
  if ((fp = fopen (tc->fn, "r")) == NULL) {
    printf ("ERROR: thermo_trace_init can not open input file [%s].\n", tc->fn);
    exit (1);
  }

   /****************************************************************/
  /* Read in the header line which holds:                         */
  /*    extra information                                         */
   /****************************************************************/
  fgets (tc->TC_header, MAX_STRING_LEN, fp);
  /*fscanf(fp, "%s\n", tc->TC_header); */
  fprintf (stderr, "TC Header info: %s\n", tc->TC_header);
   /****************************************************************/
  /* Read in the next line which holds:                         */
  /*    number of lines expected, timestep                        */
   /****************************************************************/
  fscanf (fp, "%d %g", &(tc->Nlines), &(tc->Tstep));
  fprintf (stderr, "TC Nlines, Tstep info: %i %f\n", tc->Nlines, tc->Tstep);

   /*************************************************/
  /* Malloc the arrays to hold time and temp values */
   /*************************************************/
  /* check for a bad number */
  if (tc->Nlines <= 1) {
    fprintf (stderr, "ERROR:thermo_trace_init: Not enough data lines specified! %i\n", tc->Nlines);
    exit (0);
  }
  /* malloc array of time values  */

  if (!(tc->Time = (float *) malloc (tc->Nlines * sizeof (float)))) {
    fprintf (stderr, "ERROR: thermo_trace_init  time array  malloc failed\n");
    return (2);
  }

  /* malloc array of temperatures */
  if (!(tc->Temp = (float *) malloc (tc->Nlines * sizeof (float)))) {
    fprintf (stderr, "ERROR: thermo_trace_init  temp array  malloc failed\n");
    return (3);
  }

   /************************************************/
  /* Now read in all the rest of the file.        */
   /************************************************/
  /* loop through reading the r locations         */
  for (i = 0; i < tc->Nlines; i++) {
    if (fscanf (fp, "%g %g", &tc->Time[i], &tc->Temp[i]) == NULL) {
      fprintf (stderr, "ERROR: thermo_trace_init ran out of data! %i\n", i);
      errors++;
    }
  }
  fclose (fp);
  return (errors);
}

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_thermo_trace_init_c ()
{
  static char const rcsid[] = "$Id: thermo_trace_init.c 1341 2008-07-23 15:23:30Z  $";

  return (rcsid);
}

/* end of rcs_id_thermo_trace_init_c subroutine */
/*
*/
