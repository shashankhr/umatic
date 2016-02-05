
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

#include <math.h>
#include <stdio.h>
#include "machine.h"

/*
   BOX-MULLER routine for returning a random number normally distributed
   about a mean of zero, with a standard deviation of 1. Adapted from Numerical Recipes in C by Press and Flannery.
*/

CA_FLOAT box_muller ()
{
  static int flag = 0;
  static CA_FLOAT gset = 0;
  CA_FLOAT v1, v2, fac, gdev, r;

  if (!flag) {
    do {
      v1 = 2 * drand48 () - 1;
      v2 = 2 * drand48 () - 1;
      r = v1 * v1 + v2 * v2;
    } while (r >= 1 || r == 0);
    fac = sqrt ((-2 * LOG (r) / r));
    gset = v1 * fac;
    gdev = v2 * fac;
    flag = 1;
  } else {
    gdev = gset;
    flag = 0;
  }
  return (gdev);
}                               /* end of box_muller */

/* return gaussian devates for a given mean and standard deviation */
/* by calling the box_muller routine. */

CA_FLOAT gaussdev (CA_FLOAT * params)
{
  CA_FLOAT stdev;
  CA_FLOAT mean;

  if (params == NULL || (params + 1) == NULL) {
    fprintf (stderr, "ERROR: gaussdev: You messed up, NULL pointer here.\n");
    exit (1);
  }

  mean = params[0];
  stdev = params[1];
  return ((stdev * box_muller ()) + mean);
}

#ifdef TEST_BOX_MULLER
void main ()
{
  int i;
  CA_FLOAT stdev = 2;
  CA_FLOAT mean = 10;
  CA_FLOAT params[2];

  params[0] = mean;
  params[1] = stdev;
  for (i = 0; i < 10000; i++) {
    printf ("%.5g\n", gaussdev (params));
  }

}
#endif

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file */
char const *rcs_id_gaussdev_c ()
{
  static char const rcsid[] = "$Id: gaussdev.c 1339 2008-07-23 13:58:29Z  $";

  return (rcsid);
}

/* end of rcs_id_subroutine */
