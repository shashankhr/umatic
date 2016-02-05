
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

/*RCS Id:$Id: recr_extra.c 1341 2008-07-23 15:23:30Z  $*/
/*
*/
#include <math.h>
#include "machine.h"

CA_FLOAT gro_recr (CA_FLOAT gg_c, CA_FLOAT temp, CA_FLOAT lnZ, CA_FLOAT dt)
{
  CA_FLOAT grorate;

  /* G = g (1/T) exp(-1/T)(lnZ-1)^2 */

  grorate = gg_c * (1 / temp) * (exp (-1 / temp)) * (lnZ - 1) * (lnZ - 1);
  return (dt * grorate);
}                               /* end gro_recr */

CA_FLOAT nuc_rate_recr (CA_FLOAT nuc_c, CA_FLOAT temp, CA_FLOAT lnZ, CA_FLOAT eps)
{
  CA_FLOAT nucrate;

  nucrate = nuc_c * lnZ * (exp (eps) + exp (-eps) + 1) * exp (-1 / temp);
  return (nucrate);
}

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_recr_extra_c ()
{
  static char const rcsid[] = "$Id: recr_extra.c 1341 2008-07-23 15:23:30Z  $";

  return (rcsid);
}

/* end of rcs_id_recr_extra_c subroutine */
