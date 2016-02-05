
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

/*$Id: setup_mould_src_function.c 1339 2008-07-23 13:58:29Z  $*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* include header files requred by subroutines */
#include "machine.h"
#include "read_ctrl.h"
#include "blocks.h"
#include "mould_sources.h"
void setup_mould_src_function (Ctrl_str * cp, BB_struct * bp, Solute_props * sp)
{

       /*******************************************/
  /* set up the mould source function */
       /*******************************************/
  fprintf (stderr, "SETUP_MOULD_SRC_FUNCTION: setting up the mould source function...\n");
  switch (sp->mould_src) {
  case (MouldSourceNONE):
    fprintf (stderr, "SETUP_MOULD_SRC_FUNCTION: Zero source function defined.\n");
    sp->mould_src_func = none_mould_src;
    break;
  case (MouldSourceCONSTANT):
    /* constant */
    fprintf (stderr, "SETUP_MOULD_SRC_FUNCTION: mould source function is CONSTANT\n");
    sp->mould_src_func = const_mould_src;
    break;
  case (MouldSourceDIFF):
    /* difference */
    fprintf (stderr, "SETUP_MOULD_SRC_FUNCTION: mould source function is DIFF\n");
    sp->mould_src_func = diff_mould_src;
    break;
  case (MouldSourcePLIN):
    /* piecewise linear */
    fprintf (stderr, "SETUP_MOULD_SRC_FUNCTION: mould source function is PLIN (piecewise linear)\n");
    sp->mould_src_func = plin_mould_src;
    break;
  case (MouldSourceFLUX):
    /* piecewise linear */
    fprintf (stderr, "SETUP_MOULD_SRC_FUNCTION: mould source function is FLUX \n");
    sp->mould_src_func = flux_mould_src;
    break;
  default:
    fprintf (stderr, "ERROR: SETUP_MOULD_SRC_FUNCTION: Invalid function defined, %i.\n", sp->mould_src);
    exit (0);
    break;
  }
  if (sp->mould_src_pert) {
    /* perturbed function */
    fprintf (stderr, "SETUP_MOULD_SRC_FUNCTION: mould source function is PERTURBED\n");
    sp->mould_src_func = perturb_mould_src;
  }
       /*******************************************/
  return;
}

/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *setup_mould_src_function_c ()
{
  static char const rcsid[] = "$Id: setup_mould_src_function.c 1339 2008-07-23 13:58:29Z  $";

  return (rcsid);
}
