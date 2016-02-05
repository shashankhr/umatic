
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

/*$Id: setup_temp_func.c 1342 2008-07-23 15:45:00Z  $*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* include header files required by subroutines */
#include "machine.h"
#include "read_ctrl.h"
#include "blocks.h"
#include "temp_calc.h"
void setup_temp_function (Ctrl_str * cp, BB_struct * bp)
{

       /*******************************************/
  /* set up the temperature function */
       /*******************************************/
  fprintf (stderr, "CA_SOLID: setting up the temperature function...\n");
  if ((cp->fgrid_input + cp->con_cast + cp->external > 1) || (cp->fgrid_input + cp->con_cast + cp->external < 0)) {
    fprintf (stderr, "ERROR: SETUP_TEMP_FUNCTION: Invalid option.\n");
    fprintf (stderr, "SETUP_TEMP_FUNCTION: cannot select multiple options.\n");
    fprintf (stderr, "                     fgrid:   %i\n", cp->fgrid_input);
    fprintf (stderr, "                     concast: %i\n", cp->con_cast);
    fprintf (stderr, "                     external: %i\n", cp->external);
    fprintf (stderr, "exiting...\n");
    exit (32);
  }

       /**********************************************************/
  /* Set up the appropriate temperature function            */
  /* by copying the pointer to the chosen function          */
  /* into the bigblock temperature function pointer variable */
       /**********************************************************/
  if (cp->fgrid_input) {
    /* finite element grid input */
    fprintf (stderr, "CA_SOLID: temperature function is FGRID_INPUT\n");
    bp->cell_temp_func = fg_temp_calc;
  } else if (cp->con_cast) {
    /* calculated directional gradient */
    fprintf (stderr, "CA_SOLID: temperature function is (concast) DIR_TEMP\n");
    bp->cell_temp_func = dir_temp_calc;
  } else if (cp->external) {
    /* external temperature already calculated */
    /**  \todo  check order of external temperature update - external - coupled */
    fprintf (stderr, "CA_SOLID: temperature function is EXTERNAL_TEMP\n");
    bp->cell_temp_func = external_temp_calc;
  } else {
    /* isothermal (const. cooling, thermocouple input, */
    fprintf (stderr, "CA_SOLID: temperature function is CONST_TEMP\n");
    bp->cell_temp_func = const_temp_calc;
  }
       /*******************************************/
}

/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *setup_temp_func_c ()
{
  static char const rcsid[] = "$Id: setup_temp_func.c 1342 2008-07-23 15:45:00Z  $";

  return (rcsid);
}
