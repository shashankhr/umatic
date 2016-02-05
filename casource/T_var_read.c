
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
/* T_var_read.c:                                                */
/* Subroutine to readin the temperature field at each step      */
/* in a transient state solution of the heat, mass and          */
/* momentum transfer in VAR ingots.                             */
/****************************************************************/
/* Written by X. Xu                          Imperial College   */
/* Feb.18, 2000                                                 */
/****************************************************************/
/*      MODIFIED by:                                            */
/*                                                              */
/*                                                              */
/****************************************************************/
/****** To Do List **********************************************/
/*General:                                                      */
/* 1)                                                           */
/****************************************************************/
#include <stdio.h>
#include <math.h>

/* include header files requred by subroutines */
#include "machine.h"
#include "blocks.h"
#include "fidap.h"              /* included for def. of FGrid_struct */
/****************************************************************/
/* T_var_read.c:                                                */
/* Subroutine to read in values.                                */
/****************************************************************/
/* Input Variables:                                             */
/*      		*Note: currently FIDAP geom. is trans.  */
/*      		to [mm] from [m] since [mm] used in CA  */
/*   *fg:		ptr to the FGrid_str structure that     */
/*       		holds all FIDAP geometry and T's data.  */
/*   r:		        radial dist. in macro model (m)         */
/*   z:		 	height in macro model (m)               */
/*   time:		time in the CA model             [s]    */
/*      		                                        */
/* Output Variables:    NONE                                    */
/*      		                                        */
/* Returned Value:        0                                     */
/****************************************************************/

void T_var_read (FGrid_str * fg, CA_FLOAT time)
{
  int i, j;
  CA_FLOAT T;

  fscanf (fg->fp, "" SCAN_F " " SCAN_F " ", &(fg->time_var), &(fg->h_ingot));

/* loop through reading the z locations         */
  for (i = 0; i < fg->nz; i++) {
    fscanf (fg->fp, "" SCAN_F "", &fg->z[i]);
  }
/* loop through reading the r locations         */
  for (j = 0; j < fg->nr; j++) {
    fscanf (fg->fp, "" SCAN_F "", &fg->r[j]);
  }
/* loop through reading the temperatures        */
  for (j = 0; j < fg->nr; j++) {
    for (i = 0; i < fg->nz; i++) {
      fscanf (fg->fp, "" SCAN_F "", &T);

/* Calculate the T array index from the i and j values */
      fg->Fidap_T[(fg->nr) * i + j] = T;
      if (i == 0 && j == 0) {
        fg->Fidap_T[(fg->nr) * i + j] = fg->Tmax;
      }
    }                           /* i = 0, ...... */
  }                             /* j = 0, ...... */

  /*   return(0); */
}

/* end of T_var_read(FGrid_str *fg, time). */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file */
char const *rcs_id_T_var_read_c ()
{
  static char const rcsid[] = "$Id: T_var_read.c 1341 2008-07-23 15:23:30Z  $";

  return (rcsid);
}

/* end of rcs_id_fidap_interp_calc_c subroutine */

/*RCS Id:$Id: T_var_read.c 1341 2008-07-23 15:23:30Z  $*/
/*
*/
