
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

/*RCS Id:$Id: sb_diffuse_alloy_multi.c 1373 2008-08-27 20:51:52Z  $*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "blocks.h"
#include "matprops.h"
#include "multi_diff_props.h"
#include "machine.h"
#include "umat_matrix.h"
#include "sb_diffuse.h"

extern CA_FLOAT get_dl (CA_FLOAT temp);
extern CA_FLOAT get_ds (CA_FLOAT temp);
extern CA_FLOAT getav_d (CA_FLOAT dl, CA_FLOAT ds, CA_FLOAT fs);
extern int analythic (BB_struct * bp, int kk, int jj, int ii);
extern int interpolate (Ctrl_str * cp, BB_struct * bp, int cell_index, CA_FLOAT * cell_temp_in, int sbnum, int inter_flag,
                        CA_FLOAT cell_temperature);
extern int fcopy_surf (CA_FLOAT * big, CA_FLOAT * lit, Frame * cubeptr_local, int fcode, int type_local, int dir);

int sb_diffuse_alloy_multi (BB_struct * bp, int sbnum)
{

}                               /* end of sb_diffuse */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_sb_diffuse_alloy_multi_c ()
{
  static char const rcsid[] = "$Id: sb_diffuse_alloy_multi.c 1373 2008-08-27 20:51:52Z  $";

  return (rcsid);
}

/* end of rcs_id_sb_diffuse_alloy_c subroutine */
/*
*/
