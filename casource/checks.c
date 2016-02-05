
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

/*RCS ID: $Id: checks.c 1339 2008-07-23 13:58:29Z  $*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "machine.h"
#include "blocks.h"

int find_ngrow (BB_struct * bp, int curr)
{
  SB_struct *sp;
  CA_FLOAT *fs, *fsend;
  int ngrow = 0;

  sp = bp->sb[curr];
  fsend = sp->c_fs + bp->ncsb;
  for (fs = sp->c_fs; fs < fsend; fs++) {
    ngrow += ((*fs > 0) && (*fs < 1));
  }
  return (ngrow);
}

CA_FLOAT find_fs_av (BB_struct * bp, int curr)
{
  SB_struct *sp;
  CA_FLOAT *fs, *fsend, sumfs = 0, fsav;

  sp = bp->sb[curr];
  fsend = sp->c_fs + bp->ncsb;
  for (fs = sp->c_fs; fs < fsend; fs++) {
    sumfs += *fs;
  }
  fsav = sumfs / (bp->ncsb - sp->nmould);
  return (fsav);
}

/* find the av. number of grwoing/solid cells in a grain */
CA_FLOAT grain_ncells (BB_struct * bp)
{
  int i;
  int ncells = 0;
  CA_FLOAT ncav = 0;

  if (bp->nprops.ngr == 0)
    return (0);
  for (i = 1; i < bp->nprops.ngr; i++) {
    ncells += bp->gr[i]->ncells;
  }
  ncav = (CA_FLOAT) ncells / bp->nprops.ngr;
  return (ncav);
}                               /*end of grain_ncells */

CA_FLOAT grain_ngrow (BB_struct * bp)
{
  int i;
  int ncells = 0;
  CA_FLOAT ncav = 0;

  if (bp->nprops.ngr == 0)
    return (0);
  for (i = 1; i < bp->nprops.ngr; i++) {
    ncells += bp->gr[i]->ngrow;
  }
  ncav = (CA_FLOAT) ncells / bp->nprops.ngr;
  return (ncav);
}                               /*end of grain_ngrow */

/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *checks_c ()
{
  static char const rcsid[] = "$Id: checks.c 1339 2008-07-23 13:58:29Z  $";

  return (rcsid);
}
