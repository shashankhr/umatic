
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

/*RCS Id:$Id: fg_read.c 1402 2008-11-20 15:36:41Z  $*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include "safeopen.h"
#include "machine.h"
#include "constants.h"
#include "fidap.h"
#include "read_list/read_fg_list.h"
#include "read_list/readwrite_fg.h"
#include "read_list/freecsv.h"

/***********************************************************************/
/* subroutine to input fg temperature files.                           */
/*   Inputs: file name where the list is stored,                       */
/*            pointer to the FGrid location,                           */
/*           flag to indicate first/last/transient step                */
/*                                                                     */
/*   Output:  Non zero if the list is finished                         */
/*                                                                     */
/*   Result:  First time - calles read_listfile to store the list of time and filename */
/*                         Reads the first set of data                 */
/*           Other times - Read the next set fo data                   */
/*           Last time - free the locally allocated arrays             */
/**            \todo  decide how to handle non-transient (only one file) data  -- coupled */
/***********************************************************************/
int fg_read (const char *listfilename, FGrid_str * fg, int fg_flag)
{
  static FGrid_list_str fgl;
  static int n_names = 0, this_idx = 0;

  CA_FLOAT thistime = 0, nexttime = 0;

  int i;
  char thisname[255];

  fgl.nheaders_list = 1;

  switch (fg_flag) {
  case FG_FIRST_READ:
    n_names = read_listfile (listfilename, &fgl);
    /* fall through */
  case FG_TRANS_READ:
    sprintf (thisname, "%s.fgb", fgl.rows[this_idx]->filename);
    thistime = fgl.rows[this_idx]->time;

    if (n_names > this_idx + 1) {
      nexttime = fgl.rows[this_idx + 1]->time;
    } else {
      fprintf (stderr, "ERROR:read_fg_list: Ran out of file names! \n");
      fprintf (stderr, "ERROR:read_fg_list: n_names: %i this_idx: %i \n",n_names,this_idx);
      fprintf (stderr, "ERROR:read_fg_list: thisname: %s \n",thisname);
      fprintf (stderr, "Trying to finish ... \n");
      raise (SIGUSR1);
      exit (0);
    }
    break;

  case FG_CLEANUP:
    free_fg_list (&fgl);
    return (0);
    break;

  default:
    break;
  }
  read_fg_bin (thisname, fg, fg_flag);
  fg->tstart = thistime;
  fg->tnext = nexttime;
  this_idx++;
  if (this_idx >= n_names)
    return (1);
  return (0);

}

/*
*/
/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *fg_read_c ()
{
  static char const rcsid[] = "$Id: fg_read.c 1402 2008-11-20 15:36:41Z  $";

  return (rcsid);
}
