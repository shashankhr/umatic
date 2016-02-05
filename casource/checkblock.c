
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

#include <stdlib.h>
#include <stdio.h>
#include "machine.h"
#include "blocks.h"
#include "pore.h"

int check_pore (BB_struct * bp, PORE_str * porelist, int npores, int nt);
int check_sb (BB_struct * bp, int i);

/* realloc arrays to allow SSMALLOC_FASTCHK to catch heap corruption*/
/* use with caution since it could be very wasteful of memory! */
/* and cause everyone elses programs to grind to a halt, or crash */
int check_bb (BB_struct * bp)
{

  int i;
  int n, nco, nc;

  n = bp->ntsb;
  nc = bp->ncsb;
  nco = (bp->nc[0] + 2) * (bp->nc[1] + 2) * (bp->nc[2] + 2);

  for (i = 0; i < bp->ntsb; i++) {

    /*       fprintf(stderr,"Checking subblock %i\n",i); */
    check_sb (bp, i);
    if (realloc (bp->sb[i], sizeof (SB_struct)) == NULL)
      fprintf (stderr, "checkblock error 1");
  }

  realloc (bp->sb, n * sizeof (SB_struct *));
  realloc (bp->c_fs_values->block_array, n * sizeof (CA_FLOAT *));
  realloc (bp->c_fs_values, sizeof (Value_struct));

  realloc (bp->sch_fs_values->block_array, n * sizeof (CA_FLOAT *));
  realloc (bp->sch_fs_values, sizeof (Value_struct));
  realloc (bp->c_sol_values->block_array, n * sizeof (CA_FLOAT *));
  realloc (bp->c_sol_values, sizeof (Value_struct));
  realloc (bp->c_sol_alloy_values->block_array, n * sizeof (CA_FLOAT *));
  realloc (bp->c_sol_alloy_values, sizeof (Value_struct));

  realloc (bp->gr_array, n * sizeof (int *));

  realloc (bp->c_elm_array, n * sizeof (int *));

  realloc (bp->ftmp_one, nco * sizeof (CA_FLOAT));
  realloc (bp->ftmp_two, nco * sizeof (CA_FLOAT));
  realloc (bp->ftmp_three, nco * sizeof (CA_FLOAT));

#ifdef OLD_TUNDER
  realloc (bp->old_Tunder, nco * sizeof (CA_FLOAT));
#endif /*OLD_TUNDER */

  realloc (bp->itmp_one, nco * sizeof (int));

  /*check the grains */
  /*   fprintf(stderr,"Checking the grains\n"); */
  for (i = 0; i < bp->nprops.ngr; i++)
    realloc (bp->gr[i], sizeof (Ind_grain));

  realloc (bp->gr, bp->nprops.gd_max_total * sizeof (Ind_grain *));

  realloc (bp->ctrl->rgbp, sizeof (RGB_struct));

  return (0);
}
int check_sb (BB_struct * bp, int i)
{

  int nc, npores;
  SB_struct *sb;

  sb = bp->sb[i];
  nc = bp->ncsb;
  npores = sb->Npores;

  if (sb->open) {
    realloc (sb->gr, nc * sizeof (int));
    realloc (sb->c_elm, nc * sizeof (int));
    realloc (sb->c_fs, nc * sizeof (CA_FLOAT));
    realloc (sb->c_sol, nc * sizeof (CA_FLOAT));
    if (bp->ctrl->scheil)
      realloc (sb->sch_fs, nc * sizeof (CA_FLOAT));
    realloc (sb->c_sol_alloy, nc * sizeof (CA_FLOAT));

    /*check all the pores */
    check_pore (bp, sb->porelist, sb->Npores, bp->pprops.P_ntrad);

    realloc (sb->porelist, npores * sizeof (PORE_str));
  }

  return (0);

}
int check_pore (BB_struct * bp, PORE_str * porelist, int npores, int nt)
{
  int i, j;

  /*   fprintf(stderr,"Checking the pores\n"); */
  for (i = 0; i < npores; i++) {
    if ((porelist[i].State != PORE_NONE) && (porelist[i].State != PORE_LATENT)) {
      for (j = 0; j < N_T_LISTS; j++) {
        realloc (porelist[i].t_lists[j], nt * sizeof (CA_FLOAT));
      }
    }
    realloc (porelist[i].t_lists, N_T_LISTS * sizeof (CA_FLOAT *));
  }
  return (0);
}

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file */
char const *rcs_id_checkblock_c ()
{
  static char const rcsid[] = "$Id: checkblock.c 1339 2008-07-23 13:58:29Z  $";

  return (rcsid);
}

/* end of rcs_id_subroutine */
