
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

/*$Id: step_output.c 1356 2008-08-18 13:41:15Z  $*/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include "machine.h"
#include "blocks.h"
#include "writeblocks.h"
#include "interp.h"

extern void write_bin_blocks (BB_struct * bp);
extern CA_FLOAT cell_temp_calc_cc (BB_struct * bp, int sbnum, int x, int y);

/* write out a binary block data file if requested */
/* and any other output for the whole block */

void make_restart_file(){
   FILE * fp;
   fp = fopen("umat_step_restart.in", "w");
   fprintf(fp,"BlockRestartFileName step.%s",BL_EXT);
   fclose(fp);
}

void step_blk_out (BB_struct * bp, int step)
{
/* perform any file output for the current step */
  write_bin_blocks(bp);
   if (bp->ctrl->pore != 0) pore_write(bp);
  make_restart_file();

  /*                                                */
  /** \todo  fix the block output routines so they work multi block -- multiblock */
  /* or just make a new output stage ......         */
}                               /* end of step_blk_out */

          /*******************************************/
          /* perform screen print for current step   */
          /*******************************************/

void step_screen_out (BB_struct * bp, int step)
{
#ifdef _DEBUG_MALLOC_INC
  unsigned long hist0, size0;
#endif /*_DEBUG_MALLOC_INC*/
  Ctrl_str *cp = bp->ctrl;
  static CA_FLOAT oldfs = 0;
  static int oldtime = 0;       /*used for speed calc, old real time */
  CA_FLOAT Tout;                /*temp. to print to screen, average or con_cast */

  Tout = bp->sb[0]->Tvals.Tavg;
  /* print the elapsed time */
  fprintf (stderr, "\nElapsed calculation time: %i seconds\n", bp->realtime);
  fprintf (stderr, "\nASSUMING one subblock ONLY\n");
  fprintf (stderr, "Calculation constant: %.5e sec/cell/step\n                      %.5e this step\n",
           ((CA_FLOAT) bp->realtime) / ((CA_FLOAT) bp->ntsb * (CA_FLOAT) bp->ncsb * (CA_FLOAT) step),
           ((CA_FLOAT) bp->realtime -
            (CA_FLOAT) oldtime) / ((CA_FLOAT) bp->ntsb * (CA_FLOAT) bp->ncsb * (CA_FLOAT) (cp->scr_dmp_freq)));
  oldtime = bp->realtime;
#ifdef CLOCK
  {
    struct tms new_times;
    static CA_FLOAT old_time = 0;
    CA_FLOAT cputime = 0;
    CA_FLOAT cpu_const;

    times (&new_times);

    cputime = ((double) new_times.tms_utime - old_time) / ((double) (THE_CLOCK));
    old_time = (double) new_times.tms_utime;

    cpu_const = cputime / ((double) bp->ntsb * (double) bp->ncsb * (double) (cp->scr_dmp_freq));

    fprintf (stderr, "                      %.5e cpu clock\n", cpu_const);
    fprintf (stderr, "                      %.5e clock tick\n", (double) (THE_CLOCK));
  }
#endif /*CLOCK*/
#ifdef OLD_CLOCK
  {
    clock_t new_clocktime, clocktime;
    static clock_t old_clocktime = 0;
    CA_FLOAT cpu_const;

    new_clocktime = clock ();
    clocktime = new_clocktime - old_clocktime;
    old_clocktime = new_clocktime;
    cpu_const = (CA_FLOAT) clocktime / ((CA_FLOAT) bp->ntsb * (CA_FLOAT) bp->ncsb * (CA_FLOAT) (cp->scr_dmp_freq));
    fprintf (stderr, "                      %.5e cpu clock\n", cpu_const / 1e6);
  }
#endif /*CLOCK*/
#ifdef MIN_OUT
    fprintf (stderr, "STEP_OUTPUT: SSmin, SSmax not calculated for MIN_OUT mode\n");
#endif /*MIN_OUT */
  fprintf (stderr, "STEP_OUTPUT: ts: %d, bstep (total micro steps): %d sim_time: %4g \n", step, bp->step, bp->sim_time);
  fprintf (stderr, "STEP_OUTPUT: Tavg:%10f, fs:%4g, nsol:%i, npor:%i, ngr: %d \n",
           Tout, bp->fs_active, bp->sb[0]->ncsolid, bp->bb_npores, bp->nprops.ngr);
  fprintf (stderr, "STEP_OUTPUT: Tmin:%10f, Tmax %.10f\n", bp->sb[0]->Tvals.Tmin, bp->sb[0]->Tvals.Tmax);
#ifdef BUG_MAY_03
  if (bp->fs_active < oldfs) {
    int dumb;

    dumb = 0;
  }
  oldfs = bp->fs_active;

#endif

#ifndef MIN_OUT
  fprintf (stderr, "STEP_OUTPUT: SATmax[%.4g], SSmin[%.4g], SSmax[%.4g], n_active_sb[%i]\n",
           bp->c_sol_values->SATmax, bp->c_sol_values->SSmin, bp->c_sol_values->SSmax, bp->n_a_sb);
#endif /*MIN_OUT */
  fprintf (stderr, "STEP_OUTPUT: W[%i], C[%i], P[%i]\n", bp->dfs_warn, bp->dfs_cap, bp->pore_err);
#ifdef _DEBUG_MALLOC_INC
  size0 = malloc_inuse (&hist0);
  fprintf (stderr, "STEP_OUTPUT: DB_MALLOC in use : %i\n", size0);
#endif /*_DEBUG_MALLOC_INC*/
}                               /* end of step_screen_out */

          /*******************************************/
          /* perform file output  for current step   */
          /*******************************************/
void step_file_out (BB_struct * bp, int step)
{
  Ctrl_str *cp = bp->ctrl;
  int i;
  static int chkflag = 0;

  if (chkflag)
    return;

  for (i = 0; i < cp->nsbslice; i++) {
    if (cp->slice[i][0] >= bp->ntsb || cp->slice[i][1] >= bp->nc[2]) {
      chkflag = 1;
      fprintf (stderr, "ERROR: step_file_out: chosen picture is out of range. sb %i slice %i \n", cp->slice[i][0], cp->slice[i][1]);
      fprintf (stderr, "ERROR: step_file_out: Deactivating picture output! \n");
      return;
    }

    if (bp->nprops.ngr > 0) {
      write_slice (bp, cp->slice[i][0], cp->slice[i][1]);
    }
#ifndef CYGWIN
    if (cp->floatdump) {
      if (cp->diffuse == TRUE) {
        write_slice_conc (bp, bp->c_sol_values, cp->slice[i][0], cp->slice[i][1]);
        write_slice_sat (bp, cp->slice[i][0], cp->slice[i][1]);

      }
      if (cp->diffuse_alloy == TRUE || cp->particle == TRUE) {
        write_slice_conc (bp, bp->c_sol_alloy_values, cp->slice[i][0], cp->slice[i][1]);
        write_slice_undercool (bp, cp->slice[i][0], cp->slice[i][1]);
      }
    }
#endif /*CYGWIN*/
      if (cp->floatdump) {
      write_shortint_slice (bp, cp->slice[i][0], cp->slice[i][1]);
      write_slice_fs (bp, cp->slice[i][0], cp->slice[i][1]);
/**  \todo  only write this if the option is selected -- merge-xly */

      if (cp->decentred_octahedron == TRUE) {
        write_slice_curv (bp, cp->slice[i][0], cp->slice[i][1]);
      }
      /* write_slice_diff(bp,cp->slice[i][0],cp->slice[i][1]); */
      /* write_slice_sat(bp,cp->slice[i][0],cp->slice[i][1]); */
    }
    if (cp->tempslice)
      write_sb_slice_temp (bp, cp->slice[i][0], cp->slice[i][1]);
  }

  for (i = 0; i < cp->nbbslice; i++) {
    if (cp->tempslice)
      write_bb_T_slice (bp, cp->bbslice[i]);
  }

  if (bp->nprops.ngr > 0) {
    for (i = 0; i < cp->nbbslice; i++) {
#ifndef CYGWIN
      write_bb_slice (bp, cp->bbslice[i]);
#endif /*CYGWIN*/
    }
  }
}                               /* end of current step file output */

/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *step_output_c ()
{
  static char const rcsid[] = "$Id: step_output.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}
