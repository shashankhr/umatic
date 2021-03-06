
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

/*
 Main subroutine for a portable Cellular Automata code.	
 It uses a multi-block approach to solve large problems, and 
 can either be run by itself (using CA_WRAPPER),              
 as a postprocessor to fem programmes, or                     
 as a user subroutine within a fem programme (e.g. with ProCast)  
 The superblock is called: bigblock, and is divided into      
 subblocks.                                                   
 The current version simulates grain nucleation and growth    
 for eutectic equiaxed growth.                               
 Written by Peter D. Lee, Imperial College    		
 Wed Jul  1 18:38:31 bst 1998                 		
*/

/*RCS ID:$Id: ca_solid.c 1405 2008-12-04 13:47:23Z  $*/
/* Includes and defines.*/
/* include system headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <signal.h>

/* include header files requred by subroutines */
#include "machine.h"
#include "read_ctrl.h"
#include "blocks.h"
#include "readfiles.h"
#include "writeblocks.h"
#include "interp.h"
#include "find_max.h"
/* headers needed for the package for reading finite-element grid data */
/* eg from Sandia code */
#include "read_list/read_fg_list.h"
#include "read_list/readwrite_fg.h"

/* different temperature calculation option funcitons */
#include "temp_calc.h"
int init_output_img (RGB_struct * rgbp);

extern int pore_setup (BB_struct * bp, int sbnum);

/* Function Prototypes.*/
extern int init_output_img (RGB_struct * rgbp);
extern int pore_setup (BB_struct * bp, int sbnum);
extern void free_fg (FGrid_str * fg);

#ifdef TEST_CLOSE
extern int close_sb (BB_struct * bp, int sbnum);
#endif /*TEST_CLOSE */
extern void write_bin_blocks (BB_struct * bp);
extern void read_bin_blocks (BB_struct * bp, const char *fname);

extern int free_bb (BB_struct * bp);
extern int sb_boundary (BB_struct * bp, int sbnum);

/* functions used from bigblock.c */
extern int init_bb (BB_struct * bp, Ctrl_str * cp);

#ifdef EXTERNAL_CA
#ifdef VERBOSE_EXTERNAL
/* access to Procast variables, only for debugging */
/* in normal mode no access is needed in this routine */
#include "common.h"
extern float *t0;
extern float *t1;
#endif /* VERBOSE_EXTERNAL */
#endif /*EXTERNAL_CA*/

/* functions used from subblock.c */
extern int init_sb (BB_struct *, int);
extern int calc_sb (BB_struct *, int);
extern int open_sb (BB_struct *, int);

/* functions used from sb\_nuc.c */
extern int sb_nuc_area (BB_struct * bp, int sbnum, NucArea_struct * nap);

/* functions used from init\_sb\_neigh.c */
extern int init_sb_neigh (BB_struct *, int);

/* from sb_temp_calc.c */
extern void cell_temp_calc_extern (BB_struct * bp, int index_ca);
extern void find_cell_pressure (BB_struct * bp, int index_ca);
extern float min_temp_calc_extern (BB_struct * bp);

/* from step_output.c */
extern void step_screen_out (BB_struct * bp, int step);
extern void step_file_out (BB_struct * bp, int step);
extern void step_blk_out (BB_struct * bp, int step);

/* from umat_wrapper.c -- to test the signal */
extern int jflg;
extern int signal_change_freq;

/* from fg_read.c */
extern int fg_read (const char *listfilename, FGrid_str * fg, int fg_flag);

/* to set up the chosen mould boundary condition functions */
extern void setup_mould_src_function (Ctrl_str * cp, BB_struct * bp, Solute_props * sp);
extern void setup_temp_function (Ctrl_str * cp, BB_struct * bp);

/* Routines to do any housekeeping necessary to allow control options to be */
/* changed upon restarting */

void restart_change_options (Ctrl_str * cp, BB_struct * bp)
{
  int sbnum = 0, i;
  int errors = 0;
  SB_struct *sp;
  CA_FLOAT k_gas, c_gas_liq;

  k_gas = bp->mprops.gasprops.part_coef[0];

  for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {
    sp = bp->sb[sbnum];
    if (sp->open != SB_OPEN)
      continue;
         /**********************************************/
         /**    Turn gas diffusion on when restarting **/
         /**********************************************/
    /* if we are turning on the gas calcualtion after a restart */
    /* we need to allocate and set the arrays */
    if (cp->restart_gas_on) {

      fprintf (stderr, "RESTARTING: turning on the gas diffusion.\n");
      sp->Svals[GAS].Cinit = bp->mprops.gasprops.Cinit; /*all the same for now! */
      /* lever rule partitioning */
      c_gas_liq = bp->mprops.gasprops.Cinit / (1 - (1 - k_gas) * sp->Tvals.fsavg);

      if (!(sp->c_sol = (CA_FLOAT *) calloc (bp->ncsb, sizeof (CA_FLOAT)))) {
        fprintf (stderr, "ERROR:restart_change_options: SB cell sol array malloc failed\n");
        exit (423);
      }
      /* set concentraion to initial value */
      /* according to the fraction solid in the cell */
      fprintf (stderr, "RESTARTING: setting initial gas concetnration to equilibrium partitioning.\n");
      for (i = 0; i < bp->ncsb; i++)
        *(sp->c_sol + i) = (1 - (1 - k_gas) * *(sp->c_fs + i)) * c_gas_liq;
      /* set value pointer array */
      bp->c_sol_values->block_array[sbnum] = sp->c_sol;
    }

    /* end restart gas */
 /************************************/
    /* Turn porosity on when restarting */
 /************************************/
    /* if we are turning on the pore calcualtion after a restart */
    /* we need to allocate and set the arrays,initialize threshold */
    if (cp->restart_pore_on) {
      errors += pore_setup (bp, sbnum);
      fprintf (stderr, "RESTARTING: turning on the pore calculation.\n");
      if (errors != 0) {
        fprintf (stderr, "ERROR:restart_change_options: pore setup failed\n");
        exit (errors);
      }
    }                           /* end restart pore */
  }                             /* end loop through subblocks */

    /*******************/
  /* reset the flags */
  /* if it is stopped and then restarted again, we do not want */
  /* to reallocate the arrays! */
    /*****************************/
  cp->restart_gas_on = 0;
  cp->restart_pore_on = 0;

}

/* decide how often to reload a finite-element grid */
void fg_find_steps (BB_struct * bp, FGrid_str * fg) {
  CA_FLOAT timedif;
  int stepdif;

  timedif = fg->tnext - fg->tstart;
  stepdif = (int) (FLOOR (timedif / bp->delt));
  fg->s_next = fg->s_start + stepdif;
  fg->stepdif = stepdif;
  fg->wf_add = (double) 1 / ((double) (stepdif));
  if (ABS (fg->wf_add) <= MINVAL) {
    fprintf (stderr, "ERROR: fg_find_steps: Incorrect wf_add or not enough precision %.5g\n", fg->wf_add);
    exit (0);
  }
}

/****************************************************************/
/* Beginning of the umat_solid program!				*/
/****************************************************************/
/* Main part.*/
/* Main part*/
/* Main part.*/

/**
 * umat_solid: function to
 * run a full micro step consisting of several micro timesteps 
 * or an entire stand-alone simulation.
 * @callgraph
 * @callergraph
*/
int umat_solid (int stat_flag, CA_FLOAT i_time, CA_FLOAT delt, Ctrl_str * cp, BB_struct * bp) {
  /* declare counters and output variables for main */
  int i, j, k, ii, sbnum;
  int fsprintflg = 0, tcprintflg = 0, tsprintflg = 0, tbprintflg = 0;   /* flags for testing output condition */
  int new_fg = FALSE;
  int h;
  static int fspchk[MAX_CTRL] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  int err = 0;
  CA_FLOAT min_temp, start_temp;
  int step, pr_step, olddmp, oldscrdmp;
  int done = FALSE;
  CA_FLOAT ptime, del_ptime = 1.0;

  /* declare the main CA structure that holds all the rest */
  time_t thetime;
  sigset_t sigmask, newsig;

#ifdef  VERBOSE_EXTERNAL
  /* only used for exhaustive debugging of external */
  /* by dumping out all node numbers etc */
  int node_number;
  FILE *test6, *test7;
#endif /*VERBOSE_EXTERNAL */
#ifdef _DEBUG_MALLOC_INC
  unsigned long oldsize;
  unsigned long size1, size2;
  unsigned long hist0, hist1, hist2;

  oldsize = malloc_inuse (&hist0);
#endif /*_DEBUG_MALLOC_INC*/

  bp->n_a_sb = 0;

#ifdef  VERBOSE_EXTERNAL
  test6 = fopen ("external.node", "w+");
  test7 = fopen ("external.step", "w+");
#endif /* VERBOSE_EXTERNAL */

/*   test11=fopen("external.current","w+"); */

/****************************************************************/
/* First check the stat\_flag to determine the operation mode.   */
/****************************************************************/
  switch (stat_flag) {
      /*******************************************/
    /* Restart: read in the bin block file     */
      /*******************************************/
  case RESTART_BB:
    {
      /* NOTE at this point cp and bp->ctrl are NOT the same */
      /* cp is the dummy ctrl structure for restart */
      /* whereas bp->ctrl is the structure that was saved in */
      /* the big block file */

      read_bin_blocks (bp, cp->fn_block_restart);
      bp->ctrl->restart = 1;    /* pass a flag for restart options */

      /* read any values that need to be changed */
      read_ctrl_vals (cp->fn_ctrl, bp->ctrl);
      read_matprop_vals (bp->ctrl, bp);
      read_geoplus_vals (bp->ctrl, bp);
      setup_temp_function (bp->ctrl, bp);
      setup_mould_src_function (bp->ctrl, bp, &(bp->mprops.gasprops));

      for(i=0;i<bp->mprops.ele_1;i++){
         setup_mould_src_function (bp->ctrl, bp, &(bp->mprops.alloyprops[i]));
      }

      /** \todo  get rgbp out of control structure  -- general */
      /*should be ONLY used in output image routines */
      init_output_img (bp->ctrl->rgbp);

      /* apply any re-initialization in order to change options */
      /* NOT ALL OPTIONS MAY BE MEANINGFULLY CHANGED until      */
      /* a routine is written to set up the changes needed.     */
      restart_change_options (bp->ctrl, bp);
      /* reset the real time counter */
      bp->starttime = time (0);
      bp->realtime = 0;

    }

#ifndef MIN_OUT
    if (bp->ctrl->excel)
      write_bb_excel (RESTART_CALL, bp);
#endif /* MIN_OUT */
    break;
/*******************************************/
/* end of RESTART_BB                       */
/*******************************************/

/****************************************************************/
/* Initialise all: move most of this into subroutines later...  */
/****************************************************************/
  case INIT_BB:                /* first call to umat_solid, do all initilisation */
    /* not restarting, the initial control structure may be used */
    /* as the real control structure */
    bp->ctrl = cp;
    if (init_bb (bp, cp) != 0) {
      fprintf (stderr, "exiting due to init_bb failure");
      exit (2);
    }
    /* fallthrough on this case, EXTERNAL initialized the big block */
    /* when generating the cells */
  case EXTERNAL_INIT_BB:
    bp->starttime = time (0);
    bp->realtime = 0;
    fprintf (stderr, "Starting to create and initialise subblocks \n");
    fprintf (stderr, "ntsb %d nsb: %d, %d, %d\n\n", bp->ntsb, bp->nsb[0], bp->nsb[1], bp->nsb[2]);
    /* loop through all required subblock for now creating them */
    fprintf (stderr, "umat_solid: B...\n");
    setup_temp_function (cp, bp);

    for(i=0;i<bp->mprops.ele_1;i++){
       setup_mould_src_function (bp->ctrl, bp, &(bp->mprops.alloyprops[i]));
    }

    setup_mould_src_function (cp, bp, &(bp->mprops.gasprops));

    /* macromodel from file --set up procedure once per run */
    if (cp->fgrid_input == TRUE) {
      /*******************************************/
      /* If the FE grid info is used, then       */
      /* Initialise the R direction interpolation */
      /* factors -- only once per run            */
      /* assuming that the R grid does not change */
      /*******************************************/
      fprintf (stderr, "Reading in fgrid structure data from %s ... \n", cp->fn_fgrid);
      /* read the first and second finit grids to set up the interpolation */
      fg_read (cp->fn_fgrid, bp->fg, FG_FIRST_READ);
      bp->fg->state = cp->fg_tr;
      fg_read (cp->fn_fgrid, bp->fg_next, FG_TRANS_READ);
      bp->fg_next->state = cp->fg_tr;
      /* set up the step to read the next file */
      bp->fg->s_start = bp->step;
      fg_find_steps (bp, bp->fg);
      bp->fg_next->s_start = bp->fg->s_next;
      fg_find_steps (bp, bp->fg_next);

      bp->fg->wfact = 0;

    }

    /* create the offset for the subblock neighbours */
    init_sb_neigh (bp, WRAP);

         /*****perform the temp and pres***/
         /*****allocation for each cell*************************/
    if (cp->external) {
      for (h = 0; h < bp->total_cell_number; h++) {
        cell_temp_calc_extern (bp, h);
        if (bp->ctrl->flow_on == TRUE) {
          find_cell_pressure (bp, h);
        }
      }
      bp->first_micro = 0;
    }
         /**********************************************/

    start_temp = bp->mprops.Tstart;
    if (FALSE) {                /** \todo  test subblocks for start condition -- multiblock */
      for (i = 0; i < bp->ntsb; i++) {
        min_temp = bp->sb[i]->Tvals.Tmin;
        if ((min_temp <= start_temp) && (bp->sb[i]->open == SB_NEW)) {
          if (open_sb (bp, i) != 0) {
            fprintf (stderr, "exiting due to open_sb failure");
            exit (2);
          }
          /* Check if there are nucleation areas, and set them! */
          for (j = 0; j < bp->nprops.nareanuc; j++) {
            if (bp->nprops.nap[j].sbnum == i)
              sb_nuc_area (bp, bp->nprops.nap[j].sbnum, &(bp->nprops.nap[j]));
          }
        }                       /*end start_temp comparison */
      }                         /*end loop through subblocks */

    } else {                    /* Non-concast: just open them all for now */
      for (i = 0; i < bp->ntsb; i++) {
        if (open_sb (bp, i) != 0) {
          fprintf (stderr, "exiting due to open_sb failure");
          exit (2);
        }
      }

      /* Check if there are nucleation areas, and set them! */
      for (i = 0; i < bp->nprops.nareanuc; i++) {
        sb_nuc_area (bp, bp->nprops.nap[i].sbnum, &(bp->nprops.nap[i]));
      }
    }                           /* end of open initial subblocks */

#ifndef MIN_OUT
    if (cp->excel)
      write_bb_excel (FIRST_CALL, bp);
#endif /* MIN_OUT */

    fprintf (stderr, "CA_SOLID: Finished creating&initialising subblocks. \n");
    break;

/*******************************************/
/* End of INIT_BB                          */
/*******************************************/

/****************************************************************/
/* CALC_BB: Perform one timestep calculation.                   */
/****************************************************************/
  case CALC_BB:                /* perform a series of sub-timesteps */
    cp = bp->ctrl;              /**AARGH  \todo  straighten this out! -- general */
       /************************************************/
    /* CASE SOLO mode:                              */
    /*    If in SOLO mode, loop through the entire  */
    /*    calc. without returning to ca\_wrapper.    */
    /* CASE CAP mode:                               */
    /*    If in CAP mode, loop through sub-steps    */
    /*    and then return for another cap timestep  */
       /************************************************/
    /* loop through time calling the ca routine...  */
       /************************************************/

    fprintf (stderr, "CA_SOLID: starting subblock calculation CALC_BB \n");
    fprintf (stderr, "with flag (done) = %i  \n",done);
    fprintf (stderr, "CA_SOLID: sim,fin,delt,nsteps,bstep: %f, %f, %f, %i, %i\n", bp->sim_time, bp->finish_time, bp->delt, bp->nsteps,
             bp->step);

       /********at the begining of each macro time step reset the counter**/

    bp->first_micro = 0;
    step = 0;
    start_temp = bp->mprops.Tstart;

    /* loop through all required subblocks performing required steps */
    while (!done) {
         /**********************************************************************/
      /* mask the USR1 and USR2 signal, so that requested out put will be consistant */
      /* note you can still force ending with TERM signal ...               */
         /**********************************************************************/
      /* the ALL_SIGS preprocesser definition will prevent the signal masking */
      /* this is useful for debugging */

#ifndef ALL_SIGS
      sigaddset (&newsig, SIGUSR1);
      sigaddset (&newsig, SIGUSR2);
      sigprocmask (SIG_BLOCK, &newsig, &sigmask);
#endif
         /**********************************************************************/

             /*****perform for the current step the temp and pres***/
             /*****allocation for each cell*************************/
      if (cp->external) {
        for (h = 0; h < bp->total_cell_number; h++) {
          cell_temp_calc_extern (bp, h);
          if (bp->ctrl->flow_on == TRUE) {
            find_cell_pressure (bp, h);
          }
        }
      }
             /**********************************************/

      bp->sim_time += bp->delt;
      thetime = time (0);
      bp->realtime = (int) (thetime) - (int) (bp->starttime);
      fsprintflg = 0;
      new_fg = FALSE;
      step++;
      bp->micro_step = step;
      /**********************************************/
      /* macromodel from file: set up per macro timestep */
      /* -- first test if this microstep is a new macrostep */
      /* Read in the temperature interpolation file */
      /* once per timestep unless steady state      */
       /**********************************************/
      if (cp->fgrid_input == TRUE) {
        if (bp->fg->state == TRANSIENT) {
          /* flip the next fg into the current one */
          /* get rid of the old fg data */
          if (bp->step == bp->fg->s_next) {
            new_fg = TRUE;
            free_fg (bp->fg);
            memcpy (bp->fg, bp->fg_next, sizeof (FGrid_str));
            fg_read (cp->fn_fgrid, bp->fg_next, FG_TRANS_READ);
            bp->fg_next->s_start = bp->fg->s_next;
            bp->fg_next->state = bp->fg->state;
            fg_find_steps (bp, bp->fg_next);
            bp->fg->wfact = 0;
          } else {
            /* update the linear weight factor for time interpolation */
            bp->fg->wfact += bp->fg->wf_add;
          }
        }
      }
      bp->n_a_sb = 0;
      bp->fs_active = 0;
      for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {
        /* loop trhough sb's checking if open */
        if (bp->sb[sbnum]->open == FALSE) {
          /* check if sb should be opened */
          /** \todo  figure out how to set this for the finite grid method -- general*/
          min_temp = bp->sb[sbnum]->Tvals.Tmin;
          if (min_temp <= start_temp) {
            if (open_sb (bp, sbnum) != 0) {
              fprintf (stderr, "CA_SOLID: Exiting due to open_sb %i failure\n", sbnum);
              exit (2);
            }
            /* open sb */
            /* Check if there are nucleation areas, and set them! */
            for (j = 0; j < bp->nprops.nareanuc; j++) {
              if (bp->nprops.nap[j].sbnum == i)
                sb_nuc_area (bp, bp->nprops.nap[j].sbnum, &(bp->nprops.nap[j]));
            }
          }                     /*end start_temp comparison */
        }

        /* end if sb not open */
 /************************/
        /*  Call the |calc_sb|  */
 /************************/
        else if (bp->sb[sbnum]->open == TRUE) { /* if open, calc it */
          /* if a new fg has been read in, update the base temeprature arrays */
          /* don't need to test fgrid_input cause new_fg is always false otherwise */
          if (new_fg) {
            fg_temp_upd (bp, sbnum);
          }
          calc_sb (bp, sbnum);

          /* track the fraction solid in active sb's */
          bp->n_a_sb++;
          bp->fs_active += bp->sb[sbnum]->Tvals.fsavg;

        }
             /************************/
        /* end calc sb          */
             /************************/
      }                         /* end loop through subblocks */

      if (bp->n_a_sb > 0)
        bp->fs_active /= bp->n_a_sb;
      else
        bp->fs_active = 0;

      /* if no sub block is open return to the higher level routine */
      if (bp->n_a_sb == 0) {

               /************************************************/
        /* close all open input/output files and exit   */
               /************************************************/
#ifdef  VERBOSE_EXTERNAL
        fclose (test6);
        fclose (test7);
#endif /*VERBOSE_EXTERNAL */
        return (1);
      }
         /**********************************************************************/
      /* unmask the USR1,USR2 signal                                             */
         /**********************************************************************/
#ifndef ALL_SIGS
      sigprocmask (SIG_UNBLOCK, &newsig, &sigmask);
#endif
         /**********************************************************************/
      /* if the USR1 signal is received, it waits until here to process it. */
         /**********************************************************************/

          /*******************************************/
      /* Check for end conditions                */
          /*******************************************/
      if (cp->solo && bp->auto_fin) {
        if (bp->sim_time > bp->autofin_time) {  /* failsafe stop of programme... */
          done = TRUE;
          fprintf (stderr, "CA_SOLID: Reached autofin_time [%g].\n", bp->autofin_time);
        }
      } else if (bp->sim_time > bp->finish_time) {      /* failsafe stop of programme... */
        done = TRUE;
        fprintf (stderr, "CA_SOLID: Reached finish_time.\n");
      }
      /*stop at specified fraction solid */
      /*needs update for MULTIBLOCK */
      if (bp->sb[0]->Tvals.fsavg >= cp->fs_finish) {
        done = TRUE;
        if ((cp->external) && (bp->last_micro == 0)) {
          bp->last_micro = 1;
        }
        fprintf (stderr, "CA_SOLID: Reached fs_finish fraction solid %g.\n", bp->sb[0]->Tvals.fsavg);
      }
#ifdef REALTIME_LIMIT
      /* stop at limit of real time */
      if (bp->realtime > bp->realtime_limit) {
        done = TRUE;
        fprintf (stderr, "CA_SOLID: Reached elapsed realtime limit %i.\n", bp->realtime);
      }
#endif

          /*******************************************************/
      /* check if fraction solid has reached a print stage   */
          /*******************************************************/
      for (j = 0; j < cp->nfsprint; j++) {
        if ((bp->sb[0]->Tvals.fsavg >= cp->fsprint[j]) && (fspchk[j] == 0)) {
          fspchk[j] = 1;
          fsprintflg = 1;
          fprintf (stderr, "CA_SOLID: Reached fsprint fraction solid %g.\n", bp->sb[0]->Tvals.fsavg);
        }
      }                         /* end check fsprint */
      if (cp->time_dump) {
        tcprintflg = 0;
        tsprintflg = 0;
        tbprintflg = 0;
        if (bp->sim_time >= ((bp->scr_dump_num + 1) * (cp->time_unit * cp->scr_dmp_freq))) {
          bp->scr_dump_num++;
          tcprintflg = 1;
        }
        if (bp->sim_time >= ((bp->slice_dump_num + 1) * (cp->time_unit * cp->slice_dmp_freq))) {
          bp->slice_dump_num++;
          tsprintflg = 1;
        }
        if (bp->sim_time >= ((bp->blk_dump_num + 1) * (cp->time_unit * cp->blk_dmp_freq))) {
          bp->blk_dump_num++;
          tbprintflg = 1;
        }
      }

          /******************************************/
      /* check if signal USR2 has been received */
          /******************************************/
      if (jflg == JFLG_WRITEBIN) {
          /************************************************/
        /* write the entire bb to a file              */
          /************************************************/
        fprintf (stderr, "writing binary bigblock, bstep %i\n", bp->step);
        write_bin_blocks (bp);
        jflg = 0;
      }

          /*******************************************/
      /* perform screen print for current step   */
          /*******************************************/
      /* protect from % operation if zero */
      if (cp->scr_dmp_freq != 0) {
        if ((bp->step % cp->scr_dmp_freq == 0 && cp->time_dump == 0) || (fsprintflg == 1) || (tcprintflg == 1)) {
          step_screen_out (bp, step);
        }                       /* end of screen dump printout */
      }

      /* protect from % operation if zero */
      if (cp->blk_dmp_freq != 0) {
        if ((bp->step % cp->blk_dmp_freq == 0 && cp->time_dump == 0) || (fsprintflg == 1) || (tbprintflg == 1)) {
          step_blk_out (bp, step);
        }                       /* end of block dump */
      }

      /* not minimum output mode */
#ifndef MIN_OUT
          /************************************************/
      /* perform any file output for the current step */
          /************************************************/
      /* protect from % operation if zero */
      if (cp->slice_dmp_freq != 0) {
        if (((cp->nsbslice > 0)
             || (cp->nbbslice > 0))
            && bp->step != 0 && (((bp->step % cp->slice_dmp_freq == 0) && (cp->time_dump == 0))
                                 || (fsprintflg == 1) || (tsprintflg == 1))) {
#ifdef IMAGEOUT
          step_file_out (bp, step);
#else
          fprintf (stderr, "CA_SOLID: image output deactivated!\n");
#endif
        }                       /* end of current step file output */
      }
#endif /*MIN_OUT */

#ifdef VERBOSE_EXTERNAL
/* external feature, not sure why needed? */
      for (j = 0; j < bp->total_cell_number; j++) {
        if (bp->cell_element_array[j] >= 0) {
          node_number = bp->cell_node_array[j];
          fprintf (test6, "%d \t %d \t %.12f \t %.12f \t %f \t %f \n", j, node_number, bp->cell_temp_extern[j] - 273.16,
                   bp->cell_temp_change_extern[j], t0[node_number] - 273.16, t1[node_number] - 273.16);
          break;
        }
      }
#endif /* VERBOSE_EXTERNAL */
#ifndef MIN_OUT
      if ((cp->excel) && ((step % cp->scr_dmp_freq == 0) || (fsprintflg == 1)))
        write_bb_excel (STEP_CALL, bp);
#endif /*MIN_OUT */

#ifdef MANYDUMPS
          /********************************************/
      /* Debug mode to                            */
      /* zero-in on interesting behaviour         */
      /* start producing output at defined higher */
      /* frequency                                */
          /********************************************/
      if (bp->step == MANY_DUMP_START) {
        olddmp = cp->slice_dmp_freq;
        oldscrdmp = cp->scr_dmp_freq;
        cp->slice_dmp_freq = MANY_FREQ;
        cp->scr_dmp_freq = MANY_FREQ;
      }
      if (bp->step == MANY_DUMP_STOP) {
        cp->slice_dmp_freq = olddmp;
        cp->scr_dmp_freq = oldscrdmp;
      }
#endif /*MANYDUMPS*/
  /*********************************************************************/
      /* apply any requested real-time change */
      if (signal_change_freq != 0 ){
          printf("%s: Output step change requested.\n",__func__);
          printf("New step: singal_change_freq= %i\n",signal_change_freq);
          cp->blk_dmp_freq = signal_change_freq;
          signal_change_freq=0;
      }
        bp->first_micro += 1;
        bp->step += 1;

#ifdef  VERBOSE_EXTERNAL
  /*********output time step information******/
      fprintf (test7, "step=%d \t bp->step=%d \t bp->first_micro=%d \t bp->sim_time=%f \t bp->finish_time=%f \t bp->nsteps=%d \n",
               step, bp->step, bp->first_micro, bp->sim_time, bp->finish_time, bp->nsteps);
  /***********************************************/
#endif /*VERBOSE_EXTERNAL */
    }

    fprintf (stderr, "End of CALC_BB: flag (done) = %i  \n",done);
    break;

/*******************************************/
/* End of CALC_BB                          */
/*******************************************/

/****************************************************************/
/* FINISH_BB: Print out final results.                          */
/****************************************************************/
  case FINISH_BB:              /* finished all calc, preform final output and free */

#ifdef MIN_OUT
#else /*not MIN_OUT */
          /************************************************/
    /* write the entire bb to a file              */
          /************************************************/
#ifdef  VERBOSE_EXTERNAL
    fclose (test6);
    fclose (test7);
#endif /*VERBOSE_EXTERNAL */
    /* return(1); *//*??? not sure why this is here */
    fprintf (stderr, "writing binary bigblock, bstep %i\n", bp->step);
    write_bin_blocks (bp);

    step_file_out (bp, step);

    if (cp->excel) {
      write_bb_excel (LAST_CALL, bp);
      write_bb_excel_sum (LAST_CALL, bp);
    }
#endif /*not MIN_OUT */
    fprintf (stderr, "DFS Warnings: %i DFS Capped: %i\n", bp->dfs_warn, bp->dfs_cap);
/**********************/
/* free data pointers */
/**********************/
#ifdef TEST_CLOSE
    close_sb (bp, 0);
#endif /*TEST_CLOSE */
#ifdef _DEBUG_MALLOC_INC
    size1 = malloc_inuse (&hist1);
    fprintf (stderr, "umat_solid list before free_bb %i %i\n", hist0, hist1);
    malloc_list (2, bp->ctrl->hist[0], hist1);
#endif       /*_DEBUG_MALLOC_INC*/

    fprintf (stderr, "Freeing the Big Block\n");

    if (cp->fgrid_input) {
      fg_read (cp->fn_fgrid, bp->fg, FG_CLEANUP);
    }

    free_bb (bp);
#ifdef _DEBUG_MALLOC_INC
    size1 = malloc_inuse (&hist1);
    fprintf (stderr, "umat_solid list after free_bb %i %i\n", hist0, hist1);
    malloc_list (2, bp->ctrl->hist[0], hist1);
#endif       /*_DEBUG_MALLOC_INC*/
    break;

  default:
    fprintf (stderr, "ERROR:umat_solid: stat_flag value [%d] undefined.\n", stat_flag);
    break;
  }                             /* end of stat_flag switch */

/************************************************/
/* close all open input/output files and exit   */
/************************************************/
#ifdef _DEBUG_MALLOC_INC
  fprintf (stderr, "umat_solid list %i %i\n", hist0, hist1);
  malloc_list (2, hist0, hist1);
#endif /*_DEBUG_MALLOC_INC*/
  return (0);
}                               /* end of umat_solid */

/* rcs id subroutine.*/
/*Little subroutine to include the |rcs Id| in the object code.*/
/*This can also be called to print or include the Id in a file.*/

char const *rcs_id_umat_solid_c ()
{
  static char const rcsid[] = "$Id: ca_solid.c 1405 2008-12-04 13:47:23Z  $";

  return (rcsid);
}

/*
*/
