
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "machine.h"
#include "blocks.h"
#include "read_sb.h"
#include "interp.h"
#include "SurCellRoutines.h"

#ifdef JUNK
/* functions not used JUNK */
/* functions used from sb_temp_calc.c */
extern int sb_temp_calc (BB_struct * bp, int sbnum);

/* subroutines later in the file ... */
extern int init_c_elm (BB_struct *, int);
extern int init_c_elm_solo (BB_struct *, int);
#endif
extern int init_nuc_thresh (BB_struct * bp, int sbnum);
extern int pore_setup (BB_struct * bp, int sbnum);
extern int sb_nuc (BB_struct * bp, int sbnum);
extern int alloc_sb (BB_struct * bp, int sbnum);

/* from user_rtns.c -- only used in external mode */
extern void external_sb_set_cells (BB_struct * bp, int sbnum);

int open_sb (BB_struct * bp, int sbnum)
{
  int errors = 0;
  int i, nc;
  int oldstep;
  SB_struct *sp;

  sp = bp->sb[sbnum];
  nc = bp->ncsb;

  if (sp->open == SB_OPEN) {    /*check if already open */
    fprintf (stderr, "ERROR: open_sb: SB %i already open.\n", sbnum);
    errors++;
    return (errors);
  }
  if (sp->open == SB_CLOSED) {  /*check if already closed */
    fprintf (stderr, "ERROR: open_sb: SB %i already used.\n", sbnum);
    errors++;
    return (errors);
  }

  sp->open = SB_OPEN;
  bp->sb_mask[sbnum] = TRUE;

      /********************************************/
      /** call the routine to allocate the memory */
      /** and cross link the pointers             */
      /** for this subblock                       */
      /********************************************/
  init_SurCell (&(sp->surface));
  errors += alloc_sb (bp, sbnum);

  sp->nmould = 0;
  if (bp->ctrl->input) {
     if (bp->ctrl->extrudemould ==1){
        read_bin_sb_slice (bp, sbnum);      /*copy image to all slices of sb for now */
    }else{
        read_bin_sb (bp, sbnum); /* read whole 3d mould */
    }
    errors += fix_grains (bp, sbnum);
  }

  if (bp->ctrl->external) {
    external_sb_set_cells (bp, sbnum);
  }
  /* I have not considered the problems of inputting a partially complete slice at this point */
  /* so issue some warnings. Maybe it works anyways? -- Robert                                */
  if (bp->ctrl->input && bp->ctrl->external) {
    fprintf (stderr, "WARNING: open_sb: Both input and external flags are active.\n");
    fprintf (stderr, "WARNING: open_sb: This behaviour has not been fully implemented.\n");
    fprintf (stderr, "WARNING: open_sb: Continuing anyways but there may be problems.\n");
  }

  /*allocate nuc threshold array */

/*THUINET 05/05*/

/*      if ((bp->ctrl->block_nuc == TRUE )) { *//*nuc threshold array*/

  if (bp->ctrl->input ==TRUE || bp->ctrl->external == TRUE){
   sb_get_surface_cells (bp, sbnum);
  }

  if ((bp->ctrl->block_nuc == TRUE )&&(bp->ctrl->diffuse_alloy_poly == FALSE)) {
  /* set up the nucleation for this block*/
    if (! init_nuc_thresh(bp,sbnum) == 0 ){
      #ifndef ALWAYS_NUC
      fprintf(stderr,"... but in open_sb: ALWAYS_NUC is not defined so exiting due to init_nuc_thresh failure.\n");
      exit(1);
      #endif /*ALWAYS_NUC*/
    }

  }

  if ((bp->ctrl->block_nuc == TRUE )&&(bp->ctrl->diffuse_alloy_poly == TRUE)) { /*nuc threshold array*/
  /* set up the nucleation for this block*/
    if (! init_nuc_thresh_poly(bp,sbnum) == 0 ){
      #ifndef ALWAYS_NUC 
      fprintf(stderr,"... but in open_sb: ALWAYS_NUC is not defined so exiting due to init_nuc_thresh failure.\n");
      exit(1);
      #endif /*ALWAYS_NUC*/
    }
  }

  /*end allocate nuc threshold array */

#  ifdef USE_ELM
      /************************************************/
  /* Set the elm # if reading CAP CA file...      */
  /* or for testing, or reading input ......      */
      /************************************************/
  if (bp->ctrl->cap) {          /* Read elm # for all cells */
    init_c_elm (bp, sbnum);
  } else if (bp->ctrl->solo) {
    init_c_elm_solo (bp, sbnum);        /* set c_elm to cell # for test purposes */
#  else
  if (bp->ctrl->solo) {

#  endif /*USE_ELM */
    /* set initial Thermal Values... */
  }
      /***************************************************************/
  if (bp->ctrl->pore == TRUE) {
    errors += pore_setup (bp, sbnum);
  }
      /*************************************/
  /* pre-create the nuclei if FIXED... */
      /*************************************/
  if (bp->ctrl->fixed_nuc) {
    oldstep = bp->step;
    bp->step = INIT;
    sb_nuc (bp, sbnum);
    bp->step = oldstep;
  }
  if (bp->ctrl->fgrid_input) {
    wfact_r_calc (bp->fg, sp->nnd, bp, sbnum);
    wfact_z_calc (bp->fg, sp->nnd, bp, sbnum);
    wfact_r_calc (bp->fg_next, sp->nnd_next, bp, sbnum);
    wfact_z_calc (bp->fg_next, sp->nnd_next, bp, sbnum);
    sb_temp_setup (bp, sbnum);
  }
  sp->t_addsol = 0;
  return (errors);
}                               /*end of open_sb subroutine */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file */
char const *rcs_id_open_sb_c ()
{
  static char const rcsid[] = "$Id: open_sb.c 1390 2008-09-25 15:43:01Z  $";

  return (rcsid);
}

/* end of rcs_id_subroutine */

/*
*/
