
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
/* subblock.c:							*/
/* Most subroutines related to the subblock,			*/
/* including:							*/
/*   init_sb:	        initialise the subblock structure	*/
/*   calc_sb:	        perform on timestep on subblock		*/
/*   write_sb:	        write out a subblock to bb_out.dat file	*/
/*   init_c_elm:	read in CAP ca element numbers          */
/*   init_c_elm_solo:   set c_elm to cell # for test purposes   */
/* Other subblock subroutine NOT included are:                  */
/*   sb_nuc:            (sb_nuc.c) Set nuclei for a subblock    */
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Jul 1, 1998                                                  */
/****************************************************************/
/* 	MODIFIED by:						*/
/*  PDL: July 2, 1998						*/
/*  PDL: Aug 22, 1998						*/
/****************************************************************/
/****** To Do List **********************************************/
/*General:							*/
/* 1) make initial conc mor general                						*/
/* 2)                						*/
/****************************************************************/
/*RCS Id:$Id: subblock.c 1386 2008-09-25 11:47:24Z  $*/
#include <stdio.h>
#include <stdlib.h>

#include "machine.h"
#include "blocks.h"
#include "umat_matrix.h"
#include "writeblocks.h"
#include "read_sb.h"
#include "find_max.h"
#include "nearnode.h"
#include "interp.h"

/* functions used from sb_nuc.c */
extern int sb_nuc (BB_struct * bp, int sbnum);

/*from sb_umat_step.c*/
extern int sb_umat_step (BB_struct * bp, int sbnum);

/*from pore.w*/
extern int sb_pore (BB_struct * bp, int sbnum);

/* functions used from sb_temp_calc.c */
extern int sb_temp_calc (BB_struct * bp, int sbnum);

/* subroutines later in the file ... */
int init_c_elm (BB_struct *, int);
int init_c_elm_solo (BB_struct *, int);
int init_sb_interp (BB_struct * bp, NODENB_str ** nndp);

/*Subroutine to set up interpolation grid */
int init_sb_interp (BB_struct * bp, NODENB_str ** nndp)
{
  int nr, nz, errors = 0;

#ifdef VERTGRAD
  nr = bp->nc[XAXIS];
  nz = bp->nc[YAXIS];
#else
  nz = bp->nc[XAXIS];
  nr = bp->nc[YAXIS];
#endif

  /*First allocate the structure body */
  if (!(*nndp = (NODENB_str *) malloc (sizeof (NODENB_str)))) {
    fprintf (stderr, "ERROR:init_sb_interp: nnd malloc failed\n");
    errors++;
  }
  /*Allocate the node number array in RADIAL (y) direction */
  if (!((*nndp)->nl = (int *) malloc (nr * sizeof (int)))) {
    fprintf (stderr, "ERROR:init_sb_interp: nl malloc failed\n");
    errors++;
  }

  /*Allocate array for right-hand radial node weight factor */
  if (!((*nndp)->wr = (CA_FLOAT *) malloc (nr * sizeof (CA_FLOAT)))) {
    fprintf (stderr, "ERROR:init_sb_interp: nnd malloc failed\n");
    errors++;
  }
  /*Allocate array for left-hand radial node weight factor */
  if (!((*nndp)->wl = (CA_FLOAT *) malloc (nr * sizeof (CA_FLOAT)))) {
    fprintf (stderr, "ERROR:init_sb_interp: nnd malloc failed\n");
    errors++;
  }
  /*Allocate the node number array in AXIAL (x) direction */
  if (!((*nndp)->nd = (int *) malloc (nz * sizeof (int)))) {
    fprintf (stderr, "ERROR:init_sb_interp: nd malloc failed\n");
    errors++;
  }
  /*Allocate array for upper axial node weight factor */
  if (!((*nndp)->wd = (CA_FLOAT *) malloc (nz * sizeof (CA_FLOAT)))) {
    fprintf (stderr, "ERROR:init_sb_interp: nnd malloc failed\n");
    errors++;
  }
  /*Allocate array for lower axial node weight factor */
  if (!((*nndp)->wu = (CA_FLOAT *) malloc (nz * sizeof (CA_FLOAT)))) {
    fprintf (stderr, "ERROR:init_sb_interp: nnd malloc failed\n");
    errors++;
  }
  return (errors);
}                               /*end of init_sb_interp */

/********************************************************/
/********************************************************/
/*      Subroutine to initialize subblock	  	*/
/********************************************************/
/********************************************************/
/* pass in a number of sb to be inited, create	  */
/* and initialise it.				  */
int init_sb (BB_struct * bp, int sbnum)
{
  int errors = 0;
  int i, j, k;                  /* tmp counters */
  SB_struct *sp;

  fprintf (stderr, "InitSb A...\n");
   /************************************************/
  /* Malloc the SB structure                      */
   /************************************************/
  /* malloc SB structure */
  if (!(bp->sb[sbnum] = (SB_struct *) malloc (sizeof (SB_struct)))) {
    fprintf (stderr, "ERROR: SB_struct malloc failed\n");
    return (1);
  }
  sp = bp->sb[sbnum];

   /********************************************************/
  /* malloc array to hold interpolation factors         */
  /* via subroutine init_sb_interp                      */
  /* This is needed to find temperature of non-open sb's */
   /********************************************************/

  if (bp->ctrl->fgrid_input == TRUE) {  /* interpolation structure */
    if (init_sb_interp (bp, &(sp->nnd)) != 0) {
      fprintf (stderr, "ERROR: init_sb: init_sb_interp failed.\n");
      return (1);
    }
    if (init_sb_interp (bp, &(sp->nnd_next)) != 0) {
      fprintf (stderr, "ERROR: init_sb: init_sb_interp failed.\n");
      return (1);
    }

  }

  /* end init interpolation structure */
 /************************************************/
  /* Set all the general values in the SB str.    */
 /************************************************/
  sp->num = sbnum;
  sp->done = FALSE;
  sp->open = FALSE;
  /* set initial Nucleation Values... */
  sp->ngr = 0;
  sp->ncsolid = 0;              /* used as check, MUST CHANGE when working with CAP */
  sp->sbnuc.Ni_active = 0;
  sp->sbnuc.N_nuc_old = 0.0;
  sp->sbnuc.N_nuc = 0.0;
  sp->sbnuc.Ni_sum = 0.0;
  sp->sbnuc.fract_nuc = 0.0;
  /* calc the orig. of the sb */
  if (bp->ntsb == 1) {
    for (i = 0; i < 3; i++) {
      sp->orig_sb[i] = bp->orig_bb[i];
      sp->orig_bb_idx[i] = 0;
    }
  } else {
    /* multiple subblocks */
     fprintf (stderr, "DEBUG: %s Reminder -- Multiple subblocks not tested. Test this routine and remove this message. \n",__func__);
     fprintf (stderr, "DEBUG: File: subblock.c line %i \n",__LINE__);
#ifndef DOMULTIBLOCK08
     exit (2);
#endif

    k = sbnum / (bp->nsb[0] * bp->nsb[1]);
    j = (sbnum - k * (bp->nsb[0] * bp->nsb[1])) / bp->nsb[0];
    i = (sbnum - k * (bp->nsb[0] * bp->nsb[1]) - j * bp->nsb[0]);
    sp->orig_sb[0] = bp->orig_bb[0] + i * (bp->size_bb[0] / ((CA_FLOAT) bp->nsb[0]));
    sp->orig_sb[1] = bp->orig_bb[1] + j * (bp->size_bb[1] / ((CA_FLOAT) bp->nsb[1]));
    sp->orig_sb[2] = bp->orig_bb[2] + k * (bp->size_bb[2] / ((CA_FLOAT) bp->nsb[2]));

    /* set up the big-block relative index of the subblock */
    sp->orig_bb_idx[0] = i;
    sp->orig_bb_idx[1] = j;
    sp->orig_bb_idx[2] = k;
  }

  if (bp->ctrl->solo) {
    /* set initial Thermal Values... */
    sp->Tvals.Tinit = bp->Tinit;
    sp->Tvals.Tavg = bp->Tinit;
    sp->Tvals.TminReached = bp->Tinit;
    sp->Tvals.Tmin = bp->Tinit;
    sp->Tvals.Tmax = bp->Tinit;
    sp->Tvals.fsavg = 0.0;
    sp->Tvals.del_fs = 0.0;
      /**************************************************/
    /* read in a binary file of grain numbers - 8 bit */
    /* or Boundaries if in RECRYSTALLIZE mode         */
      /**************************************************/
    /*Set up the pore array for the subblock */
  }

  return (errors);
}                               /* end of init_sb subroutine */

/****************************************************************/
/****************************************************************/
/* Subroutine to read in the cap element number that a          */
/* cell belongs to from CAP CA file.                            */
/* IN:                                                          */
/*    BB_struct *bp:    bigblock structure                      */
/*    int sbnum:        index for subblock to read in           */
/* OUT:                                                         */
/*    BB_struct bp->sb[sbnum].c_elm[n]:    element number that  */
/*                      each cell is in.                        */
/****************************************************************/
/****************************************************************/
int init_c_elm (BB_struct * bp, int sbnum)
{
#ifdef JUNK
  int i, j, k, itmp;            /* tmp counters */
  long offset, ltmp;
  FILE *tfd;                    /* tmp filehandle */

  ltmp = (long) bp->sb_mask[sbnum];
  tfd = bp->ibb_fd;
  fseek (tfd, ltmp, SEEK_SET);
  /*  size_t fread (void *ptr, size_t size, size_t nitems, FILE *stream); */
  itmp = fread (bp->sb[sbnum]->c_elm, sizeof (int), (long) bp->ncsb, tfd);
  if (itmp != bp->ncsb) {
    fprintf (stderr, "ERROR: only %d [should be %d] cells read from file\n", itmp, bp->ncsb);
  }

  fseek (tfd, ltmp, 0);
#endif /*JUNK*/
    return (0);
}                               /* end of init_c_elm subroutine */

/****************************************************************/
/****************************************************************/
/* Subroutine to set cell element number if in solo mode        */
/* This is just a fudge for checking output...                  */
/* IN:                                                          */
/*    BB_struct *bp:    bigblock structure                      */
/*    int sbnum:        index for subblock to read in           */
/* OUT:                                                         */
/*    BB_struct bp->sb[sbnum].c_elm[n]:    element number that  */
/*                      each cell is in.                        */
/****************************************************************/
/****************************************************************/
int init_c_elm_solo (BB_struct * bp, int sbnum)
{
  int i, j, k, itmp;            /* tmp counters */
  int *ptr_c_elm;               /* tmp counters */

  ptr_c_elm = bp->sb[sbnum]->c_elm;
  for (i = 0; i < bp->ncsb; i++) {
    *(ptr_c_elm++) = (10000 * sbnum) + i;
  }

  return (0);
}                               /* end of init_c_elm_solo subroutine */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_subblock_c ()
{
  static char const rcsid[] = "$Id: subblock.c 1386 2008-09-25 11:47:24Z  $";

  return (rcsid);
}

/* end of rcs_id_subblock_c subroutine */
/*
*/
