
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
#include "safeopen.h"
#include "machine.h"
#include "blocks.h"
#include "read_sb.h"
extern void free_sb_arrays (BB_struct * bp, int sbnum);
extern int free_pore (PORE_str * porelist, int npores);

/* print an ascii listing of all the grains and free them up*/
/**  \todo  should seperate the write and free bits -- multiblock */
/* allowing optional writing/bin dump/maybe none */
void grain_write (BB_struct * bp, int sbnum)
{
  int i = 0, j = 0, *gp, *wf;
  Ind_grain *grainp;
  FILE *fp;
  char fname[MAX_STRING_LEN];
  SB_struct *sp;

  /* flag for already written grains - multiblock */
  /* wf: written-flag */
  wf = (int *) calloc (bp->nprops.ngr, sizeof (int));

  sp = bp->sb[sbnum];

  sprintf (fname, "GR_%s_sb%i.csv", bp->ctrl->fn_base, sbnum);
  fp = fopen (fname, "w");

  fprintf (fp, "num,state,sbnum,cell,ncells,ngrow,blocks,Tnuc,Tund,Conc,xmax,ymax,zmax,xnuc,ynuc,znuc,xmin,ymin,zmin,angle\n");
  /* loop through and write grain info */
  for (gp = sp->gr; gp < (sp->gr + bp->ncsb); gp++) {
    if (bp->gr[*gp] == NULL)
      continue;                 /*already destroyed */
    if (*(wf + *gp) != 0)
      continue;                 /* already written */
    grainp = bp->gr[*gp];

    /* write out the info */
    /**  \todo  protect with an option to allow freeing without writing -- multiblock */
    /* neede for LARGE models , stats only */
    fprintf (fp, "%i,%i,%i,%i,%i,%i,%i,",
             grainp->num, grainp->state, grainp->sbnum, grainp->cell, grainp->ncells, grainp->ngrow, grainp->blocks);
    fprintf (fp, "%.5g,%.5g,%.5g,", grainp->TNuc, grainp->TunderNuc, grainp->CellConcNuc);
    for (i = 0; i < 3; i++)
      fprintf (fp, "%i,", grainp->max[i]);
    for (i = 0; i < 3; i++)
      fprintf (fp, "%i,", grainp->nuccell[i]);
    for (i = 0; i < 3; i++)
      fprintf (fp, "%i,", grainp->min[i]);
    fprintf (fp, "%.5g", grainp->dir_angle);
    fprintf (fp, "\n");

    /*DBM flag for dbMalloc usage */
    /* try to trap grains which have not been freed */
#ifdef DBM
    fprintf (stderr, "trap,*gp,num,mem:,%i,%i,%i,%x,", j++, *gp, grainp->num, grainp);
#endif
    if (grainp->blocks == 1) {
      free (grainp);
      bp->gr[*gp] = NULL;
#ifdef DBM
      fprintf (stderr, "freed.\n");
#endif
    } else {
      /* set flag to not duplicate information for multiblock grains */
      *(wf + *gp) = 1;
#ifdef DBM
      fprintf (stderr, "notfreed.\n");
#endif
    }
  }                             /* end of loop through grain array */
  fclose (fp);
  free (wf);
}                               /* end of grain_write */

/*******************************/
/**** close_sb ****************/
/* finish off a sb that has ***/
/* frozen completely **********/
/* and recover the memory    */
/*******************************/

int close_sb (BB_struct * bp, int sbnum)
{
  int errors = 0;
  SB_struct *sp;

  sp = bp->sb[sbnum];
  if (sp->open != SB_OPEN) {    /* check if it is open */
    fprintf (stderr, "ERROR: close_sb: sb %i is not open.\n", sbnum);
    errors++;
    return (errors);
  }

  fprintf (stderr, "close_sb: closing sb %i\n", sbnum);
  sp->open = SB_CLOSED;
  bp->sb_mask[sbnum] = FALSE;

/*
            Regular block/slice output shoudl be called
            just before close_sb so we don't need to do it here.
*/

/****************** Write out the grain info and free************************/
  grain_write (bp, sbnum);

/****************** Write out the pore info and free************************/
  /* TODO need to do pore_write on close_sb  -- porosity  multiblock*/
  /* need to save out the porosity info for a block that is completely finished*/
  /* so that memory can be reclaimed */
  /* pore_write(bp,sbnum); */

  if (bp->ctrl->pore) {
    free_pore (sp->porelist, sp->Npores);
    free (sp->porelist);
  }

  /*free the arrays in sb */
  free_sb_arrays (bp, sbnum);

  return (errors);
}

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file */
char const *rcs_id_close_sb_c ()
{
  static char const rcsid[] = "$Id: close_sb.c 1339 2008-07-23 13:58:29Z  $";

  return (rcsid);
}

/* end of rcs_id_subroutine */
