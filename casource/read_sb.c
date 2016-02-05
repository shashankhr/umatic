
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

/*RCS Id:$Id: read_sb.c 1341 2008-07-23 15:23:30Z  $*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "machine.h"
#include "blocks.h"
extern int init_new_grain (BB_struct * bp, int igr, int sbnum, int xcell, int ycell, int zcell, int nc);

/* read in a binary block of memory into a subblock - 8 bit data of grain-number*/
int read_bin_sb (BB_struct * bp, int sbnum)
{
  int errors = 0, i;
  int num_data, nsolid = 0;
  unsigned char *data, *dataorig;
  CA_FLOAT *fsp;
  int *gp;
  FILE *fp;

  gp = bp->sb[sbnum]->gr;
  fsp = bp->sb[sbnum]->c_fs;
  num_data = bp->ncsb;

  fprintf (stderr, "Entering read_bin_sb, sb %i \n", sbnum);
  if ((fp = fopen (bp->ctrl->fn_inp, "r")) == NULL)
    fprintf (stderr, "ERROR:read_bin_sb: can't open input file [%s]\n", bp->ctrl->fn_inp);
  bp->ctrl->fd_inp = fp;

  if (!(data = (unsigned char *) calloc (num_data, sizeof (unsigned char))))
    fprintf (stderr, "ERROR:read_bin_sb: Couldn't calloc data space.\n");
  dataorig = data;

  if (fread (data, sizeof (unsigned char), num_data, fp) == NULL)
    fprintf (stderr, "ERROR:read_bin_sb: Couldn't read subblock binary data\n");

  for (i = 0; i < num_data; i++) {
    *gp = (int) *data;
    #ifndef RECRYSTALLIZE
       /* Leave the grain numbers as they are for RECRYSTALLIZE */
       /* otherwise assign 255 not NOT_CASTING, 0 to liquid */
       /* and set the grain number for the other grains*/
       if (*gp == 255)
         *fsp = NOT_CASTING;
       else if (*gp == 0)
         *fsp = LIQUID;
       else {
         *fsp = .9999;
         nsolid++;
       }
    #endif

    gp++;
    data++;
    fsp++;
  }

  bp->sb[sbnum]->Tvals.fsavg += (((CA_FLOAT) nsolid) / bp->ncsb);
  if (fclose (fp) != 0)
    fprintf (stderr, "ERROR:read_bin_sb:couldn't fclose\n");
  free (dataorig);
  return (errors);

}                               /*end of read_bin_sb */

/* just read in a single slice for pseudo 2D purposes */
int read_bin_sb_slice (BB_struct * bp, int sbnum)
{
  int errors = 0, i, slice = 0;
  int num_slices, num_data, nsolid = 0;
  unsigned char *data = NULL, *datap = NULL;
  CA_FLOAT *fsp, *fs;
  int *gp, *g;
  FILE *fp;

  fprintf (stderr, "Entering read_bin_sb_slice, sb %i slice%i\n", sbnum, slice);
  num_data = bp->nc[0] * bp->nc[1];
  num_slices = bp->nc[2];
  g = gp = ((bp->sb[sbnum]->gr) + (slice * num_data));
  fs = fsp = ((bp->sb[sbnum]->c_fs) + (slice * num_data));

  if ((fp = fopen (bp->ctrl->fn_inp, "r")) == NULL)
    fprintf (stderr, "ERROR:read_bin_sb_slice: can't open input file [%s]\n", bp->ctrl->fn_inp);
  bp->ctrl->fd_inp = fp;

  if (!(datap = data = (unsigned char *) malloc (num_data * sizeof (unsigned char))))
    fprintf (stderr, "ERROR:read_bin_sb_slice: Couldn't calloc data space.\n");

  if (fread (data, sizeof (unsigned char), num_data, fp) == NULL)
    fprintf (stderr, "ERROR:read_bin_sb_slice: Couldn't read subblock binary data\n");

  for (i = 0; i < num_data; i++) {
    *gp = (int) *data;

       /* assign 255 to NOT_CASTING, 0 to liquid */
       /* and set the grain number for the other grains*/
    if (*gp == 255)
      *fsp = NOT_CASTING;
    else if (*gp == 0)
      *fsp = LIQUID;
    else {
      *fsp = .9999;
      nsolid++;
    }

    gp++;
    data++;
    fsp++;

  }
  /* copy the input slice to all slices */
  for (i = 1; i < num_slices; i++) {
    memcpy ((g + i * num_data), g, num_data * sizeof (int));
    memcpy ((fs + i * num_data), fs, num_data * sizeof (CA_FLOAT));
  }
  nsolid *= num_slices;

  bp->sb[sbnum]->Tvals.fsavg += (((CA_FLOAT) nsolid) / bp->ncsb);
  fclose (fp);
  free (datap);

  return (errors);
}                               /*end of read_bin_sb_slice */

/* fix up the grain structures after reading in a subblock */
int fix_grains (BB_struct * bp, int sbnum)
{
  double dtmp;
  int errors = 0;
  int ngrains, nbound = 0;
  int i, j, k;
  int grainlist[256], ncells[256];
  int xcells[256], ycells[256], zcells[256];
  int num_data;
  int *gp;                      /*Grain Pointer */
  int *nucp;                    /*NUCleation list Pointer */
  int *elp;                     /*ELement Pointer used as recr. flag for g.b. */
  CA_FLOAT *fs, *fsp;

  for (i = 0; i < 256; i++)
    grainlist[i] = ncells[i] = xcells[i] = ycells[i] = zcells[i] = 0;

  ngrains = bp->nprops.ngr;
  gp = bp->sb[sbnum]->gr;       /*rewind */
  fsp = fs = bp->sb[sbnum]->c_fs;
  elp = bp->sb[sbnum]->c_elm;   /*rewind */
  num_data = bp->ncsb;

  for (k = 0; k < bp->nc[2]; k++) {
    for (j = 0; j < bp->nc[1]; j++) {
      for (i = 0; i < bp->nc[0]; i++) {
        {
          if (*gp <= 1) {
            *gp = 0;
          } else if (*gp == 255) {
            *gp = NOT_CASTING;
            *fsp = NOT_CASTING;
          } else if (grainlist[*gp] == 0) {
            ngrains++;
            grainlist[*gp] = ngrains;
            xcells[ngrains] = i;
            ycells[ngrains] = j;
            zcells[ngrains] = k;
          }
          fsp++;
          gp++;
          elp++;
        }
      }
    }
  }
  gp = bp->sb[sbnum]->gr;       /*rewind grain */
  elp = bp->sb[sbnum]->c_elm;   /*rewind element */
  for (i = 0; i < num_data; i++) {
    /* renumber all the grains avoiding any missing numbers */
    if (*gp != NOT_CASTING) {
      *gp = grainlist[*gp];
      ncells[*gp]++;
    }
    gp++;
    elp++;

  }
  for (i = bp->nprops.ngr; i <= ngrains; i++) {
    Ind_grain *g;
    CA_FLOAT a0, a1, a2;
    CA_FLOAT c0, s0, c1, s1, c2, s2;

    init_new_grain (bp, i, sbnum, xcells[i], ycells[i], zcells[i], ncells[i]);
    g = bp->gr[i];
    /* generate rotate matrix for grain */
    /* added by Wei WANG on 20-09-02 */
    /* copied here as quick way to choose angle for */
    /* grain which is read in from the image file */
    /**  \todo  improve the input for grain angle -- decentred - maybe obsolete */
    if (bp->ctrl->decentred_octahedron == TRUE) {
#ifdef XY_ANGLE
      a0 = XY_ANGLE;
#else
      a0 = dtmp;                /* need to be changed */
#endif /*XY_ANGLE */
      a1 = 0.0;
      a2 = 0.0;
      g->ang[0] = a0;
      g->ang[1] = a1;
      g->ang[2] = a2;

      c0 = g->cang[0] = cos (g->ang[0]);
      s0 = g->sang[0] = sin (g->ang[0]);
      c1 = g->cang[1] = cos (g->ang[1]);
      s1 = g->sang[1] = sin (g->ang[1]);
      c2 = g->cang[2] = cos (g->ang[2]);
      s2 = g->sang[2] = sin (g->ang[2]);

      /* rotation matrix */
      g->g[0][0] = c2 * c0 - s2 * s1 * s0;
      g->g[1][0] = c2 * s0 + s2 * c1 * c0;
      g->g[2][0] = s2 * s1;

      g->g[0][1] = -1. * s2 * c0 - c2 * c1 * s0;
      g->g[1][1] = -1. * s2 * s0 + c2 * c1 * c0;
      g->g[2][1] = c2 * s1;

      g->g[0][2] = s1 * s0;
      g->g[1][2] = -1. * s1 * c0;
      g->g[2][2] = c1;
    }

  }
  bp->nprops.ngr = ngrains;
  return (errors);
}                               /*endof fix_grains */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_read_sb_c ()
{
  static char const rcsid[] = "$Id: read_sb.c 1341 2008-07-23 15:23:30Z  $";

  return (rcsid);
}

/* end of rcs_id_read_sb_c subroutine */
/*
*/
