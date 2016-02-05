
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
#include <math.h>
#include <string.h>
#include "machine.h"
#include "blocks.h"
#include "umat_histo.h"
#include "castats.h"

extern int getxyz (int cellnum, int *nc, int *Cell);
void write_grain_histo (BB_struct * bp, int stat_flag)
{

  int errors = 0;
  Histo_struct histoparams;
  FILE *listfile;
  char filename[MAX_STRING_LEN];
  int *grain_list;
  int grain_ctr = 0;
  Ind_grain **gr;
  int ngr;
  int i, ii;
  int zone, gzone, nzones, cellxyz[3];
  int gzx, gzy;                 /* x and y index of current grain size zone */
  int nzx, nzy;                 /* zones for grain size histo's */
  int zonex = 0, zoney = 0;     /* zones for grain size histo's */
  CA_FLOAT *size, *sizep, thissize, thisbin, grainstat[4];

  /* total number of grains -- possibly including deactivated grains */
  /* needs to be the number of entries in the bp->gr array to be considered */
  ngr = bp->nprops.ngr;
  gr = bp->gr;
  grain_ctr = 0;
  nzx = bp->nzones[0];
  nzy = bp->nzones[1];
  nzones = nzx * nzy;

  /* set up the histogram parameters */
  histoparams.binsize = G_SIZE_BINSIZE;
  histoparams.minbin = G_SIZE_MINBIN;
  histoparams.nbins = G_SIZE_NBINS;

  sprintf (filename, "G_H_%s.csv", bp->ctrl->fn_base);
  if (stat_flag == FIRST_CALL) {
    thisbin = histoparams.minbin;
    listfile = fopen (filename, "w");
    fprintf (listfile, "Histogram of log grain size\n");
    fprintf (listfile, "binsize,%.5g\n", histoparams.binsize);
    fprintf (listfile, "minbin,%.5g\n", histoparams.minbin);
    fprintf (listfile, "nbins,%i\n", histoparams.nbins);
    fprintf (listfile, "gzonex,gzoney,ngrains,step,avgcells,min,max,sdev");
    for (i = 0; i <= histoparams.nbins; i++) {
      fprintf (listfile, ",%.5g", thisbin);
      thisbin += histoparams.binsize;
    }
    fprintf (listfile, "\n");

  } else {
    listfile = fopen (filename, "a");
  }

  if (ngr != 0) {

    size = malloc (ngr * sizeof (CA_FLOAT));

    for (zoney = 0, zone = 0; zoney < nzy; zoney++) {
      for (zonex = 0; zonex < nzx; zonex++, zone++) {
        /** \todo  do this properly, go through the grains ONCE and put into the correct zone array  -- general - output*/
        /*also three-dimensional zones? */

        memset (size, 0, ngr * sizeof (CA_FLOAT));
        sizep = size;
        grain_ctr = 0;
        init_stat_val (histoparams.stat);
        /* collate the grain data into the size list */
        for (i = 1; i < ngr; i++) {
          if (gr[i] == NULL)
            continue;

          /*
             errors += getxyz(gr[i]->cell,bp->nc,cellxyz);
           */

          for (ii = 0; ii < 3; ii++) {
            cellxyz[ii] = (int) (FLOOR ((gr[i]->max[ii] + gr[i]->min[ii]) / 2));
          }

          gzx = (cellxyz[0] * nzx) / (bp->nc[0]);
          gzy = (cellxyz[1] * nzy / bp->nc[1]);
          gzone = gzx + nzx * gzy;
          if (gzone == zone) {
            thissize = (CA_FLOAT) (gr[i]->ncells);
            add_stat_val (histoparams.stat, thissize);
            *sizep++ = LOG (thissize * bp->vol_c * 1e18);
            grain_ctr++;
          }                     /* end check zone */
        }                       /* end go through grains */

        histoparams.ndata = grain_ctr;
        calc_stat_val (histoparams.stat, grain_ctr);
        fprintf (listfile, "%i,%i,%i,", zonex, zoney, grain_ctr);
        errors += umat_histo (listfile, &histoparams, size, bp->step);
      }
    }

    free (size);

  } else {
    errors += umat_histo (listfile, &histoparams, NULL, bp->step);
  }

  fclose (listfile);
  return;
}

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file */
char const *rcs_id_write_grain_histo_c ()
{
  static char const rcsid[] = "$Id: write_grain_histo.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}

/* end of rcs_id_subroutine */
