
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

/*RCS Id:$Id: sb_diffuse_alloy_decentred.c 1356 2008-08-18 13:41:15Z  $*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "blocks.h"
#include "machine.h"
#include "umat_matrix.h"
#include "sb_diffuse.h"
#include "mould_sources.h"
extern void sb_mould_src (BB_struct * bp, Solute_props * sp, int sbnum, CA_FLOAT * osol);

extern CA_FLOAT get_dl (CA_FLOAT temp);
extern CA_FLOAT get_ds (CA_FLOAT temp);
extern CA_FLOAT getav_d (CA_FLOAT dl, CA_FLOAT ds, CA_FLOAT fs);
extern CA_FLOAT cell_temp_calc_cc (BB_struct * bp, int sbnum, int x, int y);

int sb_diffuse_alloy_decentred (BB_struct * bp, int sbnum)
{

  int errflg = 0, fileflag = 0, errors = 0;
  int nx, ny, nz, skip;
  static int wmess = 0, emess = 0, courant_messg = 0;
  int *oni, *onip, *onend;
  SB_struct *sp;
  int i, j, k, n;
  CA_FLOAT **sol_alloy_values;
  CA_FLOAT *ofs, *ocl, *oce, *nfs, *ncl, *nce;
  CA_FLOAT dtx, r, rs, rl, rsrl;
  CA_FLOAT fstot, nbfs, fs_av, conc, nbconc, nbsum;

/* set up local neighbourhood */
/* use 6cell only for now     */
  fstot = 0;
  oni = bp->nbhd.onq;           /*padded */
  onip = oni;
  onend = oni + 6;

  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  skip = 2 * (nx + 2);

/* set up local values and pointers */
  sp = bp->sb[sbnum];
  bp->cubeptr.curr = sbnum;

  nfs = sp->c_fs;
  ncl = sp->c_sol_alloy;
  nce = sp->c_eqv_alloy;

  ofs = bp->ftmp_one;
  ocl = bp->ftmp_four;
  oce = bp->ftmp_five;

  /* copy C_L, C_E to temporary buffers */
  fcopy_matrix (PAD, ocl, ncl, bp, NULL, sbnum);        /*not work for multi-blocks */
  fcopy_matrix (PAD, oce, nce, bp, NULL, sbnum);

  /* add the source amount to the interface cells */
  if (bp->ctrl->mould_src) {
    sb_mould_src (bp, &(bp->mprops.alloyprops[0]), sbnum, ocl);
  }

  ofs += bp->cubeptr.flist[0][START];
  ocl += bp->cubeptr.flist[0][START];
/*set up parameter values*/

  dtx = bp->delt / (bp->size_c[0] * bp->size_c[0]);     /* dt /dx^2 */
  rs = (bp->mprops.alloyprops[0].Dsol[0] * dtx);      /* D_S * dt /dx^2 */
  rl = (bp->mprops.alloyprops[0].Dliq * dtx);      /* D_L * dt /dx^2 */
  rsrl = rs * rl;               /* precalc to save time in loop */

/* check courant stability*/
  if (rs > COURANT_LIMIT || rl > COURANT_LIMIT) {
    if (courant_messg < MAX_WARN_MSG)
      fprintf (stderr, "SB_DIFFUSE_ALLOY: WARNING: Possible instability by Courant criterion!\n Solid: , %1.2e, Liquid: %1.2e\n", rs,
               rl);
    courant_messg++;
#ifdef ERROR_EXIT
    if (courant_messg > WARN_EXIT_LIMIT) {
      fprintf (stderr, "SB_DIFFUSE_ALLOY: ERROR_EXIT: Courant Stability warning limit exceeded. %i warnings.\n", courant_messg);
      exit (courant_messg);
    }
#endif /*ERROR_EXIT */
  }

      /************************************************/
  /* now calculate the finite difference */
      /************************************************/
  /* DIFFUSION LOOP                               */
      /************************************************/
  /* Run through all cells updating as needed.    */
      /************************************************/
  for (k = 0; k < nz; k++) {    /* loop cells in z direction */
    for (j = 0; j < ny; j++) {  /* loop cells in y direction */
      for (i = 0; i < nx; i++) {        /* loop cells in x direction */
        /* skip cells that are not in the casting */
        if (*ofs == NOT_CASTING) {
        } else {

          nbsum = 0;
          conc = *ocl;
          for (onip = oni; onip < onend; onip++) {
            nbfs = *(ofs + *onip);
            /* skip nb cells that are not in the casting */
              if (nbfs == NOT_CASTING)
                continue;
              /* averaged frac solid */
              fs_av = 0.5 * (*ofs + nbfs);
            /* Linear Averaged diff coeff */
            r = rs * fs_av + rl * (1 - fs_av);
            nbconc = *(ocl + *onip);
            nbsum += r * (nbconc - conc);
          }                     /* end of neighbour sum loop */
          *nce += nbsum;        /* increment of equivalent cocentration */
          if (*nce < 0.) {
            if (emess < MAX_ERR_MSG) {
              fprintf (stderr, "ERROR:sb_diffuse_alloy_decentred:Instability C_E=%g I=%d J=%d K=%d\n", *nce, i, j, k);
            }
            emess++;
          }
          fstot += *ofs;
        }                       /* end of NOT_CASTING test */
        ofs++;
        ocl++;
        nce++;
      }                         /*x */
      ofs += 2;
      ocl += 2;
    }                           /*y */
    ofs += skip;
    ocl += skip;
  }                             /*z */

  sp->Tvals.fsavg = fstot / bp->ncsb;   /*fix up fraction solid */
  return (errflg);
}                               /* end of sb_diffuse */

/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *sb_diffuse_alloy_decentred_c ()
{
  static char const rcsid[] = "$Id: sb_diffuse_alloy_decentred.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}
