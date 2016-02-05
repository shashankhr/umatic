
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
#include "umat_matrix.h"
#include "props.h"

int interface_normal (BB_struct * bp, int sbnum)
{
  SB_struct *sp;
  int i, j, k;                  /* tmp counter */
  int *nid;
  int nx, ny, nz, skip;
  CA_FLOAT *ofs, *ofs_e, *ofs_w, *ofs_n, *ofs_s, *ofs_ne, *ofs_se, *ofs_nw, *ofs_sw;
  CA_FLOAT fs_e, fs_w, fs_s, fs_n;
  CA_FLOAT *normal_x, *normal_y;

  /*float size_cell_x, size_cell_y; */

  sp = bp->sb[sbnum];

  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  skip = 2 * (nx + 2);

  /* size_cell_x = bp->size_c[0];
     size_cell_y = bp->size_c[1]; */

  /*set local pointer */
  normal_x = sp->norm_x;
  normal_y = sp->norm_y;
  ofs = bp->ftmp_one;
  nid = sp->index;

  /*move the ptr from outside corner to inside corner */
  ofs += bp->cubeptr.flist[0][START];

  ofs_e = ofs + 1;
  ofs_w = ofs - 1;
  ofs_n = ofs + nx + 2;
  ofs_s = ofs - nx - 2;
  ofs_ne = ofs_n + 1;
  ofs_nw = ofs_n - 1;
  ofs_se = ofs_s + 1;
  ofs_sw = ofs_s - 1;

   /***********************************************/
  /* Beginning of main loop(s)                   */
   /***********************************************/

  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {

        fs_e = 0.25 * (*ofs_se + *ofs_ne + 2 * (*ofs_e));
        fs_w = 0.25 * (*ofs_sw + *ofs_nw + 2 * (*ofs_w));
        fs_s = 0.25 * (*ofs_sw + *ofs_se + 2 * (*ofs_s));
        fs_n = 0.25 * (*ofs_ne + *ofs_nw + 2 * (*ofs_n));

        if (*nid == 2) {
          *normal_x = 0.5 * (fs_w - fs_e);
          *normal_y = 0.5 * (fs_s - fs_n);
        } else {
          *normal_x = 0.;
          *normal_y = 0.;
        }
        nid++;
        normal_x++;
        normal_y++;
        ofs++;
        ofs_e++;
        ofs_w++;
        ofs_n++;
        ofs_s++;
        ofs_ne++;
        ofs_nw++;
        ofs_se++;
        ofs_sw++;
      }                         /* end of I loop */
      ofs += 2;
      ofs_e += 2;
      ofs_w += 2;
      ofs_n += 2;
      ofs_s += 2;
      ofs_ne += 2;
      ofs_nw += 2;
      ofs_se += 2;
      ofs_sw += 2;

    }                           /* end of J loop */
    ofs += skip;
    ofs_e += skip;
    ofs_w += skip;
    ofs_n += skip;
    ofs_s += skip;
    ofs_ne += skip;
    ofs_nw += skip;
    ofs_se += skip;
    ofs_sw += skip;

  }                             /* end of K loop */
  return (1);
}

int interface_curvature (BB_struct * bp, int sbnum)
{
  SB_struct *sp;
  int i, j, k;                  /* tmp counter */
  int nx, ny, nz, skip;
  int *nid;
  CA_FLOAT *normal_x, *onx, *onx_e, *onx_w, *onx_n, *onx_s, *onx_ne, *onx_nw, *onx_se, *onx_sw;
  CA_FLOAT *normal_y, *ony, *ony_e, *ony_w, *ony_n, *ony_s, *ony_ne, *ony_nw, *ony_se, *ony_sw;
  CA_FLOAT normal, normal_e, normal_w, normal_s, normal_n, normal_ne, normal_nw, normal_se, normal_sw;
  CA_FLOAT *curv;
  CA_FLOAT k_x, k_y;

  /*float size_cell_x, size_cell_y; */

  sp = bp->sb[sbnum];
  bp->cubeptr.curr = sbnum;

  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  skip = 2 * (nx + 2);

  normal_x = sp->norm_x;
  normal_y = sp->norm_y;

  onx = bp->ftmp_nx;
  ony = bp->ftmp_ny;

  fcopy_matrix (PAD, onx, normal_x, bp, NULL, sbnum);
  fcopy_matrix (PAD, ony, normal_y, bp, NULL, sbnum);

  onx += bp->cubeptr.flist[0][START];
  ony += bp->cubeptr.flist[0][START];

  /* size_cell_x = bp->size_c[0];
     size_cell_y = bp->size_c[1]; */

  /*set local pointer */
  onx_e = onx + 1;
  onx_w = onx - 1;
  onx_n = onx + nx + 2;
  onx_s = onx - nx - 2;
  onx_ne = onx_n + 1;
  onx_nw = onx_n - 1;
  onx_se = onx_s + 1;
  onx_sw = onx_s - 1;
  ony_e = ony + 1;
  ony_w = ony - 1;
  ony_n = ony + nx + 2;
  ony_s = ony - nx - 2;
  ony_ne = ony_n + 1;
  ony_nw = ony_n - 1;
  ony_se = ony_s + 1;
  ony_sw = ony_s - 1;

  nid = sp->index;
  curv = sp->curv;

   /***********************************************/
  /* Beginning of main loop(s)                   */
   /***********************************************/

  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
        if (*nid == 2) {
          normal = sqrt (*onx * (*onx) + *ony * (*ony));
          normal_e = sqrt (*onx_e * (*onx_e) + *ony_e * (*ony_e));
          normal_w = sqrt (*onx_w * (*onx_w) + *ony_w * (*ony_w));
          normal_s = sqrt (*onx_s * (*onx_s) + *ony_s * (*ony_s));
          normal_n = sqrt (*onx_n * (*onx_n) + *ony_n * (*ony_n));
          normal_ne = sqrt (*onx_ne * (*onx_ne) + *ony_ne * (*ony_ne));
          normal_nw = sqrt (*onx_nw * (*onx_nw) + *ony_nw * (*ony_nw));
          normal_se = sqrt (*onx_se * (*onx_se) + *ony_se * (*ony_se));
          normal_sw = sqrt (*onx_sw * (*onx_sw) + *ony_sw * (*ony_sw));

          *onx = (normal > 0) ? *onx / normal : 0;
          *ony = (normal > 0) ? *ony / normal : 0;

          *onx_e = (normal_e > 0) ? *onx_e / (normal_e) : 0;
          *onx_w = (normal_w > 0) ? *onx_w / (normal_w) : 0;
          *onx_s = (normal_s > 0) ? *onx_s / (normal_s) : 0;
          *onx_n = (normal_n > 0) ? *onx_n / (normal_n) : 0;
          *onx_ne = (normal_ne > 0) ? *onx_ne / (normal_ne) : 0;
          *onx_nw = (normal_nw > 0) ? *onx_nw / (normal_nw) : 0;
          *onx_se = (normal_se > 0) ? *onx_se / (normal_se) : 0;
          *onx_sw = (normal_sw > 0) ? *onx_sw / (normal_sw) : 0;

          *ony_e = (normal_e > 0) ? *ony_e / (normal_e) : 0;
          *ony_w = (normal_w > 0) ? *ony_w / (normal_w) : 0;
          *ony_s = (normal_s > 0) ? *ony_s / (normal_s) : 0;
          *ony_n = (normal_n > 0) ? *ony_n / (normal_n) : 0;
          *ony_ne = (normal_ne > 0) ? *ony_ne / (normal_ne) : 0;
          *ony_nw = (normal_nw > 0) ? *ony_nw / (normal_nw) : 0;
          *ony_se = (normal_se > 0) ? *ony_se / (normal_se) : 0;
          *ony_sw = (normal_sw > 0) ? *ony_sw / (normal_sw) : 0;

          *onx_e = 0.25 * (*onx_se + *onx_ne + 2 * (*onx_e));
          *onx_w = 0.25 * (*onx_sw + *onx_nw + 2 * (*onx_w));
          *ony_s = 0.25 * (*ony_se + *ony_sw + 2 * (*ony_s));
          *ony_n = 0.25 * (*ony_nw + *ony_ne + 2 * (*ony_n));

          k_x = 0.5 * (*onx_e - *onx_w);
          k_y = 0.5 * (*ony_n - *ony_s);

          *curv = 1 * (k_x + k_y);
        } else {
          *curv = 0;
        }

        nid++;
        curv++;
        onx++;
        ony++;
        onx_e++;
        onx_w++;
        onx_n++;
        onx_s++;
        onx_ne++;
        onx_nw++;
        onx_se++;
        onx_sw++;
        ony_e++;
        ony_w++;
        ony_n++;
        ony_s++;
        ony_ne++;
        ony_nw++;
        ony_se++;
        ony_sw++;

      }                         /*end of I loop */
      onx += 2;
      ony += 2;
      onx_e += 2;
      onx_w += 2;
      onx_n += 2;
      onx_s += 2;
      onx_ne += 2;
      onx_nw += 2;
      onx_se += 2;
      onx_sw += 2;
      ony_e += 2;
      ony_w += 2;
      ony_n += 2;
      ony_s += 2;
      ony_ne += 2;
      ony_nw += 2;
      ony_se += 2;
      ony_sw += 2;
    }                           /*end of J loop */
    onx += skip;
    ony += skip;
    onx_e += skip;
    onx_w += skip;
    onx_n += skip;
    onx_s += skip;
    onx_ne += skip;
    onx_nw += skip;
    onx_se += skip;
    onx_sw += skip;
    ony_e += skip;
    ony_w += skip;
    ony_n += skip;
    ony_s += skip;
    ony_ne += skip;
    ony_nw += skip;
    ony_se += skip;
    ony_sw += skip;
  }                             /*end of K loop */
  return (1);
}

/***************************************************/
/* Surface Tension                                 */
/* by Wei Wang on 06-03-01                         */
/***************************************************/
/* using NASTAC 'template' rule (1999) */
int surface_curvature (BB_struct * bp, int sbnum)
{
  SB_struct *sp;
  int i, j, k;                  /* tmp counter */
  int nx, ny, nz, skip;
  int *nid;
  CA_FLOAT *ofs, *ofs_e, *ofs_w, *ofs_n, *ofs_s, *ofs_ne, *ofs_nw, *ofs_se, *ofs_sw;
  CA_FLOAT *ncv;
  CA_FLOAT fs_aver;

  sp = bp->sb[sbnum];

  nx = bp->nc[0];
  ny = bp->nc[1];
  nz = bp->nc[2];
  skip = 2 * (nx + 2);

  /* Set up local pointers */

  ofs = bp->ftmp_one;
  ncv = sp->curv;               /* array of surface curvature, size nc[i]^3 */
  nid = sp->index;

  /*move the ptr to outside corner to inside corner */
  ofs += bp->cubeptr.flist[0][START];

  ofs_e = ofs + 1;
  ofs_w = ofs - 1;
  ofs_n = ofs + nx + 2;
  ofs_s = ofs - nx - 2;
  ofs_ne = ofs_n + 1;
  ofs_nw = ofs_n - 1;
  ofs_se = ofs_s + 1;
  ofs_sw = ofs_s - 1;

   /***********************************************/
  /* Beginning of main loop(s)                   */
   /***********************************************/

  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {

        if (*nid != 0) {

          fs_aver = *ofs + *ofs_e + *ofs_w + *ofs_n + *ofs_s + *ofs_ne + *ofs_nw + *ofs_se + *ofs_sw;

          *ncv = 0.0 * (1.0 - 0.22222222 * fs_aver);    /* NASTAC 1999 */
        }
        nid++;
        ncv++;
        ofs++;
        ofs_e++;
        ofs_w++;
        ofs_n++;
        ofs_s++;
        ofs_ne++;
        ofs_nw++;
        ofs_se++;
        ofs_sw++;
      }                         /*end of I loop */
      ofs += 2;
      ofs_e += 2;
      ofs_w += 2;
      ofs_n += 2;
      ofs_s += 2;
      ofs_ne += 2;
      ofs_nw += 2;
      ofs_se += 2;
      ofs_sw += 2;
    }                           /* end of J loop */
    ofs += skip;
    ofs_e += skip;
    ofs_w += skip;
    ofs_n += skip;
    ofs_s += skip;
    ofs_ne += skip;
    ofs_nw += skip;
    ofs_se += skip;
    ofs_sw += skip;
  }                             /* end of K loop */

  return (1);
}

char const *rcs_id_curvature_c ()
{
  static char const rcsid[] = "$Id: curvature.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}
