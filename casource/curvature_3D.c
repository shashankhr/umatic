
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


int interface_normal_3D (BB_struct *bp, int sbnum)
{  
   SB_struct * sp;
   int i, j, k;                             /* tmp counter */
   int *nid;
   int nx, ny, nz, skip;
   int beta = 2;
   int gamma = 4;
   CA_FLOAT *ofs ,*ofs_e, *ofs_w, *ofs_n, *ofs_s, *ofs_ne, *ofs_se, *ofs_nw, *ofs_sw;
   CA_FLOAT *ofs_t, *ofs_b, *ofs_te, *ofs_tw, *ofs_tn, *ofs_ts, *ofs_be, *ofs_bw, *ofs_bn, *ofs_bs;
   CA_FLOAT *ofs_bse, *ofs_bsw, *ofs_bne, *ofs_bnw, *ofs_tse, *ofs_tsw, *ofs_tne, *ofs_tnw;
   CA_FLOAT fs_e, fs_w, fs_s, fs_n, fs_t, fs_b;
   CA_FLOAT *normal_x, *normal_y, *normal_z;
   /*float size_cell_x, size_cell_y;*/
   
   sp = bp->sb[sbnum];

   nx = bp->nc[0];
   ny = bp->nc[1];
   nz = bp->nc[2];
   skip = 2 * (nx + 2);
  
  /* size_cell_x = bp->size_c[0];
   size_cell_y = bp->size_c[1];*/
    
   /*set local pointer */
   normal_x = sp->norm_x;
   normal_y = sp->norm_y;
   normal_z = sp->norm_z; 
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
   
   /* PLIC extension to 3D - Dylan Ness 18/10/05 */

   ofs_t = ofs + (nx+2)*(ny+2);
   ofs_b = ofs - (nx+2)*(ny+2);
   ofs_te = ofs_t + 1;
   ofs_tw = ofs_t - 1;
   ofs_tn = ofs_t + nx + 2;
   ofs_ts = ofs_t - nx - 2;
   ofs_be = ofs_b + 1;
   ofs_bw = ofs_b - 1;
   ofs_bn = ofs_b + nx + 2;
   ofs_bs = ofs_b - nx - 2;
   ofs_bse = ofs_bs + 1;
   ofs_bsw = ofs_bs - 1;
   ofs_bne = ofs_bn + 1;
   ofs_bnw = ofs_bn - 1;
   ofs_tse = ofs_ts + 1;
   ofs_tsw = ofs_ts - 1;
   ofs_tne = ofs_tn + 1;
   ofs_tnw = ofs_tn - 1;
   
   /***********************************************/
   /* Beginning of main loop(s)                   */
   /***********************************************/

   for (k = 0; k < nz; k++) {
     for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
         
             fs_e = (*ofs_bse + *ofs_bne + *ofs_tse + *ofs_tne + beta*(*ofs_be + *ofs_se + *ofs_ne + *ofs_te) + gamma*(*ofs_e))/(4 + 4*beta + gamma);
	     fs_w = (*ofs_bsw + *ofs_bnw + *ofs_tsw + *ofs_tnw + beta*(*ofs_bw + *ofs_sw + *ofs_nw + *ofs_tw) + gamma*(*ofs_w))/(4 + 4*beta + gamma);
             fs_s = (*ofs_bsw + *ofs_bse + *ofs_tsw + *ofs_tse + beta*(*ofs_bs + *ofs_se + *ofs_sw + *ofs_ts) + gamma*(*ofs_s))/(4 + 4*beta + gamma);
             fs_n = (*ofs_bnw + *ofs_bne + *ofs_tnw + *ofs_tne + beta*(*ofs_bn + *ofs_ne + *ofs_nw + *ofs_tn) + gamma*(*ofs_n))/(4 + 4*beta + gamma);
	     fs_t = (*ofs_tne + *ofs_tnw + *ofs_tsw + *ofs_tse + beta*(*ofs_tw + *ofs_te + *ofs_tn + *ofs_ts) + gamma*(*ofs_t))/(4 + 4*beta + gamma);
	     fs_b = (*ofs_bne + *ofs_bnw + *ofs_bsw + *ofs_bse + beta*(*ofs_bw + *ofs_be + *ofs_bn + *ofs_bs) + gamma*(*ofs_b))/(4 + 4*beta + gamma);

            if (*nid == 2) { 
               *normal_x = 0.5 * (fs_w - fs_e);
               *normal_y = 0.5 * (fs_s - fs_n);
	       *normal_z = 0.5 * (fs_t - fs_b);
            }
            else {
               *normal_x = 0.;
               *normal_y = 0.;
	       *normal_z = 0.;
            }
            nid ++;
            normal_x ++;
            normal_y ++;
	    normal_z ++;
            ofs ++;
            ofs_e ++;
            ofs_w ++;
            ofs_n ++;
            ofs_s ++;
	    ofs_t ++;
	    ofs_b ++;
            ofs_ne ++;
            ofs_nw ++;
            ofs_se ++;
            ofs_sw ++;
	    ofs_te ++;
	    ofs_tw ++;
	    ofs_tn ++;
	    ofs_ts ++;
	    ofs_be ++;
	    ofs_bw ++;
	    ofs_bn ++;
	    ofs_bs ++;
	    ofs_tne ++;
	    ofs_tnw ++;
	    ofs_tse ++;
	    ofs_tsw ++;
	    ofs_bne ++;
	    ofs_bnw ++;
	    ofs_bse ++;
	    ofs_bsw ++;

	}/* end of I loop */
        ofs += 2;
        ofs_e += 2;
        ofs_w += 2;
        ofs_n += 2;
        ofs_s += 2;
        ofs_ne += 2;
        ofs_nw += 2;
        ofs_se += 2;
        ofs_sw += 2;
	ofs_te += 2;
        ofs_tw += 2;
        ofs_tn += 2;
        ofs_ts += 2;
        ofs_be += 2;
        ofs_bw += 2;
        ofs_bn += 2;
        ofs_bs += 2;
	ofs_tne += 2;
        ofs_tnw += 2;
        ofs_tse += 2;
        ofs_tsw += 2;
        ofs_bne += 2;
        ofs_bnw += 2;
        ofs_bse += 2;
        ofs_bsw += 2;


     }/* end of J loop */
     ofs += skip;
     ofs_e += skip;
     ofs_w += skip;
     ofs_n += skip;
     ofs_s += skip;
     ofs_ne += skip;
     ofs_nw += skip;
     ofs_se += skip;
     ofs_sw += skip;
     ofs_te += skip;
     ofs_tw += skip;
     ofs_tn += skip;
     ofs_ts += skip;
     ofs_be += skip;
     ofs_bw += skip;
     ofs_bn += skip;
     ofs_bs += skip;
     ofs_tne += skip;
     ofs_tnw += skip;
     ofs_tse += skip;
     ofs_tsw += skip;
     ofs_bne += skip;
     ofs_bnw += skip;
     ofs_bse += skip;
     ofs_bsw += skip;


  } /* end of K loop */
  return (1);
}   


int interface_curvature_3D (BB_struct *bp, int sbnum)
{  
   SB_struct * sp;
   int i, j, k;                            /* tmp counter */
   int nx, ny, nz, skip;
   int *nid;
   int beta = 2;
   int gamma = 4;
   CA_FLOAT *normal_x ,*onx, *onx_e, *onx_w, *onx_n, *onx_s, *onx_ne, *onx_nw, *onx_se, *onx_sw;
   CA_FLOAT *onx_t, *onx_b, *onx_te, *onx_tw, *onx_tn, *onx_ts, *onx_be, *onx_bw, *onx_bn, *onx_bs;
   CA_FLOAT *onx_tne, *onx_tnw, *onx_tsw, *onx_tse, *onx_bne, *onx_bnw, *onx_bse, *onx_bsw;
   CA_FLOAT *normal_y, *ony, *ony_e, *ony_w, *ony_n, *ony_s, *ony_ne, *ony_nw, *ony_se, *ony_sw;
   CA_FLOAT *ony_t, *ony_b, *ony_te, *ony_tw, *ony_tn, *ony_ts, *ony_be, *ony_bw, *ony_bn, *ony_bs;
   CA_FLOAT *ony_tne, *ony_tnw, *ony_tsw, *ony_tse, *ony_bne, *ony_bnw, *ony_bse, *ony_bsw;
   CA_FLOAT *normal_z, *onz, *onz_e, *onz_w, *onz_n, *onz_s, *onz_ne, *onz_nw, *onz_se, *onz_sw;
   CA_FLOAT *onz_t, *onz_b, *onz_te, *onz_tw, *onz_tn, *onz_ts, *onz_be, *onz_bw, *onz_bn, *onz_bs;
   CA_FLOAT *onz_tne, *onz_tnw, *onz_tsw, *onz_tse, *onz_bne, *onz_bnw, *onz_bse, *onz_bsw;
   CA_FLOAT normal, normal_e, normal_w, normal_s, normal_n, normal_ne, normal_nw, normal_se, normal_sw;
   CA_FLOAT normal_t, normal_b, normal_te, normal_tw, normal_tn, normal_ts, normal_be, normal_bw, normal_bn, normal_bs;
   CA_FLOAT normal_tne, normal_tnw, normal_tse, normal_tsw, normal_bne, normal_bnw, normal_bse, normal_bsw;
   CA_FLOAT *curv, curv_x, curv_y, curv_z; 
   CA_FLOAT k_xx, k_xy, k_xz, k_yx, k_yy, k_yz, k_zx, k_zy, k_zz;
  
   
   /*float size_cell_x, size_cell_y;*/
   
   sp = bp->sb[sbnum];
   bp->cubeptr.curr = sbnum;

   nx = bp->nc[0];
   ny = bp->nc[1];
   nz = bp->nc[2];
   skip = 2 * (nx + 2);
   
   normal_x = sp->norm_x;
   normal_y = sp->norm_y;
   normal_z = sp->norm_z;

   onx = bp->ftmp_nx;
   ony = bp->ftmp_ny;
   onz = bp->ftmp_nz;
  
   fcopy_matrix(PAD, onx, normal_x, bp, NULL, sbnum); 
   fcopy_matrix(PAD, ony, normal_y, bp, NULL, sbnum);
   fcopy_matrix(PAD, onz, normal_z, bp, NULL, sbnum);
   
   onx += bp->cubeptr.flist[0][START];
   ony += bp->cubeptr.flist[0][START];
   onz += bp->cubeptr.flist[0][START];

  /* size_cell_x = bp->size_c[0];
   size_cell_y = bp->size_c[1];*/
    
   /*set local pointer */
   onx_e = onx + 1;
   onx_w = onx - 1;
   onx_n = onx + nx + 2;
   onx_s = onx - (nx + 2);
   onx_t = onx + (nx+2)*(ny+2);
   onx_b = onx - (nx+2)*(ny+2);

	/*

   onx_ne = onx_n + 1;
   onx_nw = onx_n - 1;
   onx_se = onx_s + 1;
   onx_sw = onx_s - 1;
   onx_te = onx_t + 1;
   onx_tw = onx_t - 1;
   onx_tn = onx_t + (nx + 2);
   onx_ts = onx_t - (nx + 2);
   onx_be = onx_b + 1;
   onx_bw = onx_b - 1;
   onx_bn = onx_b + (nx + 2);
   onx_bs = onx_b - (nx + 2);
   onx_tne = onx_tn + 1;
   onx_tnw = onx_tn - 1;
   onx_tse = onx_ts + 1;
   onx_tsw = onx_ts - 1;
   onx_bne = onx_bn + 1;
   onx_bnw = onx_bn - 1;
   onx_bse = onx_bs + 1;
   onx_bsw = onx_bs - 1;

	*/

   ony_e = ony + 1;
   ony_w = ony - 1;
   ony_n = ony + nx + 2;
   ony_s = ony - nx - 2;
   ony_t = ony + (nx + 2);
   ony_b = ony - (nx + 2);

	/*

   ony_ne = ony_n + 1;
   ony_nw = ony_n - 1;
   ony_se = ony_s + 1;
   ony_sw = ony_s - 1;
   ony_te = ony_t + 1;
   ony_tw = ony_t - 1;
   ony_tn = ony_t + (nx + 2);
   ony_ts = ony_t - (nx + 2);
   ony_be = ony_b + 1;
   ony_bw = ony_b - 1;
   ony_bn = ony_b + (ny + 2);
   ony_bs = ony_b - (ny + 2);
   ony_tne = ony_tn + 1;
   ony_tnw = ony_tn - 1;
   ony_tse = ony_ts + 1;
   ony_tsw = ony_ts - 1; 
   ony_bne = ony_bn + 1;
   ony_bnw = ony_bn - 1;
   ony_bse = ony_bs + 1;
   ony_bsw = ony_bs - 1;

	*/

   onz_e = onz + 1;
   onz_w = onz - 1;
   onz_n = onz + (nx+2);
   onz_s = onz - (nx+2);
   onz_t = onz + (nx+2)*(ny+2);
   onz_b = onz - (nx+2)*(nx+2);
   
	/*

   onz_ne = onz_n + 1;
   onz_nw = onz_n - 1;
   onz_se = onz_s + 1;
   onz_sw = onz_s - 1;
   onz_te = onz_t + 1;
   onz_tw = onz_t - 1;
   onz_tn = onz_t + (nx + 2);
   onz_ts = onz_t - (nx + 2);
   onz_be = onz_b + 1;
   onz_bw = onz_b - 1; 
   onz_bn = onz_b + (nx + 2);
   onz_bs = onz_b - (nx + 2);
   onz_tne = onz_tn + 1;
   onz_tnw = onz_tn - 1;
   onz_tse = onz_ts + 1;
   onz_tsw = onz_ts - 1;
   onz_bne = onz_bn + 1;
   onz_bnw = onz_bn - 1;
   onz_bse = onz_bs + 1;
   onz_bsw = onz_bs - 1;

	*/

   nid = sp->index;
   curv = sp->curv;


   /***********************************************/
   /* Beginning of main loop(s)                   */
   /***********************************************/

   for (k = 0; k < nz; k++) {
     for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) { 
          if (*nid == 2){
            normal = sqrt(*onx * (*onx) + *ony * (*ony) + *onz * (*onz));  
            normal_e = sqrt(*onx_e * (*onx_e) + *ony_e * (*ony_e) + *onz_e * (*onz_e));
            normal_w = sqrt(*onx_w * (*onx_w) + *ony_w * (*ony_w) + *onz_w * (*onz_w));
            normal_s = sqrt(*onx_s * (*onx_s) + *ony_s * (*ony_s) + *onz_s * (*onz_s));
            normal_n = sqrt(*onx_n * (*onx_n) + *ony_n * (*ony_n) + *onz_n * (*onz_n));
            normal_t = sqrt(*onx_t * (*onx_t) + *ony_t * (*ony_t) + *onz_t * (*onz_t));
	    normal_b = sqrt(*onx_b * (*onx_b) + *ony_b * (*ony_b) + *onz_b * (*onz_b));

		/*

	    normal_ne = sqrt(*onx_ne * (*onx_ne) + *ony_ne * (*ony_ne) + *onz_ne * (*onz_ne));
            normal_nw = sqrt(*onx_nw * (*onx_nw) + *ony_nw * (*ony_nw) + *onz_nw * (*onz_nw)); 
            normal_se = sqrt(*onx_se * (*onx_se) + *ony_se * (*ony_se) + *onz_se * (*onz_se));
            normal_sw = sqrt(*onx_sw * (*onx_sw) + *ony_sw * (*ony_sw) + *onz_sw * (*onz_sw));
	    normal_te = sqrt(*onx_te * (*onx_te) + *ony_te * (*ony_te) + *onz_te * (*onz_te));
	    normal_tw = sqrt(*onx_tw * (*onx_tw) + *ony_tw * (*ony_tw) + *onz_tw * (*onz_tw));
	    normal_tn = sqrt(*onx_tn * (*onx_tn) + *ony_tn * (*ony_tn) + *onz_tn * (*onz_tn));
	    normal_ts = sqrt(*onx_ts * (*onx_ts) + *ony_ts * (*ony_ts) + *onz_ts * (*onz_ts));
	    normal_be = sqrt(*onx_be * (*onx_be) + *ony_be * (*ony_be) + *onz_be * (*onz_be));
	    normal_bw = sqrt(*onx_bw * (*onx_bw) + *ony_bw * (*ony_bw) + *onz_bw * (*onz_bw));
	    normal_bn = sqrt(*onx_bn * (*onx_bn) + *ony_bn * (*ony_bn) + *onz_bn * (*onz_bn));
	    normal_bs = sqrt(*onx_bs * (*onx_bs) + *ony_bs * (*ony_bs) + *onz_bs * (*onz_bs));
            normal_tne = sqrt(*onx_tne * (*onx_tne) + *ony_tne * (*ony_tne) + *onz_tne * (*onz_tne));
	    normal_tnw = sqrt(*onx_tnw * (*onx_tnw) + *ony_tnw * (*ony_tnw) + *onz_tnw * (*onz_tnw));
	    normal_tse = sqrt(*onx_tse * (*onx_tse) + *ony_tse * (*ony_tse) + *onz_tse * (*onz_tse));
	    normal_tsw = sqrt(*onx_tsw * (*onx_tsw) + *ony_tsw * (*ony_tsw) + *onz_tsw * (*onz_tsw));
	    normal_bne = sqrt(*onx_bne * (*onx_bne) + *ony_bne * (*ony_bne) + *onz_bne * (*onz_bne));
	    normal_bnw = sqrt(*onx_bnw * (*onx_bnw) + *ony_bnw * (*ony_bnw) + *onz_bnw * (*onz_bnw));
	    normal_bse = sqrt(*onx_bse * (*onx_bse) + *ony_bse * (*ony_bse) + *onz_bse * (*onz_bse));
	    normal_bsw = sqrt(*onx_bsw * (*onx_bsw) + *ony_bsw * (*ony_bsw) + *onz_bsw * (*onz_bsw));

		*/

            *onx = (normal > 0) ? *onx / normal : 0;
            *ony = (normal > 0) ? *ony / normal : 0;
	    *onz = (normal > 0) ? *onz / normal : 0;

             
	     *onx_e = (normal_e > 0) ? *onx_e / normal_e : 0;
	     *onx_w = (normal_w > 0) ? *onx_w / normal_w : 0;
	     *onx_n = (normal_n > 0) ? *onx_n / normal_n : 0;
	     *onx_s = (normal_s > 0) ? *onx_s / normal_s : 0;
	     *onx_t = (normal_t > 0) ? *onx_t / normal_t : 0;
	     *onx_b = (normal_b > 0) ? *onx_b / normal_b : 0;	    

		/*

	     *onx_ne = (normal_ne > 0) ? *onx_ne / (normal_ne) : 0;
             *onx_nw = (normal_nw > 0) ? *onx_nw / (normal_nw) : 0;
             *onx_se = (normal_se > 0) ? *onx_se / (normal_se) : 0;
             *onx_sw = (normal_sw > 0) ? *onx_sw / (normal_sw) : 0;	     
	     *onx_te = (normal_te > 0) ? *onx_te / normal_te : 0;
	     *onx_tw = (normal_tw > 0) ? *onx_tw / normal_tw : 0;
	     *onx_tn = (normal_tn > 0) ? *onx_tn / normal_tn : 0;
	     *onx_ts = (normal_ts > 0) ? *onx_ts / normal_ts : 0;
	     *onx_be = (normal_be > 0) ? *onx_be / normal_be : 0;
	     *onx_bw = (normal_bw > 0) ? *onx_bw / normal_bw : 0;
	     *onx_bn = (normal_bn > 0) ? *onx_bn / normal_bn : 0;
	     *onx_bs = (normal_bs > 0) ? *onx_bs / normal_bs : 0;
	     *onx_tne = (normal_tne > 0) ? *onx_tne / normal_tne : 0;
             *onx_tnw = (normal_tnw > 0) ? *onx_tnw / normal_tnw : 0;
             *onx_tsw = (normal_tsw > 0) ? *onx_tsw / normal_tsw : 0;
             *onx_tse = (normal_tse > 0) ? *onx_tse / normal_tse : 0;
             *onx_bne = (normal_bne > 0) ? *onx_bne / normal_bne : 0;
             *onx_bnw = (normal_bnw > 0) ? *onx_bnw / normal_bnw : 0;
             *onx_bse = (normal_bse > 0) ? *onx_bse / normal_bse : 0;
             *onx_bsw = (normal_bsw > 0) ? *onx_bsw / normal_bsw : 0;

		*/

	     *ony_e = (normal_e > 0) ? *ony_e / normal_e : 0;
	     *ony_w = (normal_w > 0) ? *ony_w / normal_w : 0;
	     *ony_n = (normal_n > 0) ? *ony_n / normal_n : 0;
	     *ony_s = (normal_s > 0) ? *ony_s / normal_s : 0;
	     *ony_t = (normal_t > 0) ? *ony_t / normal_t : 0;
	     *ony_b = (normal_b > 0) ? *ony_b / normal_b : 0;

		/*

             *ony_ne = (normal_ne > 0) ? *ony_ne / (normal_ne) : 0;
             *ony_nw = (normal_nw > 0) ? *ony_nw / (normal_nw) : 0;
             *ony_se = (normal_se > 0) ? *ony_se / (normal_se) : 0;
             *ony_sw = (normal_sw > 0) ? *ony_sw / (normal_sw) : 0;
	     *ony_te = (normal_te > 0) ? *ony_te / (normal_te) : 0;
	     *ony_tw = (normal_tw > 0) ? *ony_tw / (normal_tw) : 0;
	     *ony_tn = (normal_tn > 0) ? *ony_tn / (normal_tn) : 0;
	     *ony_ts = (normal_ts > 0) ? *ony_ts / (normal_ts) : 0;
	     *ony_be = (normal_be > 0) ? *ony_be / (normal_be) : 0;
	     *ony_bw = (normal_bw > 0) ? *ony_bw / (normal_bw) : 0;
	     *ony_bn = (normal_bn > 0) ? *ony_bn / (normal_bn) : 0;
	     *ony_bs = (normal_bs > 0) ? *ony_bs / (normal_bs) : 0;
	     *ony_tne = (normal_tne > 0) ? *ony_tne / (normal_tne) : 0;
             *ony_tnw = (normal_tnw > 0) ? *ony_tnw / normal_tnw : 0;
             *ony_tsw = (normal_tsw > 0) ? *ony_tsw / normal_tsw : 0;
             *ony_tse = (normal_tse > 0) ? *ony_tse / normal_tse : 0;
             *ony_bne = (normal_bne > 0) ? *ony_bne / normal_bne : 0;
             *ony_bnw = (normal_bnw > 0) ? *ony_bnw / normal_bnw : 0;
             *ony_bse = (normal_bse > 0) ? *ony_bse / normal_bse : 0;
             *ony_bsw = (normal_bsw > 0) ? *ony_bsw / normal_bsw : 0;
	
		*/

	     *onz_e = (normal_e > 0) ? *onz_e / normal_e : 0;
	     *onz_w = (normal_w > 0) ? *onz_w / normal_w : 0;
	     *onz_n = (normal_n > 0) ? *onz_n / normal_n : 0;
	     *onz_s = (normal_s > 0) ? *onz_s / normal_s : 0;
	     *onz_t = (normal_t > 0) ? *onz_t / normal_t : 0;                                   
	     *onz_b = (normal_b > 0) ? *onz_b / normal_b : 0;

		/*

	     *onz_ne = (normal_ne > 0) ? *onz_ne / (normal_ne) : 0;
	     *onz_nw = (normal_nw > 0) ? *onz_nw / (normal_nw) : 0;
	     *onz_se = (normal_se > 0) ? *onz_se / (normal_se) : 0;
	     *onz_sw = (normal_sw > 0) ? *onz_sw / (normal_sw) : 0;
	     *onz_te = (normal_te > 0) ? *onz_te / (normal_te) : 0;
             *onz_tw = (normal_tw > 0) ? *ony_tw / (normal_tw) : 0;
	     *onz_tn = (normal_tn > 0) ? *ony_tn / (normal_tn) : 0;
	     *onz_ts = (normal_ts > 0) ? *ony_ts / (normal_ts) : 0;
	     *onz_be = (normal_be > 0) ? *onz_be / (normal_be) : 0;
             *onz_bw = (normal_bw > 0) ? *ony_bw / (normal_bw) : 0;
             *onz_bn = (normal_bn > 0) ? *ony_bn / (normal_bn) : 0;
             *onz_bs = (normal_bs > 0) ? *ony_bs / (normal_bs) : 0;
	     *onz_tne = (normal_tne > 0) ? *onz_tne / (normal_tne) : 0;
	     *onz_tnw = (normal_tnw > 0) ? *onz_tnw / (normal_tnw) : 0;
	     *onz_tsw = (normal_tsw > 0) ? *onz_tsw / normal_tsw : 0;
	     *onz_tse = (normal_tse > 0) ? *onz_tse / normal_tse : 0;
	     *onz_bne = (normal_bne > 0) ? *onz_bne / normal_bne : 0;
	     *onz_bnw = (normal_bnw > 0) ? *onz_bnw / normal_bnw : 0;
	     *onz_bse = (normal_bse > 0) ? *onz_bse / normal_bse : 0;
	     *onz_bsw = (normal_bsw > 0) ? *onz_bsw / normal_bsw : 0;
	     
	     

             *onx_e = (*onx_bse + *onx_bne + *onx_tse + *onx_bne + beta*(*onx_be + *onx_se + *onx_ne + *onx_te)+gamma*(*onx_e))/(4 + 4*beta + gamma);
             *onx_w = (*onx_bsw + *onx_bnw + *onx_tsw + *onx_bnw + beta*(*onx_bw + *onx_sw + *onx_nw + *onx_tw)+gamma*(*onx_w))/(4 + 4*beta + gamma);
	     *ony_s = (*ony_bsw + *ony_bse + *ony_tsw + *ony_tse + beta*(*ony_bs + *ony_se + *ony_sw + *ony_ts)+gamma*(*ony_s))/(4 + 4*beta + gamma); 
             *ony_n = (*ony_bnw + *ony_bne + *ony_tnw + *ony_tne + beta*(*ony_bn + *ony_ne + *ony_nw + *ony_tn)+gamma*(*ony_n))/(4 + 4*beta + gamma);
	     *onz_t = (*onz_tne + *onz_tnw + *onz_tsw + *onz_tse + beta*(*onz_tw + *onz_te + *onz_tn + *onz_ts)+gamma*(*onz_t))/(4 + 4*beta + gamma);
	     *onz_b = (*onz_bne + *onz_bnw + *onz_bsw + *onz_bse + beta*(*onz_bw + *onz_be + *onz_bn + *onz_bs)+gamma*(*onz_b))/(4 + 4*beta + gamma);
	     */

	     k_xx = 0.5 * (*onx_e - *onx_w);
             k_xy = 0.5 * (*onx_n - *onx_s);
	     k_xz = 0.5 * (*onx_t - *onx_b);
             
	     k_yx = 0.5 * (*ony_e - *ony_w);
             k_yy = 0.5 * (*ony_n - *ony_s);
             k_yz = 0.5 * (*ony_t - *ony_b);

	     k_zx = 0.5 * (*onz_e - *onz_w);
             k_zy = 0.5 * (*onz_n - *onz_s);
             k_zz = 0.5 * (*onz_t - *onz_b);

	     curv_x = (k_xx + k_xy + k_xz);
	     curv_y = (k_yx + k_yy + k_yz);
	     curv_z = (k_zx + k_zy + k_zz);

	     if (abs(curv_x) >= abs(curv_y) && abs(curv_x) >= abs(curv_z)){
		*curv = curv_x;}
	     
	     else if (abs(curv_y) >= abs(curv_x) && abs(curv_y) >= abs(curv_z)){
                *curv = curv_y;}

	     else if (abs(curv_z) >= abs(curv_y) && abs(curv_z) >= abs(curv_x)){
                *curv = curv_z;}

	     else {*curv = 0;}

             /* *curv = (k_x + k_y + k_z); */
          }
          else {
               *curv = 0;
          }

            nid ++;
            curv ++;
            onx ++;
            ony ++;
	    onz ++;            

	    onx_e ++;
            onx_w ++;
            onx_n ++;
            onx_s ++;
	    onx_t ++;
	    onx_b ++;            

		/*

	    onx_ne ++;
            onx_nw ++;
            onx_se ++;
            onx_sw ++;
	    onx_te ++;
	    onx_tw ++;
	    onx_tn ++;
	    onx_ts ++;
	    onx_be ++;
	    onx_bw ++;
	    onx_bn ++;
	    onx_bs ++;
	    onx_bne ++;
	    onx_bnw ++;
	    onx_bse ++;
	    onx_bsw ++;
	    onx_tne ++;
	    onx_tnw ++;
	    onx_tse ++;
	    onx_tsw ++;	    
            
		*/

	    ony_e ++;
            ony_w ++;
            ony_n ++;
            ony_s ++;
            ony_t ++;
	    ony_b ++;

		/*

	    ony_ne ++;
            ony_nw ++;
            ony_se ++;
            ony_sw ++;
	    ony_te ++;
	    ony_tw ++;
	    ony_tn ++;
	    ony_ts ++;
	    ony_be ++;
	    ony_bw ++;
	    ony_bn ++;
	    ony_bs ++;
	    ony_bne ++;
            ony_bnw ++;
            ony_bse ++;
            ony_bsw ++;
            ony_tne ++;
            ony_tnw ++;
            ony_tse ++;
            ony_tsw ++;

		*/

	    onz_e ++;
            onz_w ++;
            onz_n ++;
            onz_s ++;
            onz_t ++;
            onz_b ++;

		/*

            onz_ne ++;
            onz_nw ++;
            onz_se ++;
            onz_sw ++;
            onz_te ++;
            onz_tw ++;
            onz_tn ++;
            onz_ts ++;
            onz_be ++;
            onz_bw ++;
            onz_bn ++;
            onz_bs ++;	
	    onz_bne ++;
            onz_bnw ++;
            onz_bse ++;
            onz_bsw ++;
            onz_tne ++;
            onz_tnw ++;
            onz_tse ++;
            onx_tsw ++;

*/
            /* add in onz and rest of onx, ony */

          } /*end of I loop */
            onx += 2;
            ony += 2;
            onz += 2;

	    onx_e += 2;
            onx_w += 2;
            onx_n += 2;
            onx_s += 2;
            onx_t += 2;
            onx_b += 2;

		/*

            onx_ne += 2;
            onx_nw += 2;
            onx_se += 2;
            onx_sw += 2;
            onx_te += 2;
            onx_tw += 2;
            onx_tn += 2;
            onx_ts += 2;
            onx_be += 2;
            onx_bw += 2;
            onx_bn += 2;
            onx_bs += 2;
	    onx_bne += 2;
            onx_bnw += 2;
            onx_bse += 2;
            onx_bsw += 2;
            onx_tne += 2;
            onx_tnw += 2;
            onx_tse += 2;
            onx_tsw += 2;

		*/

            ony_e += 2;
            ony_w += 2;
            ony_n += 2;
            ony_s += 2;
            ony_t += 2;
            ony_b += 2;

		/*

            ony_ne += 2;
            ony_nw += 2;
            ony_se += 2;
            ony_sw += 2;
            ony_te += 2;
            ony_tw += 2;
            ony_tn += 2;
            ony_ts += 2;
            ony_be += 2;
            ony_bw += 2;
            ony_bn += 2;
            ony_bs += 2;
	    ony_bne += 2;
            ony_bnw += 2;
            ony_bse += 2;
            ony_bsw += 2;
            ony_tne += 2;
            ony_tnw += 2;
            ony_tse += 2;
            ony_tsw += 2;

		*/

            onz_e += 2;
            onz_w += 2;
            onz_n += 2;
            onz_s += 2;
            onz_t += 2;
            onz_b += 2;

		/*

            onz_ne += 2; 
            onz_nw += 2;
            onz_se += 2;
            onz_sw += 2;
            onz_te += 2;
            onz_tw += 2;
            onz_tn += 2;
            onz_ts += 2;
            onz_be += 2;
            onz_bw += 2;
            onz_bn += 2;
            onz_bs += 2;
	    onz_bne += 2;
            onz_bnw += 2;
            onz_bse += 2;
            onz_bsw += 2;
            onz_tne += 2;
            onz_tnw += 2;
            onz_tse += 2;
            onz_tsw += 2;
		
		*/

        } /*end of J loop */
            onx += skip;
            ony += skip;
 	    onz += skip;

	    onx_e += skip;
            onx_w += skip;
            onx_n += skip;
            onx_s += skip;
            onx_t += skip;
            onx_b += skip;

		/*

            onx_ne += skip;
            onx_nw += skip;
            onx_se += skip;
            onx_sw += skip;
            onx_te += skip;
            onx_tw += skip;
            onx_tn += skip;
            onx_ts += skip;
            onx_be += skip;
            onx_bw += skip;
            onx_bn += skip;
            onx_bs += skip;
	    onx_bne += skip;
            onx_bnw += skip;
            onx_bse += skip;
            onx_bsw += skip;
            onx_tne += skip;
            onx_tnw += skip;
            onx_tse += skip;
            onx_tsw += skip;

		*/

            ony_e += skip;
            ony_w += skip;
            ony_n += skip;
            ony_s += skip;
            ony_t += skip;
            ony_b += skip;

		/*

            ony_ne += skip;
            ony_nw += skip;
            ony_se += skip;
            ony_sw += skip;
            ony_te += skip;
            ony_tw += skip;
            ony_tn += skip;
            ony_ts += skip;
            ony_be += skip;
            ony_bw += skip;
            ony_bn += skip;
            ony_bs += skip;
	    ony_bne += skip;
            ony_bnw += skip;
            ony_bse += skip;
            ony_bsw += skip;
            ony_tne += skip;
            ony_tnw += skip;
            ony_tse += skip;
            ony_tsw += skip;

		*/

            onz_e += skip;
            onz_w += skip;
            onz_n += skip;
            onz_s += skip;
            onz_t += skip;
            onz_b += skip;

		/*

            onz_ne += skip;
            onz_nw += skip;
            onz_se += skip;
            onz_sw += skip;
            onz_te += skip;
            onz_tw += skip;
            onz_tn += skip;
            onz_ts += skip;
            onz_be += skip;
            onz_bw += skip;
            onz_bn += skip;
            onz_bs += skip;
	    onz_bne += skip;
            onz_bnw += skip;
            onz_bse += skip;
            onz_bsw += skip;
            onz_tne += skip;
            onz_tnw += skip;
            onz_tse += skip;
            onz_tsw += skip;

		*/

       } /*end of K loop */
   return (1);
 }

char const *rcs_id_curvature_3D_c()
{
   static char const rcsid[] = "$Id: curvature_3D.c 1356 2008-08-18 13:41:15Z  $";
   return(rcsid);
}
