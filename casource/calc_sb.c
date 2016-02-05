
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

/**\file calc_sb.c */
/****************************************************************/
/* Subroutine to perform one timestep on a subblock             */
/****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>/* for memset prototype */

#include "machine.h"
#include "blocks.h"
#include "umat_matrix.h"
#include "writeblocks.h"
#include "read_sb.h"
#include "find_max.h"
#include "interp.h"

#include "sb_head.h"
#include "sb_diffuse.h"
#include "sb_decentred_step.h"
#include "curvature.h"

#ifdef CHECK_GAS
extern CA_FLOAT checkgas(BB_struct *bp,int callflag) ;
#endif



extern int grow_octahedron_poly (BB_struct * bp, int sbnum);

/** 
* Function calc_sb: Perform all necessary calculation for one microstep on 
*   one subblock.
*
* @callgraph
* @callergraph
*/
int calc_sb (BB_struct * bp, int sbnum)
{
   int      errors = 0;
   int      i;    /* tmp counters */
  CA_FLOAT *ofs, *nfs, *ofs_pad;
   int      *ogr, *ngr;
   CA_FLOAT coor[3];   /* tmp CA_FLOAT var.*/
   int      **gv_p; /* grain block value array pointer pointer :-/ */
   SB_struct * sp;
  int iphs, iphs_tot;
  CA_FLOAT *ofs_poly[NPHAMAX],*nfs_poly[NPHAMAX];

   #ifdef CHECK_GAS
      static CA_FLOAT diffgas0=0,diffgas1=0,change=0;

      diffgas0 = checkgas(bp,5);
      change = diffgas0 - diffgas1;
      if (ABS(change) > 1e-5){
         fprintf(stderr,"WARNING: sb_pore: gas change limit exceeded! %.5g\n",change);
      }
   #endif

  iphs_tot = bp->ctrl->NUM_PHS;

   sp = bp->sb[sbnum];
   gv_p = bp->gr_array;
   ofs = bp->ftmp_one;
   nfs = sp->c_fs;
   ogr = bp->itmp_one;
   ngr = sp->gr;

  if (bp->ctrl->diffuse_alloy_poly == TRUE){
    for (iphs = 0; iphs < iphs_tot; iphs++) {
      ofs_poly[iphs]=bp->ftmp_one_poly[iphs];
      nfs_poly[iphs]=sp->c_fs_poly[iphs];
    }
  }


   coor[0] = coor[1] = coor[2] = 1.0;

   /************************************************/
   /* Calculate the new temperature of the sb...   */
   /************************************************/
   sb_temp_calc(bp, sbnum);
   /************************************************/
   /* Copy the f_s and gr  arrays...           */
   /************************************************/
   if (bp->ctrl->scheil == TRUE) {
      /* (flag, to, from, bp, block_values,sbnum) */
      fcopy_matrix(PAD, bp->ftmp_three, sp->sch_fs, bp, bp->sch_fs_values->block_array, sbnum);
  }

  /*end of schiel test */
   /* The fraction solid is copied into a padded buffer array */
   /* (flag, to, from, bp) */
   fcopy_matrix(PAD, ofs, nfs, bp, bp->c_fs_values->block_array, sbnum);

  if (bp->ctrl->diffuse_alloy_poly == TRUE){
    for (iphs = 0; iphs < iphs_tot; iphs++) {  
      fcopy_matrix (PAD,ofs_poly[iphs],nfs_poly[iphs], bp, bp->c_fs_values[iphs].block_array, sbnum);
    }
   }
   
   /************************************************/
   /* If still liquid call nucleation and growth   */
   /************************************************/
   if (!bp->sb[sbnum]->done) {

      /*******************************************/
      /* Nucleate ...                            */
      /* This nucleation routine gets called     */
      /* for those modes where there is not      */
      /* cell-to-cell variation.                 */
      /*                                         */
      /*******************************************/

      /* this is set in combine_ctrl routine */
      if (!(bp->ctrl->use_cell_nuc)){
       sb_nuc(bp, sbnum);
      }
       /************************************************/
       /* Grow...                                      */
       /* Call the CA subroutine for one step...       */
       /*                                              */
       /************************************************/

      /* Wei's decentred square / octahedron method */
      if (bp->ctrl->decentred_octahedron == TRUE) {  /* by Wei WANG 15-07-02*/

      /* this is set in combine_ctrl routine */
      if ((bp->ctrl->use_cell_nuc) && (bp->ctrl->diffuse_alloy_poly == FALSE)) {
          cell_nucleation (bp, sbnum);                         /*cell nucleation*/
	}
        
      if ((bp->ctrl->use_cell_nuc) && (bp->ctrl->diffuse_alloy_poly == TRUE)) {
        cell_nucleation_poly (bp, sbnum);       /*cell nucleation */
	}
 
      if ((bp->ctrl->diffuse_alloy == TRUE) && (bp->ctrl->diffuse_alloy_poly == FALSE)) {
         icopy_matrix(PAD, ogr, ngr, bp, gv_p, sbnum);
         /* need to use old GR stored in temporary buffer */
         cell_index(bp, sbnum);                                 /* cell index */
 
	cell_temp(bp, sbnum);                                /* cell temperature */

         /* need to use old FS stored in temporary buffer 
            new Curvature changed in SB */
         
	 if (bp->ctrl->curvature_3D){
	    interface_normal_3D(bp,sbnum);
	    interface_curvature_3D(bp, sbnum);
	 } else if (bp->ctrl->curvature_2D == 1){
	    interface_normal(bp, sbnum); 
       interface_curvature(bp, sbnum);
	 } else if (bp->ctrl->curvature_2D == 2){
       surface_curvature(bp, sbnum); /* Nastac interface curvature */
    } 

         /* alloy diffusion */
         /* need to use old FS stored in temporary buffer
            old C_L is copied to temporary buffer
            old C_E is copied to temporary buffer, new C_E is stored in SB */ 
         sb_diffuse_alloy_decentred(bp, sbnum);                           /* diffuse alloy */

         /* need to use old C_E stored in temporary buffer and new C_E in SB
            need to use old Curvature in SB
            new FS is changed in ftmp_one
            new C_L is changed in SB
            new GR assoicated with melting back is changed in SB */
         fs_change_diffuse(bp, sbnum);                          /* fs change with diffusion */
 
         /* need to use new FS in ftmp_one
            new half diagonal of octahedron is stored in SB */
         grow_octahedron (bp, sbnum);                    /* grow octahedron */

         /* capture new cells */
         /* need to use GR, FS stoed in temoprary buffer
            new GR associated with captures are stored in SB
            new FS associated with captures are stored in ftmp_one
            old informations of octahedron are copied to temporary buffers, new ones are stored in SB */
         capture_octahedron_diffuse(bp, sbnum);                 /* octahedron capture cell with diffusion*/
      }

      else if (bp->ctrl->diffuse_alloy_poly == TRUE) {  /* wei wang's algorithm extension to multicomponent systems by thuinet */

        /* need to use old GR stored in temporary buffer */
        cell_index (bp, sbnum); /* cell index */

        cell_temp (bp, sbnum);  /* cell temperature */

        /* need to use old FS stored in temporary buffer
           new Curvature changed in SB */
        /* cannot use existing routine with polycomponent since the fraction-solid does not have the correct meaning */
        /**todo: modify curvature to use a meaningful value in the polycomponent case */
        /*surface_curvature (bp, sbnum);*/  /* interface curvature */
        if ((bp->ctrl->curvature_3D !=0) || (bp->ctrl->curvature_2D != 0 )){
            fprintf(stderr,"ERROR:%s: Curvature not yet implemented for poly-component option! \n",__func__);
            exit(0);
        }

        /* alloy diffusion */
        /* need to use old FS stored in temporary buffer
           old C_L is copied to temporary buffer
           old C_E is copied to temporary buffer, new C_E is stored in SB */
        sb_diffuse_alloy_decentred_poly (bp, sbnum);    /*diffuse alloy in multicomponent systems */

        /* need to use old C_E stored in temporary buffer and new C_E in SB
           need to use old Curvature in SB
           new FS is changed in ftmp_one
           new C_L is changed in SB
           new GR assoicated with melting back is changed in SB */
        fs_change_diffuse_poly (bp, sbnum);     /* fs change with diffusion in multicomponents systems */

        /* need to use new FS in ftmp_one
           new half diagonal of octahedron is stored in SB */
        grow_octahedron_poly (bp, sbnum);       /* grow octahedron */

        /* capture new cells */
        /* need to use GR, FS stoed in temoprary buffer
           new GR associated with captures are stored in SB
           new FS associated with captures are stored in ftmp_one
           old informations of octahedron are copied to temporary buffers, new ones are stored in SB */
        capture_octahedron_diffuse_poly (bp, sbnum);    /* octahedron capture cell with diffusion */
	} else { /*No diffusion */
         fs_change_nodiffuse (bp, sbnum);
         grow_octahedron (bp, sbnum);
         capture_octahedron (bp, sbnum);
	}

         /*******************************************/
         /* Diffuse gas                             */
         /*******************************************/
      if (bp->ctrl->diffuse == TRUE){
        if (bp->ctrl->diffuse_alloy_poly == FALSE) 
          errors += sb_diffuse_gas(bp, sbnum);
        if (bp->ctrl->diffuse_alloy_poly == TRUE)
          errors += sb_diffuse_gas_poly (bp, sbnum);
      }

      /******************************************/
      /* Call the pore calculation if needed    */
      /******************************************/

      if (bp->ctrl->pore == TRUE)
        errors += sb_pore (bp, sbnum);

      fcopy_mat_back (PAD, nfs, ofs, bp, bp->c_fs_values->block_array, sbnum);

      if (bp->ctrl->diffuse_alloy_poly == TRUE){
        for (iphs = 0; iphs < iphs_tot; iphs++) {
          fcopy_mat_back(PAD,nfs_poly[iphs],ofs_poly[iphs], bp, bp->c_fs_values[iphs].block_array, sbnum);
	       }
	  }
	  


      } else {  /* Robert's algorithm  -- 1-neighbour capture*/
       /*ensure new grains are included*/
       icopy_matrix(PAD, ogr, ngr, bp, gv_p, sbnum);

       /**********************************/
       /*   call the ca step algorithm   */
       /**********************************/

       /* fraction solid is updated in the buffer ftmp_one */

       /****************************************************/
       /*                                                  */
       /*             DO THE CA STEP                       */
       /*                                                  */
       /****************************************************/

       sb_umat_step(bp, sbnum);

       /* after sb_umat_step, the NEW fraction solid is in the subblock fs array */
       /* and the OLD fraction solid is still in the bp->ftmp_one array. This is */
       /* needed in order for the diffusion (non-decentered method and gas diffusion)*/
       /* to calculate the partition of the dissolved species */

       /* this is a potential problem for the multi block case -- where does the old frction */
       /* solid for another block go to ? possibly the whole old big block needs to be buffered */

         /*******************************************/
         /* Diffuse alloy                           */
         /*******************************************/
      if ((bp->ctrl->diffuse_alloy == TRUE)&&(bp->ctrl->diffuse_alloy_multi==FALSE)){
         errors += sb_diffuse_alloy(bp, sbnum);
      }


      /* needs to have OLD fs in subblock array, NEW fs in buffer! */
      if (bp->ctrl->diffuse == TRUE)
         errors += sb_diffuse_gas(bp, sbnum);

          /******************************************/
          /* Call the pore calculation if needed    */
          /******************************************/
      if (bp->ctrl->pore == TRUE)
        errors += sb_pore (bp, sbnum);

         /*******************************************/
         /* update scheil concentration             */
         /*******************************************/
      if (bp->ctrl->scheil == TRUE) {
        fcopy_mat_back (PAD, sp->sch_fs, bp->ftmp_three, bp, bp->sch_fs_values->block_array, sbnum);
        /* (flag, to, from, bp) */
      }

      /*end of schiel test */
      /*************************************/
      /* at this point the old f_s is lost */
      /*************************************/
      fcopy_mat_back(PAD, nfs, ofs, bp, bp->c_fs_values->block_array, sbnum);
      icopy_mat_back(PAD, ngr, ogr, bp, bp->gr_array, sbnum);

      i = (bp->nc[0] + 2) *  (bp->nc[1] + 2) *  (bp->nc[2] + 2);
      if (memset(bp->itmp_one,0,i*sizeof(int)) == NULL){
         fprintf(stderr,"ERROR: calc_sb: memset failed! Debug me! \n");
         exit(25);
      }
      
   } /* end Robert's algorithm*/

  }

  /*end of test for "done" condition */
   /* window moving by Wei Wang on 04-11-02 */
   if (bp->ctrl->window_moving == TRUE) {
     bp->window_disp += bp->window_velo * bp->delt;
/*THUINET 03/05*/
    bp->window_disp2 += bp->window_velo * bp->delt;
/*END THUINET 03/05*/
     if (bp->window_disp >= bp->size_c[0]) {
       window_move (bp,sbnum);    /* in sb_umat_step.c */
       bp->window_disp -= bp->size_c[0];
     }
   }
   return (errors);
} /* end of calc_sb subroutine */

char const *rcs_id_calc_sb_c()
{
  static char const rcsid[] = "$Id: calc_sb.c 1373 2008-08-27 20:51:52Z  $";

   return(rcsid);
}

/*
*/
