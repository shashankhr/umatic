
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
#include "machine.h"
#include "blocks.h"
#include "read_sb.h"
#include "SurCellRoutines.h"
extern int      sb_nuc(BB_struct *bp, int sbnum);

/* functions used from sb_temp_calc.c */
extern int      sb_temp_calc(BB_struct *bp, int sbnum);

/* subroutines later in the file ... */
extern int      init_c_elm(BB_struct *, int);
extern int      init_c_elm_solo(BB_struct *, int);
extern int      init_nuc_thresh(BB_struct *bp,int sbnum);

int alloc_sb (BB_struct * bp, int sbnum)
{
      int      errors = 0;
  int i, j, ele_1,npha, nc, nxny;
      int      oldstep;

/*THUINET 05/05*/
  int iphs, iphs_tot;

/*FIN THUINET 05/05*/
      SB_struct * sp;
   MultiS_struct *ms;

/*THUINET 05/05*/
  iphs_tot = (bp->ctrl->NUM_PHS);
/*FIN THUINET 05/05*/

      sp = bp->sb[sbnum];
      nc = bp->ncsb;
      nxny = bp->nc[0]*bp->nc[1];
      ms = &(bp->MultiSvals);

#  ifdef USE_ELM
      /* malloc array to hold the FE element that the cells belong to  */
      if (!(sp->c_elm = (int *) malloc(nc * sizeof(int)))) {
          fprintf(stderr, "ERROR: alloc_sb: SB cell element  array malloc failed\n");
          return(1);
      }
      bp->c_elm_array[sbnum] = sp->c_elm;
#  endif /*USE_ELM*/
      /* malloc array to hold fs of cells and set = 0 */
      if (!(sp->c_fs = (CA_FLOAT * ) calloc(nc, sizeof(CA_FLOAT)))) {
          fprintf(stderr, "ERROR: alloc_sb: SB cell fs array malloc failed\n");
          return(1);
      }
      /* set value pointer array */
      bp->c_fs_values->block_array[sbnum] = sp->c_fs;

/********************only for multi component***************/

  /* end of multi diff test */
/*******************************************************/
  if (bp->ctrl->diffuse_alloy_poly == TRUE) {
     ele_1 = (bp->ctrl->NUM_COMP - 1);
     npha = bp->ctrl->NUM_PHS;

      if (!(sp->nat_cell = (int *) calloc(nc, sizeof(int)))) {
          fprintf(stderr, "ERROR: alloc_sb: SB cell gr array malloc failed\n");
          return(1);
        }

    if (!(sp->nat_grain = (int *) calloc (nc, sizeof (int)))) {
      fprintf (stderr, "ERROR: alloc_sb: SB cell gr array malloc failed\n");
      return (1);
    }


  }else{
     ele_1 = 1;
     npha = 1;
  }




    for (iphs = 0; iphs < npha ; iphs++) {

      if (!(sp->c_fs_poly[iphs] = (CA_FLOAT *) calloc (nc, sizeof (CA_FLOAT)))) {
        fprintf (stderr, "ERROR: %s: malloc failed at %s\n",__func__,__LINE__);
        return (1);
      }
    }

  /*allocation and initialisation for multi-component system (THUINET) */
  /* now used for pure binary by setting eel_1 to the appropriate value */
    for (j = 0; j < ele_1; j++) {

      if (!(sp->c_sol_poly[j] = (CA_FLOAT *) calloc (nc, sizeof (CA_FLOAT)))) {
        fprintf (stderr, "ERROR: alloc_sb: SB cell alloysol array malloc failed for multi\n");
        return (1);
      }

      if (!(sp->c_eqv_poly[j] = (CA_FLOAT *) calloc (nc, sizeof (CA_FLOAT)))) {
        fprintf (stderr, "ERROR: alloc_sb: SB cell equiv. sol array malloc failed for multi\n");
        return (1);
      }

      /* set concentration to initial value */
      for (i = 0; i < nc; i++) {
        *(sp->c_sol_poly[j] + i) = bp->mprops.alloyprops[j].Cinit;
        *(sp->c_eqv_poly[j] + i) = bp->mprops.alloyprops[j].Cinit;
      }
      /* set value pointer array
         bp->c_sol_alloy_values->block_array[sbnum] = sp->c_sol_alloy; */
    }


  /* end of multi component test (THUINET) */

      /* malloc array to hold schiel_fs of cells and set = 0 */
      if (bp->ctrl->scheil == TRUE) {
          if (!(sp->sch_fs = (CA_FLOAT * ) calloc(nc, sizeof(CA_FLOAT)))) {
              fprintf(stderr, "ERROR: alloc_sb: SB sch_fs array malloc failed\n");
              return(1);
          }
          /* set value pointer array */
          bp->sch_fs_values->block_array[sbnum] = sp->sch_fs;
      }

      /* malloc array for cell index & set = 0 */ /*by Wei WANG 11-07-02*/
      if (!(sp->index = (int *) calloc(nc, sizeof(int)))) {
         fprintf(stderr, "ERROR: alloc_sb: SB cell index array malloc failed\n");
        return(1);
      }
      /* malloc array to hold grain # of cells and set = 0 */
      if (!(sp->gr = (int *) calloc(nc, sizeof(int)))) {
          fprintf(stderr, "ERROR: alloc_sb: SB cell gr array malloc failed\n");
          return(1);
      }
      /* set value pointer array */
      bp->gr_array[sbnum] = sp->gr;

      /* malloc array for cell temperature */ /*by Wei Wang on 11-07-02*/
      if (!(sp->c_temp = (CA_FLOAT *) calloc(nc, sizeof(CA_FLOAT)))) {
         fprintf(stderr, "ERROR: alloc_sb: SB cell temperature array malloc failed\n");
         return(1);
      }

      /* malloc array for next cell temperarue, for 2d temp. dist. */ /*by Robert Atwood*/
      if(bp->ctrl->fgrid_input){
         if (!(sp->c_fg_temp = (CA_FLOAT *) calloc(nxny, sizeof(CA_FLOAT)))) {
            fprintf(stderr, "ERROR: alloc_sb: SB cell temp_next array malloc failed\n");
            return(1);
         }
         if (!(sp->c_fg_temp_next = (CA_FLOAT *) calloc(nxny, sizeof(CA_FLOAT)))) {
            fprintf(stderr, "ERROR: alloc_sb: SB cell temp_next array malloc failed\n");
            return(1);
         }
      }

      /* malloc array for curvature */ /*by xly 2004/08/09*/
      /* protected by options rca 2007 01 15 */
      /**todo: finetune the options to avoid excessive memory allocation when the options are not used*/
      if ((bp->ctrl->curvature_3D != 0 ) || (bp->ctrl->curvature_2D != 0)) {
        if (!(sp->curv = (CA_FLOAT *) calloc(nc+2, sizeof(CA_FLOAT)))) {

           fprintf(stderr, "ERROR: alloc_sb: SB array of ptrs to curv malloc failed\n");
           return (1);
        }   
         if (!(sp->norm_x = (CA_FLOAT *) calloc(nc+2, sizeof(CA_FLOAT)))) {
           fprintf(stderr, "ERROR: alloc_sb: SB array of ptrs to norm_x malloc failed\n");
           return (1);   
        }
         if (!(sp->norm_y = (CA_FLOAT *) calloc(nc+2, sizeof(CA_FLOAT)))) {
           fprintf(stderr, "ERROR: alloc_sb: SB array of ptrs to norm_y malloc failed\n");
           return (1);   
        } 
      	if (!(sp->norm_z = (CA_FLOAT *) calloc(nc+2, sizeof(CA_FLOAT)))) {
           fprintf(stderr, "ERROR: alloc_sb: SB array of ptrs to norm_z malloc failed\n");
           return (1);
        }
      }

      /* malloc array to hold solute conc of cells and initialise */
      if (bp->ctrl->diffuse == TRUE) {/* GAS solute array */
          sp->Svals[GAS].Cinit = bp->mprops.gasprops.Cinit;/*all the same for now!*/
          if (!(sp->c_sol = (CA_FLOAT * ) calloc(nc, sizeof(CA_FLOAT)))) {
              fprintf(stderr, "ERROR: alloc_sb: SB cell sol array malloc failed\n");
              return(1);
          }
          /* set concentraion to initial value */
          for (i = 0; i < nc; i++)
/** \todo  compartmentalize the initialization seperate from allocaton -- merge-xly */
              *(sp->c_sol + i) = sp->Svals[GAS].Cinit;
          /* set value pointer array */
          bp->c_sol_values->block_array[sbnum] = sp->c_sol;
      }

  if (((bp->ctrl->diffuse_alloy == TRUE) || (bp->ctrl->particle == TRUE)) && (bp->ctrl->diffuse_alloy_multi == FALSE)&&(bp->ctrl->diffuse_alloy_poly == FALSE)) {        /*ALLOY solute array */
    sp->Svals[ALLOY].Cinit = bp->mprops.alloyprops[0].Cinit;       /*all the same for now! */
         if (!(sp->c_sol_alloy = (CA_FLOAT * ) calloc(nc, sizeof(CA_FLOAT)))) {
            fprintf(stderr, "ERROR: alloc_sb: SB cell alloysol array malloc failed\n");
            return(1);
         }
         /* malloc array for C_E*/ /*by Wei WANG on 11-07-02*/
         if (bp->ctrl->decentred_octahedron == TRUE){
            if (!(sp->c_eqv_alloy = (CA_FLOAT * ) calloc(nc, sizeof(CA_FLOAT)))) {
               fprintf(stderr, "ERROR: alloc_sb: SB cell equiv. sol array malloc failed\n");
               return(1);
            }
         }

         /* set concentraion to initial value */
         for (i = 0; i < nc; i++) {
            *(sp->c_sol_alloy + i) = sp->Svals[ALLOY].Cinit;
            if (bp->ctrl->decentred_octahedron == TRUE)
               *(sp->c_eqv_alloy + i) = sp->Svals[ALLOY].Cinit; /*by Wei WANG 11-07-02*/
         }
         /* set value pointer array */
         bp->c_sol_alloy_values->block_array[sbnum] = sp->c_sol_alloy;
      }

      /* malloc array for decented octahedron */ /*by Wei Wang 11-07-02*/
      if (bp->ctrl->decentred_octahedron == TRUE) {        
        if (!(sp->dc_d = (CA_FLOAT *) calloc(nc, sizeof(CA_FLOAT)))) {
           fprintf(stderr, "ERROR: alloc_sb: SB array of ptrs to dc_d malloc failed\n");
        }
        if (!(sp->dc_x = (CA_FLOAT *) calloc(nc, sizeof(CA_FLOAT)))) {
           fprintf(stderr, "ERROR: alloc_sb: SB array of ptrs to dc_x malloc failed\n");
        }
        if (!(sp->dc_y = (CA_FLOAT *) calloc(nc, sizeof(CA_FLOAT)))) {
           fprintf(stderr, "ERROR: alloc_sb: SB array of ptrs to dc_y malloc failed\n");
        }
        if (!(sp->dc_z = (CA_FLOAT *) calloc(nc, sizeof(CA_FLOAT)))) {
           fprintf(stderr, "ERROR: alloc_sb: SB array of ptrs to dc_z malloc failed\n");
        }
      }

   /*allocate and initialise the multi component array*/

  #ifdef CHIRAZI_MULTI
   if (bp->ctrl->diffuse_alloy_multi == TRUE ) {

      sp->c_sol_alloy_multi=(CA_FLOAT **)calloc(ele_1,sizeof(CA_FLOAT *));
      sp->c_sol_tot_multi=(CA_FLOAT **)calloc(ele_1,sizeof(CA_FLOAT *));
      sp->c_sol_tot_old_multi=(CA_FLOAT **)calloc(ele_1,sizeof(CA_FLOAT *));
      sp->c_sol_n_eut_multi=(CA_FLOAT **)calloc(ele_1,sizeof(CA_FLOAT *));
      sp->c_sol_n_eut_old_multi=(CA_FLOAT **)calloc(ele_1,sizeof(CA_FLOAT *));

      for(j=0;j<ele_1;j++){

      sp->c_sol_alloy_multi[j]=(CA_FLOAT *)calloc(nc,sizeof(CA_FLOAT));
      sp->c_sol_tot_multi[j]=(CA_FLOAT *)calloc(nc,sizeof(CA_FLOAT));
      sp->c_sol_tot_old_multi[j]=(CA_FLOAT *)calloc(nc,sizeof(CA_FLOAT));
      sp->c_sol_n_eut_multi[j]=(CA_FLOAT *)calloc(nc,sizeof(CA_FLOAT));
      sp->c_sol_n_eut_old_multi[j]=(CA_FLOAT *)calloc(nc,sizeof(CA_FLOAT));

      /* set concentraion to initial value */
      for (i = 0; i < nc; i++){
         sp->c_sol_alloy_multi[j][i]=ms->Cinit_multi[j];
         sp->c_sol_tot_multi[j][i]=ms->Cinit_multi[j];
         sp->c_sol_tot_old_multi[j][i]=ms->Cinit_multi[j];
       }
      /* set value pointer array */
      bp->multi_conc[j]->block_array[sbnum] = sp->c_sol_alloy_multi[j];
      bp->multi_conc_tot[j]->block_array[sbnum] = sp->c_sol_tot_multi[j];
      bp->multi_conc_tot_old[j]->block_array[sbnum] = sp->c_sol_tot_old_multi[j];
      bp->multi_conc_n_eut[j]->block_array[sbnum] = sp->c_sol_n_eut_multi[j];
      bp->multi_conc_n_eut_old[j]->block_array[sbnum] = sp->c_sol_n_eut_old_multi[j];
    }
  }
  #endif /* CHIRAZI_MULTI */
  /* end of multi diff test */
   /************************************************/
      /*allocate nuc threshold array*/
      if ((bp->ctrl->block_nuc == TRUE )) { /*nuc threshold array*/
          if (!(sp->c_nuc_thresh = malloc(nc * sizeof(CA_FLOAT)))) {
              fprintf(stderr, "ERROR: alloc_sb: SB cell nuc thresh. array malloc failed, %i\n",nc*sizeof(CA_FLOAT));
              return(1);
          }
/*THUINET 05/05*/
      if (!(sp->nat_sol_site = malloc(nc * sizeof(CA_FLOAT)))) {
          fprintf(stderr, "ERROR: alloc_sb: SB cell nat_sol_site array malloc failed, %i\n",nc*sizeof(CA_FLOAT));
          return(1);
        }
/*END THUINET 05/05*/

        }
      /* allocate the initial structure to hold mould/casting interface cell location */
      alloc_SurCell(&(sp->surface));
      return(errors);
}/*end of alloc_sb subroutine*/

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file */
char const *rcs_id_alloc_sb_c ()
{
  static char const rcsid[] = "$Id: alloc_sb.c 1373 2008-08-27 20:51:52Z  $";

   return(rcsid);
   }

/* end of rcs_id_subroutine */
/*
*/
