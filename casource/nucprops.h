
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
/* nucprops.h                                                   */
/* Header file defining nucleation related structures.          */
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Wed Jul  1 18:38:31 bst 1998                                 */
/****************************************************************/
/*RCS Id:$Id: nucprops.h 1341 2008-07-23 15:23:30Z  $*/

#ifndef NUCPROPS_H
#define NUCPROPS_H
#define NUCPROPSREV "nucprops.h $Revision: 1341 $"

#define MAX_NSP 50
#define MAX_NUCAREA 150
#define MAX_NUC_PERTURB 20
#define NUC_PERTURB_F 3

/****************************************************************/
/* Structure to define a single clasical nucleation curve.      */
/* One structure must be created for each particle type.        */
/****************************************************************/
typedef struct nucleus {
   int fin;             /* TRUE if all nuclei consumed         */
   CA_FLOAT CA_A;           	/* constant "A"				*/
   CA_FLOAT B;           	/* constant "B"				*/
   CA_FLOAT ftheta;        /* f(theta)             		*/
   CA_FLOAT rad;	         /* average size				*/
   CA_FLOAT Nrate;	    	/* current nucleation rate 		*/
   CA_FLOAT Nmax;	      	/* max number of nuclei of this type    */
   CA_FLOAT Ntot;	      	/* total num. nucleated 		*/
} Nucleus_str;

/**************
**************************************************/
/* Structure to define an area for fixed het. nucleation.       */
/****************************************************************/
typedef struct perturb_nuc {
   CA_FLOAT v_f[NUC_PERTURB_F];  /* 0: magnitude of perturbation         */
                        /* 1: time perturb starts               */
                        /* 2: time perturbation lasts           */
} Nuc_Perturb_str;

#ifdef NEW_NUC
   /*
      The structure for storing the list of potential nuclei for the grains.
   */
   typedef struct nuc_list {
      int Cellnum ;
      int Cellxyz[3];
      CA_FLOAT thresh;
      struct nuc_list *next;
   } NUCLIST_str;

#endif /*NEW_NUC*/

 
/****************************************************************/
/* Structure to define an area for fixed het. nucleation.       */
/****************************************************************/
typedef struct nuc_area {
   CA_FLOAT dens;         	/* nucl. dens. in #/m^3                */
   int sbnum;	        /* sb that area is in                   */
   int lleft[3];	/* lower lefthand corner of area [cells]*/
   int uright[3];	/* upper righthand corner of area[cells]*/
} NucArea_struct;

/*******************************************/
/*                                         */
/*******************************************/
typedef struct nucdiststruc {
   int nuc_dist_func_type;
   int n_U_bins;          /* Number of nUc-activity bins           */
   int n_Tund_bins;       /* Number of Tundercooling bins          */
   CA_FLOAT Tund_incr;       /* increment of undercooling             */
   CA_FLOAT Unuc_incr;       /* increment of U                        */
   CA_FLOAT Tund_factor;     /* factor to multiply Tund for bin number*/
   CA_FLOAT * Number_U;      /* Array of number in each U bin         */
   CA_FLOAT * Nuctable;      /* Array of nuc rate by undercooling     */
   char ** Nuc_func_id;   /* Store the names of the functions for output */
   CA_FLOAT (**Nuc_func_table)(CA_FLOAT Unuc,CA_FLOAT par_one,CA_FLOAT par_two,CA_FLOAT par_three);
}Nuc_dist_str;


/****************************************************************/
/* Nucleation properties structures, Nuc_str contains           */
/* constans used for calculation of nucleation and              */
/* do not vary from one subblock to the next;                   */
/****************************************************************/
typedef struct nucprops  {
   int nmodel;          /* type of nucleation model:            */
                        /*    N_RAPPAZ - rappaz model           */
                        /*    N_HETERO - heterogenous model     */
                        /* etc. see nucprops.h                  */
   int nhet;		/* number of different het. nuclei      */
   Nucleus_str nuc[MAX_NSP];    /* ptr to array of heterogeneous nucl.  */
   int nareanuc;	/* number of different area nucl. spec. */
   NucArea_struct nap[MAX_NUCAREA];    /* ptr to array of area nuc str.  */
   int ngr;		/* number of grains			*/
   CA_FLOAT gd_max;        /* max. grain density per [cm^3]        */
   CA_FLOAT gd_max_surf;        /* max. grain density per [cm^3]        */
   CA_FLOAT gd_max_beut;        /* max. grain density for binary eutectic nucleation sites per [cm^3] */
   CA_FLOAT gd_max_teut;        /* max. grain density for ternary eutectic nucleation sites per [cm^3] */

/*used in the multiphase extension of the CAFD model THUINET  04/05 */
   CA_FLOAT gd_max_poly[NPHAMAX];  /* max. grain density for the second solid phase nucleation sites per [cm^3] */
   CA_FLOAT NucParams_poly[NPHAMAX][4];
/*END THUINET*/

   int gd_max_total;    /* total max. # grains for array size   */
   int oriented;        /* Orientation Calc: True/Flase         */
   CA_FLOAT NucParams[4];   /* the nucleation parameters array      */
   CA_FLOAT NucParamsBeut[4];   /* the nucleation parameters array for binary eutectic      */
   CA_FLOAT NucParamsTeut[4];   /* the nucleation parameters array for ternary eutectic     */
   CA_FLOAT NucParamsSurf[4];   /* The nucleatio parameters array for the surface */

                        /* Perturbation things                  */
   int perturb_on;        /* Using Fidap? TRUE/FALSE              */
   int n_perturb;  /* Using Fidap? TRUE/FALSE              */
   Nuc_Perturb_str perturb[MAX_NUC_PERTURB];  /* structure to add perturbations       */
   Nuc_dist_str nucdist; 
   CA_FLOAT (*rand_function) (CA_FLOAT * params); /* random number function with required distribution */
   CA_FLOAT (*rand_function_beut) (CA_FLOAT * params); /* random number function for binary eutectic */
   CA_FLOAT (*rand_function_teut) (CA_FLOAT * params); /* random number function ternary eutectic */
   CA_FLOAT *beut_threshold;
   CA_FLOAT *teut_threshold;
} Nuc_str;
/****************************************************************/
/* SB_buc_str contains nucleation values calculated             */
/* for a particular subblock (at a particular iteration)        */
/****************************************************************/
typedef struct sbnuc  {
                       /* Nucleation stuff for each subblock  	*/
   int Ni_active;      /* number made active                   	*/
   CA_FLOAT N_nuc_old;    /* Number calculated for previous iter. */
   CA_FLOAT N_nuc;        /* Number calculated this iteration     	*/
   CA_FLOAT Ni_sum;       /* sum of Ni that should be active      	*/
   CA_FLOAT fract_nuc;    /* carry fractions of nuclei             */
} SB_nuc_str;


#endif /* NUCPROPS_H */
/*
*/
