
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
/* output.h                                                     */
/* Header file defining output structure                        */
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Aug 15, 1998                                                 */
/*RCS Id:$Id: read_ctrl.h 1382 2008-09-24 14:59:34Z  $*/
/****************************************************************/


#ifndef READCTRL_H
#define READCTRL_H
#define CTRLREV "read_ctrl.h $Revision: 1382 $"

#ifndef OUTPUT_H
  #include "output.h"
#endif /* OUTPUT_H */


/****************************************************************/
/* Ouput_Str holds things which set output options     */
/****************************************************************/
typedef struct ctrl_opt {
   char * cflags;
   int jflg;                    /* behaviour due to signal handler*/
   int restart;                 /* is this a restart? */
   int restart_pore_on;         /* reinitialize with pores turned on at restart */
   int restart_gas_on;          /* reinitialize the gas array at restart */
   int mould_nuc;
   int mould_src;              /*  mould function 0=no 1=yes */
   int mould_src_pert;          /* perturb the mould source? */

   CA_FLOAT mould_source_value;
   CA_FLOAT mould_source_pert;  /* magnitude of perturbation */
   CA_FLOAT mould_source_freq;  /* spatial lenght scale of pert. */

   int solo;                    /* flag: true running alone       */
   int cap;                     /* flag to state if CAP mode      */
   int post;                    /* flag: true - postprocessing    */
   int input;                   /* flag: true - input bigblock    */
   int extrudemould;            /* flag: 1=read just xy slice and extrude 0=read whole domain */
   int t_input;                 /* flag: true - temperature input */
   int fgrid_input;             /* flag: true - input fgrid       */
   int fg_tr;                   /* flag: true - transient       */
   int con_cast;                /* contiuous casting postpro      */
   int particle;                /* particle mode 0=no 1=yes       */
   int diffuse;                 /* gas diffusion 0=no 1=yes          */
   int pore;                    /* pores 0=no 1=yes               */
   int temp_lookup;             /* Ali version lookup for temperature */
   int pr_lookup;               /* lookup for pressure */
   int swap_axes;
   int diffuse_alloy;
   int gradtilt;
   double grad_angle[3];
   int diffuse_alloy_multi;     /* 0=no and 1=yes */

/*THUINET 07-02-05*/
   int diffuse_alloy_poly;      /* 0=no and 1=yes */
/*End THUINET 07-02-05*/

   int thermocalc;              /* 0=no and 1=yes */
   int temp_dep_diff;           /* 0=no and 1=yes */
   int interpolate;             /* 1, 2 or 3 */
   int show_eut;                /*show eutectic as special colour*/
   int phase_diag_on;           /* use phase diag. for cell Tund*/
   int das_limrad;              /* use limiting radius based on das*/
   int global_undercooling;     /* use und. wrt bulk liquidus for nucleation */
   int diffuse_step;            /*  no of diff steps per growth */
   int window_moving;           /* true window moving */ /*by Wei WANG on 11-07-02*/
   int melt_back;               /* melt-back */
   int isotherm_curv;           /* producing concave or convex isotherm */
   int init_cont;               /* 0 initial 1 continuous calculation */ /*by Wei WANG on 11-07-02 */
   int decentred_octahedron;    /* true use decentred square_octahedron algorithm */ /*by Wei WANG on 11-07-02*/
   int random_angles;           /* 0 none 1 original (2d) 2 uniform3d */
   int external;                 /* the external mode is on *******/
   int umat_feedback;             /* the fully coupled version is used*/
   int flow_on;                 /* the flow option of external is on */
   int curvature_3D;		/* PLIC algorithm for 3D curvature */
   int curvature_2D;             /* PLIC algorithm for 2D curvature */
				/* Control of CA neighbourhood  */
   int umat_method;               /* use push or pull method?     */
   int scheil;                   /* use schiel rules?            */
   int n_neigh;                 /* # of neighbours to use       */
   int32_t seed;              /* rand. num. gen. seedval      */

				/* FILES FOR INPUT              */
   char *fn_umat_extern;                /* Filename for umat_extern file */
   char *fn_ctrl;                      /* Filename the control options were read from*/ 
   char *fn_block_restart;             /* Filename for restart data */ 
   FILE *fd_umat_extern;                /* File Handle for umat_extern file */
   char *fn_phadia;             /* Filename for phase diagram   */
   FILE *fd_phadia;             /* File Handle for phase diagram */
   char *fn_cap;                /* Filename for CAP CA file     */
   FILE *fd_cap;                /* File Handle for CAP CA file  */
   char *fn_geo;                /* Filename for readgeoplus     */
   FILE *fd_geo;                /* File Handle for readgeoplus  */
   char *fn_mat;                /* Filename for readmat         */
   FILE *fd_mat;                /* File Handle for readmat      */
   char *fn_inp;                /* Filename for read Input      */
   FILE *fd_inp;                /* File Handle for read Input   */
   char *fn_base;               /* Filename base for all output */
   char *fn_fgrid;              /* filename for input of FGrid_str nodal values data */
   char *fn_t_input;
   char *fn_solprops_gas;       /* filename for input of gas solute props */
   char *fn_solprops_alloy[NSOLMAX];     /* filename for input of alloy solute props */
				/* VARIABLES RELATED TO OUTPUT  */
   int nsbslice;                /* # of SB slices to print      */
   int nbbslice;                /* # of BB slices to print      */
   int pore_dump_sb;            /* dump pore radius for sb #    */
   int write_block;
   int grainslice;              /* do line intercept on slice # */
   int slice[MAX_CTRL][2];      /* gif ouput of slice: sb#, sl# */
   int bbslice[MAX_CTRL];       /* gif ouput of BB slice:  sl#  */
   int nfsprint;                /* create output when fs is reached */
   CA_FLOAT fsprint[MAX_CTRL];
   int do_conc_prof;            /* do a conc profile> true/false */
   int conc_prof[3];            /* concentration profile sb#,sl#,row*/
   int time_dump;               /* interpret output step as time (not absolute steps */
   int time_exp;                /* unit of time 0=seconds -3 = millisecond +1 = 10sec */
   CA_FLOAT time_unit;          /* unit of time calculated from time_exp */
   int slice_dmp_freq;          /* Slice  print Freq.           */
   int blk_dmp_freq;          /* Slice  print Freq.           */
   int scr_dmp_freq;            /* Screen print Freq.           */
   int nscrdumps;               /*est. number of screen dumps */
   int tempslice;                /*print a temperature slice of bb*/
   int floatdump;               /* dump the float info for slice*/
   int excel;                   /* Printout excel files 1/0=t/f */
   FILE *fd_ex;                 /* File Handle for excel output */
   FILE *fd_ex_diff;            /* File Handle for conc. output */
   FILE *fd_ex_pore;            /* File handle for pores output */
   FILE *fd_ex_pore_stats;            /* File handle for pores output */
   FILE *fd_ex_pore_extent;     /* File handle for  pore extent output */
   RGB_struct *rgbp;             /* RGB colouring scheme         */
   int rgbmode;                 /* RGB colouring scheme          */
   int rgbgrey;
				/* VARIABLES FOR CONCENTRATION OUTPUT */
   int diff_rgbmode;            /* RGB colouring scheme for diff */
   int diff_log_disp;           /* Logarithmic Concentration Display*/
   int diff_ratio_disp;         /* Use max at each step for ratio*/
   CA_FLOAT alloy_disp_max;           /* Value for max colour          */
   CA_FLOAT gas_disp_max;           /* Value for max colour          */
   int diff_disp_cap;           /* use const colour over max, or wrap? */
				/* VARIABLES RELATED TO TESTING */
   int fixed_Pore;               /* true/false: use fixed porenucl.  */
   int nfPore;                   /* number of fixed nuclei       */
   CA_FLOAT nPsite[MAX_CTRL][5];      /* location of fixed nuclei     */

   int fixed_nuc;               /* true/false: use fixed nucl.  */
   int block_nuc;
   int nfnuc;                   /* number of fixed nuclei       */
 
/*THUINET 04/05*/
   int nfnuc_poly[NPHAMAX];
/*END THUINET 04/05*/

/*THUINET 04/05*/
  /*CA_FLOAT nsite[MAX_CTRL][7];*/      /* fixed nuclei */ /*modified by Wei WANG on 29-01-03*/
  CA_FLOAT nsite[MAX_CTRL][8];
/*END THUINET 05/05*/  

/*THUINET 04/05*/
   int nsite_poly[NPHAMAX][MAX_CTRL][7];
/*END THUINET 04/05*/

   int coolrate;               /* Const. cooling */
   CA_FLOAT fs_finish;
   CA_FLOAT ref_pres;             /* reference pressure used for pressure calculation */
   CA_FLOAT delT;                  /* cooling rate  */
   int NUM_COMP;                  /*number of component*/
   int ele_1;
/*THUINET 04/05*/
   int NUM_PHS; /*number of solid phases for multiphase solidification*/
   int NUM_EQ; /*number of equilibria for multiphase solidification*/
   int NUM_NUC_LAW;/*number of distinct gaussian nucleation law*/
   int ietat; /*define the state of a cell*/
/*END THUINET*/ 
   /* logically combined options -- precalculated to save time */

   int use_global_undercooling;
   int use_csol_alloy;
   int use_cell_nuc;
   int use_cell_temp;
   int eut_nuc_option;
   int strontium;
   /* for debugging memory with dbMalloc */
   unsigned long int hist[5];
   /*DBM - can't ifdef this without messing up the read_bin_blocks*/

} Ctrl_str;

/****************************************************************/
/* Define all of the default values below to be used            */
/* if the user does not specify them.                           */
/****************************************************************/
#define D_SOLO TRUE             /* flag: true running alone     */
#define D_PHASE_DIAG FALSE
#define D_GLOBAL_UNDERCOOLING FALSE
#define D_SHOW_EUT FALSE
#define D_DAS_LIMRAD FALSE
#define D_INP  FALSE            /* flag to state if CAP mode    */
#define D_CAP  FALSE            /* flag to state if CAP mode    */
#define D_POST FALSE            /* flag: true - postprocessing  */
#define D_EXTERNAL FALSE
#define D_RESTART FALSE
#define D_RESTART_PORE FALSE 
#define D_RESTART_GAS FALSE 
#define D_FEEDBACK FALSE
#define D_FLOW_ON FALSE
#define D_EutNucOption 0         /* default value for the eutectic nucleation threshold setting */
#define D_StrontiumOption 0      /* default value for the strontium flag */
#define D_REF_PRES 1.0           /* default for the reference pressure */
#define D_TEMPSLICE FALSE         /* flag: true -  print temp slice*/
#define D_PARTICLE FALSE
#define D_DIFFUSE FALSE         /* flag: true - diffuse gas     */
#define D_DIFFUSE_ALLOY FALSE   /* flag: true - diffuse alloy   */
#define D_DIFFUSE_ALLOY_MULTI FALSE  /* flag for multi diffusion */
#define D_TEMP_DEP_DIFF FALSE    /* flag for temperature dependent diffusion */
#define D_INTERPOLATION LINEAR1  /* flag for the choice of inter function */
#define D_DIFFUSE_STEP 1        /* no of diff steps per gro.st. */
#define D_WINDOW_MOVING FALSE   /* by Wei WANG on 11-07-02*/
#define D_INIT_CONT 0           /* by Wei WANG on 11-07-02*/
#define D_DECENTRED_OCTAHEDRON FALSE /*by Wei WANG on 11-07-02*/
#define D_ISOTHERM_CURV FALSE
#define D_MELT_BACK FALSE
#define D_CON_CAST FALSE        /* flag: true - con_cast mode   */
#define D_USE_GRAD_TILT FALSE
#define D_GRAD_ANGLE (0) 
#define D_CAP_FILE "bigblock.dat"
#define D_BLOCK_RESTART_FILE "no_such_file.blk"
#define D_GEO_FILE "umat_geoplus.in"
#define D_MAT_FILE "umat_matprop.in"
#define D_INP_FILE "umat_test.inp"
#define D_PHADIA_FILE "umat_phadia.in"
#define D_FG_FILE "umat_fg.fbn"
#define D_PORE TRUE
#define D_BASEFILE "umat_testout"
#define D_FS_FINISH 1.1
#define D_SEEDVAL 98724351	/* random num. gen. seedval     */
#define D_THERMOCALC FALSE
#define D_NUM_COMP 3           /* default num of components */
#define D_NUMTIETRI 20         /* default num of tie triangles produced by thermocalc for the ternary alloy */

/****************************************************************/
/* Define CA method defs...                                     */
/****************************************************************/
#define CA_PULL 1               /* turn node solid dep. on neigh*/
#define CA_PUSH 2               /* if node turn sol, push on n's*/
#define CA_NOSCH 0               /* if node turn sol, push on n's*/

/****************************************************************/
/* Define all of the default output values below to be used     */
/* if the user does not specify them.                           */
/****************************************************************/
#define D_CA_FLOATDUMP 0
#define D_EXCEL 1               /* Printout excel files 1/0=t/f */
#define D_SLICE_FREQ 10         /* Screen print Freq.           */
#define D_SCREEN_FREQ 1         /* Screen print Freq.           */
#define D_BLK_FREQ 0         /* block print Freq.           */
#define D_RGB_RANDOM 1          /* RGB colour scheme            */
#define D_DIFF_RGBMODE 1        /* RGB colour scheme for diff.  */
#define D_DIFF_DISP_MAX 100     /*max for concentration display */

/****************************************************************/
/* Define default DEBUGGING states                              */
/****************************************************************/
#define D_DEL_TEMP 1.0          /* fixed T drop for testing     */



#endif /* READCTRL_H */
/*
*/



