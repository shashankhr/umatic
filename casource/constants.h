
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
/* constants.h                                                  */
/* Header file defining numerical constants and similar things  */
/*                                                              */
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Tue Nov 23 16:53:55 gmt  1999                                */
/****************************************************************/

/*RCS Id:$Id: constants.h 1390 2008-09-25 15:43:01Z  $*/

#ifndef CONSTANTS_H
#define CONSTANTS_H
/* line data */
#define N_DATA 100
#define ADD_DATA 100

/* fields per line */
#define N_FIELDS 10
#define ADD_FIELDS 5

/* z and r data*/
#define INIT_ARRAY 10
#define ADD_ARRAY 10

/*Version multicomponent multiphase THUINET 11/05*/

/*Max number of solutes in the system*/
#define NSOLMAX 10
/*Max number of phases or equilibrium between phases in the system*/
#define NPHAMAX 5 

/*End THUINET*/

#define CONSTREV "constants.h $Revision: 1390 $"

#define EXTERNAL_MIN_STEPS 10
#define FIRST_CALL -1
#define STEP_CALL   0
#define LAST_CALL   1
#define RESTART_CALL (3)

#define JFLG_END 0    /* default behaviour if a signal received */
#define JFLG_WRITEXIT 1
#define JFLG_WRITEBIN 2

/* used for debugging around a known timestep */
#define MANY_DUMP_START 775000 
#define MANY_DUMP_STOP  775500 
#define MANY_FREQ 1
/**/

/* cell state definitions */
#define LIQUID 0.0
#define SOLID 1.0
#define NOT_CASTING -1
#define NUC_THRESH_NONE (-1.0) 

#define EMBRYO 0.000001
#define RAD_EMBRYO 0.000001
#define COURANT_LIMIT 0.16
#define SOLID 1.0
#define NOT_CASTING -1
#define DFS_WARNING 0.4
#define EUT_WARNING 0  
#define MAX_WARN_MSG 5
#define MAX_ERR_MSG 5
#define WARN_EXIT_LIMIT 1000
#define DFS_CAP 0.8
#define MAX_CAP_MSG 5
#define MAX_NUC_MSG 5
#define CAP_GROWTH 0.5
#define START_TEMP_OFFSET 2000
#define CAP_TEMP_OFFSET 20
#define LIQUIDUS 1
#define CONC_MULTI 2
#define CONC_MULTI_MONO 3
#define LINEAR1 1
#define LINEAR2 2
#define REGRESSION 3
#define FROZEN 400.0 /* failsafe temperature to stop at */
#define MAX_CTRL 100 /* maximum number of fixed nuclei or slices */
/***********************************************************/
/* Some numerical constants etc...                         */
/***********************************************************/
#define SQRT2 1.4142136
#define THREE_BY_4PI    0.2387324146378    /* 3/4pi */
#define FOURPI_BY_THREE    4.1887902       /* 4/3 * pi */
#define PI    3.14159265       /*  pi */
#define PI_BY_FOUR      0.7853981633974    /* pi/4 */

#define PIBY180 (0.01745329)
#define  DEGRAD  (57.2957795)
#define SIN1  (0.0174524)
#define SIN44 (0.69465837)
#define SIN44M1 (0.677205964)

#define DEG43 (0.75049158)
#define DEG1  (0.01745329)


#define LARGE 1.0e+21
#define P_AP_PA	101000.0                 /* applied pressure default value in Pa*/
#define P_AP_ATM 1.0                       /* applied pressure default value in Atm*/
#define GAS_CONST_ATM	82.056875	       /* R in $cm^3 x Atm \over (g.mol x K)$ */
#define GAS_CONST_SI	8.31441  	          /* R in $J\over K mol$*/
#define GAS_CONST_CM	8.31441e6  	          /* R in $J\over K mol$*/
#define STD_TEMP        273.16             /* Std temp in K */
#define CBRT2		1.2599210499         	 /* cbrt(2,0): 1.259921 */
#define POW_3_TO_1PT5	5.1961524227	    /* pow(3,3.0/2.0): 5.196152 */

#define MPMMCB 1.174812e-9                 /*conv ml/100g stp to mol/mm cub*/
#define MPLITRE 1.174812e-3                /*conv ml/100g stp to mol/l*/
#define MPMETERCUB 1.174812                /*conv ml/100g stp to mol/m cub*/
#define CONV_ATM_CM 9.8716683e-4           /* conv pa.m to atm.cm */


#define SB_NEW 0
#define SB_OPEN 1
#define SB_CLOSED 2
#define SURF_NUC_ALLOC_START (10000)
#define SURF_NUC_ALLOC_STEP (1000)

/***********************************************************/
/* properties of the alloy solute    (SCHEIL)              */
/*  also used in phase diagram mode for pores              */
/***********************************************************/
/* default to SILICON */
#if  !defined MAGNESIUM && !defined SILICON && !defined A356 && !defined INCONEL && !defined MULTICOMP
  #define SILICON
#endif

#ifdef MAGNESIUM
   #define DOUTRE_F 0.0170
   /* the following values from NLM - Hiromi Nagaumi */
      #define DAS_PRE 48.837e-6
      #define DAS_EXP -0.3085
   /*end NLM values*/

#elif defined INCONEL
   #define DOUTRE_F 0
   #define TRAD_END 500
   #define DAS_PRE 50.0e-6
   #define DAS_EXP -0.3333

   #ifndef ALLOY_EXPONENT
      #define ALLOY_EXPONENT 0.33
   #endif
   #ifndef DAS_COS_THETA
      #define DAS_COS_THETA 0.9848 /* about 10 degrees */
   #endif /*DAS_COS_THETA*/

#elif defined MULTICOMP
   #define MAX_B_CONC_1  12
   #define MAX_B_CONC_2  27.0
   #define T_EUT_1 570.0
   #define T_EUT_2 540.0
   #define TP_1    650.0
   #define TP_2    650.0
   #define DOUTRE_F -0.0119
   #define FS_EUT 0.65
   #define T_EUT 570.0 /*C*/
   #define T_LIQ 616.7 /*C*/
   #define T_SOL 500.0
   #define TP 650.0
   #define M -7.1212
   #define PD_SLOPE -7.1212
   #define KB 0.13
   #define CINITB 7.0 /*w%*/
   #define MAX_B_CONC 12.6 /*w%*/
   #define TRAD_END 545

   #ifndef ALLOY_EXPONENT
      #define ALLOY_EXPONENT 0.33
   #endif
   #ifndef DAS_COS_THETA
      #define DAS_COS_THETA 0.9848 /* about 10 degrees */
   #endif /*DAS_COS_THETA*/
   #define TRAD_END 449 
#elif defined SILICON
   #define DOUTRE_F -0.0119
   #define TRAD_END 545
   #define DAS_PRE 50.0e-6
   #define DAS_EXP -0.3333

   #ifndef ALLOY_EXPONENT
      #define ALLOY_EXPONENT 0.33
   #endif
   #ifndef DAS_COS_THETA
      #define DAS_COS_THETA 0.9848 /* about 10 degrees */
   #endif /*DAS_COS_THETA*/
#elif defined A356
   #define DOUTRE_F -0.0119
   #define TRAD_END 545
   #define DAS_PRE 50.0e-6
   #define DAS_EXP -0.3333

   #ifndef ALLOY_EXPONENT
      #define ALLOY_EXPONENT 0.33
   #endif
   #ifndef DAS_COS_THETA
      #define DAS_COS_THETA 0.9848 /* about 10 degrees */
   #endif /*DAS_COS_THETA*/
#elif defined VAR

#endif /*MAGNESIUM or SILICON*/

   #ifdef BUG_OCT_02
      #undef TRAD_END
      #define TRAD_END 610
      #define TRAD_START 615
   #endif
/***********************************************************/
/* end of SCHEIL parameters section  (SCHEIL)              */
/***********************************************************/
/*******************************************/
/* Nuc distribution parameters             */
/*******************************************/
   #define N_NUC_MODELS 7
      /* the names of the models */
      #define N_RAPPAZ 1	/*  rappaz model                        */
      #define N_HETERO 2	/*  heterogenous model                  */
      #define N_DIST 3   	/*  distribution heterogenous model                 */
      #define N_RAP_DEL 4  /*use delta Undercooling and Rappaz                 */ 
      #define N_OLDF_DEL 5  /*use delta Undercooling and Oldfiled              */ 
      #define N_BLOCK 6     /* use method with threshold per cell in the block */

   #define D_GD_MAX_TOTAL 100000	/* Max number of grains in total*/
   #define N_NUC_FUNCS 6
   #define N_TUND_BINS 10000
   #define TUND_MAX 50
   #define U_MAX 10
   #define N_U_BINS 100
   #define K_ONE 1e19

/*******************************************/
/* end of Nuc distribution parameters      */
/*******************************************/
#define FLAT_GROWTH 7e-5
/* test gas diff. coeff. for constant Dl and Ds*/
#define TEST_DL 3.3e-7
#define TEST_DS 6.1e-8
/* constants used in pore routines */
   #define PORE_T_EUT_OFFSET 1.0 /* print pores with tnuc > TEUT + offset */
   #define PORE_NONE 0
   #define PORE_SPHERE 1
   #define PORE_OFF 2 
   #define PORE_TUBE 3
   #define PORE_LATENT 4
   #define PORE_MULTI 5
   #define PORE_FROZEN 6
   #define PORE_NTRAD 100 /* number of temperature benchmarks for pore data */
   #define PORE_EXTRA 25  /* if PORE_REALLOC defined, how many array elements to increase */
   #define TAU 0.1 /* dwell time constant */
   #define PORE_MIN_TUBE_FS 1e-4
   #define PORE_MAX_FS 0.9
#ifndef MIN_LIMRAD             /*minimum pore limrad */
#define MIN_LIMRAD 1e-7
#endif

/* the minimum number of allowed pores */
#define PNUC_MIN (1)
/* the different pore nuc function choices */
#define PNUC_GAUSS (0)
#define PNUC_STEP (1)
#define PNUC_TRUESTEP (2)
#define PNUC_FUNCTION (3)
#define PNUC_INPUT (4)
#define PNUC_GAUSSDEV (5)
#define PNUC_SQUARE (6)

#define P_NUC_MIN_SS 1.2
#define P_NUC_SIG_MULT 3.0
/* used for preparing nuclaetion  histogram for pores*/
#define P_MINBIN 1.0
#define P_NBINS 100
#define P_BINSIZE 0.05

/* random generator function for grains */
#define G_NUC_STEP 0
#define G_NUC_SQUARE 1
#define G_NUC_GAUSSDEV 2
#define G_NUC_TWOSTEP 3

#define G_NUC_MIN_UND 1.0e-5
#define G_NUC_SIG_MULT 6.0
/* used for preparing nuclaetion  histogram for grains*/
#define G_MINBIN 0.0
#define G_NBINS  50
#define G_BINSIZE 0.2
/* used for preparing size histogram for grains*/
#define G_SIZE_MINBIN  0.0
#define G_SIZE_NBINS  50
#define G_SIZE_BINSIZE 0.5 

#endif /*CONSTANTS_H*/
/*
*/
