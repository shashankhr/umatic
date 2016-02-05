
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
/* solprops.h                                                   */
/* Header file defining concentrationsolution related structures*/
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* june 28, 1999                                                 */
/*RCS Id:$Id: solprops.h 1341 2008-07-23 15:23:30Z  $*/
/****************************************************************/

#ifndef SOLPROPS_H
#define SOLPROPS_H
#define SOLPROPREV "solprops.h $Revision: 1341 $"


struct bigblock;

/****************************************************************/
/*  the solute S_struct structure, which holds                */
/*  information about the solute and pahse diagram              */
/*  for a whole subblock (extensive values eg average, total..) */
/****************************************************************/
typedef struct sb_solprops {
   CA_FLOAT Cinit;         /* the initial concentration [mol/m^3] */

   CA_FLOAT Cmax;          /* max concentration of whole sb     */
   CA_FLOAT CLmax;          /* max concentration of whole sb     */
   CA_FLOAT CLmin;          /* max concentration of whole sb     */

   CA_FLOAT Cavg;          /* Avg. concentration of whole sb     */
   CA_FLOAT Ctot;          /* tot. concentration of whole sb     */
   CA_FLOAT Cmin;          /* min concentration of whole sb     */
   CA_FLOAT addsol;        /* keep track of the amount of solute */
   CA_FLOAT t_addsol;      /* keep tracke of added solute */
} S_struct;

/* structure for the properties associated with one component */
/* ie. gas or alloy in the basic configuration */
typedef struct solprops {
   char my_name[32]; /* store the name of the solute */
   int my_type; /* 0=gas 1=alloy */
   int my_num;  /* index in the array for alloy */
   int NUM_PHS;
   CA_FLOAT Cinit;
   CA_FLOAT Dliq;
   CA_FLOAT Dsol[NPHAMAX];
    /* coefficients for determining arrhenius eqn for diffusion coeff. */
   CA_FLOAT QaLiq;
   CA_FLOAT DoLiq;
   CA_FLOAT QaSol[NPHAMAX];
   CA_FLOAT DoSol[NPHAMAX];

   SrcFn_T mould_src;           /* choose the type of mould function 0=none */
   CA_FLOAT mould_source_value;
   int mould_src_pert;          /* perturb the mould source? */

   CA_FLOAT part_coef[NPHAMAX]; /*Solute partition coefficient for each equilibrium*/
   CA_FLOAT cs_stoechio; 
   CA_FLOAT km;
   CA_FLOAT kminv;
   CA_FLOAT c_eut;
/*   CA_FLOAT T_liq;*/  /*?*/
/*   CA_FLOAT T_pure;*/ /*?*/
   CA_FLOAT T_eut;
   CA_FLOAT Fs_eut;
   CA_FLOAT m_solute[NPHAMAX]; /*solute slope for each equilibrium surface*/
   CA_FLOAT m_inv_solute;
   CA_FLOAT surf_tens_coef;

   CA_FLOAT coef_doutre; 

   CA_FLOAT (*mould_src_func)(struct bigblock *bp, struct solprops *sp, CA_FLOAT cell_temp, CA_FLOAT conc, int i, int j, int k);
}Solute_props;

#endif /* SOLPROPS_H */
/*
*/
