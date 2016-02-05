
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
/* multi_diff_props.h                                                   */
/* Header file defining concentrationsolution related structures*/
/* for multi-component diffusion                                */
/****************************************************************/
/****************************************************************/
/* Written by A. Chirazi Imperial College */
/* October 4, 2000                                                 */
/*RCS Id:$Id: multi_diff_props.h 1339 2008-07-23 13:58:29Z  $*/
/****************************************************************/

#ifndef MULTIPROPS_H
#define MULTIPROPS_H

/****************************************************************/
/*  the solute MultiS_struct structure, which holds                */
/*  information about solute phases and pahse diagram              */
/****************************************************************/
typedef struct multisolprops {
                        /***********T AND FL VARIABLES***********/
   CA_FLOAT *Cinit_multi;  /* the initial concentration array [mol/m^3] */
   CA_FLOAT *Cmax_multi;   /* max concentration array of whole sb     */
   CA_FLOAT *Cmin_multi;   /* min concentration aaray of whole sb     */
   CA_FLOAT *Cavg_multi;   /* Avg. concentration array of whole sb     */
   CA_FLOAT *Ctot_multi;   /* tot. concentration array of whole sb     */
   CA_FLOAT *slope_multi;  /* the slope for each solute phase */
   CA_FLOAT **Clim_multi;   /* limiting concentration for each phase and cell*/

   CA_FLOAT *part_coef_multi; /* partition coefficeint for different elements */
   CA_FLOAT **part_coef_matrix; /*variation of Ki for each cell*/
   CA_FLOAT *LDiff_multi;   /* liquid diffusion array */
   CA_FLOAT *SDiff_multi;  /* solid diffusion array */
   CA_FLOAT **Diff_matrix_liq; /* diffusion parameters matrix, liquid phase */
   CA_FLOAT **Diff_matrix_sol; /* diffusion parameters matrix, solid phase */
   CA_FLOAT **temp_multi;   /* a range of temparory arrays for the*/
                        /* storage of the solute values */
                        /* old values are stored at each time step*/
   CA_FLOAT **begin_buffer; /* the begin pointer of the conc buffer */
   CA_FLOAT phaval[30000][4];     /*phase diagram values for liqidus and al si cu conc*/
   CA_FLOAT **ya1, **ya2, **ya3, **xxa; /**interpolation matrices for phase dia**/
   int counter[5];           /* counter array for phase dia ***/
   int numtietri;      /*number of tie triangles produced by thermocalc*/

   /***********phase diagram critical points (bin and ter eut points) **/
   CA_FLOAT bin_eut_max[3];
   CA_FLOAT bin_eut_temp[3];
   CA_FLOAT ter_eut_max[3];
   CA_FLOAT ter_eut_temp[3];
} MultiS_struct;
#endif /* MULTIPROPS_H */
