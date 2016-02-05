
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
/* readgeoplus.h                                                */
/* Header file defining default values via cpp operatives       */
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Aug 16, 1998                                                 */
/****************************************************************/
/*RCS Id:$Id: readgeoplus.h 1341 2008-07-23 15:23:30Z  $*/
/****************************/
/*                          */
/*  ADDING OPTIONS          */
/*                          */
/**************************************************/
/*                                                */
/* 1. add member to the structure                 */
/* 2. add default value                           */
/* 3. add lines to set_geo_defaults               */
/* 4. add parsing lines to read_geo_values        */
/* 5. modify read_bin_blocks and write_bin_blocks */
/*    as necessary (if pointer members are added) */
/* 6. guard the option functions with a flag      */
/*    as necessary                                */
/*                                                */
/**************************************************/


#ifndef READGEOPLUS_H
#define READGEOPLUS_H

/****************************************************************/
/* Define all of the default geometry values below to be used   */
/* if the user does not specify them.                           */
/****************************************************************/
#define D_FS_GROW 1.0
#define D_NSB 1                 /* # sub-blocks in x,y,z        */
#define D_SIZE_CELL 1e-4                 /* # sub-blocks in x,y,z        */
#define D_NC_SB 50              /* total # cells in x,y,z       */
#define D_SIZE_BB 7.0           /* size of big-block            */
#define D_FINISH_TIME 200.0     /* time to stop simulation      */
#define D_DELT 0.1              /* ca time step                 */
#define D_TINIT 1190            /* initial temperature          */
#define D_PINIT 1.0           /* initial pressure atm          */
#define D_DEL_TEMP 1.0          /* fixed T drop for testing     */

#define D_WINDOW_VELO 0.0       /* window velocity */ /*by Wei WANG 11-07-02*/
#define D_NZONES 1              /* dnumber of zones in subblock */
#define D_WINDOW_DISP 0.0       /* window displacement */ /*by Wei WANG 11-07-02*/
#define D_ISO_COEF1 0.0 	/* curved isotherm dong*/
#define D_ISO_COEF2 0.0		/*curved isotherm dong*/
#define D_VELO_COEF 0.0          /*varying V*/
#define D_GRAD_COEF 0.0         /*varying G */
#define D_FIDAP_VEL 0.1         /* vel. of VAR filling          */
#define D_FIDAP_TMIN 500.0      /* min. T if no FIDAP result    */
#define D_FIDAP_TMAX 2000.0     /* max. T if no FIDAP result    */
#define D_FIDAP_STATE 1    /* steady - 0/transient -1 state of temp. */
#define D_TCTRACE_FILE "trace.dat"
#define D_FIDAP_FILE "trace.dat"

#endif /* READGEOPLUS_H */
/*
*/
