
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
/* blocks.h                                                                           */
/* Header file defining superblock structure                          */
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Jul 1, 1998                                                  */
/* Jan 19 1999 RCA added fixed neighbourhood choices            */
/* June22 1999 RCA added some more diffusion stuff              */
/****************************************************************/
/****************************************************************/
/*RCS Id:$Id: blocks.h 1342 2008-07-23 15:45:00Z  $*/


#ifndef BLOCKS_H
#define BLOCKS_H


#include <time.h>

#define BLOCKREV "blocks.h $Revision: 1342 $"
#ifndef MACHINE_H
  #include "machine.h"
#endif /* MACHINE_H */
#ifndef PORE_H
  #include "pore.h"
#endif /* PORE_H */
#ifndef MATPROPS_H
  #include "matprops.h"
#endif /* MATPROPS_H */
#ifndef NUCPROPS_H
  #include "nucprops.h"
#endif /* NUCPROPS_H */
#ifndef READCTRL_H
  #include "read_ctrl.h"
#endif /* READCTRL_H */
#ifndef TEMPROPS_H
  #include "temprops.h"
#endif /* TEMPROPS_H */
#ifndef SOLROPS_H
  #include "solprops.h"
#endif /* SOLPROPS_H */
#ifndef MULTIPROPS_H
  #include "multi_diff_props.h"
#endif /* MULTIPROPS_H */
#ifndef FEM_H
  #include "fem.h"
#endif /* FEM_H */
#ifndef FIDAP_INTERP_H
  #include "fidap.h"
#endif /* FIDAP_INTERP_H */
#ifndef GRAIN_H
  #include "grain.h"
#endif /* GRAIN_H */
#ifndef CUBE_H
  #include "cube.h"
#endif /* CUBE_H */
#ifndef TC_INTERP_H
  #include "tcpl.h"
#endif /* TC_INTERP_H */
#ifndef NBHD_H
   #include "nbhd.h"
#endif /*NBHD_H*/
#ifndef NEARNODE_H
   #include "nearnode.h"
#endif /*NEARNODE_H*/

#ifndef PR_STRUCT_H
   #include "pr_struct.h"
#endif /*PR_STRUCT_H*/

#ifndef  SURCELL_H 
   #include "SurCell.h"
#endif /* SURCELL_H */

typedef struct blocksize_str{
    int fsize;
    int csize;
    int bbsize;
    int sbsize;
    int gsize;
    int psize;
    int PORE_size ;
    int pcnsize ;
    int pclsize ;
    int nucpropssize;
} bsize;
#ifndef TEMP_STRUCT_H
   #include "temp_struct.h"
#endif /*TEMP_STRUCT_H*/

typedef struct values {
   CA_FLOAT ** block_array;  /* the array of pointers into value array subblocks */
   FILE * fd_exel;
   char id_string[MAX_STRING_LEN]; /* id string pre-pended to output for this value*/

   CA_FLOAT disp_max; /*the max. value for displaying this value */
   CA_FLOAT part_coef; /* the partition coefficient for this value */
   CA_FLOAT CLmin;
   CA_FLOAT CLmax;
   CA_FLOAT Cmin;
   CA_FLOAT Cmax;
   CA_FLOAT Ctot;
   CA_FLOAT SSmax;
   CA_FLOAT SSmin;
   CA_FLOAT SATmax;
} Value_struct;


#include "subblock.h"
#include "bigblock.h"

#define EXTERNAL_INIT_BB (-3)
#define READ_BB  (-2) 
#define INIT_BB  (-1)
#define CALC_BB  (0)
#define FINISH_BB  (1)
#define WRITE_BB  (2)
#define RESTART_BB  (3)

#endif /* BLOCKS_H */
/*
*/
