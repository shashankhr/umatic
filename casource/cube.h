
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
/* CUBE.H:                                                      */
/* Structures and flags for copysurf.c and initface.c.          */
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Aug 24, 1998                                                  */
/****************************************************************/
/*      MODIFIED by:                                            */
/*  PDL: Aug 24, 1998                                           */
/*  RCA: Sept 4, 1998                                                            */
/*                                                              */
/*                                                              */
/****************************************************************/
/*RCS Id:$Id: cube.h 1341 2008-07-23 15:23:30Z  $*/
#ifndef CUBE_H
#define CUBE_H


/******************************************************/
/*  flags for copy_surf routines                      */
/******************************************************/
#define WRAP    0
#define PAD     1
#define FIX_BDY 2
#define INS     3
#define FLUX_BDY 4

#define FORWARD 0
#define BACK 1

#define TWO_D_CA   2
#define THREE_D 3
/******************************************************/
/* indices into face code lookup table                */
/******************************************************/
#define START   0
#define SKIP    1
#define JUMP    2
#define POFF    3
#define WOFF    4
#define ISTART  5
#define ISKIP   6
#define IJUMP   7
#define IPOFF   8
#define NBOFF   9
/* axis defs  -- might not use these after all as they are not useful as loop indices  */
#define XAXIS 0
#define YAXIS 1
#define ZAXIS 2

#define CM_INS  0
#define CM_FACE 1
#define CM_EDGE 2
#define CM_CORN 3
#define XY_INS 2
#define XY_EDGE 3
#define XY_CORN 4

/***********************************************/
/* a useful little structure that has x,y,z axes  */
/* this is used by Frame structure below          */
/***********************************************/
typedef struct axes {
   int x,y,z;     
} Geom;

/***********************************************/
/*  This structure stores the information      */
/*  that describes the symmetry of the (cubic) */
/*  neighbourhoods and subblock shapes used.   */
/*  elist and clist store the permutations of  */
/*  faces that adjoin each edge and corner     */
/*  so that the external neighbours of a cell  */
/*  at the edge of a subblock can be easily    */
/*  looked up.  The other members store the    */
/*  offsets in the main array of the various   */
/*  neighbour cells.                           */
/* it is a little confusing but it works!!!    */
/* and saves doing lots of cnditional tests    */
/***********************************************/
typedef struct symm { /* This assumes hexahedral subblocks! It should be*/ 
   CA_FLOAT ivalue; /*to pass the init. value for type FIX_BDY copy_surf*/
   CA_FLOAT dtbydx; /* flux coefficient deltat / deltax */
   int code;
   int nouts;
   int curr;
   int bbins[3];
   int ins[3];        /* to pass the inside block size       */
   int outs[3];       /* and not have to figure it out again */
   int flist[6][10];   /* the face offsets for pad copying    */
   int ndim;          /* the number of dimensions */
   int edge[3][4];    /* possibe to use even if not rectangular, I think!*/
   int elist[12][6];  /* the edge permutations       */
   int clist[8][6];   /* the corner permutations     */
   int dlist[8][6];   /* the flat 2-d corner permutations  */
   int neigh[6];      /* offset of normal neighbours */
   int corner[8];     /* location of corner cells    */
   int face[6];       /* offset of opposite face     */
   int facectrl[6];    /* flags to control face wrap/pad */
   Geom nn;    /* for the bigblock, # sblocks to avoid duplicate calc's   */
} Frame;

#endif /* CUBE_H */
/*
*/
