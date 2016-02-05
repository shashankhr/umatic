
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
/* grain.h                                                      */
/* Header file defining grain structure.                        */
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Wed Oct. 28 1998                                             */
/****************************************************************/
/*RCS Id:$Id: grain.h 1341 2008-07-23 15:23:30Z  $*/

#ifndef GRAIN_H
#define GRAIN_H

#define GRAINREV "grain.h $Revision: 1341 $"  
/* define cpp variables for neighbourhood definition */
#define N_NEIGH 27
#define NEIGH_6  6
#define NEIGH_8  8
#define NEIGH_10  10
#define NEIGH_26  26


/****************************************************************/
/* Structure to define a single grain.                          */
/****************************************************************/
typedef struct grain {
        int num;                /* grain number                 */
        int state;              /* current state                */
        int sbnum;              /* number of sb it centres in   */
        int cell;               /* cell at centre               */
        int ncells;             /* # cells it occupies          */
        int ngrow;              /* # of growing cells */
        int blocks;		/* how many blocks for this grain*/
        CA_FLOAT NucTh; /* the threshold for the nucleation */
        CA_FLOAT TNuc;             /* cell temp at which grain nuc. */
        CA_FLOAT TunderNuc;        /* cell undercooling at which nuc*/
        CA_FLOAT CellConcNuc;      /* cell conc. at which nuc.     */
        int max[3];             /* highest x,y,z value cell    */
        int nuccell[3];
        int min[3];             /* lowest xyz value cell      */

                                /* stuff for dir. growth..      */
/* NOTE: mdir should be moved out with a limited # of dirs! */
        CA_FLOAT mdir[N_NEIGH];         /* Factor for directionnal growth*/
        CA_FLOAT dir_fact;         /* Factor for directional growth*/
        CA_FLOAT gro_fact;         /* Factor for directional growth*/
        CA_FLOAT ph;         /* Factor for hex growth*/
        CA_FLOAT pq;         /* Factor for quad growth*/
        CA_FLOAT po;        /* factor for oct growht          */
        CA_FLOAT dir_angle;        /* Crystal orientation          */

        /*orientation used in decentred octahedron*/ /*by Wei WANG 11-07-02*/
        CA_FLOAT ang[3];
        CA_FLOAT sang[3];       /*sin(ang[])*/
        CA_FLOAT cang[3];       /*cos(ang[])*/
        CA_FLOAT g[3][3];       /*rotation matrix*/

} Ind_grain;

#endif /* GRAIN_H */
/*
*/
