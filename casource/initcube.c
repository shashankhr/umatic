
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
/*  initcube.c:                                             */
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Jul 1, 1998                                                  */
/****************************************************************/
/*      MODIFIED by:                                            */
/****************************************************************/
/****** To Do List **********************************************/
/*General:                                                      */
/****************************************************************/
/*RCS Id:$Id: initcube.c 1341 2008-07-23 15:23:30Z  $*/
#include <stdio.h>
#include "machine.h"
#include "blocks.h"

/******************************************************/
/*                                                    */
/*          Initialize the edge permutation list      */
/*                                                    */
/* Face 0 = -ve x axis   Face 1 = +ve x axis          */
/* Face 2 = -ve y axis   Face 3 = +ve y axis          */
/* Face 3 = -ve z axis   Face 4 = +ve z axis          */
/*                                                    */
/* The first two members of the list are the face     */
/* directions of the outside neighbours. The other    */
/* four are the inside neighbours.                    */
/* THe middle two are the axis parrallel to the       */
/* current edge.                                      */
/*                                                    */
/******************************************************/
void init_edge_list (Frame * cubeptr)
{
  int i, j;
  int edgelist[12][6] = {
    {2, 4, 0, 1, 3, 5},         /* 0  */
    {2, 5, 0, 1, 3, 4},         /* 1  */
    {3, 4, 0, 1, 2, 5},         /* 2  */
    {3, 5, 0, 1, 2, 4},         /* 3  */
    {0, 4, 2, 3, 1, 5},         /* 4  */
    {1, 4, 2, 3, 0, 5},         /* 5  */
    {0, 5, 2, 3, 1, 4},         /* 6  */
    {1, 5, 2, 3, 0, 4},         /* 7  */
    {0, 2, 4, 5, 1, 3},         /* 8  */
    {1, 2, 4, 5, 0, 3},         /* 9  */
    {0, 3, 4, 5, 1, 2},         /* 10 */
    {1, 3, 4, 5, 0, 2}          /* 11 */
  };

  /* copy the codes into the frame structure */
  for (i = 0; i < 12; ++i) {
    for (j = 0; j < 6; ++j) {
      cubeptr->elist[i][j] = edgelist[i][j];
    }
  }
}

/******************************************************/
/*                                                    */
/*    Initialise the corner permutation list          */
/*    Similar to the edge list except three outside   */
/*    and three inside neighbours.                    */
/*                                                    */
/*                                                    */
/******************************************************/
void init_corner_list (Frame * cubeptr)
{
  int i, j;
  int cornlist[8][6] = {
    {0, 2, 4, 1, 3, 5},         /* 0 */
    {1, 2, 4, 0, 3, 5},         /* 1 */
    {0, 3, 4, 1, 2, 5},         /* 2 */
    {1, 3, 4, 0, 2, 5},         /* 3 */
    {0, 2, 5, 1, 3, 4},         /* 4 */
    {1, 2, 5, 0, 3, 4},         /* 5 */
    {0, 3, 5, 1, 2, 4},         /* 6 */
    {1, 3, 5, 0, 2, 4}          /* 7 */
  };
  /* copy the codes in to the frame structure */
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 6; ++j) {
      cubeptr->clist[i][j] = cornlist[i][j];
    }
  }
}

/*                                                           */
/* the special cases for an XY plane multiblock (nsb[z] = 1) */
/*                                                           */
void init_ins_list (Frame * cubeptr)
{
  int i, j;

  /* always 2 outside, in the +z and -z direction */
  /* so use the same structure as 3d edge blocks  */
  int edgelist[12][6] = {
    {4, 5, 0, 1, 2, 3},         /* 0 always  */
    {2, 5, 0, 1, 3, 4},         /* 1 unused  */
    {3, 4, 0, 1, 2, 5},         /* 2 unused  */
    {3, 5, 0, 1, 2, 4},         /* 3 unused  */
    {0, 4, 2, 3, 1, 5},         /* 4 unused  */
    {1, 4, 2, 3, 0, 5},         /* 5 unused  */
    {0, 5, 2, 3, 1, 4},         /* 6 unused  */
    {1, 5, 2, 3, 0, 4},         /* 7 unused  */
    {0, 2, 4, 5, 1, 3},         /* 8 unused  */
    {1, 2, 4, 5, 0, 3},         /* 9 unused  */
    {0, 3, 4, 5, 1, 2},         /* 10 unused */
    {1, 3, 4, 5, 0, 2}          /* 11 unused */
  };

  /* copy the codes into the frame structure */
  for (i = 0; i < 12; ++i) {
    for (j = 0; j < 6; ++j) {
      cubeptr->elist[i][j] = edgelist[i][j];
    }
  }
}

/******************************************************/
/*                                                    */
/*    Initialise the XY edge permutation list          */
/*    use the structure used for 3d corner            */
/*                                                    */
/*                                    three outside   */
/*    and three inside neighbours.                    */
/*                                                    */
/*                                                    */
/******************************************************/
void init_ins_edge_list (Frame * cubeptr)
{
  int i, j;

  /*  edge has 3 outside, including the +z and -z direction */
  /* so use the same structure as 3d corner blocks  */
  int cornlist[8][6] = {
    {5, 4, 2, 0, 1, 3},         /* 0 front */
    {5, 4, 1, 0, 2, 3},         /* 1 unused */
    {5, 4, 3, 0, 1, 2},         /* 2 back */
    {1, 3, 4, 0, 2, 5},         /* 3 unused */
    {5, 4, 0, 1, 2, 3},         /* 4 right */
    {5, 4, 1, 0, 2, 3},         /* 5 left */
    {0, 3, 5, 1, 2, 4},         /* 6 unused */
    {1, 3, 5, 0, 2, 4}          /* 7 unused */
  };
  /* copy the codes in to the frame structure */
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 6; ++j) {
      cubeptr->clist[i][j] = cornlist[i][j];
    }
  }
}

/******************************************************/
/*                                                    */
/*    Initialise the XY corner permutation list          */
/*                                                    */
/*                                    four outside   */
/*    and two inside neighbours.                    */
/*                                                    */
/*                                                    */
/******************************************************/
void init_ins_corn_list (Frame * cubeptr)
{
  int i, j;

  /*  edge has 4 outside, including the +z and -z direction */
  /*  use the same structure as 3d corner blocks  */
  int cornlist[8][6] = {
    {5, 4, 0, 2, 1, 3},         /* 0 front right */
    {5, 4, 1, 2, 0, 3},         /* 1 front left */
    {5, 4, 0, 3, 1, 2},         /* 2 back right */
    {4, 5, 1, 3, 0, 2},         /* 3 back left */
    {5, 4, 1, 0, 2, 3},         /* 4 unused */
    {5, 4, 0, 1, 2, 3},         /* 5 unused */
    {0, 3, 5, 1, 2, 4},         /* 6 unused */
    {1, 3, 5, 0, 2, 4}          /* 7 unused */
  };
  /* copy the codes in to the frame structure */
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 6; ++j) {
      cubeptr->dlist[i][j] = cornlist[i][j];
    }
  }
}

/******************************************************/
/*                                                    */
/*    Dump out the cube geometry for debugging.       */
/*                                                    */
/******************************************************/

void dumpgeom (Frame * cubeptr)
{
/*dump the edge info*/
  int a, b, i;

  for (a = 0; a < 3; a++) {
    printf ("Edge %i : ", a);
    for (b = 0; b < 4; b++) {
      printf ("%i ", cubeptr->edge[a][b]);
    }
    printf ("\n");
  }

  for (i = 0; i < 8; ++i) {
    printf ("corner %i : %i\n", i, cubeptr->corner[i]);
  }

}

/******************************************************/
/*                                                    */
/*                                                    */
/*  Subroutine to initialise the offsets and sizes    */
/*  of the neighbourhoods within the subblock, using  */
/*  global values passed via big block structure.     */
/*                                                    */
/*                                                    */
/******************************************************/
int init_cube (BB_struct * bp)
{
  int i, j, k;
  int up;
  int cnum;
  Frame *cubeptr;

  cubeptr = &(bp->cubeptr);
  /* test the geometry to see if it is supported */
  if ((bp->nsb[0] == 1 || bp->nsb[1] == 1) && bp->ntsb != 1) {
    fprintf (stderr, "ERROR:init_cube: Multiblock geometry not supported, %i %i %i\n", bp->nsb[0], bp->nsb[1], bp->nsb[2]);
    fprintf (stderr, "Only -- Single subblock\n     -- X-Y plane\n     -- Fully 3-d\ngeometry allowed!\n");
    exit (0);
  }
  /* flat multiblock */
  if (bp->nsb[2] == 1 && bp->ntsb != 1) {
    fprintf (stderr, "Initializing flat multiblock .... \n");
    init_ins_list (cubeptr);
    init_ins_edge_list (cubeptr);
    init_ins_corn_list (cubeptr);

    /* bad number , zero or negative in nsb array */
  } else if (bp->ntsb < 1) {
    fprintf (stderr, "ERROR:init_cube:Invalid multiblock ntsb %i \n", bp->ntsb);
    exit (0);
    /*              */
    /* normal case  */
    /* single or 3-d */
    /*               */
  } else if (bp->ntsb >= 1) {
    fprintf (stderr, "Initializing multiblock .... \n");
    init_corner_list (cubeptr);
    init_edge_list (cubeptr);
    /* impossible case? */
  } else {
    fprintf (stderr, "ERROR:init_cube:Impossible multiblock ntsb %i \n", bp->ntsb);
    exit (0);
  }

/*Get the basic cell dimensions from big-block*/
  for (i = 0; i < 3; i++)
    cubeptr->bbins[i] = bp->nsb[i] - 1;
  up = bp->nsb[0] * bp->nsb[1]; /* offset to go up a layer  */

/* initialise the neighbour offsets */
  cubeptr->neigh[0] = -1;
  cubeptr->neigh[1] = 1;
  cubeptr->neigh[2] = -bp->nsb[0];
  cubeptr->neigh[3] = bp->nsb[0];
  cubeptr->neigh[4] = -up;
  cubeptr->neigh[5] = up;

/*initialise the face offset array */
  cubeptr->face[0] = cubeptr->bbins[0];
  cubeptr->face[1] = -cubeptr->bbins[0];
  cubeptr->face[2] = cubeptr->bbins[1] * bp->nsb[0];
  cubeptr->face[3] = -cubeptr->bbins[1] * bp->nsb[0];
  cubeptr->face[4] = cubeptr->bbins[2] * up;
  cubeptr->face[5] = -cubeptr->bbins[2] * up;

/* initialise the edge array  I donth think this is needed any more.
	cubeptr->edge[XAXIS][0] =0;
	cubeptr->edge[XAXIS][1] =cubeptr->bbins[2]*up;
	cubeptr->edge[XAXIS][2] =cubeptr->bbins[1]*bp->nsb[0];
	cubeptr->edge[XAXIS][3] =cubeptr->bbins[2]*up+cubeptr->bbins[1]*bp->nsb[0];

	cubeptr->edge[YAXIS][0] =0;
	cubeptr->edge[YAXIS][1] =cubeptr->bbins[0];
	cubeptr->edge[YAXIS][2] =cubeptr->bbins[2]*up;
	cubeptr->edge[YAXIS][3] =cubeptr->bbins[0]+cubeptr->bbins[2]*up;
	
	cubeptr->edge[ZAXIS][0] =0;
	cubeptr->edge[ZAXIS][1] =cubeptr->bbins[0];
	cubeptr->edge[ZAXIS][2] =cubeptr->bbins[1]*bp->nsb[0];
	cubeptr->edge[ZAXIS][3] =cubeptr->bbins[0]+cubeptr->bbins[1]*bp->nsb[0];
*/

/*Initialise the corner array*/
  for (i = 0; i < 2; ++i) {
    for (j = 0; j < 2; ++j) {
      for (k = 0; k < 2; ++k) {
        cnum = k * 4 + j * 2 + i;
        cubeptr->corner[cnum] =
          k * (bp->nsb[2]) * (bp->nsb[1]) * (cubeptr->bbins[2]) + j * (bp->nsb[0]) * (cubeptr->bbins[1]) + i * (cubeptr->bbins[0]);
        /* AAAAARRRRRGGGGGHHHHHH */
      }
    }
  }

/* end initialise arrays */
  /* delta t / delta x for flux calculation */
  cubeptr->dtbydx = bp->delt / bp->size_c[0];
  return 0;
}

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_initcube_c ()
{
  static char const rcsid[] = "$Id: initcube.c 1341 2008-07-23 15:23:30Z  $";

  return (rcsid);
}

/* end of rcs_id_initcube_c subroutine */
/*
*/
