
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

/*RCS Id:$Id: init_sb_neigh.c 1341 2008-07-23 15:23:30Z  $*/
#include <stdio.h>
#include "machine.h"
#include "blocks.h"
#define LP_COMMAND  fprintf(stderr,"ISBN:ii before %i"\
                                        "\n nouts: %i\n",ii,nouts);\
                    bp->sb[ii]->nouts=nouts;\
                    bp->sb[ii]->code=code;\
                    if(ii>=bp->ntsb) {fprintf(stderr,"INSB:ERROR too many sb");\
                    exit(0);}\
                    ii++;\
                    fprintf(stderr,"ISBN:Ecount %i / %i\n",ecount++,bp->ntsb);

/**************************************************/
/*                                                */
/* init_sb_neigh -- set up the pointers to the    */
/* neighbours of each subblock.                   */
/*                                                */
/**************************************************/
/* bp = pointer to bigblock                       */
/* flag = scheme for wrapping, padding, fixed     */
/*        boundary etc.                           */
/**************************************************/

int init_sb_neigh (BB_struct * bp, int flag)
{
  int ii, i, j, k;
  int ecount = 0;
  int nouts, code;
  int bbins[3];                 /* local version of inside dim. of bigblock */
  FILE *out;
  int xy_only = FALSE;

  for (i = 0; i < 3; i++)
    bbins[i] = bp->cubeptr.bbins[i];

  if (bbins[2] == 0)
    xy_only = TRUE;
  else if (bbins[0] == 0 || bbins[1] == 0) {
    fprintf (stderr, "ERROR:init_sb_neigh: Sorry, only 3-d or X-Y 2-dimensional multiblocks supported.\n");
    fprintf (stderr, "*******************************************************************************");
    fprintf (stderr, "* Please rewrite the array-copying subsystem if you require other geometries. *\n");
    fprintf (stderr, "*******************************************************************************");
    exit (0);
  }
#ifdef VERBOSE
  out = fopen ("code.table", "w");
#endif
  ii = 0;
  k = 0;
  j = 0;
  i = 0;

/**************************************************/
/*  If one singel subblock, easy out.             */
/**************************************************/
  if (bp->ntsb == 1) {

    bp->sb[ii]->nouts = 6;
    bp->sb[ii]->code = 0;
    return (0);
  } else if (xy_only) {

/**************************************************/
/*                                                */
/*                                                */
/*                                                */
/*   flat 2d multiblocks                          */
/*   XY plane only !                              */
/*                                                */
/*                                                */
/**************************************************/
/**************************************************/
/*  Front  Left corner (000)                */
/**************************************************/
    nouts = XY_CORN;
    code = 0;
    LP_COMMAND
/**************************************************/
/*  Front Edge (i00)                        */
/**************************************************/
      nouts = XY_EDGE;
    code = 0;
    for (i = 1; i < bbins[0]; i++) {
    LP_COMMAND}

/**************************************************/
/*  Front  Right corner (n00)                */
/**************************************************/
    nouts = XY_CORN;
    code = 1;
    LP_COMMAND
/**************************************************/
/* The bottom surface                             */
/**************************************************/
      for (j = 1; j < bbins[1]; j++) {
   /**************************************************/
      /* a   left edge                            */
   /**************************************************/
      nouts = XY_EDGE;
      code = 4;
      LP_COMMAND
   /**************************************************/
        /* a row of face units                     */
   /**************************************************/
        nouts = XY_INS;
      code = 0;
      for (i = 1; i < bbins[0]; i++) {
      LP_COMMAND}
   /**************************************************/
      /* a right edge                                   */
   /**************************************************/
      nouts = XY_EDGE;
      code = 5;
    LP_COMMAND}
/**************************************************/
/*  Back  Left corner                       */
/**************************************************/
    nouts = XY_CORN;
    code = 2;
    LP_COMMAND
/**************************************************/
/* Back Edge                               */
/**************************************************/
      nouts = XY_EDGE;
    code = 2;
    for (i = 1; i < bbins[0]; i++) {
    LP_COMMAND}
/**************************************************/
/*  Back Right corner                       */
/**************************************************/
    nouts = XY_CORN;
    code = 3;
    LP_COMMAND return (0);
  } else {
/**************************************************/
/*                                                */
/*                                                */
/*                                                */
/*   three d multiblocks                          */
/*                                                */
/*                                                */
/*                                                */
/**************************************************/
/**************************************************/
/* Bottom Front  Left corner (000)                */
/**************************************************/
    nouts = CM_CORN;
    code = 0;
    LP_COMMAND
/**************************************************/
/* Bottom Front Edge (00i)                        */
/**************************************************/
      nouts = CM_EDGE;
    code = 0;
    for (i = 1; i < bbins[0]; i++) {
    LP_COMMAND}
/**************************************************/
/* Bottom Front Right corner (00n)                 */
/**************************************************/
    nouts = CM_CORN;
    code = 1;
    LP_COMMAND
/**************************************************/
/* The bottom surface                             */
/**************************************************/
      for (j = 1; j < bbins[1]; j++) {
   /**************************************************/
      /* a  bottom left edge                            */
   /**************************************************/
      nouts = CM_EDGE;
      code = 4;
      LP_COMMAND
   /**************************************************/
        /* a row of bottom face units                     */
   /**************************************************/
        nouts = CM_FACE;
      code = 4;
      for (i = 1; i < bbins[0]; i++) {
      LP_COMMAND}
   /**************************************************/
      /* a right edge                                   */
   /**************************************************/
      nouts = CM_EDGE;
      code = 5;
    LP_COMMAND}
/**************************************************/
/* Bottom Back  Left corner                       */
/**************************************************/
    nouts = CM_CORN;
    code = 2;
    LP_COMMAND
/**************************************************/
/* Bottom Back Edge                               */
/**************************************************/
      nouts = CM_EDGE;
    code = 2;
    for (i = 1; i < bbins[0]; i++) {
    LP_COMMAND}
/**************************************************/
/* Bottom Back Right corner                       */
/**************************************************/
    nouts = CM_CORN;
    code = 3;
    LP_COMMAND
/**************************************************/
/* The central layers                             */
/**************************************************/
      for (k = 1; k < bbins[2]; k++) {
   /**************************************************/
      /*  The front left edge                           */
   /**************************************************/
      nouts = CM_EDGE;
      code = 8;
      LP_COMMAND
   /**************************************************/
        /*  the row in the front surface                  */
   /**************************************************/
        nouts = CM_FACE;
      code = 2;
      for (i = 1; i < bbins[0]; i++) {
      LP_COMMAND}
   /**************************************************/
      /*  The front right edge                          */
   /**************************************************/
      nouts = CM_EDGE;
      code = 8;
      LP_COMMAND
   /**************************************************/
        /*  An Internal Row                               */
   /**************************************************/
        for (j = 1; j < bbins[1]; j++) {
      /**************************************************/
        /*  The Left face                                 */
      /**************************************************/
        nouts = CM_FACE;
        code = 0;
        LP_COMMAND
         /**************************************************/
          /* The Internal blocks                            */
         /**************************************************/
          nouts = CM_INS;
        code = 0;
        for (i = 1; i < bbins[0]; i++) {
        LP_COMMAND}

      /**************************************************/
        /*  The Right face                                */
      /**************************************************/
        nouts = CM_FACE;
        code = 1;
      LP_COMMAND}
   /**************************************************/
      /*  End Internal Row                              */
   /**************************************************/
   /**************************************************/
      /*  The back  left edge                           */
   /**************************************************/
      nouts = CM_EDGE;
      code = 10;
      LP_COMMAND
   /**************************************************/
        /*  the row in the back  surface                  */
   /**************************************************/
        nouts = CM_FACE;
      code = 3;
      for (i = 1; i < bbins[0]; i++) {
      LP_COMMAND}
   /**************************************************/
      /*  The back  right edge                          */
   /**************************************************/
      nouts = CM_EDGE;
      code = 11;
    LP_COMMAND}
/**************************************************/
/*  End Central Layer                             */
/**************************************************/
/**************************************************/
/* Top Front  Left corner (000)                */
/**************************************************/
    nouts = CM_CORN;
    code = 4;
    LP_COMMAND
/**************************************************/
/* Top Front Edge (00i)                        */
/**************************************************/
      nouts = CM_EDGE;
    code = 1;
    for (i = 1; i < bbins[0]; i++) {
    LP_COMMAND}
/**************************************************/
/* Top Front Right corner (00n)                 */
/**************************************************/
    nouts = CM_CORN;
    code = 5;
    LP_COMMAND
/**************************************************/
/* The Top surface                             */
/**************************************************/
      for (j = 1; j < bbins[1]; j++) {
   /**************************************************/
      /* a  Top left edge                            */
   /**************************************************/
      nouts = CM_EDGE;
      code = 6;
      LP_COMMAND
   /**************************************************/
        /* a row of Top face units                     */
   /**************************************************/
        nouts = CM_FACE;
      code = 5;
      for (i = 1; i < bbins[0]; i++) {
      LP_COMMAND}
   /**************************************************/
      /* a right edge                                   */
   /**************************************************/
      nouts = CM_EDGE;
      code = 7;
    LP_COMMAND}
/**************************************************/
/* Top Back  Left corner                       */
/**************************************************/
    nouts = CM_CORN;
    code = 6;
    LP_COMMAND
/**************************************************/
/* Top Back Edge                               */
/**************************************************/
      nouts = CM_EDGE;
    code = 3;
    for (i = 1; i < bbins[0]; i++) {
    LP_COMMAND}
/**************************************************/
/* Top Back Right corner                       */
/**************************************************/
    nouts = CM_CORN;
    code = 7;
    LP_COMMAND
/**************************************************/
/* All done                                       */
/**************************************************/
  }
#ifdef VERBOSE
  fclose (out);
#endif
  return (0);
}

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_init_sb_neigh_c ()
{
  static char const rcsid[] = "$Id: init_sb_neigh.c 1341 2008-07-23 15:23:30Z  $";

  return (rcsid);
}

/* end of rcs_id_init_sb_neigh_c subroutine */
/*
*/
