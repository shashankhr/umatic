
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
/* READ_MAT.C:   (Part of CA code)                              */
/* Subroutine to read the material properties from a file.      */
/* The file is formated, using the style:                       */
/*    # as first character:     Comment                         */
/* and values are input in the format:                          */
/*    command value  #comments                                  */
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Jul 1, 1998                                                  */
/* Ludovic Thuinet                                              */
/* 2005                                                         */
/****************************************************************/

/****************************************************************/
/* Versions maintained with CVS                                 */
/* see log at end of file ***************************************/
/****************************************************************/
/*RCS id $Id: calc_solprops_poly.c 1339 2008-07-23 13:58:29Z  $*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "machine.h"
#include "read_ctrl.h"
#include "readmat.h"
#include "blocks.h"
#include "solprops_reader.h"
#include "solprops_writer.h"
void calc_solprops_poly(BB_struct *bp, Mat_str *mp, Solute_props *sp, int tp_flag, int ieq_tot) {
  /* subroutine to calculate derived properties for POLYCOMPONENT alloys */
  /* this routine shoud not be called at all if poly-component mode is not selected */
  /*THUINET 02/05 et 04/05 */

  int ieq;

  if (tp_flag == 1) {

    fprintf (stderr, "T_pure specified for polycomponent alloys, OVERRIDING liquidus\n");

    for (ieq = 0; ieq < ieq_tot; ieq++) {
      mp->Tliq_poly[ieq] += sp->m_solute[ieq] * sp->Cinit;
      fprintf (stderr, "Liquidus Temperature Tliq[%d] = %.6g\n", ieq, mp->Tliq_poly[ieq]);
    }

  } else {

    fprintf (stderr, "T_pure not specified for polycomponent alloys, calculating..\n");

    for (ieq = 0; ieq < ieq_tot; ieq++) {
      mp->tp_poly[ieq] -= sp->m_solute[ieq] * sp->Cinit;
      fprintf (stderr, "Reference Temperature Tp[%d] = %.6g\n", ieq, mp->tp_poly[ieq]);
    }

  }

  sp->Fs_eut = 1.0;

  /*FIN THUINET 02/05 and 04/05 */

}/* end of calc_solprops_poly */
