
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

/*$Id: histo_struct.h 1342 2008-07-23 15:45:00Z  $*/
#ifndef HISTO_STRUCT_H
#define HISTO_STRUCT_H
typedef struct H_str {
   int nbins;
   CA_FLOAT binsize;
   CA_FLOAT minbin;
   int ndata;
} Histo_struct;
#endif
/*
$Log$
Revision 11.1  2006/03/01 18:20:39  rcatwood
Merging polycomponent and gas with meltback

Revision 10.2  2005/12/01 14:38:01  rcatwood
Merged xly_05 changes into the main trunk
Primarily involving melt-back

Revision 10.1.2.2  2005/11/23 18:18:53  rcatwood
Result of merging mould_source and xly meltback+curvature 2d versions

Revision 10.1  2005/11/03 11:56:46  rcatwood
New version number -- using mould_src as base

Revision 8.1.12.2  2005/11/02 11:55:05  rcatwood
Fixing up the revision nubmer after loss of repository

Revision 9.2  2003/10/27 20:01:11  rcatwood
Harmonized header file cpp protection
Fixed filename bug for restart

Revision 9.1  2003/08/14 14:38:35  rcatwood
Working merge with decentered/porosity/external, also including
Ali Chirazi's multicomponent (not tested in this version)

Revision 8.1.6.1  2003/01/22 16:53:44  rcatwood
Almost working read_fg version

Revision 8.1  2002/10/17 17:01:01  rcatwood
New version number! for decentered/porosity merge! Alpha Version!

Revision 1.2  2001/02/19 19:28:47  rcatwood
fixed histo
for grains

and also make TcTrace mode override const. cooling rate

*/

