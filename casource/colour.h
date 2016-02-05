
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
/* colour.h                                                     */
/* Header file defining colour lookup table for gif output      */
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Aug 16, 1998                                                 */
/****************************************************************/

/*RCS Id:$Id: colour.h 1341 2008-07-23 15:23:30Z  $*/

#ifndef COLOUR_H
#define COLOUR_H

/* Define a set of 32 colours for gif output */
#define COLOUR_SIZE 32

char colours[32*3] = {
0xff, 0xff, 0xff,	/*black*/
0xff, 0x00, 0x00,	/*red*/
0x00, 0xff, 0x00,	/*green*/
0x00, 0x00, 0xff,	/*blue*/

0xff, 0xa0, 0xa0,	/*red*/
0xa0, 0xff, 0xa0,	/*green*/
0xa0, 0xa0, 0xff,	/*blue*/
0xff, 0xa0, 0xa0,	/*red*/
0xf0, 0xff, 0xa0,	/*green*/
0xa0, 0xf0, 0xff,	/*blue*/

0xff, 0xa0, 0xf0,	/*red*/
0xa0, 0xff, 0xf0,	/*green*/
0xf0, 0xa0, 0xff,	/*blue*/
0xff, 0xf0, 0x40,	/*red*/
0x40, 0xff, 0xf0,	/*green*/
0xf0, 0x40, 0xff,	/*blue*/

0xaf, 0x10, 0xa0,	/*red*/
0xa0, 0xaf, 0x10,	/*green*/
0x10, 0xa0, 0xaf,	/*blue*/
0x2f, 0xa0, 0xa0,	/*red*/
0xa0, 0x2f, 0xa0,	/*green*/
0xa0, 0xa0, 0x2f,	/*blue*/

0xf0, 0x10, 0xf0,	/*red*/
0xf0, 0x30, 0xf0,	/*green*/
0xf0, 0x70, 0xff,	/*blue*/
0x2f, 0xf0, 0xa0,	/*red*/
0xa0, 0xff, 0xa0,	/*green*/
0xa0, 0xf0, 0x2f,	/*blue*/

0xbf, 0x40, 0xa0,	/*red*/
0xbb, 0x4f, 0x10,	/*green*/
0xb0, 0x40, 0xaf,	/*blue*/
0xbf, 0x40, 0xb0,	/*red*/
} ;

/* colour for eutectic grain */
#define EUT_RED 0    /*black*/
#define EUT_GREEN 0 /*black*/
#define EUT_BLUE 0 /*black*/

#endif /* COLOUR_H */
/*
*/
