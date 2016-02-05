
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

#ifndef GD_IO_H
#define GD_IO_H 1

#include <stdio.h>
 
typedef struct gdIOCtx {
	int	(*getC)(struct gdIOCtx*);
	int	(*getBuf)(struct gdIOCtx*, void*, int);

        void     (*putC)(struct gdIOCtx*, int);
	int	(*putBuf)(struct gdIOCtx*, const void*, int);

	int	(*seek)(struct gdIOCtx*, const int);
	long	(*tell)(struct gdIOCtx*);

	void    (*free)(struct gdIOCtx*);

} gdIOCtx;

typedef struct gdIOCtx	*gdIOCtxPtr;

void Putword(int w, gdIOCtx *ctx);
void Putchar(int c, gdIOCtx *ctx);

void gdPutC(const unsigned char c, gdIOCtx *ctx);
int gdPutBuf(const void *, int, gdIOCtx*);
void gdPutWord(int w, gdIOCtx *ctx);
void gdPutInt(int w, gdIOCtx *ctx);

int gdGetC(gdIOCtx *ctx);
int gdGetBuf(void *, int, gdIOCtx*);
int gdGetByte(int *result, gdIOCtx *ctx);
int gdGetWord(int *result, gdIOCtx *ctx);
int gdGetInt(int *result, gdIOCtx *ctx);

int gdSeek(gdIOCtx *ctx, const int);
long gdTell(gdIOCtx *ctx);

#endif
