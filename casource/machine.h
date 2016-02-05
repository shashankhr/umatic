
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
/* machine.h                                                    */
/* Header file defining the cpp operatives                      */
/* that are machine and compile flag specific.                  */
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Wed Jul  1 18:38:31 bst 1998                                 */
/****************************************************************/
/*RCS Id:$Id: machine.h 1401 2008-11-13 15:32:56Z  $*/
#ifndef MACHINE_H
#define MACHINE_H
/* if specified, include the headers for dbmalloc debugging */
   #include <stdio.h>
#ifdef DBM
   #include <sys/types.h>
   #include "/usr/local/include/dbmalloc.h"
#elif DMALLOC
   #include <stdlib.h>
   #include <dmalloc.h>
#else /* no debugging of malloc */
   #include <stdlib.h>
#endif /* DBM debug malloc */

/* macro to replace the random generator using the preprocessor */
/* used for debugging random value assignment */
#ifdef REPLACERAND
   extern double replacerand(); 
   extern void replaceSrand(long int seedval); 

#  ifdef drand48
#      undef drand48
#   endif
#   define drand48 replacerand

#   ifdef srand48
#      undef sdrand48
#   endif
#   define srand48 replaceSrand

#endif
#define MAX_REWRITES (100)

#ifdef SGI
#define THE_CLOCK CLK_TCK
#else
#define THE_CLOCK CLOCKS_PER_SEC
#endif

/* the debug malloc header sets this flag, otherwise these funcitons are deactivated */
#ifndef _DEBUG_MALLOC_INC
   #define malloc_enter(func)
   #define malloc_leave(func)
   #define malloc_chain_check()
   #define malloc_dump(fd)
   #define malloc_list(a,b,c)
   #define malloc_inuse(hist)	(*(hist) = 0, 0)
#endif


#include "constants.h"
#include "safeopen.h"

#ifdef DMALLOC /* dmalloc debugging library instead (redhat default) */
#include <dmalloc.h>
#endif


/****************************************************************/
/* set to either single or double precision!                    */
/****************************************************************/ 
/* single precision! *//*by Wei WANG 31-07-02*/
#ifdef CA_DOUBLE
#define CA_FLOAT double
#define SCAN_F "%lf"
#define SQRT(A) sqrt(A)
#define SIN(a) sin(a)
#define POW(A,B) pow( (A),(B) )
#define EXP(A) exp( (A) )
#define CEIL(A) ceil(A)
#define FLOOR(A) floor(A)
#define LOGP(A) log1p(A)
#define LOG(A) log(A)
#define ABS(A) fabs(A)
#define MINVAL (1e-14)

#else /* SINGLE */
#define CA_FLOAT float
#define SCAN_F "%g"
#ifndef SGI
#define SIN(A) ((float) sin((double) A))
#define SQRT(A) ((float) sqrt( (double) A ))
#define POW(A,B) ((float) pow( (double) (A), (double) (B) ))
#define EXP(A) ((float) exp( (double) (A) ))
#define CEIL(A) ((float) ceil((double) A))
#define FLOOR(A) ((float) floor((double) A))
#define LOG(A) ((float) log((double) A))
#define LOGP(A) ((float) log1p((double)  A ))
#define ABS(A)  ((float) fabs((double) A))
#else
#define SQRT(A) sqrtf( A )
#define POW(A,B) powf( (A),(B) )
#define EXP(A) expf( (A) )
#define CEIL(A) ceilf( A )
#define FLOOR(A) floorf( A )
#define LOGP(A) (flog1p( A ))
#define LOG(A) (flog( A ))
#define ABS(A) fabsf(A)
#endif
#define MINVAL (1e-6)
#endif /* DOUBLE/SINGLE */
/* flags for diffusion routine */
#define GAS 0
#define ALLOY 1
/*flag to indicate eutectic solid*/
#define EUT_GRAIN -1
/*flags for nucleation distribution (pores)*/

/*logical flags*/
#if     !defined(TRUE) || ((TRUE) != 1)
#define TRUE    (1)
#endif

#define LINELENGTH 1200
#if     !defined(FALSE) || ((FALSE) != 0)
#define FALSE   (0)
#endif
/* useful small functions */
#ifndef MAX
#define MAX(A,B) ( (A) > (B) ? (A) : (B))
#endif /* MAX */

#ifndef MIN
#define MIN(A,B) ( (A) < (B) ? (A) : (B))
#endif /* MIN */

#define INIT (-1)


/***********************************************************/
/* Default sizes of character strings, etc...              */
/***********************************************************/
#define MAX_STRING_LEN (256)

typedef enum SourceFunctionType{
        MouldSourceNONE,
        MouldSourceCONSTANT,
        MouldSourceDIFF,
        MouldSourcePLIN,
        MouldSourceFLUX
} SrcFn_T;

/* base file name to use for autogenerated restarts and checkpoints */
/* better leave it as "step" , or alter any queue and external checkpointing scripts */

#ifndef DEFAULT_CHECKPOINT_NAME
#  define DEFAULT_CHECKPOINT_NAME ("step") 
#endif
#ifdef BL_COMPRESS
#   define READARRAY read_comp_array
#   define WRITEARRAY write_comp_array
#   define BL_EXT "blz"
#else
#   define READARRAY fread
#   define WRITEARRAY fwrite
#   define BL_EXT "blk"
#endif

#endif /* MACHINE_H */
/*
*/

