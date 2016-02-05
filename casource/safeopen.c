
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

/*$Id: safeopen.c 1339 2008-07-23 13:58:29Z  $*/
/* safeopen file opener */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
int safeclose (FILE * fp)
{
  int ret = 0;

#ifdef VERBOSE_FILE
  int fno = 0;                  /* the file descriptor number */

  fno = fileno (fp);
  fprintf (stderr, "SAFECLOSE: closing,%i\n", fno);
#endif /*VERBOSE_FILE */

  ret = fclose (fp);
  return (ret);
}

FILE *safeopen (const char *fname, const char *type)
{
  FILE *fp;
  int myerror;

#ifdef VERBOSE_FILE
  int fno = 0;                  /* the file descriptor number */
#endif /*VERBOSE_FILE */
  errno=0;
  myerror=0;
  fp = fopen (fname, type);
  myerror=errno; /* because errno may be changed by fprintf */
  /* exit if the file did not open */
  if (fp == NULL) {
    fprintf (stderr, "ERROR: safeopen: could not open file named [%s]\n", fname);
    fprintf (stderr, "       safeopen: Check for the following problems:\n");
    fprintf (stderr, "       safeopen: (1) The file may not exist,\n");
    fprintf (stderr, "       safeopen: (2) The file may be locked,\n");
    fprintf (stderr, "       safeopen: (3) The file may not have the correct permissions.\n");
    fprintf (stderr, "       safeopen: (4) The file name may be invalid (spaces, linefeeds included, etc.\n");
    fprintf (stderr, "       safeopen: The attempt was made to open the file with mode: \"%s\".\n", type);
    fprintf (stderr, "\nERROR message: %s\n",strerror(myerror));
    fprintf (stderr, "Exiting...\n");
    exit (17);
  }
#ifdef VERBOSE_FILE
  fno = fileno (fp);
  fprintf (stderr, "SAFEOPEN: %s,%i\n", fname, fno);
#endif /*VERBOSE_FILE */
  return (fp);
}

/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *safeopen_c ()
{
  static char const rcsid[] = "$Id: safeopen.c 1339 2008-07-23 13:58:29Z  $";

  return (rcsid);
}
