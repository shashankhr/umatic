
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

/*$Id: handlers.c 1356 2008-08-18 13:41:15Z  $*/
#include <stdlib.h>
#include <stdio.h>
#include <sys/signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>
#include <setjmp.h>
#include "machine.h"
#include "constants.h"
#ifdef GNU_TRAP_FPE
#define _GNU_SOURCE
#   include <fenv.h>
extern int feenableexcept (int excepts);
#endif
extern int jflg;
extern int the_signo;
extern jmp_buf env;
extern int errno;
extern int signal_change_freq;

void float_error (int signo)
{                               /* signal handler for sigfpe */
  static int count = 0;

  printf ("\nSignal handler: float_error: ERROR: fpe detected (signal %i) !!!!\n",signo);
  if (count == 0)
    printf ("\nRaising USR2 to obtain output step.\n");
    raise (SIGUSR2);
  count++;
  if (count >= 10)
    printf ("\n%s: More than %i fpe errors detected. Exiting...\n",__func__);
    exit (6);
  signal (SIGFPE, float_error);
  return;
}

void catchit (int signo)
{                               /* signal handler -- jump to finish_bb */
  printf ("\nCA_WRAPPER: catchit: Signal %d received \n", signo);
  the_signo = signo;
  fflush (NULL);
  FILE * fp;
  struct stat  status;
  
  if (signo == SIGUSR1) {
    jflg = JFLG_WRITEXIT;
  } else if (signo == SIGUSR2) {
    jflg = JFLG_WRITEBIN;
    signal (SIGUSR2, catchit);

    /********************************************************/
    /* mechanism for changing output timestep while runnign */
    /* if the file is there, then it gets read when the signal */
    /* is caught */

    if ( stat("umat_newstep.in",&status) == 0 ){
       char * line, *sep,*token;
       line = (char * ) calloc(MAX_STRING_LEN,sizeof(char));
       sep = (char *) calloc (MAX_STRING_LEN, sizeof (char));
       sprintf (sep, " ,;\t\n\r");
       fprintf(stderr,"%s: Changing the output step as requested.\n",__func__);

       /* change the time step */
       fp = fopen("umat_newstep.in","r");

         while (fgets (line, MAX_STRING_LEN, fp) != NULL) {
       /* ignore comment and blank lines */
       if (line[0] == '%' || line[0] == '#' || (token = strtok (line, sep)) == NULL) {
         continue;
       } else if (strcasecmp (token, "NewStep") == 0) {
         if ((token = strtok (NULL, sep)) != NULL)
           signal_change_freq = atoi (token);
           if (signal_change_freq <= 0 ) {
              signal_change_freq = 0 ;
           }
         }else{
           signal_change_freq = 0 ;
         }
      }
    }else{
      /* no change */
        signal_change_freq = 0 ;
        fprintf(stderr,"%s:Returning with output flag set.\n",__func__);
    }
    return;
  }
  longjmp (env, 1);
  return;
}

#ifdef GNU_TRAP_FPE
void enable_fpe_traps ()
{
  /*
   * This installs the SIGFPE signal handler and enables traps for
   *
   * debugging purposes.
   */

  fenv_t env;
  int excepts;

  /* Enable exception trapping. */

  excepts = FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW;
  feholdexcept (&env);          /* start clean */
  if (feenableexcept (excepts) == -1)
    exit (1);                   /* set traps */
}
#endif

/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *handlers_c ()
{
  static char const rcsid[] = "$Id: handlers.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}
