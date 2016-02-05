
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
/* READ_CA_EXTERNAL.C:   (Part of CA code)                       */
/* Subroutine to read the initial values for umat_extern.        */
/* The file is formated, using the style:                       */
/*    # as first character:     Comment                         */
/* and values are input in the format:                          */
/*    command value  #comments                                  */
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee, Robert C. Atwood & A. Chirazi       */
/*                                           , Imperial College */
/* Feb 25, 2002                                                 */
/****************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "blocks.h"

#define D_CELL_OPTION 2
#define D_NSB 1
#define D_CELL_SIZE 10.0E-06
#define D_NC_SB 50

int read_umat_extern (Ctrl_str * cp, BB_struct * bp)
{
  char line[MAX_STRING_LEN];
  char *token;
  int i, error;
  int rflag = 0;

   /*********************************************************/
  /* Open the file                                         */
   /*********************************************************/
  if ((cp->fd_umat_extern = fopen (cp->fn_umat_extern, "r")) == NULL) {
    fprintf (stderr, "Error: can't open input file [%s]\n", cp->fn_umat_extern);
    exit (0);
  }

  while (fgets (line, MAX_STRING_LEN, cp->fd_umat_extern) != NULL) {

    /* ignore comment and blank lines */
    if (line[0] == '%' || line[0] == '#' || (token = strtok (line, " ,;\t")) == NULL) {
      continue;

      /*********************************************************/
      /* All values in the bigblock structure                  */
      /*********************************************************/
      /* CellOption int */
    } else if (strcasecmp (token, "CellOption") == 0) {
      if ((token = strtok (NULL, " ,;\t\n")) != NULL)
        bp->cell_option = atoi (token);
      else {
        bp->cell_option = D_CELL_OPTION;
        rflag++;
        fprintf (stderr, "Error: CellOption, default value used: %d.\n", bp->cell_option);
      }
      /* BigBlockOrigin CA_FLOAT CA_FLOAT CA_FLOAT */
    } else if (strcasecmp (token, "BigBlockOrigin") == 0 && bp->cell_option == 2) {
      for (i = 0; i < 3; i++) {
        if ((token = strtok (NULL, " ,;\t\n")) != NULL)
          bp->orig_bb[i] = atof (token);
        else {
          bp->orig_bb[i] = 0.0;
          fprintf (stderr, "Error: BigBlockOrigin, default value used: %f.\n", 0.0);
        }
      }
      /* Default if command not recognised */
    } else {
      fprintf (stderr, "Warning: Unknown command: %s.\n", token);
      rflag++;
    }
  }                             /* end while */

  return rflag;

}                               /* end of read_umat_extern routine */

/***************************************************/
/* rcs id routine to include rcs id in the program */
/* generated by make_rcs_sub.sh script             */
/***************************************************/
char const *read_umat_extern_c ()
{
  static char const rcsid[] = "$Id: read_ca_procast.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}
