
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

/*RCS Id:$Id: prop_wrapper.c 1356 2008-08-18 13:41:15Z  $*/
/*
*/
/* include system headers */
#include "stdio.h"
#include "math.h"
#include "time.h"

/* include header files requred by subroutines */
#include "machine.h"
#include "read_ctrl.h"          /* included for def. of READ_CTRL subroutine */
#include "blocks.h"             /* included for def. of INIT_BB, CALC_BB, FINISH_BB */

/* functions used from umat_solid.c */
extern int umat_solid (int stat_flag, CA_FLOAT time, CA_FLOAT delt, Ctrl_str * cp);

/* functions used from read_ctrl.c */
extern int read_ctrl (char *filename, Ctrl_str * cp);

#ifdef DEBUG
FILE *efd;
#endif /* DEBUG */

void print_usage (char *prog_name)
{                               /* print the usage message on error */
  fprintf (stderr, "   %s: Use the CA method to simulate grain\n", prog_name);
  fprintf (stderr, "   nucleation and growth.\n");
  fprintf (stderr, "   %s acts as a wrapper to test the umat_solid\n", prog_name);
  fprintf (stderr, "   that can also be called from the FEM software CAP.\n");
  fprintf (stderr, "   We also plan to allow postprocessing FEM results.\n");
  fprintf (stderr, "\n   Usage: %s contro_file\n", prog_name);
  fprintf (stderr, "\tWhere the control_file if supplied, tells the\n");
  fprintf (stderr, "\tprgramme which files to read information from.\n");
  fprintf (stderr, "\tThe following command line options are allowed:\n");
  fprintf (stderr, "\t-?\t\t-> print this message\n");
  fprintf (stderr, "\t\n");
}

/****************************************************************/
/* Beginning of the MAIN program!				*/
/* This is just a wrapper to simulate the grain subroutine	*/
/* being called from within CAP. It asks a few questions	*/
/* that that CAP would know, and then goes on from there.	*/
/* It gives you the option to enter all the thermophysical	*/
/* data that would be read from a file if using CAP.		*/
/* This file must be in the current directory called umat_mp.dat	*/
/****************************************************************/
extern CA_FLOAT get_dl (CA_FLOAT temp);
extern CA_FLOAT get_ds (CA_FLOAT temp);
extern CA_FLOAT schiel (CA_FLOAT Tcell);
extern CA_FLOAT find_sch_conc (CA_FLOAT tempK, CA_FLOAT fl);
extern CA_FLOAT find_sat (CA_FLOAT cell_tempK, CA_FLOAT cell_si, CA_FLOAT cell_fs);
extern CA_FLOAT getav_d (CA_FLOAT dl, CA_FLOAT ds, CA_FLOAT fs);

void main (int argc, char *argv[])
{

  /* declare counters and output variables for main */
  int i;
  int stat_flag = -1;           /* stage in analysis: -1=before, 0=during, 1=final */

  /* variables for input */
  Ctrl_str ctrl, *cp;
  int cflg, errflg;
  int finput = FALSE;           /* inquire control info from user as default */
  char *ctrl_fname;

  int step, pr_step;
  CA_FLOAT ptime, del_ptime = 1.0;
  CA_FLOAT time = 0.0;
  CA_FLOAT delt = 1.0;
  char *time_string;            /* local date/time in string */
  FILE *propsfile;
  CA_FLOAT tempK, temp, tempdrop, tempstart, tempend;
  CA_FLOAT kH;
  CA_FLOAT d_av, ss, Hinit, Hliq, Hsol, dl, ds, fs, con, cellsat, cellsatFU;

  propsfile = fopen ("props.txt", "w");

  kH = 0.1;
  Hinit = .816;
  cp = &ctrl;
  tempstart = 650;
  tempend = 550;
  tempdrop = 1.0;
#ifdef NORMAL
  fprintf (propsfile, "temp\tdl\tds\tfs\td_av\tsi-con\tHliq-FU\tcellsat\tcellsatFU\tHss\n");

  for (temp = tempstart; temp > tempend; temp -= tempdrop) {
    tempK = temp + STD_TEMP;
    dl = get_dl (temp);
    ds = get_ds (temp);
    fs = schiel (temp);
    d_av = getav_d (dl, ds, fs);
    con = find_sch_conc ((tempK), (1 - fs));
    cellsat = find_sat (tempK, con, fs);
    cellsatFU = cellsat / MPMETERCUB;
    Hliq = Hinit / (1 - (1 - kH) * fs);
    ss = Hliq / cellsatFU;
    fprintf (propsfile, "%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\n", temp, dl, ds, fs, d_av, con, Hliq,
             cellsat, cellsatFU, ss);
  }
#endif /* NORMAL */
  fprintf (propsfile, "temp/con");
  for (con = 7; con < 12.5; con += 0.1) {
    fprintf (propsfile, "\t%.5g", con);
  }
  fprintf (propsfile, "\n");

  for (temp = tempstart; temp > tempend; temp -= tempdrop) {
    tempK = temp + STD_TEMP;
    fprintf (propsfile, "%.5g", temp);
    for (con = 7; con < 12.5; con += 0.1) {
      cellsat = find_sat (tempK, con, 0.5);
      fprintf (propsfile, "\t%.5g", cellsat);

    }
    fprintf (propsfile, "\n");
  }

  fclose (propsfile);
  fprintf (stderr, "Finshed making %s\n", "props.txt");
  exit (0);
}                               /* end of main program, umat_wrapper */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_prop_wrapper_c ()
{
  static char const rcsid[] = "$Id: prop_wrapper.c 1356 2008-08-18 13:41:15Z  $";

  return (rcsid);
}

/* end of rcs_id_prop_wrapper_c subroutine */
