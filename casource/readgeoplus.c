
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
/* READGEOPLUS.C:   (Part of CA code)                           */
/* Subroutine to read the geometry and initial conditions.      */
/* The file is formated, using the style:                       */
/*    # as first character:     Comment                         */
/* and values are input in the format:                          */
/*    command value  #comments                                  */
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Jul 1, 1998                                                  */
/* Modified:                                                    */
/* dec 7, 1998 RCA -- added thermocopule trace option           */
/****************************************************************/
/*RCS Id:$Id: readgeoplus.c 1342 2008-07-23 15:45:00Z  $*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "machine.h"
#include "readgeoplus.h"
#include "blocks.h"

void set_geo_defaults (BB_struct * bp)
{
   /*********************************************************/
  /* Set all the default values...                         */
   /*********************************************************/
  /* from blocks.h */
  int i;

  if (bp->ctrl->external == FALSE) {
    for (i = 0; i < 3; i++) {
      bp->nsb[i] = D_NSB;
      bp->nc[i] = D_NC_SB;
      bp->size_bb[i] = D_SIZE_BB;
    }
    bp->finish_time = D_FINISH_TIME;
  }
  bp->delt = D_DELT;
  bp->fsgrow = D_FS_GROW;
  bp->window_velo = D_WINDOW_VELO;      /*by Wei WANG 11-07-02 */
  bp->window_disp = D_WINDOW_DISP;      /*by Wei WANG 11-07-02 */
  bp->iso_coef1 = D_ISO_COEF1;  /*curved isotherm*dong */
  bp->iso_coef2 = D_ISO_COEF2;  /*curved isotherm dong */
  bp->velo_coef = D_VELO_COEF;  /* Varying V */
  bp->grad_coef = D_GRAD_COEF;  /*varying G */
  /* below are things for solo and testing modes...    */
  bp->Tinit = D_TINIT;
  bp->dim = THREE_D;
  /* details for FIDAP mode */
  bp->Cbdy_alloy = 0;
  bp->Cbdy_gas = 0;

#ifdef JUNK
  bp->fg->tring_on = FALSE;
  bp->fg->ntring = 0;

#endif
  bp->auto_fin = FALSE;
  /* below are things for bc's ...    */
  for (i = 0; i < 6; i++)
    bp->cubeptr.facectrl[i] = 1;        /* default is wrap */
  for (i = 0; i < 3; i++)
    bp->nzones[i] = D_NZONES;
}

int read_geoplus_vals (Ctrl_str * cp, BB_struct * bp)
{
  char line[MAX_STRING_LEN];
  char *token;
  int i, error;
  CA_FLOAT cellsize;
  double a, b, c, d;
  float e;
  int rflag = 0;
  extern CA_FLOAT global_pressure;

   /*********************************************************/
  /* Open the file                                         */
   /*********************************************************/
  if ((cp->fd_geo = fopen (cp->fn_geo, "r")) == NULL) {
    fprintf (stderr, "Error: can't open input file [%s]\n", cp->fn_geo);
    exit (0);
  }

  while (fgets (line, MAX_STRING_LEN, cp->fd_geo) != NULL) {

    /* ignore comment and blank lines */
    if (line[0] == '%' || line[0] == '#' || (token = strtok (line, " ,;\t")) == NULL) {
      continue;

      /*********************************************************/
      /* All values in the bigblock structure                  */
      /*********************************************************/
      /* NSubblocks int int int */
    } else if (strcasecmp (token, "NSubblocks") == 0) {
      for (i = 0; i < 3; i++) {
        if ((token = strtok (NULL, " ,;\t")) != NULL)
          bp->nsb[i] = atoi (token);
      }
      /* NCellsPerSB int int int */
    } else if (strcasecmp (token, "NCellsPerSB") == 0) {
      for (i = 0; i < 3; i++) {
        if ((token = strtok (NULL, " ,;\t")) != NULL)
          bp->nc[i] = atoi (token);
      }
#ifdef ORIGIN
      /* BigBlockOrigin CA_FLOAT CA_FLOAT CA_FLOAT */
    } else if (strcasecmp (token, "BigBlockOrigin") == 0) {
      origin_flag = 1;
      for (i = 0; i < 3; i++) {
        if ((token = strtok (NULL, " ,;\t")) != NULL)
          bp->orig_bb[i] = atof (token);
      }
#endif /*ORIGIN*/
#ifdef BB_SIZE
        /* BigBlockSize CA_FLOAT CA_FLOAT CA_FLOAT */
    } else if (strcasecmp (token, "BigBlockSize") == 0) {
      for (i = 0; i < 3; i++) {
        if ((token = strtok (NULL, " ,;\t")) != NULL)
          bp->size_bb[i] = atof (token);
      }
#else /*input cell size */
    } else if (strcasecmp (token, "CellSize") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        cellsize = atof (token);

      for (i = 0; i < 3; i++) {
        bp->size_bb[i] = cellsize * bp->nc[i] * bp->nsb[i];
      }
#endif /*BB or Cell size */
      /* AutoFinishTime */
    } else if (strcasecmp (token, "Nzones") == 0) {
      for (i = 0; i < 3; i++) {
        if ((token = strtok (NULL, " ,;\t")) != NULL)
          bp->nzones[i] = atof (token);
      }
    } else if (strcasecmp (token, "AutoFinishTime") == 0) {
      bp->auto_fin = TRUE;
      /* FinishTime CA_FLOAT */
    } else if (strcasecmp (token, "FinishTime") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->finish_time = atof (token);
      /* TimeStep CA_FLOAT */
    } else if (strcasecmp (token, "TimeStep") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->delt = atof (token);
      /* InitialTemperature CA_FLOAT */
      /**  \todo  Check initial temperature is below liquidus  -- general -- easy */
    } else if (strcasecmp (token, "InitialTemperature") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->Tinit = atof (token);
      /* InitialPressure CA_FLOAT */
    } else if (strcasecmp (token, "InitialPressure") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        global_pressure = atof (token);
      /* FSGrow CA_FLOAT */
    } else if (strcasecmp (token, "FSGrow") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->fsgrow = atof (token);
      /* NumDimensions */
    } else if (strcasecmp (token, "NumDimensions") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->dim = atoi (token);

      /* GradSlope */
    } else if (strcasecmp (token, "GradSlope") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->grad_slope = atof (token);

      /* Gradient */
    } else if (strcasecmp (token, "Gradient") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->gradient = atof (token);

      /*  Curved isothermal coef1   */
    } else if (strcasecmp (token, "Iso_Coef_One") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->iso_coef1 = atof (token);

      /* Curved isothermal coef2  */
    } else if (strcasecmp (token, "Iso_Coef_Two") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->iso_coef2 = atof (token);
      /* Velocity for directional (non-fidap) */

    } else if (strcasecmp (token, "Velocity") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->velocity = atof (token);

      /* variables for melt back xly 20041018 */
      /* time held for constant initial velocity */

    } else if (strcasecmp (token, "Time_Velo1") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->time_velo1 = atof (token);
      /* time stop decreasing velocity */
    } else if (strcasecmp (token, "Time_Velo2") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->time_velo2 = atof (token);
      /* time stop increasing velocity */
    } else if (strcasecmp (token, "Time_Velo4") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->time_velo4 = atof (token);
      /* time held of zero velocity */
    } else if (strcasecmp (token, "Time_Hold") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->time_hold = atof (token);
      /* velocity decrease per second, negative value */
    } else if (strcasecmp (token, "Velo_Coef1") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->velo_coef1 = atof (token);
      /* velocity increase per second, positive value */
    } else if (strcasecmp (token, "Velo_Coef2") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->velo_coef2 = atof (token);

      /* Velocity variation (direct Sol) */
    } else if (strcasecmp (token, "Velo_Coef") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->velo_coef = atof (token);
      /*Gradient Coefficient */
    } else if (strcasecmp (token, "Grad_Coef") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->grad_coef = atof (token);
      /*temperature curvature */
    } else if (strcasecmp (token, "Time_Curv") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->time_curv = atof (token);
      /*temperature transient time */
    } else if (strcasecmp (token, "Time_Tran") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->time_tran = atof (token);
      /*cell number for curved isotherm */
    } else if (strcasecmp (token, "Cell_No") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->cell_no = atof (token);
      /*trsnsient coef. for curved isotherm */
    } else if (strcasecmp (token, "Therm_Coef") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->therm_coef = atof (token);
      /* OctFactor */
    } else if (strcasecmp (token, "OctFactor") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->oct_factor = atof (token);

      /* Window_Velo CA_FLOAT *//* by Wei WANG 11-07-02 */
    } else if (strcasecmp (token, "Window_Velo") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->window_velo = atof (token);

      /* FaceCtrl bot(int) bot(int) bot(int) bot(int) bot(int) bot(int)  */
    } else if (strcasecmp (token, "FaceCtrl") == 0) {
      for (i = 0; i < 6; i++) {
        if ((token = strtok (NULL, " ,;\t")) != NULL)
          bp->cubeptr.facectrl[i] = atof (token);
      }
      /* boundary concentraiont for FIX_BDY condition  */
    } else if (strcasecmp (token, "Cbdy_alloy") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->Cbdy_alloy = atof (token);
      /* boundary concentraiont for FIX_BDY condition  */
    } else if (strcasecmp (token, "Cbdy_gas") == 0) {
      if ((token = strtok (NULL, " ,;\t")) != NULL)
        bp->Cbdy_gas = atof (token);

      /* Default if command not recognised */
    } else {
      fprintf (stderr, "Warning: Unknown command: %s.\n", token);
      rflag++;
    }
  }                             /* end while */

   /*********************************************************/
  /* Calculate subsiduary values...                        */
   /*********************************************************/
  if (bp->tc.on)
    bp->ctrl->coolrate = FALSE;
  bp->ntsb = bp->nsb[0] * bp->nsb[1] * bp->nsb[2];
  if (bp->ntsb <= 0) {
    fprintf (stderr, "ERROR:read_geoplus_vals: bad values for subbblock numbers %i %i %i \n", bp->nsb[0], bp->nsb[1], bp->nsb[2]);
    fprintf (stderr, "exiting... \n");
    exit (0);
  }
  bp->ncsb = bp->nc[0] * bp->nc[1] * bp->nc[2];
  if (bp->ncsb <= 0) {
    fprintf (stderr, "ERROR:read_geoplus_vals: bad values for subblock size %i %i %i \n", bp->nc[0], bp->nc[1], bp->nc[2]);
    fprintf (stderr, "exiting... \n");
    exit (0);
  }
#ifdef ORIGIN
  for (i = 0; i < 3; i++) {
    if (!origin_flag) {
      bp->orig_bb[i] = 0.0;     /* set origin to 0,0,0 as default */
    }
    bp->tnc[i] = bp->nc[i] * bp->nsb[i];
    bp->size_c[i] = bp->size_bb[i] / ((CA_FLOAT) bp->tnc[i]);
  }
#else /*not ORIGIN */
  for (i = 0; i < 3; i++) {
    bp->orig_bb[i] = 0.0;       /* set origin to 0,0,0 as it does not matter */
    bp->tnc[i] = bp->nc[i] * bp->nsb[i];
    bp->size_c[i] = bp->size_bb[i] / ((CA_FLOAT) bp->tnc[i]);
  }
#endif /*ORIGIN*/
    bp->total_cell_number = bp->tnc[0] * bp->tnc[1] * bp->tnc[2];
  a = (double) bp->size_c[0];
  b = (double) bp->size_c[1];
  c = (double) bp->size_c[2];
  d = a * b * c;
  e = (float) d;
  fprintf (stderr, "x,y,z,volume,single-precision-volume, %.10g,%.10g,%.10g,%.10g,%.10g\n", a, b, c, d, e);
  d = (double) ((double) (bp->size_c[0]) * (double) (bp->size_c[1]) * (double) (bp->size_c[2]));
/*   bp->vol_c = (double)((double)(bp->size_c[0]) *  (double)(bp->size_c[1]) * (double)(bp->size_c[2]));*/
  bp->vol_c = bp->size_c[0] * bp->size_c[1] * bp->size_c[2];
  bp->yinv = 1 / ((CA_FLOAT) (bp->nc[1]));

  fprintf (stderr, "total #sb: %d, # cells/sb: %d\n", bp->ntsb, bp->ncsb);
  fprintf (stderr, "#cells in x,y,z: %d, %d, %d\n", bp->nc[0], bp->nc[1], bp->nc[2]);
  fprintf (stderr, "bigblock size: %f, %f, %f,", bp->size_bb[0], bp->size_bb[1], bp->size_bb[2]);
  fprintf (stderr, "cell vol:, %.10g\n", bp->vol_c);

#ifdef DEBUG_READGEO
  /* from blocks.h */
  fprintf (stderr, "NSubblocks %d %d %d\n", bp->nsb[0], bp->nsb[1], bp->nsb[2]);
  fprintf (stderr, "NCellsPerSB %d %d %d\n", bp->nc[0], bp->nc[1], bp->nc[2]);
  fprintf (stderr, "BigBlockSize %f %f %f\n", bp->size_bb[0], bp->size_bb[1], bp->size_bb[2]);
  fprintf (stderr, "FinishTime %f\n", bp->finish_time);
  fprintf (stderr, "TimeStep %f\n", bp->delt);
  fprintf (stderr, "InitialTemperature %f\n", bp->Tinit);
#endif /* DEBUG_READGEO */

  fclose (cp->fd_geo);
  fprintf (stderr, "Exiting read_geoplus().\n");
  return rflag;                 /* return the number of errors whilst reading file */
}                               /* end of subroutine read_geoplus */

int read_geoplus (Ctrl_str * cp, BB_struct * bp)
{
  int retval = 0;

  set_geo_defaults (bp);
  retval = read_geoplus_vals (cp, bp);
  return (retval);
}

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/
char const *rcs_id_readgeoplus_c ()
{
  static char const rcsid[] = "$Id: readgeoplus.c 1342 2008-07-23 15:45:00Z  $";

  return (rcsid);
}

/* end of rcs_id_readgeoplus_c subroutine */
/*
*/
