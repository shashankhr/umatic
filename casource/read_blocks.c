
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


/*RCS Id: $Id: read_blocks.c 1402 2008-11-20 15:36:41Z  $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef BL_COMPRESS
#include <zlib.h>
#endif
#include <errno.h>
#include "machine.h"
#include "constants.h"
#include "blocks.h"
#include "read_ctrl.h"
#include "grain.h"
#include "safeopen.h"
#include "pore.h"
#include "pore_routines.h"
#include "SurCellRoutines.h"
#include "sb_tableofcontents.h"


extern int init_output_img (RGB_struct *);
extern int alloc_sb (BB_struct * bp, int sbnum);
extern int alloc_bb (BB_struct * bp);

extern CA_FLOAT global_pressure;
Sb_tableofcontents * sb_contents_list;
extern void write_sbtoc(FILE *fp, Sb_tableofcontents * stocp);

void read_bin_errmsg (const char *msg, int exflag)
{
  fprintf (stderr, "**************************************************\n");
  fprintf (stderr, "read_blocks:ERROR: %s\n", msg);
  fprintf (stderr, "That is too hard to fix. \n");
  fprintf (stderr, "quitting...\n\n\n");
  fprintf (stderr, "**************************************************\n");
  if (exflag)
    exit (1);
  return;
}

void rst (FILE * fp, char **s)
{
  char line[MAX_STRING_LEN];

  if (fgets (line, MAX_STRING_LEN, fp) == NULL) {
    fprintf (stderr, "read_blocks:ERROR: Cannot read a line\n");
    exit (1);
  }
  /*strip the \n return character off the end of the string */
  line[strlen (line) - 1] = (char) 0;
  if (strcmp (line,"NULL") == 0 ){
      /* this spot is unused */
      *s = NULL ;
  }else{
     *s = (char *) calloc (MAX_STRING_LEN, sizeof (char));
     memcpy (*s, line, sizeof (char) * strlen (line));
  }
}

void read_bin_ctrl (Ctrl_str ** cp, FILE * fp)
{
  int i;
  int nread;

  (*cp) = (Ctrl_str *) calloc (1, sizeof (Ctrl_str));
  nread = fread (*cp, sizeof (Ctrl_str), 1, fp);
  fprintf (stderr, "Read %i bytes for ctrl struct\n", nread * sizeof (Ctrl_str));

  /* get the file name strings */
  rst (fp, &((*cp)->fn_cap));
  rst (fp, &((*cp)->fn_geo));
  rst (fp, &((*cp)->fn_mat));
  rst (fp, &((*cp)->fn_inp));
  rst (fp, &((*cp)->fn_base));

  rst (fp, &((*cp)->fn_solprops_gas));
  for (i=0;i<NSOLMAX;i++){
     rst (fp, &((*cp)->fn_solprops_alloy[i]));
  }

  /** \todo  include these quantities too? test cflags? -- general   */
  (*cp)->cflags = (char *) malloc (sizeof (char));
  *((*cp)->cflags) = '\000';
  (*cp)->rgbp = (RGB_struct *) malloc (sizeof (RGB_struct));
  init_output_img ((*cp)->rgbp);
}

void read_f_bin_sb (FILE * fp, BB_struct * bp, int sbnum)
{
  fread (bp->sb[sbnum], sizeof (SB_struct), 1, fp);
  alloc_sb (bp, sbnum);
}

void read_bin_pore (FILE * fp, PORE_str * cp, long *trad_pos)
{
  int i;
  p_c_node *node;
  long this_pos;

  switch (cp->State) {
    /* cases where no nodes or arrays are allocated */
  case NOT_CASTING:
    break;

    /* cases where there are data arrays and pore nodes */
  case PORE_NONE:
  case PORE_SPHERE:
  case PORE_TUBE:
  case PORE_MULTI:
  case PORE_OFF:
  case PORE_FROZEN:
    cp->boundary = (p_c_list *) calloc (1, sizeof (p_c_list));
    cp->body = (p_c_list *) calloc (1, sizeof (p_c_list));

    node = (p_c_node *) calloc (1, sizeof (p_c_node));
    fread (node, sizeof (p_c_node), 1, fp);
    while (node->next != NULL) {
      add_to_list (cp->boundary, node);
      node = (p_c_node *) calloc (1, sizeof (p_c_node));
      fread (node, sizeof (p_c_node), 1, fp);
    }
    add_to_list (cp->boundary, node);

    /* save the location of the data (t,rad) information */
    this_pos = ftell (fp);
#ifndef LINUX
    if ((*trad_pos) != (this_pos)) {
      fprintf (stderr, "ERROR: trad_pos does not match! %i %i\n", this_pos, (*trad_pos));
    }
#endif

    cp->t_lists = (CA_FLOAT **) calloc (N_T_LISTS, sizeof (CA_FLOAT *));

    for (i = 0; i < N_T_LISTS; i++) {
      cp->t_lists[i] = (CA_FLOAT *) calloc (PORE_NTRAD, sizeof (CA_FLOAT));
      fread (cp->t_lists[i], sizeof (CA_FLOAT), cp->trad_last, fp);
    }

    break;

    /* must be an error */
  default:
    fprintf (stderr, "ERROR:write_bin_pore: This cannot be happening! Pore state %i\n", cp->State);
    exit (0);
    break;
  }                             /* end of pore state switch */

  return;
}

/* read the pore structure list */
void read_f_bin_sb_pores (FILE * fp, BB_struct * bp, int sbnum)
{
  int pnum;
  long *pore_data_pos, *pore_trad_pos, this_pos;
  SB_struct *sp = bp->sb[sbnum];

  pore_data_pos = (long *) calloc (sp->Npores, sizeof (long));
  pore_trad_pos = (long *) calloc (sp->Npores, sizeof (long));

  /* make space for the pore array */
  sp->porelist = (PORE_str *) calloc (sp->Npores, sizeof (PORE_str));

  /* read the array of pore structures */
  fread (sp->porelist, sizeof (PORE_str), sp->Npores, fp);

  /* read the data position arrays */
  fread (pore_data_pos, sizeof (long), sp->Npores, fp);
  fread (pore_trad_pos, sizeof (long), sp->Npores, fp);
  for (pnum = 0; pnum < sp->Npores; pnum++) {
    this_pos = ftell (fp);
#ifndef LINUX
    if (pore_data_pos[pnum] != this_pos) {
      fprintf (stderr, "ERROR:read_f_bin_sb_pores: data position does not match for pore %i -- this:%i file:%i\n", pnum, this_pos,
               (pore_data_pos[pnum]));
    }
#endif
    read_bin_pore (fp, sp->porelist + pnum, pore_trad_pos + pnum);
  }

  free (pore_data_pos);
  free (pore_trad_pos);
}

#ifdef BL_COMPRESS
/* use zlib to uncompress the data */
/* needs to be ifdef'd out to compile if zlib is not present */
int read_comp_array (void *datap, size_t vsize, int nmemb, FILE * fp)
{
  Bytef *comp_data;
  uLongf comp_size, datasize;
#ifdef GZFILES
  static int gnum=0;
  char fname[32];
  FILE * gzfile;

  snprintf(fname,32,"f%i.gz",gnum);
  gzfile=fopen(fname,"w");
#endif 

  datasize = (uLongf) (nmemb * vsize);

  /* read in the compressed buffer size and allocate enough space */
  fread (&comp_size, sizeof (uLongf), 1, fp);
  comp_data = (Bytef *) calloc (comp_size, sizeof (Bytef));

  /* read in the compressed data */
  fread (comp_data, sizeof (Bytef), comp_size, fp);
#ifdef GZFILES
  fwrite(comp_data,sizeof(Bytef),comp_size,gzfile);
  fclose(gzfile);
#endif

  /* call the routine from libz, this requires the -lz flag to the linker */
  uncompress ((Bytef *) datap, &datasize, comp_data, comp_size);

  free (comp_data);
  return ((int) datasize * sizeof (Bytef));
}
#endif /* BL_COMPRESS */

void read_f_bin_igrain_data (FILE * fp, BB_struct * bp, int * buf, off_t pos){
  /* grain # */
  fseeko(fp,pos,SEEK_SET);
  READARRAY (buf, sizeof (int), bp->ncsb, fp);
}

void read_f_bin_sb_data (FILE * fp, BB_struct * bp, int sbnum)
{
  long this_pos = 0;
  static int i = 0;

/*THUINET 03/05*/
  int iphs, isol, ele_num, ele_1;

/*end THUINET 03/05*/

  /* allocate the arrays in the sb, according to the ctrl information */
  /* which should determine the options. */
  /**  \todo  arrange to alter options when needed, eg running part-way without */
  /* porosity then adding it it or somethign like that ... */
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif

#ifdef USE_ELM
  READARRAY (bp->sb[sbnum]->c_elm, sizeof (int), bp->ncsb, fp);
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
#endif
  /* fraction solid */
  READARRAY (bp->sb[sbnum]->c_fs, sizeof (CA_FLOAT), bp->ncsb, fp);
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
  /* schiel */
  if (bp->ctrl->scheil == TRUE) {
    READARRAY (bp->sb[sbnum]->sch_fs, sizeof (CA_FLOAT), bp->ncsb, fp);
#ifdef TB
    this_pos = ftell (fp);
    fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
  }

  /* cell index for Wei WANG */
  READARRAY (bp->sb[sbnum]->index, sizeof (int), bp->ncsb, fp);
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif

  /* grain # */
  READARRAY (bp->sb[sbnum]->gr, sizeof (int), bp->ncsb, fp);
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif

  /* cell temperature for Wei WANG */
  READARRAY (bp->sb[sbnum]->c_temp, sizeof (CA_FLOAT), bp->ncsb, fp);
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif

  /* curvature for Wei WANG */
      /* protected by options rca 2007 01 15 */
      /**todo: finetune the options to avoid excessive memory allocation when the options are not used*/
  if ((bp->ctrl->curvature_3D != 0) || (bp->ctrl->curvature_2D != 0)) {
    READARRAY (bp->sb[sbnum]->curv, sizeof (CA_FLOAT), bp->ncsb+2, fp);
    READARRAY (bp->sb[sbnum]->norm_x, sizeof (CA_FLOAT), bp->ncsb+2, fp);
    READARRAY (bp->sb[sbnum]->norm_y, sizeof (CA_FLOAT), bp->ncsb+2, fp);
    READARRAY (bp->sb[sbnum]->norm_z, sizeof (CA_FLOAT), bp->ncsb+2, fp);
#ifdef TB
    this_pos = ftell (fp);
    fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
  }

  /* gas solute array */
  if (bp->ctrl->diffuse == TRUE) {
    READARRAY (bp->sb[sbnum]->c_sol, sizeof (CA_FLOAT), bp->ncsb, fp);
#ifdef TB
    this_pos = ftell (fp);
    fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
  }
  /*ALLOY solute array */
  if ((bp->ctrl->diffuse_alloy == TRUE) || (bp->ctrl->particle == TRUE)) {
    READARRAY (bp->sb[sbnum]->c_sol_alloy, sizeof (CA_FLOAT), bp->ncsb, fp);
#ifdef TB
    this_pos = ftell (fp);
    fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
    if (bp->ctrl->decentred_octahedron == TRUE) {
      READARRAY (bp->sb[sbnum]->c_eqv_alloy, sizeof (CA_FLOAT), bp->ncsb, fp);
#ifdef TB
      this_pos = ftell (fp);
      fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
    }
  }
  /* Wei WANG - dencentered octa arrays */
  if (bp->ctrl->decentred_octahedron == TRUE) {
    READARRAY (bp->sb[sbnum]->dc_d, sizeof (CA_FLOAT), bp->ncsb, fp);
#ifdef TB
    this_pos = ftell (fp);
    fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
    READARRAY (bp->sb[sbnum]->dc_x, sizeof (CA_FLOAT), bp->ncsb, fp);
#ifdef TB
    this_pos = ftell (fp);
    fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
    READARRAY (bp->sb[sbnum]->dc_y, sizeof (CA_FLOAT), bp->ncsb, fp);
#ifdef TB
    this_pos = ftell (fp);
    fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
    READARRAY (bp->sb[sbnum]->dc_z, sizeof (CA_FLOAT), bp->ncsb, fp);
#ifdef TB
    this_pos = ftell (fp);
    fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
  }

  /* multicomponent for L. THUINET 03/05 */
  ele_num = bp->ctrl->NUM_COMP; /* number of elements in the alloy */
  ele_1 = ele_num - 1;

  if ((bp->ctrl->decentred_octahedron == TRUE) && (bp->ctrl->diffuse_alloy_poly == TRUE)) {
    for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes by Ludovic THUINET */
      READARRAY (bp->sb[sbnum]->c_eqv_poly[isol], sizeof (CA_FLOAT), bp->ncsb, fp);
#ifdef TB
      this_pos = ftell (fp);
      fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
    }

    for (iphs = 0; iphs < bp->ctrl->NUM_PHS; iphs++) {  /*loop on solutes by Ludovic THUINET */
      READARRAY (bp->sb[sbnum]->c_fs_poly[iphs], sizeof (CA_FLOAT), bp->ncsb, fp);
#ifdef TB
      this_pos = ftell (fp);
      fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
    }

  }
  /*end multicomponent for THUINET */

  if ((bp->ctrl->block_nuc == TRUE)) {  /*nuc threshold array */
    READARRAY (bp->sb[sbnum]->c_nuc_thresh, sizeof (CA_FLOAT), bp->ncsb, fp);
#ifdef TB
    this_pos = ftell (fp);
    fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
  }

  /* the surface cell location array /structure */

  READARRAY (bp->sb[sbnum]->surface.c_surfp, sizeof (int), bp->sb[sbnum]->surface.n_alloc, fp);
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
  /* the surface cell xyz location array */

  READARRAY (bp->sb[sbnum]->surface.surf_xyz, sizeof (int[3]), bp->sb[sbnum]->surface.n_alloc, fp);
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "read_f_bin_sb_data: call %i file position %i\n", i++, this_pos);
#endif
}

/*read the grain array */
void read_bin_grain_data (BB_struct * bp, FILE * fp)
{
  int i;

  /*
     bp->gr = (Ind_grain **)calloc(bp->nprops.gd_max_total , sizeof(Ind_grain *));
   */

  for (i = 1; i <= bp->nprops.ngr; i++) {
    bp->gr[i] = (Ind_grain *) malloc (sizeof (Ind_grain));
    fread (bp->gr[i], sizeof (Ind_grain), 1, fp);
  }
  bp->gr[0] = NULL;
}

void test_bin_block (FILE * fp, bsize * bsz)
{
  const char revision[] = BLOCKREV;
  const char bbrevision[] = BIGBLOCKREV;
  const char sbrevision[] = SUBBLOCKREV;
  const char crevision[] = CTRLREV;
  const char grevision[] = GRAINREV;
  const char prevision[] = POREREV;
  const char nrevision[] = NUCPROPSREV;
  char line[MAX_STRING_LEN];
  int b;

  if (fgets (line, MAX_STRING_LEN, fp) == NULL) {
    fprintf (stderr, "read_blocks:ERROR: Cannot read a line\n");
    exit (1);
  }
  line[strlen (line) - 1] = (char) 0;
  fprintf (stderr, " read_blocks: reading block data from\n%s (file)\n%s (this program)\n", line, revision);

  /*test the revision numbers */
  /* test the cvs revision number of blocks.h */
  if (strcasecmp (line, revision) != 0) {
    fprintf (stderr, "**************************************************\n");
    fprintf (stderr, "* read_blocks: WARNING: block revision numbers do not match.\n");
    fprintf (stderr, "* In the file: %s, In this program: %s\n", line, revision);
    fprintf (stderr, "* Attempting to continue...\n");
    fprintf (stderr, "* If there are more errors below, the block revisions\n");
    fprintf (stderr, "* are probably not compatible.\n");
    fprintf (stderr, "* In this case you should checkout the correct\n");
    fprintf (stderr, "* older version of the ca program and use that,\n");
    fprintf (stderr, "* or use (or create) a filter version to convert the data.\n");
    fprintf (stderr, "**************************************************\n");
  }

  if (fgets (line, MAX_STRING_LEN, fp) == NULL) {
    fprintf (stderr, "read_blocks:ERROR: Cannot read a line\n");
    exit (1);
  }
  line[strlen (line) - 1] = (char) 0;

  /*test the revision numbers */
  /* test the cvs revision number of bigblock.h */
  if (strcasecmp (line, bbrevision) != 0) {
    fprintf (stderr, "**************************************************\n");
    fprintf (stderr, "* read_blocks: WARNING: bigblock revision numbers do not match.\n");
    fprintf (stderr, "* In the file: %s, In this program: %s\n", line, revision);
    fprintf (stderr, "* Attempting to continue...\n");
    fprintf (stderr, "* If there are more errors below, the block revisions\n");
    fprintf (stderr, "* are probably not compatible.\n");
    fprintf (stderr, "* In this case you should checkout the correct\n");
    fprintf (stderr, "* older version of the ca program and use that,\n");
    fprintf (stderr, "* or use (or create) a filter version to convert the data.\n");
    fprintf (stderr, "**************************************************\n");
  }
  if (fgets (line, MAX_STRING_LEN, fp) == NULL) {
    fprintf (stderr, "read_blocks:ERROR: Cannot read a line\n");
    exit (1);
  }
  line[strlen (line) - 1] = (char) 0;

  /*test the revision numbers */
  /* test the cvs revision number of subblock.h */
  if (strcasecmp (line, sbrevision) != 0) {
    fprintf (stderr, "**************************************************\n");
    fprintf (stderr, "* read_blocks: WARNING: block revision numbers do not match.\n");
    fprintf (stderr, "* In the file: %s, In this program: %s\n", line, revision);
    fprintf (stderr, "* Attempting to continue...\n");
    fprintf (stderr, "* If there are more errors below, the block revisions\n");
    fprintf (stderr, "* are probably not compatible.\n");
    fprintf (stderr, "* In this case you should checkout the correct\n");
    fprintf (stderr, "* older version of the ca program and use that,\n");
    fprintf (stderr, "* or use (or create) a filter version to convert the data.\n");
    fprintf (stderr, "**************************************************\n");
  }
  /* test the cvs revision number of read_ctrl.h */
  if (fgets (line, MAX_STRING_LEN, fp) == NULL) {
    fprintf (stderr, "read_blocks:ERROR: Cannot read a line\n");
    exit (1);
  }
  line[strlen (line) - 1] = (char) 0;
  fprintf (stderr, "reading ctrl revision: %s %s\n", line, crevision);
  if (strcasecmp (line, crevision) != 0) {
    fprintf (stderr, "**************************************************\n");
    fprintf (stderr, "* read_blocks: WARNING: Ctrl revision numbers do not match.\n");
    fprintf (stderr, "* In the file: %s, In this program: %s\n", line, crevision);
    fprintf (stderr, "* Attempting to continue...\n");
    fprintf (stderr, "**************************************************\n");
  }
  /* test the cvs revision number of grain.h */
  if (fgets (line, MAX_STRING_LEN, fp) == NULL) {
    fprintf (stderr, "read_blocks:ERROR: Cannot read a line\n");
    exit (1);
  }
  line[strlen (line) - 1] = (char) 0;
  fprintf (stderr, "reading grain revision: %s %s\n", line, grevision);
  if (strcasecmp (line, grevision) != 0) {
    fprintf (stderr, "**************************************************\n");
    fprintf (stderr, "* read_blocks: WARNING: grain revision numbers do not match.\n");
    fprintf (stderr, "* In the file: %s, In this program: %s\n", line, grevision);
    fprintf (stderr, "* Attempting to continue...\n");
    fprintf (stderr, "**************************************************\n");
  }

  /* test the cvs revision number of pore.h */
  if (fgets (line, MAX_STRING_LEN, fp) == NULL) {
    fprintf (stderr, "read_blocks:ERROR: Cannot read a line\n");
    exit (1);
  }
  line[strlen (line) - 1] = (char) 0;
  fprintf (stderr, "reading grain revision: %s %s\n", line, prevision);
  if (strcasecmp (line, prevision) != 0) {
    fprintf (stderr, "**************************************************\n");
    fprintf (stderr, "* read_blocks: WARNING: pore revision numbers do not match.\n");
    fprintf (stderr, "* In the file: %s, In this program: %s\n", line, prevision);
    fprintf (stderr, "* Attempting to continue...\n");
    fprintf (stderr, "**************************************************\n");
  }
  /* test the cvs revision number of nucprops.h */
  if (fgets (line, MAX_STRING_LEN, fp) == NULL) {
    fprintf (stderr, "read_blocks:ERROR: Cannot read a line\n");
    exit (1);
  }
  line[strlen (line) - 1] = (char) 0;
  fprintf (stderr, "reading nucprops revision: %s %s\n", line, nrevision);
  if (strcasecmp (line, nrevision) != 0) {
    fprintf (stderr, "**************************************************\n");
    fprintf (stderr, "* read_blocks: WARNING: nucprops revision numbers do not match.\n");
    fprintf (stderr, "* In the file: %s, In this program: %s\n", line, nrevision);
    fprintf (stderr, "* Attempting to continue...\n");
    fprintf (stderr, "**************************************************\n");
  }

  /* read in the size data -- first check the size of the size structure! */
  fread (&b, sizeof (int), 1, fp);
  if (b != sizeof (bsize)) {
#ifdef TB
    fprintf (stderr, "%i %i\n", b, sizeof (bsize));
#endif
    read_bin_errmsg ("bsize struct does not match.", 1);
  }
  fread (bsz, sizeof (bsize), 1, fp);

  /* test the size of CA_FLOAT - compiled precision */
  fprintf (stderr, "CA_FLOAT size: %i %i\n", bsz->fsize, sizeof (CA_FLOAT));
  if (bsz->fsize != sizeof (CA_FLOAT)) {
    read_bin_errmsg ("CA_FLOAT precision does not match.", 1);
  }

  /* test the size of big block struct */
  fprintf (stderr, "bigblock size: %i %i\n", bsz->bbsize, sizeof (BB_struct));
  if (bsz->bbsize != sizeof (BB_struct)) {
    read_bin_errmsg ("bigblock size does not match.", 1);
  }

  /* test the size of sub block struct */
  fprintf (stderr, "subblock size: %i %i\n", bsz->sbsize, sizeof (SB_struct));
  if (bsz->sbsize != sizeof (SB_struct)) {
    read_bin_errmsg ("subblock size does not match.", 1);
  }

  /* test the size of control struct */
  fprintf (stderr, "control block size: %i %i\n", bsz->csize, sizeof (Ctrl_str));
  if (bsz->csize != sizeof (Ctrl_str)) {
    read_bin_errmsg ("control block size does not match.", 1);
  }

  /* test the size of grain struct */
  fprintf (stderr, "grain struct size: %i %i\n", bsz->gsize, sizeof (Ind_grain));
  if (bsz->gsize != sizeof (Ind_grain)) {
    read_bin_errmsg ("grain struct size does not match.", 1);
  }
  /* test the size of poreprops struct */
  fprintf (stderr, "poreprops struct size: %i %i\n", bsz->psize, sizeof (P_str));
  if (bsz->psize != sizeof (P_str)) {
    read_bin_errmsg ("poreprops struct size does not match.", 1);
  }
  /* test the size of pore struct */
  fprintf (stderr, "pore struct size: %i %i\n", bsz->PORE_size, sizeof (PORE_str));
  if (bsz->PORE_size != sizeof (PORE_str)) {
    read_bin_errmsg ("pore struct size does not match.", 1);
  }
  /* test the size of pore list struct */
  fprintf (stderr, "pore list struct size: %i %i\n", bsz->pclsize, sizeof (p_c_list));
  if (bsz->pclsize != sizeof (p_c_list)) {
    read_bin_errmsg ("pore list struct size does not match.", 1);
  }
  /* test the size of pore node struct */
  fprintf (stderr, "pore node struct size: %i %i\n", bsz->pcnsize, sizeof (p_c_node));
  if (bsz->pcnsize != sizeof (p_c_node)) {
    read_bin_errmsg ("pore node struct size does not match.", 1);
  }
  /* test the size of nuc prop struct */
  fprintf (stderr, "nucprops struct size: %i %i\n", bsz->nucpropssize, sizeof (Nuc_str));
  if (bsz->nucpropssize != sizeof (Nuc_str)) {
    read_bin_errmsg ("nucprops struct size does not match.", 1);
  }
}                               /* end of test section */

/************************************************/
/* This is the MAIN entry for reading a whole big block */
/* into a bigblock structure pointed at by bp   */
/* using the ctrl options located in the target file */
/************************************************/

void read_bin_blocks (BB_struct * bp, const char *fname)
{
  bsize bsz;
  int nread, sbnum;
  FILE *fp;
  long **sb_data_pos, this_pos, grain_pos;
  char head[256];


  fp = fopen (fname, "r");
  test_bin_block (fp, &bsz);

  /* if we made it this far, read in the BB_struct data */
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "start of bb struct data: %i \n", (this_pos));
#endif /*TB*/
    clearerr (fp);
  nread = fread (bp, bsz.bbsize, 1, fp);
  if (nread == 0)
    fprintf (stderr, "FERROR: %i\n", ferror (fp));

  fprintf (stderr, "Read %i bytes for bb struct\n", nread * bsz.bbsize);
  fprintf (stderr, "Bigblock Header: %s\n", bp->header);
  fprintf (stderr, "Bigblock Tailer: %s\n", bp->tailer);

  /*then the ctrl data */
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "start of ctrl struct data: %i \n", (this_pos));
#endif /*TB*/
    read_bin_ctrl (&(bp->ctrl), fp);
  /* read the global pressure value */

  /** \todo  change pressure to structue member */
  fread (&global_pressure, sizeof (CA_FLOAT), 1, fp);

  /* read in the location of the grain data */
  fread (&grain_pos, sizeof (long), 1, fp);

  /* now we know how big the bb is & how many subblocks so */
  /* create and zero an array of file position pointers */
  sb_data_pos = (long **) calloc (bp->ntsb, sizeof (long *));
  sb_contents_list = (Sb_tableofcontents *)calloc(bp->ntsb,sizeof(Sb_tableofcontents));
  for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {
    sb_data_pos[sbnum] = (long *) calloc (1, sizeof (long));
  }
  /* allocate the memory space */
  alloc_bb (bp);

  /* now read in the sb mask */
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "start of mask data: %i \n", (this_pos));
#endif /*TB*/
    fread (bp->sb_mask, sizeof (int), bp->ntsb, fp);

#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "start of table of contents : %i \n", (this_pos));
#endif /*TB*/
    fread (sb_contents_list, sizeof (Sb_tableofcontents), bp->ntsb, fp);
    for (sbnum=0;sbnum<bp->ntsb;sbnum++){
      write_sbtoc(stdout, (sb_contents_list+sbnum));
    }
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "end of table of contents : %i \n", (this_pos));
#endif /*TB*/

#ifdef TOC
 {
     int * buf;
     FILE * ofp;
  buf=(int *) calloc(bp->ncsb,sizeof(int));
  read_f_bin_igrain_data(fp,bp,buf,sb_contents_list[0].gr_pos);
  ofp=fopen("RAW_GR.bin","w");
  fwrite(buf,sizeof(int),bp->ncsb,ofp);
  free(buf);
  fclose(ofp);
  exit(0);
}
#endif TOC

  /* read the sub block member data */
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "start of sb struct data: %i \n", (this_pos));
#endif /*TB*/
    for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {
    if (bp->sb_mask[sbnum] == TRUE) {
      read_f_bin_sb (fp, bp, sbnum);
      fread (sb_data_pos[sbnum], sizeof (long), 1, fp);
    }
  }
  /* read the subblock (pointer) data arrays */
  for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {
    if (bp->sb_mask[sbnum] == TRUE) {
      this_pos = ftell (fp);
#ifndef LINUX
      if (this_pos != *(sb_data_pos[sbnum])) {
        fprintf (stderr, "ERROR:read_bin_blocks: data position does not match for subblock %i -- %i %i\n", sbnum, (int) this_pos,
                 (int) (*(sb_data_pos[sbnum])));
      }
#endif
#ifdef TB
      fprintf (stderr, "read_bin_blocks: data position for subblock %i -- %i %i\n", sbnum, (int) this_pos,
               (int) (*(sb_data_pos[sbnum])));
#endif /*TB*/
        read_f_bin_sb_data (fp, bp, sbnum);
#ifdef TB
      this_pos = ftell (fp);
      fprintf (stderr, "read_bin_blocks: pore position for subblock %i -- %i \n", sbnum, (this_pos));
#endif /*TB*/
        if (bp->ctrl->pore == TRUE)
        read_f_bin_sb_pores (fp, bp, sbnum);

    }
  } /* end of reading subblock data (pointer) arrays */

#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "grain pos %i\n", (this_pos));
#endif /*TB*/
    read_bin_grain_data (bp, fp);
  /* copy the partition coefficient values */
  /**  \todo  get rid of this awkward duplication of partcoef values */
  /** this is only used to create a scaled output slice but could be done in postrpocessing instead */

  bp->c_sol_values->part_coef = bp->mprops.gasprops.part_coef[0];
  bp->c_sol_alloy_values->part_coef = bp->mprops.alloyprops[0].part_coef[0];

  /* then free the pointers to the position data */
  for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {
    free (sb_data_pos[sbnum]);
  }
  free (sb_data_pos);
  free (sb_contents_list);

  fclose (fp);
}

/* rcs id subroutine.*/
/*Little subroutine to include the |rcs Id| in the object code.*/
/*This can also be called to print or include the Id in a file.*/

char const *rcs_id_read_blocks_c ()
{
  static char const rcsid[] = "$Id: read_blocks.c 1402 2008-11-20 15:36:41Z  $";

  return (rcsid);
}

/*
*/
