
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef BL_COMPRESS
#   include <zlib.h>
#endif

#include "machine.h"
#include "constants.h"
#include "blocks.h"
#include "read_ctrl.h"
#include "grain.h"
#include "pore.h"
#include "safeopen.h"
#include "sb_tableofcontents.h"


extern CA_FLOAT global_pressure;
Sb_tableofcontents * sb_contents_list;

void write_sbtoc(FILE * fp,Sb_tableofcontents * stocp ){
   int i;
   fprintf(fp," my_addresse: %li\n",(long) stocp-> my_addresse);
   fprintf(fp," c_elm_pos: %li\n",(long) stocp-> c_elm_pos);
   fprintf(fp," c_fs_pos: %li\n",(long) stocp-> c_fs_pos);
   fprintf(fp," sch_fs_pos: %li\n",(long) stocp-> sch_fs_pos);
   fprintf(fp," index_pos: %li\n",(long) stocp-> index_pos);
   fprintf(fp," gr_pos: %li\n",(long) stocp-> gr_pos);
   fprintf(fp," c_temp_pos: %li\n",(long) stocp-> c_temp_pos);
   fprintf(fp," curv_pos: %li\n",(long) stocp-> curv_pos);
   fprintf(fp," norm_x_pos: %li\n",(long) stocp-> norm_x_pos);
   fprintf(fp," norm_y_pos: %li\n",(long) stocp-> norm_y_pos);
   fprintf(fp," norm_z_pos: %li\n",(long) stocp-> norm_z_pos);
   fprintf(fp," c_sol_pos: %li\n",(long) stocp-> c_sol_pos);
   fprintf(fp," c_sol_alloy_pos: %li\n",(long) stocp-> c_sol_alloy_pos);
   fprintf(fp," c_eqv_alloy_pos: %li\n",(long) stocp-> c_eqv_alloy_pos);
   fprintf(fp," dc_d_pos: %li\n",(long) stocp-> dc_d_pos);
   fprintf(fp," dc_x_pos: %li\n",(long) stocp-> dc_x_pos);
   fprintf(fp," dc_y_pos: %li\n",(long) stocp-> dc_y_pos);
   fprintf(fp," dc_z_pos: %li\n",(long) stocp-> dc_z_pos);
   for (i=0;i<NSOLMAX;i++){
     fprintf(fp," c_eqv_poly_pos[%i]: %li\n",i,(long) stocp-> c_eqv_poly_pos[i]);
   }
   for (i=0;i<NPHAMAX;i++){
      fprintf(fp," c_fs_poly_pos[%i]: %li\n",i,(long) stocp-> c_fs_poly_pos[i]);
   }
   fprintf(fp," c_nuc_thresh_pos: %li\n",(long) stocp-> c_nuc_thresh_pos);
   fprintf(fp," c_surfp_pos: %li\n",(long) stocp-> c_surfp_pos);
   fprintf(fp," surf_xyz_pos: %li\n",(long) stocp-> surf_xyz_pos);
}

/* write a string tag to the file and place a carriage return after */
int wst (FILE * fp, const char *s)
{
  int nwrite = 0;
  if (s == NULL ){
      nwrite += fwrite("NULL",sizeof(char),4,fp);
  }else{
      nwrite += fwrite (s, sizeof (char), strlen (s), fp);
  }
  nwrite += fwrite ("\n", sizeof (char), 1, fp);
  return (nwrite);
}                               /* end of wst  -- write a string tag */

/****************************************************/
/* check integrity of sb masking                    */
/* see if the mask matches the actual subblock data */
/****************************************************/
void checkmask (BB_struct * bp)
{
  int i;

  for (i = 0; i < bp->ntsb; i++) {
    if (bp->sb_mask[i] != bp->sb[i]->open) {
      fprintf (stderr, "ERROR: checkmask: %i does not match\n", i);
      bp->sb_mask[i] = bp->sb[i]->open;
    }
  }
}                               /* end of checkmask */

/*********************************/
/* write the control information */
/*********************************/
int write_bin_ctrl (FILE * fp, Ctrl_str * cp)
{
  int i;
  int nwrite = 0;

  nwrite += fwrite (cp, sizeof (Ctrl_str), 1, fp);      /* ctrl struct data */

  /* Write the file names used as in ctrl file */
  nwrite += wst (fp, cp->fn_cap);
  nwrite += wst (fp, cp->fn_geo);
  nwrite += wst (fp, cp->fn_mat);
  nwrite += wst (fp, cp->fn_inp);
  nwrite += wst (fp, cp->fn_base);

  nwrite += wst (fp, cp->fn_solprops_gas);

  for (i=0;i<NSOLMAX;i++){
    nwrite += wst (fp,cp->fn_solprops_alloy[i]);
  }

  return (nwrite);
}

/*************************************************************/
/* write out the information stored directly in the subblock */
/*************************************************************/
void write_bin_sb (FILE * fp, BB_struct * bp, int sbnum)
{
  fwrite (bp->sb[sbnum], sizeof (SB_struct), 1, fp);
}

/*******************************************************/
/* write out all the information for a particular pore */
/*******************************************************/

void write_bin_pore (FILE * fp, PORE_str * cp, long *trad_pos)
{
  int i;
  p_c_list *bdy;
  p_c_node *node;

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
    bdy = cp->boundary;
    /*fwrite(bdy,sizeof(p_c_list),1,fp); */
    for (node = bdy->first; node != NULL; node = node->next) {
      fwrite (node, sizeof (p_c_node), 1, fp);
    }
    /* save the location of the data (t,rad) information */
    *trad_pos = ftell (fp);
    for (i = 0; i < N_T_LISTS; i++) {
      fwrite (cp->t_lists[i], sizeof (CA_FLOAT), cp->trad_last, fp);
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

/*************************************************************/
/* write out the pore list information                       */
/*************************************************************/
void write_bin_sb_pore (FILE * fp, BB_struct * bp, int sbnum)
{
  int pnum;
  long *pore_data_pos, *pore_trad_pos, info_start_pos, pore_end_pos;

  SB_struct *sp = bp->sb[sbnum];
  pore_data_pos = (long *) calloc (sp->Npores, sizeof (long));
  pore_trad_pos = (long *) calloc (sp->Npores, sizeof (long));

  /* write out the array of pore structures */
  fwrite (sp->porelist, sizeof (PORE_str), sp->Npores, fp);
  /* write an (empty) list of the data locations */
  info_start_pos = ftell (fp);
  fwrite (pore_data_pos, sizeof (long), sp->Npores, fp);
  /* and the t,rad info locations */
  fwrite (pore_trad_pos, sizeof (long), sp->Npores, fp);

  for (pnum = 0; pnum < sp->Npores; pnum++) {
    *(pore_data_pos + pnum) = ftell (fp);
    write_bin_pore (fp, (sp->porelist + pnum), (pore_trad_pos + pnum));
  }
  pore_end_pos = ftell (fp);
  /* rewind to the data location info position */
  fseek (fp, info_start_pos, SEEK_SET);
  /* write a list of the data locations */
  fwrite (pore_data_pos, sizeof (long), sp->Npores, fp);
  fwrite (pore_trad_pos, sizeof (long), sp->Npores, fp);
  /* reset to the end of the file Ihope */
  fseek (fp, 0, SEEK_END);
  free (pore_data_pos);
  free (pore_trad_pos);
}

#ifdef BL_COMPRESS
/* compress one array into a compressed buffer and write to the file */
/* return the size in bytes of the compressed data */
/* needs to be ifdef'd out to compile if zlib is not present */
int write_comp_array (void *datap, size_t vsize, int nmemb, FILE * fp)
{
  Bytef *comp_data;
  uLongf comp_size, datasize;

  datasize = (uLongf) (nmemb * vsize);

  /* calculate the maximum compressed buffer size */
  comp_size = (uLongf) (ceil (1.001 * (double) (datasize)) + 13);
  comp_data = (Bytef *) calloc (comp_size, sizeof (Bytef));

  /* call the routine from libz, this  requires -lz linker option */
  compress (comp_data, &comp_size, (Bytef *) datap, datasize);

  /* store the size of the compressed data */
  fwrite (&comp_size, sizeof (uLongf), 1, fp);

  /* then store the data */
  fwrite (comp_data, sizeof (Bytef), comp_size, fp);
  free (comp_data);
  return ((int) comp_size * sizeof (Bytef));
}
#endif /*BL_COMPRESS*/

/*************************************************************/
/* write out the information stored in pointer-arrays in the */
/* subblock                                                  */
/*************************************************************/

void write_bin_sb_data (FILE * fp, BB_struct * bp, int sbnum)
{
/*THUINET 03/05*/
  int iphs, isol, ele_num, ele_1;

/*end THUINET 03/05*/
#ifdef USE_ELM
  sb_contents_list[sbnum].c_elm_pos=ftello(fp);
  WRITEARRAY (bp->sb[sbnum]->c_elm, sizeof (int), bp->ncsb, fp);
#endif
  /* fraction solid */
  sb_contents_list[sbnum].c_fs_pos=ftello(fp);
  WRITEARRAY (bp->sb[sbnum]->c_fs, sizeof (CA_FLOAT), bp->ncsb, fp);
  /* schiel */
  if (bp->ctrl->scheil == TRUE) {
    sb_contents_list[sbnum].sch_fs_pos=ftello(fp);
    WRITEARRAY (bp->sb[sbnum]->sch_fs, sizeof (CA_FLOAT), bp->ncsb, fp);
  }

  /* cell index for Wei WANG */
  sb_contents_list[sbnum].index_pos=ftello(fp);
  WRITEARRAY (bp->sb[sbnum]->index, sizeof (int), bp->ncsb, fp);

  /* grain # */
  sb_contents_list[sbnum].gr_pos=ftello(fp);
  WRITEARRAY (bp->sb[sbnum]->gr, sizeof (int), bp->ncsb, fp);

  /* cell temperature -- used for all modes now! */
  sb_contents_list[sbnum].c_temp_pos=ftello(fp);
  WRITEARRAY (bp->sb[sbnum]->c_temp, sizeof (CA_FLOAT), bp->ncsb, fp);

  /* curvature for Wei WANG; Ness */
      /* protected by options rca 2007 01 15 */
      /**todo: finetune the options to avoid excessive memory allocation when the options are not used*/
  if ((bp->ctrl->curvature_3D != 0) || (bp->ctrl->curvature_2D != 0)) {
    sb_contents_list[sbnum].curv_pos=ftello(fp);
    WRITEARRAY (bp->sb[sbnum]->curv, sizeof (CA_FLOAT), bp->ncsb+2, fp);
    sb_contents_list[sbnum].norm_x_pos=ftello(fp);
    WRITEARRAY (bp->sb[sbnum]->norm_x, sizeof (CA_FLOAT), bp->ncsb+2, fp);
    sb_contents_list[sbnum].norm_y_pos=ftello(fp);
    WRITEARRAY (bp->sb[sbnum]->norm_y, sizeof (CA_FLOAT), bp->ncsb+2, fp);
    sb_contents_list[sbnum].norm_z_pos=ftello(fp);
    WRITEARRAY (bp->sb[sbnum]->norm_z, sizeof (CA_FLOAT), bp->ncsb+2, fp);
  }

  /* gas solute array */
  if (bp->ctrl->diffuse == TRUE) {
    sb_contents_list[sbnum].c_sol_pos=ftello(fp);
    WRITEARRAY (bp->sb[sbnum]->c_sol, sizeof (CA_FLOAT), bp->ncsb, fp);
  }
  /*ALLOY solute array */
  if ((bp->ctrl->diffuse_alloy == TRUE) || (bp->ctrl->particle == TRUE)) {
    sb_contents_list[sbnum].c_sol_alloy_pos=ftello(fp);
    WRITEARRAY (bp->sb[sbnum]->c_sol_alloy, sizeof (CA_FLOAT), bp->ncsb, fp);
    if (bp->ctrl->decentred_octahedron == TRUE) {
      sb_contents_list[sbnum].c_eqv_alloy_pos=ftello(fp);
      WRITEARRAY (bp->sb[sbnum]->c_eqv_alloy, sizeof (CA_FLOAT), bp->ncsb, fp);
    }
  }
  /* Wei WANG - dencentered octa arrays */
  if (bp->ctrl->decentred_octahedron == TRUE) {
    sb_contents_list[sbnum].dc_d_pos=ftello(fp);
    WRITEARRAY (bp->sb[sbnum]->dc_d, sizeof (CA_FLOAT), bp->ncsb, fp);
    sb_contents_list[sbnum].dc_x_pos=ftello(fp);
    WRITEARRAY (bp->sb[sbnum]->dc_x, sizeof (CA_FLOAT), bp->ncsb, fp);
    sb_contents_list[sbnum].dc_y_pos=ftello(fp);
    WRITEARRAY (bp->sb[sbnum]->dc_y, sizeof (CA_FLOAT), bp->ncsb, fp);
    sb_contents_list[sbnum].dc_z_pos=ftello(fp);
    WRITEARRAY (bp->sb[sbnum]->dc_z, sizeof (CA_FLOAT), bp->ncsb, fp);
  }

  /* multicomponent for L. THUINET 03/05 */
  ele_num = bp->ctrl->NUM_COMP; /* number of elements in the alloy */
  ele_1 = ele_num - 1;

  if ((bp->ctrl->decentred_octahedron == TRUE) && (bp->ctrl->diffuse_alloy_poly == TRUE)) {
    for (isol = 0; isol < ele_1; isol++) {      /*loop on solutes by Ludovic THUINET */
      sb_contents_list[sbnum].c_eqv_poly_pos[isol]=ftello(fp);
      WRITEARRAY (bp->sb[sbnum]->c_eqv_poly[isol], sizeof (CA_FLOAT), bp->ncsb, fp);
    }
    for (iphs = 0; iphs < bp->ctrl->NUM_PHS; iphs++) {  /*loop on solutes by Ludovic THUINET */
      sb_contents_list[sbnum].c_fs_poly_pos[iphs]=ftello(fp);
      WRITEARRAY (bp->sb[sbnum]->c_fs_poly[iphs], sizeof (CA_FLOAT), bp->ncsb, fp);
    }
  }
  /*end multicomponent for THUINET */

  if ((bp->ctrl->block_nuc == TRUE)) {  /*nuc threshold array */
    sb_contents_list[sbnum].c_nuc_thresh_pos=ftello(fp);
    WRITEARRAY (bp->sb[sbnum]->c_nuc_thresh, sizeof (CA_FLOAT), bp->ncsb, fp);
  }

  /* the surface cell location pointers */
  sb_contents_list[sbnum].c_surfp_pos=ftello(fp);
  WRITEARRAY (bp->sb[sbnum]->surface.c_surfp, sizeof (int), bp->sb[sbnum]->surface.n_alloc, fp);
  sb_contents_list[sbnum].surf_xyz_pos=ftello(fp);
  WRITEARRAY (bp->sb[sbnum]->surface.surf_xyz, sizeof (int[3]), bp->sb[sbnum]->surface.n_alloc, fp);

}

void write_bin_grain_data (BB_struct * bp, FILE * fp)
{
  int i;

  for (i = 1; i <= bp->nprops.ngr; i++) {
    fwrite (bp->gr[i], sizeof (Ind_grain), 1, fp);
  }
}

/**************************************************************/
/* this is the MAIN entry point for this group of routines    */
/* that is, if you want to write a full big-block binary file */
/* call this routine                                          */
/**************************************************************/

void write_bin_blocks (BB_struct * bp)
{
  int sbnum;
  const char revision[] = BLOCKREV;
  const char bbrevision[] = BIGBLOCKREV;
  const char sbrevision[] = SUBBLOCKREV;
  const char crevision[] = CTRLREV;
  const char grevision[] = GRAINREV;
  const char prevision[] = POREREV;
  const char nrevision[] = NUCPROPSREV;
  const char endcomment[39] = "End of CA-FD big-block restart file.\n\0";
  char fname[MAX_STRING_LEN];
  char checkname[MAX_STRING_LEN];
  /*
  char command[ ( MAX_STRING_LEN + MAX_STRING_LEN + 10 )];
  */
  int bszsize, bbsize, fsize, sbsize, csize, gsize, revsize;
  int nwrite = 0, ntot = 0;
  bsize bsz;                    /* store all the structure sizes to check compatability */
  FILE *fp;
  struct stat filestat;
  int statreturn=0,try=0,myerror=0;
  long **sb_data_pos, **sb_info_pos, this_pos, grain_pos, grain_info_pos,sb_toc_pos;

    /***********************************************/
    /* save the size of the various structures     */
    /* to determine if there is an incompatibility */
    /***********************************************/

    /*************************************************************/
    /* the precision of the CA_FLOAT variable is compiled in and is */
    /* independant of the version.                               */
    /*************************************************************/

  /* create and zero an array of file position pointers */
  sb_info_pos = (long **) calloc (bp->ntsb, sizeof (long *));
  sb_data_pos = (long **) calloc (bp->ntsb, sizeof (long *));
  for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {
    sb_data_pos[sbnum] = (long *) calloc (1, sizeof (long));
    sb_info_pos[sbnum] = (long *) calloc (1, sizeof (long));
  }
  sb_contents_list = (Sb_tableofcontents *)calloc(bp->ntsb,sizeof(Sb_tableofcontents));

  sprintf (fname, "BLOCK_%s_t%06d.%s", bp->ctrl->fn_base, bp->step, BL_EXT);
  bsz.fsize = sizeof (CA_FLOAT);
  bsz.bbsize = sizeof (BB_struct);
  bsz.sbsize = sizeof (SB_struct);
  bsz.csize = sizeof (Ctrl_str);
  bsz.gsize = sizeof (Ind_grain);
  bsz.psize = sizeof (P_str);
  bsz.PORE_size = sizeof (PORE_str);
  bsz.pcnsize = sizeof (p_c_node);
  bsz.pclsize = sizeof (p_c_list);
  bsz.nucpropssize = sizeof (Nuc_str);

    /**********************************************/
    /* also use the revision number from cvs      */
    /* to help with changed versions              */
    /**********************************************/
#ifdef TB
  fprintf (stderr, "Verbose info for write_blocks\n************************\n");
#endif /*TB*/
  /* check the file, if it exists use a different filename */
    errno=0;
    try=0;
    while ( ( statreturn = stat(fname,&filestat)) == 0){;
       if (try >= MAX_REWRITES){
	  fprintf(stderr,"ERROR:%s: Maximum restart file limit exceeded! Exiting to avoid filling up the disk!\n",__func__);
	  exit(100);
       }
       fprintf(stderr,"%s exists! Saving to a new restart file name\n",fname);
       try++;
       snprintf (fname,MAX_STRING_LEN, "BLOCK_%s_r%02d_t%06d.%s", bp->ctrl->fn_base,try, bp->step, BL_EXT);
    } 

  fp = fopen (fname, "w");

  revsize = strlen (revision);
  fwrite (revision, sizeof (char), revsize, fp);        /*blocks.h revision ID string */
  fwrite ("\n", sizeof (char), 1, fp);

  revsize = strlen (bbrevision);
  fwrite (bbrevision, sizeof (char), revsize, fp);      /*bigblock.h revision ID string */
  fwrite ("\n", sizeof (char), 1, fp);

  revsize = strlen (sbrevision);
  fwrite (sbrevision, sizeof (char), revsize, fp);      /*subblock.h revision ID string */
  fwrite ("\n", sizeof (char), 1, fp);

  revsize = strlen (crevision);
  fwrite (crevision, sizeof (char), revsize, fp);       /*ctrl revision ID string */
  fwrite ("\n", sizeof (char), 1, fp);

  revsize = strlen (grevision);
  fwrite (grevision, sizeof (char), revsize, fp);       /*grain revision ID string */
  fwrite ("\n", sizeof (char), 1, fp);

  revsize = strlen (prevision);
  fwrite (prevision, sizeof (char), revsize, fp);       /*pore revision ID string */
  fwrite ("\n", sizeof (char), 1, fp);

  revsize = strlen (nrevision);
  fwrite (nrevision, sizeof (char), revsize, fp);       /*nucprops revision ID string */
  fwrite ("\n", sizeof (char), 1, fp);

  /* write out the size values */
  bszsize = sizeof (bsize);
  fwrite (&bszsize, sizeof (int), 1, fp);       /* size of the block size struct */
  fwrite (&bsz, sizeof (bsize), 1, fp); /* block size structure   */

  /*write out the big block data */
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "start of bb struct data: %i \n", (this_pos));
#endif /*TB*/
    nwrite = fwrite (bp, bsz.bbsize, 1, fp);    /* bb struct data */
#ifdef TB
  fprintf (stderr, "nwrite =%i\n", nwrite);
  this_pos = ftell (fp);
  fprintf (stderr, "start of ctrl struct data: %i \n", (this_pos));
#endif /*TB*/
    /* write the control structure data */
    nwrite = write_bin_ctrl (fp, bp->ctrl);     /* ctrl struct data, strings */
  /* write the global pressure */

  /**  \todo  change pressure to structure member -- general - porosity */

  nwrite = fwrite (&global_pressure, sizeof (CA_FLOAT), 1, fp);

  /* make some space for the location of grain data */
  grain_info_pos = ftell (fp);
  fwrite (&grain_pos, sizeof (long), 1, fp);

  /* test and write out the mask of open subblocks */
#ifdef TB
  fprintf (stderr, "nwrite =%i\n", nwrite);
  this_pos = ftell (fp);
  fprintf (stderr, "start of mask data: %i \n", (this_pos));
#endif /*TB*/
    checkmask (bp);
  nwrite = fwrite (bp->sb_mask, sizeof (int), bp->ntsb, fp);
#ifdef TB
  fprintf (stderr, "nwrite =%i\n", nwrite);
  this_pos = ftell (fp);
  fprintf (stderr, "end of mask data: %i \n", (this_pos));
#endif /*TB*/

  this_pos = ftell (fp);
  sb_toc_pos = this_pos;

#ifdef TB
  fprintf (stderr, "start of sb table of contents : %i \n", (this_pos));
#endif /*TB*/

  nwrite = fwrite(sb_contents_list,sizeof(Sb_tableofcontents),bp->ntsb,fp);

#ifdef TB
  fprintf (stderr, "nwrite =%i\n", nwrite);
  this_pos = ftell (fp);
  fprintf (stderr, "end of sb table of contents : %i \n", (this_pos));
#endif /*TB*/


  /* write the subblock data for the specified subblocks */
#ifdef TB
  fprintf (stderr, "nwrite =%i\n", nwrite);
  this_pos = ftell (fp);
  fprintf (stderr, "start of sb struct data: %i \n", (this_pos));
#endif /*TB*/
    for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {
    if (bp->sb_mask[sbnum] == TRUE) {
      write_bin_sb (fp, bp, sbnum);
      /* store this location to write the data pointer */
      /* into once we know where it is! */
      *(sb_info_pos[sbnum]) = ftell (fp);
#ifdef TB
      fprintf (stderr, "info_pos: sb %i %i\n", sbnum, (int) (*(sb_info_pos[sbnum])));
#endif /*TB*/
        /* write some zeros to make space for the data pointer */
        fwrite (sb_data_pos[sbnum], sizeof (long), 1, fp);
    }
  }

  /* write the data arrays for the specified subblocks   */
  for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {
    if (bp->sb_mask[sbnum] == TRUE) {
      /* store the pointer to this location */
      *(sb_data_pos[sbnum]) = ftell (fp);
#ifdef TB
      fprintf (stderr, "data_pos: sb %i %i\n", sbnum, (*(sb_data_pos[sbnum])));
#endif /*TB*/
        /* write the data arrays */
        write_bin_sb_data (fp, bp, sbnum);
      /* write the pore data */
#ifdef TB
      this_pos = ftell (fp);
      fprintf (stderr, "pore pos: sb %i %i\n", sbnum, (this_pos));
#endif /*TB*/
        if (bp->ctrl->pore == TRUE)
        write_bin_sb_pore (fp, bp, sbnum);
#ifdef TB
      this_pos = ftell (fp);
      fprintf (stderr, "end sb pos: sb %i %i\n", sbnum, (this_pos));
#endif /*TB*/
    }
  }

  /* write the grain data */
  grain_info_pos = ftell (fp);
#ifdef TB
  this_pos = ftell (fp);
  fprintf (stderr, "grain pos %i\n", (this_pos));
#endif /*TB*/
    write_bin_grain_data (bp, fp);

  /* now go back and write the data pointer information */
  /* into the stored locations. */
  for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {
    if (bp->sb_mask[sbnum] == TRUE) {
      fseek (fp, *(sb_info_pos[sbnum]), SEEK_SET);
      fwrite (sb_data_pos[sbnum], sizeof (long), 1, fp);
    }
  }
  /* write the grain location */
  fseek (fp, grain_info_pos, SEEK_SET);
  fwrite (&grain_pos, sizeof (long), 1, fp);


#ifdef TB

  fprintf (stderr, "start of sb table of contents : %i \n", (sb_toc_pos));
#endif /*TB*/

  fseek(fp,sb_toc_pos,SEEK_SET);
  nwrite = fwrite(sb_contents_list,sizeof(Sb_tableofcontents),bp->ntsb,fp);

#ifdef TB
  fprintf (stderr, "nwrite =%i blocks, %i bytes\n", nwrite,nwrite*sizeof(Sb_tableofcontents));
  this_pos = ftell (fp);
  fprintf (stderr, "end of sb table of contents : %i \n", (this_pos));
#endif /*TB*/

  /* then free the pointers to the position data */
  for (sbnum = 0; sbnum < bp->ntsb; sbnum++) {
    free (sb_data_pos[sbnum]);
    free (sb_info_pos[sbnum]);
  }
  free (sb_data_pos);
  free (sb_info_pos);
  for (sbnum=0;sbnum<bp->ntsb;sbnum++){
   write_sbtoc(stdout,(sb_contents_list+sbnum) );
  }
  free (sb_contents_list);

  /* write an end comment to the file */
  fflush (fp);
  wst (fp, bp->tailer);
  wst (fp, endcomment);
  /* and close the file */
  fclose (fp);
  /* create the latest step's checkpoint block file */
  /* it's a hard link */
  /* delete the existing checkpoint if any */
  fprintf (stderr, "write_bin_blocks: Linking checkpoint file.\n");
  snprintf (checkname,MAX_STRING_LEN, "%s.%s",DEFAULT_CHECKPOINT_NAME, BL_EXT);

  statreturn=0;
  statreturn = unlink(checkname);
  myerror=errno;
  if ( myerror != ENOENT &&  statreturn == -1 ) {
     static int errct=0;
     fprintf(stderr,"ERROR:__func__: unexpected file problem with the checkpoint file!\n");
     fprintf(stderr,"%s: %s\n",checkname,strerror(myerror));
     if (errct >= 1000){
          fprintf(stderr,"ERROR:__func__: error limit exceeded for this error at line %s!\n",__func__, __LINE__);
	  exit(200);
     }
     errct++;
  }

  /* create a new checkpoint link to the latest block file */
  statreturn=0;
  statreturn = link(fname,checkname);
  myerror=errno;
  if (statreturn == -1 ) {
     static int errct=0;
     fprintf(stderr,"ERROR:__func__: unexpected file problem with the checkpoint or block file!\n");
     fprintf(stderr,"%s: %s\n",checkname,strerror(myerror));
     if (errct >= 1000){
          fprintf(stderr,"ERROR:__func__: error limit exceeded for this error at line %s!\n",__func__, __LINE__);
	  exit(200);
     }
     errct++;
  }
  fprintf (stderr, "write_bin_blocks: Finished writing restart bigblock file.\n");
}                               /* end of write_bin_blocks */

/* subroutine to return rcs-id string */
char const *rcs_id_write_blocks_c ()
{
  static char const rcsid[] = "$Id: write_blocks.c 1402 2008-11-20 15:36:41Z  $";

  return (rcsid);
}

/*
*/
