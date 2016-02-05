/*$Id: writephasediag.c 892 2006-03-10 15:24:59Z rcatwood $*/
/****************************************************************/
/*   Copyright (c) 1998 - 2004 Dept. of Materials, ICSTM        */
/*   All Rights Reserved                                        */
/*   THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE OF ICSTM       */
/*   The copyright notice above does not evidence any           */
/*   actual or intended publication of such source code,        */
/*   and is an unpublished work by Dept. of Materials, ICSTM.   */
/*   This material contains CONFIDENTIAL INFORMATION that       */
/*   is the property of Imperial College. Any use,              */
/*   duplication or disclosure not specifically authorized      */
/*   by Imperial College is strictly prohibited.                */
/****************************************************************/
/* This code is part of the umats routines developed at in the  */
/* Materials Processing Group, Dept. of Materials, ICSTM.       */
/*      email p.d.lee or r.atwood @ic.ac.uk for details         */
/****************************************************************/
/*This file was created at Wed Jul 14 15:57:05 BST 2004 by rcatwood on hive.beowulf.cluster */
#define WRITEPHASEDIAG_C_REV "writephasediag.c $Id: writephasediag.c 892 2006-03-10 15:24:59Z rcatwood $"
/* Include the common header for this project */
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "../machine.h"
#include "phasediag.h"

#ifdef BL_COMPRESS
#   define WRITEARRAY write_comp_array
#   define PD_EXT "pdz"
#else
#   define WRITEARRAY fwrite
#   define PD_EXT "pdg"
#endif
const char wpd_c_rev[] = WRITEPHASEDIAG_C_REV ;

extern int write_comp_array(void * datap,size_t vsize,int nmemb,FILE * fp);

/* write the data referred to by pointers in the region struct*/
size_t write_region(Phdiag * pdp, int num,FILE * fp){
   int i;
   int written;
   Pd_num_t phnum;
   size_t nwrite = 0;
   /* save the phase number, not the pointer address */
   for (i=0;i<pdp->regions[num].np;i++){
       phnum = pdp->regions[num].phases[i].my_num;
       nwrite += fwrite(&phnum,sizeof(Pd_num_t),1,fp) * sizeof(Pd_num_t);
   }
   
   /* write the array of tie triangle data */
   if (pdp->regions[num].ntie > 0 ){
      written = WRITEARRAY(pdp->regions[num].tiedata,sizeof(Pd_tri),pdp->regions[num].ntie,fp);
     #ifndef BL_COMPRESS
      written *= sizeof(Pd_tri);
     #endif
      nwrite += written; 
   }
return (nwrite);
}

/* write the data referred to by pointers in the component struct */
size_t write_comp(Phdiag * pdp, int num,FILE * fp){
   /*
   int i;
   Pd_num_t cnum;
   */
   size_t nwrite = 0;
   /*except there is no info stored here just now*/
   /* just a stub */
return (nwrite);
}

/* write the data referred to by pointers in the phase struct */
size_t write_phase(Phdiag * pdp, int num,FILE * fp){
   /*
   int i;
   Pd_num_t phnum;
   */
   size_t nwrite = 0;
   /*except there is no info stored here just now*/
   /* just a stub */
return (nwrite);
}

/* write the data referred to by pointers in the tslice struct */
size_t write_tslice(Phdiag * pdp, int num,FILE * fp){
   int i;
   size_t nwrite = 0,written=0;
   CA_FLOAT * out_array,*outp;
   size_t n_array;
   out_array = (CA_FLOAT *) calloc(pdp->tslices[num].ndata * PHASEDIAG_NPERPOINT, sizeof(CA_FLOAT));
   outp=out_array;
   n_array=0;


   written = WRITEARRAY(pdp->tslices[num].regdata,sizeof(Pd_num_t),pdp->tslices[num].ndata,fp);
  #ifndef BL_COMPRESS
   written *= sizeof(Pd_num_t);
  #endif
   nwrite += written; 

   for (i=0;i<pdp->tslices[num].ndata;i++){
      if ( pdp->regions[pdp->tslices[num].regdata[i]].np == 2 ){
         memcpy(outp,pdp->tslices[num].data[i],PHASEDIAG_NPERPOINT * sizeof(CA_FLOAT));
         outp += PHASEDIAG_NPERPOINT ;
         n_array += PHASEDIAG_NPERPOINT ;
      }
   }

   written = fwrite(&(n_array),sizeof(size_t),1,fp);

   nwrite += written * sizeof(size_t);

   written = WRITEARRAY(out_array,sizeof(CA_FLOAT),n_array,fp);
  #ifndef BL_COMPRESS
   written *= sizeof(CA_FLOAT);
  #endif
   nwrite += written; 

free(out_array);
return (nwrite);
}

/* write out a phase diagram data structrue in binary form */
size_t write_phasediag(Phdiag * pdp, const char * filename){
   FILE * fp;
   char fname[255];
   int i;
   size_t nwrite = 0;
   size_t written = 0;

   snprintf(fname,255,"%s.%s",filename,PD_EXT);


   fp = fopen(fname,"w");

   /* first the header string */
   /* included in the data struct */
   /* second, the phase diag structure data */
   nwrite += fwrite(pdp,sizeof(Phdiag),1,fp) * sizeof(Phdiag);

   /* third, the phase data */
   nwrite += fwrite(pdp->comp,sizeof(Comp),pdp->ncomp,fp) * sizeof(Comp);
      for (i=0;i<pdp->ncomp;i++){
          nwrite += write_comp(pdp,i,fp);
      }

   nwrite += fwrite(pdp->phases,sizeof(Pha),pdp->n_phases,fp) * sizeof(Pha);
      /* followed by the phase info data for all regions */
      for (i=0;i<pdp->n_phases;i++){
          nwrite += write_phase(pdp,i,fp);
      }
   /* fourth, the region data */
   nwrite += fwrite(pdp->regions,sizeof(Reg),pdp->n_regions,fp) * sizeof(Reg);

      /* followed by the region-phase array data for all regions */
      for (i=0;i<pdp->n_regions;i++){
          nwrite += write_region(pdp,i,fp);
      }

   /* Fifth, the isothermal sections data */
   written = WRITEARRAY(pdp->tslices,sizeof(Tslice),pdp->n_tslice,fp);
  #ifndef BL_COMPRESS
   written *= sizeof(Tslice);
  #endif
   nwrite += written;


   /* followed by its sub-data */
      for (i=0;i<pdp->n_tslice;i++){
          written = write_tslice(pdp,i,fp);

          nwrite += written;
      }
   /* last, the tailer string */
   /* included in the data struct */
   written = fwrite(pdp->tailer,sizeof(char),255,fp);
   nwrite += written;
   fclose(fp);
return (nwrite);
}


/* create a volume image file of the phase diagram */
/* TODO: account for non-constant boundaries and number of points */
void two_d_image(Phdiag * pdp, FILE * fp){
    int j;
    for(j=0;j<pdp->n_tslice;j++){
        fwrite(pdp->tslices[j].regdata,sizeof(Pd_num_t),pdp->tslices[j].n_csteps[0],fp);
    }

}
void three_d_image(Phdiag * pdp, FILE * fp){
    int j;
    for(j=0;j<pdp->n_tslice;j++){
        fwrite(pdp->tslices[j].regdata,sizeof(Pd_num_t),pdp->tslices[j].ndata,fp);
    }
}

/*write the phase diagram region information as a volume graphic file */
/* include a header with the size and origin information */
void vol_phasediag(Phdiag * pdp, const char * filename){
     char raw_filename[MAX_STRING_LEN],hdr_filename[MAX_STRING_LEN];
     int xs,ys,zs;
     double xo,yo,zo;
     FILE * fp;

      xs = pdp->tslices[0].n_csteps[0];
      ys = pdp->tslices[0].n_csteps[1];
      zs = pdp->n_tslice;
      xo = pdp->tslices[0].cmin[0];
      yo = pdp->tslices[0].cmin[1];
      zo = pdp->tmin;

     snprintf(raw_filename,MAX_STRING_LEN,"%s.raw",filename);
     snprintf(hdr_filename,MAX_STRING_LEN,"%s_hdr.txt",filename);


      fp = fopen(raw_filename,"w");
      if (pdp->ncomp == 2){
         two_d_image(pdp,fp);
         ys = 1;
      }else if (pdp->ncomp == 3){
         three_d_image(pdp,fp);
      }else{
         fprintf(stderr,"vol_phasediag: Sorry, cannot handle %i dimensional diagrams\n",pdp->ncomp);
      }
      fclose (fp);
      /* write a header file for the raw volume */
      fp= fopen(hdr_filename,"w");

      fprintf(fp,"infile %s\n",raw_filename);
      fprintf(fp,"ncomp %i\n",pdp->ncomp);
      fprintf(fp,"datasize %i %i %i\n",xs,ys,zs);
      fprintf(fp,"origin %g %g %g \n",xo,yo,zo);
      fclose(fp);
return;
}










