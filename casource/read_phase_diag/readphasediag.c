/*$Id: readphasediag.c 892 2006-03-10 15:24:59Z rcatwood $*/
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
/*This file was created at Wed Jul 14 18:26:00 BST 2004 by rcatwood on hive.beowulf.cluster */
#define READPHASEDIAG_C_REV "readphasediag.c $Id: readphasediag.c 892 2006-03-10 15:24:59Z rcatwood $"
/* Include the common header for this project */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../machine.h"
#include "phasediag.h"

#ifdef BL_COMPRESS
#   define READARRAY read_comp_array
#   define BL_EXT "blz"
#else
#   define READARRAY fread
#   define BL_EXT "blk"
#endif

extern int read_comp_array(void * datap,size_t vsize,int nmemb,FILE * fp);
extern void alloc_phasediag(Phdiag * pdp);
extern void alloc_tiedata(Reg * regp);


size_t checkread(void * dest,size_t size,size_t nmemb,FILE * fp){
     size_t nr;

     nr = fread(dest,size,nmemb,fp);
    if (nr != nmemb){
       fprintf(stderr,"ERROR:read_phasediag: Read error, expecting %i got %i\n",nmemb,nr);
       exit(0);
    }

    return(nr);
}

size_t read_region(Phdiag * pdp,int num,FILE * fp){
   int i;
   Pd_num_t phnum;
   size_t nr=0;

   pdp->regions[num].phases = (Pha *) calloc(pdp->regions[num].np,sizeof(Pha));

   for(i=0;i<pdp->regions[num].np;i++){
      nr += checkread(&phnum,sizeof(Pd_num_t),1,fp);
      pdp->regions[num].phases[i] = pdp->phases[phnum];
   }
   if (pdp->regions[num].ntie > 0 ){
      alloc_tiedata( &(pdp->regions[num]) );
      nr += READARRAY(pdp->regions[num].tiedata,sizeof(Pd_tri),pdp->regions[num].ntie,fp);
   }
   return(nr);
}

size_t read_comp(Phdiag * pdp,int num,FILE * fp){

return(0);
}

size_t read_phase(Phdiag * pdp,int num,FILE * fp){

return(0);
}
   
/* read the data referred to by pointers in the tslice struct */
size_t read_tslice(Phdiag * pdp, int num,FILE * fp){
   int i;
   size_t nr = 0;
   CA_FLOAT * in_array,*inp;
   int n_array;
   in_array = (CA_FLOAT *) calloc(pdp->tslices[num].ndata * PHASEDIAG_NPERPOINT, sizeof(CA_FLOAT));
   inp=in_array;
   n_array=0;


   pdp->tslices[num].regdata = (Pd_num_t *) calloc(pdp->tslices[num].ndata,sizeof(Pd_num_t));
   pdp->tslices[num].data = (CA_FLOAT **) calloc(pdp->tslices[num].ndata,sizeof(CA_FLOAT *));

   nr  = READARRAY(pdp->tslices[num].regdata,sizeof(Pd_num_t),pdp->tslices[num].ndata,fp);
   /*TODO: check the size */
   nr = checkread(&n_array,sizeof(size_t),1,fp);

   nr = READARRAY(in_array,sizeof(CA_FLOAT),n_array,fp);
   /*TODO: check the size */

   for (i=0;i<pdp->tslices[num].ndata;i++){
      if ( pdp->regions[pdp->tslices[num].regdata[i]].np == 2 ){
         pdp->tslices[num].data[i] = (CA_FLOAT *) calloc(PHASEDIAG_NPERPOINT,sizeof(CA_FLOAT));
         memcpy(pdp->tslices[num].data[i],inp,PHASEDIAG_NPERPOINT * sizeof(CA_FLOAT));
         inp += PHASEDIAG_NPERPOINT ;
      }
   }

free(in_array);
return (nr);
}
size_t read_phasediag(Phdiag * pdp, const char * filename){
    FILE * fp;
    int i;
    size_t nr=0;
    size_t nread=0;


    fp=fopen(filename,"r");


    nr = checkread(pdp,sizeof(Phdiag),1,fp);

    alloc_phasediag(pdp);

    nr = checkread(pdp->comp, sizeof(Comp),pdp->ncomp,fp);
       for(i=0;i<pdp->ncomp;i++){
           nr = read_comp(pdp,i,fp);
       }

    nr = checkread(pdp->phases, sizeof(Pha),pdp->n_phases,fp);
       for(i=0;i<pdp->n_phases;i++){
           nr = read_phase(pdp,i,fp);
       }

    nr = checkread(pdp->regions,sizeof(Reg),pdp->n_regions,fp);
       for (i=0;i<pdp->n_regions;i++){
          nr = read_region(pdp,i,fp);
       }
    nr = READARRAY(pdp->tslices,sizeof(Tslice),pdp->n_tslice,fp);
       #ifndef BL_COMPRESS
          nr *= sizeof(Tslice);
       #endif
       if (nr != pdp->n_tslice * sizeof(Tslice) ) {
          fprintf(stderr,"ERROR:read_phasediag: Read error, expecting %li got %li\n",pdp->n_tslice,nr);
          exit(0);
       }
       for (i=0;i<pdp->n_tslice;i++){
          nr = read_tslice(pdp,i,fp);
       }


       

return(nread);
}
    


