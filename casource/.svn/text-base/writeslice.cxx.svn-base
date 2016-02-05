#include <stdio.h>
/*another test rsybnc */
/*another test cvs - ssh*/
#include <stdlib.h>
#include <math.h>
#include "gd.h"
#include <string.h>

#include "machine.h"
#include "blocks.h"
#include "colour.h"
#include "interp.h"
#include "pore_routines.h"

class block:public BB_struct {
public:
};

char make_char_mod(void * info){
   return (char) (10 + *( (int *) info)%(245));
}

char make_char_float (void * info){
   float imax;
   imax = 100;
   return (char) (( ( *(float *)info) / imax) * 245 ); 
}

int      write_slice(BB_struct *bp,Value_struct *val, char (*make_char) (void *) ) {
   char      command[MAX_STRING_LEN];
   int      nsbx, nsby, nsbz;    /* tmp counters */
   int      bz, cz, by, cy, bx, cx;    /* tmp counters */
   int      sbnum, cnum, ncx, ncy, ncz;    /* tmp counters */
   FLOAT tmp;   /* tmp FLOAT var.*/
   FILE * fp;   /* tmp filehandle */
   FLOAT * info;
   char      *ip, *imgdata;
   char bsize[3]; /* x y and z size of block */
   /************************************************/
   /* open output file to write binary image to   */
   /* fixed filename for now         */
   /************************************************/
   sprintf(command, "BL_%s%st%09d.bin",val->id_string, bp->ctrl->fn_base, bp->step);
   if ((fp = fopen(command, "w")) == NULL) {
      fprintf(stderr, "Error: can't open output file [%s]\n", command);
      return (1);
   }
   nsbx = bp->nsb[0]; 
   nsby = bp->nsb[1]; 
   nsbz = bp->nsb[2];
   ncx = bp->nc[0]; 
   ncy = bp->nc[1]; 
   ncz = bp->nc[2];

   /* limitation for AVS : max 255 cells on a side */
   bsize[0] = (char) ncx;
   bsize[1] = (char) ncy;
   bsize[2] = (char) ncz;
   if (!(imgdata = (char *) calloc(ncx+1, sizeof(char)))) { 
      fprintf(stderr, "image: pdata malloc failed\n"); 
      return 1; 
   }

   fwrite(bsize, sizeof(char), 3, fp);/* prepend volume size as used by avs */
                                      /*one byte each for x,y,and z information */
   sbnum = 0;
   for (bz = 0; bz < nsbz; bz++) {      /* loop sb's in z direction */
      for (cz = 0; cz < ncz; cz++) {   /* loop cellss in z direction */
       for (by = 0; by < nsby; by++) {      /* loop sb's in y direction */
          for (cy = 0; cy < ncy; cy++) {   /* loop cells in y direction */
             for (bx = 0; bx < nsbx; bx++) {   /* loop sb's in x direction */

              /* compose and write out a line of this subblock */
              sbnum = bx + by * nsbx + bz * (nsbx * nsby);
              ip = imgdata;
              cnum = cy * ncx + cz * (ncx * ncy);
              info = (val->block_array[sbnum])+cnum ;

              for (cx = 0; cx < ncx; cx++) {   /* loop cells in x */
                  *(ip++) = (*make_char)((void *)(info++));
              }

              fwrite(imgdata, sizeof(char), ncx, fp);
              /* finished this line */

             } /* sb's in x */
          } /* cell's in y */
       } /* sb's in y */
      } /* cells's in z */
   } /* sb's in z */

   free (imgdata);
   fclose(fp);
   return (0);
} /* end of write_bb_conc subroutine */
