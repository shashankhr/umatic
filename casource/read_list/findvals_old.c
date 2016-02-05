/*RCS Id:$Id: findvals_old.c 887 2006-03-01 18:21:01Z rcatwood $*/
#include <stdio.h>
#include <string.h>
#include "../machine.h"
#include "../constants.h"
#include "csv.h"
#include "../fidap.h"

/* copy the selected column into the FGrid_str array T */
/* and find the min,max temperature */
void find_temp(CsvFloat * the_float,FGrid_str * fg,int tfield){
int i;
fg->Tmax = -400;
fg->Tmin = 10000;
    for (i=0;i<the_float->line_count;i++){
       fg->Fidap_T[i] = the_float->data[i][tfield];
       if (fg->Fidap_T[i] > fg->Tmax) fg->Tmax = fg->Fidap_T[i];
       if (fg->Fidap_T[i] < fg->Tmin) fg->Tmin = fg->Fidap_T[i];
    }
}

/* Find out how many points in each of z and r directions, and copy the */
/* values into the FGrid_str member arrays z and r */
/* hence determine the number of nodes (should equal number of lines input)*/
void find_z_r(CsvFloat * the_float,FGrid_str * fg,int zfield,int rfield){
   int i;
   int rflag=0;
   CA_FLOAT this_r, this_z;
   CA_FLOAT *z,*r;
   int nr=0,nz=0,maxr=INIT_ARRAY,maxz=INIT_ARRAY,add=ADD_ARRAY;
   

   z=(CA_FLOAT *)malloc(maxz*sizeof(CA_FLOAT));
   r=(CA_FLOAT *)malloc(maxr*sizeof(CA_FLOAT));
   this_z = *z = the_float->data[0][zfield];
   this_r = *r = the_float->data[0][rfield];
   for (i=1;i<the_float->line_count;i++){
      this_r = the_float->data[i][rfield];
      if ( r[nr] != this_r){
         rflag = 1;
         nr++;
         if (nr >= maxr) {
            maxr += add;
            r = (CA_FLOAT *) realloc(r,maxr*sizeof(CA_FLOAT));
         }
         r[nr]=this_r;
      }
      if (rflag == 0 ){
         this_z = the_float->data[i][zfield];
         if ( z[nz] != this_z){
            nz++;
            if (nz >= maxz) {
               maxz += add;
               z = (CA_FLOAT *) realloc(z,maxz*sizeof(CA_FLOAT));
            }
            z[nz]=this_z;
         }
      }
   }
   fg->r = r;
   fg->z = z;
   fg->nr = nr+1;
   fg->nz = nz+1;
   fg->nnodes = fg->nr*fg->nz;
}


/*
RCS Log:$Log$
RCS Log:Revision 11.1  2006/03/01 18:21:00  rcatwood
RCS Log:Merging polycomponent and gas with meltback
RCS Log:
RCS Log:Revision 10.2  2005/12/01 14:38:02  rcatwood
RCS Log:Merged xly_05 changes into the main trunk
RCS Log:Primarily involving melt-back
RCS Log:
RCS Log:Revision 10.1.2.2  2005/11/23 18:19:10  rcatwood
RCS Log:Result of merging mould_source and xly meltback+curvature 2d versions
RCS Log:
RCS Log:Revision 10.1  2005/11/03 11:56:48  rcatwood
RCS Log:New version number -- using mould_src as base
RCS Log:
RCS Log:Revision 1.1.4.2  2005/11/02 11:50:56  rcatwood
RCS Log:Read list files
RCS Log:
RCS Log:Revision 9.2  2003/12/10 17:36:23  kerman
RCS Log:added unit conversions & fixed for counting number of fields from the header
RCS Log:
RCS Log:Revision 1.1.2.1  2003/02/05 14:26:06  rcatwood
RCS Log:Added read_list package to ca cvs tree
RCS Log:
RCS Log:Revision 2.2.4.1  2003/01/15 19:44:40  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 2.2  2002/12/13 17:04:25  rcatwood
RCS Log:Changed to a part of ca source treeh
RCS Log:
RCS Log:Revision 2.1  2002/10/17 16:56:01  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 1.3  2002/09/06 14:55:47  rcatwood
RCS Log:Removed all lint warnings
RCS Log:
RCS Log:Revision 1.2  2002/09/05 18:05:44  rcatwood
RCS Log:Included convert, read and write fg binary, and tested. It seems to work.
RCS Log:
RCS Log:Revision 1.1  2002/09/04 18:40:06  rcatwood
RCS Log:included FGrid_str structure from ca code , added routine to find number of nodes from input file
RCS Log:
RCS Log:
*/
