/*RCS Id:$Id: unit_conv.c 887 2006-03-01 18:21:01Z rcatwood $*/
/* to change units from BS to cgs */
#include <stdio.h>
#include <string.h>
#include "csv.h"

void unit_conv_line(CA_FLOAT * the_line,int r,int z,int t, int flag){

    switch (flag){
	case 1:
      	  the_line[r]*=0.01; /* radius from cm to m */
     	  the_line[z]*=0.01; /* height from cm to m */
          break;
 	case 2: /* CGS to SI */
      	  the_line[r]*=0.01; /* radius from cm to m */
     	  the_line[z]*=0.01; /* height from cm to m */
          the_line[t]=(5./9.)*(the_line[t]-32.); /* temp. from deg.F to C */
          break;
	case 3: /* BS to SI */
      	  the_line[r]*=0.0254; /* radius from in to m */
     	  the_line[z]*=0.0254; /* height from in to m */
          the_line[t]=(5./9.)*(the_line[t]-32.); /* temp. from deg.F to C */
          break; 
    }
}

void unit_conv(CsvFloat * the_float,int rfield,int zfield,int tfield, int flag){

   int i;

   for (i=0;i<the_float->line_count;i++){
      unit_conv_line((the_float->data[i]),rfield,zfield,tfield,flag);
   }

} /* end of unit conversion */
