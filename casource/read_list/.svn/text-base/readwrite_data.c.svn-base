/* readwrite_data.c */

#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include "machine.h"
#include "constants.h"

#include "csv.h"
#include "read_fg_list.h"
#include "fidap.h"
#include "readwrite.h"
#include "readwrite_fg.h"
#include "init.h"
#include "convert.h"
#include "findvals.h"
#include "freecsv.h"
#include "qsort_csv.h"
#include "debin.h"

extern void unit_conv(CsvFloat * the_float,int rfield,int zfield,int tfield, int flag);
void readwrite_data(const Flags flags, Fg_row * fgr,  int rfield, int zfield, int tfield){

   char * vfilename;
   char * outfilename;
   char * infilename;
   char * basefilename;
   CA_FLOAT maxht =0;
   CsvData * the_data;
   CsvFloat * the_float;
   FGrid_str * fg;
   int i;


   if (flags.V) fprintf(stderr,"readwrite_data: reading %s\n",fgr->filename);
   the_data = (CsvData *)malloc(sizeof(CsvData));
   the_float = (CsvFloat *)malloc(sizeof(CsvFloat));
   the_data->nheaders = flags.nh;

   fg = (FGrid_str *) calloc(1,sizeof(FGrid_str));
   outfilename = (char *) calloc(255,sizeof(char));
   infilename = (char *) calloc(255,sizeof(char));
   basefilename = fgr->filename;
   sprintf(infilename,"%s%s",basefilename,flags.ext);
if (flags.V) fprintf(stderr,"readwrite_data: infile %s\n",infilename); 
   init_csv_data(the_data);
   read_csv(infilename, the_data);

   convert_csv(the_float, the_data);
   if ( rfield > the_float->nfields ||
       zfield > the_float->nfields ||
       tfield > the_float->nfields){

      fprintf(stderr,"ERROR:readwrite_data:field overflow, n %i z %i r %i t %i\n",the_float->nfields,zfield,rfield,tfield);
      exit(0);
   }

   if(flags.u){ /* to change the units */
      unit_conv(the_float, rfield, zfield, tfield, flags.u); 
   }

   /* reverse the height as in SMPC data */
   if(flags.f){
		for (i=0;i<the_float->line_count;i++){
			if (the_float->data[i][zfield] >= maxht){
			 maxht = the_float->data[i][zfield];
			}
		}
		for (i=0;i<the_float->line_count;i++){
			the_float->data[i][zfield]  = maxht - the_float->data[i][zfield];
		}
   }


   qsort_csv(the_float, zfield, rfield); 

   if (flags.V){ /* verbose output */
      vfilename = (char *) calloc(255,sizeof(char));
      /* write the data as received into <base>_out.csv */
      sprintf(vfilename,"AR_%s.csv",basefilename);
      write_csv(vfilename,the_data);

      /* write the data after conversion to float into f<base>_out.dat */
      sprintf(vfilename,"F_%s.csv",basefilename);
      write_float(vfilename,the_float);
      free (vfilename);
   }

   find_z_r(the_float, fg, zfield, rfield);
   /* quick check on number of nodes */
   /* versus data lines actually read */
   if (fg->nnodes != the_float->line_count){
      fprintf(stderr,"ERROR:read_csv: number of lines != number of nodes, %i, %i\n",fg->nnodes,the_float->line_count);
      #ifdef ERROR_EXIT 
      exit(0);
      #else
      fprintf(stderr,"Continuing anyways ... \n");
      #endif /* ERROR_EXIT */
   }
   malloc_nodes(fg);

   /* transfer the data to the fg structure */
   find_temp(the_float,fg,tfield);
   fg->time_var = fgr->time;
   fg->z_off = fgr->z_offset;
   fg->h_ingot = fg->z[0];

   /* flip the z values as in weishan's code */
	for (i=0;i<fg->nz;i++){
      		fg->z[i] = fg->h_ingot - fg->z[i];
   } 

   sprintf(outfilename,"%s_out.csv",basefilename);
   write_fg_csv(outfilename,fg);

   /* write a tecplot file of the float-converted (sorted) data */
   if (flags.t){
      sprintf(outfilename,"%s.dat",basefilename);
      write_tec(outfilename,the_float,fg->nz,fg->nr);
   }

   sprintf(outfilename,"%s.fgb",basefilename);
   write_fg_bin(outfilename,fg);
   
   free_csv_data(the_data);
   free_float_data(the_float);
   free(the_data);
   free(the_float);
   free_fg(fg);
   free(fg);
   free(outfilename);
   free(infilename);

} /* end of readwrite_data.c file */
