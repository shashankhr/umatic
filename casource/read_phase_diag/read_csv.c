/* MAIN program file for read_csv project */
/*RCS Id:$Id: read_csv.c 892 2006-03-10 15:24:59Z rcatwood $*/


#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "safeopen.h"
#include "debin.h"

#include "../machine.h"
#include "../constants.h"


#include "csv.h"
#include "../fidap.h"
#include "readwrite.h"
#include "readwrite_fg.h"
#include "read_fg_list.h"
#include "init.h"
#include "convert.h"
#include "findvals.h"
#include "freecsv.h"
#include "../write_fidap_struct.h"
#include "qsort_csv.h"

void print_usage(char * prog_name)      /* print the usage message on error */
{

   fprintf(stderr, "\n\n*************************************************\n");
   fprintf (stderr,"%s \n",prog_name);
   fprintf(stderr, "*\tFilter to read comma-separated file and\n");
   fprintf(stderr, "*\textract the necessary information\n");
   fprintf(stderr, "\n*\tThe following command line options are allowed:\n");
   fprintf(stderr, "*\t-i\t\t-> specify input file\n");
   fprintf(stderr, "*\t-o\t\t-> specify output file\n");
   fprintf(stderr, "*\t-z\t\t-> specify column for z data\n");
   fprintf(stderr, "*\t-r\t\t-> specify column for r data\n");
   fprintf(stderr, "*\t-c\t\t-> specify column for data\n");
   fprintf(stderr, "*\t-t\t\t-> write a tecplot file of original data\n");
   fprintf(stderr, "*\t-V\t\t-> Verbose mode (lots of messages and output files)\n");
   fprintf(stderr, "*\t-h\t\t-> print this message\n");
   fprintf(stderr, "*\t\n");
   fprintf(stderr, "*************************************************\n");
}

void main (int argc, char* argv[]){

   CsvData data;
   CsvFloat floatdata;
   FGrid_str fg,fg_in;
   int errflg=0,cflg;
   int iflag=0,oflag=0,Vflag=0,tflag=0,fflag=0; /* flags for the options */
   int nh=1; /* default number of header lines */
   char * infilename;
   char * outfilename;
   char * ffilename;
   char * basefilename;
   int zfield = 2,rfield=1,tfield=7;
   FILE * fstructp;


#ifdef _DEBUG_MALLOC_INC
  	unsigned long		  oldsize;
	unsigned long		  size1,size2;
   unsigned long  hist0, hist1,hist2;
   oldsize = malloc_inuse( &(hist0));
#endif /*_DEBUG_MALLOC_INC*/

   basefilename = (char *) calloc(255,sizeof(char));
   infilename = (char *) calloc(255,sizeof(char));
   outfilename = (char *) calloc(255,sizeof(char));
   ffilename = (char *) calloc(255,sizeof(char));

   sprintf(infilename,"test.csv");
   sprintf(basefilename,"testout");
   sprintf(ffilename,"f%s.csv",basefilename);

   /*******************************************/
   /* get the command line optinos            */
   /*******************************************/
   errflg=0;
   if (argc < 1) errflg++;
   while ((cflg = getopt(argc, argv, "fai:o:z:r:c:n:htV")) != -1) {
      switch (cflg) {
      case 'a':
	      fprintf(stderr,"%s: sorry, flag 'a' is not defined\n",argv[0]);
         break;
      /* only convert an fg data file to csv  */
      case'f':
          fflag=1;
          break;
      /*get input file name*/
      case 'i':       
         sprintf(infilename,"%s",optarg);
         iflag = 1;
         break;
      /*get output file name*/
      case 'o':       
         sprintf(basefilename,"%s",optarg);
         oflag = 1;
         break;
      /* get the column numbers for processing or use default */
      case 'z':
         zfield = atoi(optarg);
         break;
      case 'r':
         rfield = atoi(optarg);
         break;
      case 'c':
         tfield = atoi(optarg);
         break;
      /* number of header lines */
      case 'n':
         nh = atoi(optarg);
         break;
      /* Tecplot file mode */
      case 't':
         tflag = 1;
         break;
      /* Verbose mode */
      case 'V':
         Vflag = 1;
         break;
      /* print help message */
      case 'h':
      default:
         errflg++;
      break;
      }
   }
   if (errflg) {
      print_usage(argv[0]);
      exit(0);
   }

   /* generate the base file from input file if unspecified */
   if (iflag && !(oflag)){
      debin(basefilename,infilename);
   }


   if (!(fflag)){

      if (Vflag){
         fprintf(stderr,"input filename: %s\n",infilename);
         fprintf(stderr,"output base filename: %s\n",basefilename);
         fprintf(stderr,"z column %i, r column %i, data column %i\n",zfield,rfield,tfield);
      }




      data.nheaders = nh;
      init_csv_data(&data);
      read_csv(infilename,&data);
      convert_csv(&floatdata,&data);
      qsort_csv(&floatdata,rfield,zfield);

      if (Vflag){
         /* write the data as received into <base>_out.csv */
         sprintf(outfilename,"%s_out.csv",basefilename);
         write_csv(outfilename,&data);


         /* write the data after conversion to float into f<base>_out.dat */

         sprintf(outfilename,"f%s_out.csv",basefilename);
         write_float(outfilename,&floatdata);
      }


      /* count the number of nodes, this needs the data to be properly sorted first */
      find_z_r(&floatdata,&fg,zfield,rfield);
      if(tflag){
         /* write the data after conversion to float into f<base>_out.dat */

         sprintf(outfilename,"f%s_out.dat",basefilename);
         write_tec(outfilename,&floatdata,fg.nz,fg.nr);
       }

      /* quick check on number of nodes */
      /* versus data lines actually read */

      if (fg.nnodes != floatdata.line_count){
         fprintf(stderr,"ERROR:read_csv: number of lines != number of nodes, %i, %i\n",fg.nnodes,floatdata.line_count);
         #ifdef ERROR_EXIT
         fprintf(stderr,"Exiting .... \n");
         exit(0);
         #else
         fprintf(stderr,"This will probably cause a segmentation fault! \n");
         fprintf(stderr,"Continuing anyways ... \n");
         #endif /* ERROR_EXIT */
      }
      malloc_nodes(&fg);

      /* put the required data numbers into the fg structure */
      find_temp(&floatdata,&fg,tfield);


      /* free the data structures for csv and float */
      free_csv_data(&data);
      free_float_data(&floatdata);

      if (Vflag){
         /* write the data which has been put in the fg structure (for testing) */
         sprintf(outfilename,"fg_test_%s.csv",basefilename);
         write_fg_csv(outfilename,&fg);
      }

      /**************************************/
      /* write the binary fg structure file */
      /**************************************/
      sprintf(outfilename,"fg_%s.fgb",basefilename);
      write_fg_bin(outfilename,&fg);
      if (Vflag) sprintf(infilename,"%s",outfilename);
   }/* end of (not fflag) */

   if (Vflag||fflag){
      /* read in the file to test , and write again as csv */
      read_fg_bin(infilename,&fg_in,FG_FIRST_READ);

      sprintf(outfilename,"r_%s.csv",basefilename);
      write_fg_csv(outfilename,&fg_in);

      /* write otu the fg structure values as well */
      sprintf(outfilename,"fg_struct_%s.txt",basefilename);
      fstructp=safeopen(outfilename,"w");
      write_fidap_geo_values(fstructp,&fg_in);
      fclose (fstructp);
      free_fg(&fg_in);
   }

   free_fg(&fg);

/* free the string pointers */ 
   free(basefilename);
   free(infilename);
   free(outfilename);
   free(ffilename);

#ifdef _DEBUG_MALLOC_INC
  	size2 = malloc_inuse(&hist2 );
	if( size2 != (oldsize))
	{
		fprintf(stderr,"ERROR: ca_wrapper: dbMalloc test of size of allocated memory\n");
		fprintf(stderr,"\toldsize = %lu, size = %lu - should be %lu\n",
			oldsize, size2, oldsize);
			fprintf(stderr,"Malloc list \n");
			malloc_list(2,hist0,hist2);
			fprintf(stderr,"Finshed dbMalloc size check \n");
	}
	else
	{
		fprintf(stderr,"dbMalloc size check: OK\n");
	}
#endif /*_DEBUG_MALLOC_INC*/

}/* end of MAIN procedure */


/*
RCS Log:$Log$
RCS Log:Revision 1.1  2006/03/10 15:24:59  rcatwood
RCS Log:Added read_phase_diag for phasediagram lookup table
RCS Log:
RCS Log:Revision 9.2  2003/12/10 17:36:23  kerman
RCS Log:added unit conversions & fixed for counting number of fields from the header
RCS Log:
RCS Log:Revision 2.4.4.1  2003/01/15 20:05:41  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 2.4.2.1  2003/01/09 16:26:30  rcatwood
RCS Log:Sorted out several memory prblems.
RCS Log:Modified to handle allvac's data as well as fluent data
RCS Log:
RCS Log:Revision 2.4  2003/01/08 15:56:25  rcatwood
RCS Log:Changes to use Allvac - tecplot output, variable header, alternate seperators
RCS Log:
RCS Log:Revision 2.3  2002/12/16 17:34:10  rcatwood
RCS Log:Improved gnu Makefile Gmakefile, and test of write_fidap_struct.c
RCS Log:
RCS Log:Revision 2.2  2002/12/13 17:04:25  rcatwood
RCS Log:Changed to a part of ca source treeh
RCS Log:
RCS Log:Revision 2.1  2002/10/17 16:56:01  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 1.8  2002/09/12 11:34:19  rcatwood
RCS Log:added profile loop to repeat routine for speed check
RCS Log:
RCS Log:Revision 1.7  2002/09/06 14:55:47  rcatwood
RCS Log:Removed all lint warnings
RCS Log:
RCS Log:Revision 1.6  2002/09/06 12:52:36  rcatwood
RCS Log:Cured all memory allocation problems .. I hope.
RCS Log:
RCS Log:Revision 1.5  2002/09/05 18:05:44  rcatwood
RCS Log:Included convert, read and write fg binary, and tested. It seems to work.
RCS Log:
RCS Log:Revision 1.4  2002/09/04 18:40:06  rcatwood
RCS Log:included FGrid_str structure from ca code , added routine to find number of nodes from input file
RCS Log:
RCS Log:Revision 1.3  2002/09/04 14:58:33  rcatwood
RCS Log:First working version -- reads and writes CSV files (no conversion)
RCS Log:
RCS Log:Revision 1.2  2002/09/04 10:41:07  rcatwood
RCS Log:Added id and log msg
RCS Log:
*/
