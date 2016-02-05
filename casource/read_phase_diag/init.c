/*RCS Id:$Id: init.c 892 2006-03-10 15:24:59Z rcatwood $*/
#include <stdio.h>
#include "../machine.h"
#include "../constants.h"
#include "csv.h"
#include "convert.h"
#include "../fidap.h"
LineData * LineCreate(void){
   return (LineData *) calloc(1,sizeof(LineData));
}
void init_csv_data(CsvData * the_data){
   the_data->nlines = N_DATA;
   the_data->line_count = 0;
   the_data->headers = (LineData **)calloc(the_data->nheaders,sizeof(LineData *));
   the_data->lines = (LineData **)calloc(the_data->nlines,sizeof(LineData *));
}
void init_line_data(LineData * the_line,int n){
   #ifdef PC
   static int count = 0;
   the_line->ct = count;
   if(count == 109){
      int dumb;
      dumb = 0;
      dumb = 1;
   }
   fprintf(stderr,"callocing %i\n",count++);
   #endif
   the_line->nfields = n;
   the_line->fields = (char **)calloc(the_line->nfields,sizeof(char *));
}
   
void Float_init(CsvFloat * the_float,CsvData * the_data){
     int i;
     the_float->nheaders = the_data->nheaders;
     the_float->nlines = the_data->nlines;
     the_float->line_count = the_data->line_count;
     the_float->nfields = the_data->headers[0]->nfields;
     the_float->headers = (LineData **)calloc(the_float->nheaders,sizeof(LineData *));

     for (i=0;i<the_data->nheaders;i++){
        the_float->headers[i] = (LineData *)calloc(1,sizeof(LineData));
        the_float->headers[i]->nfields = the_data->headers[i]->nfields;
        init_line_data(the_float->headers[i],the_float->headers[i]->nfields);
     }

     the_float->data = (CA_FLOAT **) calloc(the_data->line_count,sizeof(CA_FLOAT *));

     for (i=0;i<the_data->line_count;i++){
        (*(the_float->data + i)) = (CA_FLOAT *) calloc(the_data->lines[i]->nfields,sizeof(CA_FLOAT));
     }

}

void malloc_nodes(FGrid_str *fg){
     fg->Fidap_T = (CA_FLOAT *)malloc(fg->nnodes * sizeof(CA_FLOAT));
     if (fg->Fidap_T == NULL){
        fprintf(stderr,"Malloc failed: fg->Fidap_T");
        exit(0);
     }
}

/*
RCS Log:$Log$
RCS Log:Revision 1.1  2006/03/10 15:24:59  rcatwood
RCS Log:Added read_phase_diag for phasediagram lookup table
RCS Log:
RCS Log:Revision 1.1.2.1  2004/07/28 10:39:06  rcatwood
RCS Log:Added copies of read-list files suitably modified
RCS Log:
RCS Log:Revision 9.2  2003/12/10 17:36:23  kerman
RCS Log:added unit conversions & fixed for counting number of fields from the header
RCS Log:
RCS Log:Revision 1.1.2.3  2003/02/25 19:43:58  kerman
RCS Log:change the unit system from BS to SI as it's necessary for decentred square
RCS Log:
RCS Log:Revision 2.2.4.2  2003/01/22 16:52:14  rcatwood
RCS Log:Changed the sort order to conform with xuehua and weishan input routines.
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
RCS Log:Revision 1.6  2002/09/06 16:10:55  rcatwood
RCS Log:Changed line structure to double-dereferencet to improve reallocation efficiency
RCS Log:
RCS Log:Revision 1.5  2002/09/06 14:55:47  rcatwood
RCS Log:Removed all lint warnings
RCS Log:
RCS Log:Revision 1.4  2002/09/06 12:52:36  rcatwood
RCS Log:Cured all memory allocation problems .. I hope.
RCS Log:
RCS Log:Revision 1.3  2002/09/05 18:05:44  rcatwood
RCS Log:Included convert, read and write fg binary, and tested. It seems to work.
RCS Log:
RCS Log:Revision 1.2  2002/09/04 18:40:06  rcatwood
RCS Log:included FGrid_str structure from ca code , added routine to find number of nodes from input file
RCS Log:
RCS Log:Revision 1.1  2002/09/04 14:58:33  rcatwood
RCS Log:First working version -- reads and writes CSV files (no conversion)
RCS Log:
*/
