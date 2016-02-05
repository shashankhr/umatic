/*RCS Id:$Id: read_fg_list.c 887 2006-03-01 18:21:01Z rcatwood $*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include "safeopen.h"
#include "machine.h"
#include "constants.h"
#include "fidap.h"
#include "csv.h"
#include "read_fg_list.h"
#include "freecsv.h"
#include "readwrite.h"
#include "readwrite_fg.h"
#include "init.h"

    
/***********************************************************************/
/* subroutine to read the control list for input fg temperature files. */
/*   Inputs: file name where the list is stored */
/*   Output: The nubmer of entries in the list  */
/*   Result: Stores the necessary information in the fg_list structure fgl */
/* THis should include the actual filename and the appropriate time for that data */
/* and may include other information as necessary */
/***********************************************************************/
extern  void copy_line (LineData * toline,LineData * fromline);

void fill_list (FGrid_list_str * fgl, CsvData * list_ptr){
   int i,j;
   /* keep track of the number of lines and headers */
   fgl->nrows = list_ptr->line_count;
   fgl->nheaders_list=list_ptr->nheaders;

   /* Put the data into the correct fields of the list structure */
   fgl->rows = (Fg_row **) calloc(fgl->nrows,sizeof(Fg_row *));
   for (i=0;i<list_ptr->line_count;i++){
      fgl->rows[i] = (Fg_row *) calloc (1,sizeof(Fg_row));
      fgl->rows[i]->filename = strdup(list_ptr->lines[i]->fields[0]);
      fgl->rows[i]->time = atof(list_ptr->lines[i]->fields[1]);
      fgl->rows[i]->z_offset = atof(list_ptr->lines[i]->fields[2]);
   }

   /* copy the header information */
   fgl->headers= (LineData **) calloc(fgl->nheaders_list,sizeof(LineData *));
   for (i=0;i<list_ptr->nheaders;i++){
        fgl->headers[i] = (LineData *)calloc(1,sizeof(LineData));
        fgl->headers[i]->nfields = list_ptr->headers[i]->nfields;
        init_line_data(fgl->headers[i],list_ptr->headers[i]->nfields);
       copy_line(fgl->headers[i],list_ptr->headers[i]);
   }

}/* end of fill_list */

int read_listfile(const char * listfilename, FGrid_list_str * fgl){
   int n_names=0;
   CsvData list_data;
   int nfields_list;
   int i;
   int done=0;
    
   list_data.nheaders=fgl->nheaders_list;
   init_csv_data(&list_data);
   read_csv(listfilename,&list_data);
    /* now we know how many rows */
   fill_list(fgl,&list_data);

   fprintf(stderr,"Done reading %i Rows from the list.\n",fgl->nrows);
   n_names=fgl->nrows;
   free_csv_data(&list_data);
   return(n_names);
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
RCS Log:Revision 1.1.2.3  2003/02/25 19:43:58  kerman
RCS Log:change the unit system from BS to SI as it's necessary for decentred square
RCS Log:
RCS Log:Revision 1.1.2.1  2003/01/16 19:26:22  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 1.1.6.1  2003/01/15 20:05:41  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 1.1.4.2  2003/01/09 16:26:30  rcatwood
RCS Log:Sorted out several memory prblems.
RCS Log:Modified to handle allvac's data as well as fluent data
RCS Log:
RCS Log:Revision 1.1.4.1  2003/01/08 17:38:57  rcatwood
RCS Log:adding new files neede for readin gthe list of files
RCS Log:
RCS Log:Revision 1.1.2.1  2003/01/08 16:50:15  rcatwood
RCS Log:added files for list processing
RCS Log:
RCS Log:Revision 8.2  2002/12/13 14:04:43  rcatwood
RCS Log:Compile without errors
RCS Log:
RCS Log:Revision 8.1  2002/12/13 13:42:25  rcatwood
RCS Log:Read and organize the data files for inputting
RCS Log:finite-element data
RCS Log:
RCS Log:Revision 2.1  2002/10/17 16:56:01  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 1.2  2002/09/06 13:06:42  rcatwood
RCS Log:improved for compatiblity with CA code
RCS Log:
RCS Log:Revision 1.1  2002/09/06 12:52:36  rcatwood
RCS Log:Cured all memory allocation problems .. I hope.
RCS Log:
RCS Log:Revision 1.2  2002/09/05 18:05:44  rcatwood
RCS Log:Included convert, read and write fg binary, and tested. It seems to work.
RCS Log:
RCS Log:Revision 1.1  2002/09/04 14:58:33  rcatwood
RCS Log:First working version -- reads and writes CSV files (no conversion)
RCS Log:
*/
