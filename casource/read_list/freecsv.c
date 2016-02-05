/*RCS Id:$Id: freecsv.c 887 2006-03-01 18:21:01Z rcatwood $*/
#include <stdio.h>
#include <string.h>
#include "../machine.h"
#include "../constants.h"
#include "csv.h"
#include "../fidap.h"
#include "read_fg_list.h"

void free_line(LineData * lp){
   int i;
   #ifdef PC
   fprintf(stderr,"Freeing %i \n",lp->ct);
   #endif
   if (lp->nfields > 0){
      for (i=0;i<lp->nfields;i++){
         free(lp->fields[i]);
   #ifdef PC
         fprintf(stderr,"%i... ",i);
   #endif
      }
   }
   #ifdef PC
   fprintf(stderr,"\n");
   #endif
   lp->nfields = 0;
   free (lp->fields);
   free(lp);
}
void free_fg_row( Fg_row * row){

   free (row->filename);
   free (row);

}

void free_csv_data (CsvData * cp){
   int i;
   for (i=0;i<cp->line_count;i++){
      free_line((cp->lines[i]));
   }
   for (i=0;i<cp->nheaders;i++){
      free_line((cp->headers[i]));
   }
   free(cp->lines);
   free(cp->headers);
   cp->nheaders=cp->nlines=cp->line_count=0;
}

void free_float_data(CsvFloat * cp){
   int i;
   for (i=0;i<cp->nheaders;i++){
      free_line((cp->headers[i]));
   }
   for (i=0;i<cp->line_count;i++){
      free(cp->data[i]);
   }
   free (cp->data);
   free (cp->headers);
   cp->nheaders = 0;
   cp->nlines=0;
   cp->line_count=0;
   cp->nfields=0;
}
void free_fg(FGrid_str *fg){
   free(fg->r);
   free(fg->z);
   free(fg->Fidap_T);
   fg->nr = 0;
   fg->nz = 0;
   fg->nnodes = 0;
}

void free_fg_list(FGrid_list_str *fgl){
   int i;
   for (i=0;i<fgl->nheaders_list;i++){
      free_line((fgl->headers[i]));
   }
   free(fgl->headers);
   fgl->nheaders_list = 0;
   for (i=0;i<fgl->nrows;i++){
      free_fg_row((fgl->rows[i]));
   }
   free(fgl->rows);
   fgl->nrows=0;
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
RCS Log:Revision 2.2.4.2  2003/01/15 20:05:41  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 2.2.2.1  2003/01/09 16:26:29  rcatwood
RCS Log:Sorted out several memory prblems.
RCS Log:Modified to handle allvac's data as well as fluent data
RCS Log:
RCS Log:Revision 2.2  2002/12/13 17:04:25  rcatwood
RCS Log:Changed to a part of ca source treeh
RCS Log:
RCS Log:Revision 2.1  2002/10/17 16:56:01  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 1.2  2002/09/06 16:10:55  rcatwood
RCS Log:Changed line structure to double-dereferencet to improve reallocation efficiency
RCS Log:
RCS Log:Revision 1.1  2002/09/06 12:52:36  rcatwood
RCS Log:Cured all memory allocation problems .. I hope.
RCS Log:
*/
