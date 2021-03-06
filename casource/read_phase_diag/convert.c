/*RCS Id:$Id: convert.c 892 2006-03-10 15:24:59Z rcatwood $*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../machine.h"
#include "../constants.h"
#include "csv.h"
#include "init.h"

void copy_line (LineData * toline,LineData * fromline){
     int i;
     toline->nperline = fromline->nperline;
     toline->nfields = fromline->nfields;
     #ifdef PC
     toline->ct = fromline->ct;
     #endif
     for (i=0;i<fromline->nfields;i++){
        toline->fields[i] = strdup(fromline->fields[i]);
     }
}



void convert_line (CA_FLOAT * fdat,LineData * the_line){
   int i;
   for (i=0;i<the_line->nfields;i++){
      if(the_line->fields[i] == NULL) *fdat++ = 0;
      else *fdat++ = (CA_FLOAT)atof(the_line->fields[i]);
   }
}

void convert_csv(CsvFloat * the_float,CsvData * the_data){
   int i;

   Float_init(the_float,the_data);
  for (i=0;i<the_data->nheaders;i++){
     copy_line((the_float->headers[i]),(the_data->headers[i]));
  }

   for (i=0;i<the_data->line_count;i++){
      convert_line(*(the_float->data + i),(the_data->lines[i]));
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
RCS Log:Revision 2.3.4.2  2003/01/22 16:52:14  rcatwood
RCS Log:Changed the sort order to conform with xuehua and weishan input routines.
RCS Log:
RCS Log:Revision 2.3.4.1  2003/01/15 20:05:41  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 2.3.2.1  2003/01/08 17:38:06  rcatwood
RCS Log:merge robert and ahmad read + sort + list programs
RCS Log:
RCS Log:Revision 2.3  2003/01/08 15:56:25  rcatwood
RCS Log:Changes to use Allvac - tecplot output, variable header, alternate seperators
RCS Log:
RCS Log:Revision 2.2  2002/12/13 17:04:25  rcatwood
RCS Log:Changed to a part of ca source treeh
RCS Log:
RCS Log:Revision 2.1  2002/10/17 16:56:01  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 1.4  2002/09/06 16:10:55  rcatwood
RCS Log:Changed line structure to double-dereferencet to improve reallocation efficiency
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
