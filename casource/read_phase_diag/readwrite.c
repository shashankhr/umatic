/*RCS Id:$Id: readwrite.c 892 2006-03-10 15:24:59Z rcatwood $*/
#include <stdio.h>
#include <string.h>
#include "../safeopen.h"
#include "../machine.h"
#include "../constants.h"
#include "csv.h"
#include "init.h"

/* define the seperator characters */
/* to be used by all strtok parsing commands */
const char * sep = ", \n\t\r";

/* read a header line from the CSV file */
/* and figure out how many fields per line */
void read_csv_header(FILE * infile,LineData * the_header){
	char * token;
	char * line;
	int nwords=0;
	line = (char *)malloc(LINELENGTH * sizeof(char));

	fgets(line,LINELENGTH,infile);
	token = strtok(line,sep);
	if (token == NULL) return;
	the_header->fields[nwords] = strdup(token);
	nwords ++;
	while ((token = strtok(NULL,sep)) != NULL){
         the_header->fields[nwords] = strdup(token);
			nwords++;
			/* TODO: there is a problem when the header is not strictly seperated */
			/* by the correct characters eg " time (s) vel (m/s) " is read as 4 fields */
			if(nwords >= the_header->nperline){
				the_header->nperline += ADD_FIELDS;
				the_header->fields = (char **) realloc(the_header->fields,the_header->nperline*sizeof(char *));
			}
	}
	the_header->nfields = nwords;
	free(line);
	free(token);
    
}

/* read a line from the csv file, */
/* knowing already how many fields there are */
int read_csv_line(FILE * infile, LineData ** the_line,int nfields){
	char * token;
	char * line;
	int nwords=0;
	line = (char *)malloc((LINELENGTH+1) * sizeof(char));

	/*try to read a line, if the end of file is encountered, */
	/*then we are done! */
	if ( fgets(line,LINELENGTH,infile) == NULL ){
	  free(line);
	  return(1);
	}

	/* if the line does not parse, we are also done! */
	token = strtok(line,sep);
	if (token == NULL){
      free(token);
      free(line);
      return(1);
	}

   *the_line = LineCreate();
   init_line_data (*the_line,nfields);
	(*the_line)->fields[nwords] = strdup(token);
	nwords ++;
	while ((token = strtok(NULL,sep)) != NULL){
         (*the_line)->fields[nwords] = strdup(token);
		nwords++;
			if(nwords > (*the_line)->nfields){
            fprintf(stderr,"ERROR: read_csv_line: Number of fields per line exceeded!\n");
			}
	}
	free(line);
	free(token);
	return(0);
}

/* read in the csv file, first find out how many fields per line, */
/* then figure how many lines, dynamically allocating the space   */
void read_csv(const char * infilename,CsvData * the_data){
   int i;
   int done=0;
   int nfields;
   FILE * infile;
   infile = safeopen(infilename,"r");
   the_data->line_count = 0;

   /* read the header lines */
   for(i=0;i<the_data->nheaders;i++){
      the_data->headers[i] = LineCreate();
      init_line_data(the_data->headers[i],N_FIELDS);
      read_csv_header(infile,(the_data->headers[i]));
   }
   /* set the number of fields from the last header line*/
   nfields = the_data->headers[the_data->nheaders-1]->nfields;
/*   nfields = the_data->headers[0]->nfields; */
   /* read the data lines */
   while (done==0){
      done = read_csv_line(infile,(the_data->lines + the_data->line_count),nfields);
      (the_data->line_count)++;
      if (the_data->line_count >= the_data->nlines){
         #ifdef VERBOSE
         fprintf(stderr,"Read %i Lines...\n",the_data->line_count);
         #endif /*VERBOSE*/
         the_data->nlines += ADD_DATA;
         the_data->lines = (LineData **)realloc(the_data->lines,the_data->nlines * sizeof(LineData *));
      }
   }
   the_data->line_count--;
   fprintf(stderr,"Done reading %i Lines.\n",the_data->line_count);
   fclose(infile);
}

/***************************************/
/* routines to write out the csv file  */
/* from the string data structure      */
/***************************************/

/* write a line */
void write_csv_line(FILE * ofile,LineData * the_line){
   int i;
   int nm=the_line->nfields-1;

   for (i=0;i<nm;i++){
      fprintf(ofile,"%s,",the_line->fields[i]);
   }
   fprintf(ofile,"%s",the_line->fields[i]);
   fprintf(ofile,"\n");
}

/*write a line of float data*/
void write_float_line(FILE * ofile,CA_FLOAT * the_line,int nm){
   int i;

   for (i=0;i<nm-1;i++){
      fprintf(ofile,"%.10g,",the_line[i]);
   }
   fprintf(ofile,"%.10g",the_line[i]);
   fprintf(ofile,"\n");
}

/* write the file */
void write_csv(const char * outfilename,CsvData * the_data){
   int i;
   FILE * ofile;

   ofile = safeopen(outfilename,"w");
   for (i=0;i<the_data->nheaders;i++){
      #ifdef NUMLINE
      fprintf(ofile,"%i,",i);
      #endif
      write_csv_line(ofile,(the_data->headers[i]));
   }

   for (i=0;i<the_data->line_count;i++){
      #ifdef NUMLINE
      fprintf(ofile,"%i,",i);
      #endif
       write_csv_line(ofile,(the_data->lines[i]));
   }
   fclose(ofile);
}

/************************/
/* write the CA_FLOAT data */
/* to a csv file        */
/************************/

void write_float(const char * outfilename,CsvFloat * the_data){
   int i;
   FILE * ofile;
   ofile = safeopen(outfilename,"w");
   for (i=0;i<the_data->nheaders;i++){
      #ifdef NUMLINE
      fprintf(ofile,"%i,",i);
      #endif
      write_csv_line(ofile,(the_data->headers[i]));
   }

   for (i=0;i<the_data->line_count;i++){
      #ifdef NUMLINE
      fprintf(ofile,"%i,",i);
      #endif
       write_float_line(ofile,(the_data->data[i]),the_data->nfields);
   }
   fclose(ofile);
}

/*write the CA_FLOAT data to a tecplot DAT file */
void write_tec(const char * outfilename,CsvFloat * the_float, int nz, int nr){
   int i;
   FILE * ofile;
   LineData * the_head;
   int nm;
   the_head = the_float->headers[the_float->nheaders-2];
   nm=the_head->nfields;

   ofile = safeopen(outfilename,"w");
   fprintf(ofile,"variables=\"");

   for (i=0;i<nm;i++){
      fprintf(ofile,"%s\",\"",the_head->fields[i]);
   }
/*   fprintf(ofile,"%s\"",the_head->fields[i]); */
   fprintf(ofile,"\n");
   fprintf(ofile,"zone i=%i j=%i f=point\n",nr,nz);

   for (i=0;i<the_float->line_count;i++){
       write_float_line(ofile,(the_float->data[i]), nm);
   }
   fclose(ofile);
}


/*
RCS Log:$Log$
RCS Log:Revision 1.1  2006/03/10 15:24:59  rcatwood
RCS Log:Added read_phase_diag for phasediagram lookup table
RCS Log:
RCS Log:Revision 1.1.2.1  2004/07/28 10:39:06  rcatwood
RCS Log:Added copies of read-list files suitably modified
RCS Log:
RCS Log:Revision 9.3  2003/12/10 17:36:23  kerman
RCS Log:added unit conversions & fixed for counting number of fields from the header
RCS Log:
RCS Log:Revision 2.4.4.4  2003/05/16 16:13:41  kerman
RCS Log:Change the unit conversion routine and the location where the number of fields
RCS Log:are being read from the data header files
RCS Log:
RCS Log:Revision 2.4.4.3  2003/01/23 17:47:33  rcatwood
RCS Log:finite grid applied to decentered square,
RCS Log:works, but not checked for correct results.
RCS Log:
RCS Log:Revision 2.4.4.2  2003/01/22 16:52:15  rcatwood
RCS Log:Changed the sort order to conform with xuehua and weishan input routines.
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
RCS Log:Revision 2.3  2002/12/13 18:09:18  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 2.2  2002/12/13 17:04:25  rcatwood
RCS Log:Changed to a part of ca source treeh
RCS Log:
RCS Log:Revision 2.1  2002/10/17 16:56:01  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 1.5  2002/09/06 16:10:55  rcatwood
RCS Log:Changed line structure to double-dereferencet to improve reallocation efficiency
RCS Log:
RCS Log:Revision 1.4  2002/09/06 14:55:47  rcatwood
RCS Log:Removed all lint warnings
RCS Log:
RCS Log:Revision 1.3  2002/09/06 12:52:36  rcatwood
RCS Log:Cured all memory allocation problems .. I hope.
RCS Log:
RCS Log:Revision 1.2  2002/09/05 18:05:44  rcatwood
RCS Log:Included convert, read and write fg binary, and tested. It seems to work.
RCS Log:
RCS Log:Revision 1.1  2002/09/04 14:58:33  rcatwood
RCS Log:First working version -- reads and writes CSV files (no conversion)
RCS Log:
*/
