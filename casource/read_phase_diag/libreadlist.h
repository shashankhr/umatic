/*RCS Id:$Id: libreadlist.h 892 2006-03-10 15:24:59Z rcatwood $*/
#ifndef LIBREADLIST_H
#define LIBREADLIST_H

void copy_line (LineData * toline,LineData * fromline);
void convert_line (CA_FLOAT * fdat,LineData * the_line);
void convert_csv(CsvFloat * the_float,CsvData * the_data);

void debin(char * result, char * fname);


extern void free_csv_data (CsvData * cp);
extern void free_float_data(CsvFloat * cp);

void init_csv_data(CsvData * the_data);
void init_line_data(LineData * the_line,int nfields);
void Float_init(CsvFloat * the_float,CsvData * the_data);

LineData * LineCreate(void);

void read_csv_header(FILE * infile,LineData * the_header);
int read_csv_line(FILE * infile, LineData ** the_line,int nfields);
void read_csv(const char * infilename,CsvData * the_data);
void write_csv(const char * outfilename,CsvData * the_data);
void write_float(const char * outfilename,CsvFloat * the_data);
void write_tec(const char * outfilename,CsvFloat * the_float, int nz, int nr);


#ifndef SAFEOPEN_H
#define SAFEOPEN_H
extern FILE * safeopen(const char * fname, const char * type);
#ifdef fopen
#undef fopen
#endif
#define fopen safeopen
#endif

#endif /*LIBREADLIST_H*/

/*
RCS Log:$Log$
RCS Log:Revision 1.1  2006/03/10 15:24:59  rcatwood
RCS Log:Added read_phase_diag for phasediagram lookup table
RCS Log:
RCS Log:Revision 9.2  2003/12/10 17:36:23  kerman
RCS Log:added unit conversions & fixed for counting number of fields from the header
RCS Log:
RCS Log:Revision 1.1.2.3  2003/02/25 19:43:58  kerman
RCS Log:change the unit system from BS to SI as it's necessary for decentred square
RCS Log:
RCS Log:Revision 2.1.6.1  2003/01/15 20:05:41  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 2.1.4.1  2003/01/09 16:26:30  rcatwood
RCS Log:Sorted out several memory prblems.
RCS Log:Modified to handle allvac's data as well as fluent data
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
/*
RCS Log:$Log$
RCS Log:Revision 1.1  2006/03/10 15:24:59  rcatwood
RCS Log:Added read_phase_diag for phasediagram lookup table
RCS Log:
RCS Log:Revision 9.2  2003/12/10 17:36:23  kerman
RCS Log:added unit conversions & fixed for counting number of fields from the header
RCS Log:
RCS Log:Revision 1.1.2.3  2003/02/25 19:43:58  kerman
RCS Log:change the unit system from BS to SI as it's necessary for decentred square
RCS Log:
RCS Log:Revision 2.1.6.1  2003/01/15 20:05:41  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 2.1.4.2  2003/01/09 16:26:30  rcatwood
RCS Log:Sorted out several memory prblems.
RCS Log:Modified to handle allvac's data as well as fluent data
RCS Log:
RCS Log:Revision 2.1.4.1  2003/01/08 17:38:07  rcatwood
RCS Log:merge robert and ahmad read + sort + list programs
RCS Log:
RCS Log:Revision 2.1.2.1  2003/01/08 16:00:34  rcatwood
RCS Log:Ahmad version with sorting and input list of files to process.
RCS Log:
RCS Log:Revision 2.1  2002/10/17 16:56:01  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 1.5  2002/09/06 16:10:55  rcatwood
RCS Log:Changed line structure to double-dereferencet to improve reallocation efficiency
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
RCS Log:Revision 2.1.6.1  2003/01/15 20:05:41  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 2.1.4.1  2003/01/09 16:26:29  rcatwood
RCS Log:Sorted out several memory prblems.
RCS Log:Modified to handle allvac's data as well as fluent data
RCS Log:
RCS Log:Revision 2.1  2002/10/17 16:56:01  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 1.1  2002/09/06 12:52:36  rcatwood
RCS Log:Cured all memory allocation problems .. I hope.
RCS Log:
*/
/*
RCS Log:$Log$
RCS Log:Revision 1.1  2006/03/10 15:24:59  rcatwood
RCS Log:Added read_phase_diag for phasediagram lookup table
RCS Log:
RCS Log:Revision 9.2  2003/12/10 17:36:23  kerman
RCS Log:added unit conversions & fixed for counting number of fields from the header
RCS Log:
RCS Log:Revision 1.1.2.3  2003/02/25 19:43:58  kerman
RCS Log:change the unit system from BS to SI as it's necessary for decentred square
RCS Log:
RCS Log:Revision 2.1.6.1  2003/01/22 16:52:14  rcatwood
RCS Log:Changed the sort order to conform with xuehua and weishan input routines.
RCS Log:
RCS Log:Revision 2.1  2002/10/17 16:56:01  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 1.1  2002/09/04 18:40:54  rcatwood
RCS Log:
RCS Log:header files needed for converting and finding number of nodes, and fgrid structure
RCS Log:
*/
