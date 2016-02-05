/*RCS Id:$Id: freecsv.h 892 2006-03-10 15:24:59Z rcatwood $*/

extern void free_csv_data (CsvData * cp);
extern void free_float_data(CsvFloat * cp);
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
