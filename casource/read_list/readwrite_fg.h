/*RCS Id:$Id: readwrite_fg.h 1007 2007-04-27 18:57:20Z  $*/
#ifndef READWRITE_FG_H

#define READWRITE_FG_H 1

#define FG_FIRST_READ 0
#define FG_TRANS_READ 1
#define FG_CLEANUP 2

#include "read_fg_list.h"
void write_fg_csv(const char * outfilename,FGrid_str *fg);
void read_fg_bin(const char * infilename,FGrid_str *fg,int fg_flag);
void write_fg_bin(const char * outfilename,FGrid_str *fg);
int read_listfile(const char * listfilename, FGrid_list_str * fgl);
#endif /*READWRITE_FG_H*/
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
RCS Log:Revision 1.1.2.4  2003/02/27 23:04:50  rcatwood
RCS Log:Removed use of old temperature routines , all temperatures shoudl
RCS Log:be determined by checking the array c_temp in teh subblock, if the
RCS Log:subblock is open
RCS Log:
RCS Log:Revision 1.1.2.3  2003/02/25 19:43:58  kerman
RCS Log:change the unit system from BS to SI as it's necessary for decentred square
RCS Log:
RCS Log:Revision 2.2  2002/12/13 17:04:25  rcatwood
RCS Log:Changed to a part of ca source treeh
RCS Log:
RCS Log:Revision 2.1  2002/10/17 16:56:01  rcatwood
RCS Log:*** empty log message ***
RCS Log:
RCS Log:Revision 1.2  2002/09/06 14:55:47  rcatwood
RCS Log:Removed all lint warnings
RCS Log:
RCS Log:Revision 1.1  2002/09/06 12:52:36  rcatwood
RCS Log:Cured all memory allocation problems .. I hope.
RCS Log:
*/
