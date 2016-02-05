/* read_fg_list.h */

#include "machine.h"
#include "constants.h"
#include "csv.h"

typedef struct fg_list_row{
   char * filename;
   CA_FLOAT time;
   CA_FLOAT z_offset;   
}Fg_row;

typedef struct fg_list_str{
   int nheaders_list;
   int nrows;
   LineData ** headers;
   Fg_row ** rows;
}FGrid_list_str;
