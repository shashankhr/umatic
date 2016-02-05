#include <stdio.h>
#include <string.h>
/* remove the file extension from a string representing a filename */
/* assuming only a single . in the filename ... */

void debin(char * result, char * fname){



   sscanf(fname,"%[^.].",result);
}

#ifdef PIPE
void main(){
   char name[256];
   char res[256];

   gets(name);
   debin (res,name);
   printf ("%s\n",res);

}
#endif

