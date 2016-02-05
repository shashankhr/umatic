#include <stdio.h>
#include <string.h>
/* remove the file extension from a string representing a filename */
/* assuming only a single . in the filename ... */

void debin(char * result, char * fname){



   sscanf(fname,"%[^.].",result);
}

void deprefix(char * result, char * fname){
   sscanf(fname,"F_[GA]_%s",result);
}

#ifdef TEST_MAIN
void main(){
   char name[256];
   char res[256];

   sprintf(name,"bl.ah.txt");
   printf ("name = %s\n",name);
   debin (res,name);
   printf ("res = %s\n",res);

}
#endif

