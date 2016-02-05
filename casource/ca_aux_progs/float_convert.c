#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

float * float_convert(double * inptr,int ncsb){
    float * outptr;
    int i;
    outptr = (float *) calloc (ncsb,sizeof(float));
      for (i=0;i<ncsb;i++){
         outptr[i]=(float)(inptr[i]);
      }
return(outptr);
}

void print_usage(const char * pname){
   printf("%s <filename> <number> <skip>: Convert double precision raw binary file to float\n",pname);
   printf("Output is placed in fl_<filename>\n");
   printf("<number> is the total number of data points.\n");
   printf("<skip> is the number of header bytes.\n");
}

int main (int argc, char * argv[]){
   char * fname,*ofname,*header;
   FILE * fp;
   float * outptr;
   double * inptr;
   int nread;
   int ndata,skip;

   if (argc != 4){
      print_usage(argv[0]);
      exit(0);
   }
   fname = strdup(argv[1]);
   ndata = atoi(argv[2]);
   skip = atoi(argv[3]);

  inptr=(double *)calloc (ndata,sizeof(double));
  header=(char *) calloc(skip,sizeof(char));

   fp = fopen(fname,"r");
   fread(header,sizeof(char),skip,fp);
   nread=fread(inptr,sizeof(double),ndata,fp);
   fclose(fp);

   if (nread != ndata){
      fprintf(stderr,"ERROR: %s: Wrong number %i of values read;\n possibly the input file is not double precision,\nor the wrong number of values were requested %i.\n\n",argv[0],nread,ndata);
      if (feof (fp)) fprintf(stderr,"End-of-file encountered.\n");
      if (ferror (fp)){
          fprintf(stderr,"File error encountered.\n");
          fprintf(stderr,"%i  %s\n",errno,strerror(errno));
      }
      
      print_usage(argv[0]);
      exit(0);
   }

   outptr = float_convert(inptr,ndata);

   ofname = (char *) calloc(255,sizeof(char));

   snprintf(ofname,255,"fl_%s",fname);
   fp = fopen(ofname,"w");
   fwrite(outptr,sizeof(float),ndata,fp);
   fclose(fp);
   free(outptr);
   free(inptr);
   free(fname);
   free(ofname);
   return(0);
}








   
