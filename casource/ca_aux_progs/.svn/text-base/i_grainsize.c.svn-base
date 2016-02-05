#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "../machine.h"
#include "../ca_histo.h"
#include "debin.h"

#ifndef TRUE
#define TRUE 1 
#endif
#ifndef INIT_NUM
#define INIT_NUM 1000
#endif
#ifndef ADD_NUM
#define ADD_NUM 100
#endif


void main(int argc, char * argv[]){
   
	FILE *infile,*outfile,*geofile;
   char *token,*line,*in_fname,*out_fname,*base_fname,*cut_fname,*label;
	int nxny,linelength,nlines,i=0,j=0,k=0;
	int xstart,xend,ystart,yend,cutx,cuty;
	int finput = 0,foutput = 0,flabel=0,cflg,errflg,ctrlflg;
   int xsize=0,ysize=0,nc[3]; 
   int randflag=FALSE,copyflag=FALSE,datflag=FALSE,cutflag=FALSE;
   int maxnum,nremap,ndat=0,idx=0;
   int ** sizelist;
   int * maxnum,*num;
   int nzones=4;

   extern int optind;
   extern char *optarg;
   short int *data, *remap,*datap;
   char * cutdata;


    out_fname = (char *)calloc(255,1);
    cut_fname = (char *)calloc(255,1);
    in_fname = (char *)calloc(255,1);
    base_fname = (char *)calloc(255,1);
    label = (char *)calloc(255,1);
    token = (char *) calloc(MAX_STRING_LEN ,sizeof(char));
    line = (char *) calloc(MAX_STRING_LEN ,sizeof(char));
    errflg=0;
   if (argc == 1) errflg=TRUE;
   fprintf(stderr,"\n");
   #ifdef DEBUG
   fprintf(stderr,"argc = %i\n",argc);
   #endif/*DEBUG*/


   while ((cflg = getopt(argc, argv, "CDai:o:d:hl:r:c:")) != -1) {

      
      switch (cflg) {
      case 'D':
          datflag = TRUE;
          break;
      case 'a':
	      fprintf(stderr,"%s: sorry, flag 'a' is not defined\n",argv[0]);
         break;
      case 'g':
	      fprintf(stderr,"Using control file ca_geoplus.in\n");
	      ctrlflg = 1;
         break;
      case 'i':       /* get input filename */
         finput = TRUE;
         in_fname=strdup(optarg);
         debin (base_fname,in_fname);
          #ifdef DEBUG
          fprintf(stderr,"in_fname: %s\n",in_fname);
          #endif/*DEBUG*/
         break;
      case 'd':
          #ifdef DEBUG
          fprintf(stderr,"d optarg =%s %s \n",optarg,argv[optind]);
          #endif /*DEBUG*/
         if (!(xsize = atoi(optarg))){
            fprintf(stderr,"ERROR with d argument: %s\n",optarg);
            errflg = TRUE; 
         }
         if (!( ysize = atoi(argv[optind]))){
            fprintf(stderr,"ERROR with d argument: %s\n",argv[optind]);
            errflg = TRUE; 
         }
         optind++;
         break;
      case 'c':
         xstart = atoi(optarg);
         ystart = atoi(argv[optind]);
         optind++;
         xend = atoi(argv[optind]);
         optind++;
         yend = atoi(argv[optind]);
         optind++;
         cutflag =TRUE;
         break;
      case 'C': /* copy to another binary file, maybe transformed */
         copyflag = TRUE;
         break;
      case 'o':       /* get output filename */
         foutput = TRUE;
         out_fname = strdup(optarg);
         break;
      case 'l':       /* get label string */
         #ifdef DEBUG
         fprintf(stderr,"l optarg =%s \n",optarg);
         #endif /*DEBUG*/

         flabel = TRUE;
         label = strdup(optarg);
         break;
      case 'r':
         randflag=TRUE;
         nremap = atoi(optarg);
         break;
      case 'h':
         errflg=TRUE;
         break;
      default:
         errflg=TRUE;
         break;
      }
   }


   /* set up rectangle cut */
   if (cutflag){
      if (xstart < 0 || ystart < 0 || xend > xsize || yend > ysize){
         errflg++;
         fprintf(stderr," out of bounds: xstart %i ystart %i xend %i yend %i xsize %i ysize %i \n",xstart,ystart,xend,yend,xsize, ysize);
      }
      cutx = xend - xstart;
      cuty = yend - ystart;
      cutdata = (char *) calloc((cutx*cuty), sizeof(char));
      sprintf(cut_fname,"%s_cut.bin",base_fname);
      fprintf(stderr,"Default cut filename used: %s\n",cut_fname);

   }

   /* handl errors */
   if (errflg ) {
      fprintf(stderr,"i_res: Process I short int (16-bit) grain number output from CA model.\n"\
                     "i_res -i infile [-o outfile] [-r (number) = remap] [-D = make DAT file] [-C = save binary file\n"\
                     "-g Use ca_geoplus.in file for size\n"\
                     "-d xsize ysize\n"\
                     "-c xtart ystart xend yend (cut a rectangle out)\n"\
                      "Default out is [infile].dat\n");
      exit(0);
   }

   if (!foutput){
      sprintf(out_fname,"%s.dat",base_fname);
      fprintf(stderr,"Default output filename used: %s\n",out_fname);
   }
   if (!finput){
      fprintf(stderr,"ERROR: no input file specified! \n");
      exit(0);
   }

    if (ctrlflg){
       if ( (geofile = fopen("ca_geoplus.in","r"))== NULL ){
          fprintf(stderr,"ERROR, could not open ca_geoplus.in\n");
          exit(0);
       }else{
       while (fgets(line,MAX_STRING_LEN,geofile)!= NULL){
            if(line[0]=='%' || line[0]=='#'||(token = strtok(line," ,;\t"))==NULL) {
	         continue;
           }else if (strcasecmp(token,"NCellsPerSB") == 0){
             for (i=0; i<3; i++) {
                if ((token = strtok(NULL, " ,;\t")) != NULL)
                   nc[i]  = atoi(token);
                else {
                   fprintf(stderr,"Error: NCellsPerSB,not found\n");
                   exit(0);
                }
             }
           } 
        }/*end of while loop for getting lines*/
        fclose(geofile);
        xsize = nc[0];
        ysize = nc[1];
     }
     }/* end of read geo file */

    fprintf(stderr,"Processing short int results from %s into %s\n",in_fname,out_fname);


    nxny = xsize * ysize;
    data = (short int *) malloc ( nxny * sizeof(short int));

   /* set up zone size lists now that we know nzones*/
   nz = (xnz * ynz); /* x and y zones */
   sizelist = (int **) calloc((nz),sizeof(int *));
   num = (int *) calloc((nz),sizeof(int));
   maxnum = (int *) calloc((nz),sizeof(int));

   for (i=0;i<nz;i++){
      sizelist[i]=(int *)calloc(INIT_NUM,sizeof(int));
   }
   /*
   zonesize = (int)((float)nxny / (float)nz);
   */


    /* read the data */
    if (!(infile = fopen(in_fname,"r"))){
      fprintf(stderr,"ERROR: could not open  infile %s\n",in_fname);
      exit(0);
    }
    ndat=fread(data,sizeof(short int),nxny,infile);
    fclose(infile);
    fprintf(stdout,"Number of data read: %i\n",ndat);
    if (ndat != nxny){
       fprintf(stderr,"ERROR: i_res: Wrong number of data ! \n");
       exit(0);
    }

    /* process the grain size into the zone histograms */
    for (j=0;j<ynz;j++){
       for(i=0;i<xnz;i++){


    /*remap the data if desired change oriented to random */
    if (randflag){
      fprintf(stderr,"Remapping to %i values\n",nremap);

       maxnum=0;
       for (i=0;i<nxny;i++){
          if ( *(data + i) > maxnum ){
             maxnum = *(data+i);
          }
       }
       remap = (short int *) malloc (maxnum * sizeof(short int));
       *remap = 0;
       for (i=1;i<maxnum;i++){
          (*(remap +i)) = (short int)floor( (( nremap * drand48())+1)); 
       }
       for (i=0;i<nxny;i++){
          data[i] = remap[data[i]];
       }
    }

    
    /*DAT file output for tecplot --just a text file with the header*/
    if (datflag){
       datap=data;
       if (!(outfile = fopen(out_fname,"w"))){
         fprintf(stderr,"ERROR: could not open  outfile %s\n",out_fname);
         exit(0);
       }
       if (flabel) fprintf(outfile,"VARIABLES = \"%s\"\n",label);
       else fprintf (outfile,"VARIABLES = \"CA RESULTS\"\n");

       fprintf (outfile,"ZONE I=%i,J=%i\n",xsize,ysize);

       for (i=0;i<nxny;i++){
          fprintf(outfile,"%i\n",*datap++);
       }
       fclose(outfile);
    }

    /* Binary copy file (possibly remapped) */
    if(copyflag){
        sprintf(out_fname,"C_%s.bin",base_fname);
        fprintf(stdout, "copying transformed data to TWO BYTE binary file ...%s\n",out_fname);
       if (!(outfile = fopen(out_fname,"w"))){
         fprintf(stderr,"ERROR: could not open  outfile %s\n",out_fname);
         exit(0);
       }
       ndat = fwrite(data,sizeof(short int),nxny,outfile);
       fprintf(stdout,"Wrote %i data \n",ndat);
       fclose(outfile);
    }

    /* save the selected rectangle as one-byte data */
    if (cutflag){
       fprintf(stdout, "copying rectangle to ONE BYTE binary file ... %s\n",cut_fname);
       for (j=0;j<cuty;j++){
           for(i=0;i<cutx;i++){
              idx = (j+ystart)*xsize+(i+xstart);
              if (idx >= nxny){
                 fprintf(stderr, "ERROR: i_res: index outside of image ! x= %i y= %i\n",i,j);
                 exit(0);
              }
              if (data[idx]==0) cutdata[j+cutx+i]=(char)(0);
              else cutdata[j*cutx + i] = (char) (data[idx] % 254 + 1);
           }
        }
        if (!(outfile = fopen(cut_fname,"w"))){
          fprintf(stderr,"ERROR: could not open  outfile %s\n",out_fname);
          exit(0);
        }
        fprintf(stderr,"Cut image size is %i x %i\n",cutx,cuty);
        fwrite(cutdata,sizeof(char),(cutx*cuty),outfile);
        fclose(outfile);
    }
    fprintf(stderr,"Total lines processed %i\n",i);
}
