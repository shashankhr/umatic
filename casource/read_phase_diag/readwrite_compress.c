#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>

/* compress one array into a compressed buffer and write to the file */
/* return the size in bytes of the compressed data written */
size_t write_comp_array(void * datap,size_t vsize,int nmemb,FILE * fp){
    Bytef * comp_data;
    uLongf comp_size,datasize;
    datasize =  (uLongf)(nmemb*vsize);
    /* minimum compressed buffer size */
    comp_size = (uLongf)(ceil(1.001 * (double)(datasize)) + 13);
    comp_data = (Bytef *) calloc(comp_size,sizeof(Bytef));
    compress(comp_data,&comp_size,(Bytef * )datap,datasize);
    fwrite(&comp_size,sizeof(uLongf),1,fp);
    fwrite(comp_data,sizeof(Bytef),comp_size,fp);
    free(comp_data);
    return ( (size_t) comp_size * sizeof(Bytef) + sizeof(uLongf));
}

/* read teh comressed data, return the size of the uncompressed data */
int read_comp_array(void * datap,size_t vsize,int nmemb,FILE * fp){
    Bytef * comp_data;
    uLongf comp_size,datasize;
    size_t nread=0;
    datasize =  (uLongf)(nmemb*vsize);
    /* minimum compressed buffer size */

    nread = fread(&comp_size,sizeof(uLongf),1,fp);
    if (nread != 1){
       fprintf(stderr,"ERROR:read_comp_array: Read error, expected 1 got %li\n",nread);
    }
    comp_data = (Bytef *) calloc(comp_size,sizeof(Bytef));
    nread = fread(comp_data,sizeof(Bytef),comp_size,fp);
    if (nread != comp_size){
       fprintf(stderr,"ERROR:read_comp_array: Read error, expected %li got %li\n",comp_size,nread);
    }
    uncompress((Bytef *)datap,&datasize,comp_data,comp_size);
    free(comp_data);
    return ( (int) datasize * sizeof(Bytef) );
}
