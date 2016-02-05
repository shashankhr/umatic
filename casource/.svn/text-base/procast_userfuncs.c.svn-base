/* */
/* This is a stub for the user functions called by the commercial software package PROCAST 
   distributed by ESI */
/*END of LICENSE NOTICE*/

#include <stdio.h>
#include <stdlib.h>

#define real double
#ifdef EXTERNAL_HEATFLUX
real func_heatflux(char*, int, real, real, real, real, real, real, int);
#endif
#ifdef EXTERNAL_HEATTRANSFER
real func_convehtransfer(char*, int, real, real, real, real, real, real, int);
real func_texternal(char*, int, real, real, real, real, real, real, int);
#endif

void func_externalcompute(char*,int,int,real);
extern real usertemp1(int);
extern real userfs1(int);
extern real uservx1(int);
extern real uservy1(int);
extern real uservz1(int);
extern int nodNum (real,real,real,int,real* ,real* , real* );


/*
 *    function called at the beginning, the end of the calculation 
 *    and at the beginning of each timestep
 *    loop = 0 : function is called at the start of the computation
 *    loop = 1 : function is called at the start of the current timestep
 *    loop = 2 : function is called at the end of the current timestep
 *    loop = 3 : function is called at the end of the computation
 *
 */
void func_externalcompute(
   char  prefix[],     /* case name  */
   int  loop,          /* loop value : 0/1/2/3  */
   int  timestep,      /* current timestep      */
   real  time)         /* current time          */
{
/* ------------- Do not change anything above this line ------------- *
 * ------------- Program your function below this line  ------------- */

printf("Running the  uMatIc external subroutine\n");
printf("prefix:%s\n",prefix);
printf("Loop value:%i\n",loop);
printf("timestep:%i\n",timestep);
printf("time:%g\n",time);
/* if loop ==0 mignt mean the initial step, and could be used for memory allocation */
	if (loop == 2){
	   umat_extern_step();
	}
/* if loop == 3 might be useful for freeing the pointers closing output files */	
printf("Finished the uMatIc external subroutine\n");
}/* end of func_externalcompute */

/*
 *    convective heat transfer coefficient (applied on external surfaces)
 */

#ifdef EXTERNAL_HEATTRANSFER
real func_convehtransfer(
   char  prefix[],     /* case name  */
   int   dimension,    /* 2 = 2D ; 3 = 3D */
   real  temp,         /* current temperature */
   real  fs,           /* current fraction of solid */
   real  time,         /* current time */
   real  x_coor,       /* local coordinates: x */
   real  y_coor,       /* local coordinates: y */
   real  z_coor,       /* local coordinates: z */
   int   numBC)        /* number of boundary condition */
{
/* ------------- Do not change anything above this line ------------- *
 * ------------- Program your function below this line  ------------- */





   return (50.0);
}

/*
 *    ambient temperature (applied on external surfaces)
 */

real func_texternal(
   char  prefix[],     /* case name  */
   int   dimension,    /* 2 = 2D ; 3 = 3D */
   real  temp,         /* current temperature */
   real  fs,           /* current fraction of solid */
   real  time,         /* current time */
   real  x_coor,       /* local coordinates: x */
   real  y_coor,       /* local coordinates: y */
   real  z_coor,       /* local coordinates: z */
   int   numBC)        /* number of boundary condition */
{
/* ------------- Do not change anything above this line ------------- *
 * ------------- Program your function below this line  ------------- */



/* ------------ Do not forget to remove the call to exit ------------ *
 * ------------ hereafter before running the calculation ------------ */

   return (1550.0);
}

#endif

/*
 *    heat flux coefficient (applied on external surfaces)
 */

#ifdef EXTERNAL_HEATFLUX
real func_heatflux(
   char  prefix[],     /* case name  */
   int   dimension,    /* 2 = 2D ; 3 = 3D */
   real  temp,         /* current temperature */
   real  fs,           /* current fraction of solid */
   real  time,         /* current time */
   real  x_coor,       /* local coordinates: x */
   real  y_coor,       /* local coordinates: y */
   real  z_coor,       /* local coordinates: z */
   int   numBC)        /* number of boundary condition */
{
/* ------------- Do not change anything above this line ------------- *
 * ------------- Program your function below this line  ------------- */



/* ------------ Do not forget to remove the call to exit ------------ *
 * ------------ hereafter before running the calculation ------------ */
   return 0;
}

/***************** functions missing from exported external functions */
#endif

