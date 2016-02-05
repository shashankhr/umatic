/*$Id: phasediag.h 892 2006-03-10 15:24:59Z rcatwood $*/
/****************************************************************/
/*   Copyright (c) 1998 - 2004 Dept. of Materials, ICSTM        */
/*   All Rights Reserved                                        */
/*   THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE OF ICSTM       */
/*   The copyright notice above does not evidence any           */
/*   actual or intended publication of such source code,        */
/*   and is an unpublished work by Dept. of Materials, ICSTM.   */
/*   This material contains CONFIDENTIAL INFORMATION that       */
/*   is the property of Imperial College. Any use,              */
/*   duplication or disclosure not specifically authorized      */
/*   by Imperial College is strictly prohibited.                */
/****************************************************************/
/* This code is part of the umats routines developed at in the  */
/* Materials Processing Group, Dept. of Materials, ICSTM.       */
/*      email p.d.lee or r.atwood @ic.ac.uk for details         */
/****************************************************************/
/*This file was created at Mon Jul 12 14:47:05 BST 2004 by rcatwood on hive.beowulf.cluster */
#ifndef PHASEDIAG_H
#define PHASEDIAG_H
#define PHASEDIAG_H_REV "phasediag.h $Id: phasediag.h 892 2006-03-10 15:24:59Z rcatwood $"
#define PHASEDIAG_ADD (10) 
#define PHASEDIAG_NPERPOINT (5)
#define MAX_NAME_LEN (32)

/* the max can be changed if the datatype is changed*/
/* for Pd_num_t */
#define MAX_PD_NUM (255)
typedef unsigned char Pd_num_t ;

typedef enum {
P1C1,P2C1,P1C2,P2C2,P3C1,P3C2
} Pha_Info;

/* structure for transferring the information at a phase diagram point*/
typedef struct pd_point_struct{
   CA_FLOAT temp;
   CA_FLOAT conc1;
   CA_FLOAT conc2;
   Pd_num_t region;
   CA_FLOAT p1c1;
   CA_FLOAT p2c1;
   CA_FLOAT p1c2;
   CA_FLOAT p2c2;
   /* only used if the point is in a tie triangle */
   CA_FLOAT p3c1;
   CA_FLOAT p3c2;
}Pd_point;

typedef struct tie_triangle_struct{
   Pd_num_t region;
   CA_FLOAT temp;
   CA_FLOAT p1c1;
   CA_FLOAT p2c1;
   CA_FLOAT p1c2;
   CA_FLOAT p2c2;
   CA_FLOAT p3c1;
   CA_FLOAT p3c2;
}Pd_tri;

typedef struct component_struct{
   Pd_num_t my_num;
   char name[MAX_NAME_LEN];
}Comp;
/*information concerning a specific phase */
typedef struct phase_struct{
   Pd_num_t my_num;
   char name[MAX_NAME_LEN];
}Pha;

/* information concerning a phase diagram region */
typedef struct region_struct{
   Pd_num_t my_num;
   int np;        /* number of phases */
   Pha * phases;  /* dynamically sized list of pointers */

   /* parameters describing the array of tie triangle data */
   int ntie;
   CA_FLOAT tmin;
   CA_FLOAT tmax;
   Pd_tri * tiedata; /* tie triangle at each temperature */
   
}Reg;


/* An isothermal section of the phase diagram */
typedef struct tslice_struct{
   int ncomp;
   CA_FLOAT my_temp;

   CA_FLOAT cmax[2];
   CA_FLOAT cmin[2];
   CA_FLOAT cstep[2]; /* only use first element if 2-comp */
   int n_csteps[2];

   int ndata;
   int n_alloc;

   Pd_num_t * regdata; /* large array of region number at every point*/

   /* array of tie triangle data */
   /*
      Pd_tri * tiedata; 
   */

   CA_FLOAT ** data;   /* only allocated when needed */
}Tslice;


/* the parent structure holding the whole phase diagram */
typedef struct phdiag_struct {
   char header[255];
   int ncomp;
   CA_FLOAT tmin;
   CA_FLOAT tmax;
   CA_FLOAT tstep;
   int n_regions;
   int n_c_alloc;
   int n_r_alloc;
   int n_phases;
   int n_p_alloc;
   int n_tslice; /* number of isothermal sections actuall stored */
   int n_alloc;  /* number of isothermal pointers allocated */

   Reg * regions;
   Pha * phases;
   Comp * comp; /* identifiers for the components */
   Tslice * tslices;
   char tailer[255];
}Phdiag;
#endif /*PHASEDIAG_H*/
