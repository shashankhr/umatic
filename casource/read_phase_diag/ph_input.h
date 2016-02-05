/*$Id: ph_input.h 892 2006-03-10 15:24:59Z rcatwood $*/
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
/*This file was created at Tue Jul 27 12:08:20 BST 2004 by rcatwood on hive.beowulf.cluster */
#ifndef PH_INPUT_H
#define PH_INPUT_H
#define PH_INPUT_H_REV "ph_input.h $Id: ph_input.h 892 2006-03-10 15:24:59Z rcatwood $"
     typedef struct ph_input_struct{
     char * filename;
     int ncomp;
     double Tstart;
     double      Tmin;
     double      Tend;
     double      Tmax;
     double      Tstep;

     double      Climits[6];
       
     int      nphases ;
     char **  phase;
     char ** components;
     int      nregions; 
      /* region number nphases phase phase phase */ 
     Pd_num_t   (*  region) [6]; 
}PhInp;

#endif /*PH_INPUT_H*/
