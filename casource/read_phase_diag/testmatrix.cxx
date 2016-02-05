/*$Id: testmatrix.cxx 1086 2007-10-10 10:55:39Z  $*/
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
/*This file was created at Thu Jul 15 17:09:31 BST 2004 by rcatwood on hive.beowulf.cluster */
#define MAKETERNARY_C_REV "maketernary.c $Id: testmatrix.cxx 1086 2007-10-10 10:55:39Z  $"
/* Include the common header for this project */
#include <stdio.h>
#include <math.h>
#include "../machine.h"
#include "vec.h"

int main (int argc, char * argv[]){
   vec p1(12,0,550);
   vec p2(0,6,525);
   vec p3(14,8,500);
   vec alpha,beta;
   Mat   m1(3,2);
   Mat  m2 (2,3);
   Mat  m3(3);
   FLOAT indata[9];

   FLOAT Tpure=650,m_a=-12,m_b=-6,k_a=0.1,k_b=0.15;
   FLOAT Tliq,Tsol;
   FLOAT Te,Tmin=400,Tmax=660,Tstep=1;
   FLOAT Ca,CaMin=0,CaMax=10,CaStep=.1;
   FLOAT Cb,CbMin=0,CbMax=10,CbStep=.1;
   FLOAT Ceut_a=12,Ceut_b=6;
   FLOAT Teut_a,Teut_b;
   FLOAT maxeut;
   FLOAT S_a,S_b;


   unsigned char *data,*datap;
   int ndata = 1;
   int nsteps,ntsteps,ncasteps,ncbsteps;
   FLOAT i,j,k;
   FILE * fp;
   int slice;

   S_a = Ceut_a*k_a;
   S_b = Ceut_b*k_b;
   Teut_a = Tpure + m_a * Ceut_a;
   Teut_b = Tpure + m_b * Ceut_b;
   maxeut = MAX(Teut_a,Teut_b);



   ntsteps = (Tmax-Tmin)/Tstep;
   ndata *= ntsteps;

   ncasteps = (CaMax - CaMin)/CaStep;
   ndata *= ncasteps;

   ncbsteps = (CbMax - CbMin)/CbStep;
   ndata *= ncbsteps;
   slice = ncasteps*ncbsteps;
   
   data = (unsigned char *)calloc(ndata,sizeof(unsigned char *));
   datap = data;

   printf("A,B,T: %i %i %i\n",ncasteps,ncbsteps,ntsteps);

   Te=Tmin;
   for (k=0;k<ntsteps;k++){
      Cb=CbMin;
      for(j=0;j<ncbsteps;j++){
         Ca = CaMin;
         for(i=0;i<ncasteps;i++){
         Tliq = Tpure + m_a * Ca + m_b * Cb;
         Tsol = Tpure + (m_a/k_a) * Ca + (m_b/k_b) * Cb;

         if (Te > Tliq) *datap = 0;
         else if (Te> Tsol) *datap = 1;
         else *datap=2;

         datap++;
         Ca += CaStep;
         }
      Cb += CbStep;
      }
   Te += Tstep;
   }

      int n_isteps = 1000;
      FLOAT Ca_istep = (Ceut_a - S_a) / n_isteps;
      FLOAT Cb_istep = (Ceut_b - S_b) / n_isteps;
      FLOAT ia,ib,Ta,Tb;
      FLOAT Cpa,Cpb;
      FLOAT dT = Teut_a - Teut_b;
      FLOAT lp,ratio,Lab,Dca;
      int ctr=0,idx;

      int Ts;
      int Cas;
      int Cbs;
      Cb=S_b;
      Ca=S_a;
      for(j=0;j<n_isteps;j++){

         Lab = sqrt(Ca*Ca+Cb*Cb);
         for(i=0;i<n_isteps;i++){
            Dca = i * (Ca/n_isteps);
            Cpa = Ca - Dca;
            Cpb = i*(Cb/n_isteps);
            lp = sqrt(Dca * Dca + Cpb * Cpb);

            ratio = lp/Lab;

            Te = Teut_a - ratio * (Teut_a - Teut_b);
            Ts = (int)(floor( (Te-Tmin) / Tstep));
            //printf("%i Lab,ratio,lp,Cpa,Cpb,Te: %g,%g,%g,%g,%g,%g\n",ctr++,Lab,ratio,lp,Cpa,Cpb,Te);
            Cas = (int)(floor( ((Cpa - CaMin )/ CaStep)));
            Cbs = (int)(floor( ((Cpb - CbMin)/ CbStep)));
            idx= Ts * slice + Cbs * ncbsteps + Cas;
            if (idx > 0 && idx < ndata && Cas <ncasteps && Cbs < ncbsteps && Ts < ntsteps ){
               if (*(data + idx) == 1){
               *(data + idx) = 5;
               }
            }

         }
      Ca += Ca_istep;
      Cb += Cb_istep;
      }

   fp = fopen("pd.raw","w");
   fwrite(data,sizeof(unsigned char),ndata,fp);
   fclose (fp);










   indata[0] = -1;
   indata[1] = 1;
   indata[2] = -2;
   indata[3] = 2;
   indata[4] = 1;
   indata[5] = 3;
   indata[6] = 2;
   indata[7] = 0;
   indata[8] = -1;


   m1.set_Mat(indata);
   m1.printMat(stdout);

   m2.set_Mat(indata);
   m3 = (m1 * m2);
   m1.printMat(stdout);
   m2.printMat(stdout);
   m3.printMat(stdout);

   




   p1.list(stdout);
   printf("\n");

   p2.list(stdout);
   printf("\n");

   alpha = p2 - p3;
   beta = p1 - p3;

   alpha.list(stdout);
   printf("\n");


   printf("A,B,T: %i %i %i\n",ncasteps,ncbsteps,ntsteps);
}


