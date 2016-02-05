
/****************************************************************/
/*      Copyright (c) 1993 Peter D Lee                          */
/*      Copyright (c) 1998 Dept. of Materials, ICSTM            */
/*      All Rights Reserved                                     */
/*      The copyright notice above does not evidence any        */
/*      actual or intended publication of such source code,     */
/*      and is an unpublished work by Dept. of Materials, ICSTM.*/
/*      continuing D Phil work from University of Oxford        */
/****************************************************************/
/* This code is part of the umats routines developed at in the  */
/* Materials Processing Group, Dept. of Materials, ICSTM.       */
/*      email p.d.lee or r.atwood @imperial.ac.uk for details   */
/****************************************************************/

/********************************************************************************/
/*  This version is distributed under a BSD style public license, as follows:   */
/*                                                                              */
/*  Copyright (c) 2007, Dept. of Materials, Imperial College London             */
/*  All rights reserved.                                                        */
/*  Redistribution and use in source and binary forms, with or without          */
/*  modification, are permitted provided that the following conditions          */
/*  are met:                                                                    */
/*                                                                              */
/*  * Redistributions of source code must retain the above copyright            */
/*  notice, this list of conditions and the following disclaimer.               */
/*                                                                              */
/*  * Redistributions in binary form must reproduce the above                   */
/*  copyright notice, this list of conditions and the following                 */
/*  disclaimer in the documentation and/or other materials provided             */
/*  with the distribution.                                                      */
/*                                                                              */
/*  * Neither the name of the Dept. of Materials, Imperial College London, nor  */
/*  the names of its contributors may be used to endorse or promote products    */
/*  derived from this software without specific prior written permission.       */
/*                                                                              */
/*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS         */
/*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT           */
/*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR       */
/*  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT        */
/*  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,       */
/*  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    */
/*  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR      */
/*  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF      */
/*  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING        */
/*  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS          */
/*  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                */
/********************************************************************************/
/*END of LICENSE NOTICE*/

#ifndef PROPS_H
#define PROPS_H
CA_FLOAT get_d_arr(CA_FLOAT temp,CA_FLOAT Qa, CA_FLOAT D_o);
CA_FLOAT get_ds_ti(CA_FLOAT temp) ;
CA_FLOAT get_dl_ti(CA_FLOAT temp) ;
CA_FLOAT get_reac_rate(CA_FLOAT temp) ;
CA_FLOAT get_dg_ti(CA_FLOAT temp) ;
#ifdef CELL_DIFF_ARR
CA_FLOAT get_ds(CA_FLOAT temp,Solute_props *sp);
CA_FLOAT get_dl(CA_FLOAT temp,Solute_props *sp);
#else
CA_FLOAT get_ds(CA_FLOAT temp);
CA_FLOAT get_dl(CA_FLOAT temp);
#endif
CA_FLOAT getav_d(CA_FLOAT dl,CA_FLOAT ds,CA_FLOAT fs);
CA_FLOAT schiel(CA_FLOAT Tcell);
CA_FLOAT find_sch_conc(CA_FLOAT tempK,CA_FLOAT fl);
CA_FLOAT growth(CA_FLOAT gg_const,CA_FLOAT gg_cub, CA_FLOAT dt, CA_FLOAT Tunder);
CA_FLOAT cell_liq_calc(CA_FLOAT cell_conc,Mat_str *mp);
CA_FLOAT cell_liq_calc_multi(CA_FLOAT **conc_multi, int index, BB_struct *bp);
CA_FLOAT cell_liq_calc_poly(CA_FLOAT cell_conc, Mat_str * mp, int ele_1);
CA_FLOAT find_sat(Mat_str *mp,CA_FLOAT cell_tempK,CA_FLOAT cell_si,CA_FLOAT cell_fs);
CA_FLOAT find_sat_poly (BB_struct * bp, int sbnum, CA_FLOAT cell_tempK);

#ifdef SENS
CA_FLOAT find_sat_multi(BB_struct *bp,CA_FLOAT cell_tempK,CA_FLOAT **conc_multi,CA_FLOAT cell_fs,CA_FLOAT sens_P,int index);
#else
CA_FLOAT find_sat_multi(BB_struct *bp,CA_FLOAT cell_tempK,CA_FLOAT **conc_multi,CA_FLOAT cell_fs,int index);
#endif
#ifdef SENS
CA_FLOAT r_from_ideal_gas (CA_FLOAT temperature, \
                        CA_FLOAT n_gmol,CA_FLOAT gammaSI, CA_FLOAT old_rSI,CA_FLOAT sens_PRESS);
#else
CA_FLOAT r_from_ideal_gas (CA_FLOAT temperature, CA_FLOAT n_gmol,CA_FLOAT gammaSI, CA_FLOAT old_rSI);
#endif
CA_FLOAT get_mat(CA_FLOAT thresh,CA_FLOAT delta_t, CA_FLOAT ss, CA_FLOAT oldmat);
#endif /*PROPS_H*/
