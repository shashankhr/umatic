
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

/*RCS Id:$Id: rcs_id.h 1356 2008-08-18 13:41:15Z  $*/

#ifndef RCS_ID_H
#define RCS_ID_H

extern char const *alloc_multicomp_c();
extern char const *checkgas_c();
extern char const *checks_c();
extern char const *fg_read_c();
extern char const *grain_functions_c();
extern char const *handlers_c();
extern char const *init_nuc_mould_c();
extern char const *rcs_id_bigblock_c();
extern char const *rcs_id_c();
extern char const *rcs_id_umat_histo_c();
extern char const *rcs_id_calc_sb_c();
extern char const *rcs_id_umat_solid_c();
extern char const *rcs_id_castats_c();
extern char const *rcs_id_umat_wrapper_c();
/*extern char const *rcs_id_checkblock_c();*/
extern char const *rcs_id_curvature_c();
extern char const *rcs_id_curvature_3D_c();
extern char const *rcs_id_find_max_c();
extern char const *rcs_id_find_nmols_c();
extern char const *rcs_id_findroot_c();
extern char const *rcs_id_freeblock_c();
extern char const *rcs_id_gaussdev_c();
extern char const *rcs_id_getxyz_c();
extern char const *rcs_id_initcube_c();
extern char const *rcs_id_initface_c();
extern char const *rcs_id_init_sb_neigh_c();
extern char const *rcs_id_nbhd_def_c();
extern char const *rcs_id_nuc_functions_c();
extern char const *rcs_id_nuc_lookup_c();
extern char const *rcs_id_open_sb_c();
extern char const *rcs_id_output_ex_c();
extern char const *rcs_id_output_img_c();
extern char const *rcs_id_pore_c();
extern char const *rcs_id_pore_routines_c();
extern char const *rcs_id_pore_setup_c();
extern char const *rcs_id_props_c();
extern char const *rcs_id_rand_square_c();
extern char const *rcs_id_rand_step_c();
extern char const *rcs_id_read_blocks_c();
extern char const *rcs_id_read_cap_umat_c();
extern char const *rcs_id_read_ctrl_c();
extern char const *rcs_id_readgeoplus_c();
extern char const *rcs_id_readmat_c();
extern char const *rcs_id_read_sb_c();
extern char const *rcs_id_recr_extra_c();
extern char const *rcs_id_sb_boundary_c();
extern char const *rcs_id_sb_umat_step_c();
extern char const *rcs_id_sb_diffuse_alloy_c();
extern char const *rcs_id_sb_diffuse_gas_c();
extern char const *rcs_id_sb_line_int_c();
extern char const *rcs_id_sb_nuc_c();
extern char const *rcs_id_sb_temp_calc_c();
extern char const *rcs_id_subblock_c();
extern char const *rcs_id_SurCellRoutines_c();
extern char const *rcs_id_thermo_trace_calc_c();
extern char const *rcs_id_thermo_trace_init_c();
extern char const *rcs_id_trans_interp_calc_c();
extern char const *rcs_id_wfact_r_calc_c();
extern char const *rcs_id_wfact_z_calc_c();
extern char const *rcs_id_write_blocks_c();
extern char const *read_umat_extern_c();
extern char const *safeopen_c();
extern char const *sb_decentred_step_c();
extern char const *sb_diffuse_alloy_decentred_c();
extern char const *setup_mould_src_function_c();
extern char const *setup_temp_func_c();
extern char const *step_output_c();
extern char const *window_move_c();
extern char const *write_bb_values_c();
extern char const *write_block_values_c();
extern char const *write_ctrl_values_c();
extern char const *write_matprop_values_c();
extern char const *write_nprop_values_c();
extern char const *write_POREprop_values_c();
extern char const *write_pprop_values_c();

extern char const *rcs_id_icopymat_c();
extern char const *rcs_id_fcopymat_c();

#ifdef EXTERNAL_CA
	extern char const *umat_feedback_c();
	extern char const *user_rtn_cell_c();
#endif /* EXTERNAL_CA */

#endif
/*
*/

