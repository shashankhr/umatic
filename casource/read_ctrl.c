
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

/****************************************************************/
/* READ_CTRL.C:   (Part of CA code)                             */
/* Subroutine to read the main controlling info from a file.    */
/* The file is formated, using the style:                       */
/*    # as first character:     Comment                         */
/* and values are input in the format:                          */
/*    command value  #comments                                  */
/****************************************************************/
/****************************************************************/
/* Written by Peter D. Lee & Robert C. Atwood, Imperial College */
/* Jul 1, 1998                                                  */
/****************************************************************/

/****************************************************************/
/* Versions maintained with RCS                                 */
/* Version 1.0: Aug. 13, pdl                                    */
/****************************************************************/
/*RCS Id:$Id: read_ctrl.c 1382 2008-09-24 14:59:34Z  $*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include "getcflags.h"

#include "machine.h"
#include "blocks.h"

int combine_ctrl(Ctrl_str * cp);
int read_ctrl(char *filename, Ctrl_str * cp);
int read_ctrl_vals(char *filename, Ctrl_str * cp);
void set_ctrl_defaults(Ctrl_str * cp);

void set_ctrl_defaults(Ctrl_str * cp) {
    /* copy the compilation flags into a variable -- requires */
    char *Cflags;
    int i;

    Cflags = strdup(GETCFLAGS);

    /*********************************************************/
    /* Set all the default values...                         */
    /*********************************************************/
    /* from read_ctrl.h */
    cp->solo = D_SOLO;
    cp->cap = D_CAP;
    cp->post = D_POST;
    cp->do_conc_prof = FALSE;
    cp->phase_diag_on = D_PHASE_DIAG;
    cp->global_undercooling = D_GLOBAL_UNDERCOOLING;
    cp->particle = D_PARTICLE;
    cp->restart = D_RESTART;
    cp->restart_gas_on = D_RESTART_PORE;
    cp->restart_gas_on = D_RESTART_GAS;
    cp->diffuse = D_DIFFUSE;
    cp->das_limrad = D_DAS_LIMRAD;
    cp->diffuse_alloy = D_DIFFUSE_ALLOY;
    cp->diffuse_alloy_multi = D_DIFFUSE_ALLOY_MULTI;
    cp->temp_dep_diff = D_TEMP_DEP_DIFF;
    cp->interpolate = D_INTERPOLATION;
    cp->thermocalc = D_THERMOCALC;
    cp->diffuse_step = D_DIFFUSE_STEP;
    cp->temp_lookup = 0;
    cp->window_moving = D_WINDOW_MOVING; /*by Wei WANG on 11-07-02 */
    cp->init_cont = D_INIT_CONT; /*by Wei WANG on 11-07-02 */
    cp->decentred_octahedron = D_DECENTRED_OCTAHEDRON; /*by Wei WANG on 11-07-02 */
    cp->random_angles = 1; /* choose wei original option as default */
    cp->melt_back = D_MELT_BACK; /* melt back option, xly 20040802 */
    cp->isotherm_curv = D_ISOTHERM_CURV; /* option for producing transient isotherm curvature, xly 20040802 */
    cp->pr_lookup = 0;
    cp->diff_rgbmode = D_DIFF_RGBMODE;
    cp->rgbgrey = FALSE;
    cp->input = D_INP;
    cp->extrudemould = 0;
    cp->t_input = FALSE;
    cp->time_dump = FALSE;
    cp->time_exp = 0;
    cp->time_unit = 1;
    cp->fgrid_input = FALSE;
    cp->con_cast = D_CON_CAST;
    cp->n_neigh = NEIGH_6;
    cp->umat_method = CA_PULL;
    cp->seed = D_SEEDVAL;
    cp->scheil = FALSE;
    cp->swap_axes = 0;
    cp->gradtilt = D_USE_GRAD_TILT;

    for (i = 0; i < 3; i++) {
        cp->grad_angle[i] = D_GRAD_ANGLE;
    }

    cp->curvature_2D = 0;
    cp->curvature_3D = 0;

    cp->fn_cap = strdup(D_CAP_FILE);

    /* default file names and allocation of strings */
    cp->fn_block_restart = (char *) calloc(MAX_STRING_LEN, sizeof (char));
    sprintf(cp->fn_block_restart, "%s\0", (D_BLOCK_RESTART_FILE));
    cp->fn_geo = (char *) calloc(MAX_STRING_LEN, sizeof (char));
    sprintf(cp->fn_geo, "%s\0", (D_GEO_FILE));
    cp->fn_mat = (char *) calloc(MAX_STRING_LEN, sizeof (char));
    sprintf(cp->fn_mat, "%s\0", (D_MAT_FILE));
    cp->fn_base = (char *) calloc(MAX_STRING_LEN, sizeof (char));
    sprintf(cp->fn_base, "%s\0", "CHOOSE_a_better_file_name");
    cp->fn_phadia = (char *) calloc(MAX_STRING_LEN, sizeof (char));
    sprintf(cp->fn_phadia, "%s\0", (D_PHADIA_FILE));
    cp->fn_inp = (char *) calloc(MAX_STRING_LEN, sizeof (char));
    sprintf(cp->fn_inp, "%s\0", (D_INP_FILE));
    cp->fn_fgrid = (char *) calloc(MAX_STRING_LEN, sizeof (char));
    sprintf(cp->fn_fgrid, "%s\0", (D_FG_FILE));

    cp->fg_tr = TRANSIENT;

    /* output variables */
    cp->write_block = TRUE;
    cp->excel = D_EXCEL;
    cp->nsbslice = 0;
    cp->tempslice = D_TEMPSLICE;
    cp->nbbslice = 0;
    cp->grainslice = 0;
    cp->slice_dmp_freq = D_SLICE_FREQ;
    cp->blk_dmp_freq = D_BLK_FREQ;
    cp->floatdump = D_CA_FLOATDUMP;
    cp->scr_dmp_freq = D_SCREEN_FREQ;
    cp->rgbmode = D_RGB_RANDOM;

    /* test mode variables */
    cp->fixed_Pore = FALSE;
    cp->nfPore = 0;
    cp->fixed_nuc = FALSE;
    cp->nfnuc = 0;
    cp->coolrate = FALSE;
    cp->delT = D_DEL_TEMP;
    cp->fs_finish = D_FS_FINISH;
}

int read_ctrl_vals(char *filename, Ctrl_str * cp) {
    char *line;
    char *token;
    FILE *fp;
    char *sep;
    int i, error;
    int rflag = 0;
    line = (char *) calloc(MAX_STRING_LEN, sizeof (char));
    sep = (char *) calloc(MAX_STRING_LEN, sizeof (char));
    sprintf(sep, " ,;\t\n\r");

    /*********************************************************/
    /* Open the control file                                 */
    /*********************************************************/
    if ((fp = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Error: read_ctrl: can't open input file [%s]\n", filename);
        return (1);
    }

    while (fgets(line, MAX_STRING_LEN, fp) != NULL) {
        /* ignore comment and blank lines */
        if (line[0] == '%' || line[0] == '#' || (token = strtok(line, sep)) == NULL) {
            continue;

            /*********************************************************/
            /* All values in the Ctrl_str structure                  */
            /*********************************************************/
            /* SoloMode */
        } else if (strcasecmp(token, "SoloMode") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->solo = atoi(token);
            /* CAPMode */
        } else if (strcasecmp(token, "CAPMode") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->cap = atoi(token);
            /* PostProcessingMode */
        } else if (strcasecmp(token, "PostProcessingMode") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->post = atoi(token);
            /* CA_FEEDBACK */
        } else if (strcasecmp(token, "CA_FEEDBACK") == 0) {
            if ((token = strtok(NULL, " ,;\t")) != NULL)
                cp->umat_feedback = atoi(token);

            /*Nucleation on the mould surface */
        } else if (strcasecmp(token, "MouldNuc") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->mould_nuc = atoi(token);

            /* Source of solute at the mould */
            /* uses enum variable for the type of function */
            /* typedef in machine.h */
        } else if (strcasecmp(token, "MouldSrc") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->mould_src = (atoi(token));

            /*Perturb the mould source ? */
        } else if (strcasecmp(token, "MouldSrcPert") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->mould_src_pert = atoi(token);

            /* Source magnitude of solute at the mould */
        } else if (strcasecmp(token, "MouldSourceValue") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->mould_source_value = atof(token);

            /* Source magnitude of solute prerturbation at the mould */
        } else if (strcasecmp(token, "MouldSourcePertValue") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->mould_source_pert = atof(token);

            /* Source spatial frequency of solute at the mould */
        } else if (strcasecmp(token, "MouldSourceFreq") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->mould_source_freq = atof(token);

            /* ConCastMode -- directional solidification */
        } else if (strcasecmp(token, "ConCastMode") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->con_cast = atoi(token);
        } else if (strcasecmp(token, "GradAngle") == 0) {
            for (i = 0; i < 3; i++) {
                if ((token = strtok(NULL, " ,;\t")) != NULL)
                    cp->grad_angle[i] = (double) (atof(token));
            }
            /* Particle */
        } else if (strcasecmp(token, "Particle") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->particle = atoi(token);
            /* das_limrad */
        } else if (strcasecmp(token, "das_limrad") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->das_limrad = atoi(token);
            /* Diffuse gas */
        } else if ((strcasecmp(token, "diffuse") == 0) ||
                /* alternate key for gas diffusion */
                (strcasecmp(token, "diffuseGas") == 0)) {
            int old_dif;

            if (cp->restart)
                old_dif = cp->diffuse;
            if ((token = strtok(NULL, sep)) != NULL)
                cp->diffuse = atoi(token);
            /* check for restart with gas turned on */
            if (cp->restart && (old_dif == 0) && (cp->diffuse == 1)) {
                cp->restart_gas_on = 1;
            }

            /* Diffuse Alloy */
        } else if (strcasecmp(token, "diffuseAlloy") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->diffuse_alloy = atoi(token);

            /* Number of alloying elements *//* L. THUINET 07/02/05 */
        } else if (strcasecmp(token, "NUMCOMP") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->NUM_COMP = atoi(token);

            /* Number of solid phases for multiphase solidification *//* L. THUINET 04/05 */
        } else if (strcasecmp(token, "NUMPHS") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->NUM_PHS = atoi(token);

            /* Number of equilibria for multiphase solidification *//* L. THUINET 04/05 */
        } else if (strcasecmp(token, "NUMEQ") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->NUM_EQ = atoi(token);

            /* Number of distinct gaussian nucleation laws *//* L. THUINET 10/05 */
        } else if (strcasecmp(token, "NUM_NUC_LAW") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->NUM_NUC_LAW = atoi(token);

            /* Diffuse Alloy Poly *//* L. THUINET 04/02/05 */
        } else if (strcasecmp(token, "diffuseAlloyPoly") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->diffuse_alloy_poly = atoi(token);

            /* Phase Diagram */
        } else if (strcasecmp(token, "PhaseDiagram") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->phase_diag_on = atoi(token);

            /* Global Undercooking */
        } else if (strcasecmp(token, "GlobalUndercooling") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->global_undercooling = atoi(token);
            /* Show Eutectic */
        } else if (strcasecmp(token, "ShowEutectic") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->show_eut = atoi(token);
            /* fraction solid stop */
        } else if (strcasecmp(token, "FsFinish") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fs_finish = (CA_FLOAT) atof(token);

            /* swap-xy axis for directional */
        } else if (strcasecmp(token, "SwapXy") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->swap_axes = (CA_FLOAT) atof(token);
        } else if (strcasecmp(token, "SwapAxes") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->swap_axes = (CA_FLOAT) atof(token);
            /* use a tilted gradient in directional */
        } else if (strcasecmp(token, "GradTilt") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->gradtilt = atoi(token);
            /* diffuse_step */
        } else if (strcasecmp(token, "diffuse_step") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->diffuse_step = atoi(token);

            /* Window Moving. by Wei WANG on 11-07-02 */
        } else if (strcasecmp(token, "Window_moving") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->window_moving = atoi(token);

            /* Initial or continuing calculation. by Wei WANG on 11-07-02 */
        } else if (strcasecmp(token, "Init_cont") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->init_cont = atoi(token);

            /* Decentred Octahedron. by Wei Wang on 11-07-02 */
        } else if (strcasecmp(token, "Decentred_octahedron") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->decentred_octahedron = atoi(token);
        } else if (strcasecmp(token, "Random_Angles") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->random_angles = atoi(token);
            /* option for change pulling velocity, for melt-back only, xly 20040802 */
        } else if (strcasecmp(token, "Melt_Back") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->melt_back = atoi(token);

            /* option for producing transient isotherm curvature, xly 20040802 */
        } else if (strcasecmp(token, "Isotherm_Curv") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->isotherm_curv = atoi(token);
            /* schiel */
        } else if (strcasecmp(token, "Scheil") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->scheil = atoi(token);

            /* N_Neighbours */
        } else if (strcasecmp(token, "N_Neighbours") == 0) {
            if ((token = strtok(NULL, sep)) != NULL) {
                cp->n_neigh = atoi(token);
                switch (cp->n_neigh) {
                    case NEIGH_6:
                        fprintf(stderr, "N_Neighbours: set equal to 6.\n");
                        break;
                    case NEIGH_8:
                        fprintf(stderr, "N_Neighbours: set equal to 8.\n");
                        break;
                    case NEIGH_10:
                        fprintf(stderr, "N_Neighbours: set equal to 10.\n");
                        break;
                    case NEIGH_26:
                        fprintf(stderr, "N_Neighbours: set equal to 27.\n");
                        break;
                    default:
                        fprintf(stderr, "ERROR: unknown neighbourhood [%d], set to 6.\n", cp->n_neigh);
                        cp->n_neigh = 6;
                        break;
                } /* end of Flag switch */

            }
            /* CAPFileName */
        } else if (strcasecmp(token, "CAPFileName") == 0) {
            free(cp->fn_cap);
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fn_cap = strdup(token);
            /* Restart Block file name */
        } else if (strcasecmp(token, "BlockRestartFileName") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                sprintf(cp->fn_block_restart, "%s\0", token);
            /* GeoFileName */
        } else if (strcasecmp(token, "GeoFileName") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                sprintf(cp->fn_geo, "%s\0", token);
            /* MatFileName */
        } else if (strcasecmp(token, "MatFileName") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                sprintf(cp->fn_mat, "%s\0", token);
            /* InpFileName */
        } else if (strcasecmp(token, "InpFileName") == 0) {
            if ((token = strtok(NULL, sep)) != NULL) {
                sprintf(cp->fn_inp, "%s\0", token);
                cp->input = TRUE;
            }
            /* ExtrudeMould */
        } else if ((strcasecmp(token, "ExtrudeMould") == 0) || (strcasecmp(token, "ExtrudeMold") == 0)) {
            if ((token = strtok(NULL, sep)) != NULL) {
                sprintf(cp->fn_inp, "%s\0", token);
                cp->input = TRUE;
            }
            /* temperatureinput filename(string) */
        } else if (strcasecmp(token, "TemperatureInputFile") == 0) {
            cp->t_input = TRUE;
            if ((token = strtok(NULL, " ,;\t")) != NULL) {
                cp->fn_t_input = strdup(token);
                fprintf(stderr, "readgeoplus: TemperatureInputFile set to [%s].\n", cp->fn_t_input);
            }
            /* BaseFileName */
        } else if (strcasecmp(token, "BaseFileName") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                sprintf(cp->fn_base, "%s\0", token);
            /* FgridFileName */
        } else if (strcasecmp(token, "FgridFileName") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                sprintf(cp->fn_fgrid, "%s\0", token);
            cp->fgrid_input = TRUE;
            /* FgridFileName */
        } else if (strcasecmp(token, "FgridTransient") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fg_tr = atoi(token);
            /* gas properties filename */
        } else if (strcasecmp(token, "GasPropsFile") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fn_solprops_gas = strdup(token);
            fprintf(stderr, "read_ctrl:gas solute properties file is %s\n", cp->fn_solprops_gas);

            /* alloy properties filename */

        } else if (strcasecmp(token, "AlloyPropsFile0") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fn_solprops_alloy[0] = strdup(token);
            fprintf(stderr, "read_ctrl:alloy solute properties file 0  is %s\n", cp->fn_solprops_alloy[0]);
        } else if (strcasecmp(token, "AlloyPropsFile1") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fn_solprops_alloy[1] = strdup(token);
            fprintf(stderr, "read_ctrl:alloy solute properties file 1 is %s\n", cp->fn_solprops_alloy[1]);
        } else if (strcasecmp(token, "AlloyPropsFile2") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fn_solprops_alloy[2] = strdup(token);
            fprintf(stderr, "read_ctrl:alloy solute properties file 2 is %s\n", cp->fn_solprops_alloy[2]);
        } else if (strcasecmp(token, "AlloyPropsFile3") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fn_solprops_alloy[3] = strdup(token);
            fprintf(stderr, "read_ctrl:alloy solute properties file 3 is %s\n", cp->fn_solprops_alloy[3]);
        } else if (strcasecmp(token, "AlloyPropsFile4") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fn_solprops_alloy[4] = strdup(token);
            fprintf(stderr, "read_ctrl:alloy solute properties file 4 is %s\n", cp->fn_solprops_alloy[4]);
        } else if (strcasecmp(token, "AlloyPropsFile5") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fn_solprops_alloy[5] = strdup(token);
            fprintf(stderr, "read_ctrl:alloy solute properties file 5 is %s\n", cp->fn_solprops_alloy[5]);
        } else if (strcasecmp(token, "AlloyPropsFile6") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fn_solprops_alloy[6] = strdup(token);
            fprintf(stderr, "read_ctrl:alloy solute properties file 6 is %s\n", cp->fn_solprops_alloy[6]);
        } else if (strcasecmp(token, "AlloyPropsFile7") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fn_solprops_alloy[7] = strdup(token);
            fprintf(stderr, "read_ctrl:alloy solute properties file 7 is %s\n", cp->fn_solprops_alloy[7]);
        } else if (strcasecmp(token, "AlloyPropsFile8") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fn_solprops_alloy[8] = strdup(token);
            fprintf(stderr, "read_ctrl:alloy solute properties file 8 is %s\n", cp->fn_solprops_alloy[8]);
        } else if (strcasecmp(token, "AlloyPropsFile9") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->fn_solprops_alloy[9] = strdup(token);
            fprintf(stderr, "read_ctrl:alloy solute properties file 9 is %s\n", cp->fn_solprops_alloy[9]);


            /* pressure lookup? */
        } else if (strcasecmp(token, "PrLookup") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->pr_lookup = atoi(token);
            /* Porosity */
        } else if (strcasecmp(token, "Pore") == 0) {
            int old_pore;

            if (cp->restart)
                old_pore = cp->pore;
            if ((token = strtok(NULL, sep)) != NULL)
                cp->pore = atoi(token);
            /* check for restart with pore turned on */
            if (cp->restart && (old_pore == 0) && (cp->pore == 1)) {
                cp->restart_pore_on = 1;
            }

            /* RandomSeed Number */
        } else if (strcasecmp(token, "RandSeedVal") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->seed = (int32_t) (atoi(token));
            if (cp->seed == -2) {
                int fd;

                fd = open("/dev/random", O_RDONLY);
                read(fd, &(cp->seed), sizeof (int32_t));
                fprintf(stderr, "Kernel random seed requested: using %i\n", cp->seed);
                close(fd);
            }
            /* PLIC Curvature in 3D option */
        } else if (strcasecmp(token, "Curvature_3D") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->curvature_3D = atoi(token);
            /* PLIC Curvature in 2D option */
        } else if (strcasecmp(token, "Curvature_2D") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->curvature_2D = atoi(token);



#ifdef JUNK
            /* ?? */
        } else if (strcasecmp(token, "??") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp-> ? ? = atoi(token);
#endif /* JUNK */

            /*********************************************************/
            /* All Output settings                                   */
            /*********************************************************/
            /* PrintExcel 1/0 (true/false) */
        } else if (strcasecmp(token, "PrintExcel") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->excel = atof(token);

            /* create 3d blocks for AVS */
        } else if (strcasecmp(token, "WriteBlock") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->write_block = atof(token);

            /* Dump pore infor for subblock */
        } else if (strcasecmp(token, "PoreDumpSb") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->pore_dump_sb = atoi(token);

            /* Dump pore infor for subblock */
        } else if (strcasecmp(token, "GrainSlice") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->grainslice = atoi(token);

            /* conc_prof sb# zslice# row# */
        } else if (strcasecmp(token, "ConcProfile") == 0) {
            error = FALSE;
            for (i = 0; i < 3; i++) {
                if ((token = strtok(NULL, sep)) != NULL) {
                    cp->do_conc_prof = TRUE;
                    cp->conc_prof[i] = atoi(token);
                } else {
                    rflag++;
                    error = TRUE;
                    fprintf(stderr, "Error: ConcProfile, must have three params: sb#, slice#,row#.\n");
                }
            }

            /* PrintSlice sb# zslice# */
        } else if (strcasecmp(token, "PrintSlice") == 0) {
            if (cp->nsbslice < MAX_CTRL) {
                error = FALSE;
                for (i = 0; i < 2; i++) {
                    if ((token = strtok(NULL, sep)) != NULL) {
                        cp->slice[cp->nsbslice][i] = atoi(token);
                    } else {
                        rflag++;
                        error = TRUE;
                        fprintf(stderr, "Error: PrintSlice, must have two params: sb#, slice#.\n");
                    }
                }
                if (!error)
                    cp->nsbslice++;
            } else { /* exceeded MAX_CTRL print commands */
                fprintf(stderr, "Error: PrintSlice, exceed Max of %d slices.\n", MAX_CTRL);
            }

            /* PrintBBSlice zslice# */
        } else if (strcasecmp(token, "PrintBBSlice") == 0) {
            if (cp->nbbslice < MAX_CTRL) {
                if ((token = strtok(NULL, sep)) != NULL) {
                    cp->bbslice[cp->nbbslice] = atoi(token);
                    cp->nbbslice++;
                } else {
                    fprintf(stderr, "Error: PrintBBSlice, must have two params: sb#, slice#.\n");
                }
            } else { /* exceeded MAX_CTRL */
                fprintf(stderr, "Error: PrintSlice, exceed Max of %d slices.\n", MAX_CTRL);
            }

            /*TempSlice int */
        } else if (strcasecmp(token, "TempSlice") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->tempslice = atoi(token);
            /*FloatDump int */
        } else if (strcasecmp(token, "BinDump") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->floatdump = atoi(token);
            /*TimeDump int */
        } else if (strcasecmp(token, "Time_Dump") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->time_dump = atoi(token);
            /*TimeExp int */
        } else if (strcasecmp(token, "Time_Exp") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->time_exp = atoi(token);
            cp->time_unit = pow(10.0, (double) cp->time_exp);
            /* SlicePFreq int */
        } else if (strcasecmp(token, "SlicePFreq") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->slice_dmp_freq = atoi(token);
            /* BlockPFreq int */
        } else if (strcasecmp(token, "BlockPFreq") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->blk_dmp_freq = atoi(token);
            /* SreenPFreq int -- here for backward compatibility with old ctrl files */
        } else if (strcasecmp(token, "SreenPFreq") == 0) {
            if ((token = strtok(NULL, sep)) != NULL) {
                fprintf(stderr, "WARNING: read_ctrl: Old erroneous keyword (SreenPFreq) used. \n");
                cp->scr_dmp_freq = atoi(token);
            }

            /* ScreenPFreq int */
        } else if (strcasecmp(token, "ScreenPFreq") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->scr_dmp_freq = atoi(token);

            /* RGBmode int */
        } else if (strcasecmp(token, "RGBmode") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->rgbmode = atoi(token);

            /* RGBmode int */
        } else if (strcasecmp(token, "RGBgrey") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->rgbgrey = atoi(token);

            /* RGBmode diffusion */
        } else if (strcasecmp(token, "diff_rgbmode") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->diff_rgbmode = atoi(token);
            /* LogMode diffusion */
        } else if (strcasecmp(token, "diff_log_disp") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->diff_log_disp = atoi(token);
            /* Ratio Disp diffusion */
        } else if (strcasecmp(token, "diff_ratio_disp") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->diff_ratio_disp = atoi(token);
            /* display max for gas */
        } else if (strcasecmp(token, "gas_disp_max") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->gas_disp_max = atof(token);
            /* display max for alloy */
        } else if (strcasecmp(token, "alloy_disp_max") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->alloy_disp_max = atof(token);
            /* cap display at max for diffusion */
        } else if (strcasecmp(token, "diff_disp_cap") == 0) {
            if ((token = strtok(NULL, sep)) != NULL)
                cp->diff_disp_cap = atoi(token);
            /*********************************************************/
            /* Testing modes...                                      */
            /*********************************************************/
            /* put in only a few fixed nuclei!                       */
            /*********************************************************/
            /* FixedPore nx ny nz  */
        } else if (strcasecmp(token, "FixedPore") == 0) {
            error = FALSE;
            if (cp->nfPore < MAX_CTRL) {
                for (i = 0; i < 4; i++) {
                    if ((token = strtok(NULL, sep)) != NULL) {
                        cp->nPsite[cp->nfPore][i] = atof(token);
                    } else {
                        rflag++;
                        error = TRUE;
                        fprintf(stderr, "Error: FixedPore, must have four params (cell index): nx,ny,nz,thresh.\n");
                    }
                }
                if (!error) {
                    cp->fixed_Pore = TRUE;
                    cp->nfPore++;
                }
            } else { /* exceeded MAX_CTRL nucleation points */
                fprintf(stderr, "Error: FixedPore, exceed Max: %d nucleation points.\n", MAX_CTRL);
            }
            /* FsPrint float */
        } else if (strcasecmp(token, "FsPrint") == 0) {
            error = FALSE;
            if (cp->nfsprint < MAX_CTRL) {
                if ((token = strtok(NULL, sep)) != NULL) {
                    cp->fsprint[cp->nfsprint] = (CA_FLOAT) atof(token);
                } else {
                    rflag++;
                    error = TRUE;
                    fprintf(stderr, "Error: FsPrint: No value \n");
                }

                if (!error) {
                    cp->nfsprint++;
                }
            } else { /* exceeded MAX_CTRL nucleation points */
                fprintf(stderr, "Error: FsPrint, exceed Max: %d fsprint points.\n", MAX_CTRL);
            }

            /* FixedNuc nx ny nz ang0 ang1 ang2 NucThreshold *//*modified by Wei WANG on 29-01-03 */
        } else if (strcasecmp(token, "FixedNuc") == 0) {
            error = FALSE;
            if (cp->nfnuc < MAX_CTRL) {
                /*THUINET 05/05*/
                /*for (i=0; i<7; i++) { */
                for (i = 0; i < 8; i++) {
                    /*FIN THUINET 05/05*/
                    if ((token = strtok(NULL, sep)) != NULL) {
                        cp->nsite[cp->nfnuc][i] = atoi(token);
                    } else {
                        rflag++;
                        error = TRUE;
                        /*THUINET 05/05*/
                        /*fprintf(stderr,"Error: FixedNuc, must have 7 params (cell index): nx,ny,nz,ang0,ang1,ang2, NucThreshold.\n"); */
                        fprintf(stderr, "Error: FixedNuc, must have 8 params (cell index): nx,ny,nz,ang0,ang1,ang2, NucThreshold,nat_sol.\n");
                        /*FIN THUINET 05/05*/
                    }
                }
                if (!error) {
                    cp->fixed_nuc = TRUE;
                    cp->nfnuc++;
                }
            } else { /* exceeded MAX_CTRL nucleation points */
                fprintf(stderr, "Error: FixedNuc, exceed Max: %d nucleation points.\n", MAX_CTRL);
            }
            /*********************************************************/
            /* Set the cooling rate to be constant everywhere!       */
            /*********************************************************/
            /* cooling rate CA_FLOAT */
        } else if (strcasecmp(token, "CoolingRate") == 0) {
            if ((token = strtok(NULL, sep)) != NULL) {
                cp->coolrate = TRUE;
                /* allow positive or negative numbers, always cool down! */
                cp->delT = -(ABS(atof(token)));
            }

            /*********************************************************/
            /* Unknown command                                       */
            /*********************************************************/
        } else {
            fprintf(stderr, "Warning: Unknown command: %s.\n", token);
        }
    } /* end while */

    /***********************************************************/
    /*  check for illogical combination                        */
    /***********************************************************/
    if (cp->phase_diag_on && !cp->diffuse_alloy) {
        fprintf(stderr, "ERROR: read_ctrl: cannot use phase diagram without alloy diff\n");
        fprintf(stderr, "ERROR: read_ctrl: Turning off phase diagram\n");
        cp->phase_diag_on = 0;
    }
    if (cp->pore && !cp->diffuse) {
        fprintf(stderr, "ERROR: read_ctrl: cannot use pore without gas diff\n");
        fprintf(stderr, "ERROR: read_ctrl: Turning off pore mode.\n");
        cp->pore = 0;
    }
    /*   if (cp->decentred_octahedron  && !cp->fixed_nuc){
       fprintf(stderr,"ERROR: read_ctrl: Decentered mode needs fixed nuclei\n");
       exit(0);
       } *//* removed by Wei WANG on 17-09-02 */
    /* combine the logical values to create some useful flags */
    rflag += combine_ctrl(cp);

    /* creat an array for output colours */
    if (!(cp->rgbp = (RGB_struct *) malloc(sizeof (RGB_struct)))) {
        fprintf(stderr, "ERROR: RGB_struct malloc failed\n");
        return (1);
    }
    fclose(fp);
    free(line);
    free(sep);
    return rflag; /* return the number of errors whilst reading file */
} /* end of subroutine read_ctrl */

/*****************************************************************/
/* subroutine to combine the logical values to create some useful flags */

/*****************************************************************/

int combine_ctrl(Ctrl_str * cp) {

    cp->use_global_undercooling = (cp->phase_diag_on && cp->global_undercooling);
    cp->use_csol_alloy = (cp->phase_diag_on || cp->particle);
#ifdef CELL_NUC_OFF
    cp->use_cell_nuc = ((cp->diffuse_alloy_poly || cp->con_cast || cp->phase_diag_on || cp->fgrid_input) /*&& !(cp->fixed_nuc) */);
#else
    cp->use_cell_nuc = 1;
#endif

    cp->use_cell_temp = (cp->con_cast || cp->fgrid_input);
    if (cp->diffuse_alloy_poly == 0) {
        cp->ele_1 = 1; /* force correct value for binary-only case */
        cp->NUM_PHS = 1; /* pure binary/oversimplified eutectic model */
        cp->NUM_EQ = 1; /* gibbs phase rule */
        cp->NUM_NUC_LAW = 1; /* default case -- not nucleating eutectic */
    } else {
        cp->ele_1 = cp->NUM_COMP - 1; /* for polycomponent THUINET model */
    }
    /** \todo  find incorrect use of con_cast and replace with use_cell_temp -- general - easy */
    return (0);
}

/*******************************************************/
/*    Entry point for reading the control file        **/
/*    Set the defaults, then read the values          **/
/*    On restart, the values are read without         **/
/*    setting the defaults -- after the restart       **/
/*    values are read in.                             **/

/*******************************************************/
int read_ctrl(char *filename, Ctrl_str * cp) {
    int errors = 0;

    set_ctrl_defaults(cp);
    errors = read_ctrl_vals(filename, cp);
    return (errors);
}

/********************************************************/
/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file*/

/********************************************************/
char const *rcs_id_read_ctrl_c() {
    static char const rcsid[] = "$Id: read_ctrl.c 1382 2008-09-24 14:59:34Z  $";

    return (rcsid);
}

/* end of rcs_id_read_ctrl_c subroutine */

