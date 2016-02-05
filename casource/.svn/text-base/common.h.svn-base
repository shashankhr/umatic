#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>

#define M_ENTH_P       200
#define M_TEMP_P       200
#define M_QUAD_P       10
#define M_TIME_P       200
#define M_PRES_P       50
#define M_COOL_P       50
#define M_CONC_P       50

#define U 0
#define V 1
#define W 2
#define P 3
#define T 4
#define F 5
#define K 6
#define E 7
#define H 8
#define D 9
#define X 10
#define Y 11
#define Z 12
#define DD 13
#define Q 14
#define A 15
#define C 16
#define PHI_E 17
#define FF 18
#define REMOVE 19

extern int dim;
extern int nel;
extern int ntel;
extern int nvel;
extern int nnod;
extern int ntnod;
extern int nvnod;
extern int n_link;
extern int n_link_mpc;
extern int n_neu_data;
extern int n_coin_face;
extern int n_coin_data;
extern int n_dir_data;
extern int n_dir_info;
extern int n_encl_face;
extern int n_encl_node;
extern int n_encl_data;
extern int n_vel_data;
extern int n_sym;
extern int nradl;
extern int mirror;
extern int n_neu_face;
extern int n_total_face;
extern int n_mpc;
extern int n_id;
extern int n_mat;
extern int n_foam;
extern int n_enth_f;
extern int n_enth_p;
extern int n_temp_f;
extern int n_temp_p;
extern int n_quad_f;
extern int n_quad_p;
extern int n_time_f;
extern int n_time_p;
extern int n_pres_f;
extern int n_pres_p;
extern int n_cool_f;
extern int n_cool_p;
extern int n_conc_f;
extern int n_conc_p;
extern int n_fs_f;
extern int n_fs_p;
extern int n_q_vol_data;
extern int l_rad;
extern int NLINEAR;
extern int istep;
extern int iter;
extern int l_per;
extern int n_solid;
extern int n_vents;
extern int n_inject;
extern int n_aniso;
extern int n_sym_face;
extern int n_drag_face;
extern int n_nucl_face;
extern int n_periodic;
extern int n_mic_nod;
extern int n_mom_data;

extern int mcoef;

extern int npe[];
extern int fpe[];
extern int n_nodes[][6];
extern int local_node[][6][6];

extern int *(*ncon);
extern int *node_list;
extern int *el_type;
extern int *mat_id;
extern int *neu_status;
extern int *q_vol_status;
extern int *free_face;
extern int *el_st;
extern int *adv_st;
extern int *wall_elem;
extern float *wall_dist;
extern float *w_dist;
extern float *k_turb;
extern int *zero_dir;
extern int *nonz_dir;
extern int *temp_dir;
extern int *or_dir;
extern int *dir_node;
extern int *dir_index;
extern int *dir_bit;
extern int *ld;
extern int *i_per;
extern int *icoef;
extern int (*i_link)[5];
extern int *vent_node;
extern int *inject_node;
extern int *pvent_time;
extern int *inj_time;
extern int *inj_pres;
extern int *sym_el_face;
extern int *drag_el_face;
extern int *trans_coef;
extern int *inside_node;
extern int *t_form;
extern int *node_use;
extern int *node_use_ptr;
extern int *node_use_freq;
extern int *node_map;
extern int *mom_status;
extern int w_shear;
extern char *no_slip_node;

extern int g_time_func;
extern float grav_vec[6];
extern double gx;
extern double gy;
extern double gz;
extern double grav_mag;
extern int omega_time_fnc;
extern int *omega_node;
extern float *omega_data;
extern int n_rot_node;
extern int n_coin_mpc;

extern double *d_coef;
extern double *pd_coef;
extern double *d_fac_coef;
extern float *ucoef;
extern float *lcoef;
extern float *pcoef;
extern float *fac_coef;
extern float *ufac_coef;
extern float *lfac_coef;
extern double *dd_fac_coef;
extern double *id_coef;
extern float *iucoef;

extern double *rhs_v;
extern double *rhs_x;
extern double *rhs_y;
extern double *rhs_z;
extern double *rhs_s;
extern double *rhs_store;
extern double *dpdx;
extern double *dpdy;
extern double *dpdz;
extern double *diverge;
extern float *el_theta;
extern double *x_cord;
extern double *y_cord;
extern double *z_cord;
extern float *q_rad;
extern float *wdv;
extern float *wdv_frac;
extern float *aii_x;
extern float *aii_y;
extern float *f_link;
extern float sym_axis[3][9];
extern float (*per_data)[10];
extern double (*per_rot)[9];

extern double detj[8];
extern double a11[8];
extern double a12[8];
extern double a13[8];
extern double a21[8];
extern double a22[8];
extern double a23[8];
extern double a31[8];
extern double a32[8];
extern double a33[8];

extern double bshape[8][8];
extern double wshape[6][6];
extern double tshape[4][4];
extern double sshape[5][6][4];
extern double hotshape[10][4];
extern double bderv[3][8][8];
extern double wderv[3][6][6];
extern double tderv[3][4][4];
extern double sderv[5][2][6][4];
extern double hotderv[3][10][4];
extern double t6w[4];
extern double gp_temp[2][8];
extern double gp_fs[2][8];
extern double gp_freeze[8];
extern double gp_conc[8];

extern float (*enth_fnc)[ M_ENTH_P ][2];
extern float (*temp_fnc)[ M_TEMP_P ][2];
extern float (*quad_fnc)[ M_QUAD_P ][4];
extern float (*time_fnc)[ M_TIME_P ][2];
extern float (*pres_fnc)[ M_PRES_P ][2];
extern float (*cool_fnc)[ M_COOL_P ][2];
extern float (*conc_fnc)[ M_CONC_P ][2];
extern float (*fs_fnc)[ M_TEMP_P ][2];
extern int *enth_pts; 
extern int *temp_pts; 
extern int *quad_pts;
extern int *time_pts; 
extern int *pres_pts; 
extern int *cool_pts; 
extern int *conc_pts; 
extern int *fs_pts; 

extern int *enth_ptr;
extern int *i_dens;
extern int *i_cp;
extern int *i_fs;
extern int *i_enth;
extern int *i_cond;
extern int *i_vis;
extern int *i_fl_perm;
extern int *i_st;
extern int *i_nc;
extern int *i_elastic;
extern int *i_pois;
extern int *i_yield;
extern int *i_ultimate;
extern int *i_th_exp;
extern int *i_hard;
extern int *i_hard_exp;
extern int (*i_kin_hard)[2];
extern int *i_vis_pow;
extern int *i_st_vis;

extern int *mat_num;
extern int *solid_vel;
extern int *eq_store;
extern int *mold;
extern int n_filters;
extern int *fluid_state;
extern int *stress_type;
extern int (*i_filter)[2];
extern float (*f_filter)[4];
extern int *filter_nodes;
extern int *filter_num;
extern int f_inject;
extern int *f_inj_node;
extern int *foam_m_id;
extern int *ff_ptr;
extern int surf_bc;

extern float *foam_inj;
extern float *area_ratio;
extern float *density;
extern float *spec_heat;
extern float *cond;
extern float *aniso;
extern float *vis;
extern float *liquidus;
extern float *solidus;
extern float *elastic_mod;
extern float *pois_ratio;
extern float *yield_stress;
extern float *ultimate_stress;
extern float *th_expansion;
extern float *hardening;
extern float *hardening_exp;
extern float (*kin_hard)[2];
extern float *visco_power;
extern float *st_vis;

extern float *mat_theta;
extern float *gas_perm;
extern float *fl_perm;
extern float *surf_tension;
extern float *lat_heat;
extern double rho[8];
extern double cp[8];
extern double conductivity[8];
extern double viscosity[8];
extern double elast_loc[8];
extern double pois_loc[8];
extern double yield_loc[8];
extern double ultimate_loc[8];
extern double th_exp_loc[8][2];
extern double hard_loc[8];
extern double hard_exp_loc[8];
extern double kin_hard_loc[8][2];
extern double vis_pow_loc[8];
extern double st_vis_loc[8];
extern double mag_p[8];
extern double mag_c[8];
extern double mass_diffusivity[8];
extern double perm[8];
extern double frac[8];

extern int (*neu_info)[3];
extern int (*i_neu_data)[6];
extern int (*nucl_info)[2];
extern float (*f_neu_data)[4];
extern int (*i_dir_data);
extern float (*f_dir_data);
extern int (*i_dir_info)[2];
extern float *f_dir_info;
extern int (*i_non_newt)[5];
extern float (*f_non_newt)[5];
extern int (*i_mic_data)[6];
extern float (*f_mic_data)[6];
extern int (*mic_ptr)[20];
extern int *i_mic_array;
extern float *f_mic_array;
extern int (*i_coin_mpc)[4];
extern float (*f_coin_mpc)[3];
extern int *i_mom_data;
extern float (*f_mom_data)[3];
extern int *cover_up;

extern int *coin_info;
extern int (*i_coin_data)[2];
extern float *f_coin_data;
extern int (*i_q_vol_data)[2];
extern float *f_q_vol_data;

extern int (*i_rad)[6];
extern int (*k_rad)[4];
extern float (*f_rad)[2];
extern int *j_rad;
extern int *s_rad;
extern int *i_encl_vel;
extern float (*f_encl_vel)[3];
extern int (*i_encl_prop)[2];
extern float (*f_encl_prop)[2];
extern int (*i_mpc)[5];
extern float (*f_mpc)[4];
extern int *vel_node;
extern int *elf_mat;
extern int *node_mat;

extern float *t0, *t1;
extern float *fs0, *fs1;
extern float *u0, *u1;
extern float *v0, *v1;
extern float *w0, *w1;
extern float *p0, *p1;
extern float *k0, *k1;
extern float *e0, *e1;
extern float *h0, *h1;
extern float *rho0, *rho1;
extern float *vf0, *vf1, *vf2;
extern float *bubble_radius;
extern float *start_freezing;
extern float *axr, *axi;
extern float *ayr, *ayi;
extern float *azr, *azi;

extern float *dtdh;
extern float *mu_t, *mu_t0;
extern float *el_vol;
extern float *rho_bar;
extern float *t_start;
extern float *f_vol0, *f_vol1, *f_recent;
extern float *k_fac_x;
extern float *k_fac_y;
extern float *k_fac_z;
extern double *u_hat;
extern double *v_hat;
extern double *w_hat;
extern double *u_hat_store;
extern double *v_hat_store;
extern double *w_hat_store;
extern float *t_old;
extern float *u_old;
extern float *v_old;
extern float *w_old;
extern float *u_varia;
extern float *v_varia;
extern float *w_varia;
extern float *t_varia;
extern float *dp_tension;
extern double *doo;
extern float *nn_visc;
extern float *nn_shear;
extern float *fse, *fsd, *fsp, *fsm;
extern float *dqdt;

extern float *p_vent;
extern float *vent_dia;
extern float *vent_roughness;
extern float *vent_length;
extern float *m_inject;

extern float *fic_temp1;
extern float *fic_temp0;
extern float *fic_rho0;
extern float *fic_rho1;
extern float *fic_mu0;
extern float *fic_mu1;

extern int n_source;
extern int (*i_source_data)[5];
extern float (*f_source_data)[5];
extern int *source_node;

extern int (*i_vm_data)[5];
extern float *f_vm_data;

extern int n_magnetic;
extern int *i_mp;
extern float *mag_perm;
extern int *i_mc;
extern float *mag_cond;
extern int *i_current;
extern float *current_density;
extern float current_freq;

extern int n_species;
extern int i_sp_off;
extern int i_specie;
extern char (*species_name)[31];
extern int *i_mass_diff;
extern float *f_mass_diff;
extern float *conc1, *conc0;
extern int n_cdir_data;
extern int *cdir_specie;
extern int *cdir_node;
extern int *cdir_temp;
extern float *cdir_value;
extern int *cdir_ptr;

extern struct md { float val; int node; } max_del[ A + 6 ];
extern int solver_type[ A + 6 ];
extern int solver_iter[ A + 6 ];
extern double solver_del[ A + 6 ];

extern char prefix[60];
extern char p_out[60];
extern char d_comments[3][72];
extern char p_comments[3][72];
extern int tunits;
extern int vunits;
extern int punits;
extern int qunits;
extern int prnlev;
extern int conlev;
extern int debug;
extern int rdebug;
extern int sdebug;
extern int nstep;
extern int maxcor;
extern int ncoru;
extern int ncorl;
extern int nrstar;
extern int inilev;
extern int nprfr;
extern int encl_id;
extern int cover;
extern int start_elem;
extern int edge;
extern int n_cycle;
extern int cycle_num;
extern int lufac;
extern int ffreq;
extern int pp_node;
extern int vv_node;
extern int switch_off;
extern int UVW;
#ifdef NT4
extern unsigned long int offset;
#else
extern long long int offset;
#endif

extern float t_cycle;
extern float p_relax;
extern float t_relax;
extern float c_relax;
extern float m_relax;
extern float tb_relax;
extern float conv_crit[5];
extern float convt;
extern float convv;
extern float dt;
extern float dtmax;
extern float t_final;
extern float vf_time;
extern float vf_disp;
extern float beta;
extern float tmods;
extern float tmodr;
extern float mlump;  
extern float clump;  
extern float accel[3];
extern float current_time;
extern float total_time;
extern float dt_old;
extern float stefan;
extern float courant;
extern float lev_surf;
extern float adv_weight;
extern float p_ref;
extern double inlet_flow;
extern double max_velocity;
extern double c_lim2;
extern double datum;
extern double percent_filled;
extern double di_mass;
extern double v_rel[3];
extern double new_mass;
extern double old_mass;
extern int mass_imbalance;
extern int i_rel;

extern float c_mu;
extern float sigma_k;
extern float sigma_e;
extern float c_one;
extern float c_two;
extern float kappa;
extern float intense;
extern float char_len;

extern int FLOW;
extern int THERMAL;
extern int TWO_D;
extern int AXISYM;
extern int FREE_SURFACE;
extern int TURB;
extern int COMPRESS;
extern int STRESS;
extern int MICRO;
extern int EM;
extern int POROSITY;
extern int md_ptr;
extern int i_pool;
extern int f_time;
extern int hi_visc;
extern int coupled;
extern int aveprop;
extern int g_perm;
extern int GAS;
extern int GRAVITY;
extern int SHRINK_IT;
extern int turn_off;
extern int surf_ten;
extern int n_fic;
extern int rad_update;
extern int non_newtonian;
extern int split;
extern int rel_vel;
extern float vflim;
extern float eptol;
extern float angtol;
extern float p_limit;
extern float conv_tol;
extern float flow3_delay;
extern float mobile;

extern char data_vers[6];
extern char data_time[20];
extern char data_date[20];

extern int INT;
extern int FLOAT;
extern int DOUBLE;

#ifndef MYALLOC
extern int *int_alloc();
extern long *long_int_alloc();
extern float *float_alloc();
extern double *double_alloc();
extern char *char_alloc();
extern int *int_realloc();
extern float *float_realloc();
extern double *double_realloc();
extern char *char_realloc();
extern void Free();
#endif

extern int n_crnt_data;
extern int n_crnt_face;
extern int *crnt_status;
extern int (*crnt_info)[3];

extern float *phie;
extern float *bxr, *bxi, *byr, *byi, *bzr, *bzi;

extern int i_frequency;
extern float current_frequency;

extern float *vf_s1,*vf_s0; /* volume fraction of porosity (shrinkage) */

