#Cmake project file for uMat IC microstructure simulation software 
#By Robert Atwood
#$Id: CMakeLists.txt 1405 2008-12-04 13:47:23Z  $
PROJECT (uMatIC)
cmake_minimum_required(VERSION 2.4)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Ensure that the source directory is referenced for #include files
include_directories(${uMatIC_BINARY_DIR})
include_directories(${uMatIC_SOURCE_DIR})
# Give the location of additional project files
# Functions for reading the external data files
add_subdirectory(./read_list)
link_directories(./read_list)
#add the default c flags

#Create the custom options. 
#Using UMAT prefix ensures that options are grouped together (alphabetically)
#and clarifies which options are custom for this project.

#If NEITHER of COMPRESS_BLOCKS or IMAGE_OUTPUT is selected then the compression 
#library zlib is not included. 
option(UMAT_COMPRESS_BLOCKS " Compile with block compression" ON)
if(UMAT_COMPRESS_BLOCKS)
   set(COMPRESS_DEFS -DBL_COMPRESS)
endif(UMAT_COMPRESS_BLOCKS)

option(UMAT_32 " Compile a 32 bit version using -m32 flag" OFF)
   if(UMAT_32)
   add_definitions(-m32)
   set(CMAKE_EXE_LINKER_FLAGS -m32)
   set(CMAKE_SHARED_LINKER_FLAGS -m32)
   endif(UMAT_32)

option(UMAT_LIST_ALL_NUC " List all nuclei instead of just a summary" OFF)
if (UMAT_LIST_ALL_NUC)
   add_definitions(-DLIST_ALL_NUC)
endif(UMAT_LIST_ALL_NUC)

option(UMAT_IMAGE_OUTPUT " Compile with png image output" ON)
if(UMAT_IMAGE_OUTPUT )
   set(IMAGE_DEFS -DIMAGEOUT)
endif(UMAT_IMAGE_OUTPUT )

option(UMAT_BUILD_READER " Compile umat_read utility" OFF)
option(UMAT_BUILD_MAIN " Compile the main program" ON)
option(UMAT_BUILD_UTIL " Compile additional utilities"  OFF)
option(UMAT_EXTERNAL_LIB " Compile shared library for external macromodel" OFF)
option(UMAT_EXTERNAL_STUB " Debugging stub subroutines for external macromodel" OFF)
option(UMAT_DEBUG_MANYDUMPS " set the MANYDUMPS flag (for extra output)" OFF)
option(ELECTRIC_FENCE "use electric fence for debugging malloc" OFF)

if(UMAT_EXTERNAL_LIB)
   set(UMAT_BUILD_MAIN OFF)
   set(UMAT_BUILD_READER OFF)
   set(PCFILES  user_rtns.c
           user_rtn_cell.c
           procast_userfuncs.c
           ca_feedback.c
   )
   set(N_PRO_LCDEFS 
      -DEXTERNAL_HEATFLUX 
      -DEXTERNAL_CA
      -DEXTERNAL_3D
      -DUSER_2
      -DMYALLOC
#      -DALL_SIGS
#      -DVERBOSE_EXTERNAL 
#      -DDMALLOC 
#      -DNOFS_EXTERNAL  
#      -DQ_FEEDBACK_ROUTINE 
#      -DPROBE_CELL=1275 
   )
   if (UMAT_EXTERNAL_STUB)
       set(PRO_LCDEFS ${N_PRO_LCDEFS} -DEXTERNAL_STUBS)
   else (UMAT_EXTERNAL_STUB)
       set(PRO_LCDEFS ${N_PRO_LCDEFS})
   endif (UMAT_EXTERNAL_STUB)
endif(UMAT_EXTERNAL_LIB)

set(REG_LCDEFS  -DLIST_ALL_NUC 
   -DPRINT_AV_PROFILE
   -DGLOBAL_UND_GRO
   -DAV_DIFF_COEFF
   -DNOSOLID
   -DCELL_DIFF_ARR
   -DHISTO_NUC
   -DALWAYS_NUC
   -DRANSLEY
   -DTEST_LINE_INT
   -DWRITE_BLOCK
   -DINCONEL
   -DNLM_FORMULA
   -DCA_DOUBLE
   -DCLOCK
   -DPNG_OUTPUT
   -fPIC  CACHE INTERNAL "internal" FORCE )



#list of all the c files for this project
#this will need to use slightly more complicated syntax to define
#different projects e.g. ca_procast, sharing many of the same files.

set(MAINFILE ca_wrapper.c)
set(RMAINFILE ca_read_wrapper.c)


set(CFILES
   alloc_bb.c
   alloc_multicomp.c
   alloc_sb.c
   bigblock.c
   ca_histo.c
   calc_sb.c
   ca_solid.c
   castats.c
   curvature.c
   curvature_3D.c
   checkblock.c
   checkgas.c
   checks.c
   close_sb.c
   fcopymat.c
   fg_read.c
   fidap_interp_calc.c
   find_max.c
   find_nmols.c
   findroot.c
   freeblock.c
   gaussdev.c
   getxyz.c
   grain_functions.c
   handlers.c
   icopymat.c
   initcube.c
   initface.c
   init_sb_neigh.c
   init_nuc_mould.c
   mould_sources.c
   nbhd_def.c
   nuc_functions.c
   nuc_lookup.c
   open_sb.c
   output_ex.c
   output_img.c
   p_growth.c
   pnuc_function.c
   pore_multi.c
   pore_routines.c
   pore_setup_multi.c
   pr_lookup.c
   props.c
   rand_square.c
   rand_step.c
   rcs_id.c
   read_blocks.c
   read_cap_ca.c
   read_ca_procast.c
   read_ctrl.c
   readgeoplus.c
   readmat.c
   readphase.c
   read_sb.c
   recr_extra.c
   safeopen.c
   sb_boundary.c
   sb_ca_step.c
   sb_diffuse_alloy.c
   sb_diffuse_gas.c
   sb_diffuse_gas_poly.c
   sb_get_surface_cells.c
   sb_line_int.c
   sb_nuc.c
   sb_temp_calc.c
   setup_mould_src_function.c
   setup_temp_func.c
   solprops_reader.c
   solprops_writer.c
   step_output.c
   subblock.c
   SurCellRoutines.c
   thermo_trace_calc.c
   temp_calc_melt.c
   thermo_trace_init.c
   trans_interp_calc.c
   wfact_r_calc.c
   wfact_z_calc.c
   write_bb_values.c
   write_blocks.c
   write_block_values.c
   write_ctrl_values.c
   write_grain_histo.c
   write_header_ctrl.c
   write_matprop_values.c
   write_nprop_values.c
   write_POREprop_values.c
   write_pprop_values.c
   calc_solprops_poly.c
   sb_decentred_step.c
   sb_diffuse_alloy_decentred.c
   window_move.c
   sb_decentred_step_poly.c
   sb_diffuse_alloy_decentred_poly.c
   replacerand.c
   olderror.c
)

#conditionally include the library files needed for uMat_run, without the 'lib' prefix
# i.e. libz, libm, libpng, libgd and the library produced by the above subprojects

#Choose the system libraries needed depending upon the options
#We always need math library (but it may have a different name on some other systems)
set(SLIBS m )

#gather all the options into a single list variable
set(ALL_DEFS ${COMPRESS_DEFS} ${IMAGE_DEFS} ${REG_LCDEFS} $ENV{CFLAGS}  )

add_definitions(${ALL_DEFS})
# set the MANYDUMPS flag (for extra output)
if(UMAT_DEBUG_MANYDUMPS)
   add_definitions(-DMANYDUMPS)
endif(UMAT_DEBUG_MANYDUMPS)


if(ELECTRIC_FENCE)
   set(CMAKE_EXE_LINKER_FLAGS_DEBUG -lefence)
endif(ELECTRIC_FENCE)


   if(UMAT_COMPRESS_BLOCKS OR UMAT_IMAGE_OUTPUT)
#png depends upon zlib, gd depends upon png
	set(CLIBS  z )
   endif(UMAT_COMPRESS_BLOCKS OR UMAT_IMAGE_OUTPUT)
   if(UMAT_IMAGE_OUTPUT)
      set(IMLIBS png gd)
   endif(UMAT_IMAGE_OUTPUT)
   set(SLIBS ${SLIBS} ${CLIBS} ${IMLIBS} )
#end of library setting section

   if (UMAT_EXTERNAL_LIB)
      add_library(ca_umat SHARED ${PCFILES})   
      ADD_DEPENDENCIES(ca_umat includer)
      target_link_libraries(ca_umat ${SLIBS} readlist umatic)   
      add_definitions(${PRO_LCDEFS})
   endif(UMAT_EXTERNAL_LIB)

#ensure that the generated readlist library is linked in the project
   if(UMAT_BUILD_MAIN OR UMAT_BUILD_READER OR UMAT_EXTERNAL_LIB)
      add_library(umatic ${CFILES})
      ADD_DEPENDENCIES(umatic includer)
      target_link_libraries(umatic readlist)
   endif(UMAT_BUILD_MAIN OR UMAT_BUILD_READER OR UMAT_EXTERNAL_LIB)


   if (UMAT_BUILD_MAIN)
      add_executable(umat_run ${MAINFILE} )
      ADD_DEPENDENCIES(umat_run includer)
      target_link_libraries(umat_run ${SLIBS} readlist umatic)
   endif(UMAT_BUILD_MAIN)

   if(UMAT_BUILD_READER)
       add_executable(umat_read ${RMAINFILE} )
       ADD_DEPENDENCIES(umat_read includer)
       target_link_libraries(umat_read ${SLIBS} readlist umatic)
   endif(UMAT_BUILD_READER)

   if (UMAT_BUILD_UTIL)
      add_executable(read_list_util safeopen.c read_list/read_list.c)
      target_link_libraries(read_list_util readlist)
   endif(UMAT_BUILD_UTIL)




#create an include-file for incorporating the build details in the executable
#this should also somehow include a dependency for getcflags.inc ? 
#
set(INCFILE getcflags.h)

ADD_CUSTOM_TARGET(includer
         COMMAND echo  '\#define GETCFLAGS \(\"\\' > ${INCFILE}
         COMMAND echo  'Build date: \\' >> ${INCFILE}
	 COMMAND echo `date` '\\' >> ${INCFILE}
         COMMAND echo  '\\nConfigure Username:\\t'$ENV{USER}'\\' >> ${INCFILE}
         COMMAND echo  '\\nConfigure Host:\\t'$ENV{HOSTNAME}'\\' >> ${INCFILE}
         COMMAND echo  '\\nCompiler Info:\\t' `${CMAKE_C_COMPILER} --version | head -1` '\\' >> ${INCFILE}
         COMMAND echo  '\\nConfigure CFLAGS:\\t'${ALL_DEFS}'\\' >> ${INCFILE}
         COMMAND echo  '\\nSVN Version:\\t' `svnversion ${uMatIC_SOURCE_DIR}` '\\'  >> ${INCFILE}
         COMMAND echo  '\\n\"\)' >> ${INCFILE}
         COMMAND echo  '\#define GETSVNINFO \(\"\\' >> ${INCFILE}
         COMMAND echo  '\\nsvn status:\\n' `cd ${uMatIC_SOURCE_DIR} && svn status  | while read a\; do echo  $$a'\\n'\;done` '\\'  >> ${INCFILE}
         COMMAND echo  '\\nsvn info:\\n' `cd ${uMatIC_SOURCE_DIR} && svn info | while read a\; do echo $$a'\\n'\;done` '\\' >> ${INCFILE}
         COMMAND echo  '\\n\"\)' >> ${INCFILE}
)

MESSAGE("Info: The CFLAGS environment variable is listed on the following
      line. If the line is blank then no CFLAGS are set.\nCFLAGS= $ENV{CFLAGS}\n  Press 'e' to continue.")
