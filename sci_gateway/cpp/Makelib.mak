# ------------------------------------------------------
# generated by builder.sce : Please do not edit this file
# see TEMPLATE makefile for Visual Studio
# see SCI/modules/dynamic_link/src/scripts/TEMPLATE_MAKEFILE.VC
# ------------------------------------------------------
SCIDIR = C:/PROGRA~1/SCILAB~1.2
# ------------------------------------------------------
# default include options 
INCLUDES = -I"$(SCIDIR)/libs/MALLOC/includes" \
-I"$(SCIDIR)/libs/f2c" \
-I"$(SCIDIR)/libs/hashtable" \
-I"$(SCIDIR)/libs/intl" \
-I"$(SCIDIR)/modules/core/includes" \
-I"$(SCIDIR)/modules/api_scilab/includes" \
-I"$(SCIDIR)/modules/call_scilab/includes" \
-I"$(SCIDIR)/modules/output_stream/includes" \
-I"$(SCIDIR)/modules/jvm/includes" \
-I"$(SCIDIR)/modules/localization/includes" \
-I"$(SCIDIR)/modules/dynamic_link/includes" \
-I"$(SCIDIR)/modules/mexlib/includes" \
-I"$(SCIDIR)/modules/time/includes" \
-I"$(SCIDIR)/modules/windows_tools/includes"
# ------------------------------------------------------
# SCILAB_LIBS is used by the binary version of Scilab for linking external codes
SCILAB_LIBS = "$(SCIDIR)/bin/blasplus.lib" \
"$(SCIDIR)/bin/libf2c.lib" \
"$(SCIDIR)/bin/core.lib" \
"$(SCIDIR)/bin/core_f.lib" \
"$(SCIDIR)/bin/lapack.lib" \
"$(SCIDIR)/bin/libintl.lib" \
"$(SCIDIR)/bin/intersci.lib" \
"$(SCIDIR)/bin/output_stream.lib" \
"$(SCIDIR)/bin/dynamic_link.lib" \
"$(SCIDIR)/bin/integer.lib" \
"$(SCIDIR)/bin/optimization_f.lib" \
"$(SCIDIR)/bin/libjvm.lib" \
"$(SCIDIR)/bin/scilocalization.lib" \
"$(SCIDIR)/bin/linpack_f.lib" \
"$(SCIDIR)/bin/call_scilab.lib" \
"$(SCIDIR)/bin/time.lib" \
"$(SCIDIR)/bin/api_scilab.lib" \
"$(SCIDIR)/bin/libintl.lib" \
"$(SCIDIR)/bin/scilab_windows.lib"
# ------------------------------------------------------
# name of the dll to be built
LIBRARY = FOSSEE_Optimization_Toolbox
# ------------------------------------------------------
# list of files
FILES_SRC = FOSSEE_Optimization_Toolbox.c globals.cpp read_mps.cpp sci_LinProg.cpp sci_QuadNLP.cpp sci_iofunc.cpp sci_ipopt.cpp sci_ipoptfminbnd.cpp sci_ipoptfmincon.cpp sci_ipoptfminunc.cpp sci_minbndNLP.cpp sci_minconNLP.cpp sci_minuncNLP.cpp sci_solver_status_query_functions.cpp sci_sym_addrowcol.cpp sci_sym_get_dbl_arr.cpp sci_sym_get_iteration_count.cpp sci_sym_get_matrix.cpp sci_sym_get_num_int.cpp sci_sym_getinfinity.cpp sci_sym_getobjsense.cpp sci_sym_getrowact.cpp sci_sym_isenvactive.cpp sci_sym_load_mps.cpp sci_sym_loadproblem.cpp sci_sym_openclose.cpp sci_sym_primalbound.cpp sci_sym_remove.cpp sci_sym_rowmod.cpp sci_sym_set_indices.cpp sci_sym_set_variables.cpp sci_sym_setcolsoln.cpp sci_sym_setobj.cpp sci_sym_solution.cpp sci_sym_solve.cpp sci_sym_varbounds.cpp sci_vartype.cpp
# ------------------------------------------------------
# list of objects file
OBJS = FOSSEE_Optimization_Toolbox.obj globals.obj read_mps.obj sci_LinProg.obj sci_QuadNLP.obj sci_iofunc.obj sci_ipopt.obj sci_ipoptfminbnd.obj sci_ipoptfmincon.obj sci_ipoptfminunc.obj sci_minbndNLP.obj sci_minconNLP.obj sci_minuncNLP.obj sci_solver_status_query_functions.obj sci_sym_addrowcol.obj sci_sym_get_dbl_arr.obj sci_sym_get_iteration_count.obj sci_sym_get_matrix.obj sci_sym_get_num_int.obj sci_sym_getinfinity.obj sci_sym_getobjsense.obj sci_sym_getrowact.obj sci_sym_isenvactive.obj sci_sym_load_mps.obj sci_sym_loadproblem.obj sci_sym_openclose.obj sci_sym_primalbound.obj sci_sym_remove.obj sci_sym_rowmod.obj sci_sym_set_indices.obj sci_sym_set_variables.obj sci_sym_setcolsoln.obj sci_sym_setobj.obj sci_sym_solution.obj sci_sym_solve.obj sci_sym_varbounds.obj sci_vartype.obj
OBJS_WITH_PATH = Release/FOSSEE_Optimization_Toolbox.obj Release/globals.obj Release/read_mps.obj Release/sci_LinProg.obj Release/sci_QuadNLP.obj Release/sci_iofunc.obj Release/sci_ipopt.obj Release/sci_ipoptfminbnd.obj Release/sci_ipoptfmincon.obj Release/sci_ipoptfminunc.obj Release/sci_minbndNLP.obj Release/sci_minconNLP.obj Release/sci_minuncNLP.obj Release/sci_solver_status_query_functions.obj Release/sci_sym_addrowcol.obj Release/sci_sym_get_dbl_arr.obj Release/sci_sym_get_iteration_count.obj Release/sci_sym_get_matrix.obj Release/sci_sym_get_num_int.obj Release/sci_sym_getinfinity.obj Release/sci_sym_getobjsense.obj Release/sci_sym_getrowact.obj Release/sci_sym_isenvactive.obj Release/sci_sym_load_mps.obj Release/sci_sym_loadproblem.obj Release/sci_sym_openclose.obj Release/sci_sym_primalbound.obj Release/sci_sym_remove.obj Release/sci_sym_rowmod.obj Release/sci_sym_set_indices.obj Release/sci_sym_set_variables.obj Release/sci_sym_setcolsoln.obj Release/sci_sym_setobj.obj Release/sci_sym_solution.obj Release/sci_sym_solve.obj Release/sci_sym_varbounds.obj Release/sci_vartype.obj
# ------------------------------------------------------
# added libraries
FORTRAN_RUNTIME_LIBRARIES = 
OTHERLIBS = 
# ------------------------------------------------------
!include $(SCIDIR)\modules\dynamic_link\src\scripts\Makefile.incl.mak
# ------------------------------------------------------
#CC = 
# ------------------------------------------------------
CFLAGS = $(CC_OPTIONS) -D__SCILAB_TOOLBOX__ -DFORDLL -D__USE_DEPRECATED_STACK_FUNCTIONS__ -w -I C:\Users\FOSSEE\Desktop\build\FOSSEE-Optimization-Toolbox-src\sci_gateway\cpp\ -I C:\Users\FOSSEE\Desktop\build\FOSSEE-Optimization-Toolbox-src\sci_gateway\cpp\\..\..\thirdparty\windows\include\coin   
# ------------------------------------------------------
FFLAGS = $(FC_OPTIONS) -DFORDLL  
# ------------------------------------------------------
EXTRA_LDFLAGS = C:\Users\FOSSEE\Desktop\build\FOSSEE-Optimization-Toolbox-src\sci_gateway\cpp\\..\..\thirdparty\windows\lib\x86\libClp.lib C:\Users\FOSSEE\Desktop\build\FOSSEE-Optimization-Toolbox-src\sci_gateway\cpp\\..\..\thirdparty\windows\lib\x86\libCgl.lib C:\Users\FOSSEE\Desktop\build\FOSSEE-Optimization-Toolbox-src\sci_gateway\cpp\\..\..\thirdparty\windows\lib\x86\libOsi.lib C:\Users\FOSSEE\Desktop\build\FOSSEE-Optimization-Toolbox-src\sci_gateway\cpp\\..\..\thirdparty\windows\lib\x86\libOsiClp.lib C:\Users\FOSSEE\Desktop\build\FOSSEE-Optimization-Toolbox-src\sci_gateway\cpp\\..\..\thirdparty\windows\lib\x86\libCoinUtils.lib C:\Users\FOSSEE\Desktop\build\FOSSEE-Optimization-Toolbox-src\sci_gateway\cpp\\..\..\thirdparty\windows\lib\x86\libSymphony.lib C:\Users\FOSSEE\Desktop\build\FOSSEE-Optimization-Toolbox-src\sci_gateway\cpp\\..\..\thirdparty\windows\lib\x86\IpOptFSS.lib C:\Users\FOSSEE\Desktop\build\FOSSEE-Optimization-Toolbox-src\sci_gateway\cpp\\..\..\thirdparty\windows\lib\x86\IpOpt-vc10.lib 
# ------------------------------------------------------
!include $(SCIDIR)\modules\dynamic_link\src\scripts\Makedll.incl
# ------------------------------------------------------
