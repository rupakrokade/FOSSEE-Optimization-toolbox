#ifdef __cplusplus
extern "C" {
#endif
#include <mex.h> 
#include <sci_gateway.h>
#include <api_scilab.h>
#include <MALLOC.h>
static int direct_gateway(char *fname,void F(void)) { F();return 0;};
extern Gatefunc sci_sym_open;
extern Gatefunc sci_sym_close;
extern Gatefunc sci_sym_set_defaults;
extern Gatefunc sci_sym_set_int_param;
extern Gatefunc sci_sym_get_int_param;
extern Gatefunc sci_sym_set_dbl_param;
extern Gatefunc sci_sym_get_dbl_param;
extern Gatefunc sci_sym_set_str_param;
extern Gatefunc sci_sym_get_str_param;
extern Gatefunc sci_sym_loadProblemBasic;
extern Gatefunc sci_sym_loadProblem;
extern Gatefunc sci_sym_solve;
extern Gatefunc sci_sym_get_status;
extern Gatefunc sci_sym_get_solver_status;
extern Gatefunc sci_sym_get_solver_status;
extern Gatefunc sci_sym_get_solver_status;
extern Gatefunc sci_sym_get_solver_status;
extern Gatefunc sci_sym_get_solver_status;
extern Gatefunc sci_sym_get_solver_status;
extern Gatefunc sci_sym_getVarSoln;
extern Gatefunc sci_sym_getObjVal;
extern Gatefunc sci_sym_get_iteration_count;
extern Gatefunc sci_linearprog;
extern Gatefunc sci_rmps;
extern Gatefunc sci_solveqp;
extern Gatefunc sci_solveminuncp;
extern Gatefunc sci_solveminbndp;
extern Gatefunc sci_solveminconp;
extern Gatefunc cpp_intfminunc;
extern Gatefunc cpp_intfminbnd;
extern Gatefunc cpp_intfmincon;
extern Gatefunc cpp_intqpipopt;
extern Gatefunc matrix_cppintlinprog;
extern Gatefunc mps_cppintlinprog;
extern Gatefunc sci_ecos;
extern Gatefunc sci_fotversion;
static GenericTable Tab[]={
  {(Myinterfun)sci_gateway,sci_sym_open,"sym_open"},
  {(Myinterfun)sci_gateway,sci_sym_close,"sym_close"},
  {(Myinterfun)sci_gateway,sci_sym_set_defaults,"sym_resetParams"},
  {(Myinterfun)sci_gateway,sci_sym_set_int_param,"sym_setIntParam"},
  {(Myinterfun)sci_gateway,sci_sym_get_int_param,"sym_getIntParam"},
  {(Myinterfun)sci_gateway,sci_sym_set_dbl_param,"sym_setDblParam"},
  {(Myinterfun)sci_gateway,sci_sym_get_dbl_param,"sym_getDblParam"},
  {(Myinterfun)sci_gateway,sci_sym_set_str_param,"sym_setStrParam"},
  {(Myinterfun)sci_gateway,sci_sym_get_str_param,"sym_getStrParam"},
  {(Myinterfun)sci_gateway,sci_sym_loadProblemBasic,"sym_loadProblemBasic"},
  {(Myinterfun)sci_gateway,sci_sym_loadProblem,"sym_loadProblem"},
  {(Myinterfun)sci_gateway,sci_sym_solve,"sym_solve"},
  {(Myinterfun)sci_gateway,sci_sym_get_status,"sym_getStatus"},
  {(Myinterfun)sci_gateway,sci_sym_get_solver_status,"sym_isOptimal"},
  {(Myinterfun)sci_gateway,sci_sym_get_solver_status,"sym_isInfeasible"},
  {(Myinterfun)sci_gateway,sci_sym_get_solver_status,"sym_isAbandoned"},
  {(Myinterfun)sci_gateway,sci_sym_get_solver_status,"sym_isIterLimitReached"},
  {(Myinterfun)sci_gateway,sci_sym_get_solver_status,"sym_isTimeLimitReached"},
  {(Myinterfun)sci_gateway,sci_sym_get_solver_status,"sym_isTargetGapAchieved"},
  {(Myinterfun)sci_gateway,sci_sym_getVarSoln,"sym_getVarSoln"},
  {(Myinterfun)sci_gateway,sci_sym_getObjVal,"sym_getObjVal"},
  {(Myinterfun)sci_gateway,sci_sym_get_iteration_count,"sym_getIterCount"},
  {(Myinterfun)sci_gateway,sci_linearprog,"linearprog"},
  {(Myinterfun)sci_gateway,sci_rmps,"rmps"},
  {(Myinterfun)sci_gateway,sci_solveqp,"solveqp"},
  {(Myinterfun)sci_gateway,sci_solveminuncp,"solveminuncp"},
  {(Myinterfun)sci_gateway,sci_solveminbndp,"solveminbndp"},
  {(Myinterfun)sci_gateway,sci_solveminconp,"solveminconp"},
  {(Myinterfun)sci_gateway,cpp_intfminunc,"inter_fminunc"},
  {(Myinterfun)sci_gateway,cpp_intfminbnd,"inter_fminbnd"},
  {(Myinterfun)sci_gateway,cpp_intfmincon,"inter_fmincon"},
  {(Myinterfun)sci_gateway,cpp_intqpipopt,"sci_intqpipopt"},
  {(Myinterfun)sci_gateway,matrix_cppintlinprog,"sci_matrix_intlinprog"},
  {(Myinterfun)sci_gateway,mps_cppintlinprog,"sci_mps_intlinprog"},
  {(Myinterfun)sci_gateway,sci_ecos,"solveecos"},
  {(Myinterfun)sci_gateway,sci_fotversion,"fotversion"},
};
 
int C2F(libFOSSEE_Optimization_Toolbox)()
{
  Rhs = Max(0, Rhs);
  if (*(Tab[Fin-1].f) != NULL) 
  {
     if(pvApiCtx == NULL)
     {
       pvApiCtx = (StrCtx*)MALLOC(sizeof(StrCtx));
     }
     pvApiCtx->pstName = (char*)Tab[Fin-1].name;
    (*(Tab[Fin-1].f))(Tab[Fin-1].name,Tab[Fin-1].F);
  }
  return 0;
}
#ifdef __cplusplus
}
#endif
