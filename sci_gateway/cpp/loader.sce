// This file is released under the 3-clause BSD license. See COPYING-BSD.
// Generated by builder.sce : Please, do not edit this file
// ----------------------------------------------------------------------------
//
libsymphonytools_path = get_absolute_file_path('loader.sce');
//
// ulink previous function with same name
[bOK, ilib] = c_link('libsymphonytools');
if bOK then
  ulink(ilib);
end
//
list_functions = [ 'sym_open';
                   'sym_close';
                   'sym_isEnvActive';
                   'sym_resetParams';
                   'sym_setIntParam';
                   'sym_getIntParam';
                   'sym_setDblParam';
                   'sym_getDblParam';
                   'sym_setStrParam';
                   'sym_getStrParam';
                   'sym_getInfinity';
                   'sym_loadProblemBasic';
                   'sym_loadProblem';
                   'sym_loadMPS';
                   'sym_getNumConstr';
                   'sym_getNumVar';
                   'sym_getNumElements';
                   'sym_isContinuous';
                   'sym_isBinary';
                   'sym_isInteger';
                   'sym_setContinuous';
                   'sym_setInteger';
                   'sym_getVarLower';
                   'sym_getVarUpper';
                   'sym_setVarLower';
                   'sym_setVarUpper';
                   'sym_getObjCoeff';
                   'sym_setObjCoeff';
                   'sym_getObjSense';
                   'sym_setObjSense';
                   'sym_getRhs';
                   'sym_getConstrRange';
                   'sym_getConstrLower';
                   'sym_getConstrUpper';
                   'sym_setConstrLower';
                   'sym_setConstrUpper';
                   'sym_setConstrType';
                   'sym_getMatrix';
                   'sym_getConstrSense';
                   'sym_addConstr';
                   'sym_addVar';
                   'sym_deleteVars';
                   'sym_deleteConstrs';
                   'sym_getPrimalBound';
                   'sym_setPrimalBound';
                   'sym_setVarSoln';
                   'sym_solve';
                   'sym_getStatus';
                   'sym_isOptimal';
                   'sym_isInfeasible';
                   'sym_isAbandoned';
                   'sym_isIterLimitReached';
                   'sym_isTimeLimitReached';
                   'sym_isTargetGapAchieved';
                   'sym_getVarSoln';
                   'sym_getObjVal';
                   'sym_getIterCount';
                   'sym_getConstrActivity';
];
addinter(libsymphonytools_path + filesep() + 'libsymphonytools' + getdynlibext(), 'libsymphonytools', list_functions);
// remove temp. variables on stack
clear libsymphonytools_path;
clear bOK;
clear ilib;
clear list_functions;
// ----------------------------------------------------------------------------
