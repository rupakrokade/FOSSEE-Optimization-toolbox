// Copyright (C) 2016 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Harpreet Singh
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include "CoinPragma.hpp"
#include "CoinTime.hpp"
#include "CoinError.hpp"

#include "BonOsiTMINLPInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "minbndTMINLP.hpp"
#include "BonCbc.hpp"
#include "BonBonminSetup.hpp"

#include "BonOACutGenerator2.hpp"
#include "BonEcpCuts.hpp"
#include "BonOaNlpOptim.hpp"


extern  "C"
{
#include <wchar.h>
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>

const char fname[] = "cpp_intfminbnd";
int cpp_intfminbnd(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out) 
{
  	using namespace Ipopt;
  	using namespace Bonmin;

	
	// Input arguments

	static unsigned int nVars = 0;
	int x0_rows, x0_cols,intconSize, intconSize2=0; //changed temp2 from unsigned int to int
	Number *intcon = NULL,*options=NULL, *ifval=NULL;
	double* lb=NULL;
	double* ub=NULL;

	
	// Output arguments
	Number ObjVal=0,iteration=0,cpuTime=0,fobj_eval=0;	//Number *fX = NULL, ObjVal=0,iteration=0,cpuTime=0,fobj_eval=0;
	Number dual_inf, constr_viol, complementarity, kkt_error;
	int rstatus = 0;
	
	const double *fX = NULL;	//changed fX from Ipopt::Number* to const double*

	if (nin !=8)  //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname, 8);
        	return STATUS_ERROR; 
	}
	
	if (nout !=3) //Checking the output arguments

	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname, 3);
		return 1;
	} 


	//Lb
	if (scilab_isDouble(env, in[3]) == 0 || scilab_isMatrix2d(env, in[3]) == 0)
	{
		Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 4);
		return 1;
	}
	
	scilab_getDoubleArray(env, in[3], &lb);
	int size1 = scilab_getDim2d(env, in[3], &x0_rows, &x0_cols);
	
	//Ub
	if (scilab_isDouble(env, in[4]) == 0 || scilab_isMatrix2d(env, in[4]) == 0)
	{
		Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 4);
		return 1;
	}
	
	scilab_getDoubleArray(env, in[4], &ub);

	// Getting intcon
	

	if (scilab_isDouble(env, in[5]) == 0 || scilab_isMatrix2d(env, in[5]) == 0)
	{
		Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 6);
		return 1;
	}
	
	scilab_getDoubleArray(env, in[5], &intcon);
	size1 = scilab_getDim2d(env, in[5], &intconSize, &intconSize2);

    //Initialization of parameters
	nVars=x0_rows;

	//Getting parameters
	if (scilab_isList(env, in[6]) == 0)
    {
        Scierror(999, "%s: Wrong type for input argument #%d: A list expected.\n", fname, 6);
        return 1;
    }

	scilabVar temp1 = scilab_getListItem( env, in[6], 1);
	scilabVar temp2 = scilab_getListItem( env, in[6], 3);
	scilabVar temp3 = scilab_getListItem( env, in[6], 5);
	scilabVar temp4 = scilab_getListItem( env, in[6], 7);
	scilabVar temp5 = scilab_getListItem( env, in[6], 9);

	double integertolerance=0, allowable_gap=0, maxnodes =0,  cputime=0, maxiter=0;

	scilab_getDouble(env, temp1, &integertolerance);
	scilab_getDouble(env, temp2, &maxnodes);
	scilab_getDouble(env, temp3, &cpuTime);
	scilab_getDouble(env, temp4, &allowable_gap);
	scilab_getDouble(env, temp5, &maxiter);


	int max_nodes = (int)maxnodes;
	int cpu_time = (int)cpuTime;
	int iterLim = (int)maxiter;


	SmartPtr<minbndTMINLP> tminlp = new minbndTMINLP(env, in, nVars,lb,ub,intconSize,intcon);

	BonminSetup bonmin;
	bonmin.initializeOptionsAndJournalist();

	bonmin.options()->SetStringValue("mu_oracle","loqo");
    bonmin.options()->SetNumericValue("bonmin.integer_tolerance", integertolerance);
    bonmin.options()->SetIntegerValue("bonmin.node_limit", max_nodes);
    bonmin.options()->SetNumericValue("bonmin.time_limit", cpu_time);
    bonmin.options()->SetNumericValue("bonmin.allowable_gap", allowable_gap);
    bonmin.options()->SetIntegerValue("bonmin.iteration_limit", iterLim);

	//Now initialize from tminlp
	bonmin.initialize(GetRawPtr(tminlp));

	//Set up done, now let's branch and bound
	try {
	Bab bb;
	bb(bonmin);//process parameter file using Ipopt and do branch and bound using Cbc
	}
	catch(TNLPSolver::UnsolvedError *E) {
	Scierror(999, "\nIpopt has failed to solve the problem!\n");
	}
	catch(OsiTMINLPInterface::SimpleError &E) {
	  Scierror(999, "\nFailed to solve a problem!\n");
	}
	catch(CoinError &E) {
	  Scierror(999, "\nFailed to solve a problem!\n");
	}
	rstatus=tminlp->returnStatus();



	if(rstatus==0 ||rstatus== 3)
	{

		fX = tminlp->getX();
		ObjVal = tminlp->getObjVal();

		out[0] = scilab_createDoubleMatrix2d(env, 1, nVars, 0);
		scilab_setDoubleArray(env, out[0], fX);

		out[1] = scilab_createDouble(env, ObjVal);

		out[2] = scilab_createDouble(env, rstatus);
		
	}
	else
	{

		out[0] = scilab_createDoubleMatrix2d(env, 0, 0, 0);
		scilab_setDoubleArray(env, out[0], fX);

		out[1] = scilab_createDouble(env, ObjVal);

		out[2] = scilab_createDouble(env, rstatus);
  }

	return 0;
	}
}

