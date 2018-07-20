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
#include "minconTMINLP.hpp"
#include "BonCbc.hpp"
#include "BonBonminSetup.hpp"

#include "BonOACutGenerator2.hpp"
#include "BonEcpCuts.hpp"
#include "BonOaNlpOptim.hpp"

#include "sci_iofunc.hpp"
extern  "C"
{
#include "call_scilab.h"
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>

int cpp_intfmincon(char *fname)
{
  	using namespace Ipopt;
  	using namespace Bonmin;

	CheckInputArgument(pvApiCtx, 13, 13);
	CheckOutputArgument(pvApiCtx, 3, 3);

	// Input arguments
	Number *integertolerance=NULL, *maxnodes=NULL, *allowablegap=NULL, *cputime=NULL,*max_iter=NULL;
	Number *x0 = NULL, *lb = NULL, *ub = NULL,*conLb = NULL, *conUb = NULL,*LC = NULL;
	static unsigned int nVars = 0,nCons = 0;
	unsigned int temp1 = 0,temp2 = 0, iret = 0;
	int x0_rows, x0_cols,intconSize;
	Number *intcon = NULL,*options=NULL, *ifval=NULL;
	
	// Output arguments
<<<<<<< HEAD
	Number *fX = NULL, ObjVal=0,iteration=0,cpuTime=0,fobj_eval=0;
	Number dual_inf, constr_viol, complementarity, kkt_error;
	int rstatus = 0;
	
	if(getDoubleMatrixFromScilab(6, &nVars, &x0_cols, &x0))
=======
	Number  ObjVal=0,iteration=0,cpuTime=0,fobj_eval=0;	//Number *fX = NULL, ObjVal=0,iteration=0,cpuTime=0,fobj_eval=0;
	Number dual_inf, constr_viol, complementarity, kkt_error;
	int rstatus = 0;

	const double *fX = NULL;	//changed fX from Ipopt::Number* to const double*
	
	if(getDoubleMatrixFromScilab(6, (int*)&nVars, &x0_cols, &x0))	//typecast nVars from unsigned int* to *int
>>>>>>> 611a2eae153e83b49d73e0277def7c3f865b4eb3
	{
		return 1;
	}
	
	if(getDoubleMatrixFromScilab(7, &x0_rows, &x0_cols, &lb))
	{
		return 1;
	}

	if(getDoubleMatrixFromScilab(8, &x0_rows, &x0_cols, &ub))
	{
		return 1;
	}

<<<<<<< HEAD
	if(getDoubleMatrixFromScilab(9, &nCons, &x0_cols, &conLb))
=======
	if(getDoubleMatrixFromScilab(9, (int*)&nCons, &x0_cols, &conLb))	//typecast nCons from unsigned int* to int* 
>>>>>>> 611a2eae153e83b49d73e0277def7c3f865b4eb3
	{
		return 1;
	}

	if(getDoubleMatrixFromScilab(10, &x0_rows, &x0_cols, &conUb))
	{
		return 1;
	}

	// Getting intcon
<<<<<<< HEAD
	if (getDoubleMatrixFromScilab(11,&intconSize,&temp2,&intcon))
=======
	if (getDoubleMatrixFromScilab(11,&intconSize,(int*)&temp2,&intcon))	//typecast temp2 from unsigned int* to int* 
>>>>>>> 611a2eae153e83b49d73e0277def7c3f865b4eb3
	{
		return 1;
	}

<<<<<<< HEAD
	if (getDoubleMatrixFromScilab(13,&temp1,&temp2,&LC))
=======
	if (getDoubleMatrixFromScilab(13,(int*)&temp1,(int*)&temp2,&LC))	//typecast temp1 and temp2 from unsigned int* to int* 
>>>>>>> 611a2eae153e83b49d73e0277def7c3f865b4eb3
	{
		return 1;
	}

    //Initialization of parameters
	temp1 = 1;
	temp2 = 1;

	//Getting parameters
	if (getFixedSizeDoubleMatrixInList(12,2,temp1,temp2,&integertolerance))
	{
		return 1;
	}
	if (getFixedSizeDoubleMatrixInList(12,4,temp1,temp2,&maxnodes))
	{
		return 1;
	}
	if (getFixedSizeDoubleMatrixInList(12,6,temp1,temp2,&cputime))
	{
		return 1;
	}
	if (getFixedSizeDoubleMatrixInList(12,8,temp1,temp2,&allowablegap))
	{
		return 1;
	}
	if (getFixedSizeDoubleMatrixInList(12,10,temp1,temp2,&max_iter))
	{
		return 1;
	}
	
	SmartPtr<minconTMINLP> tminlp = new minconTMINLP(nVars,x0,lb,ub,(unsigned int)LC,nCons,conLb,conUb,intconSize,intcon);
	
	BonminSetup bonmin;
	bonmin.initializeOptionsAndJournalist();
	bonmin.options()->SetStringValue("mu_oracle","loqo");
	bonmin.options()->SetIntegerValue("bonmin.print_level",5);
    bonmin.options()->SetNumericValue("bonmin.integer_tolerance", *integertolerance);
    bonmin.options()->SetIntegerValue("bonmin.node_limit", (int)*maxnodes);
    bonmin.options()->SetNumericValue("bonmin.time_limit", *cputime);
    bonmin.options()->SetNumericValue("bonmin.allowable_gap", *allowablegap);
    bonmin.options()->SetIntegerValue("bonmin.iteration_limit", (int)*max_iter);

	//Now initialize from tminlp
	bonmin.initialize(GetRawPtr(tminlp));

	//Set up done, now let's branch and bound
	try {
	Bab bb;
	bb(bonmin);//process parameter file using Ipopt and do branch and bound using Cbc
	}
	catch(TNLPSolver::UnsolvedError *E) {
	}
	catch(OsiTMINLPInterface::SimpleError &E) {
	}
	catch(CoinError &E) {
	}
	rstatus=tminlp->returnStatus();

	if(rstatus==0 ||rstatus== 3)
	{
		fX = tminlp->getX();
		ObjVal = tminlp->getObjVal();
		if (returnDoubleMatrixToScilab(1, nVars, 1, fX))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(2, 1, 1, &ObjVal))
		{
			return 1;
		}

		if (returnIntegerMatrixToScilab(3, 1, 1, &rstatus))
		{
			return 1;
		}
	
	}
	else
	{
		if (returnDoubleMatrixToScilab(1, 0, 0, fX))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(2, 1, 1, &ObjVal))
		{
			return 1;
		}

		if (returnIntegerMatrixToScilab(3, 1, 1, &rstatus))
		{
			return 1;
		}
	
	}

	return 0;
	}
}

