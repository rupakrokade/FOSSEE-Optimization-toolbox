// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: R.Vidyadhar & Vignesh Kannan
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include "sci_iofunc.hpp"
#include "IpIpoptApplication.hpp"
#include "minconNLP.hpp"
#include <IpSolveStatistics.hpp>

extern "C"
{
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>

int sci_solveminconp(char *fname)
{
	using namespace Ipopt;

	CheckInputArgument(pvApiCtx, 16, 16);
	CheckOutputArgument(pvApiCtx, 12, 12);
	
	// Error management variable
	SciErr sciErr;

	//Function pointers, input matrix(Starting point) pointer, flag variable 
	double *x0ptr=NULL, *lbptr=NULL, *ubptr=NULL,*Aptr=NULL, *bptr=NULL, *Aeqptr=NULL, *beqptr=NULL;
	double nonlinCon=0,nonlinIneqCon=0;
        

    	// Input arguments
	double *cpu_time=NULL,*max_iter=NULL;
	static unsigned int nVars = 0,nCons = 0;
	unsigned int temp1 = 0,temp2 = 0, iret = 0;
	int x0_rows=0, x0_cols=0, lb_rows=0, lb_cols=0, ub_rows=0, ub_cols=0, A_rows=0, A_cols=0, b_rows=0, b_cols=0, Aeq_rows=0, Aeq_cols=0, beq_rows=0, beq_cols=0;
	
	// Output arguments	
	double  ObjVal=0,iteration=0,cpuTime=0,fobj_eval=0;
	double dual_inf, constr_viol, complementarity, kkt_error;
	const double *fX = NULL, *fGrad =  NULL;
	const double *fHess =  NULL;
	const double *fLambda = NULL;
	const double *fZl=NULL;
	const double *fZu=NULL;
	int rstatus = 0;
	int int_fobj_eval, int_constr_eval, int_fobj_grad_eval, int_constr_jac_eval, int_hess_eval;

	////////// Manage the input argument //////////

	
	//Getting matrix representing linear inequality constraints 
	if(getDoubleMatrixFromScilab(2, &A_rows, &A_cols, &Aptr))
	{
		return 1;
	}

	//Getting matrix representing bounds of linear inequality constraints 
	if(getDoubleMatrixFromScilab(3, &b_rows, &b_cols, &bptr))
	{
		return 1;
	}
	
	//Getting matrix representing linear equality constraints 
	if(getDoubleMatrixFromScilab(4, &Aeq_rows, &Aeq_cols, &Aeqptr))
	{
		return 1;
	}

	//Getting matrix representing bounds of linear inequality constraints 
	if(getDoubleMatrixFromScilab(5, &beq_rows, &beq_cols, &beqptr))
	{
		return 1;
	}

	//Getting matrix representing linear inequality constraints 
	if(getDoubleMatrixFromScilab(6, &lb_rows, &lb_cols, &lbptr))
	{
		return 1;
	}

	//Getting matrix representing linear inequality constraints 
	if(getDoubleMatrixFromScilab(7, &ub_rows, &ub_cols, &ubptr))
	{
		return 1;
	}

	//Number of non-linear constraints
	if(getDoubleFromScilab(8, &nonlinCon))
	{
		return 1;
	}

	//Number of non-linear inequality constraints
	if(getDoubleFromScilab(9, &nonlinIneqCon))
	{
		return 1;
	}
	
	//x0(starting point) matrix from scilab
	if(getDoubleMatrixFromScilab(14, &x0_rows, &x0_cols, &x0ptr))
	{
		return 1;
	}

    //Getting number of iterations
    if(getFixedSizeDoubleMatrixInList(15,2,temp1,temp2,&max_iter))
	{
		return 1;
	}

	//Getting Cpu Time
	if(getFixedSizeDoubleMatrixInList(15,4,temp1,temp2,&cpu_time))
	{
		return 1;
	}

	//Number of variables and constraints
	nVars = x0_cols;
	nCons = A_rows + Aeq_rows + nonlinCon;


        
        // Starting Ipopt

	SmartPtr<minconNLP> Prob = new minconNLP(nVars, nCons, x0ptr, Aptr, bptr, Aeqptr, beqptr, A_rows, A_cols, b_rows, b_cols, Aeq_rows, Aeq_cols, beq_rows, beq_cols, lbptr, ubptr, nonlinCon, nonlinIneqCon);
	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

	////////// Managing the parameters //////////

	app->Options()->SetNumericValue("tol", 1e-6);
	app->Options()->SetIntegerValue("max_iter", (int)*max_iter);
	app->Options()->SetNumericValue("max_cpu_time", *cpu_time);
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");

	///////// Initialize the IpoptApplication and process the options /////////
	ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Solve_Succeeded) 
	{
	  	sciprint("\n*** Error during initialization!\n");
   		return (int) status;
 	 }

	 // Ask Ipopt to solve the problem
	 status = app->OptimizeTNLP((SmartPtr<TNLP>&)Prob);
	 //Get the solve statistics
	 cpuTime = app->Statistics()->TotalCPUTime();
	 app->Statistics()->NumberOfEvaluations(int_fobj_eval, int_constr_eval, int_fobj_grad_eval, int_constr_jac_eval, int_hess_eval);
	 app->Statistics()->Infeasibilities(dual_inf, constr_viol, complementarity, kkt_error);
	 rstatus = Prob->returnStatus();
	 fobj_eval=(double)int_fobj_eval;

	////////// Manage the output argument //////////

	fX = Prob->getX();
	fGrad = Prob->getGrad();
	fHess = Prob->getHess();
	fLambda = Prob->getLambda();
	fZl = Prob->getZl();
	fZu = Prob->getZu();
	ObjVal = Prob->getObjVal();
	iteration = (double)app->Statistics()->IterationCount();

	if (returnDoubleMatrixToScilab(1, 1, nVars, fX))
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

	if (returnDoubleMatrixToScilab(4, 1, 1, &iteration))
	{
		return 1;
	}
	
	if (returnDoubleMatrixToScilab(5, 1, 1, &cpuTime))
	{
		return 1;
	}
	
	if (returnDoubleMatrixToScilab(6, 1, 1, &fobj_eval))
	{
		return 1;
	}
	
	if (returnDoubleMatrixToScilab(7, 1, 1, &dual_inf))
	{
		return 1;
	}

	if (returnDoubleMatrixToScilab(8, 1, nCons, fLambda))
	{
		return 1;
	}

	if (returnDoubleMatrixToScilab(9, 1, nVars, fZl))
	{
		return 1;
	}

	if (returnDoubleMatrixToScilab(10, 1, nVars, fZu))
	{
		return 1;
	}
		
	if (returnDoubleMatrixToScilab(11, 1, nVars, fGrad))
	{
		return 1;
	}

	if (returnDoubleMatrixToScilab(12, 1, nVars*nVars, fHess))
	{
		return 1;
	}
	
	// As the SmartPtrs go out of scope, the reference count
	// will be decremented and the objects will automatically
	// be deleted.*/
	
	return 0;
}
}
