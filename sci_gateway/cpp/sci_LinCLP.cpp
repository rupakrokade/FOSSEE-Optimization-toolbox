/*
 * Linear Solver Toolbox for Scilab using CLP library
 * Authors :
	Guru Pradeep Reddy
    Bhanu Priya Sayal
*/

#include "LinCLP.hpp"
extern "C"{
#include "api_scilab.h"
#include "Scierror.h"
#include "localization.h"
#include "sciprint.h"
#include "sci_iofunc.hpp"

//creating a problem pointer using Base class of OsiSolverInterface and
//Instantiate the object using specific derived class of ClpSolver
OsiSolverInterface* si = new OsiClpSolverInterface();

LinCLP::~LinCLP()
		 {
		 	free(objMatrix_);
		 	free(conMatrix_);
		 	free(conlb_);
		 	free(conub_);
		 	free(lb_);
		 	free(ub_);
			free(xValue_);
			free(Zl_);
			free(dual_);}

//Clp Solver function definition
LinCLP::LinCLP(int numVars_ , int numCons_ ,double objMatrix_[] , double conMatrix_[] , double conlb_[] , double conub_[] ,double lb_[] , double ub_[], double options_[])
{
   
   //Defining the constraint matrix
   CoinPackedMatrix *matrix =  new CoinPackedMatrix(false , 0 , 0);
   matrix->setDimensions(0 , numVars_);
   for(int i=0 ; i<numCons_ ; i++)
   	{
    	CoinPackedVector row;
	 	for(int j=0 ; j<numVars_ ; j++)
	 		{
   				row.insert(j, conMatrix_[i+j*numCons_]);
   	 		}
			
        matrix->appendRow(row);
   	}

   //setting options for maximum iterations
   si->setIntParam(OsiMaxNumIteration,options_[0]);

   //Load the problem to OSI
   si->loadProblem(*matrix , lb_ , ub_, objMatrix_ , conlb_ , conub_);

   //Solve the problem
   si->initialSolve();
  
}

   //Output the solution to Scilab
   //get solution for x
   const double* LinCLP::getX()
	{ 
        xValue_ = si->getColSolution();
		return xValue_;
	}

   //get objective value
   double LinCLP::getObjVal()
	{
		objValue_ = si->getObjValue();
		return objValue_;
	}
   
   //get exit status 
   int LinCLP::returnStatus()
	{
   		status_;
   		if(si->isProvenOptimal())
    			status_=0;
   		else if(si->isProvenPrimalInfeasible())
        		status_=1;
   		else if(si->isProvenDualInfeasible())
        		status_=2;
   		else if(si->isIterationLimitReached())
        		status_=3;
   		else if(si->isAbandoned())
        		status_=4;
   		else if(si->isPrimalObjectiveLimitReached())
        		status_=5;
   		else if(si->isDualObjectiveLimitReached())
        		status_=6;
		return status_;
	}

   //get number of iterations
   double LinCLP::iterCount()
	{
		iterations_ = si->getIterationCount(); 
		return iterations_;
	}

   //get lower vector
   const double* LinCLP::getReducedCost()
	{
		Zl_ = si->getReducedCost();
		return Zl_;
	}

   //get dual vector
   double* LinCLP::getDual()
	{
		dual_ = si->getRowPrice();
        return dual_;
	}

}
   
