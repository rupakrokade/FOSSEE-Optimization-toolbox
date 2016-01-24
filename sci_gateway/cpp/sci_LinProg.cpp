/*
 * Linear Solver Toolbox for Scilab using CLP library
 * Authors :
	Guru Pradeep Reddy
    Bhanu Priya Sayal
*/

#include "sci_iofunc.hpp"
#include "LinCLP.hpp"

extern "C"{
#include <api_scilab.h>
#include <Scierror.h>
#include <localization.h>
#include <sciprint.h>

//Solver function
int sci_linearprog(char *fname) 
{
	//Objective function
	double* obj;  
	//Constraint matrix coefficients
	double* conMatrix;  
	//Constraints upper bound
	double* conlb;
	//Constraints lower bound
	double* conub;
	//Lower bounds for variables
	double* lb;  
	//Upper bounds for variables
	double* ub;
	//options for maximum iterations and writing mps
	double* options;
	//Flag for Mps
	double flagMps;
	//mps file path
	char * mpsFile;
	//Error structure in Scilab  
	SciErr sciErr;
	//Number of rows and columns in objective function
	int nVars=0, nCons=0,temp1=0,temp2=0;
	
	CheckInputArgument(pvApiCtx , 9 , 9);             //Checking the input arguments
	CheckOutputArgument(pvApiCtx , 6, 6);               //Checking the output arguments

	////////// Manage the input argument //////////
	
	//Number of Variables
	if(getIntFromScilab(1,&nVars))
	{
		return 1;
	}

	//Number of Constraints
	if (getIntFromScilab(2,&nCons))
	{
		return 1;
	}

	//Objective function from Scilab
	temp1 = nVars;
	temp2 = nCons;
	if (getFixedSizeDoubleMatrixFromScilab(3,1,temp1,&obj))
	{
		return 1;
	}

	if (nCons!=0)
	{
		//conMatrix matrix from scilab
		temp1 = nCons;
		temp2 = nVars;

		if (getFixedSizeDoubleMatrixFromScilab(4,temp1,temp2,&conMatrix))
		{
			return 1;
		}

		//conLB matrix from scilab
		temp1 = nCons;
		temp2 = 1;
		if (getFixedSizeDoubleMatrixFromScilab(5,temp1,temp2,&conlb))
		{
			return 1;
		}

		//conUB matrix from scilab
		if (getFixedSizeDoubleMatrixFromScilab(6,temp1,temp2,&conub))
		{
			return 1;
		}

	}

	//lb matrix from scilab
	temp1 = 1;
	temp2 = nVars;
	if (getFixedSizeDoubleMatrixFromScilab(7,temp1,temp2,&lb))
	{
		return 1;
	}


	//ub matrix from scilab
	if (getFixedSizeDoubleMatrixFromScilab(8,temp1,temp2,&ub))
	{
		return 1;
	}

	//get options from scilab
	if(getFixedSizeDoubleMatrixInList(9 , 2 , 1 , 1 , &options))
	{
		return 1;      
	}

	//Call to the Clp Solver
	LinCLP* Prob = new LinCLP(nVars,nCons,obj,conMatrix,conlb,conub,lb,ub,options);

	//Output the solution to Scilab
	//get solution for x
	double* xValue = Prob->getX();

	//get objective value
	double objValue = Prob->getObjVal();

	//get Status value
	double status = Prob->returnStatus();

	//get number of iterations
	double iterations = Prob->iterCount();

	//get reduced cost
	double* reducedCost = Prob->getReducedCost();

	//get dual vector
	double* dual = Prob->getDual();

	returnDoubleMatrixToScilab(1 , 1 , nVars , xValue);
	returnDoubleMatrixToScilab(2 , 1 , 1 , &objValue);
	returnDoubleMatrixToScilab(3 , 1 , 1 , &status);
	returnDoubleMatrixToScilab(4 , 1 , 1 , &iterations);
	returnDoubleMatrixToScilab(5 , 1 , nVars , reducedCost);
	returnDoubleMatrixToScilab(6 , 1 , nCons , dual);

	}
}


  	
   
