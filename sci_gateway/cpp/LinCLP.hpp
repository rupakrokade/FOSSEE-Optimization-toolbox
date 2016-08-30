// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Guru Pradeep Reddy Bhanu Priya Sayal
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#ifndef __LinCLP_HPP__
#define __LinCLP_HPP__

#include"OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"

class LinCLP
{
	private:

		int numVars_;				//Number of variables

		int numCons_;          		//Number of inequality constraints

		double* objMatrix_[];   		//Objective function vector 
 
		double* conMatrix_[];   		//Inequality constraint matrix
  
		double* conlb_[];     		//Inequality constraint vector

		double* conub_[];   		//Equality constraint vector
 
		double* lb_[];       		//Lower bounds for all variables

		double* ub_[];       		//Upper bounds for all variables

        double options_[];          //options for setting maximum iterations and writing mps and lp files
      
		double* xValue_ = NULL;		//Optimal value of variables 

		double objValue_ =0 ;		//Optimal values of objective

		double status_ = 0; 		//Return Status

		double iterations_ = 0;		//Number of iteration 

		double* Zl_ = NULL;			//Reduced cost

		double* Zu_ = NULL;			//Reduced cost

		double* dual_ = NULL;		// Dual of the solution


	public:
/*
 * Constructor 
*/
		LinCLP(int numVars_ , int numCons_ ,double objMatrix_[] , double conMatrix_[] , double conlb_[] , double conub_[] ,double lb_[] , double ub_[], double options_[]);
		

		virtual ~LinCLP();			//Destructor to free memory

		const double* getX();   	//Returns a pointer to matrix of size 
									//1*numVars with final values for the objective variables

	    double getObjVal();     	//Returns the output of the final value of the objective

		int returnStatus();     	//Returns the status of the problem

	    double iterCount();     	//Returns the iteration count

        const double* getReducedCost();   //Returns a pointer to matrix of size 
									      //1*numVars with values for lower dual vector

	    double* getDual();   	//Returns a pointer to matrix of size
								//1*numCons with values for dual vector

};

#endif __LinCLP_HPP__

