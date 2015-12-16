// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// Author: Harpreet Singh
// Organization: FOSSEE, IIT Bombay
// Email: harpreet.mertia@gmail.com
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt


function [xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg (varargin)
	// Solves nonnegative least-squares curve fitting problems.
	//
	//   Calling Sequence
	//   x = lsqnonneg(C,d)
	//   x = lsqnonneg(C,d,param)
	//   [xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg( ... )
	//   
	//   Parameters
	//   C : a matrix of doubles, represents the multiplier of the solution x in the expression C*x - d. C is M-by-N, where M is the number of equations, and N is the number of elements of x.
	//   d : a vector of doubles, represents the additive constant term in the expression C*x - d. d is M-by-1, where M is the number of equations.
	//   xopt : a vector of doubles, the computed solution of the optimization problem.
	//   resnorm : a double, objective value returned as the scalar value norm(C*x-d)^2.
	//   residual : a vector of doubles, solution residuals returned as the vector C*x-d.
	//   exitflag : Integer identifying the reason the algorithm terminated.
	//   output : Structure containing information about the optimization.
	//   lambda : Structure containing the Lagrange multipliers at the solution x (separated by constraint type).
	//   
	//   Description
	//   Solves nonnegative least-squares curve fitting problems specified by :
	//
	//   <latex>
	//    \begin{eqnarray}
	//    &\mbox{min}_{x}
	//    & 1/2||C*x - d||_2^2  \\
	//    & & x \geq 0 \\
	//    \end{eqnarray}
	//   </latex>
	//   
	//   We are calling IPOpt for solving the nonnegative least-squares curve fitting problems, IPOpt is a library written in C++. The code has been written by ​Andreas Wächter and ​Carl Laird.
	//    
	// Examples 
	// A basic lsqnonneg problem
	//	C = [
	//		0.0372    0.2869
	//		0.6861    0.7071
	//		0.6233    0.6245
	//		0.6344    0.6170];
	//	d = [
	//    	0.8587
	//    	0.1781
	//   	0.0747
	//	    0.8405];
	// [xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg(C,d)
	//
	// Authors
	// Harpreet Singh


	//To check the number of input and output argument
	[lhs , rhs] = argn();

	//To check the number of argument given by user
	if ( rhs < 2 | rhs > 3 ) then
		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set of [2 3]"), "lsqlin", rhs);
		error(errmsg)
	end

	C = varargin(1);
	d = varargin(2);
	nbVar = size(C,2);
	if ( rhs<3 | size(varargin(3)) ==0 ) then
		param = list();
	else
		param =varargin(10);
	end

	if (type(param) ~= 15) then
		errmsg = msprintf(gettext("%s: param should be a list "), "lsqlin");
		error(errmsg);
	end


	if (modulo(size(param),2)) then
		errmsg = msprintf(gettext("%s: Size of parameters should be even"), "lsqlin");
		error(errmsg);
	end

	options = list(	"MaxIter"     , [3000], ...
					"CpuTime"   , [600] ...
	);

	for i = 1:(size(param))/2

		select param(2*i-1)
			case "MaxIter" then
				options(2*i) = param(2*i);
			case "CpuTime" then
				options(2*i) = param(2*i);
			else
				errmsg = msprintf(gettext("%s: Unrecognized parameter name ''%s''."), "lsqlin", param(2*i-1));
				error(errmsg)
		end
	end

	// Check if the user gives row vector 
	// and Changing it to a column matrix


	if (size(d,2)== [nbVar]) then
		d=d';
	end

	//Check the size of f which should equal to the number of variable
	if ( size(d,1) ~= size(C,1)) then
		errmsg = msprintf(gettext("%s: The number of rows in C must be equal the number of elements of d"), "lsqlin");
		error(errmsg);
	end

	//Converting it into Quadratic Programming Problem

	Q = C'*C;
	p = [-C'*d]';
	op_add = d'*d;
	LB = repmat(0,1,nbVar);
	UB = repmat(%inf,1,nbVar);	
	x0 = repmat(0,1,nbVar);;
	conMatrix = [];
	nbCon = size(conMatrix,1);
	conLB = [];
	conUB = [] ; 
	[xopt,fopt,status,iter,Zl,Zu,lmbda] = solveqp(nbVar,nbCon,Q,p,conMatrix,conLB,conUB,LB,UB,x0,options);

	xopt = xopt';
	residual = -1*(C*xopt-d);
	resnorm = residual'*residual;
	exitflag = status;
	output = struct("Iterations"      , []);
	output.Iterations = iter;
   lambda = struct("lower"           , [], ..
                   "upper"           , []);
   
   lambda.lower = Zl;
   lambda.upper = Zu;

	select status 
		case 0 then
			printf("\nOptimal Solution Found.\n");
		case 1 then
			printf("\nMaximum Number of Iterations Exceeded. Output may not be optimal.\n");
		case 2 then
			printf("\nMaximum CPU Time exceeded. Output may not be optimal.\n");
		case 3 then
			printf("\nStop at Tiny Step\n");
		case 4 then
			printf("\nSolved To Acceptable Level\n");
		case 5 then
			printf("\nConverged to a point of local infeasibility.\n");
		case 6 then
			printf("\nStopping optimization at current point as requested by user.\n");
		case 7 then
			printf("\nFeasible point for square problem found.\n");
		case 8 then 
			printf("\nIterates diverging; problem might be unbounded.\n");
		case 9 then
			printf("\nRestoration Failed!\n");
		case 10 then
			printf("\nError in step computation (regularization becomes too large?)!\n");
		case 12 then
			printf("\nProblem has too few degrees of freedom.\n");
		case 13 then
			printf("\nInvalid option thrown back by IPOpt\n");
		case 14 then
			printf("\nNot enough memory.\n");
		case 15 then
			printf("\nINTERNAL ERROR: Unknown SolverReturn value - Notify IPOPT Authors.\n");
		else
			printf("\nInvalid status returned. Notify the Toolbox authors\n");
		break;
	end

endfunction
