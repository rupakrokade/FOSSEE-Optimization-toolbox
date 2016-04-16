// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Salman Anis, Harpreet Singh
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

function [xopt,fopt,exitflag,output,lambda,gradient,hessian] = fseminf (varargin)
  	// Solves a multi-variable constrainted optimization problem
  	//
  	//   Calling Sequence
  	//   xopt = fseminf(fun,x0,ntheta,seminfcon)
  	//   xopt = fseminf(fun,x0,ntheta,seminfcon,A,b)
  	//   xopt = fseminf(fun,x0,ntheta,seminfcon,A,b,Aeq,beq)
  	//   xopt = fseminf(fun,x0,ntheta,seminfcon,A,b,Aeq,beq,lb,ub)
  	//   xopt = fseminf(fun,x0,ntheta,seminfcon,A,b,Aeq,beq,lb,ub,options)
  	//   [xopt,fopt] = fseminf(.....)
  	//   [xopt,fopt,exitflag]= fseminf(.....)
  	//   [xopt,fopt,exitflag,output]= fseminf(.....)
  	//   [xopt,fopt,exitflag,output,lambda]=fseminf(.....)
  	//   [xopt,fopt,exitflag,output,lambda]=fseminf(.....)
  	//   [xopt,fopt,exitflag,output,lambda]=fseminf(.....)
  	//
  	//   Parameters
  	//   fun : a function, representing the objective function of the problem 
  	//   x0 : a vector of doubles, containing the starting values of variables of size (1 X n) or (n X 1) where 'n' is the number of Variables
	//	 ntheta : The number of semi-infinite constraints.
	//	 seminfcon : a function that calculates the vector of nonlinear inequality constraints c, a vector of nonlinear equality constraints ceq, and ntheta semi-infinite constraints. See below for details.
	//   A : a matrix of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b. 
	//   b : a vector of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b.
	//   Aeq : a matrix of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   beq : a vector of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   lb : a vector of double, contains lower bounds of the variables.
	//   ub : a vector of double, contains upper bounds of the variables.
  	//   options : a list, containing the option for user to specify. See below for details. 

	//   xopt : a vector of double, the computed solution of the optimization problem.
	//   fopt : a double, the value of the function at x.
	//   exitflag : The exit status. See below for details.
	//   output : The structure consist of statistics about the optimization. See below for details.
	//   lambda : The structure consist of the Lagrange multipliers at the solution of problem. See below for details.
  	//
  	//   Description
  	//   Search the minimum of a constrained optimization problem specified by :
  	//   Find the minimum of f(x) such that 
  	//
  	//   <latex>
  	//    \begin{eqnarray}
  	//    &\mbox{min}_{x}
  	//    & f(x) \\
  	//    & \text{subject to} & A*x \leq b \\
  	//    & & Aeq*x \ = beq\\
  	//    & & lb \leq x \leq ub \\
  	//	  & & c(x) \leq  0\\
  	//    & & ceq(x) \ =  0\\
	//	  & & K_i(x,w_i) \leq 0, 1 \leq i \leq n.
  	//    \end{eqnarray}
  	//   </latex>
  	//
  	//   The routine calls Ipopt for solving the Constrained Optimization problem, Ipopt is a library written in C++.
  	//
  	// The options allows the user to set various parameters of the Optimization problem. 
  	// It should be defined as type "list" and contains the following fields.
	// <itemizedlist>
	//   <listitem>Syntax : options= list("MaxIter", [---], "CpuTime", [---], "GradObj", ---);</listitem>
	//   <listitem>MaxIter : a Scalar, containing the Maximum Number of Iteration that the solver should take.</listitem>
	//   <listitem>CpuTime : a Scalar, containing the Maximum amount of CPU Time that the solver should take.</listitem>
	//   <listitem>GradObj : a function, representing the gradient function of the Objective in Vector Form.</listitem>
	//   <listitem>Default Values : options = list("MaxIter", [3000], "CpuTime", [600]);</listitem>
	// </itemizedlist>
	//
	// The exitflag allows to know the status of the optimization which is given back by Ipopt.
	// <itemizedlist>
	//   <listitem>exitflag=0 : Optimal Solution Found </listitem>
	//   <listitem>exitflag=1 : Maximum Number of Iterations Exceeded. Output may not be optimal.</listitem>
	//   <listitem>exitflag=2 : Maximum amount of CPU Time exceeded. Output may not be optimal.</listitem>
	//   <listitem>exitflag=3 : Stop at Tiny Step.</listitem>
	//   <listitem>exitflag=4 : Solved To Acceptable Level.</listitem>
	//   <listitem>exitflag=5 : Converged to a point of local infeasibility.</listitem>
	// </itemizedlist>
	//
	// For more details on exitflag see the ipopt documentation, go to http://www.coin-or.org/Ipopt/documentation/
	//
	// The output data structure contains detailed informations about the optimization process. 
	// It has type "struct" and contains the following fields.
	// <itemizedlist>
	//   <listitem>output.Iterations: The number of iterations performed during the search</listitem>
	//   <listitem>output.Cpu_Time: The total cpu-time spend during the search</listitem>
	//   <listitem>output.Objective_Evaluation: The number of Objective Evaluations performed during the search</listitem>
	//   <listitem>output.Dual_Infeasibility: The Dual Infeasiblity of the final soution</listitem>
	// </itemizedlist>
	//
	// The lambda data structure contains the Lagrange multipliers at the end 
	// of optimization. In the current version the values are returned only when the the solution is optimal. 
	// It has type "struct" and contains the following fields.
	// <itemizedlist>
	//   <listitem>lambda.lower: The Lagrange multipliers for the lower bound constraints.</listitem>
	//   <listitem>lambda.upper: The Lagrange multipliers for the upper bound constraints.</listitem>
	//   <listitem>lambda.eqlin: The Lagrange multipliers for the linear equality constraints.</listitem>
	//   <listitem>lambda.ineqlin: The Lagrange multipliers for the linear inequality constraints.</listitem>
	//   <listitem>lambda.eqnonlin: The Lagrange multipliers for the non-linear equality constraints.</listitem>
	//   <listitem>lambda.ineqnonlin: The Lagrange multipliers for the non-linear inequality constraints.</listitem>
	// </itemizedlist>
	//
  	// Examples    
	//    function [y] = obj(x)
	//        y = (x-1)^2;
	//    endfunction
	//    function [c, ceq, K1, s] = seminfcon(x,s)
	//    // No finite nonlinear inequality and equality constraints
	//    c = [];
	//    ceq = [];
	//    // Sample set
	//    if isnan(s)
	//        // Initial sampling interval
	//        s = [0.01 0];
	//    end
	//    t = 0:s(1):1;
	//    // Evaluate the semi-infinite constraint
	//    K1 = (x - 0.5) - (t - 0.5).^2;
	//    endfunction
	//	  x = fseminf(obj,0.2,1,seminfcon)
	//    
  	// Authors
  	// 	Salman Anis, Harpreet Singh
 	
	
	//To check the number of input and output arguments
   	[lhs , rhs] = argn();
	
	//To check the number of arguments given by the user
   	if ( rhs<4 | rhs==5 | rhs==7 | rhs>13 ) then
    		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while it should be 4,6,8,10,11"), "fseminf", rhs);
    		error(errmsg)
   	end
 
	//Storing the Input Parameters  
   	_fun   	  = varargin(1);
   	x0   	  = varargin(2);
	ntheta	  = varargin(3);
	seminfcon = varargin(4); 
    nbVar = size(x0,'*');

	if(nbVar == 0) then
		errmsg = msprintf(gettext("%s: Cannot determine the number of variables because input initial guess is empty"), "lsqcurvefit");
		error(errmsg);
	end

   	A    	 = [];
   	b    	 = [];
   	Aeq  	 = [];
   	beq  	 = [];
   	lb       = [];
   	ub       = [];
   	nlc      = [];
    param = list();
   	
   	if (rhs>4) then
   		A	 = varargin(5);
   		b  	 = varargin(6);
   	end

   	if (rhs>6) then
   		Aeq       = varargin(7);
   		beq       = varargin(8);
   	end

   	if (rhs>8) then
   		lb      = varargin(9);
		ub      = varargin(10);
	end

	if (rhs>10) then
   		param = varargin(10);
	end

	if (size(lb,2)==0) then
		lb = repmat(-%inf,nbVar,1);
	end

	if (size(ub,2)==0) then
		ub = repmat(%inf,nbVar,1);
	end

	//Check type of variables
	Checktype("fseminf", _fun, "fun", 1, "function")
	Checktype("fseminf", x0, "x0", 2, "constant")
	Checktype("fseminf", ntheta, "ntheta", 3, "constant")
	Checktype("fseminf", seminfcon, "seminfcon", 4, ["function","constant"])
	Checktype("fseminf", A, "A", 5, "constant")
	Checktype("fseminf", b, "b", 6, "constant")
	Checktype("fseminf", Aeq, "Aeq", 7, "constant")
	Checktype("fseminf", beq, "beq", 8, "constant")
	Checktype("fseminf", lb, "lb", 9, "constant")
	Checktype("fseminf", ub, "ub", 10, "constant")
	Checktype("fseminf", param, "param", 10, "list")
	
	//To check the user entry for options and storing it
   	for i = 1:(size(param))/2
       	select convstr(param(2*i-1),'l')
          	case "maxiter" then
      			Checktype("fseminf", param(2*i), "maxiter", 10, "constant")
  				options(2) = param(2*i);    //Setting the maximum number of iterations as per user entry
       		case "cputime" then
      			Checktype("fseminf", param(2*i), "cputime", 10, "constant")
  				options(4) = param(2*i);    //Setting the maximum CPU time as per user entry
        	case "gradobj" then
      			Checktype("fseminf", param(2*i), "gradobj", 10, "string")
      			if(convstr(param(2*i),'l') == "on") then
      				function dy = graObj(x)
      					
      				endfunction
      			end
       		else
    	     	 	errmsg = msprintf(gettext("%s: Unrecognized parameter name %s."), "fmincon", param(2*i-1));
    	      		error(errmsg);
	     end        					
     end 	      
	
  	//To check and convert the 2nd Input argument (x0) to a row vector 
   	if((size(x0,1)~=1) & (size(x0,2)~=1)) then
   		errmsg = msprintf(gettext("%s: Expected Vector for initial guess"), "fseminf");
   		error(errmsg);
    end

   	if(size(x0,2)==1) then
   		x0=x0(:);
   	end
   	  	
  	//To check the match between fun (1st Parameter) and x0 (2nd Parameter)
   	if(execstr('init=_fun(x0)','errcatch')==21) then
		errmsg = msprintf(gettext("%s: Objective function and x0 did not match"), "fseminf");
   		error(errmsg);
	end
	   	
	//Check the size of inequality constraint which should be equal to the number of variables
	if ( size(A,2) ~= nbVar & size(A,2) ~= 0) then
		errmsg = msprintf(gettext("%s: The number of columns in A must be the same as the number of elements of x0"), "fseminf");
		error(errmsg);
	end
  	
	nbConInEq = size(A,"r");
	  	
	//Check the size of equality constraint which should be equal to the number of variables
	if ( size(Aeq,2) ~= nbVar & size(Aeq,2) ~= 0 ) then
		errmsg = msprintf(gettext("%s: The number of columns in Aeq must be the same as the number of elements of f"), "fseminf");
		error(errmsg);
	end
	
	b = b(:);
	beq = beq(:);  	
	lb = lb(:);
	ub = ub(:);
   	
    //To check the contents of lb & ub
    for i = 1:nbVar
		if(ub(i)<lb(i)) then
			errmsg = msprintf(gettext("%s: Problem has inconsistent variable bounds"), "fseminf");
			error(errmsg);
    	end
	end

	if(typeof(seminfcon) == "function")
		sample_S = %nan
		if(execstr('[sample_c,sample_ceq,sample_K,sample_S] = seminfcon(x0,sample_S)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Semi-Infinite Constraint function & x0 did not match"), "fseminf");
   			error(errmsg);
		end
  		[sample_c, sample_ceq, sample_K, sample_S] = seminfcon(x0,sample_S);

		if (size(sample_c,1)~=1 & size(sample_c,1)~=0) then	
  			errmsg = msprintf(gettext("%s: c in seminfcon should be a row vector or empty matrix"), "fseminf");
    		error(errmsg)
    	end
  		
  		if (size(sample_ceq,1)~=1 & size(sample_ceq,1)~=0) then
  			errmsg = msprintf(gettext("%s: ceq in seminfcon should be a row vector or empty matrix"), "fseminf");
    		error(errmsg)
    	end
		
  		if (size(sample_K,"r")~=ntheta) then
  			errmsg = msprintf(gettext("%s: Number of rows in K should be equal to ntheta"), "fseminf");
    		error(errmsg)
    	end
		
  		if (size(sample_S,"c")~=2) then
  			errmsg = msprintf(gettext("%s: Number of columns in sampling interval should be equal to 2"), "fseminf");
    		error(errmsg)
    	end
	end
	
    ierr = execstr('init=_fun(x0)', "errcatch")
    if ierr <> 0 then
    lamsg = lasterror();
    lclmsg = "%s: Error while evaluating the function: ""%s""\n";
    error(msprintf(gettext(lclmsg), "fseminf", lamsg));
    end

	S = %nan;
    
    ierr = execstr('init=seminfcon(x0,S)', "errcatch")
    if ierr <> 0 then
    lamsg = lasterror();
    lclmsg = "%s: Error while evaluating the function: ""%s""\n";
    error(msprintf(gettext(lclmsg), "fseminf", lamsg));
    end

	function [c, ceq] = _seminfcon(x)
		[c, ceq, K, S] = seminfcon(x,S)
		K_max = max(K,"c");
		c= [c;K_max];
		ceq = ceq;
	endfunction

   	//Calling the fmincon function for solving the above problem
	[xopt,fopt,exitflag,output,lambda,gradient] = fmincon(_fun,x0,A,b,Aeq,beq,lb,ub,_seminfcon,param)
    		
endfunction
