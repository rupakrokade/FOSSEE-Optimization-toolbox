mode(1)
//
// Demo of intfminunc.sci
//

//Find x in R^2 such that it minimizes the Rosenbrock function
//f = 100*(x2 - x1^2)^2 + (1-x1)^2
//Objective function to be minimised
function y= f(x)
y= 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
endfunction
//Starting point
x0=[-1,2];
intcon = [2]
//Options
options=list("MaxIter", [1500], "CpuTime", [500]);
//Calling
[xopt,fopt,exitflag,gradient,hessian]=intfminunc(f,x0,intcon,options)
// Press ENTER to continue
halt()   // Press return to continue
 
//Find x in R^2 such that the below function is minimum
//f = x1^2 + x2^2
//Objective function to be minimised
function y= f(x)
y= x(1)^2 + x(2)^2;
endfunction
//Starting point
x0=[2,1];
intcon = [1];
[xopt,fopt]=intfminunc(f,x0,intcon)
// Press ENTER to continue
halt()   // Press return to continue
 
//The below problem is an unbounded problem:
//Find x in R^2 such that the below function is minimum
//f = - x1^2 - x2^2
//Objective function to be minimised
function [y,g,h] = f(x)
y = -x(1)^2 - x(2)^2;
g = [-2*x(1),-2*x(2)];
h = [-2,0;0,-2];
endfunction
//Starting point
x0=[2,1];
intcon = [1]
options = list("gradobj","ON","hessian","on");
[xopt,fopt,exitflag,gradient,hessian]=intfminunc(f,x0,intcon,options)
//========= E N D === O F === D E M O =========//
