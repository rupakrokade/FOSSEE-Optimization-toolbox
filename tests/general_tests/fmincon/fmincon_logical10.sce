//Example where user provides hessian

function y=f(x)
y=x(1)*x(2)+x(2)*x(3);
endfunction
//Starting point, linear constraints and variable bounds
x0=[0.1 , 0.1 , 0.1];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[];
ub=[];
//Nonlinear constraints
function [c,ceq]=nlc(x)
c = [x(1)^2 - x(2)^2 + x(3)^2 - 2 , x(1)^2 + x(2)^2 + x(3)^2 - 10];
ceq = [];
endfunction

//Hessian of the Lagrange Function
function y= lHess(x,obj,lambda)
y= obj*[0,1,0;1,0,1;0,1,0] + lambda(1)*[2,0,0;0,-2,0;0,0,2] + lambda(2)*[2,0,0;0,2,0;0,0,2]
endfunction

//Options
options=list("MaxIter", [1500], "CpuTime", [500], "Hessian", lHess);

//Output
//Optimal Solution Found.
// hessian  =
// 
//    1.4142136    1.72D-322    2.12D-314  
//    1.           0.           5.82D+252  
//    1.           1.4142136    1.         
// gradient  =
// 
//    2.236068  - 3.1622776    2.236068  
// lambda  =
// 
//   lower: [0,0,0]
//   upper: [0,0,0]
//   ineqlin: [0x0 constant]
//   eqlin: [0x0 constant]
//   ineqnonlin: [4.545D-09,0.7071068]
//   eqnonlin: [0x0 constant]
// output  =
// 
//   Iterations: 23
//   Cpu_Time: 0.164
//   Objective_Evaluation: 24
//   Dual_Infeasibility: 6.124D-08
//   Message: "Optimal Solution Found"
// exitflag  =
// 
//  0  
// fopt  =
// 
//  - 7.0710678  
// xopt  =
// 
//  - 1.5811388  
//    2.236068   
//  - 1.5811388 

//Calling Ipopt
[xopt,fopt,exitflag,output,lambda,gradient,hessian] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options)
