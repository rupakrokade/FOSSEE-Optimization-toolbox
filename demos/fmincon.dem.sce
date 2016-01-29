mode(1)
//
// Demo of fmincon.sci
//

//Find x in R^2 such that it minimizes:
//f(x)= -x1 -x2/3
//x0=[0,0]
//constraint-1 (c1): x1 + x2 <= 2
//constraint-2 (c2): x1 + x2/4 <= 1
//constraint-3 (c3): x1 - x2 <= 2
//constraint-4 (c4): -x1/4 - x2 <= 1
//constraint-5 (c5): -x1 - x2 <= -1
//constraint-6 (c6): -x1 + x2 <= 2
//constraint-7 (c7): x1 + x2 = 2
//Objective function to be minimised
function y=f(x)
y=-x(1)-x(2)/3;
endfunction
//Starting point, linear constraints and variable bounds
x0=[0 , 0];
A=[1,1 ; 1,1/4 ; 1,-1 ; -1/4,-1 ; -1,-1 ; -1,1];
b=[2;1;2;1;-1;2];
Aeq=[1,1];
beq=[2];
lb=[];
ub=[];
nlc=[];
//Gradient of objective function
function y= fGrad(x)
y= [-1,-1/3];
endfunction
//Hessian of lagrangian
function y= lHess(x,obj,lambda)
y= obj*[0,0;0,0]
endfunction
//Options
options=list("GradObj", fGrad, "Hessian", lHess);
//Calling Ipopt
[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options)
// Press ENTER to continue
halt()   // Press return to continue
 
//Find x in R^3 such that it minimizes:
//f(x)= x1*x2 + x2*x3
//x0=[0.1 , 0.1 , 0.1]
//constraint-1 (c1): x1^2 - x2^2 + x3^2 <= 2
//constraint-2 (c2): x1^2 + x2^2 + x3^2 <= 10
//Objective function to be minimised
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
//Gradient of objective function
function y= fGrad(x)
y= [x(2),x(1)+x(3),x(2)];
endfunction
//Hessian of the Lagrange Function
function y= lHess(x,obj,lambda)
y= obj*[0,1,0;1,0,1;0,1,0] + lambda(1)*[2,0,0;0,-2,0;0,0,2] + lambda(2)*[2,0,0;0,2,0;0,0,2]
endfunction
//Gradient of Non-Linear Constraints
function [cg,ceqg] = cGrad(x)
cg=[2*x(1) , -2*x(2) , 2*x(3) ; 2*x(1) , 2*x(2) , 2*x(3)];
ceqg=[];
endfunction
//Options
options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", fGrad, "Hessian", lHess,"GradCon", cGrad);
//Calling Ipopt
[x,fval,exitflag,output] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options)
// Press ENTER to continue
halt()   // Press return to continue
 
//The below problem is an unbounded problem:
//Find x in R^3 such that it minimizes:
//f(x)= -(x1^2 + x2^2 + x3^2)
//x0=[0.1 , 0.1 , 0.1]
//  x1 <= 0
//  x2 <= 0
//  x3 <= 0
//Objective function to be minimised
function y=f(x)
y=-(x(1)^2+x(2)^2+x(3)^2);
endfunction
//Starting point, linear constraints and variable bounds
x0=[0.1 , 0.1 , 0.1];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[];
ub=[0,0,0];
//Options
options=list("MaxIter", [1500], "CpuTime", [500]);
//Calling Ipopt
[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,[],options)
// Press ENTER to continue
halt()   // Press return to continue
 
//The below problem is an infeasible problem:
//Find x in R^3 such that in minimizes:
//f(x)=x1*x2 + x2*x3
//x0=[1,1,1]
//constraint-1 (c1): x1^2 <= 1
//constraint-2 (c2): x1^2 + x2^2 <= 1
//constraint-3 (c3): x3^2 <= 1
//constraint-4 (c4): x1^3 = 0.5
//constraint-5 (c5): x2^2 + x3^2 = 0.75
// 0 <= x1 <=0.6
// 0.2 <= x2 <= inf
// -inf <= x3 <= 1
//Objective function to be minimised
function y=f(x)
y=x(1)*x(2)+x(2)*x(3);
endfunction
//Starting point, linear constraints and variable bounds
x0=[1,1,1];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[0 0.2,-%inf];
ub=[0.6 %inf,1];
//Nonlinear constraints
function [c,ceq]=nlc(x)
c=[x(1)^2-1,x(1)^2+x(2)^2-1,x(3)^2-1];
ceq=[x(1)^3-0.5,x(2)^2+x(3)^2-0.75];
endfunction
//Gradient of objective function
function y= fGrad(x)
y= [x(2),x(1)+x(3),x(2)];
endfunction
//Hessian of the Lagrange Function
function y= lHess(x,obj,lambda)
y= obj*[0,1,0;1,0,1;0,1,0] + lambda(1)*[2,0,0;0,0,0;0,0,0] + lambda(2)*[2,0,0;0,2,0;0,0,0] +lambda(3)*[0,0,0;0,0,0;0,0,2] + lambda(4)*[6*x(1    ),0,0;0,0,0;0,0,0] + lambda(5)*[0,0,0;0,2,0;0,0,2];
endfunction
//Gradient of Non-Linear Constraints
function [cg,ceqg] = cGrad(x)
cg = [2*x(1),0,0;2*x(1),2*x(2),0;0,0,2*x(3)];
ceqg = [3*x(1)^2,0,0;0,2*x(2),2*x(3)];
endfunction
//Options
options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", fGrad, "Hessian", lHess,"GradCon", cGrad);
//Calling Ipopt
[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options)
//========= E N D === O F === D E M O =========//
