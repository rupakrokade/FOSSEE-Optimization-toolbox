clear
clc
cd /home/siddharth/FOSSEE-Optimization-toolbox
exec cleaner.sce;
ulink;

exec builder.sce;
exec("loader.sce");

///Example 2:
//Objective function to be minimised
function y=f(x)
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
