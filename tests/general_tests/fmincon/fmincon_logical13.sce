// Example with objective function and inequality constraints

function y=fun(x)
    y=0
    for i = 1:20
       y = y + exp(x(i)) 
    end
endfunction

x0 = repmat(1,1,20);
lb = repmat(0,1,20);

//Nonlinear constraints
function [c,ceq]=nlc(x)
c = [- 0.5*exp(x(1)) - 2.5*exp(x(2)) - 1.5*exp(x(3)) + 50];
ceq = [];
endfunction

[xopt,fopt,exitflag,output,lambda,gradient,hessian] = fmincon (fun, x0,[],[],[],[],lb,[],nlc)
