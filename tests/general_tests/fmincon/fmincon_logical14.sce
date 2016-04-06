// Example with objective function and inequality constraints

function y=fun(x)
    y=0
    for i = 1:20
       y = y + exp(x(i)) 
    end
endfunction

x0 = repmat(1,1,20);
lb = repmat(0,1,20);

A=[-1,-5,-3 repmat(0,1,17); -0.5,-2.5 -1.5 repmat(0,1,17);];
b=[-100 -50]';

//Nonlinear constraints
function [c,ceq]=nlc(x)
    cfor = 0;
    for i = 1:20
        cfor = cfor + 2*exp(x(i)) 
    end
    c = [ cfor + 1];
    ceq = [];
endfunction

[xopt,fopt,exitflag,output,lambda,gradient,hessian] = fmincon (fun, x0,A,b,[],[],lb,[],nlc)
