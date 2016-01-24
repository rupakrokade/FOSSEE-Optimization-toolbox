mode(1)
//
// Demo of fminbnd.sci
//

//Find x in R^6 such that it minimizes:
//f(x)= sin(x1) + sin(x2) + sin(x3) + sin(x4) + sin(x5) + sin(x6)
//-2 <= x1,x2,x3,x4,x5,x6 <= 2
//Objective function to be minimised
function y=f(x)
y=0
for i =1:6
y=y+sin(x(i));
end
endfunction
//Variable bounds
x1 = [-2, -2, -2, -2, -2, -2];
x2 = [2, 2, 2, 2, 2, 2];
//Options
options=list("MaxIter",[1500],"CpuTime", [100],"TolX",[1e-6])
//Calling Ipopt
[x,fval] =fminbnd(f, x1, x2, options)
halt()   // Press return to continue
 
//Find x in R such that it minimizes:
//f(x)= 1/x^2
//0 <= x <= 1000
//Objective function to be minimised
function y=f(x)
y=1/x^2
endfunction
//Variable bounds
x1 = [0];
x2 = [1000];
//Calling Ipopt
[x,fval,exitflag,output,lambda] =fminbnd(f, x1, x2)
halt()   // Press return to continue
 
//The below problem is an unbounded problem:
//Find x in R^2 such that it minimizes:
//f(x)= -[(x1-1)^2 + (x2-1)^2]
//-inf <= x1,x2 <= inf
//Objective function to be minimised
function y=f(x)
y=-((x(1)-1)^2+(x(2)-1)^2);
endfunction
//Variable bounds
x1 = [-%inf , -%inf];
x2 = [];
//Options
options=list("MaxIter",[1500],"CpuTime", [100],"TolX",[1e-6])
//Calling Ipopt
[x,fval,exitflag,output,lambda] =fminbnd(f, x1, x2, options)
//========= E N D === O F === D E M O =========//
