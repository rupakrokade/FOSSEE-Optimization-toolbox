mode(1)
//
// Demo of lsqnonlin.sci
//

//A simple non-linear least square example taken from leastsq default present in scilab
function y=yth(t, x)
y  = x(1)*exp(-x(2)*t)
endfunction
// we have the m measures (ti, yi):
m = 10;
tm = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5]';
ym = [0.79, 0.59, 0.47, 0.36, 0.29, 0.23, 0.17, 0.15, 0.12, 0.08]';
// measure weights (here all equal to 1...)
wm = ones(m,1);
// and we want to find the parameters x such that the model fits the given
// data in the least square sense:
//
//  minimize  f(x) = sum_i  wm(i)^2 ( yth(tm(i),x) - ym(i) )^2
// initial parameters guess
x0 = [1.5 ; 0.8];
// in the first examples, we define the function fun and dfun
// in scilab language
function y=myfun(x, tm, ym, wm)
y = wm.*( yth(tm, x) - ym )
endfunction
// the simplest call
[xopt,resnorm,residual,exitflag,output,lambda,gradient] = lsqnonlin(myfun,x0)
// Press ENTER to continue
halt()   // Press return to continue
 
//A basic example taken from leastsq default present in scilab with gradient
function y=yth(t, x)
y  = x(1)*exp(-x(2)*t)
endfunction
// we have the m measures (ti, yi):
m = 10;
tm = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5]';
ym = [0.79, 0.59, 0.47, 0.36, 0.29, 0.23, 0.17, 0.15, 0.12, 0.08]';
// measure weights (here all equal to 1...)
wm = ones(m,1);
// and we want to find the parameters x such that the model fits the given
// data in the least square sense:
//
//  minimize  f(x) = sum_i  wm(i)^2 ( yth(tm(i),x) - ym(i) )^2
// initial parameters guess
x0 = [1.5 ; 0.8];
// in the first examples, we define the function fun and dfun
// in scilab language
function [y,dy]=myfun(x, tm, ym, wm)
y = wm.*( yth(tm, x) - ym )
v = wm.*exp(-x(2)*tm)
dy = [v , -x(1)*tm.*v]
endfunction
options = list("GradObj", "on")
[xopt,resnorm,residual,exitflag,output,lambda,gradient] = lsqnonlin(myfun,x0,[],[],options)
//========= E N D === O F === D E M O =========//
