// Check for the size of parameters
C = [
	0.0372    0.2869
	0.6861    0.7071
	0.6233    0.6245
	0.6344    0.6170];
d = [
	0.8587
	0.1781
	0.0747
    0.8405];

param = list("MaxIter");

//Error
//lsqlin: Size of parameters should be even
//at line      75 of function lsqnonneg called by :  
//[xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg(C,d,param)

[xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg(C,d,param)

