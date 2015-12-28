// Check for the input arguments
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

// Error
//lsqlin: Unexpected number of input arguments : 1 provided while should be in the set of [2 3]
//at line      55 of function lsqnonneg called by :  
// [xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg(C)

 [xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg(C)

