// Check for the size of C and d
C = [
	0.0372    0.2869
	0.6861    0.7071
	0.6233    0.6245
	0.6344    0.6170];
d = [
	0.8587
	0.1781
	0.0747
    0.8405
	0.2356];

// Error
//lsqlin: The number of rows in C must be equal the number of elements of d
//at line     106 of function lsqnonneg called by :  
// [xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg(C,d)

 [xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg(C,d)

