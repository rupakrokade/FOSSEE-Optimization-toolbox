// An example with C and d as input
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

//Output
//Optimal Solution Found.
// lambda  =
// 
//   lower: [0.1506118,1.441D-11]
//   upper: [0,0]
// output  =
// 
//   Iterations: 5
// exitflag  =
// 
//  0  
// residual  =
// 
//    0.6598971  
//  - 0.3118739  
//  - 0.3580375  
//    0.4129595  
// resnorm  =
// 
//    0.8314560  
// xopt  =
// 
//    0.         
//    0.6929344  

[xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg(C,d)

