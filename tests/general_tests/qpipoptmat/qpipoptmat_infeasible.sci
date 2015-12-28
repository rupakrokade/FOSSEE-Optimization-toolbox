// qpipopt infeasibility test

H = [2 0;0 8];
f = [0 -32];
A = [-1 0; 0, -1; 1 1];
b = [-6 -6 11];
ub = [%inf %inf];
lb = -1*ub;
[xopt,fopt,exitflag,output,lamda] = qpipoptmat(H,f,A,b,[],[],lb,ub)

//Output
//Converged to a point of local infeasibility.
// lamda  =
// 
//   lower: [0x0 constant]
//   upper: [0x0 constant]
//   eqlin: [0x0 constant]
//   ineqlin: [0x0 constant]
// output  =
// 
//   Iterations: 0
// exitflag  =
// 
//  5  
// fopt  =
// 
//    0.  
// xopt  =
// 
//     []
