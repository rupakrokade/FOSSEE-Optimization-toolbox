// qpipopt infeasibility test

H = [2 0;0 8];
f = [0 -32];
A = [-1 0; 0, -1; 1 1];
conUB = [-6 -6 11];
conLB = -1*[%inf %inf %inf];
nbVar = 2;
nbCon = 3;
ub = [%inf %inf];
lb = -1*ub;

// Output
//Converged to a point of local infeasibility.
// lamda  =
// 
//   lower: [0x0 constant]
//   upper: [0x0 constant]
//   constraint: [0x0 constant]
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

[xopt,fopt,exitflag,output,lamda] = qpipopt(nbVar,nbCon,H,f,lb,ub,A,conLB,conUB)
