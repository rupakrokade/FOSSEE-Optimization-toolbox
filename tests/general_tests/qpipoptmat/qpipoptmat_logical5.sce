//Find x in R^6 such that:
// A simple example with constraints

Aeq= [1,-1,1,0,3,1;
-1,0,-3,-4,5,6;
2,5,3,0,1,0];
beq=[1; 2; 3];
A= [0,1,0,1,2,-1;
-1,0,2,1,1,0];
b = [-1; 2.5];
lb=[-1000; -10000; 0; -1000; -1000; -1000];
ub=[10000; 100; 1.5; 100; 100; 1000];
param = list("MaxIter", 300, "CpuTime",100);
//and minimize 0.5*x'*H*x + f'*x with
f=[1; 2; 3; 4; 5; 6]; H=eye(6,6);
x0 = repmat(0,6,1);
//Output
//
//Optimal Solution Found.
// lambda  =
// 
//   lower: [0,0,0,0,0,0]
//   upper: [0,0,0,0,0,0]
//   eqlin: [0x0 constant]
//   ineqlin: [1.232D-12,6.884D-13]
// output  =
// 
//   Iterations: 5
//   ConstrViolation: 0
// exitflag  =
// 
//  0  
// fopt  =
// 
//  - 45.5  
// xopt  =
// 
//  - 1.  
//  - 2.  
//  - 3.  
//  - 4.  
//  - 5.  
//  - 6. 

[xopt,fopt,exitflag,output,lambda]=qpipoptmat(H,f,A,b)
   
