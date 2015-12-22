
// A simple example without constraints
conMatrix= [];
conLB=[];
conUB = [];
lb=[];
ub=[];
p=[2 -35 -47]'; 
Q =[5   -2   -1; -2   4   3; -1   3   5];
nbVar = 3;
nbCon = 0;
[xopt,fopt,exitflag,output,lambda]=qpipopt(nbVar,nbCon,Q,p,lb,ub,conMatrix,conLB,conUB)
disp("xopt",xopt,"fopt",fopt,"exitflag",exitflag,"output",output,"lambda",lambda)


//Output
//
//Optimal Solution Found.
// 
//   lower: [0,0,0]
//   upper: [0,0,0]
//   constraint: [0x0 constant]
// 
// lambda   
// 
//   Iterations: 1
// 
// output   
// 
//  0  
// 
// exitflag   
// 
//  - 249.  
// 
// fopt   
// 
//    3.  
//    5.  
//    7.  
// 
// xopt 
