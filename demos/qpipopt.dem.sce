mode(1)
//
// Demo of qpipopt.sci
//

//Find x in R^6 such that:
A= [1,-1,1,0,3,1;
-1,0,-3,-4,5,6;
2,5,3,0,1,0
0,1,0,1,2,-1;
-1,0,2,1,1,0];
conLB=[1;2;3;-%inf;-%inf];
conUB = [1;2;3;-1;2.5];
lb=[-1000;-10000; 0; -1000; -1000; -1000];
ub=[10000; 100; 1.5; 100; 100; 1000];
//and minimize 0.5*x'⋅H⋅x + f'⋅x with
f=[1; 2; 3; 4; 5; 6]; H=eye(6,6);
nbVar = 6;
nbCon = 5;
x0 = repmat(0,nbVar,1);
param = list("MaxIter", 300, "CpuTime", 100);
[xopt,fopt,exitflag,output,lambda]=qpipopt(nbVar,nbCon,H,f,lb,ub,A,conLB,conUB,x0,param)
// Press ENTER to continue
halt()   // Press return to continue
 
//Find the value of x that minimize following function
// f(x) = 0.5*x1^2 + x2^2 - x1*x2 - 2*x1 - 6*x2
// Subject to:
// x1 + x2 ≤ 2
// –x1 + 2x2 ≤ 2
// 2x1 + x2 ≤ 3
// 0 ≤ x1, 0 ≤ x2.
H = [1 -1; -1 2];
f = [-2; -6];
A = [1 1; -1 2; 2 1];
conUB = [2; 2; 3];
conLB = [-%inf; -%inf; -%inf];
lb = [0; 0];
ub = [%inf; %inf];
nbVar = 2;
nbCon = 3;
[xopt,fopt,exitflag,output,lambda] = qpipopt(nbVar,nbCon,H,f,lb,ub,A,conLB,conUB)
//========= E N D === O F === D E M O =========//
