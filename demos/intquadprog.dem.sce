mode(1)
//
// Demo of intquadprog.sci
//

H = [1 -1; -1 2];
f = [-2; -6];
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];
lb=[0,0];
ub=[%inf, %inf];
intcon = [1 2];
[xopt,fopt,status,output]=intquadprog(H,f,intcon,A,b,[],[],lb,ub)
halt()   // Press return to continue
 
//Find x in R^6 such that:
Aeq= [1,-1,1,0,3,1;
-1,0,-3,-4,5,6;
2,5,3,0,1,0];
beq=[1; 2; 3];
A= [0,1,0,1,2,-1;
-1,0,2,1,1,0];
b = [-1; 2.5];
lb=[-1000; -10000; 0; -1000; -1000; -1000];
ub=[10000; 100; 1.5; 100; 100; 1000];
x0 = repmat(0,6,1);
param = list("MaxIter", 300, "CpuTime", 100);
//and minimize 0.5*x'*H*x + f'*x with
f=[1; 2; 3; 4; 5; 6]; H=eye(6,6);
intcon = [2 4];
[xopt,fopt,exitflag,output]=intquadprog(H,f,intcon,A,b,Aeq,beq,lb,ub,x0,param)
//========= E N D === O F === D E M O =========//
