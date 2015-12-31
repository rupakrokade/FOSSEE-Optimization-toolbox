mode(1)
//
// Demo of lsqlin.sci
//

//A simple linear least square example
C = [ 2 0;
-1 1;
0 2]
d = [1
0
-1];
A = [10 -2;
-2 10];
b = [4
-4];
[xopt,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,A,b)
// Press ENTER to continue
halt()   // Press return to continue
 
//A basic example for equality, inequality constraints and variable bounds
C = [1 1 1;
1 1 0;
0 1 1;
1 0 0;
0 0 1]
d = [89;
67;
53;
35;
20;]
A = [3 2 1;
2 3 4;
1 2 3];
b = [191
209
162];
Aeq = [1 2 1];
beq = 10;
lb = repmat(0.1,3,1);
ub = repmat(4,3,1);
[xopt,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,A,b,Aeq,beq,lb,ub)
//========= E N D === O F === D E M O =========//
