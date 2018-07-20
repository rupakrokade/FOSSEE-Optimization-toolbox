mode(1)
//
// Demo of lsqnonneg.sci
//

// A basic lsqnonneg problem
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
[xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg(C,d)
//========= E N D === O F === D E M O =========//
