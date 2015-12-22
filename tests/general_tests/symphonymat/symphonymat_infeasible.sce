// Infeasible problem
C = -1 * [1 1]'
A = [-1 0; 0, -1; 1 1]
b = [-6 -6 11]
[xopt, fopt, exitflag, output] = symphonymat(C,1,A,b);
