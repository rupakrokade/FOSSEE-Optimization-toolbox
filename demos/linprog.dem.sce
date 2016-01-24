mode(1)
//
// Demo of linprog.sci
//

//Optimal problems
//Linear program, linear inequality constraints
c=[-1,-1/3]'
A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1]
b=[2,1,2,1,-1,2]
[xopt,fopt,exitflag,output,lambda]=linprog(c, A, b)
// Press ENTER to continue
halt()   // Press return to continue
 
//Linear program with Linear Inequalities and Equalities`
c=[-1,-1/3]'
A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1]
b=[2,1,2,1,-1,2]
Aeq=[1,1/4]
beq=[1/2]
[xopt,fopt,exitflag,output,lambda]=linprog(c, A, b, Aeq, beq)
// Press ENTER to continue
halt()   // Press return to continue
 
//Linear program with all constraint types
c=[-1,-1/3]'
A=[1,1;1,1/4;1,-1;-1/4,-1;-1,-1;-1,1]
b=[2,1,2,1,-1,2]
Aeq=[1,1/4]
beq=[1/2]
lb=[-1,-0.5]
ub=[1.5,1.25]
[xopt,fopt,exitflag,output,lambda]=linprog(c, A, b, Aeq, beq, lb, ub)
// Press ENTER to continue
halt()   // Press return to continue
 
//Primal Infeasible Problem
c=[-1,-1,-1]'
A=[1,2,-1]
b=[-4]
Aeq=[1,5,3;1,1,0]
beq=[10,100]
lb=[0,0,0]
ub=[%inf,%inf,%inf]
[xopt,fopt,exitflag,output,lambda]= linprog(c,A,b,Aeq,beq,lb,ub)
// Press ENTER to continue
halt()   // Press return to continue
 
//Dual Infeasible Problem
c=[3,5,-7]'
A=[-1,-1,4;1,1,4]
b=[-8,5]
Aeq=[]
beq=[]
lb=[-%inf,-%inf,-%inf]
ub=[%inf,%inf,%inf]
[xopt,fopt,exitflag,output,lambda]= linprog(c,A,b,Aeq,beq,lb,ub)
// Press ENTER to continue
halt()   // Press return to continue
 
filepath = get_absolute_file_path('linprog.dem.sce');
filepath = filepath + "exmip1.mps"
[xopt,fopt,exitflag,output,lambda] =linprog(filepath);
//========= E N D === O F === D E M O =========//
