//Find x in R^6 such that:
// Check if a user gives unequal number of constraints as given by him

A= [1,-1,1,0,3,1;
-1,0,-3,-4,5,6;
2,5,3,0,1,0
0,1,0,1,2,-1;];
b = [1;2;3;-1;2.5];
//and minimize 0.5*x'*H*x + f'*x with
f=[1; 2; 3; 4; 5; 6]; H=eye(6,6);
intcon=[1,2];

//  !--error 10000 
// intqpipopt: The number of rows in A must be the same as the number of elements of b
// at line     291 of function intqpipopt called by :  
// [xopt,fopt,exitflag,output,lambda]=intqpipopt(H,f,intcon,A,b)
// at line      24 of exec file called by :    
// exec intqpipopt_A1.sce

[xopt,fopt,exitflag,output,lambda]=intqpipopt(H,f,intcon,A,b)

