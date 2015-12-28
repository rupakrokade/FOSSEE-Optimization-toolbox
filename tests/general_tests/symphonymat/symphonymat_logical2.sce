// An example with equality constraints and variable bounds

// A basic case :

// Objective function
c = -1*[20,10,15]';

// Lower Bound of variable
lb = repmat(0,3,1);

// Upper Bound of variables
ub = repmat(%inf,3,1);

// Constraint Matrix
A = [3,2,5;
     2,1,1;
     1,1,3;
     5,2,4]

// Upper Bound of constrains
b = [ 55;26;30;57]

// Row Matrix for telling symphony that the is integer or not
intcon = [];

// Output
//Problem loaded into environment.
//
//Note: There is no limit on time.
//
//An optimal solution has been found.
// 
//    0.  
// 
//   Iterations: 1
// 
// output   
// 
//    227.  
// 
// status   
// 
//  - 268.  
// 
// f   
// 
//    1.8   
//    20.8  
//    1.6   
// 
// x  

// Calling Symphony
[x,f,status,output] = symphonymat(c,intcon,A,b,[],[],lb,ub)
disp("x",x,"f",f,"status",status,"output",output);

