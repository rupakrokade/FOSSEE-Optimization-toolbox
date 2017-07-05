//
// Demo of cbcintlinprog.sci
//


// Objective function
    c = [-750 -1000];
    // Constraint inequality matrix LHS
    G = [
      1 1;
      1 2;
      4 3;
        ];
    // Constraint inequality matrix RHS
    h = [10 15 25]';
    // Constraint equality matrix LHS
    A = [
      0.5 1
    ];
    // Constraint equality matrix RHS
    b=[7.5];
    // Dimension of positive orthant
    l = [3];
    q = [];
    e = [];
    dims=list("l",l,"q",q,"e",e)
     //Calling ecos
    [x,y,s,z,info,status] =ecos(c,G,h,dims,A,b);


halt()   // Press return to continue

    // Objective function
    c = [0 0 0 0 1];
     //Constraint inequality matrix LHS
    G = [
      0.4167578    0.0562668    0.           0.           0.  
        2.1361961   -1.6402708    0.           0.           0.  
        1.7934356    0.8417474    0.           0.           0.  
        0.           0.           0.4167578    0.0562668    0.  
        0.           0.           2.1361961   -1.6402708    0.  
        0.           0.           1.7934356    0.8417474    0.  
        0.           0.           0.           0.          -1.  
       -1.           0.           0.           0.           0.  
        0.          -1.           0.           0.           0.  
        0.           0.          -1.           0.           0.  
        0.           0.           0.          -1.           0.
        ];
     //Constraint inequality matrix RHS
    h = [0 0  0  0  0  0  0  0  0  0  0]';
    // Dimension of positive orthant
    l = [6];
    q = [5];
    e = [0]
    dims=list("l",l,"q",q,"e",e)
    [x,y,s,z,info,status] =ecos(c,G,h,dims);