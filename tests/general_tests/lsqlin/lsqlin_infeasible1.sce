// Check for the infeasible problem
C = [0.9501    0.7620    0.6153    0.4057
	 0.2311    0.4564    0.7919    0.9354
	 0.6068    0.0185    0.9218    0.9169
	 0.4859    0.8214    0.7382    0.4102
	 0.8912    0.4447    0.1762    0.8936];
d = [0.0578
	 0.3528
	 0.8131
	 0.0098
	 0.1388];
Aeq = [-1 0 0 0; 0 -1 0 0; 1 1 0 0 ];
beq = [-6 -6 11];

//Converged to a point of local infeasibility.
// lambda  =
// 
//   lower: [0x0 constant]
//   upper: [0x0 constant]
//   eqlin: [0x0 constant]
//   ineqlin: [0x0 constant]
// output  =
// 
//   Iterations: 0
// exitflag  =
// 
//  5  
// residual  =
// 
//    0.0578  
//    0.3528  
//    0.8131  
//    0.0098  
//    0.1388  
// resnorm  =
// 
//    0.8083018  
// xopt  =
// 
//     []

[xopt,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,[],[],Aeq,beq)

