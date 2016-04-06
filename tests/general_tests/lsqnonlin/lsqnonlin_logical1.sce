function retF = testmyfun(x)
	km = [1:10]';
	retF = 2 + 2*km-exp(km*x(1))-exp(km*x(2));
endfunction

x0 = [0.3 0.4]'
[x,resnorm] = lsqnonlin(testmyfun,x0) 
