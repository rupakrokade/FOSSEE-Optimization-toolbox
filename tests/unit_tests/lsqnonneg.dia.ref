// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// Author: Harpreet Singh
// Organization: FOSSEE, IIT Bombay
// Email: harpreet.mertia@gmail.com
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

// <-- JVM NOT MANDATORY -->
// <-- ENGLISH IMPOSED -->


//
// assert_close --
//   Returns 1 if the two real matrices computed and expected are close,
//   i.e. if the relative distance between computed and expected is lesser than epsilon.
// Arguments
//   computed, expected : the two matrices to compare
//   epsilon : a small number
//
function flag = assert_close ( computed, expected, epsilon )
  if expected==0.0 then
    shift = norm(computed-expected);
  else
    shift = norm(computed-expected)/norm(expected);
  end
//  if shift < epsilon then
//    flag = 1;
//  else
//    flag = 0;
//  end
//  if flag <> 1 then pause,end
    flag = assert_checktrue ( shift < epsilon );
endfunction
//
// assert_equal --
//   Returns 1 if the two real matrices computed and expected are equal.
// Arguments
//   computed, expected : the two matrices to compare
//   epsilon : a small number
//
//function flag = assert_equal ( computed , expected )
//  if computed==expected then
//    flag = 1;
//  else
//    flag = 0;
//  end
//  if flag <> 1 then pause,end
//endfunction

// A basic lsqnonneg problem
C = [
	0.0372    0.2869
	0.6861    0.7071
	0.6233    0.6245
	0.6344    0.6170];
d = [
	0.8587
	0.1781
	0.0747
	0.8405];
[xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg(C,d)

assert_close ( xopt , [ 0 0.6929344 ]' , 0.0005 );
assert_close ( residual , [0.6598971 -0.3118739 -0.3580375 0.4129595]' , 0.0005 );
assert_close ( resnorm , [ 0.8314560] , 0.0005 );
assert_checkequal( exitflag , int32(0) );
printf("Test Successful");
