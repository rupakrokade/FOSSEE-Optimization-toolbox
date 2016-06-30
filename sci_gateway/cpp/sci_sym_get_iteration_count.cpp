// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Sai Kiran
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include <symphony.h>
#include <sci_iofunc.hpp>
extern sym_environment* global_sym_env;//defined in globals.cpp

extern "C" {
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <stdlib.h>
#include <malloc.h>
#include <localization.h>
#include <sciprint.h>

#include <string.h>

/*
 * This function is used to get iteration count after solving a problem
*/
int sci_sym_get_iteration_count(char *fname, unsigned long fname_len){
	
	//check whether we have no input and one output argument or not
	CheckInputArgument(pvApiCtx, 0, 0) ; //no input argument
	CheckOutputArgument(pvApiCtx, 1, 1) ; //one output argument
	
	int iteration_count=0; // return value to the caller
	if(global_sym_env==NULL) //There is no environment opened.
		sciprint("Error: Symphony environment is not initialized.\n");
	else { //There is an environment opened
		 //Call symphony function
		int status=sym_get_iteration_count(global_sym_env,&iteration_count);
		if (status == FUNCTION_TERMINATED_ABNORMALLY) {
			sciprint("\nHave you solved a problem ?\n");
			iteration_count = 0;
			}
		}
	// Write the result to scilab
	return returnDoubleToScilab(iteration_count);
	}

}
