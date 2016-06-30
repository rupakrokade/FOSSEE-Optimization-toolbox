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
 * Generelized function for sym_getNumVars,
 * sym_getNumConstrs,sym_get_NumElements
*/
int sci_sym_get_num_int(char *fname, unsigned long fname_len){

	int result=-1;/* Result of the callar */
  
	//check whether we have no input and one output argument or not
	CheckInputArgument(pvApiCtx, 0, 0) ; //no input argument
	CheckOutputArgument(pvApiCtx, 1, 1) ; //one output argument

	/* Array of possible callers of this function */
	char* arr_caller[]={"sym_getNumConstr","sym_getNumVar","sym_getNumElements"};

	/* Array of functions to be called */
	int (*fun[])(sym_environment*,int*)= { sym_get_num_rows,
											sym_get_num_cols,
											sym_get_num_elements
											};
	
	if(global_sym_env==NULL) //There is no environment opened.
		sciprint("Error: Symphony environment is not initialized.\n");
	else {
		//There is an environment opened
		int iter=0,length=sizeof(arr_caller)/sizeof(char*),found_at= -1;
		for (;iter < length ;++iter){
			if (!strcmp(fname,arr_caller[iter]))
				found_at=iter;
			}
		if (found_at != -1) {
			int ret_val=fun[found_at](global_sym_env,&result);
			sciprint("\nFunction invoked unsuccessfully.\n");
			if (ret_val == FUNCTION_TERMINATED_ABNORMALLY)
				result=-1;
			}
		else //very rare case
			sciprint("\nError in function mapping in scilab script\n");
		}
	
	//Copy the result to scilab. Location is position next to input arguments.
	return returnDoubleToScilab(result);
	}
}