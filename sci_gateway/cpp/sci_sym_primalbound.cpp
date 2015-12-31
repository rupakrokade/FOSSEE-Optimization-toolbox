// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Keyur Joshi
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include "symphony.h"
#include "sci_iofunc.hpp"

extern sym_environment* global_sym_env; //defined in globals.cpp

extern "C" {
#include "api_scilab.h"
#include "Scierror.h"
#include "sciprint.h"
#include "BOOL.h"
#include <localization.h>

int sci_sym_getPrimalBound(char *fname){
	
	//error management variable
	SciErr sciErr;
	int iRet;
	
	//data declarations
	double retVal;
	
	//ensure that environment is active
	if(global_sym_env==NULL){
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
		return 1;
	}
	
	//code to check arguments and get them
	CheckInputArgument(pvApiCtx,0,0) ;
	CheckOutputArgument(pvApiCtx,1,1) ;
	
	//code to process input
	iRet=sym_get_primal_bound(global_sym_env,&retVal);
	if(iRet==FUNCTION_TERMINATED_ABNORMALLY){
		Scierror(999, "An error occured. Has the problem been loaded?\n");
		return 1;
	}
	
	//code to give output
	if(returnDoubleToScilab(retVal))
		return 1;
	
	return 0;
}

int sci_sym_setPrimalBound(char *fname){
	
	//error management variable
	SciErr sciErr;
	int iRet;
	
	//data declarations
	double bound;
	
	//ensure that environment is active
	if(global_sym_env==NULL){
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
		return 1;
	}
	
	//code to check arguments and get them
	CheckInputArgument(pvApiCtx,1,1) ;
	CheckOutputArgument(pvApiCtx,1,1) ;
	
	//code to process input
	if(getDoubleFromScilab(1,&bound))
		return 1;

	iRet=sym_set_primal_bound(global_sym_env,bound);
	if(iRet==FUNCTION_TERMINATED_ABNORMALLY){
		Scierror(999, "An error occured. Has the problem been loaded?\n");
		return 1;
	}
	
	//code to give output
	if(return0toScilab())
		return 1;
	
	return 0;
}

}
