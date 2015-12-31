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

int sci_sym_isEnvActive(char *fname){
	
	//error management variable
	SciErr sciErr;
	int iRet;
	
	//data declarations
	double returnVal;
	
	//code to check arguments and get them
	CheckInputArgument(pvApiCtx,0,0) ;
	CheckOutputArgument(pvApiCtx,1,1) ;
	
	//code to process input
	if(global_sym_env==NULL){
		sciprint("Symphony environment is not initialized. Please run 'sym_open()' first.\n");
		returnVal=0.0;
	}else{
		sciprint("Symphony environment is active and ready for use.\n");
		returnVal=1.0;
	}
	
	//code to give output
	if(returnDoubleToScilab(returnVal))
		return 1;

	return 0;
}

}
