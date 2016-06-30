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

/* This is generelized function for 
 * sym_getVarLower,sym_getVarUpper,sym_getRhs,sym_getConstrRange,sym_getConstrLower,
 * sym_getConstrUpper and sym_getObjCoeff .
 * (Functions taking symphony env and pointer to array of doubles as arguments)
*/
int sci_sym_get_dbl_arr(char *fname, unsigned long fname_len){

	int result_len=0;/* Length of the output double array */
	double *result=NULL;/* Pointer to output double array */
  
	//check whether we have no input and one output argument or not
	CheckInputArgument(pvApiCtx, 0, 0) ; //no input argument
	CheckOutputArgument(pvApiCtx, 1, 1) ; //one output argument
	
	/* Array of possible callers of this function */
	char* arr_caller[]={"sym_getVarLower","sym_getVarUpper",
						"sym_getRhs","sym_getConstrRange",
						"sym_getConstrLower","sym_getConstrUpper",
						"sym_getObjCoeff"};

	/* Array of functions to be called */
	int (*fun[])(sym_environment*,double*)= {sym_get_col_lower,sym_get_col_upper,
											 sym_get_rhs,sym_get_row_range,
											 sym_get_row_lower,sym_get_row_upper,
											sym_get_obj_coeff };
	
	/* Array of functions the above functions depend on */
	int (*fun_depends[])(sym_environment*,int*) = {sym_get_num_cols,sym_get_num_cols,
													sym_get_num_rows,sym_get_num_rows,
													sym_get_num_rows,sym_get_num_rows,
													sym_get_num_cols };

	/* We want to ouput row-matrix if we are dealing with column data .
	 * column matrix if we are dealing with row data .
	 * 0 - output a row matrix.
	 * 1 - output a column matrix.
	 */
	int representation = 0; //output a row matrix

	/* Array of representations of output depending on the above functions.
	 * It's length is same as above arrays.
	*/
	int matrix_representation[] = { 0 ,0 , 1, 1, 1, 1, 0};
	if(global_sym_env==NULL) //There is no environment opened.
		sciprint("Error: Symphony environment is not initialized.\n");
	else {
		//There is an environment opened
		int iter=0,length=sizeof(arr_caller)/sizeof(char*),found_at= -1;
		for (;iter < length ;++iter){
			if (!strcmp(fname,arr_caller[iter]))
				found_at=iter;
			}
		if (found_at != -1){
			int status1=fun_depends[found_at](global_sym_env,&result_len);
			if ( status1 == FUNCTION_TERMINATED_NORMALLY && result_len ) {
				result=(double*)malloc( sizeof(double) * result_len );
				int ret_val=fun[found_at](global_sym_env,result);
				sciprint("\nFunction invoked unsuccessfully.\n");
				if (ret_val == FUNCTION_TERMINATED_ABNORMALLY)
					result_len=0;
				else {
					if (found_at == 6) {//if called function is sym_getObjCoeff
						int iter=0,sense=0,status2 = sym_get_obj_sense(global_sym_env,&sense);
						if (sense == -1) // Multiply with -1 while showing 
							for (;iter < result_len;++iter) result[iter] *= -1; 						
						}
					representation = matrix_representation[found_at];
					}				
				}
			else
				sciprint("\n Is a problem loaded ? \n");
			}
		else //very rare case
			sciprint("\nError in function mapping in scilab script\n");
		}

	//Copy the result to scilab. Location is position next to input arguments.
	SciErr err;
	if (representation) // output a column-matrix
		err=createMatrixOfDouble(pvApiCtx,nbInputArgument(pvApiCtx)+1,result_len,1,result);
	else // output a row-matrix
		err=createMatrixOfDouble(pvApiCtx,nbInputArgument(pvApiCtx)+1,1,result_len,result);
	free(result); //Free the allocated space
	result=NULL; //Set to NULL
	if (err.iErr){ //Process error
		AssignOutputVariable(pvApiCtx, 1) = 0;
        printError(&err, 0);
        return 1;
    	}

	//assign result position to output argument
	AssignOutputVariable(pvApiCtx, 1) = nbInputArgument(pvApiCtx) + 1;
	//ReturnArguments(pvApiCtx);
	return 0;
	}
}
