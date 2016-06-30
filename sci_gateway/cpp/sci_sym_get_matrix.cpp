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
 * Proto-type of function that converts column-major (sparse) representation
 * to row-major (sparse) representation .
*/
void column_major_to_row_major(int,int,int,double *,int *,int *,double *,int *,int *);


/* This function is to retrieve the problem's constraint matrix (sparse) .
 * Symphony uses column-major (sparse) representation.
 * Scilab uses row-major (sparse) representation.
 * So, This function takes column-major (sparse) representation from symphony ,
 * converts that to row-major (sparse) representation and writes to scilab's memory.
 *
 **/
int sci_sym_get_matrix(char *fname, unsigned long fname_len){
	int nz_ele=0;// No.of non-zero elements of the matrix
	int rows=0; //No. of rows in constraint matrix
	int columns=0; //No. of columns in constraint matrix

	/* Variables that store column-major representation of matrix.
	 * These variables will be filled by symphony
	*/
	int *column_start=NULL;// Starting index(in elements array) of each column 
	int *row_indices=NULL;// Row indices corresponding to each non-zero element
	double *elements=NULL;// Non-zero elements of matrix
	
	/* Variables that store row-major representation of matrix.
	 * Filled by a function column_major_to_row_major.
	*/
	double *new_list=NULL; // Non-zero elements of row-major representation
	int *count_per_row=NULL; //Count of non-zero elements in earch row
	int *column_position=NULL; //Column of each non-zero element

	//check whether we have no input and one output argument or not
	CheckInputArgument(pvApiCtx, 0, 0) ; //no input argument
	CheckOutputArgument(pvApiCtx, 1, 1) ; //one output argument

	if(global_sym_env==NULL) //There is no environment opened.
		sciprint("Error: Symphony environment is not initialized.\n");
	else { //There is an environment opened
		int status1=sym_get_num_elements(global_sym_env,&nz_ele); //No. of non-zero elements
		int status2=sym_get_num_cols(global_sym_env , &columns); //Columns
		int status3=sym_get_num_rows(global_sym_env , &rows); //Rows
		int status4=FUNCTION_TERMINATED_ABNORMALLY;
		
		//Make sure functions terminated normally
		if (status1 == status2 && status1 == status3 && status1 == FUNCTION_TERMINATED_NORMALLY){ 
			//Allocate memory for column-major representation			
			column_start=(int*)malloc(sizeof(int) * (columns+1));
			row_indices=(int*)malloc(sizeof(int) * nz_ele);
			elements=(double*)malloc(sizeof(double) * nz_ele);
			
			//Take column-major representation from symphony
			status4=sym_get_matrix(global_sym_env,&nz_ele,column_start,row_indices,elements); 
			if (status1 == status4) { //Check termination status of function, if normal
				
				//Allocate memory for row-major representation 
				new_list=(double*) calloc( nz_ele , sizeof(double));
				count_per_row=(int*) calloc(  rows, sizeof(int ) );
				column_position=(int*) calloc( nz_ele, sizeof(int));
				
				//Convert column-major representation to row-major representation
				column_major_to_row_major(rows,columns,nz_ele,elements,row_indices,column_start,new_list,count_per_row,column_position);
				
				/*
				(Important)Scilab considers indices from 1 , But we have column indices starting from 0 in column_position.
				 Hence add 1 to each index
				*/				
				int iter=0;
				for (;iter < nz_ele ; ++iter) column_position[iter]++;

				}
			else { //If termination status is abnormal
				sciprint("\nFunction invoked unsuccessfully.\n");
				sciprint("\n Is a problem loaded ? \n");
				}				
			}
		else //If termination status of any of functions is abnormal
			sciprint("\nFunction invoked unsuccessfully.\n");
			
		}
	
	//Copy the result to scilab. Location is position next to input arguments.
	SciErr err=createSparseMatrix(pvApiCtx,nbInputArgument(pvApiCtx)+1,rows,columns,nz_ele,count_per_row,column_position,new_list);

	/*
	 *Free allocated memory before exit
	*/
	free(row_indices);
	free(column_start);
	free(elements);
	free(new_list);
	free(count_per_row);
	free(column_position);

	if (err.iErr){ //Process error
        printError(&err, 0);
		AssignOutputVariable(pvApiCtx, 1) = 0;
        return 1;
    	}

	//assign result position to output argument
	AssignOutputVariable(pvApiCtx, 1) = nbInputArgument(pvApiCtx) + 1;
	//ReturnArguments(pvApiCtx);
	return 0;
	}
}

/*
 * It converts column-major representation to row-major representation
 * :: ARGUMENTS ::
 * rows - No. of rows IN
 * columns - No. of columns IN
 * nz_ele - No. of non-zero elements IN 
 * elements - Non-zero elements in column-major representation IN
 * row_indices - Row index( starts from 0 : symphony) of each non-zero element	IN
 * column_start - Starting index in elements of each column 	IN
 * new_list - Non-zero elements in row-major representation	 OUT 
 * count_per_row - Count of non-zero elements in each row	OUT
 * column_position - Column index ( starts from 0 (we'll add 1 to each index later)) of each non-zero element	OUT
*/
void column_major_to_row_major(int rows,int columns,int nz_ele,double *elements,int *row_indices,int *column_start,double *new_list,int *count_per_row,int *column_position) {

	int iter=0,iter2,iter3=0,index=0;
	for (iter=0;iter < rows;++iter) {
		for (iter2=0;iter2 < nz_ele;++iter2) {
			if (row_indices[iter2] == iter) {
				count_per_row[iter]++; //Count of non-zero elements per row.
				new_list[index]=elements[iter2];
				for (iter3=0; iter3 < columns+1 ; ++iter3) {
					if (iter2 < column_start[iter3])
						break;
					}
				column_position[index] = iter3 - 1;
				index++ ;
				}
			}
		}
	}