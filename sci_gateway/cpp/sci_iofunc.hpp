// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Keyur Joshi, Harpreet Singh
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#ifndef SCI_IOFUNCHEADER
#define SCI_IOFUNCHEADER

//input
int getDoubleFromScilab(int argNum, double *dest);
int getUIntFromScilab(int argNum, int *dest);
int getIntFromScilab(int argNum, int *dest);
int getFixedSizeDoubleMatrixFromScilab(int argNum, int rows, int cols, double **dest);
int getDoubleMatrixFromScilab(int argNum, int *rows, int *cols, double **dest);
int getFixedSizeDoubleMatrixInList(int argNum, int itemPos, int rows, int cols, double **dest);


//output
int return0toScilab();
int returnDoubleToScilab(double retVal);
int returnDoubleMatrixToScilab(int itemPos, int rows, int cols, double *dest);
int returnIntegerMatrixToScilab(int itemPos, int rows, int cols, int *dest);

#endif //SCI_IOFUNCHEADER
