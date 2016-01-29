// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Harpreet Singh
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

demopath = get_absolute_file_path("sci_FOSSEE_Optimization_Toolbox.dem.gateway.sce");

subdemolist = ["Linprog","linprog.dem.sce";"Symphony", "symphony.dem.sce"; "SymphonyMat", "symphonymat.dem.sce"; "Qpipopt", "qpipopt.dem.sce"; "QpipoptMat", "qpipoptmat.dem.sce";"Lsqlin","lsqlin.dem.sce";"Lsqnonneg","lsqnonneg.dem.sce";"Fminunc","fminunc.dem.sce";"Fminbnd","fminbnd.dem.sce";"Fmincon","fmincon.dem.sce"];

subdemolist(:,2) = demopath + subdemolist(:,2);
