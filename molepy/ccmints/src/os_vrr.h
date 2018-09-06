/*Copyright 2018 Oinam Romesh Meitei. All Rights Reserved.

Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/
#ifndef OS_VRR_H
#define OS_VRR_H
#define MAXAM 7


static std::vector<double> Fn(6*MAXAM+1);

double pVRR(double A[3], double norma, int la, 
	    int ma, int na, double a1,
	    double B[3], double normb, double a2,
	    double C[3], double normc, int lc,
	    int mc, int nc, double a3,
	    double D[3], double normd, double a4);

double cVRR(double A[3], std::vector<double> &norma, std::vector<double> &aexps,
	    std::vector<double> &acoefs, int la, int ma, int na,
	    double B[3], std::vector<double> &normb, std::vector<double> &bexps,
	    std::vector<double> &bcoefs,
	    double C[3], std::vector<double> &normc, std::vector<double> &cexps,
	    std::vector<double> &ccoefs, int lc, int mc, int nc,
	    double D[3], std::vector<double> &normd, std::vector<double> &dexps,
	    std::vector<double> &dcoefs);
#endif
