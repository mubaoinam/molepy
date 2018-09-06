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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <python2.7/Python.h>
#include <math.h>
#define MAXAM 7

namespace py = pybind11;

static double os_overlap(double x1, double y1, double z1,
			 std::vector<double> &norm1, std::vector<double> &exps1,
			 std::vector<double> &coef1, int l1, int m1, int n1,
			 double x2, double y2, double z2,
			 std::vector<double> &norm2, std::vector<double> &exps2,
			 std::vector<double> &coef2, int l2, int m2, int n2){

  double SIJ = 0.0;
  double AB2 = pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2);

  int p1,p2,i,j;
  
  const int nbas1 = coef1.size();
  const int nbas2 = coef2.size();

  static std::vector<double> x(MAXAM*MAXAM);
  static std::vector<double> y(MAXAM*MAXAM);
  static std::vector<double> z(MAXAM*MAXAM);

  double a1,a2,c1,c2,gamma,gamma_half,pre_fac,et2;
  double P[3],PA[3],PB[3];
  
  for (p1=0; p1<nbas1; p1++){
    a1 = exps1[p1];
    c1 = coef1[p1];

    for(p2=0;p2<nbas2; p2++){
      a2 = exps2[p2];
      c2 = coef2[p2];
      
      gamma = a1+a2;
      gamma_half = 1.0/gamma;

      P[0] = (a1*x1+a2*x2)*gamma_half;
      P[1] = (a1*y1+a2*y2)*gamma_half;
      P[2] = (a1*z1+a2*z2)*gamma_half;

      PA[0] = P[0] - x1;
      PA[1] = P[1] - y1;
      PA[2] = P[2] - z1;
      PB[0] = P[0] - x2;
      PB[1] = P[1] - y2;
      PB[2] = P[2] - z2;

      pre_fac = exp(-a1*a2*AB2*gamma_half)*
	sqrt(M_PI*gamma_half)*M_PI*gamma_half*
	c1*c2*norm1[p1]*norm2[p2];
      x[0] = 1.0;
      y[0] = 1.0;
      z[0] = 1.0;
      
      et2 = 1/(2*gamma);
      x[1*MAXAM+0] = PA[0];
      for(i=1;i<l1;i++){
	x[(i+1)*MAXAM+0] = PA[0]*x[i*MAXAM+0]+et2*i*x[(i-1)*MAXAM+0];
      }
      for(j=0;j<=l2;j++){
	for(i=0;i<=l1;i++){
	  //std::cout<<i*l1+j+1<<std::endl;
	  x[i*MAXAM+j+1] = PB[0]*x[i*MAXAM+j];
	  if (i>0){
	    x[i*MAXAM+j+1] += et2*(i)*x[(i-1)*MAXAM+j];
	  }
	  if (j>0){
	    x[i*MAXAM+j+1] += et2*(j)*x[i*MAXAM+j-1];
	  }
	}
      }
      y[1*MAXAM+0] = PA[1];
      for (i=1;i<m1;i++){
	y[(i+1)*MAXAM+0] = PA[1]*y[i*MAXAM+0]+et2*i*y[(i-1)*MAXAM+0];
      }
      for (j=0;j<=m2;j++){
	for (i=0;i<=m1;i++){
	  y[i*MAXAM+j+1] = PB[1]*y[i*MAXAM+j];
	  if (i > 0){
	    y[i*MAXAM+j+1] += et2*(i)*y[(i-1)*MAXAM+j];
	  }
	  if( j > 0){
	    y[i*MAXAM+j+1] += et2*(j)*y[i*MAXAM+j-1];
	  }
	}
      }                
      z[1*MAXAM+0] = PA[2];
      for (i=1;i<n1;i++){
	z[(i+1)*MAXAM+0] = PA[2]*z[i*MAXAM+0]+et2*i*z[(i-1)*MAXAM+0];
      }
      for(j=0;j<=n2; j++){
	for (i=0;i<=n1;i++){
	  z[i*MAXAM+j+1] = PB[2]*z[i*MAXAM+j];
	  if (i > 0){
	    z[i*MAXAM+j+1] += et2*(i)*z[(i-1)*MAXAM+j];
	  }
	  if( j > 0){
	    z[i*MAXAM+j+1] += et2*(j)*z[i*MAXAM+j-1];
	  }
	}
      }
      SIJ += pre_fac*x[l1*MAXAM+l2]*y[m1*MAXAM+m2]*z[n1*MAXAM+n2];
    }
  }
  return SIJ;
}


	   
PYBIND11_MODULE(overlap,m) {
  m.def("os_overlap", &os_overlap, "Overlap");
}
