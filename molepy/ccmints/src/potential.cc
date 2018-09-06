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
#include <pybind11/stl.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "boys.h"
#define MAXAM 7

namespace py = pybind11;

static int tsize = MAXAM*MAXAM*6*MAXAM;
static std::vector<double> tir_x(tsize);
static std::vector<double> tir_y(tsize);
static std::vector<double> tir_z(tsize);
static std::vector<double> Fn(6*MAXAM+1);

int iindex(int a1, int a2, int a3){
  int ival=0;
  ival = a1*MAXAM*(6*MAXAM)+a2*(6*MAXAM)+a3;
  return ival;
}

static int onefunc( int x){
  return (x<0) ? 0 : x;
}


static double os_Vint(double x1, double y1, double z1,
		      std::vector<double> &norm1, std::vector<double> &exps1,
		      std::vector<double> &coef1, int l1, int m1, int n1,
		      double x2, double y2, double z2,
		      std::vector<double> &norm2, std::vector<double> &exps2,
		      std::vector<double> &coef2, int l2, int m2, int n2,
		      std::vector<double> &coords, std::vector<int> &anum){
  double VIJ = 0.0;
  double AB2 = pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2);
  

  int mtot = l1+m1+n1+l2+m2+n2+0;

  const int nbas1 = coef1.size();
  const int nbas2 = coef2.size();
  const int natom = anum.size();

  double C[3],P[3],PA[3],PB[3],PC[3];
  double a1,c1,a2,c2,gamma,Kab,pre_fac,rpc2;

  int p1,p2,m,i,j,atom;

  

  for (p1=0; p1<nbas1; p1++){
    a1 = exps1[p1];
    c1 = coef1[p1];

    for(p2=0;p2<nbas2; p2++){
      a2 = exps2[p2];
      c2 = coef2[p2];
      
      gamma = a1+a2;
      Kab = exp(-a1*a2*AB2/gamma);
      pre_fac = (2.*M_PI/gamma)*Kab*c1*c2*norm1[p1]*norm2[p2];


      P[0] = (a1*x1+a2*x2)/gamma;
      P[1] = (a1*y1+a2*y2)/gamma;
      P[2] = (a1*z1+a2*z2)/gamma;

      PA[0] = P[0] - x1;
      PA[1] = P[1] - y1;
      PA[2] = P[2] - z1;
      PB[0] = P[0] - x2;
      PB[1] = P[1] - y2;
      PB[2] = P[2] - z2;

      for (atom = 0;atom<natom;atom++){
	C[0] = coords[atom*3+0];
	C[1] = coords[atom*3+1];
	C[2] = coords[atom*3+2];
	
	PC[0] = P[0] - C[0];
	PC[1] = P[1] - C[1];
	PC[2] = P[2] - C[2];
	rpc2 = pow(PC[0],2)+pow(PC[1],2)+pow(PC[2],2);
	
	Fn[mtot] = Fm(mtot,gamma*rpc2);
	for(m=mtot-1;m>=0;m--){
	  Fn[m] = (2.*gamma*rpc2*Fn[m+1]+exp(-gamma*rpc2))/(2.*m+1);
	}

	// OS recurrence begin
	for(m=0;m<=mtot;m++){
	  tir_x[iindex(0,0,m)] = pre_fac*Fn[m];
	}
       	if (l1>0||l2>0){
	  for (i=0;i<=l1;i++){
	    if (i>0){
	      for (m=0;m<mtot-(i-1);m++){
		tir_x[iindex(i,0,m)] = PA[0]*tir_x[iindex(i-1,0,m)]-PC[0]*
		  tir_x[iindex(i-1,0,m+1)];
		if (i-1>0){
		  tir_x[iindex(i,0,m)] += (i-1)/2./gamma*   
		    (tir_x[iindex(i-2,0,m)]-
		     tir_x[iindex(i-2,0,m+1)]);
		}
	      }
	    }
	    for (j=0;j<onefunc(i+l2-l1);j++){
	      for (m=0;m<mtot-j-i;m++){
		tir_x[iindex(i,j+1,m)] = PB[0]*tir_x[iindex(i,j,m)] -PC[0]* 
		  tir_x[iindex(i,j,m+1)];
		if (j>0){
		  tir_x[iindex(i,j+1,m)] += j/2./gamma*
		    (tir_x[iindex(i,j-1,m)]- 
		     tir_x[iindex(i,j-1,m+1)]);
		}
		if (i>0){
		  tir_x[iindex(i,j+1,m)] += i/2./gamma*
		    (tir_x[iindex(i-1,j,m)]-tir_x[iindex(i-1,j,m+1)]);
		}
	      }
	    }
	  }   
	}
	for(m=0;m<=mtot-l1-l2;m++){
	  tir_y[iindex(0,0,m)] = tir_x[iindex(l1,l2,m)];
	}
	if (m1>0|| m2>0){
	  for (i=0;i<=m1;i++){
	    if (i>0){
	      for(m=0; m<mtot-(i-1)-l1-l2;m++){
		tir_y[iindex(i,0,m)] = PA[1]*tir_y[iindex(i-1,0,m)] -PC[1]*	
		  tir_y[iindex(i-1,0,m+1)];
		if (i-1>0){
		  tir_y[iindex(i,0,m)] += (i-1)/2./gamma*
		    (tir_y[iindex(i-2,0,m)]-
		     tir_y[iindex(i-2,0,m+1)]);
		}
	      }
	    }
	    for(j=0; j<onefunc(i+m2-m1);j++){
	      for (m=0;m<mtot-j-i-l1-l2;m++){
		tir_y[iindex(i,j+1,m)] = PB[1]*tir_y[iindex(i,j,m)] -PC[1]*
		  tir_y[iindex(i,j,m+1)];
		if (j>0){
		  tir_y[iindex(i,j+1,m)] += j/2./gamma*
		    (tir_y[iindex(i,j-1,m)]-
		     tir_y[iindex(i,j-1,m+1)]);
		}
		if (i>0){
		  tir_y[iindex(i,j+1,m)] += i/2./gamma*
		    (tir_y[iindex(i-1,j,m)]-tir_y[iindex(i-1,j,m+1)]);
		}
	      }
	    }
	  }
	}
	for (m=0;m<=mtot-m1-m2;m++){
	  tir_z[iindex(0,0,m)] = tir_y[iindex(m1,m2,m)];
	}
	if (n1>0|| n2>0){
	  for (i=0;i<=n1;i++){
	    if (i>0){
	      for( m=0;m<mtot-(i-1)-m1-m2;m++){
		tir_z[iindex(i,0,m)] = PA[2]*tir_z[iindex(i-1,0,m)] -PC[2]*	
		  tir_z[iindex(i-1,0,m+1)];
		if (i-1>0){
		  tir_z[iindex(i,0,m)] += (i-1)/2./gamma*	
		    (tir_z[iindex(i-2,0,m)]-
		     tir_z[iindex(i-2,0,m+1)]);
		}
	      }
	    }
	    for (j=0;j<onefunc(i+n2-n1);j++){
	      for (m=0;m <mtot-j-i-m1-m2;m++){
		tir_z[iindex(i,j+1,m)] = PB[2]*tir_z[iindex(i,j,m)] -PC[2]*	
		  tir_z[iindex(i,j,m+1)];
		if (j>0){
		  tir_z[iindex(i,j+1,m)] += j/2./gamma*
		    (tir_z[iindex(i,j-1,m)]-
		     tir_z[iindex(i,j-1,m+1)]);
		}
		if( i>0){
		  tir_z[iindex(i,j+1,m)] += i/2./gamma*
		    (tir_z[iindex(i-1,j,m)]-tir_z[iindex(i-1,j,m+1)]);
		}
	      }
	    }
	  }
	}
	VIJ += anum[atom]*-tir_z[iindex(n1,n2,0)];
      }
    }
  }
  return VIJ;
}
	
	   
PYBIND11_MODULE(potential,m) {
  m.def("os_Vint", &os_Vint, "Potentials");
}
