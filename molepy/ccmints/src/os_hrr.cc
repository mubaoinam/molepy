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
#include "os_vrr.h"

namespace py = pybind11;

static int ksize = (1-pow(MAXAM,7))/(1-MAXAM);

static std::vector<double> hrr_(ksize);
static std::vector<double> tmp_(ksize);


int kindex(int a1, int a2, int a3, int a4, int a5, int a6){
   int jval = 0;
   jval = a1*pow(MAXAM,5)+a2*pow(MAXAM,4)+a3*pow(MAXAM,3)+a4*pow(MAXAM,2)+a5*pow(MAXAM,1)+a6;
   //jval = a6+a5*pow(MAXAM,1)+a4*pow(MAXAM,2)+a3*pow(MAXAM,3)+a2*pow(MAXAM,4)+a1*pow(MAXAM,5);
   return jval;
}

double twofunc_(double A[3], std::vector<double> &norma, int I, int J, int K,
		std::vector<double> &aexps, std::vector<double> &acoefs,
		double B[3],std::vector<double> &normb,
		std::vector<double> &bexps, std::vector<double> &bcoefs,
		double C[3], std::vector<double> &normc, int lc, int mc, int nc,
		std::vector<double> &cexps, std::vector<double> &ccoefs,
		double D[3], std::vector<double> &normd, int ld, int md, int nd, 
		std::vector<double> &dexps, std::vector<double> &dcoefs){  
  int i,j,k,l;

  for(i=0;i<nd+2;i++){
    for(j=lc+ld;j>lc-1;j--){
      for(k=mc+md;k>mc-1;k--){
	for(l=nc+nd-i;l>nc-1;l--){
	  if(i==0){
	    tmp_[kindex(0,0,0,j,k,l)] = cVRR(A, norma, aexps, acoefs, I, J, K,
					     B, normb, bexps, bcoefs,
					     C, normc, cexps, ccoefs, j, k, l,
					     D, normd, dexps, dcoefs);
	    
	  } else {
	    tmp_[kindex(i,0,0,j,k,l)] = tmp_[kindex(i-1,0,0,j,k,l+1)]+(C[2]-D[2])*tmp_[kindex(i-1,0,0,j,k,l)];
	  }
           
	}
      }
    }
  }
  
  for(i=1;i<md+1;i++){
    for(j=lc+ld;j>lc-1;j--){
      for(k=mc+md-i;k>mc-1;k--){
	tmp_[kindex(nd,i,0,j,k,nc)]=tmp_[kindex(nd,i-1,0,j,k+1,nc)]+(C[1]-D[1])*tmp_[kindex(nd,i-1,0,j,k,nc)];
      }
    }
  }
  for(i=1;i<ld+1;i++){
    for(j=lc+ld-i;j>lc-1;j--){
      tmp_[kindex(nd,md,i,j,mc,nc)] = tmp_[kindex(nd,md,i-1,j+1,mc,nc)]+(C[0]-D[0])*tmp_[kindex(nd,md,i-1,j,mc,nc)];
    }
  }
  double result = tmp_[kindex(nd,md,ld,lc,mc,nc)];
  return result;
}
	      

double ccHRR(double xa, double ya, double za ,
	     std::vector<double> &norma, int la, int ma, int na,
	     std::vector<double> &aexps, std::vector<double> &acoefs,
	     double xb, double yb, double zb,
	     std::vector<double> &normb, int lb, int mb, int nb, 
	     std::vector<double> &bexps, std::vector<double> &bcoefs,
	     double xc, double yc, double zc,
	     std::vector<double> &normc, int lc, int mc, int nc,
	     std::vector<double> &cexps, std::vector<double> &ccoefs,
	     double xd, double yd, double zd,
	     std::vector<double> &normd, int ld, int md, int nd, 
	     std::vector<double> &dexps, std::vector<double> &dcoefs){
  
  int i,j,k,l;
  double A[3] = {xa,ya,za};
  double B[3] = {xb,yb,zb};
  double C[3] = {xc,yc,zc};
  double D[3] = {xd,yd,zd};  
  
  for(i=la;i<la+lb+1;i++){
    for(j=ma;j<ma+mb+1;j++){
      for(k=na;k<na+nb+1;k++){
	
	hrr_[kindex(i,j,k,0,0,0)]=twofunc_(A,norma,i,j,k,aexps,acoefs,
					   B,normb,bexps,bcoefs,
					   C,normc,lc,mc,nc,cexps,ccoefs,
					   D,normd,ld,md,nd,dexps,dcoefs);
	
      }
    }
  }
  
  for(i=1;i<nb+1;i++){
    for(j=la;j<la+lb+1;j++){
      for(k=ma;k<ma+mb+1;k++){
	for(l=na;l<na+nb+1-i;l++){
	  hrr_[kindex(j,k,l,0,0,i)]=hrr_[kindex(j,k,l+1,0,0,i-1)]+(A[2]-B[2])*hrr_[kindex(j,k,l,0,0,i-1)];
	}
      }
    }
  }
  for(i=1;i<mb+1;i++){
    for(j=la;j<la+lb+1;j++){
      for(k=ma;k<ma+mb+1-i;k++){
	hrr_[kindex(j,k,na,0,i,nb)]=hrr_[kindex(j,k+1,na,0,i-1,nb)]+(A[1]-B[1])*hrr_[kindex(j,k,na,0,i-1,nb)];
      }
    }
  }
  for(i=1;i<lb+1;i++){
    for(j=la;j<la+lb+1-i;j++){
      hrr_[kindex(j,ma,na,i,mb,nb)]=hrr_[kindex(j+1,ma,na,i-1,mb,nb)]+(A[0]-B[0])*hrr_[kindex(j,ma,na,i-1,mb,nb)];
    }
  }
  double result= hrr_[kindex(la,ma,na,lb,mb,nb)];
  return result;
}
	  
	   
PYBIND11_MODULE(os_hrr,m) {
  m.def("ccHRR", &ccHRR, "Contracted HRR");
}
