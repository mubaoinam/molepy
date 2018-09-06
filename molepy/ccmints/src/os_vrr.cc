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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "boys.h"
#include "os_vrr.h"

static int tsize = MAXAM*MAXAM*6*MAXAM;
static std::vector<double> tir_x(tsize);
static std::vector<double> tir_y(tsize);
static std::vector<double> tir_z(tsize);

int iindex(int a1, int a2, int a3){
  int ival=0;
  //ival = a1*MAXAM*(6*MAXAM)+a2*(6*MAXAM)+a3;
  ival = a3+a2*6*MAXAM+a1*MAXAM*6*MAXAM;
  return ival;
}

int func1(int x1){
  return (x1<0) ? 0: x1;
}

double pVRR(double A[3], double norma, int la, 
	    int ma, int na, double a1,
	    double B[3], double normb, double a2,
	    double C[3], double normc, int lc,
	    int mc, int nc, double a3,
	    double D[3], double normd, double a4){
  double zeta,eta,AB2,CD2,PQ2,Kab,Kcd,Theta,pre_fac;
  int Ntot,i,m,j;
  zeta = a1+a2;
  eta = a3+a4;
  double P[3], Q[3], W[3];


  P[0] = (a1*A[0] + a2*B[0])/zeta;
  P[1] = (a1*A[1] + a2*B[1])/zeta;
  P[2] = (a1*A[2] + a2*B[2])/zeta;

  Q[0] = (a3*C[0] + a4*D[0])/eta;
  Q[1] = (a3*C[1] + a4*D[1])/eta;
  Q[2] = (a3*C[2] + a4*D[2])/eta;
    
  W[0] = (zeta*P[0] + eta*Q[0])/(zeta+eta);
  W[1] = (zeta*P[1] + eta*Q[1])/(zeta+eta);
  W[2] = (zeta*P[2] + eta*Q[2])/(zeta+eta);

  AB2 = pow(A[0]-B[0],2) + pow(A[1]-B[1],2) + pow(A[2]-B[2],2);
  CD2 = pow(C[0]-D[0],2) + pow(C[1]-D[1],2) + pow(C[2]-D[2],2);
  PQ2 = pow(P[0]-Q[0],2) + pow(P[1]-Q[1],2) + pow(P[2]-Q[2],2);

  Kab = exp(-a1*a2/zeta*AB2);
  Kcd = exp(-a3*a4/eta *CD2);

  Theta = zeta*eta/(zeta+eta);

  Ntot = la+ma+na+lc+mc+nc;

  Fn[Ntot] = Fm(Ntot,Theta*PQ2);
  
  for (i=Ntot-1; i>=0; i--){
    Fn[i] = (2.0*Theta*PQ2*Fn[i+1]+exp(-Theta*PQ2))/(2.0*i+1);
  }
  pre_fac = 2.0*pow(M_PI,2.5)*Kab*Kcd/(zeta*eta*sqrt(zeta+eta))*norma*normb*normc*normd;
  
  double PA[3],QC[3],Xp[3],Xq[3];

  PA[0] = P[0] - A[0];
  PA[1] = P[1] - A[1];
  PA[2] = P[2] - A[2];
  QC[0] = Q[0] - C[0];
  QC[1] = Q[1] - C[1];
  QC[2] = Q[2] - C[2];
  Xp[0] = W[0] - P[0];
  Xp[1] = W[1] - P[1];
  Xp[2] = W[2] - P[2];
  Xq[0] = W[0] - Q[0];
  Xq[1] = W[1] - Q[1];
  Xq[2] = W[2] - Q[2];

  for(m=0;m<=Ntot;m++){
    tir_x[iindex(0,0,m)] = pre_fac*Fn[m];
  }
  if (la>0||lc>0){
    for (i=0;i<=la;i++){
      if (i>0){
	for (m=0;m<Ntot-(i-1);m++){
	  tir_x[iindex(i,0,m)] = PA[0]*tir_x[iindex(i-1,0,m)]+
	    Xp[0]*tir_x[iindex(i-1,0,m+1)];
	  if (i-1>0){
	    tir_x[iindex(i,0,m)] += (i-1)/2./zeta*   
	      (tir_x[iindex(i-2,0,m)]-
	       eta/(zeta+eta)*tir_x[iindex(i-2,0,m+1)]);
	  }
	}
      }
      for (j=0;j<func1(i+lc-la);j++){
	for (m=0;m<Ntot-j-i;m++){
	  tir_x[iindex(i,j+1,m)] = QC[0]*tir_x[iindex(i,j,m)] + 
	    Xq[0]*tir_x[iindex(i,j,m+1)];
	  if (j>0){
	    tir_x[iindex(i,j+1,m)] += j/2./eta*
	      (tir_x[iindex(i,j-1,m)]- 
	       zeta/(zeta+eta)*tir_x[iindex(i,j-1,m+1)]);
	  }
	  if (i>0){
	    tir_x[iindex(i,j+1,m)] += i/2./(zeta+eta)*tir_x[iindex(i-1,j,m+1)];
	  }
	}
      }
    }   
  }
  Ntot -= la-lc;
  for(m=0;m<=Ntot;m++){
    tir_y[iindex(0,0,m)] = tir_x[iindex(la,lc,m)];
  }



  if (ma>0|| mc>0){
    for (i=0;i<=ma;i++){
      if (i>0){
	for(m=0; m<Ntot-(i-1);m++){
	  tir_y[iindex(i,0,m)] = PA[1]*tir_y[iindex(i-1,0,m)] +	
	    Xp[1]*tir_y[iindex(i-1,0,m+1)];
	  if (i-1>0){
	    tir_y[iindex(i,0,m)] += (i-1)/2./zeta*
	      (tir_y[iindex(i-2,0,m)]-
	       eta/(zeta+eta)*tir_y[iindex(i-2,0,m+1)]);
	  }
	}
      }
      for(j=0; j<func1(i+mc-ma);j++){
	for (m=0;m<Ntot-j-i;m++){
	  tir_y[iindex(i,j+1,m)] = QC[1]*tir_y[iindex(i,j,m)] +
	    Xq[1]*tir_y[iindex(i,j,m+1)];
	  if (j>0){
	    tir_y[iindex(i,j+1,m)] += j/2./eta*
	      (tir_y[iindex(i,j-1,m)]-
	       zeta/(zeta+eta)*tir_y[iindex(i,j-1,m+1)]);
	  }
	  if (i>0){
	    tir_y[iindex(i,j+1,m)] += i/2./(zeta+eta)*tir_y[iindex(i-1,j,m+1)];
	  }
	}
      }
    }
  }
  Ntot -= ma-mc;
  for (m=0;m<=Ntot;m++){
    tir_z[iindex(0,0,m)] = tir_y[iindex(ma,mc,m)];
  }
  if (na>0|| nc>0){
    for (i=0;i<=na;i++){
      if (i>0){
	for( m=0;m<Ntot-(i-1);m++){
	  tir_z[iindex(i,0,m)] = PA[2]*tir_z[iindex(i-1,0,m)] +	
	    Xp[2]*tir_z[iindex(i-1,0,m+1)];
	  if (i-1>0){
	    tir_z[iindex(i,0,m)] += (i-1)/2./zeta*	
	      (tir_z[iindex(i-2,0,m)]-
	       eta/(zeta+eta)*tir_z[iindex(i-2,0,m+1)]);
	  }
	}
      }
      for (j=0;j<func1(i+nc-na);j++){
	for (m=0;m <Ntot-j-i;m++){
	  tir_z[iindex(i,j+1,m)] = QC[2]*tir_z[iindex(i,j,m)] +	
	    Xq[2]*tir_z[iindex(i,j,m+1)];
	  if (j>0){
	    tir_z[iindex(i,j+1,m)] += j/2./eta*
	      (tir_z[iindex(i,j-1,m)]-
	       zeta/(zeta+eta)*tir_z[iindex(i,j-1,m+1)]);
	  }
	  if( i>0){
	    tir_z[iindex(i,j+1,m)] += i/2./(zeta+eta)*tir_z[iindex(i-1,j,m+1)];
	  }
	}
      }
    }
  }
  double result = tir_z[iindex(na,nc,0)];
  return result;
}


double cVRR(double A[3], std::vector<double> &norma, std::vector<double> &aexps,
	    std::vector<double> &acoefs, int la, int ma, int na,
	    double B[3], std::vector<double> &normb, std::vector<double> &bexps,
	    std::vector<double> &bcoefs,
	    double C[3], std::vector<double> &normc, std::vector<double> &cexps,
	    std::vector<double> &ccoefs, int lc, int mc, int nc,
	    double D[3], std::vector<double> &normd, std::vector<double> &dexps,
	    std::vector<double> &dcoefs){
  int lena, lenb, lenc, lend;
  int i,j,k,l;


  lena = aexps.size();
  lenb = bexps.size();
  lenc = cexps.size();
  lend = dexps.size();
  
  double val = 0.;
  for(i=0;i<lena;i++){
    for(j=0;j<lenb;j++){
      for(k=0;k<lenc;k++){
	for(l=0;l<lend;l++){
	  val += acoefs[i]*bcoefs[j]*ccoefs[k]*dcoefs[l]*
	    pVRR(A,norma[i],la,ma,na,aexps[i],
		 B,normb[j],bexps[j],
		 C,normc[k],lc,mc,nc,cexps[k],
		 D,normd[l],dexps[l]);
	}
      }
    }
  }  
  return val;
}
