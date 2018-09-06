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
#include <assert.h>
#include <iostream>
#include "boys.h"
#define SMALL 0.00000001
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30


static const double lgam_[12] = {0.572364942924700,-0.120782237635245,
				 0.284682870472919,1.200973602347074,
				 2.453736570842442,3.957813967618717,
				 5.662562059857141,7.534364236758734,
				 9.549267257300997,11.689333420797269,
				 13.940625219403762,16.292000476567239};

static const double lgam_e[12] = {1.772453850905517,
				  0.886226925452759,
				  1.329340388179137,
				  3.323350970447843,
				  11.631728396567436,
				  52.342777784553483,
				  287.885277815043651,
				  1871.254305797787765,
				  14034.407293483451213,
				  119292.461994609009707,
				  1133278.388948787236586,
				  11899423.0839622188359500};

double gammln(double xx){
  double x,y,ser,tmp;
  int j;
  static double cof[14] ={57.1562356658629235,-59.5979603554754912,
			  14.1360979747417471,-0.491913816097620199,
			  .339946499848118887e-4,.465236289270485756e-4,
			  -.983744753048795646e-4,.158088703224912494e-3,
			  -.210264441724104883e-3,.217439618115212643e-3,
			  -.164318106536763890e-3,.844182239838527433e-4,
			  -.261908384015814087e-4,.368991826595316234e-5};
  
  y = x = xx;
  tmp = x+5.24218750000000000;
  tmp = (x+0.5)*log(tmp)-tmp;
  ser = 0.999999999999997092;
  for (j=0;j<14;j++){
    y += 1;
    ser += cof[j]/y;
  }
  return tmp+log(2.5066282746310005*ser/x);
}


double gser( int a, double x){
  
  //double gln = gammln(a);
  double sum,ap,del_;
  int n;
  ap = a+0.5;
  del_ = sum = 1.0/(a+0.5);
  for (n=1;n<=100;n++){
    ++ap;
    del_ *= x/ap;
    sum += del_;
    if (fabs(del_) < fabs(sum)*EPS){
      return sum*exp(-x+(a+0.5)*log(x)-lgam_[a]);
    }
  }
  return sum*exp(-x+(a+0.5)*log(x)-lgam_[a]);
}

double gcf(int a, double x){
  
  //double gln = gammln(a);
  double b,c,d,h,an,del_;
  int i;

  b=x+1.0-(a+0.5);
  c = 1.0/FPMIN;
  d = 1.0/b;
  h=d;
  for (i=1;i<100;i++){
    an = -i*(i-(a+0.5));
    b+=2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del_=d*c;
    h *= del_;
    if (fabs(del_-1.0) < EPS) return exp(-x+(a+0.5)*log(x)-lgam_[a])*h;
  }
  return exp(-x+(a+0.5)*log(x)-lgam_[a])*h;
}

double gammp(int a, double x){
  assert(x>0.);
  assert((a+0.5)>=0.);
  
  if (x == 0.0){
    return 0.0;
  }
  if (x<(a+0.5)+1.0){
    return gser(a,x);
  }
  else{
    return 1.0-gcf(a,x);
  }
}
      
double Fm(int m, double x){
  //double incomplete;
  if (fabs(x) < 0.00000001) x = 0.00000001;
  //incomplete = gammp(m,x);
  //gamma = exp(gammln(m+0.5));
  return lgam_e[m]*gammp(m,x)*0.5*pow(x,-m-0.5);
}  

