#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <python2.7/Python.h>
#include <pybind11/numpy.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

namespace py = pybind11;

struct LCart{
  int lx, ly, lz;
  LCart(int lx, int ly, int lz)
    : lx(lx), ly(ly), lz(lz){}
};

static std::vector<std::vector<LCart>> lcart = {
  {{0,0,0}},{{1,0,0},{0,1,0},{0,0,1}},
  {{2,0,0},{1,1,0},{1,0,1},{0,2,0},{0,1,1},{0,0,2}},
  {{3,0,0},{2,1,0},{2,0,1},{1,2,0},{1,1,1},{1,0,2},
   {0,3,0},{0,2,1},{0,1,2}, {0,0,3}}};

static std::vector<std::vector<int>> mquant = {
  {0},{-1,0,1},{-2,-1,0,1,2},{-3,-2,-1,0,1,2,3}};


// Coefficients

static const double coeff_0[1] = { 1.000000000000000};
static const double coeff_1[3*3] = {0.000000000000000,1.000000000000000,0.000000000000000,
				    0.000000000000000,0.000000000000000,1.000000000000000,
				    1.000000000000000,0.000000000000000,0.000000000000000};
static const double coeff_2[5*6] = {0.000000000000000,1.000000000000000,0.000000000000000,
				    0.000000000000000,0.000000000000000,0.000000000000000,
				    0.000000000000000,0.000000000000000,0.000000000000000,
				    0.000000000000000,1.000000000000000,0.000000000000000,
				    -0.500000000000000,0.000000000000000,0.000000000000000,
				    -0.500000000000000,0.000000000000000,1.000000000000000,
				    0.000000000000000,0.000000000000000,1.000000000000000,
				    0.000000000000000,0.000000000000000,0.000000000000000,
				    0.866025403784439,0.000000000000000,0.000000000000000,
				    -0.866025403784439,0.000000000000000,0.000000000000000};
static const double coeff_3[7*10] = {0.000000000000000,1.060660171779821,0.000000000000000,
				     0.000000000000000,0.000000000000000,0.000000000000000,
				     -0.790569415042095,0.000000000000000,0.000000000000000,
				     0.000000000000000,
				     0.000000000000000,0.000000000000000,0.000000000000000,
				     0.000000000000000,1.000000000000000,0.000000000000000,
				     0.000000000000000,0.000000000000000,0.000000000000000,
				     0.000000000000000,
				     0.000000000000000,-0.273861278752583,0.000000000000000,
				     0.000000000000000,0.000000000000000,0.000000000000000,
				     -0.612372435695795,0.000000000000000,1.095445115010332,
				     0.000000000000000,
				     0.000000000000000,0.000000000000000,-0.670820393249937,
				     0.000000000000000,0.000000000000000,0.000000000000000,
				     0.000000000000000,-0.670820393249937,0.000000000000000,
				     1.000000000000000,
				     -0.612372435695795,0.000000000000000,0.000000000000000,
				     -0.273861278752583,0.000000000000000,1.095445115010332,
				     0.000000000000000,0.000000000000000,0.000000000000000,
				     0.000000000000000,
				     0.000000000000000,0.000000000000000,0.866025403784439,
				     0.000000000000000,0.000000000000000,0.000000000000000,
				     0.000000000000000,-0.866025403784439,0.000000000000000,
				     0.000000000000000,
				     0.790569415042095,0.000000000000000,0.000000000000000,
				     -1.060660171779821,0.000000000000000,0.000000000000000,
				     0.000000000000000,0.000000000000000,0.000000000000000,
				     0.000000000000000};

static const double* CScoeff[] = { coeff_0, coeff_1, coeff_2, coeff_3};
				     

static int Index(int i, int j, int k, int l){
  int tmp,ij,kl;
  if(i<j){
    tmp = i;
    i = j;
    j = tmp;
  }
  if (k<l){
    tmp = k;
    k = l;
    l = tmp;
  }
  ij = i*(i+1)/2+j;
  kl = k*(k+1)/2+l;
  if (ij < kl){
    tmp = ij;
    ij = kl;
    kl = tmp;
  }
  return ij*(ij+1)/2+kl;
}

int factorial(int n){
  return (n>1) ? n*factorial(n-1) : 1;
}

int dfact(int n){
  return (n<=1) ? 1:n*dfact(n-2);
}

int binomial(int x,int y){
  return factorial(x)/factorial(y)/factorial(x-y);
}
 

double coeff_compute(int l,int m,int lx,int ly,int lz){
  int u;  
  int t = l-lz-fabs(m);
  int p = (m>=0) ? (lx-m) : (lx+m-1);
  
  if ((t%2 != 0) || (p%2 != 0)){
    return 0.0E+0;
  }

  t /= 2;
  p /= 2;

  double pfac = double(sqrt(factorial(l+fabs(m))*factorial(l-fabs(m))));
  pfac *= (m==0) ? 1.000000E+00 : 1.4142135623730950488016887E+00;
  pfac /= double(pow(2,fabs(m))*factorial(l));
  pfac *= double(binomial(l,t)*binomial(l-t,fabs(m)+t))/pow(4,t);
  double sum =0.0000000000E+00;
  for (u=0;u<=t;u++){
    double sum1 = double(binomial(t,u)*binomial(fabs(m),ly-u-u));
    if ((ly-u+p)%2 !=0){
      sum -= sum1;
    } else {
      sum += sum1;
    }
  }
  sum *= pfac;
  sum *= sqrt(double(dfact(2*lx-1)*dfact(2*ly-1)*dfact(2*lz-1))/
	      dfact(2*l-1));
  
  return sum;
}


std::vector<double> Cart2Spher2(std::vector<int> shells, std::vector<double> &CartInt, int ncart){

  int nsph = 0;
  for(auto const& i : shells){
    nsph += 2*i+1;
  }
  std::vector<double> SphInt(nsph*nsph);

  int c1ind,c2ind;
  double coeff1,coeff2;
  int sph2,s2ind;
  int sph1 = 0;
  int s1ind = 0;
  for(auto const& l1: shells){
    for(auto const& m1: mquant[l1]){
      sph2 = 0;
      s2ind = 0;
      for(auto const& l2: shells){
	for(auto const& m2: mquant[l2]){
	  c1ind = sph1;
	  for(auto const& c1: lcart[l1]){
	    c2ind  = sph2;
	    coeff1 = (l1>3) ? coeff_compute(l1,m1,c1.lx,c1.ly,c1.lz):
	      CScoeff[l1][(&m1-&mquant[l1][0])*lcart[l1].size()+(&c1-&lcart[l1][0])];
	    for(auto const& c2: lcart[l2]){
	      coeff2 = (l2>3) ? coeff_compute(l2,m2,c2.lx,c2.ly,c2.lz):
		CScoeff[l2][(&m2-&mquant[l2][0])*lcart[l2].size()+(&c2-&lcart[l2][0])];	      
	      SphInt[s1ind*nsph+s2ind] += coeff1*coeff2*
		CartInt[c1ind*ncart+c2ind];
	      c2ind += 1;
	    }
	    c1ind += 1;
	  }
	  s2ind += 1;
	}
	sph2 += lcart[l2].size();
      }
      s1ind += 1;
    }
    sph1 += lcart[l1].size();
  }
  return SphInt;
}





std::vector<double> Cart2Spher4(std::vector<int> shells, std::vector<double> &CartInt){
  
  int i,c1ind,c2ind,c3ind,c4ind;
  const int shells_size = shells.size();
  
  int nsph = 0;
  for(i=0;i<shells_size;i++){
    nsph += 2*shells[i]+1;
  }
  const int tbas = nsph*(nsph+1)*(nsph*nsph+nsph+2)/8;
  static std::vector<double> SphInt(tbas);
  
  int sph2,s2ind,m1m2,sph3,s3ind,shell3,sph4,s4ind,m3m4;
  double coeff1,coeff2,coeff3,coeff4;

  int sph1 = 0;
  int s1ind = 0;
  int shell1 = 0;
  for (auto const& l1: shells){
    for (auto const& m1: mquant[l1]){   
      
      sph2 = 0;
      s2ind = 0;
      for (auto const& l2 : std::vector<int> (shells.begin(),shells.begin()+shell1+1)){
	
	for (auto const& m2: mquant[l2]){

	  if (s2ind > s1ind)
	    continue;
	  m1m2 =  s1ind*(s1ind+1)/2+s2ind;
                    
	  sph3 = 0;
	  s3ind = 0;
	  shell3 = 0;
	  for (auto const& l3: shells){
	    for ( auto const& m3: mquant[l3]){

	      sph4 = 0;
	      s4ind = 0;
	      for (auto const& l4 : std::vector<int> (shells.begin(),shells.begin()+shell3+1)){
		for (auto const& m4: mquant[l4]){
		  if (s4ind > s3ind)
		    continue;
		  m3m4 = s3ind*(s3ind+1)/2+s4ind;
		  if (m1m2 <= m3m4){
                                    
		    c1ind = sph1;
		    for (auto const& c1: lcart[l1]){
                                        
		      c2ind = sph2;
		      coeff1 = (l1 > 3) ? coeff_compute(l1,m1,c1.lx,c1.ly,c1.lz): 
			CScoeff[l1][(&m1-&mquant[l1][0])*lcart[l1].size()+(&c1-&lcart[l1][0])];
		      for (auto const& c2: lcart[l2]){
                                        
			c3ind = sph3;
			coeff2 = (l2>3) ? coeff_compute(l2,m2,c2.lx,c2.ly,c2.lz):
			  CScoeff[l2][(&m2-&mquant[l2][0])*lcart[l2].size()+(&c2-&lcart[l2][0])];
			for (auto const& c3: lcart[l3]){
                                        
			  c4ind = sph4;
			  coeff3 = (l3>3) ? coeff_compute(l3,m3,c3.lx,c3.ly,c3.lz):
			    CScoeff[l3][(&m3-&mquant[l3][0])*lcart[l3].size()+(&c3-&lcart[l3][0])];
			  for (auto const& c4 : lcart[l4]){
                                        
			    coeff4 = (l4>3) ? coeff_compute(l4,m4,c4.lx,c4.ly,c4.lz):
			      CScoeff[l4][(&m4-&mquant[l4][0])*lcart[l4].size()+(&c4-&lcart[l4][0])];
			    SphInt[Index(s1ind,s2ind,s3ind,s4ind)] += coeff1*coeff2*coeff3*coeff4* 
			      CartInt[Index(c1ind,c2ind,c3ind,c4ind)];
                                                        
			    c4ind += 1;
			  }
			  c3ind += 1;
			}
			c2ind += 1;
		      }
		      c1ind += 1;
		    }
		  }
		  s4ind += 1;
		}
		sph4 += lcart[l4].size() ;
	      }    
	      s3ind += 1 ;
	    }
	    sph3 += lcart[l3].size();
	    shell3 += 1;
	  }
	  s2ind += 1 ;
	}
	sph2 += lcart[l2].size();
      }
      s1ind += 1;
    }
    sph1 += lcart[l1].size();
    shell1 += 1;
    
  }
  return SphInt;
}

            

	   
PYBIND11_MODULE(CartSpher,m) {
  m.def("Cart2Spher4", &Cart2Spher4, "Cart to Spher");
  m.def("Cart2Spher2", &Cart2Spher2, "Cart to Spher");
}
