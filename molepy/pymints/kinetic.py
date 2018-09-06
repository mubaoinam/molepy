# Copyright 2018 Oinam Romesh Meitei. All Rights Reserved.
# 
# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
# 
#   http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
import numpy
from math import sqrt,pi,exp,pow

# kinetic ints with OS RR


def os_kinetic(bas1,bas2):
    
    TIJ = 0.0
    
    A = bas1.center
    B = bas2.center
    AB2 = pow(A[0]-B[0],2)+pow(A[1]-B[1],2)+pow(A[2]-B[2],2)
    
    l1,m1,n1 = bas1.cartam
    l2,m2,n2 = bas2.cartam
    
    for p1 in range(bas1.nbas):
        a1 = bas1.exps[p1]
        c1 = bas1.coefs[p1]
        
        for p2 in range(bas2.nbas):
            a2 = bas2.exps[p2]
            c2 = bas2.coefs[p2]
            
            gamma = a1 + a2
            gamma_half = 1.0/gamma

            P = numpy.zeros(3,'d')
            PA = numpy.zeros(3,'d')
            PB = numpy.zeros(3,'d')

            P[0] = (a1*A[0] + a2*B[0])*gamma_half
            P[1] = (a1*A[1] + a2*B[1])*gamma_half
            P[2] = (a1*A[2] + a2*B[2])*gamma_half

            PA = P - A
            PB = P - B

            pre_fac = exp(-a1*a2*AB2*gamma_half) * \
                      sqrt(pi*gamma_half)*pi*\
                      gamma_half*c1*c2*\
                      bas1.norms[p1] * bas2.norms[p2]

            x = {}
            y = {}
            z = {}

            x[0,0] = y[0,0] = z[0,0] = 1.0

            et2 = 1/(2*gamma)

            x[1,0] = PA[0]
            for i in range(1,l1+1):
                x[i+1,0] = PA[0]*x[i,0]+et2*i*x[i-1,0]  
            for j in range(l2+2):
                for i in range(l1+2):
                    x[i,j+1] = PB[0]*x[i,j]
                    if i > 0:
                        x[i,j+1] += et2*(i)*x[i-1,j]
                    if j > 0:
                        x[i,j+1] += et2*(j)*x[i,j-1]
            y[1,0] = PA[1]
            for i in range(1,m1+1):
                y[i+1,0] = PA[1]*y[i,0]+et2*i*y[i-1,0]
            for j in range(m2+2):
                for i in range(m1+2):
                    y[i,j+1] = PB[1]*y[i,j]
                    if i > 0:
                        y[i,j+1] += et2*(i)*y[i-1,j]
                    if j > 0:
                        y[i,j+1] += et2*(j)*y[i,j-1]
                
            z[1,0] = PA[2]
            for i in range(1,n1+1):
                z[i+1,0] = PA[2]*z[i,0]+et2*i*z[i-1,0]
            for j in range(n2+2):
                for i in range(n1+2):
                    z[i,j+1] = PB[2]*z[i,j]
                    if i > 0:
                        z[i,j+1] += et2*(i)*z[i-1,j]
                    if j > 0:
                        z[i,j+1] += et2*(j)*z[i,j-1]


            I1 = 0.0 if (l1==0 or l2==0) else\
                 0.5*l1*l2*x[l1-1,l2-1]*y[m1,m2]*z[n1,n2]*pre_fac
            I2 = 2.0*a1*a2*x[l1+1,l2+1]*y[m1,m2]*z[n1,n2]*pre_fac
            I3 = 0.0 if l2==0 else\
                 a1*l2*x[l1+1,l2-1]*y[m1,m2]*z[n1,n2]*pre_fac
            I4 = 0.0 if l1==0 else\
                 a2*l1*x[l1-1,l2+1]*y[m1,m2]*z[n1,n2]*pre_fac
            Ix = I1 + I2 -I3 - I4

            I1 = 0.0 if (m1==0 or m2==0) else\
                 0.5*m1*m2*x[l1,l2]*y[m1-1,m2-1]*z[n1,n2]*pre_fac
            I2 = 2.0*a1*a2*x[l1,l2]*y[m1+1,m2+1]*z[n1,n2]*pre_fac
            I3 = 0.0 if m2==0 else\
                 a1*m2*x[l1,l2]*y[m1+1,m2-1]*z[n1,n2]*pre_fac
            I4 = 0.0 if m1==0 else\
                 a2*m1*x[l1,l2]*y[m1-1,m2+1]*z[n1,n2]*pre_fac
            Iy = I1 + I2 -I3 - I4

            I1 = 0.0 if (n1==0 or n2==0) else\
                 0.5*n1*n2*x[l1,l2]*y[m1,m2]*z[n1-1,n2-1]*pre_fac
            I2 = 2.0*a1*a2*x[l1,l2]*y[m1,m2]*z[n1+1,n2+1]*pre_fac
            I3 = 0.0 if n2==0 else\
                 a1*n2*x[l1,l2]*y[m1,m2]*z[n1+1,n2-1]*pre_fac
            I4 = 0.0 if n1==0 else\
                 a2*n1*x[l1,l2]*y[m1,m2]*z[n1-1,n2+1]*pre_fac
            Iz = I1 + I2 -I3 - I4
            
            TIJ += Ix+Iy+Iz
    return bas1.norm*bas2.norm*TIJ
    
