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
from math import exp,pi,pow,sqrt
from molepy.pymints.boys import Fm

def cVRR(A,norma,Aexps,Acoefs,amA,
         B,normb,Bexps,Bcoefs,
         C,normc,Cexps,Ccoefs,amC,
         D,normd,Dexps,Dcoefs):
    " Bra and ket contraction "
    " (m0|n0) from [m0|n0](0) "
    
    val = 0.0
    
    for i in range(len(Acoefs)):
        for j in range(len(Bcoefs)):
            for k in range(len(Ccoefs)):
                for l in range(len(Dcoefs)):
                    val += Acoefs[i]*Bcoefs[j]*\
                           Ccoefs[k]*Dcoefs[l]*\
                           pVRR(A,norma[i],amA,Aexps[i],
                                B,normb[j],Bexps[j],
                                C,normc[k],amC,Cexps[k],
                                D,normd[l],Dexps[l])
    return val
    
 
def pVRR(A,norma,(la,ma,na),a1,
         B,normb,a2,
         C,normc,(lc,mc,nc),a3,
         D,normd,a4):
    " Obara Saika Recurrence "
    " comput [m0|n0](0) through vertical RR"
    
    def func1(x):
        return 0 if x<0 else x

    zeta = a1+a2
    eta = a3+a4
        
    
    P = numpy.zeros(3,'d')
    Q = numpy.zeros(3,'d')
    W = numpy.zeros(3,'d')

    P[0] = (a1*A[0] + a2*B[0])/zeta
    P[1] = (a1*A[1] + a2*B[1])/zeta
    P[2] = (a1*A[2] + a2*B[2])/zeta

    Q[0] = (a3*C[0] + a4*D[0])/eta
    Q[1] = (a3*C[1] + a4*D[1])/eta
    Q[2] = (a3*C[2] + a4*D[2])/eta
    
    W[0] = (zeta*P[0] + eta*Q[0])/(zeta+eta)
    W[1] = (zeta*P[1] + eta*Q[1])/(zeta+eta)
    W[2] = (zeta*P[2] + eta*Q[2])/(zeta+eta)
    
    AB2 = pow(A[0]-B[0],2) + pow(A[1]-B[1],2) + pow(A[2]-B[2],2)
    CD2 = pow(C[0]-D[0],2) + pow(C[1]-D[1],2) + pow(C[2]-D[2],2)
    PQ2 = pow(P[0]-Q[0],2) + pow(P[1]-Q[1],2) + pow(P[2]-Q[2],2)

    Kab = exp(-a1*a2/zeta*AB2)
    Kcd = exp(-a3*a4/eta *CD2)

    Theta = zeta*eta/(zeta+eta)

    Ntot = la+ma+na+lc+mc+nc
    Fn = numpy.zeros(Ntot+1,'d')
    Fn[Ntot] = Fm(Ntot,Theta*PQ2)
    "Downward recursion"
    for i in range(Ntot-1,-1,-1):
        Fn[i] = (2.*Theta*PQ2*Fn[i+1]+exp(-Theta*PQ2))/(2.*i+1)

    pre_fac = 2.*pow(pi,2.5)*Kab*Kcd/(zeta*eta*sqrt(zeta+eta))*\
              norma*normb*normc*normd

    PA = P-A
    QC = Q-C
    Xp = W-P
    Xq = W-Q

    tir_ = {}
    for m in range(Ntot+1):
        tir_[0,0,0,m] = (pre_fac*Fn[m])
        
    if la or lc:
        for i in range(la+1):
            if i:
                " in x, [m+1,0|0,0] "
                for m in range(Ntot-(i-1)):
                    tir_[0,i,0,m] = (PA[0]*tir_[0,i-1,0,m] + \
                                     Xp[0]*tir_[0,i-1,0,m+1])
                    if i-1>0:
                        tir_[0,i,0,m] += ((i-1)/2./zeta*\
                                          (tir_[0,i-2,0,m]- \
                                           eta/(zeta+eta)*tir_[0,i-2,0,m+1]))
            " in x, [m,0|n+1,0] "
            for j in range(func1(i+lc-la)):
                for m in range(Ntot-j-i):
                    tir_[0,i,j+1,m] = (QC[0]*tir_[0,i,j,m] + \
                                       Xq[0]*tir_[0,i,j,m+1])
                    if j:
                        tir_[0,i,j+1,m] += (j/2./eta*\
                                             (tir_[0,i,j-1,m]-\
                                              zeta/(zeta+eta)*tir_[0,i,j-1,m+1]))
                    if i:
                        tir_[0,i,j+1,m] += (i/2./(zeta+eta)*tir_[0,i-1,j,m+1])
                        
    for m in range(Ntot-la-lc+1):
        tir_[1,0,0,m] = tir_[0,la,lc,m]
        
    if ma or mc:
        for i in range(ma+1):
            if i:
                " in y, [m+1,0|0,0] "
                for m in range(Ntot-la-lc-(i-1)):
                    tir_[1,i,0,m] = (PA[1]*tir_[1,i-1,0,m] + \
                                     Xp[1]*tir_[1,i-1,0,m+1])
                    if i-1>0:
                        tir_[1,i,0,m] += ((i-1)/2./zeta*\
                                          (tir_[1,i-2,0,m]-\
                                           eta/(zeta+eta)*tir_[1,i-2,0,m+1]))
            " in y, [m,0|n+1,0] "
            for j in range(func1(i+mc-ma)):
                for m in range(Ntot-la-lc-j-i):
                    
                    tir_[1,i,j+1,m] = (QC[1]*tir_[1,i,j,m] + \
                                       Xq[1]*tir_[1,i,j,m+1])
                    if j:
                        tir_[1,i,j+1,m] += (j/2./eta*\
                                            (tir_[1,i,j-1,m]-\
                                             zeta/(zeta+eta)*tir_[1,i,j-1,m+1]))
                    if i:
                        tir_[1,i,j+1,m] += (i/2./(zeta+eta)*tir_[1,i-1,j,m+1])
                        
    for m in range(Ntot-la-lc-ma-mc+1):
        tir_[2,0,0,m] = tir_[1,ma,mc,m]
        
    if na or nc:
        for i in range(na+1):
            if i:
                " in z, [m+1,0|0,0] "
                for m in xrange(Ntot-la-lc-ma-mc-(i-1)):
                    tir_[2,i,0,m] = (PA[2]*tir_[2,i-1,0,m] + \
                                     Xp[2]*tir_[2,i-1,0,m+1])
                    if i-1>0:
                        tir_[2,i,0,m] += ((i-1)/2./zeta*\
                                          (tir_[2,i-2,0,m]-\
                                           eta/(zeta+eta)*tir_[2,i-2,0,m+1]))
            " in z, [m,0|n+1,0] "
            for j in xrange(func1(i+nc-na)):
                for m in xrange(Ntot-la-lc-ma-mc-j-i):
                    tir_[2,i,j+1,m] = (QC[2]*tir_[2,i,j,m] + \
                                       Xq[2]*tir_[2,i,j,m+1])
                    if j:
                        tir_[2,i,j+1,m] += (j/2./eta*\
                                             (tir_[2,i,j-1,m]-\
                                              zeta/(zeta+eta)*tir_[2,i,j-1,m+1]))
                    if i:
                        tir_[2,i,j+1,m] += (i/2./(zeta+eta)*tir_[2,i-1,j,m+1])
    
    return tir_[2,na,nc,0]
                        
