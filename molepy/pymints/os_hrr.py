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
import numpy,time
from math import exp,pi,pow,sqrt
from molepy.pymints.os_vrr import cVRR
from molepy.ccmints.os_hrr import ccHRR
from molepy.lib.mole_param import hfparam


def twofunc_(A,norma,Aexps,Acoefs,I,J,K,
             B,normb,Bexps,Bcoefs,
             C,normc,Cexps,Ccoefs,lc,mc,nc,
             D,normd,Dexps,Dcoefs,ld,md,nd):
    " Part of Horizontal RR "
    " obtain (a0|cd) "

    tmp_ = {}
    for i in range(0,nd+2):
        for j in range(lc+ld,lc-1,-1):
            for k in range(mc+md,mc-1,-1):
                for l in range(nc+nd-i,nc-1,-1):
                    if i==0:
                        tmp_[j,k,l,0,0,0] = (
                            cVRR(A,norma,Aexps,Acoefs,(I,J,K),
                                 B,normb,Bexps,Bcoefs,
                                 C,normc,Cexps,Ccoefs,(j,k,l),
                                 D,normd,Dexps,Dcoefs))
                    else:
                        tmp_[j,k,l,0,0,i] = (tmp_[j,k,l+1,0,0,i-1]+\
                                             (C[2]-D[2])*tmp_[j,k,l,0,0,i-1])
    for i in range(1,md+1):
        for j in range(lc+ld,lc-1,-1):
            for k in range(mc+md-i,mc-1,-1):
                tmp_[j,k,nc,0,i,nd] = (tmp_[j,k+1,nc,0,i-1,nd]+\
                                       (C[1]-D[1])*tmp_[j,k,nc,0,i-1,nd])
    for i in range(1,ld+1):
        for j in range(lc+ld-i,lc-1,-1):
            tmp_[j,mc,nc,i,md,nd] = (tmp_[j+1,mc,nc,i-1,md,nd]+\
                                     (C[0]-D[0])*tmp_[j,mc,nc,i-1,md,nd])
          
    return tmp_[lc,mc,nc,ld,md,nd]

def cHRR(a,b,c,d):
    if hfparam.ccmints:
        return ccHRR(a.center[0],a.center[1],a.center[2],a.norms,
                     a.cartam[0],a.cartam[1],a.cartam[2],a.exps,a.coefs,
                     b.center[0],b.center[1],b.center[2],b.norms,
                     b.cartam[0],b.cartam[1],b.cartam[2],b.exps,b.coefs,
                     c.center[0],c.center[1],c.center[2],c.norms,
                     c.cartam[0],c.cartam[1],c.cartam[2],c.exps,c.coefs,
                     d.center[0],d.center[1],d.center[2],d.norms,
                     d.cartam[0],d.cartam[1],d.cartam[2],d.exps,d.coefs)
    else:
        return cHRR_py(a,b,c,d)

     
def cHRR_py(bas1,bas2,bas3,bas4):
    " Obara Saika Recurrence "
    " comput (ab|cd) through horizontal RR "
    " using (a0|c0) from vertical RR "

    A = bas1.center
    B = bas2.center
    C = bas3.center
    D = bas4.center

    norma = bas1.norms
    normb = bas2.norms
    normc = bas3.norms
    normd = bas4.norms

    Aexps = bas1.exps
    Bexps = bas2.exps
    Cexps = bas3.exps
    Dexps = bas4.exps

    Acoefs = bas1.coefs
    Bcoefs = bas2.coefs
    Ccoefs = bas3.coefs
    Dcoefs = bas4.coefs

    (la,ma,na) = bas1.cartam
    (lb,mb,nb) = bas2.cartam
    (lc,mc,nc) = bas3.cartam
    (ld,md,nd) = bas4.cartam
    
    hrr_ = {}
    for i in range(la,la+lb+1):
        for j in range(ma,ma+mb+1):
            for k in range(na,na+nb+1):
                hrr_[i,j,k,0,0,0] = (
                    twofunc_(A,norma,Aexps,Acoefs,i,j,k,
                             B,normb,Bexps,Bcoefs,
                             C,normc,Cexps,Ccoefs,lc,mc,nc,
                             D,normd,Dexps,Dcoefs,ld,md,nd,))
    for i in xrange(1,nb+1):
        for j in xrange(la,la+lb+1):
            for k in xrange(ma,ma+mb+1):
                for l in xrange(na,na+nb+1-i):
                    hrr_[j,k,l,0,0,i] = (
                        hrr_[j,k,l+1,0,0,i-1]+(A[2]-B[2])*hrr_[j,k,l  ,0,0,i-1])
    for i in xrange(1,mb+1):
        for j in xrange(la,la+lb+1):
            for k in xrange(ma,ma+mb+1-i):
                hrr_[j,k,na,0,i,nb] = (
                    hrr_[j,k+1,na,0,i-1,nb]+(A[1]-B[1])*hrr_[j,k  ,na,0,i-1,nb])
    for i in xrange(1,lb+1):
        for j in xrange(la,la+lb+1-i):
            hrr_[j,ma,na,i,mb,nb] = (
                hrr_[j+1,ma,na,i-1,mb,nb]+(A[0]-B[0])*hrr_[j  ,ma,na,i-1,mb,nb])
    
    return hrr_[la,ma,na,lb,mb,nb]

    
