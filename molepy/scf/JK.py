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
from molepy.pymints.getints import Index

"JK"

def JK(Nbas,dens,J,K):
    G = numpy.empty((Nbas,Nbas),'d')
    dens1 = numpy.ravel(dens)
    for i in range(Nbas):
        for j in range(i+1):
            J2K = 2*J[i*(i+1)/2+j]-K[i*(i+1)/2+j]
            G[i,j] = numpy.dot(J2K,dens1)
            G[j,i] = G[i,j]
    return G

def Jints(eri,I,J,Nbas):
    jmat = numpy.empty(Nbas*Nbas,'d')
    ij = 0
    for k in range(Nbas):
        for l in range(Nbas):
            jmat[ij] = eri[Index(I,J,k,l)]
            ij += 1
    return jmat

def Kints(eri,I,J,Nbas):
    kmat = numpy.empty(Nbas*Nbas,'d')
    ij = 0
    for k in range(Nbas):
        for l in range(Nbas):
            kmat[ij] = eri[Index(I,k,J,l)]
            ij += 1
    return kmat
    
def makeJK(Nbas,eri):
    " J K matrices " 

    J = [[] for i in range(Nbas*(Nbas+1)/2)]
    K = [[] for i in range(Nbas*(Nbas+1)/2)]

    for i in range(Nbas):
        for j in range(i+1):
            J[i*(i+1)/2+j] = Jints(eri,i,j,Nbas)
            K[i*(i+1)/2+j] = Kints(eri,i,j,Nbas)
    return J,K
