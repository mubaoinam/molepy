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
from molepy.lib.mole_param import intparam
from molepy.pymints.os_hrr import cHRR
def Index(i,j,k,l):
    if i<j:
        tmp = i
        i = j
        j = tmp
    if k<l:
        tmp = k
        k = l
        l = tmp

    ij = i*(i+1)/2+j
    kl = k*(k+1)/2+l
    if ij < kl:
        tmp = ij
        ij = kl
        kl = tmp
    return ij*(ij+1)/2+kl

def flatit(coord):
    natom = len(coord)
    Coords = [0.0 for i in range(natom*3)]
    Anums = [0 for i in range(natom)]
    for i in range(natom):
        Coords[i*3+0] = coord[i].center[0]
        Coords[i*3+1] = coord[i].center[1]
        Coords[i*3+2] = coord[i].center[2]
        Anums[i] = coord[i].anum
    return Coords,Anums
    
def getints1e(basis,coord):
    " S and Hcore "
    
    nbas = len(basis)
    S = numpy.empty((nbas,nbas),dtype=numpy.float64)
    Hcore = numpy.empty((nbas,nbas),dtype=numpy.float64)

    #--tmp
    Coords,Anums = flatit(coord)
    for i in range(nbas):
        for j in range(nbas):
            #S[i,j] = os_overlap(basis[i],basis[j])
            S[i,j] = intparam.os_overlap(basis[i].center[0],basis[i].center[1],basis[i].center[2],
                                  basis[i].norms,basis[i].exps,basis[i].coefs,
                                  basis[i].cartam[0],basis[i].cartam[1],basis[i].cartam[2],
                                  basis[j].center[0],basis[j].center[1],basis[j].center[2],
                                  basis[j].norms,basis[j].exps,basis[j].coefs,
                                  basis[j].cartam[0],basis[j].cartam[1],basis[j].cartam[2])
            
            #Hcore[i,j] = os_kinetic(basis[i],basis[j])
            Hcore[i,j] = intparam.os_kinetic(basis[i].center[0],basis[i].center[1],basis[i].center[2],
                                  basis[i].norms,basis[i].exps,basis[i].coefs,
                                  basis[i].cartam[0],basis[i].cartam[1],basis[i].cartam[2],
                                  basis[j].center[0],basis[j].center[1],basis[j].center[2],
                                  basis[j].norms,basis[j].exps,basis[j].coefs,
                                  basis[j].cartam[0],basis[j].cartam[1],basis[j].cartam[2])
            #Hcore[i,j] += os_Vint(basis[i],basis[j],coord)
            
            Hcore[i,j] += intparam.os_Vint(basis[i].center[0],basis[i].center[1],basis[i].center[2],
                                  basis[i].norms,basis[i].exps,basis[i].coefs,
                                  basis[i].cartam[0],basis[i].cartam[1],basis[i].cartam[2],
                                  basis[j].center[0],basis[j].center[1],basis[j].center[2],
                                  basis[j].norms,basis[j].exps,basis[j].coefs,
                                  basis[j].cartam[0],basis[j].cartam[1],basis[j].cartam[2],
                                  Coords,Anums)
                                
            
            
    return(S,Hcore)


def getints2e(basis):

    nbas = len(basis)
    tbas = nbas*(nbas+1)*(nbas*nbas+nbas+2)/8    
    eri = numpy.empty(tbas,'d')
    ttt = 0.0
    for i in range(nbas):
        for j in range(i+1):
            ij = i*(i+1)/2+j
            for k in range(nbas):
                for l in range(k+1):
                    kl = k*(k+1)/2+l
                    if ij<=kl:
                        eri[Index(i,j,k,l)] = cHRR(basis[i],basis[j],
                                                   basis[k],basis[l])
    return eri
