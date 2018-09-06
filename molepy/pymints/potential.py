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
from molepy.pymints.boys import Fm

    

def os_Vint(bas1,bas2,coord):
    
    def onefunc(x):
        return 0 if x<0 else x
            
    VIJ = 0.0
    A = bas1.center
    B = bas2.center
    AB2 = pow(A[0]-B[0],2)+pow(A[1]-B[1],2)+pow(A[2]-B[2],2)

    l1,m1,n1 = bas1.cartam
    l2,m2,n2 = bas2.cartam
    mtot = l1+m1+n1+l2+m2+n2+0
    
    for p1 in range(bas1.nbas):
        a1 = bas1.exps[p1]
        c1 = bas1.coefs[p1]
        
        for p2 in range(bas2.nbas):
            a2 = bas2.exps[p2]
            c2 = bas2.coefs[p2]

            gamma = a1+a2
            Kab = exp(-a1*a2*AB2/gamma)
            pre_fac = (2.*pi/gamma)*Kab*c1*c2*bas1.norms[p1]*bas2.norms[p2]

            P = numpy.zeros(3,'d')
            PA = numpy.zeros(3,'d')
            PB = numpy.zeros(3,'d')

            P[0] = (a1*A[0] + a2*B[0])/gamma
            P[1] = (a1*A[1] + a2*B[1])/gamma
            P[2] = (a1*A[2] + a2*B[2])/gamma
            
            PA[0] = P[0] - A[0]
            PA[1] = P[1] - A[1]
            PA[2] = P[2] - A[2]
            PB[0] = P[0] - B[0]
            PB[1] = P[1] - B[1]
            PB[2] = P[2] - B[2]
            
            for atom in coord:
                C = atom.center

                PC = numpy.zeros(3,'d')
                
                PC[0] = P[0] - C[0]
                PC[1] = P[1] - C[1]
                PC[2] = P[2] - C[2]
                rpc2 = pow(PC[0],2)+pow(PC[1],2)+pow(PC[2],2)

                Fn = numpy.zeros(mtot+1,'d')
                Fn[mtot] = Fm(mtot,gamma*rpc2)
                for im in range(mtot-1,-1,-1):
                    Fn[im] = (2.*gamma*rpc2*Fn[im+1]+exp(-gamma*rpc2))/(2.*im+1)

                
                " OS Recurrence Begin "
                vi = {}
                for m in range(mtot+1): 
                    vi[0,0,0,m] = (pre_fac*Fn[m])
                #-----x
                if l1 or l2:
                    for i in range(l1+1):
                        if i:
                            # (1|0) (2|0) ...
                            for m in range(mtot-(i-1)):
                                vi[0,i,0,m] = (PA[0]*vi[0,i-1,0,m] - PC[0]*\
                                               vi[0,i-1,0,m+1])
                                if i-1>0:
                                    vi[0,i,0,m] += ((i-1)/2./gamma*\
                                                    (vi[0,i-2,0,m]-vi[0,i-2,0,m+1]))
                        # (i|1) (i|2) ...
                        for j in range(onefunc(i+l2-l1)):
                            for m in range(mtot-j-i):
                                vi[0,i,j+1,m] = (PB[0]*vi[0,i,j,m] - PC[0]*\
                                                 vi[0,i,j,m+1])
                                if j:
                                    vi[0,i,j+1,m] += (j/2./gamma*\
                                                      (vi[0,i,j-1,m]-vi[0,i,j-1,m+1]))
                                if i:
                                    vi[0,i,j+1,m] += (i/2./gamma*\
                                                      (vi[0,i-1,j,m]-vi[0,i-1,j,m+1]))
                #-----y
                for m in range(mtot-l1-l2+1):
                    vi[1,0,0,m] = vi[0,l1,l2,m]
                if m1 or m2:
                    for i in range(m1+1):
                        if i:
                            # (1|0) (2|0) ...
                            for m in range(mtot-l1-l2-(i-1)):
                                vi[1,i,0,m] = (PA[1]*vi[1,i-1,0,m] - PC[1]*vi[1,i-1,0,m+1])
                                if i-1>0:
                                    vi[1,i,0,m] += ((i-1)/2./gamma*\
                                                    (vi[1,i-2,0,m]-vi[1,i-2,0,m+1]))
                        # (i|1) (i|2) ...
                        for j in range(onefunc(i+m2-m1)):
                            for m in range(mtot-l1-l2-j-i):
                                vi[1,i,j+1,m] = (PB[1]*vi[1,i,j,m] - PC[1]*vi[1,i,j,m+1])
                                if j:
                                    vi[1,i,j+1,m] += (j/2./gamma*\
                                                      (vi[1,i,j-1,m]-vi[1,i,j-1,m+1]))
                                if i:
                                    vi[1,i,j+1,m] += (i/2./gamma*\
                                                      (vi[1,i-1,j,m]-vi[1,i-1,j,m+1]))
                #-----z
                for m in range(mtot-l1-l2-m1-m2+1):
                    vi[2,0,0,m] = vi[1,m1,m2,m]
                if n1 or n2:
                    for i in range(n1+1):
                        if i:
                            # (1|0) (2|0) ...
                            for m in range(mtot-l1-l2-m1-m2-(i-1)):
                                vi[2,i,0,m] = (PA[2]*vi[2,i-1,0,m] - PC[2]*vi[2,i-1,0,m+1])
                                if i-1>0:
                                    vi[2,i,0,m] += ((i-1)/2./gamma*\
                                                    (vi[2,i-2,0,m]-vi[2,i-2,0,m+1]))
                        # (i|1) (i|2) ...
                        for j in range(onefunc(i+n2-n1)):
                            for m in range(mtot-l1-l2-m1-m2-j-i):
                                vi[2,i,j+1,m] = (PB[2]*vi[2,i,j,m] - PC[2]*vi[2,i,j,m+1])
                                if j:
                                    vi[2,i,j+1,m] += (j/2./gamma*\
                                                      (vi[2,i,j-1,m]-vi[2,i,j-1,m+1]))
                                if i:
                                    vi[2,i,j+1,m] += (i/2./gamma*\
                                                      (vi[2,i-1,j,m]-vi[2,i-1,j,m+1]))
                " OS Recurrence End "
                VIJ += atom.anum*-vi[2,n1,n2,0]
    return bas1.norm*bas2.norm*VIJ
               
