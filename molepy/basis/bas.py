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
from math import sqrt,pi,pow

lquant = {'S':0,'P':1,'D':2,'F':3,'G':4}
# lexical ordering
cartAM_list = [[[0,0,0]],#s
               [[1,0,0],[0,1,0],[0,0,1]],#p
               [[2,0,0],[1,1,0],[1,0,1],[0,2,0],[0,1,1],[0,0,2]],#d
               [[3,0,0],[2,1,0],[2,0,1],[1,2,0],[1,1,0],[1,0,2],
                [0,3,0],[0,2,1],[0,2,1],[0,1,2],[0,0,3]]]#f

def get_basis(coord,bas_name):
    data = __import__(bas_name,globals(),locals(),['BasisSet'])
    bas_data = getattr(data,'BasisSet')


    BASIS = Get_Basis()
    
    basis = []
    for centers in coord:
        bfs = bas_data[centers.sym]
        for bf in bfs:
            L = lquant[bf[0]]
            for cartam in cartAM_list[L]:
                basis_set = GausBas(centers.center,cartam,centers.atom)
                for i in range(len(bf[1])):
                    basis_set.primitives(bf[1][i],bf[2][i])  
                basis_set.normalize()
                basis.append(basis_set)
                BASIS.AddBas(basis_set)
            BASIS.AddShells(L)
    return(basis,BASIS)


class Get_Basis:

    def __init__(self):
        self.CartBasis = []
        self.Shells = []
    def AddBas(self,BFS_):
        self.CartBasis.append(BFS_)
    def AddShells(self,shell):
        " Shells contains only l quantum number\
          in the same sequence as would be in\
          cartesian basis functions if they would\
          have the same number of functions"
        self.Shells.append(shell)

class GausBas:

    def __init__(self,center,cartam,atom):
 
        self.center = center
        self.cartam = cartam
        self.atom = atom
        self.exps = []
        self.coefs = [] 
        self.norms = []
        self.norm = 1.0
        self.nbas = 0

    def primitives(self,expn,coef): 
        l,m,n = self.cartam
        normz = pow(2./pi,0.75)*pow(2,l+m+n)*\
               pow(expn,(2*l+2*m+2*n+3)/4.)/\
               sqrt(dfact(2*l-1)*dfact(2*m-1)*dfact(2*n-1))
        self.exps.append(expn)
        self.coefs.append(coef)
        self.norms.append(normz)
        self.nbas += 1
        

    def normalize(self):
        
        l,m,n = self.cartam
        pre_fac = pow(pi,1.5)*dfact(2*l-1)*dfact(2*m-1)*\
                  dfact(2*n-1)/pow(2,l+m+n)
        sum = 0.0
        for i in range(self.nbas):
            for j in range(self.nbas):
                sum += self.coefs[i]*self.coefs[j]*\
                       self.norms[i]*self.norms[j]/\
                       pow(self.exps[i]+self.exps[j],(l+m+n)+1.5)

        tmp_norm = 1.0/sqrt(pre_fac * sum)
        for i in range(self.nbas):
            self.coefs[i] *= tmp_norm
        
def dfact(i):    
    return 1 if (i<=1) else i*dfact(i-2)
