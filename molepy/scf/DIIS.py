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

" Pulay's DIIS "

class DIIS:
    def __init__(self):
        self.Fvec = []
        self.Evec = []
        self.nim = 0
        
    def diis(self, Fock, Dens ,S ):
        
        " Error vectors "
        Err = reduce(numpy.dot, (Fock,Dens,S)) - \
              reduce(numpy.dot, (S,Dens,Fock))
        Err = numpy.ravel(Err)

        self.Fvec.append(Fock)
        self.Evec.append(Err)
        self.nim += 1

        if self.nim > 1:
            i,j = Fock.shape
            Nvec = numpy.zeros((i,j),dtype=numpy.float64)
            
            B = numpy.zeros((self.nim+1,self.nim+1),
                            dtype=numpy.float64)
            A = numpy.zeros(self.nim+1,dtype=numpy.float64)
            for i in range(self.nim):
                for j in range(self.nim):
                    B[i,j] = numpy.dot(self.Evec[i],self.Evec[j])
                B[self.nim,i] = B[i,self.nim] = -1.0
            A[self.nim] = -1.0
            C = numpy.dot(numpy.linalg.pinv(B),A)
            for i in range(self.nim):
                Nvec += C[i]*self.Fvec[i]
            return Nvec
            
        else:
            return Fock
