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
from math import pow,sqrt,factorial

" Transform Integrals over cartesian gaussian basis functions "
" to spherical harmonics "

lcart = [[[0,0,0]],#s
         [[1,0,0],[0,1,0],[0,0,1]],#p
         [[2,0,0],[1,1,0],[1,0,1],[0,2,0],[0,1,1],[0,0,2]],#d
         [[3,0,0],[2,1,0],[2,0,1],[1,2,0],[1,1,0],[1,0,2],
          [0,3,0],[0,2,1],[0,2,1],[0,1,2],[0,0,3]]]#f

mquant = [[0],[-1,0,1],[-2,-1,0,1,2],
          [-3,-2,-1,0,1,2,3],
          [-4,-3,-2,-1,0,1,2,3,4]]


# pre-computed coeffs comes from c++ counterpart, python generates
# slightly different numbers.

coeff_0 = [ 1.000000000000000 ]
coeff_1 = [0.000000000000000,1.000000000000000,0.000000000000000,
           0.000000000000000,0.000000000000000,1.000000000000000,
           1.000000000000000,0.000000000000000,0.000000000000000]
coeff_2 = [0.000000000000000,1.000000000000000,0.000000000000000,
           0.000000000000000,0.000000000000000,0.000000000000000,
           0.000000000000000,0.000000000000000,0.000000000000000,
           0.000000000000000,1.000000000000000,0.000000000000000,
           -0.500000000000000,0.000000000000000,0.000000000000000,
           -0.500000000000000,0.000000000000000,1.000000000000000,
           0.000000000000000,0.000000000000000,1.000000000000000,
           0.000000000000000,0.000000000000000,0.000000000000000,
           0.866025403784439,0.000000000000000,0.000000000000000,
           -0.866025403784439,0.000000000000000,0.000000000000000]
coeff_3 = [0.000000000000000,1.060660171779821,0.000000000000000,
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
           0.000000000000000]

CScoeff = [ coeff_0, coeff_1, coeff_2, coeff_3]





# two center
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

def dfact(n):
    return 1 if (n<=1) else n*dfact(n-2)

def binomial(x,y):
    return factorial(x)/factorial(y)/factorial(x-y)

def coeff_compute(l,m,lx,ly,lz):

    t = l-lz-abs(m)
    p = lx-m if m>=0 else lx+m-1

    if t%2 or p%2:
        return 0.000000e0

    t /= 2
    p /= 2


    pfac = numpy.longdouble(sqrt( factorial(l+abs(m))*factorial(l-abs(m))))
    pfac *= 1.000000e0 if m==0 else \
            1.4142135623730950488016887e0
    pfac /= numpy.longdouble(pow(2,abs(m))*factorial(l))
    pfac *= numpy.longdouble(binomial(l,t)*binomial(l-t,abs(m)+t))/pow(4,t)

    sum =0.0000000000e0
    for u in range(t+1):
        sum1 = numpy.longdouble(binomial(t,u)*binomial(abs(m),ly-u-u))
        if (ly-u+p)%2:
            sum -= sum1
        else:
            sum += sum1
    sum *= pfac
    sum *= sqrt(numpy.longdouble(dfact(2*lx-1)*dfact(2*ly-1)*dfact(2*lz-1))/\
                dfact(2*l-1))
    return sum

def Cart2Spher2(shells,CartInt,ncart):
    nsph = 0
    for i in shells:
        nsph += 2*i+1
    SphInt = numpy.zeros(nsph*nsph,'d')
    sph1 = 0
    s1ind = 0
    for l1 in shells:
        for im1,m1 in enumerate(mquant[l1]):
            sph2 = 0
            s2ind = 0
            for l2 in shells:
                for im2,m2 in enumerate(mquant[l2]):
                    
                    " Cartesian contributions"
                    c1ind = sph1
                    for ic1,c1 in enumerate(lcart[l1]):
                        c2ind = sph2
                        coeff1 = coeff_compute(l1,m1,c1[0],c1[1],c1[2]) if l1>3 else\
                                 CScoeff[l1][(im1)*len(lcart[l1])+(ic1)]
                        for ic2,c2 in enumerate(lcart[l2]):
                            coeff2 = coeff_compute(l2,m2,c2[0],c2[1],c2[2]) if l2>3 else \
                                     CScoeff[l2][(im2)*len(lcart[l2])+(ic2)]
                            SphInt[s1ind*nsph+s2ind] += coeff1*coeff2* \
                                                    CartInt[c1ind*ncart+c2ind]   
                            c2ind += 1
                        c1ind += 1
                    s2ind += 1
                sph2 += len(lcart[l2])
            s1ind += 1
        sph1 += len(lcart[l1])
    return SphInt

def Cart2Spher4_save(shells,CartInt):
    print ' Starting 4 index transformation '
    start_time = time.time()
    nsph = 0
    for i in shells:
        nsph += 2*i+1
    tbas = nsph*(nsph+1)*(nsph*nsph+nsph+2)/8
    SphInt = numpy.zeros(tbas,'d')

    sph1 = 0
    s1ind = 0
    for l1 in shells:
        for m1 in mquant[l1]:

            sph2 = 0
            s2ind = 0
            for l2 in shells:
                for m2 in mquant[l2]:

                    sph3 = 0
                    s3ind = 0
                    for l3 in shells:
                        for m3 in mquant[l3]:

                            sph4 = 0
                            s4ind = 0
                            for l4 in shells:
                                for m4 in mquant[l4]:

                                    " Cartesian contributions "
                                    c1ind = sph1
                                    for c1 in lcart[l1]:

                                        c2ind = sph2
                                        coeff1 = coeff_compute(l1,m1,c1[0],c1[1],c1[2])
                                        for c2 in lcart[l2]:

                                            c1c2 = c1ind*(c1ind+1)/2+c2ind
                                            c3ind = sph3
                                            coeff2 = coeff_compute(l2,m2,c2[0],c2[1],c2[2])
                                            for c3 in lcart[l3]:

                                                c4ind = sph4
                                                coeff3 = coeff_compute(l3,m3,c3[0],c3[1],c3[2])
                                                for c4 in lcart[l4]:

                                                    c3c4 = c3ind*(c3ind+1)/2+c4ind
                                                    if (c1ind >= c2ind and c3ind >= c4ind and c1c2 <= c3c4):
                                                        coeff4 = coeff_compute(l4,m4,c4[0],c4[1],c4[2])
                                                        SphInt[Index(s1ind,s2ind,s3ind,s4ind)] = coeff1*coeff2*coeff3*coeff4*\
                                                                                                 CartInt[Index(c1ind,c2ind,
                                                                                                               c3ind,c4ind)]
                                                    c4ind += 1
                                                c3ind += 1
                                            c2ind += 1
                                        c1ind += 1
                                    s4ind += 1                
                                sph4 += len(lcart[l4])
                            s3ind += 1
                        sph3 += len(lcart[l3])
                    s2ind += 1
                sph2 += len(lcart[l2])
            s1ind += 1
        sph1 += len(lcart[l1])
    print("--- %s seconds ---" % (time.time()-start_time))
    print len(SphInt)
    return SphInt


def Cart2Spher4(shells,CartInt):
    start_time = time.time()
    nsph = 0
    for i in shells:
        nsph += 2*i+1
    tbas = nsph*(nsph+1)*(nsph*nsph+nsph+2)/8
    SphInt = numpy.zeros(tbas,'d')    
    sph1 = 0
    s1ind = 0 # index of 1st iter
    shell1 = 0
    " 1st loop for index 'i' "
    for l1 in shells:
        for im1,m1 in enumerate(mquant[l1]):    # 1st iter on basis -sph similr to for i in nbas
            

            sph2 = 0
            s2ind = 0
            "2nd loop for index 'j' "
            for l2 in shells[:shell1+1]:
                for im2,m2 in enumerate(mquant[l2]):   # 2nd iter similar to for j in i+1

                    if s2ind > s1ind: continue
                    m1m2 =  s1ind*(s1ind+1)/2+s2ind # similar to ij
                    
                    sph3 = 0
                    s3ind = 0
                    shell3 = 0
                    "3rd loop for index 'k' "
                    for l3 in shells: 
                        for im3,m3 in enumerate(mquant[l3]):  # 3rd iter similar to for k in nbas

                            sph4 = 0
                            s4ind = 0
                            " 4th loop for index 'l' "
                            for l4 in shells[:shell3+1]:
                                for im4,m4 in enumerate(mquant[l4]):  # 4th iter similar to for l in k+1:
                                    if s4ind > s3ind: continue
                                    m3m4 = s3ind*(s3ind+1)/2+s4ind # similar to kl
                                    if m1m2 <= m3m4:
                                        """
                                        loop c1 in lcart[l1] ---loopA
                                        loop c2 in lcart[l2] ---loopB
                                        loop c3 in lcart[l3] ---loopC
                                        loop c4 in lcart[l4] ---loopD
                                        """
                                        
                                        " Cartesian contributions "
                                        c1ind = sph1
                                        for ic1,c1 in enumerate(lcart[l1]):
                                        
                                            c2ind = sph2
                                            coeff1 = coeff_compute(l1,m1,c1[0],c1[1],c1[2]) if l1>3 else \
                                                     CScoeff[l1][(im1)*len(lcart[l1])+(ic1)]
                                            for ic2,c2 in enumerate(lcart[l2]):
                                        
                                                c3ind = sph3
                                                coeff2 = coeff_compute(l2,m2,c2[0],c2[1],c2[2]) if l2>3 else \
                                                         CScoeff[l2][(im2)*len(lcart[l2])+(ic2)]
                                                for ic3,c3 in enumerate(lcart[l3]):
                                        
                                                    c4ind = sph4
                                                    coeff3 = coeff_compute(l3,m3,c3[0],c3[1],c3[2]) if l3>3 else \
                                                             CScoeff[l3][(im3)*len(lcart[l3])+(ic3)]
                                                    
                                                    for ic4,c4 in enumerate(lcart[l4]):
                                        
                                                        coeff4 = coeff_compute(l4,m4,c4[0],c4[1],c4[2]) if l4 > 3 else \
                                                                 CScoeff[l4][im4*len(lcart[l4])+ic4]
                                                        SphInt[Index(s1ind,s2ind,s3ind,s4ind)] += coeff1*coeff2*coeff3*coeff4*\
                                                                                            CartInt[Index(c1ind,c2ind,c3ind,c4ind)]
                                                        c4ind += 1
                                                    c3ind += 1
                                                c2ind += 1
                                            c1ind += 1
                                    s4ind += 1 # index 'l' from 4th iter -sph
                                sph4 += len(lcart[l4])   # starting index 'l' in cartesian    
                            s3ind += 1 # index 'k' from 3rd iter -sph
                        sph3 += len(lcart[l3])  # starting index 'k' in cartesian
                        shell3 += 1
                    s2ind += 1 # index 'j' from 2nd iter -sph
                sph2 += len(lcart[l2]) # starting index 'j' in cartesian
            s1ind += 1    # index 'i' from 1st iter -sph
        sph1 += len(lcart[l1]) # starting index 'i' in cartesian
        shell1 += 1
    return SphInt

            

