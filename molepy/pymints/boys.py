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
from math import pow,exp,log


# c++ counterpart has pre-computed Euler gammas


def Fm(m,x):
    " Boys function computed from "
    " Euler gamma function and "
    " Incomplete gamma function "
    x = max(abs(x),0.00000001)
    incomplete = gammp(m+0.5,x)
    gamma = exp(gammln(m+0.5))

    boys = gamma*incomplete*0.5*pow(x,-m-0.5)

    return boys

" gammp gser gammln gcf from numerical recipe "
    
def gammp(a,x):
    
    assert (x>0. and a>=0.), "Exit from incomplete gamma!"
    
    if x == 0.0:
        return 0.0
    if x < a+1.0:
        return gser(a,x)
    else:
        return 1.0-gcf(a,x)

def gser(a,x):
    gln = gammln(a)
    ap = a
    del_ = sum = 1.0/a
    EPS = 3.e-7
    for i in range(100):
        ap += 1
        del_ *= x/ap
        sum += del_
        if abs(del_) < abs(sum)*EPS:
            break
    return sum*exp(-x+a*log(x)-gln)
 
    
def gammln(xx):
    cof = [57.1562356658629235,-59.5979603554754912,
           14.1360979747417471,-0.491913816097620199,
           .339946499848118887e-4,.465236289270485756e-4,
           -.983744753048795646e-4,.158088703224912494e-3,
           -.210264441724104883e-3,.217439618115212643e-3,
           -.164318106536763890e-3,.844182239838527433e-4,
           -.261908384015814087e-4,.368991826595316234e-5]
    y=x=xx
    tmp = x+5.24218750000000000
    tmp = (x+0.5)*log(tmp)-tmp
    ser = 0.999999999999997092
    for j in range(14):
        y += 1
        ser += cof[j]/y
    return tmp+log(2.5066282746310005*ser/x)

def gcf(a,x):

    EPS=3.e-7
    FPMIN=1.e-30    
    gln=gammln(a)
    b=x+1.0-a
    c=1.0/FPMIN
    d=1.0/b
    h=d
    for i in range(1,100):
        an = -i*(i-a)
        b += 2.0
        d=an*d+b
        if (abs(d) < FPMIN): d=FPMIN
        c=b+an/c
        if (abs(c) < FPMIN): c=FPMIN
        d=1.0/d
        del_=d*c
        h *= del_
        if (abs(del_-1.0) <= EPS): break
                 
    return exp(-x+a*log(x)-gln)*h
    

