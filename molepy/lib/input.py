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
from molepy.lib.mole_param import hfparam

def input_scan(inp):
    r=open(inp,'r')
    geom = False
    bohr = True
    cart_coord = []
    c1 = 1.0
    basis = 'vdz'
    for line in r.readlines():
    
        if 'geometry' in line and not '!' in line:
            geom = True
            if line.split()[1] == 'angstrom':
                bohr = False
            continue
    
        if geom:
            if line.isspace():
                geom=False
                continue
            s=line.split()
            if not bohr:
                c1=1.889725989
            cart_coord.append((s[0],c1*float(s[1]),c1*float(s[2]),c1*float(s[3])))
    
        if 'hf' in line and not '!' in line:
            s=line.split(',')
            for i in s[1:]:
                s1 = i.split('=')
                if s1[0] == 'basis':
                    basis = s1[1]
                elif s1[0] == 'energy':
                    hfparam.tolE = float(s1[1])
                elif s1[0] == 'density':
                    hfparam.tolD = float(s1[1])
                else:
                    print 'Parameter not recognised, proceeding without it!'
    r.close()
    return(cart_coord,basis)
                
