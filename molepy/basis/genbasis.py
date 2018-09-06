#!/usr/bin/env python
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
import sys,os,glob
from math import *
from numpy import *
r = open(sys.argv[1],'r')

w = open(sys.argv[2],'w')
atoms = ['HYDROGEN', 'HELIUM', 'LITHIUM', 'BERYLLIUM', 'BORON', 'CARBON', 'NITROGEN', 'OXYGEN', 'FLUORINE', 'NEON', 'SODIUM', 'MAGNESIUM', 'ALUMINIUM', 'SILICON', 'PHOSPHORUS' ,'SULFUR', 'CHLORINE', 'ARGON']

atm = {'HYDROGEN':'\'H\'', 'HELIUM':'\'He\'', 'LITHIUM': '\'Li\'', 'BERYLLIUM': '\'Be\'', 'BORON': '\'B\'', 'CARBON':'\'C\'', 'NITROGEN':'\'N\'', 'OXYGEN':'\'O\'', 'FLUORINE':'\'F\'', 'NEON':'\'Ne\'', 'SODIUM':'\'Na\'', 'MAGNESIUM':'\'Mg\'', 'ALUMINIUM':'\'Al\'', 'SILICON':'\'Si\'', 'PHOSPHORUS':'\'P\'', 'SULFUR':'\'S\'', 'CHLORINE':'\'Cl\'', 'ARGON':'\'Ar\''}

#atoms = ['HYDROGEN','CARBON','HELIUM','OXYGEN']
#atm = {'HYDROGEN':'\'H\'','CARBON':'\'C\'','HELIUM':'\'He\'','OXYGEN':'\'O\''}
start = False
roll = 0
cout2 = True
w.write('BasisSet = {\n')
for line in r.readlines():
    if not line.isspace():
        if line.replace('\n','') in atoms:
            if not cout2:
                w.write(',\n')
            cout2 = False
            start = True
            w.write('{:<4}{:>3}'.format(atm[line.replace('\n','')],':['))
            cout1 = 0
            continue
    else:
        if start:
            start = False
            w.write(']')
            continue

    if start:
        if roll == 0:
            if '$END' in line:
                w.write(']')
                break
            if not cout1 == 0:
                w.write(',\n')
                w.write('       ')
            s = line.split()
            roll = int(s[-1])
            sym = s[0]
            exp = []
            coef = []
            w.write('(\''+sym+'\',')
            continue
        if roll>0:
            s = line.split()
            exp.append(float(s[1].replace('D','e')))
            coef.append(float(s[2].replace('D','e')))
            roll -= 1
        if roll == 0:
            w.write(str(exp)+',')
            w.write(str(coef)+')')
        cout1 = 1
            
w.write('}\n')
