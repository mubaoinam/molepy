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

atomic_num = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4,
              'B': 5, 'C': 6, 'N': 7, 'O': 8,
              'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12,
              'Al': 13, 'Si': 14, 'P': 15, 'S': 16,
              'Cl': 17, 'Ar': 18}

def Mole(centers):
    coord = []
    cout = 1
    for s,x,y,z in centers:
        atms = Atoms(s,x,y,z,cout)
        coord.append(atms)
    return coord

class Atoms:
    def __init__(self,s,x,y,z,atom):
        self.sym = s
        self.center = (x,y,z)
        self.anum = atomic_num[s]
        self.atom = atom
