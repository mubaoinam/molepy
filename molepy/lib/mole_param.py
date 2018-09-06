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


class hf_param_:
    
    def __init__(self):
        ccmints = None
        Spheri = None
        Cart2Spher4 = None
        Cart2Spher2 = None
        tolE = 0.0
        tolD = 0.0

class mints_param_:

    def __init__(self):
        os_Vint = None
        os_overlap = None
        os_kinetic = None


hfparam = hf_param_()
hfparam.ccmints = True
hfparam.Spheri = True
hfparam.tolE = 1e-07
hfparam.tolD = 1e-07
intparam = mints_param_()


if hfparam.ccmints:
    from molepy.ccmints.CartSpher import Cart2Spher4,Cart2Spher2
    from molepy.ccmints.potential import os_Vint
    from molepy.ccmints.overlap import os_overlap
    from molepy.ccmints.kinetic import os_kinetic
    hfparam.Cart2Spher4 = Cart2Spher4
    hfparam.Cart2Spher2 = Cart2Spher2
    intparam.os_Vint = os_Vint
    intparam.os_overlap = os_overlap
    intparam.os_kinetic = os_kinetic
else:
    from molepy.pymints.CartSpher import Cart2Spher4,Cart2Spher2
    from molepy.pymints.potential import os_Vint
    from molepy.pymints.overlap import os_overlap
    from molepy.pymints.kinetic import os_kinetic
    hfparam.Cart2Spher4 = Cart2Spher4
    hfparam.Cart2Spher2 = Cart2Spher2
    intparam.os_Vint = os_Vint
    intparam.os_overlap = os_overlap
    intparam.os_kinetic = os_kinetic
