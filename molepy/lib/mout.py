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

class molog:
    filname = ''
    def __init__(self):
        apple = 0

    def set_name(self,string):
        molog.filname = string

    def initiate(self):
        out = open(molog.filname,'w')
        out.write('               ******************\n')
        out.write('               ***** MOLEPY *****\n')
        out.write('               ******************\n')
        out.write('\n')
        out.close()

    def writer(self,str1,*args):
        out = open(molog.filname,'a')
        out.write(str1.format(*args))
        out.close()
