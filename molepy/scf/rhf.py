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
from math import sqrt,pow
from numpy.core.umath_tests import inner1d
from molepy.basis.bas import get_basis
from molepy.pymints.getints import getints1e,getints2e
from molepy.scf.orth import orthog
from molepy.scf.DIIS import DIIS
from molepy.scf.JK import JK,makeJK
from molepy.lib import mout
from molepy.lib.mole_param import hfparam

def flatit(Aa,Nn):
    Bb = [0.0 for i in range(Nn*Nn)]
    for i in range(Nn):
        for j in range(Nn):
            Bb[i*Nn+j] = Aa[i][j]
    return Bb

def deflatit(Aa,Nn):
    Bb = numpy.empty((Nn,Nn),'d')
    for i in range(Nn):
        for j in range(Nn):
            Bb[i][j] = Aa[i*Nn+j]
    return Bb

def RHF(coord,bfs):
    " Restricted Hartree Fock "

    out = mout.molog()
    out.initiate()
    out.writer('{:} \n',' Program : Restricted Hartree Fock (closed shell)')
    out.writer('{:} \n','           O.R. Meitei')
    out.writer('\n','')
    out.writer('{:} \n',' Atomic Coordinates')
    out.writer('\n','')
    for i in coord:
        out.writer(' {:<3}   {:>10.6f}   {:>10.6f}   {:>10.6f}\n',i.sym,*i.center)
    out.writer('\n','')
    out.writer(' Basis Set                 : {:>10}\n',bfs)
    
    
    basis,BASIS = get_basis(coord,bfs)
    #-tmp
    NSPH = 0
    for i in BASIS.Shells:
        NSPH += 2*i+1    
    aos = 0
    nelec = 0
    for i in coord:
        nelec += i.anum
    out.writer(' Number of electrons       : {:>7} \n',nelec)
    if hfparam.Spheri:
        n_basis_ = NSPH
    else:
        n_basis_ = len(basis)
    out.writer(' Number of Basis functions : {:>7} \n',n_basis_)
    out.writer('\n','')
    
    e_nuclear = nuclear(coord)
    S,h = getints1e(basis,coord)
    eri = getints2e(basis)
    
    sval,svec = numpy.linalg.eigh(S)

    if hfparam.Spheri:
        S = flatit(S,len(basis))
        h = flatit(h,len(basis))
        S = hfparam.Cart2Spher2(BASIS.Shells,S,len(basis))
        h = hfparam.Cart2Spher2(BASIS.Shells,h,len(basis))
        
        S = deflatit(S,NSPH)
        h = deflatit(h,NSPH)
        
        eri = hfparam.Cart2Spher4(BASIS.Shells,eri)

    
    out.writer(' Nuclear repulsion energy : {:>10.7f} \n',e_nuclear)
    out.writer('\n','')
    out.writer(' Minimum eigenvalue in the overlap matrix : {:>.4e} \n',sval[0])
    out.writer('\n','')
    
    X = orthog(S)
    es1,vs1 = numpy.linalg.eigh(S)
    orb_e,orb_c = m_eigh(h,X)
    Nocc = nclosed(coord) 

    Maxiter = 50
    diis = DIIS()
    out.writer(' Initial guess is core Hamiltonian (one electron) \n')
    out.writer(' Interpolation using DIIS \n')
    out.writer('\n','')
    out.writer(' Iteration       Energy          Ediff          Ddiff  \n')
    out.writer('\n','')
    
    prevE = 0.0
    prevD = 0.0
    J,K = makeJK(n_basis_,eri)
        
    for iter in range(Maxiter):

        Coefs = orb_c[:,0:Nocc] 
        Dens = numpy.dot(Coefs,Coefs.T)
        
        Gmat = JK(n_basis_,Dens,J,K) 
        Fock = h+Gmat
        Fock = diis.diis(Fock,Dens,S)

        orb_e,orb_c = m_eigh(Fock,X)

        etot,e1e,e2e = HFenergy(h,Fock,Dens,e_nuclear)

        
        rmsD = Dens-prevD
        rmsD = numpy.linalg.norm(rmsD)
        delE = abs(etot-prevE)
        out.writer('  {:>3}      {:>15.10f}     {:>.4e}     {:>.4e} \n',\
                  iter+1,etot,delE,rmsD)
        if iter:
            if delE < hfparam.tolE and rmsD < hfparam.tolD:
                " Call it converged "
                break

        prevE = etot
        prevD = Dens
    out.writer('\n','')
    out.writer(' SCF converged \n')
    out.writer('\n','')
    out.writer(' Final Energy        : {:>15.10f} \n',etot)
    out.writer(' Nuclear Energy      : {:>15.10f} \n',e_nuclear)
    out.writer(' One electron Energy : {:>15.10f} \n',e1e)
    out.writer(' Two electron Energy : {:>15.10f} \n',e2e)
    out.writer('\n','')
    return(etot,orb_e,orb_c)


def m_eigh(A,B):
    
    BtAB = reduce(numpy.dot, (B.T,A,B))
    eval1,evec1 = numpy.linalg.eigh(BtAB)
    evec2 = numpy.dot(B,evec1)

    return(eval1,evec2)

def nuclear(coord):
    enuclear = 0.0
    for i in range(len(coord)):
        A = coord[i].center
        for j in range(i):
            B = coord[j].center
            AB = sqrt(pow(A[0]-B[0],2)+pow(A[1]-B[1],2)+pow(A[2]-B[2],2))
            enuclear += coord[i].anum*coord[j].anum/AB
    return enuclear

def nclosed(coord):
    Nocc = 0
    for i in coord:
        Nocc += i.anum
    return Nocc//2

def HFenergy(H,Fock,dens,enuclear):
    from numpy.core.umath_tests import inner1d
    E1e = sum(inner1d(dens,H.T))
    E2e = sum(inner1d(dens,Fock.T))
    Etot = sum(inner1d(dens,H.T))+sum(inner1d(dens,Fock.T))+enuclear
    return Etot,2*E1e,E2e-E1e

    
