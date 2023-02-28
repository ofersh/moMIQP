# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 16:20:45 2022
@author: ofers
"""

import numpy as np
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.variable import Real, Integer #, Choice, Binary
from scipy.linalg import hadamard

N = 32;
c1=np.array([7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7])
c2=np.array([-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4])
c3=np.array([7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7,7,-7])
c4=np.array([-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4,-4,4])
LOWER=-10
UPPER=10
class MixedVarsSpheres(ElementwiseProblem):

    def __init__(self, **kwargs):

        variables = dict()
        
        for k in range(1, N+1):
            variables[f"x{k:02}"] = Real(bounds=(1.0*LOWER, 1.0*UPPER))

        for k in range(N+1, 2*N+1):
            variables[f"x{k:02}"] = Integer(bounds=(LOWER, UPPER))

        super().__init__(vars=variables, n_obj=2, **kwargs)

    def _evaluate(self, y, out, *args, **kwargs):
        x = np.array([y[f"x{k:02}"] for k in range(1, N+1)])
        z = np.array([y[f"x{k:02}"] for k in range(N+1, 2*N+1)])
        f1 = np.array(x-c1).dot(np.array(x-c1)) + np.array(z-c1).dot(np.array(z-c1))
        f2 = np.array(x-c2).dot(np.array(x-c2)) + np.array(z-c2).dot(np.array(z-c2))
        out["F"] = [f1, f2]
#
class MixedVarsEllipsoid(ElementwiseProblem):

    def __init__(self, H, c, **kwargs):
        self.H = np.copy(H)
        self.c = c
        variables = dict()
        for k in range(1, N+1):
            variables[f"x{k:02}"] = Real(bounds=(1.0*LOWER, 1.0*UPPER))
        for k in range(N+1, 2*N+1):
            variables[f"x{k:02}"] = Integer(bounds=(LOWER, UPPER))
        super().__init__(vars=variables, n_obj=2, **kwargs)

    def _evaluate(self, y, out, *args, **kwargs):
        x = np.array([y[f"x{k:02}"] for k in range(1, N+1)])
        z = np.array([y[f"x{k:02}"] for k in range(N+1, 2*N+1)])

        f1 = (np.array(x-c1).dot(self.H).dot(np.array(x-c1)) + np.array(z-c1).dot(self.H).dot(np.array(z-c1)))/self.c
        f2 = (np.array(x-c2).dot(self.H).dot(np.array(x-c2)) + np.array(z-c2).dot(self.H).dot(np.array(z-c2)))/self.c
        out["F"] = [f1, f2]
class MixedVarsRotEllipsoid(ElementwiseProblem):

    def __init__(self, H, c, **kwargs):
        self.H = np.copy(H)
        self.c = c
        variables = dict()
        for k in range(1, N+1):
            variables[f"x{k:02}"] = Real(bounds=(1.0*LOWER, 1.0*UPPER))
        for k in range(N+1, 2*N+1):
            variables[f"x{k:02}"] = Integer(bounds=(LOWER, UPPER))
        super().__init__(vars=variables, n_obj=2, **kwargs)

    def _evaluate(self, y, out, *args, **kwargs):
        x = np.array([y[f"x{k:02}"] for k in range(1, 2*N+1)])
        
        f1 = (np.array(x-c3).dot(self.H).dot(np.array(x-c3)))/self.c
        f2 = (np.array(x-c4).dot(self.H).dot(np.array(x-c4)))/self.c
        out["F"] = [f1, f2]

def genHsphere(N,c=1) :
    return c*np.eye(N)

def genHellipse(N,c) :
	H = np.zeros((N,N))
	for i in range(N) :
		H[i,i] = c**(i/(N-1))
	return H

def genHellipse2(N,c) :
	H = np.zeros((N,N))
	for i in range(N) :
		H[i,i] = 1+(i*(c-1)/(N-1))
	return H
    
def genHadamardHellipse (N,c) :
    H = genHellipse2(N,c)
    R = hadamard(N) / np.sqrt(N)
    H = R.dot(H).dot(np.transpose(R))
    return H

def getRotation(N,theta) :
    v = np.ones((N,1))
    u = np.ones((N,1))
    for i in range(N) :
        if np.mod(i,2)==0 :
            u[i] = 0
        else :
            v[i] = 0
    v=v/np.linalg.norm(v)
    u=u/np.linalg.norm(u)
    R = np.eye(N) + np.sin(theta)*(u.dot(np.transpose(v)) - v.dot(np.transpose(u))) + (np.cos(theta)-1.0)*(v.dot(np.transpose(v)) + u.dot(np.transpose(u)))
    return R
    
def genRotatedHellipse (N,c,theta=.333*np.pi) :
    H = genHellipse(N,c)
    R = getRotation(N,theta)
    H = R.dot(H).dot(np.transpose(R))
    return H

def genHcigar(N,c) :
    H = c*np.eye(N)
    H[-1][-1] = 1.0
    return H

def genHdiscus(N,c) :
    H = np.eye(N)
    H[-1][-1] = c
    return H    