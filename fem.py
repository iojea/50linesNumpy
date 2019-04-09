#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 13:27:39 2019

Finite Element Library

@author: io
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from math import factorial

np.set_printoptions(threshold=np.inf)

# Rigidity matrix. 
    # Parameters: triang - Delaunay triangulation
    #             file_name - file for saving
    # Return: R - rigidity matrix.
def rigidity_matrix(P,T,file_name=""):
    N = P.shape[0]       #Cantidad de nodos
    d = P.shape[1]       #Dimensión
    unos = np.ones((1,d+1))
    rigi = np.vstack((np.zeros((1,d)),np.eye(d)))
    R = lil_matrix((N,N))
    for t in T:
         nodos = P[t,:]
         H = np.vstack((unos,nodos.T))
         G = np.linalg.solve(H,rigi)
         detH = np.absolute(np.linalg.det(H))
         R_loc = detH*G@G.T/factorial(d)
         R[t[:,None],t] = R[t[:,None],t] + R_loc
    if len(file_name)>0:
        archivo=open(file_name+".csv",'w')
        print(R,file=archivo)
        archivo.close()
    return R
         
# Solver for Poisson equation. 
    # Paramenters: triang -  Delaunay triangulation
    #              interior - list of interior nodes
def solver_poisson(P,T,interior):
    R = rigidity_matrix(P,T)
    N = P.shape[0]
    b = np.zeros((N,1))
    b[0] = 1
    R=R[interior[:,None],interior].tocsr()
    u = np.zeros([len(P),1])
    sol = spsolve(R,b[interior])
    u[interior] = sol.reshape(len(interior),1)
    return u
    

# Grafica soluciones 2d
def plot_solution(P,T,u):
    fig = plt.figure()
    ax =fig.gca(projection='3d')
    ax.plot_trisurf(P[:,0],P[:,1],T,u[:,0])
    plt.show()

