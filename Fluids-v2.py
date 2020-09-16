# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 23:08:51 2020

@author: Guillem
"""
import numpy as np
from matplotlib import pyplot as plt
from time import time
from numba import jit

@jit

def velocity(vel0,diff,p,D,T):
        
    """Esta función calcula la variacion del campo de velocidades en función del tiempo:
        
    Se resuelven las EDP's de Navier-Stokes, mediante la discretización de estas
    Las condiciones de contorno se mantienen siempre constantes:
        0 para la velocidad
        1 para la densidad
    """

    U, V = vel0
    dx, dy, dt = diff

    
    
    row, col = U.shape
    for t in T:
        for j in range(1, row-1):
            for i in range(1, col-1):
                
                dudx   = 1/dx*(U[i,j]-U[i-1,j])
                dudy   = 1/dy*(U[i,j]-U[i,j-1])
                dvdx   = 1/dx*(V[i,j]-V[i-1,j])
                dvdy   = 1/dy*(V[i,j]-V[i,j-1])
                d2udx2 = 1/dx**2*(U[i+1,j]-2*U[i,j]+U[i-1,j])
                d2udy2 = 1/dy**2*(U[i,j+1]-2*U[i,j]+U[i,j-1])
                d2vdx2 = 1/dx**2*(V[i+1,j]-2*V[i,j]+V[i-1,j])
                d2vdy2 = 1/dy**2*(V[i,j+1]-2*V[i,j]+V[i,j-1])
                
                dudt   = D*(d2udx2+d2udy2) - U[i,j]*(dudx+dudy)
                dvdt   = D*(d2vdx2+d2vdy2) - V[i,j]*(dvdx+dvdy)
                dpdt   = -p[i,j]*(dudx+dudy)
            
                U[i,j] += dt*dudt 
                V[i,j] += dt*dvdt 
                p[i,j] += dt*dpdt 
        
        U[0, :] = 0
        U[-1,:] = 0
        U[:, 0] = 0
        U[:,-1] = 0
        V[0, :] = 0
        V[-1,:] = 0
        V[:, 0] = 0
        V[:,-1] = 0
        p[0, :] = 1
        p[-1,:] = 1
        p[:, 0] = 1
        p[:,-1] = 1
    


def vorticity(vel0, diff, w):
    

    U, V = vel0
    dx, dy, dt = diff
    
    row, col = w.shape
    for j in range(1, row-1):        
        for i in range(1, col-1):    
            w[j,i] = (V[j,i] - V[j-1,i])/dx - (U[j,i] - U[j,i-1])/dy
    
    w[-1, :] = w[-2,:]		
    w[:, -1] = w[:,-2]


        
def plot2D(vel0,par,a):
   
    U, V = vel0
    x, y = par
    
    X, Y = np.meshgrid(x, y)
    fig = plt.figure(figsize=(11,7), dpi=100)
    plt.contourf(X, Y, a, alpha=0.5, cmap=plt.cm.winter)  
    plt.colorbar()
    plt.contour(X, Y, a, cmap=plt.cm.winter)  
    plt.quiver(X[::2, ::2], Y[::2, ::2], U[::2, ::2], V[::2, ::2]) 
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

def circle(U,V,resolution,par):
  
    x, y = par
    for j in range(0, resolution):
        for i in range(0, resolution):
            if (x[i]-3)**2 +(y[j]-3)**2 <= 4:
                U[j,i] = 2
                V[j,i] = 2

    return U, V

def main():

    resolution = 41    		                     # number of positions in axis
    Lx, Ly = 10, 10				                 # domain size along axis
    dx, dy = Lx/(resolution-1), Ly/(resolution-1)			      
    	      
    x, y  = np.linspace(0, Lx, resolution), np.linspace(0, Ly, resolution)	
    
    
    D  = 0.01
    dt = 0.01
    n  = 4500
    
    
    T = np.linspace(0,n*dt,n)
    
    
    (row,col) = (resolution,resolution)
    
    U0,V0 = np.zeros((row,col)), np.zeros((row,col))
    w = np.zeros((row,col))
    p = np.ones((row,col))
    
    par  = (x,y)
    vel0 = circle(U0,V0,resolution,par)
    diff = (dx,dy,dt)
    
    
    
    
                   
    velocity(vel0,diff,p,D,T)
    vorticity(vel0, diff, w)
    plot2D(vel0, par, w)

if __name__ == "__main__":
    main()