# -*- coding: utf-8 -*-
"""
@author: Guillem
"""

import numpy as np
from matplotlib import pyplot as plt



class Fluids:
    """Esta clase crea un objeto dada una distribución inicial:
        
    Crea un objeto con comportamiento de fluido
    Permite calcular y representar como varia la velocidad y vorticidad de este
    Representa en dos dimensiones estos resultados
    Los calculos se realizan en función del tiempo
    """
    
    
    def __init__(self,vel0,diff,par,p,w,T,D):
        
        self.U  = vel0[0]
        self.V  = vel0[1]
        self.dx = diff[0]
        self.dy = diff[1]
        self.dt = diff[2]
        self.x  = par[0]
        self.y  = par[1]
        self.p  = p
        self.w  = w
        self.T  = T
        self.D  = D


    def velocity(self):
        
        """Esta función calcula la variacion del campo de velocidades en función del tiempo:
            
        Se resuelven las EDP's de Navier-Stokes, mediante la discretización de estas
        Las condiciones de contorno se mantienen siempre contantes:
                0 para la velocidad
                1 para la densidad
        """
        
        U  = self.U
        V  = self.V
        p  = self.p
        dx = self.dx
        dy = self.dy
        dt = self.dt
        D  = self.D
        T  = self.T
        
        
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
                
                    U[i,j] = dt*dudt + U[i,j]
                    V[i,j] = dt*dvdt + V[i,j]
                    p[i,j] = dt*dpdt + p[i,j]
            
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
          
                
    def vorticity(self):
        
        U  = self.U
        V  = self.V
        w  = self.w
        dx = self.dx
        dy = self.dy

        
        
        row, col = w.shape
        for j in range(1, row-1):        
            for i in range(1, col-1):    
                w[j,i] = (V[j,i] - V[j-1,i])/dx - (U[j,i] - U[j,i-1])/dy
        
        w[-1, :] = w[-2,:]		
        w[:, -1] = w[:,-2]
      
        
    def plot2D(self,a):
        
        U  = self.U
        V  = self.V

        
        X, Y = np.meshgrid(x, y)
        fig = plt.figure(figsize=(11,7), dpi=100)
        plt.contourf(X, Y, a, alpha=0.5, cmap=plt.cm.winter)  
        plt.colorbar()
        plt.contour(X, Y, a, cmap=plt.cm.winter)  
        plt.quiver(X[::2, ::2], Y[::2, ::2], U[::2, ::2], V[::2, ::2]) 
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.show()



def circle(U,V):
    for j in range(0, resolution):
        for i in range(0, resolution):
            if (x[i]-3)**2 +(y[j]-3)**2 <= 4:
                U[j,i] = 2
                V[j,i] = 2
    return U, V



resolution = 41    		                     # number of positions in axis
Lx, Ly = 10, 10				                 # domain size along axis
dx, dy = Lx/(resolution-1), Ly/(resolution-1)			      
	      
x, y  = np.linspace(0, Lx, resolution), np.linspace(0, Ly, resolution)	


D  = 0.01
dt = 0.01
n  = 300


T = np.linspace(0,n*dt,n)


(row,col) = (resolution,resolution)

U0,V0 = np.zeros((row,col)), np.zeros((row,col))
w = np.zeros((row,col))
p = np.ones((row,col))


vel0 = circle(U0,V0)
diff = (dx,dy,dt)
par  = (x,y)


fluid = Fluids(vel0,diff,par,p,w,T,D)
               
fluid.velocity()
fluid.vorticity()
fluid.plot2D(w)
