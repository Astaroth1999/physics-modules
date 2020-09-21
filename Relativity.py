# -*- coding: utf-8 -*-

"""
Created on Fri May  1 23:01:34 2020

@author: Guillem
"""


import numpy as np
import matplotlib.pyplot as plt       
from matplotlib import rcParams                                     
from matplotlib import animation      
from scipy.optimize import brentq 
from numba import jit    

plt.ioff()


    
"""Animación de una barra que gira con consideraciones relativistas

Esta clase toma como parametros:
    
la longitud de la barra (L)
la velocidad de la luz (c)
la velocidad de la barra (V)
la velocidad angular de esta (Omega)
el objeto a actualizar (line, en este caso la barra)
"""
   



def animate(i,L,c,V,Omega,line):
    
    """Aquí generamos la animación de la barra
    
    Definimos las variables necesarias en función de los parametros
    Definimos una función f(t) según las transformaciones de Lorentz tal que f(t) = 0 
    Mediante brentq resolvemos esta ecuación (f(t) = 0) debido a que esta es no lineal
    Añadimos los resultados al objeto line, y lo devolvemos
    """

    gamma = 1/np.sqrt(1-(V/c)**2)           
    
    
    tp = i*0.05              # Tiempo propio para cada itreración
    r = 0.01*L                 # Inicializamos r para subdividir la barra
    x, y = [], []              # Inicializamos las listas de valores X e Y
   
    @jit
    def f(t): 
        
        """Definimos f(t',t) = 0
        
        Obtenida de aplicar las transformaciones de de Lorentz a X propia, Y propia y T propia (tp)
        """
        
        return t - V/(c**2)*r*np.cos(Omega*t)-tp/gamma
    
    
    # =========================================================================
    # A continuación actualizamos el valor de line, para cada iteración de la animación.
    # =========================================================================
    
    while r <= 1.01*L:
        
          """
          NOTA: Aunque siendo ti = brentq(f,a,b). Si b es = 0.19 ya se 
          obtiene una diferencia de signo entre f(a) y f(b), pero por 
          razones que desconozco, cuando menor es el valor de b, menor 
          es el tiempo de animación, llegando a ser inexistente para b < 1
          (puede que por que al actualizar el tiempo propio haga falta 
          aumentar el intervalo para obtener la condición antes 
          especificada). Para evitar problemas se ha decidido dejar b
          como 1000, ya que visto como evoluciona la animación para valores
          de b >> 1, es coherente pensar que no existe otro cambio de 
          signo que pueda llevar a desencadenar un error en el cálculo  
          siguiente
          """
          
          
          ti = brentq(f,0,100)                   
          
                                                  
          auxx = gamma*(r*np.cos(Omega*ti)-V*ti)  # Coordenada X en S' (ti) 
          auxy = r*np.sin(Omega                            # obtenidos a las listas
                                                  # de posiciones antes
                                                  # inicializadas
          *ti)               # Coordenada Y en S' (ti) 
          
          x.append(auxx), y.append(auxy)          # Añadimos los valores 
                      
          r = r + 0.01*L                          # Actualizamos r (hacerlo 
                                                  # con un valor 0.01 
                                                  # es equivalente a 
                                                  # subdividir L' en 100 
                                                  # trozos) 
                                                  
    line.set_data(x, y)                           # Añadimos a line, las 
                                                  # coordenadas de cada 
                                                  # subdivision de la barra
                                                  # para el instante ti en
                                                  # S', tp en S                            

    return line,                                  # Devolvemos line,   






# =============================================================================
# Modificamos la fuente por defecto de las etiquetas para un eventual uso de
# LaTeX en estas
# =============================================================================

rcParams['mathtext.fontset'] = 'cm'
rcParams['mathtext.rm'] = 'serif'





fig, ax = plt.subplots(1,1,facecolor='darkslategray',constrained_layout = True) 

ax.set_xlim(-17,3), ax.set_ylim(-10,10)    
ax.set_xlabel('$X$',fontsize=10,rotation=0, color='lightgreen')
ax.set_ylabel('$Y$',fontsize=10,rotation=0, color='lightgreen')
ax.tick_params(labelcolor='gold'), ax.set_facecolor('whitesmoke')
ax.grid(True, color='yellowgreen')

line, = ax.plot([], [], lw=2)



L = 2    
c = 1                              
V, Omega = 0.4*c, 0.9*c/L    



anim = animation.FuncAnimation(fig, animate, fargs=[L,c,V,Omega,line], frames=1000, interval=17, blit=True)

anim.save('prueba.mp4', writer = 'ffmpeg', fps = 120)
