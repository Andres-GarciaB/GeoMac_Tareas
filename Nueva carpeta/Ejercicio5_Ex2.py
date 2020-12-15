# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 14:42:23 2020

@author: andy_
"""
import numpy as np
import matplotlib.pyplot as plt

L = 1 # Longitud del dominio
N = 20 # Numero de incognitas internas
Tmax = 1.0 # Tiempo maximo de simulacion
ht = 0.05 # Paso de tiempo
alpha = 2 # Dato fisico
h = L / (N +1) # Tamanio de la malla espacial
Nt = int ( Tmax / ht) # Numero total de pasos
lamb = alpha * ht / h # Parametro lambda
Tmax = Nt * ht # Tiempo total de simulacion

def f(x):
    """
    Función que calcula la condición inicial

    Parameters
    ----------
    x : float
        Coordenadas de la malla.

    Returns
    -------
    float
        Condición inicial.

    """
    return np. sin (np.pi * x)


def g(x):
    """
    Función que calcula la velocidad inicial de la onda

    Parameters
    ----------
    x : float
        Coordenadas de la malla.

    Returns
    -------
    float
        velocidad inicial de la onda.

    """
    return 0


def solExacta (x, t):
    """
    Esta función calcula la solución exacta

    Parameters
    ----------
    x : float
        coordenadas de la malla.
    t : float
        Tiempo.

    Returns
    -------
    float
        Solución de la ecuación de onda.

    """
    return np. sin (np.pi * x) * np.cos (2 * np.pi * t)


def calcError (sol_n , sol_e ):
    """
    Función que calcula el error

    Parameters
    ----------
    sol_n : float
        Solución numérica.
    sol_e : float
        solución exacta.

    Returns
    -------
    float
        El valor absoluto de la diferencia entre soluciones.

    """
    return np. abs (sol_n - sol_e )

def condicionesIniciales (l, ht , u, x, op =1):
    """
    Esta función genera las condiciones iniciales

    Parameters
    ----------
    l : float
        Parámetro auxiliar lambda.
    ht : float
        Paso de tiempo.
    u : float
        Condición inicial.
    x : float
        Coordenadas de la malla.
    op : float, optional
        Tipo de solución numérica. The default is 1.

    Returns
    -------
    w : float
        Solución numérica.

    """
    N = len (u)
    w = np. zeros (N)
    for i in range (1,N -1):
         if op == 1:
              w[i] = u[i] + ht * g(x[i])
         else :
             w[i] = (1 - l **2) * u[i] + 0.5 * l **2 * (u[i +1] + u[i -1]) + ht * g(x[i])
    return w


x = np. linspace (0,L,N +2) # Coordenadas de la malla
u = f(x) # Condicion inicial
w = condicionesIniciales (lamb , ht , u, x, op =2) # Euler :op = 1
plt . suptitle ('Ecuacion  de  onda ', fontsize =14)
plt . plot (x, u,'ro--', lw = 1, label="$u_{i ,0}$")
plt . plot (x, w,'kx--', lw = 1, label = " $u_ {i ,1} $")
plt . title (' Condiciones   iniciales  $\mathcal{O}( h_t )$', color ='blue', fontsize =12)
plt . ylabel ('$u(x,t)$')
plt . xlabel ('$x$ ')
plt.legend (loc='upper left', ncol =1, framealpha =0.75, fancybox =True , fontsize =12)
plt.grid ( color ='w')
plt.subplots_adjust ( hspace =0.35)
plt . savefig ('condicion_O3.pdf')
plt.show ()

def solver (u, w, N, x, Nt , l):
    s = np.zeros(N+2)
    for n in range (1, Nt ):
        for i in range (1,N+1 ):
            s[i] = 2 * (1 - l **2) * w[i] + l**2 * (w[i +1] + w[i -1]) - u[i]
        u = w.copy ()
        w = s.copy ()
        
        plt . plot (x,s,'--')
       
    return s
    #return (s/1000000)

t = np. linspace (0,Tmax,N +2) # Coordenadas de la malla
w1 = condicionesIniciales (lamb , ht , u, x, op = 1) # Euler :op = 1
w2 = condicionesIniciales (lamb , ht , u, x, op =2) # Euler :op = 1
yE= solExacta(x,t)
s = solver (u, w2, N, x, Nt , lamb )


plt . suptitle ('Ecuacion  de  onda ', fontsize =14)
plt . plot (x, u,'ro--', lw = 1, label="$u_{i ,0}$")
plt . plot (x, w,'kx--', lw = 1, label = " $u_ {i ,1} $")
plt . plot (x, s,'C7o--', lw = 1, label = " Numérica")
plt . plot (x, yE,'C6^--', lw = 1, label = " Exacta")
plt . title (' Condiciones   iniciales  $\mathcal{O}(h^3_t )$', color ='blue', fontsize =12)
plt . ylabel ('$u(x,t)$')
plt . xlabel ('$x$ ')
plt.legend (loc='upper left', ncol =1, framealpha =0.5, fancybox =True , fontsize =8)
plt.grid ( color ='w')
plt.subplots_adjust ( hspace =0.35)
plt . savefig ('condicion_O3.pdf')
plt.show ()