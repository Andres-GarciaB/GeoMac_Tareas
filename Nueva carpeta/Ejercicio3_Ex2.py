# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 13:51:22 2020

@author: andy_
"""

import numpy as np
import matplotlib.pyplot as plt

def mesh (a, b, Nt ):
    """
    Esta función genera la malla del problema

    Parameters
    ----------
    a : float
        Valor al inicio del dominio.
    b : float
        valor al final del dominio.
    Nt : float
        Número de pasos en el tiempo.

    Returns
    -------
    ht : float
        Tamaño del paso en el tiempo.

    """
    ht = (b-a)/Nt
    return ht

def exactSolution (t, y0 , lam ):
    """
    Esta función calcula la solución exacta del problema

    Parameters
    ----------
    t : float
        tiempo.
    y0 : float
        Condición inicial. Cantidad de sustancia inicial.
    lam : float
        Constante de decaimiento.

    Returns
    -------
    float
        Cantidad de sustancia en el tiempo t.

    """
    return y0 / (y0 + ((1 - y0 ) * np.exp(- lam * t)))

def forwardEuler (y, ht , lam ):
    """
    Esta función calcula la solución numérica por el método forward euler

    Parameters
    ----------
    y : float
        Vector con las condiciones a la frontera.
    ht : float
        Tamaño del paso en el tiempo.
    lam : float
        Constante de decaimiento.

    Returns
    -------
    An : float
        vector que contiene la solución con el método Forward Euler.

    """
    A=1+ht*lam
    An = [A]
    for i, val in enumerate (y [0: -1]):
        y[i +1] = ( A * y[i] ) - ( ht*lam*(y[i]**2))
        An. append (An[i] * A)
    return An


Nt_4 = 4
Tmax = 1
ht = mesh (0, Tmax , Nt_4)
y0 = 0.01
lam = 10
t = np. linspace (0, Tmax , Nt_4 +1)
yf = np. zeros (Nt_4 +1)
yf [0] = y0
An = forwardEuler (yf , ht , lam )
tl = np. linspace (0, Tmax , 100)
y_exacta = exactSolution (tl , y0 , lam )
y_exac_p = exactSolution (t, y0 , lam)
norma_error_f = np.linalg.norm(yf - y_exac_p ,2)


Nt_16 = 16
ht_16 = mesh (0, Tmax , Nt_16)
t_16 = np. linspace (0, Tmax , Nt_16 +1)
yf_16 = np. zeros (Nt_16 +1)
yf_16 [0] = y0
An_16 = forwardEuler (yf_16 , ht_16 , lam )
t_16 = np. linspace (0, Tmax , Nt_16 +1)
y_exac_16 = exactSolution (t_16, y0 , lam)
norma_error_f_16 = np.linalg.norm(yf_16 - y_exac_16 ,2)


Nt_64 = 64
ht_64 = mesh (0, Tmax , Nt_64)
t_64 = np. linspace (0, Tmax , Nt_64 +1)
yf_64 = np. zeros (Nt_64 +1)
yf_64 [0] = y0
An_64 = forwardEuler (yf_64 , ht_64 , lam )
t_64 = np. linspace (0, Tmax , Nt_64 +1)
y_exac_64 = exactSolution (t_64, y0 , lam)
norma_error_f_64 = np.linalg.norm(yf_64 - y_exac_64 ,2)


Ecuacion = ' $y(t) = y_0 / y_0 + (1-y_0) (e^{-\lambda t})$,'
plt.style.use(['Solarize_Light2'])
plt.suptitle (' Equación Logística ', fontsize =14)
plt.plot (tl , y_exacta , 'm-', lw =3, label ='Sol.  Exacta ')
plt.plot (t, yf , 'C7o--', label ='FE. Nt = 4, Error = {:10.5f}'.format(norma_error_f ))
plt.plot (t_16, yf_16 , 'C6o--', label ='FE. Nt = 16, Error = {:10.5f}'.format(norma_error_f_16))
plt.plot (t_64, yf_64 , 'C8*-', label ='FE. Nt = 64, Error = {:10.5f}'.format(norma_error_f_64 ))
plt.title ( Ecuacion, fontsize =12 , color ='blue')
plt.xlim ( 0, 1)
plt.ylim ( 0 ,1)
plt.xlabel ('$t$ ')
plt.ylabel ('$y(t)$')
plt.legend (loc='upper left', ncol =1, framealpha =0.75 , fancybox =True , fontsize =10)
plt.grid ( color ='w')

plt.subplots_adjust ( hspace =0.35)
plt.savefig('decaimiento_Nt_{}.pdf'.format(Nt_4))
plt.show ()
