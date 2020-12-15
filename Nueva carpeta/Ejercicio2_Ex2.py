# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 00:43:10 2020

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
        Valor al final del dominio.
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
        Tiempo.
    y0 : float
        Condición inicial. Cantidad de sustancia inicial.
    lam : float
        Constante de decaimiento.

    Returns
    -------
    float
        Cantidad de sustancia en el tiempo t.

    """
    return y0 * np. exp (- lam * t)

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
    A=1-ht*lam
    An = [A]
    for i, val in enumerate (y [0: -1]):
        y[i +1] = A * y[i]
        An. append (An[i] * A)
    return An

def backwardEuler (y, ht , lam ):
    """
    Esta función calcula la solución numérica por el método Backward Euler

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
    Bn : float
        vector que contiene la solución con el método Backward Euler.

    """
    B = 1 /(1 + ht*lam)
    Bn = [B]
    for i, val in enumerate (y [0: -1]):
        y[i +1] = B * y[i]
        Bn. append (Bn[i] * B)
    return Bn


Nt = 7
Tmax = 10
ht = mesh (0, Tmax , Nt)
y0 = 20
lam = 1.5
t = np. linspace (0, Tmax , Nt +1)
yf = np. zeros (Nt +1)
yb = np. zeros (Nt +1)
yf [0] = y0
yb [0] = y0
An = forwardEuler (yf , ht , lam )
Bn = backwardEuler (yb , ht , lam )
tl = np. linspace (0, Tmax , 100)
y_exacta = exactSolution (tl , y0 , lam )
y_exac_p = exactSolution (t, y0 , lam)
norma_error_f = np.linalg.norm(yf - y_exac_p ,2)
norma_error_b = np.linalg.norm(yb - y_exac_p ,2)

Ecuacion = ' $y(t) = y_0 e^{\lambda t}$,' + '  $N_t$ ' + '= {} '. format (Nt) + ',  $h_t$ ' + '=  {:03.2f}'. format (ht)
Error = ',  Error : FE =  {:10.9f}, BE =  {:10.9f}'. format ( norma_error_f , norma_error_b )

plt.style.use(['Solarize_Light2'])
fig,(ax1 , ax2 ) = plt.subplots(2 ,1)
fig.suptitle (' Decaimiento   Radioactivo ', fontsize =14)
ax1.plot (tl , y_exacta , 'g-', lw =3, label ='Sol.  Exacta ')
ax1.plot (t, yf , 'C7o--', label ='Forward   Euler ')
ax1.plot (t, yb , 'C6o--', label ='Backward   Euler ')
ax1.set_title ( Ecuacion + Error , fontsize =12 , color ='blue')
ax1.set_xlim ( -0.5 ,t [ -1]+0.5)
ax1.set_ylim ( -30 ,30)
ax1.set_xlabel ('$t$ ')
ax1.set_ylabel ('$y(t)$')
ax1.legend (loc='upper left', ncol =1, framealpha =0.75 , fancybox =True , fontsize =10)
ax1.grid ( color ='w')
nticks = np.arange(1, Nt +1 ,1)
ax2.plot( nticks , An[: -1] , 'C7v-', label ='$A^n$ ')
ax2.plot( nticks , Bn [: -1] , 'C6^-', label ='$B^n$ ')
ax2.set_xlim( -0.5 , Nt +0.5)
ax2.set_xticks ( nticks )
ax2.set_xlabel ('$n$ ')
ax2.set_ylabel ('$A^n$ , $B^n$ ')
ax2.legend (loc='upper right', ncol =1, framealpha =0.75 , fancybox =True , fontsize =10)
ax2.grid ( color ='w')
plt.subplots_adjust ( hspace =0.35)
plt.savefig('decaimiento_Nt_{}.pdf'.format(Nt))
plt.show ()
