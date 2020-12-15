# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 16:04:03 2020

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
    ht = (b-a) / Nt
    return ht


def f(t,y):
    """
    Esta función calcula la derivada de la función con un paso antes en el tiempo

    Parameters
    ----------
    t : float
        paso del tiempo.
    y : float
        paso anterior en el timepo.

    Returns
    -------
    float
        derivada de la función.

    """
    return y - t **2 + 1


def Exacta (t):
    """
    Esta función calcula la solución exacta del problema

    Parameters
    ----------
    t : float
        tiempo.

    Returns
    -------
    float
        Solución exacta del problema.

    """
    return (t +1)**2 - 0.5 * np. exp (t)


def Euler (f, t, w, ht ):
    """
     Esta función calcula la solución numérica por el método forward euler

    Parameters
    ----------
    f : float
        función que calcula la derivada de la función.
    t : float
        tiempo.
    w : float
        paso anterior en la función.
    ht : float
        Tamaño del paso en el tiempo.

    Returns
    -------
    None.

    """
    for i, val in enumerate (w [0: -1]):
        w[i +1] = w[i] + ht * f(t[i], w[i])
        t[i +1] = t [0] + (i +1) * ht
        
        
        
def RK2 (f, t, w, ht ):
    """
    Esta función calcula la solución numérica por el método Runge-Kutta de orden dos

    Parameters
    ----------
    f : float
        función que calcula la derivada de la función.
    t : float
        tiempo.
    w : float
        paso anterior en la función.
    ht : float
        Tamaño del paso en el tiempo.

    Returns
    -------
    None.

    """
    for i, val in enumerate (w [0: -1]):
        k1 = ht * f(t[i], w[i])
        w[i +1] = w[i] + ht * f(t[i] + ht * 0.5 , w[i] + k1 * 0.5)
        t[i +1] = a + (i +1) * ht
        
        
def RK3 (f, t, w, ht ):
    """
    Esta función calcula la solución numérica por el método Runge-Kutta de orden tres

    Parameters
    ----------
     f : float
        función que calcula la derivada de la función.
    t : float
        tiempo.
    w : float
        paso anterior en la función.
    ht : float
        Tamaño del paso en el tiempo.

    Returns
    -------
    None.

    """
    for i, val in enumerate (w [0: -1]):
        k1 = ht * f(t[i], w[i])
        k2 = ht * f(t[i] + ht /3, w[i] + k1 / 3)
        k3 = ht * f(t[i] + 2 * ht / 3, w[i] + 2 * k2 / 3)
        w[i +1] = w[i] + (k1 + 3 * k3) / 4
        t[i +1] = a + (i +1) * ht
        
        
def RK4 (f, t, w, ht ):
    """
    Esta función calcula la solución numérica por el método Runge-Kutta de orden cuatro

    Parameters
    ----------
    f : float
        función que calcula la derivada de la función.
    t : float
        tiempo.
    w : float
        paso anterior en la función.
    ht : float
        Tamaño del paso en el tiempo.

    Returns
    -------
    None.

    """
    for i, val in enumerate (w [0: -1]):
        k1 = ht * f(t[i], w[i])
        k2 = ht * f(t[i] + ht /2, w[i] + k1 / 2)
        k3 = ht * f(t[i] + ht /2, w[i] + k2 / 2)
        k4 = ht * f(t[i] + ht , w[i] + k3) 
        w[i +1] = w[i] + (k1 + 2* k2 + 2* k3 + k4) / 6
        t[i +1] = a + (i +1) * ht
        

Nt =   4#8, 16, 32
a = 0
b = 4
ht = mesh (a, b, Nt)
y0 = 0.5

t = np.linspace (a, b, Nt +1)
tl = np. linspace (a, b, 100)
y_eul = np. zeros (Nt +1);
y_rk2 = np. zeros (Nt +1)
y_rk3 = np. zeros (Nt +1)
y_rk4 = np. zeros (Nt +1)

y_eul [0]= y0
y_rk2 [0]= y0
y_rk3 [0]= y0
y_rk4 [0]= y0

Euler (f, t, y_eul , ht)
RK2 (f, t, y_rk2 , ht)
RK3 (f, t, y_rk3 , ht)
RK4 (f, t, y_rk4 , ht)

yp = Exacta (tl)
yp_= Exacta (t)
e_eul = np. abs (yp_ - y_eul )
e_rk2 = np. abs (yp_ - y_rk2 )
e_rk3 = np. abs (yp_ - y_rk3 )
e_rk4 = np. abs (yp_ - y_rk4 )

n_error_eul = np. linalg . norm (e_eul , 2)
n_error_rk2 = np. linalg . norm (e_rk2 , 2)
n_error_rk3 = np. linalg . norm (e_rk3 , 2)
n_error_rk4 = np. linalg . norm (e_rk4 , 2)


Error='E Eul={:10.4f}, E RK2={:10.4f}, E RK3={:10.4f}, E RK4={:10.4f}'.format(n_error_eul,n_error_rk2,n_error_rk3,n_error_rk4)
plt.style.use(['Solarize_Light2'])
plt.suptitle (' Solución y Aproximación Nt = {:10.0f}'.format(Nt), fontsize =14)
plt.plot (tl , yp , 'm-', lw =3, label ='Sol.  Exacta ')
plt.plot (t, y_eul , 'C7o--', label ='Euler')
plt.plot (t, y_rk2 , 'C6^--', label ='RK2')
plt.plot (t, y_rk3, 'C8*--', label ='RK3')
plt.plot (t, y_rk4, 'C9*--', label ='RK4')
plt.title ( Error, fontsize =12 , color ='blue')
plt.xlim ( 0, 4)
plt.ylim ( -3 ,7)
plt.xlabel ('$t$ ')
plt.ylabel ('$y(t)$')
plt.legend (loc='upper left', ncol =1, framealpha =0.75, fancybox =True , fontsize =10)
plt.grid ( color ='w')
plt.subplots_adjust ( hspace =0.35)
plt.savefig('decaimiento_Nt_{}.pdf'.format(Nt))
plt.show ()



Error='E Eul={:10.4f}, E RK2={:10.4f}, E RK3={:10.4f}, E RK4={:10.4f}'.format(n_error_eul,n_error_rk2,n_error_rk3,n_error_rk4)
plt.style.use(['Solarize_Light2'])
plt.suptitle (' Errores Nt = {:10.0f} '.format(Nt), fontsize =14)
plt.plot (t, e_eul, 'C7o--', label ='Euler')
plt.plot (t, e_rk2, 'C6^--', label ='RK2')
plt.plot (t, e_rk3, 'C8*--', label ='RK3')
plt.plot (t, e_rk4, 'C9*--', label ='RK4')
plt.title ( Error, fontsize =12 , color ='blue')
plt.xlim ( 0, 4)
plt.ylim ( 0 ,1)
plt.xlabel ('$t$ ')
plt.ylabel ('$y(t)$')
plt.legend (loc='upper left', ncol =1, framealpha =0.75, fancybox =True , fontsize =10)
plt.grid ( color ='w')
plt.subplots_adjust ( hspace =0.35)
plt.show ()


