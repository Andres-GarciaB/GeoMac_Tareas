import numpy as np

print()
print("+----------------------------------------------------+")
print("|      Solucion de la transferencia de calor         |")
print("+----------------------------------------------------+")
print("| Autores: Hernández Terán Oscar                     |")
print("|          Ángeles Yosselin                          |")
print("|          García Andrés                             |")
print("+----------------------------------------------------+")
print("|                 Datos de entrada                   |")
print("+----------------------------------------------------+")


a = float(input("Ingrese el comienzo de la barra.                a="))
b = float(input("Ingrese el fin de la barra.                     b="))
K = float(input("Ingresa la conductividad térmica del material   k="))
N = int(input("Ingresa el número de nodos que desea            N="))
Ta = float(input("Ingrese la temperaruta al inicio.               Ta="))
Tb = float(input("Ingrese la temperaruta al final.                Tb="))

h = (b-a)/(N+1)
r = K/(h**2)

def Arreglos(N,Ta,Tb):
    A = np.zeros((N,N))
    T = np.zeros(N)    
    b = np.zeros(N)
    b[0] = Ta
    b[-1] = Tb
    return A, T, b

A, T, b = Arreglos(N,Ta,Tb)




A = np.zeros((N, N))
diagonal=-2
A[0,0] = diagonal; A[0,1] = 1
for i in range(1,N-1):
     A[i,i] = diagonal
     A[i,i+1] = 1
     A[i,i-1] = 1
A[N-1,N-2] = 1; A[N-1,N-1] = diagonal

print("\n Matriz del sistema : \n", A)
print("\n Lado derecho del sistema : \n", b)