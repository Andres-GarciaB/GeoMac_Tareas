import  numpy  as  np

print ()
print ( "+ ---------------------------------------------- ------ + " )
print ( "| Solucion de la transferencia de calor |" )
print ( "+ ---------------------------------------------- ------ + " )
print ( "| Autores: Hernández Terán Oscar |" )
print ( "| Ángeles Yosselin |" )
print ( "| García Andrés |" )
print ( "+ ---------------------------------------------- ------ + " )
print ( "| Datos de entrada |" )
print ( "+ ---------------------------------------------- ------ + " )


a  =  float ( input ( "Ingrese el comienzo de la barra. a =" ))
b  =  float ( input ( "Ingrese el fin de la barra. b =" ))
K  =  float ( input ( "Ingresa la conductividad térmica del material k =" ))
Nodos  =  int ( input ( "Ingresa el número de nodos que desea N =" ))
Ta  =  float ( input ( "Ingrese la temperaruta al inicio. Ta =" ))
Tb  =  float ( input ( "Ingrese la temperaruta al final. Tb =" ))

h  = ( b - a ) / Nodos
r  =  K / ( h ** 2 )

def  Arreglos ( Nodos , Ta , Tb ):
    A  =  np . zeros (( Nodos , Nodos ))
    T  =  np . zeros ( Nodos )    
    b  =  np . zeros ( Nodos )
    b [ 0 ] =  Ta
    b [ - 1 ] =  Tb
    return (A , T , b)

A , T , b  =  Arreglos ( Nodos , Ta , Tb )
print ( A )
print ( T )
print ( b )
