import numpy as np
import matplotlib.pyplot as plt
import sys


if(len(sys.argv)!=8):
    print "introduzca 6 parametros para correr el programa\n"
    exit()
    #print "parametros dados="+ len(sys.argv)
   

lectura=str(sys.argv[1])
c1=int(sys.argv[2])
c2=int(sys.argv[3])
c3=int(sys.argv[4])
c4=int(sys.argv[5])

#existencia del archivo de entrada
try:
    archivo=open(lectura).read().split("\n")
except IOError:
    
   print "Error: el archivo "+lectura+ " no existe"
   exit()

x1=[]
y1=[]
x2=[]
y2=[]

for i in range(int(len(archivo)-1)):
    a= archivo[i].split()
    #print a
    x1.append(a[c1])
    y1.append(a[c2])
    x2.append(a[c3])
    y2.append(a[c4])
    

X="$ "+ str(sys.argv[-2]) + " $"
Y="$ "+ str(sys.argv[-1]) + " $"
title1="$ "+"Rungekutta" + " $"
title2="$ "+"Simplectic" + " $"


plt.figure(figsize=[20,8])

plt.subplot(121)
plt.scatter(x1,y1,s=0.01)
plt.xlim([-2,2])
plt.ylim([-2,2])
plt.xlabel(X,size=20)
plt.ylabel(Y,size=20)
plt.title(title1)

    
plt.subplot(122)
plt.scatter(x2,y2,s=0.01)
plt.xlim([-2,2])
plt.ylim([-2,2])
plt.xlabel(X,size=20)
plt.ylabel(Y,size=20)
plt.title(title2)

plt.show()
   
