import numpy as np
import matplotlib.pyplot as plt
import sys


if(len(sys.argv)!=4):
    print("introduzca 3 parametros para correr el programa\n")
    exit()

lectura=str(sys.argv[1])
c1=int(sys.argv[2])
c2=int(sys.argv[3])

#existencia del archivo de entrada
try:
    archivo=open(lectura).read().split("\n")
except IOError:
   print "Error: el archivo "+lectura+ " no existe"
   exit()

x1=[]
y1=[]


for i in range(len(archivo)-1):
    a= archivo[i].split()
    #print a
    x1.append(float(a[c1]))
    y1.append(float(a[c2]))

#solucionando la funcion analitica
inicio=1
final=6
x2=np.linspace(inicio,final,len(x1))
y2=np.zeros(len(x2))

for i in range(len(x2)):
    y2[i]=x2[i]**2
div=[]
resta=[]
for i in range(len(x2)):
    div.append(abs((y1[i]-y2[i])/y1[i]))
    resta.append(y1[i]-y2[i])


X="$ x  $"
Y="$ y $"
title1="$ "+"Numerical\ solution" + " $"
title2="$ "+"Analitical\ solution" + " $"
title3="$ "+"Comparison" + " $"
title4="$ "+"Percentual\ error" + " $"


plt.figure(figsize=[10,5])
plt.subplot(221)
plt.plot(x1,y1,label=title1)
plt.xlabel(X,size=20)
plt.ylabel(Y,size=20)
plt.legend()
  
plt.subplot(222)
plt.plot(x2,y2,label=title2)
plt.xlabel(X,size=20)
plt.ylabel(Y,size=20)
plt.legend()

plt.subplot(223)
plt.plot(x1,resta,label=title3)
plt.xlabel(X,size=20)
plt.ylabel(Y,size=20)
plt.legend()

plt.subplot(224)
plt.plot(x1,div,label=title4)
plt.xlabel(X,size=20)
plt.ylabel(Y,size=20)
plt.legend()

plt.show()
   
