import numpy as np
import matplotlib.pyplot as plt
import sys


if(len(sys.argv)!=3):
    print("introduzca 3 parametros para correr el programa\n")
    exit()

lectura=str(sys.argv[1])
lectura2=str(sys.argv[2])


#existencia del archivo de entrada
try:
    archivo=open(lectura).read().split("\n")
except IOError:
   print "Error: el archivo "+lectura+ " no existe"
   exit()

#existencia del archivo de entrada
try:
    archivo2=open(lectura2).read().split("\n")
except IOError:
   print "Error: el archivo "+lectura2+ " no existe"
   exit()

#abrimos el primer archivo
x1=[]
P1=[]
rho1=[]
e1=[]
u1=[]

for i in range(len(archivo)-1):
    a= archivo[i].split(',')
    #print a
    x1.append(float(a[0]))
    P1.append(float(a[1]))
    rho1.append(float(a[2]))
    u1.append(float(a[3]))
    e1.append(float(a[4]))

#abrimos elsegundo archivo
x2=[]
P2=[]
rho2=[]
e2=[]
u2=[]

for i in range(len(archivo2)-1):
    a2= archivo2[i].split('\t')
    #print a2
    #print a
    x2.append(float(a2[0]))
    rho2.append(float(a2[1]))
    u2.append(float(a2[2]))
    e2.append(float(a2[3]))
    P2.append(float(a2[4]))

""" solucionando la funcion analitica
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
"""

X="$ x  $"
Y="$ y $"
title2="$ "+"Numerical\ solution" + " $"
title1="$ "+"Analitical\ solution" + " $"
title3="$ "+"Comparison" + " $"
title4="$ "+"Percentual\ error" + " $"


plt.figure(figsize=[15,10])
plt.subplot(221)
plt.ylim([-0.01,1.5]) 
plt.plot(x1,P1,label=title1)
plt.plot(x2,P2,label=title2)
plt.xlabel(X,size=20)
plt.ylabel("$ Pressure $",size=20)
plt.legend(loc=2)
  
plt.subplot(222)
plt.ylim([-0.01,0.3]) 
plt.plot(x1,u1,label=title1)
plt.plot(x2,u2,label=title2)
plt.xlabel(X,size=20)
plt.ylabel("$ U $",size=20)
plt.legend(loc=1)

plt.subplot(223)
plt.ylim([-0.01,1.1]) 
plt.plot(x1,rho1,label=title1)
plt.plot(x2,rho2,label=title2)
plt.xlabel(X,size=20)
plt.ylabel("$\rho$",size=20)
plt.legend(loc=3)

plt.subplot(224)
plt.ylim([-0.01,5]) 
plt.plot(x1,e1,label=title1)
plt.plot(x2,e2,label=title2)
plt.xlabel(X,size=20)
plt.ylabel('$ E $',size=20)
plt.legend(loc=2)

plt.show()
   
