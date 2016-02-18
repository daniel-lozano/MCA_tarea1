import numpy as np
import matplotlib.pyplot as plt
import sys


lectura=str(sys.argv[1])
c1=int(sys.argv[2])
c2=int(sys.argv[3])
c3=int(sys.argv[4])
c4=int(sys.argv[5])
archivo=open(lectura).read().split("\n")

x1=[]
y1=[]
x2=[]
y2=[]

for i in range(len(archivo)-1):
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


plt.figure(figsize=[10,5])
plt.subplot(221)
plt.plot(x1,y1,".")
plt.xlabel(X,size=20)
plt.ylabel(Y,size=20)
plt.title(title1)
    
plt.subplot(222)
plt.plot(x2,y2,".")
plt.xlabel(X,size=20)
plt.ylabel(Y,size=20)
plt.title(title2)

plt.show()
   
