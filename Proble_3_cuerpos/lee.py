import numpy as np
import matplotlib.pyplot as plt
import sys


if(len(sys.argv)!=9):
    print "introduzca 6 parametros para correr el programa\n"
    exit()
    #print "parametros dados="+ len(sys.argv)
   

lectura=str(sys.argv[1])
c1=int(sys.argv[2])
c2=int(sys.argv[3])
c3=int(sys.argv[4])
c4=int(sys.argv[5])
Zoom=str(sys.argv[-1])
print "c1=",c1, "c2=", c2, "c3=", c3, "c4=", c4
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
    

if(Zoom=="si"):
    print "con zoom\n"
    X="$ "+ str(sys.argv[-3]) + " $"
    Y="$ "+ str(sys.argv[-2]) + " $"
    title1="$ "+"Rungekutta" + " $"
    title2="$ "+"Simplectic" + " $"
    
    
    plt.figure(figsize=[20,8])

    plt.subplot(221)
    plt.scatter(x1,y1,s=0.01)
    #plt.xlim([-3,3])
    plt.ylim([-2.5,2.5])
    plt.xlabel(X,size=20)
    plt.ylabel(Y,size=20)
    plt.title(title1)
    
    
    plt.subplot(222)
    plt.scatter(x2,y2,s=0.01)
    #plt.xlim([-3,3])
    plt.ylim([-2.5,2.5])
    plt.xlabel(X,size=20)
    plt.ylabel(Y,size=20)
    plt.title(title2)

    plt.subplot(223)
    plt.scatter(x1,y1,s=0.01)
    plt.xlim([-2,0.5])
    plt.ylim([-1,0.5])
    plt.xlabel(X,size=20)
    plt.ylabel(Y,size=20)
    plt.title(title1)
    
    
    plt.subplot(224)
    plt.scatter(x2,y2,s=0.01)
    plt.xlim([-2,0.5])
    plt.ylim([-1,0.5])
    plt.xlabel(X,size=20)
    plt.ylabel(Y,size=20)
    plt.title(title2)

    plt.show()

if(Zoom=="no"):
    print "sin zoom\n"

    X="$ "+ str(sys.argv[-3]) + " $"
    Y="$ "+ str(sys.argv[-2]) + " $"
    title1="$ "+"Rungekutta" + " $"
    title2="$ "+"Simplectic" + " $"
    
    
    plt.figure(figsize=[20,8])

    plt.subplot(121)
    plt.scatter(x1,y1,s=0.01)
    #plt.xlim([-3,3])
    plt.ylim([-2.5,2.5])
    plt.xlabel(X,size=20)
    plt.ylabel(Y,size=20)
    plt.title(title1)
    
    
    plt.subplot(122)
    plt.scatter(x2,y2,s=0.01)
    #plt.xlim([-3,3])
    plt.ylim([-2.5,2.5])
    plt.xlabel(X,size=20)
    plt.ylabel(Y,size=20)
    plt.title(title2)

    plt.show()

   
