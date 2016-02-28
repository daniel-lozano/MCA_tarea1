import numpy as np
import matplotlib.pyplot as plt
import sys


if(len(sys.argv)!=8):
    print "introduzca 6 parametros para correr el programa\n"
    print "parametros dados=", len(sys.argv)
    
    exit()
    

lectura1=str(sys.argv[1])
lectura2=str(sys.argv[2])
c1=int(sys.argv[3])
c2=int(sys.argv[4])
Zoom=str(sys.argv[-1])
print "c1=",c1, "c2=", c2

#existencia del archivo de entrada
try:
    archivo1=open(lectura1).read().split("\n")
except IOError:
    
   print "Error: el archivo "+lectura+ " no existe"
   exit()
try:
    archivo2=open(lectura2).read().split("\n")
except IOError:
    
   print "Error: el archivo "+lectura+ " no existe"
   exit()

t1=[]
t2=[]
x1=[]
y1=[]
x2=[]
y2=[]
E1=[]
E2=[]

for i in range(int(len(archivo1)-1)):
    a= archivo1[i].split()
    t1.append(a[0])
    x1.append(a[c1])
    y1.append(a[c2])
    E1.append(0.5*float(a[c1-1])**2-0.5*(4*float(a[c2-1])**2+0.5**2)**(-0.5))

for i in range(int(len(archivo2)-1)):
    a= archivo2[i].split()
    t2.append(a[0])
    x2.append(a[c1])
    y2.append(a[c2])
    E2.append(0.5*float(a[c1-1])**2-0.5*(4*float(a[c2-1])**2+0.5**2)**(-0.5))
      

 
    

if(Zoom=="si"):
    print "con zoom\n"
    X="$ "+ str(sys.argv[-3]) + " $"
    Y="$ "+ str(sys.argv[-2]) + " $"
    title1="$ "+"Rungekutta" + " $"
    title2="$ "+"Simplectic" + " $"
    
    
    plt.figure(figsize=[20,8])

    plt.subplot(221)
    plt.scatter(x1,y1,s=0.01)
    plt.xlim([-3,3])
    plt.ylim([-2.5,2.5])
    plt.xlabel(X,size=20)
    plt.ylabel(Y,size=20)
    plt.title(title1)
    
    plt.subplot(222)
    plt.scatter(x2,y2,s=0.01)
    plt.xlim([-3,3])
    plt.ylim([-2.5,2.5])
    plt.xlabel(X,size=20)
    plt.ylabel(Y,size=20)
    plt.title(title2)

     
    plt.subplot(223)
    plt.scatter(x1,y1,s=0.01)
    plt.xlim([-2.5,-1])
    plt.ylim([-1,1])
    plt.xlabel(X,size=20)
    plt.ylabel(Y,size=20)
   
    
    
    plt.subplot(224)
    plt.scatter(x2,y2,s=0.01)
    plt.xlim([-2.5,-1])
    plt.ylim([-1,1])
    plt.xlabel(X,size=20)
    plt.ylabel(Y,size=20)
    
    
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
    plt.xlim([-3,3])
    plt.ylim([-2.5,2.5])
    plt.xlabel(X,size=20)
    plt.ylabel(Y,size=20)
    plt.title(title1)
    
    
    plt.subplot(122)
    plt.scatter(x2,y2,s=0.01)
    plt.xlim([-3,3])
    plt.ylim([-2.5,2.5])
    plt.xlabel(X,size=20)
    plt.ylabel(Y,size=20)
    plt.title(title2)

    plt.show()

   
if(Zoom=="energia"):
    print "Energia\n"
    X="$ "+ str(sys.argv[-3]) + " $"
    Y="$ "+ str(sys.argv[-2]) + " $"
    title1="$ "+"Rungekutta" + " $"
    title2="$ "+"Simplectic" + " $"
    
    
    plt.figure(figsize=[20,8])

    plt.subplot(211)
    plt.plot(t1,E1)
    plt.xlabel("$ time $",size=20)
    plt.ylabel("$ Energy\ RK\ $",size=20)
    
    plt.subplot(212)
    plt.plot(t2,E2)
    plt.xlabel("$ time $",size=20)
    plt.ylabel("$ Energy\ Simp\ $",size=20)

    plt.show()
