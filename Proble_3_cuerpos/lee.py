import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import sys


lectura=str(sys.argv[1])
coor1=int(sys.argv[2])
coor2=int(sys.argv[3])
archivo=open(lectura).read().split("\n")

x=[]
y=[]

for i in range(len(archivo)-1):
    a= archivo[i].split()
    #print a
    x.append(a[coor1])
    y.append(a[coor2])
    

print len(x),len(y)

X="$ "+ str(sys.argv[-2]) + " $"
Y="$ "+ str(sys.argv[-1]) + " $"

plt.plot(x,y,"k.")
plt.xlabel(X,size=20)
plt.ylabel(Y,size=20)
plt.show()
plt.savefig("grafica.png")
plt.close()
