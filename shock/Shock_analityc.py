
# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

#%pylab inline 

# <codecell>

import numpy as np
import matplotlib.pyplot 
import sys
import math
from scipy.optimize import fsolve

# <codecell>


#referencia: "http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html"

#tiempo de la solucion
t= 1;


#condiciones iniciales
x0 = 0;
rho_l = 1;
P_l = 1;
u_l = 0;

rho_r = 0.125;
P_r = 0.1;
u_r = 0;

gamma = 1.4;
mu = np.sqrt( (gamma-1)/(gamma+1) );

#vel sonido
c_l = math.pow( (gamma*P_l/rho_l),0.5);
c_r = math.pow( (gamma*P_r/rho_r),0.5);

print c_l
print c_r



def sod_func(P):
                 
   return  (P - P_r)*((1-mu*mu)*(rho_r*(P + mu*mu*P_r))**-1)**-(0.5)-(math.pow(P_l,(gamma-1)/(2*gamma))-math.pow(P,(gamma-1)/(2*gamma)))*(((1-mu*mu*mu*mu)*P_l**(1/gamma))*(mu*mu*mu*mu*rho_l)**-1)**(0.5);

# <codecell>

P_post = fsolve(sod_func,0.31);

print P_post

# <codecell>

v_post = 2*(np.sqrt(gamma)/(gamma - 1))*(1 - math.pow(P_post, (gamma - 1)/(2*gamma)));
rho_post = rho_r*(( (P_post/P_r) + mu**2 )/(1 + mu*mu*(P_post/P_r)));
v_shock = v_post*((rho_post/rho_r)/( (rho_post/rho_r) - 1));
rho_middle = (rho_l)*math.pow((P_post/P_l),1/gamma);

#posiciones importantes

x1 = x0 - c_l*t;
x3 = x0 + v_post*t;
x4 = x0 + v_shock*t;


c_2 = c_l - ((gamma - 1)/2)*v_post;
x2 = x0 + (v_post - c_2)*t;

#graficamos
n_points = 5000;    

x_min = -0.5;
x_max = 0.5;

x = np.linspace(x_min,x_max,n_points);
rho = np.zeros(n_points);  
P = np.zeros(n_points); 
u = np.zeros(n_points); 
e = np.zeros(n_points); 

for i in range(n_points):
    if (x[i] < x1):
       
        rho[i]=(rho_l)
        P[i]=(P_l)
        u[i]=(u_l)
    if (x1 <= x[i]  and x[i] <= x2):
       
        c = mu*mu*((x0 - x[i])/t) + (1 - mu*mu)*c_l 
        rho[i]=(rho_l*math.pow((c/c_l),2/(gamma - 1)))
        P[i]=(P_l*math.pow((rho[i]/rho_l),gamma))
        u[i]=((1 - mu*mu)*( (-(x0-x[i])/t) + c_l))
    if (x2 <=x[i] and x[i] <= x3):
        
        rho[i]=rho_middle;
        P[i] = P_post;
        u[i] = v_post;
    if (x3 <= x[i] and x[i] <= x4):
      
        rho[i] = rho_post;
        P[i] = P_post;
        u[i] = v_post;
    if (x4 < x[i]):
      
        rho[i] = rho_r;
        P[i] = P_r;
        u[i] = u_r;
    
    e[i] = P[i]/((gamma - 1)*rho[i]);

# <codecell>


# <codecell>

f= open('datos.txt', 'w')

for i in range(n_points):
    coso = str(i) +','+str(P[i])+','+str(rho[i])+','+str(u[i])+','+str(e[i])+'\n'
    f.write(coso)
    

# <codecell>


