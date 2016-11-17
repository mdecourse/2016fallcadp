"""
Kinematic analysis of four-bar linkage

Author: E. Pennestri'  pennestri@mec.uniroma2.it
Dipartimento Ingegneria Industriale
Universita' di Roma Tor Vergata
"""
from scipy.optimize import fsolve
from scipy import *
import math
import numpy as np
import scipy
from pylab import *

import fourbarfunc
import convert
close()
"""
 l1,l2,l3,l4 : Link lengths
"""
l1=1
l2=5
l3=4
l4=4
#  Initial and final angles of driving link
th1min=0
th1max=360
guess=[0,0]
# Angular speed of the crank
om1=1
#
th4=math.radians(180)
# Declare arrays
th1p=array([])
th2=array([])
th3=array([])
om2=array([])
om3=array([])
alpha2=array([])
alpha3=array([])

for th1d in range(th1min,th1max):
    th1=math.radians(th1d)
    c=fsolve(fourbarfunc.residuals,guess,args=(th1,th4,l1,l2,l3,l4))
    c[0]=convert.convert(c[0])
    c[1]=convert.convert(c[1])
    guess=[c[0],c[1]]  # Adopt current solution as guess solution
    # Angular velocity analysis
    Jac=array([[-l2*math.sin(c[0]),-l3*math.sin(c[1])],\
               [l2*math.cos(c[0]),l3*math.cos(c[1])]])   
    b=array([-l1*sin(th1)*om1,l1*cos(th1)*om1])
    om = linalg.solve(Jac,b)  
    om2=append(om2,om[0])    
    om3=append(om3,om[1])
    # Angular acceleration analysis
    b1=array([-l1*cos(th1)*pow(om1,2)-l2*cos(c[0])*pow(om[0],2)-l3*cos(c[1])*pow(om[1],2),\
    -l1*sin(th1)*pow(om1,2)-l2*sin(c[0])*pow(om[0],2)-l3*sin(c[1])*pow(om[1],2)])
    alpha = linalg.solve(Jac,b1)  
    alpha2=append(alpha2,alpha[0])    
    alpha3=append(alpha3,alpha[1])    
    
    th1p=append(th1p,th1d)
    th2=append(th2,math.degrees(c[0]))
    th3=append(th3,math.degrees(c[1]))

# end of loop 

figure(1)
p1,=plot(th1p, th2, color="blue", linewidth=1.0, linestyle="-",label="theta2")
p2,=plot(th1p, th3, color="black", linewidth=1.0, linestyle="-",label="theta3")
xlabel('Theta_1 (deg)')
ylabel('Theta_2 and Theta_3 (deg)')
legend([p2, p1], ["Theta_3", "Theta2"])
grid(True)
title('Four-bar angular positions ')   # subplot 


figure(2)
p3,=plot(th1p, om2, color="blue", linewidth=1.0, linestyle="-",label="om2")
p4,=plot(th1p, om3, color="black", linewidth=1.0, linestyle="-",label="om3")
xlabel('Theta_1 (deg)')
ylabel('Omega_2 and Omega_3 (deg)')
legend([p3, p4], ["Omega2", "Omega3"])
grid(True)
title('Four-bar angular velocities')   # subplot 

figure(3)
p5,=plot(th1p, alpha2, color="blue", linewidth=1.0, linestyle="-",label="alpha2")
p6,=plot(th1p, alpha3, color="black", linewidth=1.0, linestyle="-",label="alpha3")
xlabel('Theta_1 (deg)')
ylabel('Alpha2_2 and Alpha_3 (deg)')
legend([p5, p6], ["Alpha2", "Alpha3"])
grid(True)
title('Four-bar angular accelerations')   # subplot 
show()

