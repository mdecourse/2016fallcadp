from scipy import *
import math
def residuals(x,th1,th4,l1,l2,l3,l4):
    f=zeros(2)
    th2=x[0]
    th3=x[1]
    f[0] = l1*math.cos(th1)+l2*math.cos(th2)+l3*math.cos(th3)+l4*math.cos(th4)
    f[1] = l1*math.sin(th1)+l2*math.sin(th2)+l3*math.sin(th3)+l4*math.sin(th4)
    return f
  
