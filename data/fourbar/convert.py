""" Normalize the value of an angle 
"""
from scipy import *
import math

def convert(x):
    pi2=2*math.radians(180)
    while abs(x) > pi2:
        if x<0:
            x = x+pi2
        elif x>0:
            x=x-pi2
    return x            
