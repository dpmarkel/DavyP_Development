from __future__ import division
import numpy as np
from scipy.integrate import quad
from scipy.integrate import dblquad
from scipy.integrate import tplquad


def integrand(x):
    return x**2

def double_integrand(y, x):
    'y must be the first argument, and x the second.'
    return y * np.sin(x) + x * np.cos(y)
    


def triple_integrand(z, y, x):
    return y * np.sin(x) + z * np.cos(x)    

def singleIntegral(f, xmin, xmax, dx=0.005):
    """ My python integrationg function"""    
    # f is my function that the integral object will operate on.    
    # An integral is a mathematical object that can be interpreted as an area or a generalization of area.
    val = 0.0    
    # this takes the integration limits and slices them into dx wide chunks.  
    xs = np.arange(xmin , xmax+dx, dx)    
    # now loop through the slices calculating the line integral    
    for x in xs:
        calc = f(x)*dx
        val = (val+calc)        
    return val

def doubleIntegral(f, xmin, xmax, ymin, ymax, dx=0.005, dy=0.005):
    """ My python integrationg function""" 
    # f is my function that the integral object will operate on.    
    # An integral is a mathematical object that can be interpreted as an area or a generalization of area.
    val = 0.0    
    # this takes the integration limits and slices them into dx wide chunks.  
    xs = np.arange(xmin , xmax+dx, dx)
    ys = np.arange(ymin , ymax+dy, dy)    
    # now loop through the slices calculating the line integral    
    for x in xs:
        for y in ys:
            calc = f(y,x)*dy*dx
            val = (val+calc)        
    return val

def tripleIntegral(f, xmin, xmax, ymin, ymax, zmin, zmax, dx=0, dy=0.01, dz=0.01):
    """ My python integrationg function"""
    # f is my function that the integral object will operate on.    
    # An integral is a mathematical object that can be interpreted as an area or a generalization of area.
    val = 0.0    
    # this takes the integration limits and slices them into dx wide chunks.  
    xs = np.arange(xmin, xmax+dx, dx)
    ys = np.arange(ymin, ymax+dy, dy)
    zs = np.arange(zmin, zmax+dz, dz)    
    # now loop through the slices calculating the line integral    
    for x in xs:
        for y in ys:
            for z in zs:
                calc = f(z,y,x)*dy*dx*dz
                val = (val+calc)        
    return val


def main():
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    zmin = -1
    zmax = 1
  
    
    
    # Integration loops
    
    retVal1 = singleIntegral(integrand, xmin, xmax)
    ans1, err1 = quad(integrand, xmin, xmax)
    
    xmin=np.pi
    xmax=2*np.pi
    ymin=0
    ymax=np.pi       
    retVal2 = doubleIntegral(double_integrand, xmin, xmax, ymin, ymax)
    ans2, err2 = dblquad(double_integrand, np.pi, 2*np.pi, lambda x: 0, lambda x: np.pi)
    
    
    xmin=0
    xmax=np.pi
    ymin=0
    ymax=1
    retVal3 = tripleIntegral(triple_integrand, xmin, xmax, ymin, ymax, zmin, zmax)
    # x, y, z limits:  0 <= x <= pi  ;   0 <= y <= 1 ; -1 <= z <= 1
    ans3, err3 = tplquad(triple_integrand, 0, np.pi, lambda x: 0, lambda x: 1, lambda x,y: -1, lambda x,y: 1)
        
    print '***********************************'
    print 'Single integral: %s' % retVal1 
    print 'Quad answer: %s' % ans1    
    
    print 'Double integral: %s' % retVal2
    print 'Double quad answer: %s' % ans2    
        
    print 'Triple integral: %s' % retVal3
    print 'Triple quad answer: %s' % ans3
    print '***********************************'


if __name__== "__main__":
    main()