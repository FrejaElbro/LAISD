import collections
import random
from random import  uniform as ru
from math import *
import scipy.optimize as opt
import math 
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import minimize
from copy import copy
import numpy as np
import warnings
from scipy.special import binom
warnings.filterwarnings('ignore',category=RuntimeWarning)

def log2(x):
    """
    Binary logarithm
    """
    if x <= 0:
        return -1000
    else:
        return log(x,2)

def Hi(v):
    """
    Inverse of the binary entropy function
    """
    if 1 <= v:
        return 0.5
    if 0 >= v :
        return 0
    return fsolve(lambda x:v -(-x*log2(x)-(1-x)*log2(1-x)),0.0000001)[0]

def Hq(x,q):
    """
    q-ary entropy function   
    """
    if x>=1:
        return log(q-1,q)
    if 0 >= x :
        return 0
    return x*log(q-1,q)-x*log(x,q)-(1-x)*log(1-x,q)

def Hqi(v,q):
    """
    Inverse of the q-ary entropy function
    """
    if v>=1:
        return 1-1/q
    if 0 >= v :
        return 0
    return fsolve(lambda x:v -(x*log(q-1,q)-x*log(x,q)-(1-x)*log(1-x,q)),0.0000001)[0]

def H(c):
    """
    Binary Entropy function
    """
    if c <= 0. or c >= 1.:
        return 0.
    
    return -(c * log2(c) + (1 - c) * log2(1 - c))

def binomH(n,k):
    """
    binomial coefficient entropy approximation
    """
    if (n<=0) or (k<=0) or (k>=n):
        return 0.
    return n * H(k/n)

def multiH(n,c):
    """
    multinomial coefficient entropy approximation
    """
    if sum(c)>n:
        return 0
    tot=0
    val=n
    for i in c:
        tot+=binomH(n,i)
        n-=i
    return tot

def wrap(f,g) :
    """
    Useful for working with named tuples
    """
    def inner(x):
        return f(g(*x))
    return inner

def r(x,y,z):
    """
    Shorthand for setting up random starting points
    """
    return [(ru(x,y)) for i in range(z)]


def check_constraints(constraints, solution) : 
    """
    Helper function to check if a solution satisfies the constraints.
    """
    return [ (constraint['type'], constraint['fun'](solution)) for constraint in constraints ]


def validity(mycons, result):
    """
    Check if a solution satisfies the constraints.
    """
    const= check_constraints(mycons, result)
    for i in const:
        if i[0]=="eq":
            if abs(i[1])>1e-7:
                return False
        elif i[1]<-1e-7:
            return False
    return True



def MeurerR(kl,p,d,i,q):
    """
    Calculate the number of representations e' = e1 + e2 over F_q.
    
    Parameters:
    - kl: vector length
    - p: weight of e'
    - d: wt(e1) - p/2 = wt(e2) - p/2 
    - i: Index of largest summand
    - q: field size
    
    Reference: Meurer's thesis, Lemma 7.1.2
    """
    return binomH(p-2*i,p/2 - i) + binomH(p,2*i) + 2*i*log2(q-2) + binomH(kl-p,d-i) + (d-i)*log2(q-1)
