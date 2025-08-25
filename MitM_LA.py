import collections
import scipy.optimize as opt
import warnings
warnings.filterwarnings('ignore',category=RuntimeWarning)
from background import *

set_vars = collections.namedtuple('Our', 'wp d ell i')
def our(f) : return wrap(f, set_vars)

perms = lambda x: binomH(1., w) - binomH(k+x.ell, x.wp) - binomH(1-k-x.ell, w-x.wp)

w0 = lambda x: x.wp/2+x.d

L0 = lambda x: binomH(k+x.ell,w0(x))
T1 = lambda x: 2*L0(x) - (x.ell-w0(x))*log2(q)
P  = lambda x: min(0,r1(x)-w0(x)*log2(q))

def r1(x):
    if q == 2:
        return binomH(x.wp,x.wp/2) + binomH(k+x.ell-x.wp,x.d)
    else:
        return binomH(x.wp-2*x.i,x.wp/2 - x.i) + binomH(x.wp,2*x.i) + 2*x.i*log2(q-2) + binomH(k+x.ell-x.wp,x.d-x.i) + (x.d-x.i)*log2(q-1)



constraints = [
    { 'type' : 'ineq',   'fun' : our(lambda x : k + x.ell - x.wp - x.d)},
    { 'type' : 'ineq',   'fun' : our(lambda x : 1-k-x.ell- (w-x.wp))},
    { 'type' : 'ineq',   'fun' : our(lambda x : 1-k-x.ell)},
    { 'type' : 'ineq',   'fun' : our(lambda x : w-x.wp)},
    { 'type' : 'ineq',   'fun' : our(lambda x : x.wp/2-x.i)},
    { 'type' : 'ineq',   'fun' : our(lambda x : x.d-x.i)},
    { 'type' : 'ineq',   'fun' : our(lambda x : x.ell-w0(x))}
]

def memory(x):
    return L0(x)

def time(x):
    x = set_vars(*x)  
    return -P(x) + perms(x) + max(L0(x),T1(x))

def optimize(verb=False):
    variables=4
    start = r(0,0.008,variables)
    bounds = [(0, 1)]*variables
    
    result = opt.minimize(time, start, 
            bounds= bounds, tol=1e-10, 
            constraints=constraints, options={'maxiter':1000})  
    return result


def ISD(x,y,z):
    global k
    global w
    global q
    k = x
    w = y
    q = z
    mini=1
    succeed = 0
    tries = 0
    while succeed < 50: #for i in range(500):
        x=optimize()
        tries +=1
        #if tries % 100 == 0:
            #print("Number of tries", tries)
        if x.success:
            succeed += 1
            #print("successes", succeed)
        if x.success and x.fun<mini:
            mini=x.fun
            result = x       
    return result, set_vars( *result.x)


