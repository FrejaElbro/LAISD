import collections
import scipy.optimize as opt
import warnings
warnings.filterwarnings('ignore',category=RuntimeWarning)
from background import *


set_vars = collections.namedtuple('Our', 'p d ell i')
def our(f) : return wrap(f, set_vars)


perms = lambda x: binomH(1., w) - binomH(k+x.ell, x.p) - binomH(1-k-x.ell, w-x.p)

p1 = lambda x: x.p/2+x.d

L0 = lambda x: binomH((k+x.ell)/2,p1(x)/2)+ (p1(x)/2)*log2(q-1)
L1 = lambda x: 2*L0(x) - r1(x)
T2 = lambda x: 2*L1(x) - (log2(q)*x.ell-r1(x))

def r1(x):
    if q == 2:
        return binomH(x.p,x.p/2) + binomH(k+x.ell-x.p,x.d)
    else:
        return binomH(x.p-2*x.i,x.p/2 - x.i) + binomH(x.p,2*x.i) + 2*x.i*log2(q-2) + binomH(k+x.ell-x.p,x.d-x.i) + (x.d-x.i)*log2(q-1)


constraints = [
    { 'type' : 'ineq',   'fun' : our(lambda x : k + x.ell - x.p - x.d)},
    { 'type' : 'ineq',   'fun' : our(lambda x : 1-k-x.ell- (w-x.p))},
    { 'type' : 'ineq',   'fun' : our(lambda x : 1-k-x.ell)},
    { 'type' : 'ineq',   'fun' : our(lambda x : w-x.p)},
    { 'type' : 'ineq',   'fun' : our(lambda x : x.p/2-x.i)},
    { 'type' : 'ineq',   'fun' : our(lambda x : x.d-x.i)},
    { 'type' : 'ineq',   'fun' : our(lambda x : x.ell*log2(q)-r1(x))}
]

def memory(x):
    return max(L0(x),L1(x))

def time(x):
    x = set_vars(*x)  
    return perms(x) + max(L0(x),L1(x),T2(x))

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
    mini=100000
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
