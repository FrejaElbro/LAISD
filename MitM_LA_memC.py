import collections
import scipy.optimize as opt
import warnings
warnings.filterwarnings('ignore',category=RuntimeWarning)
from background import *

set_vars = collections.namedtuple('Our', 'wr wc d ell i kr')
def our(f) : return wrap(f, set_vars)

kc = lambda x: k+x.ell - x.kr

perms = lambda x: binomH(1., w) - binomH(x.kr, x.wr) - binomH(kc(x),x.wc)- binomH(1-k-x.ell, w-x.wc-x.wr)

wr0 = lambda x: x.wr/2+x.d
w0 = lambda x: wr0(x)+x.wc/2
L0 = lambda x: binomH(x.kr,wr0(x)) + binomH(kc(x)/2,x.wc/2)
T1 = lambda x: 2*L0(x) - (x.ell-w0(x))*log2(q)

def r1(x):
    if q == 2:
        return binomH(x.wr,x.wr/2) + binomH(x.kr-x.wr,x.d)
    else:
        return MeurerR(x.kr,x.wr,x.d,x.i,q)

P  = lambda x: min(0,r1(x)-w0(x)*log2(q))


constraints = [
    { 'type' : 'ineq',   'fun' : our(lambda x : k + x.ell - x.kr)},
    { 'type' : 'ineq',   'fun' : our(lambda x : x.kr - x.wr - x.d - 0.00001)},
    { 'type' : 'ineq',   'fun' : our(lambda x : kc(x) - x.wc)},
    { 'type' : 'ineq',   'fun' : our(lambda x : 1-k-x.ell- (w-x.wr-x.wc))},
    { 'type' : 'ineq',   'fun' : our(lambda x : 1-k-x.ell)},
    { 'type' : 'ineq',   'fun' : our(lambda x : w-x.wr-x.wc)},
    { 'type' : 'ineq',   'fun' : our(lambda x : x.wr/2-x.i)},
    { 'type' : 'ineq',   'fun' : our(lambda x : x.d-x.i)},
    { 'type' : 'ineq',   'fun' : our(lambda x : x.ell-w0(x))},
    { 'type' : 'ineq',   'fun' : our(lambda x : memcons-L0(x))}
]

def memory(x):
    return L0(x)

def time(x):
    x = set_vars(*x)  
    return -P(x) + perms(x) + max(L0(x),T1(x))

def optimize(verb=False):
    variables=6
    start = r(0,0.008,variables)
    bounds = [(0, 1)]*variables
    
    result = opt.minimize(time, start, 
            bounds= bounds, tol=1e-10, 
            constraints=constraints, options={'maxiter':1000})
   
    return result


def ISD(x,y,z,memory_constraint = 10000):
    global k
    global w
    global q
    global memcons
    k = x
    w = y
    q = z
    memcons = memory_constraint
    mini=1
    succeed = 0
    tries = 0
    while succeed < 50: #for i in range(500):
        x=optimize()
        tries +=1
        #if tries % 100 == 0:
            #print("Number of tries", tries)
        if x.success and 0<x.fun:
            succeed += 1
            #print("successes", succeed)
        if x.success and 0<x.fun<mini:
            mini=x.fun
            result = x       
    return result, set_vars( *result.x)

