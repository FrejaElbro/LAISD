import collections
import scipy.optimize as opt
import warnings
warnings.filterwarnings('ignore',category=RuntimeWarning)
from background import *

set_vars = collections.namedtuple('Our', 'wp ell')
def our(f) : return wrap(f, set_vars)

perms = lambda x: binomH(1., w) - binomH(k+x.ell, x.wp) - binomH(1-k-x.ell, w-x.wp)

L0 = lambda x: binomH((k+x.ell)/2,x.wp/2)+ (x.wp/2)*log2(q-1)
L1 = lambda x: 2*L0(x) - x.ell*log2(q)

constraints = [
    { 'type' : 'ineq',   'fun' : our(lambda x : k + x.ell - x.wp)},
    { 'type' : 'ineq',   'fun' : our(lambda x : 1-k-x.ell- (w-x.wp))},
    { 'type' : 'ineq',   'fun' : our(lambda x : 1-k-x.ell)},
    { 'type' : 'ineq',   'fun' : our(lambda x : w-x.wp)},
]

def memory(x):
    return L0(x)

def time(x):
    x = set_vars(*x)  
    return perms(x) + max(L0(x),L1(x))

def optimize(verb=False):
    variables=2
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

