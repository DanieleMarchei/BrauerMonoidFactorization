import sage.all
import sage.combinat.diagram_algebras as da
from tqdm import tqdm

import pickle as pkl
from sys import argv

def n_T_primes(f):
    s = 0
    for fact in f:
        if "T" in fact:
            s += 1
    return s

def are_intersecting_edges(tangle, edge1, edge2):
    def phi(x):
        i = abs(x)
        if x < 0: return i

        return 2*tangle.order() - i + 1

    _e1 = sorted((phi(edge1[0]), phi(edge1[1])))
    _e2 = sorted((phi(edge2[0]), phi(edge2[1])))

    b1 = _e1[0] < _e2[0] < _e1[1] < _e2[1]
    b2 = _e2[0] < _e1[0] < _e2[1] < _e1[1]

    return b1 or b2

def n_intersecting_edges(tangle, edge):
    n = 0
    for e in tangle._base_diagram:
        if e == edge: continue
        if are_intersecting_edges(tangle, edge, e):
            n += 1
    return n

def n_crossings(tangle):
    s = 0
    for edge in tangle._base_diagram:
        s += n_intersecting_edges(tangle, edge)
    
    return s // 2

_, n = argv

bd = da.BrauerDiagrams(int(n))

with open(f"../../BN/b{bd.order}.pkl", "rb") as f:
    print(f"Loading BN/b{n}.pkl")
    bn = pkl.load(f)
    counterexample_found = False
    for tangle_edges, factorization in tqdm(bn.items()):
        tangle = bd(tangle_edges)
        nc = n_crossings(tangle)
        nt = n_T_primes(factorization)
        if nc != nt:
            counterexample_found = True
            print("Found counterexample")
            print(tangle)
            print(f"Number of crossings = {nc}")
            print(f"Number of T-primes = {nt}")

if not counterexample_found:
    print(f"No counterexamples found in B{n}")

