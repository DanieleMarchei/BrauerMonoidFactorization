#walks the Cayley Graph of Bn and generates a dictionaty tangle:factorization
#where the factorization has the least amount of T-primes possible


import sage.all
import sage.combinat.diagram_algebras as da
from tqdm import tqdm

from collections import deque
import pickle as pkl
from sys import argv

def double_fact(N):
    if N == 0 or N == 1:
        return 1
    return double_fact(N-2) * N

def generators(bd):
    N = bd.order
    gens = []
    
    for i in range(1,N):
        p = [{i,i+1}, {-i, -(i+1)}]
        rest = [{j, -j} for j in range(1,N+1) if j != i and j != i+1]
        ui = bd(p + rest)
        name = f"U{i}"
        gens.append((ui, name))
    
    for i in range(1,N):
        p = [{i,-(i+1)}, {i+1, -i}]
        rest = [{j, -j} for j in range(1,N+1) if j != i and j != i+1]
        ti = bd(p + rest)
        name = f"T{i}"
        gens.append((ti, name))
    
    identity = None
    for i in range(1,N):
        rest = [{j, -j} for j in range(1,N+1)]
        identity = bd(rest)
    
    return gens, identity

_, n = argv

bd = da.BrauerDiagrams(int(n))
gens, identity = generators(bd)

queue = deque()
queue.appendleft((identity, []))

black_nodes = {identity.base_diagram() : []}

limit = double_fact(2*bd.order-1)

with tqdm(total=int(limit)) as pbar:
  while len(queue) > 0:
        
        tangle, factorization = queue.pop()

        pbar.update(1)

        for prime_tangle, prime_name in gens:
            neigh, loops_removed = tangle.compose(prime_tangle, check = False)
            neigh_factorization = factorization + [prime_name]
            
            if neigh.base_diagram() not in black_nodes:
                black_nodes[neigh.base_diagram()] = neigh_factorization
                queue.appendleft((neigh, neigh_factorization))

with open(f"BN/b{bd.order}.pkl", "wb") as f:
    pkl.dump(black_nodes, f)
