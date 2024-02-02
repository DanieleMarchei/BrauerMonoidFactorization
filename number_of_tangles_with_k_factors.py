#Suggestion: use PyPy for performance improvements


import pickle as pkl
from tqdm import tqdm
from brauermonoid import Tangle

from pmatch import brauer_diagrams

def n_components(inv):
    not_covered = set(range(1,len(inv)))
    for edge in inv:
        a,b = sorted([abs(x) for x in edge])
        for i in range(a, b):
            if i in not_covered: not_covered.remove(i)
    
    return len(not_covered) + 1


def normalize(t):
    x,y = sorted(t, key = lambda x : abs(x))
    if x < 0 and y > 0:
        x,y = y,x
    
    return (x,y)

def double_fact(N):
  if N == 0 or N == 1:
    return 1
  return double_fact(N-2) * N
        

n = 8
size_bn = double_fact(2*n - 1)

print(n)

one_component = [0 for _ in range(n*(n-1)//2 + 1)]
all_tangles = [0 for _ in range(n*(n-1)//2 + 1)]

if n <= 8:
    with open(f"BN/b{n}.pkl", "rb") as f:
        print(f"Loading BN/b{n}.pkl")
        bn = pkl.load(f)
        for inv, fact in tqdm(bn.items()):
            l = len(fact)
            all_tangles[l] += 1
            if n_components(inv) == 1:
                one_component[l] += 1
else:

    bd = brauer_diagrams(n)

    for b in tqdm(bd, total = size_bn):
        inv = [normalize(t) for t in b]
        
        tangle = Tangle(inv)
        l = len(tangle)
        all_tangles[l] += 1
        if n_components(inv) == 1:
            one_component[l] += 1

print("Tangles with one component")
print(one_component)

print()

print("All Tangles")
print(all_tangles)


