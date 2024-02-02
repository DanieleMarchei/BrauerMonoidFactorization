#Suggestion: use PyPy for performance improvements

from pmatch import brauer_diagrams
from brauermonoid import *
from tqdm import tqdm

def get_hooks_of_size_one(b):
    hs = []
    for edge in b:
        if edge_type(edge) == EdgeType.upper_hook and size(edge) == 1:
            hs.append(edge)
    
    return hs
        
def normalize(t):
    x,y = sorted(t, key = lambda x : abs(x))
    if x < 0 and y > 0:
        x,y = y,x
    
    return (x,y)

def double_fact(N):
  if N == 0 or N == 1:
    return 1
  return double_fact(N-2) * N

for n in range(3,12):
    print(n)
    size_bn = double_fact(2*n - 1)

    max_found = [-1, set()]

    for b in tqdm(brauer_diagrams(n), total = size_bn):
        
        hs = get_hooks_of_size_one(b)
        if len(hs) == 0: continue

        inv = [normalize(t) for t in b]

        t = Tangle(inv)
        l = len(t)
        max_merges_of_t = -1
        for h in hs:
            h = normalize(h)

            n_possible_merges_of_h = 0

            for edge in t.inv:
                if edge == h: continue
                _t = t.copy()
                new_edges = merge(_t, h, edge)
                if new_edges is None: continue # merge is not defined

                e1,e2 = new_edges

                for d in _t.inv:
                    if d == e1: continue
                    if _t.are_intersecting_edges(d,e1):
                        _t.n_crossings[e1] += 1
                        _t.n_crossings[d] += 1
                
                for d in _t.inv:
                    if d == e2: continue
                    if _t.are_intersecting_edges(d,e2):
                        _t.n_crossings[e2] += 1
                        _t.n_crossings[d] += 1
                
                l_new = len(_t)
                if l == l_new + 1:
                    n_possible_merges_of_h += 1
            
            if max_merges_of_t < n_possible_merges_of_h:
                max_merges_of_t = n_possible_merges_of_h
            
        if max_found[0] == max_merges_of_t:
            max_found[1].add(str(inv))
        elif max_found[0] < max_merges_of_t:
            max_found = [max_merges_of_t, {str(inv)}]

    n_max_merges = max_found[0]
    print(f"There are {len(max_found[1])} tangles with {max_found[0]} possible merges.")

# 2
# There are 1 tangles with 1 possible merges.
# 3
# There are 6 tangles with 1 possible merges.
# 4
# There are 2 tangles with 2 possible merges.
# 5
# There are 48 tangles with 2 possible merges.
# 6
# There are 18 tangles with 3 possible merges.
# 7
# There are 936 tangles with 3 possible merges.
# 8
# There are 360 tangles with 4 possible merges.
# 9
# There are 32400 tangles with 4 possible merges.
# 10
# There are 12600 tangles with 5 possible merges.
