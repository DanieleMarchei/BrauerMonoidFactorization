from enum import Enum
from functools import partial

class EdgeType(Enum):
    upper_hook = 1
    lower_hook = 2
    zero_transversal = 3
    positive_transversal = 4
    negative_transversal = 5

def edge_type(edge):
    x,y = edge
    if x < 0 and y < 0: return EdgeType.lower_hook
    if x > 0 and y > 0: return EdgeType.upper_hook
    if x > 0 and y < 0 and x > -y: return EdgeType.positive_transversal
    if x > 0 and y < 0 and x == -y: return EdgeType.zero_transversal
    if x > 0 and y < 0 and x < -y: return EdgeType.negative_transversal

def size(edge):
   return abs(abs(edge[0]) - abs(edge[1]))

class Tangle:
    def __init__(self, inv, calculate_crossings = True):
        self.N = len(inv)
        self.inv = set(inv)
        self.inv_dict = {}
        for x,y in inv:
           self.inv_dict[x] = (x,y)
           self.inv_dict[y] = (x,y)
        
        if calculate_crossings:
            self._calculate_crossings()
        
    def get_edge_from(self, node):
       return self.inv_dict[node]

    def get_n_crossings(self, edge):
       return self.n_crossings[edge]

    def _calculate_crossings(self):
        self.n_crossings = {}
        for edge in self.inv:
           self.n_crossings[tuple(edge)] = self.n_intersecting_edges(edge)

    def __len__(self):
        l = 0
        for edge in self.inv:
            l += max(size(edge), self.n_crossings[edge])
        
        return l // 2

    def __contains__(self, edge):
        return edge in self.inv or edge[::-1] in self.inv

    def are_intersecting_edges(self, e1,e2):
        
        def phi(x):
            i = abs(x)
            if x < 0: return i

            return 2*self.N - i + 1

        _e1 = sorted((phi(e1[0]), phi(e1[1])))
        _e2 = sorted((phi(e2[0]), phi(e2[1])))

        b1 = _e1[0] < _e2[0] < _e1[1] < _e2[1]
        b2 = _e2[0] < _e1[0] < _e2[1] < _e1[1]

        return b1 or b2

    def n_intersecting_edges(self, edge):
        n = 0
        for e in self.inv:
            if e == edge:
                continue
            if self.are_intersecting_edges(edge, e):
                n += 1
        return n

    def __repr__(self) -> str:
        return str(self.inv)
    
    def copy(self):
        X_new = Tangle(self.inv, calculate_crossings=False)
        X_new.n_crossings = self.n_crossings.copy()
        return X_new

    def add_edge(self, edge):
        self.inv.add(edge)
        x,y = edge
        self.inv_dict[x] = (x,y)
        self.inv_dict[y] = (x,y)
        self.n_crossings[edge] = 0

    def delete_edge(self, edge):
        self.inv.discard(edge)
        del self.inv_dict[edge[0]]
        del self.inv_dict[edge[1]]
        del self.n_crossings[edge]

def text_to_tangle(text):
    inv = []
    pairs = text.split(",")
    for pair in pairs:
        a,b = pair.split(":")
        if "'" in a:
            a = -int(a.replace("'",""))
        else:
            a = int(a)
        if "'" in b:
            b = -int(b.replace("'",""))
        else:
            b = int(b)
        
        inv.append((a,b))
    
    return Tangle(inv)





def merge(X : Tangle, h, e):
    i = h[0]
    x,y = e
    t = edge_type(e)
    if t == EdgeType.upper_hook and x < i and i+1 < y: 
        e1 = (x,i)
        e2 = (i+1, y)
    elif t == EdgeType.lower_hook and x <= i and i+1 <= abs(y):
        e1 = (i,x)
        e2 = (i+1, y)
    elif t == EdgeType.negative_transversal and x < i and i+1 <= abs(y):
        e1 = (x,i)
        e2 = (i+1, y)
    elif t == EdgeType.positive_transversal and x > i+1 and i >= abs(y):
        e1 = (i,y)
        e2 = (i+1, x)
    else: 
        return None
    
    X.delete_edge(h)
    X.delete_edge(e)

    X.add_edge(e1)
    X.add_edge(e2)

    return e1, e2
    
def compose_with_T(i, X : Tangle):
    a = X.get_edge_from(i)
    b = X.get_edge_from(i+1)

    n_crossings_a = X.n_crossings[a]
    n_crossings_b = X.n_crossings[b]

    X.delete_edge(a)
    X.delete_edge(b)

    idx_i = a.index(i)
    idx_i_1 = b.index(i+1)

    a = list(a)
    a[idx_i] = i+1
    a = tuple(a)

    b = list(b)
    b[idx_i_1] = i
    b = tuple(b)

    ta,tb = edge_type(a), edge_type(b)
    is_transversal = lambda et: et in (EdgeType.positive_transversal, EdgeType.zero_transversal, EdgeType.negative_transversal)

    if ta == EdgeType.upper_hook: a = tuple(sorted(a))
    if is_transversal(ta): a = tuple(reversed(sorted(a)))  

    if tb == EdgeType.upper_hook: b = tuple(sorted(b))
    if is_transversal(tb): b = tuple(reversed(sorted(b)))

    X.add_edge(a)
    X.add_edge(b)

    X.n_crossings[a] = n_crossings_a
    X.n_crossings[b] = n_crossings_b



def node_polarity(X : Tangle, node):
    edge = X.get_edge_from(node)
    t = edge_type(edge)
    is_left_node = edge[0] == node
    if t == EdgeType.upper_hook: return "-" if is_left_node else "+"
    if t == EdgeType.lower_hook: return "+" if is_left_node else "-"
    if t == EdgeType.positive_transversal: return "+"
    if t == EdgeType.zero_transversal: return "0"
    if t == EdgeType.negative_transversal: return "-"

def label_node_polarity(X : Tangle):
    rho = partial(node_polarity, X)
    S1 = []
    S2 = []
    for i in range(1, X.N+1):
        edge = X.get_edge_from(i)
        if X.get_n_crossings(edge) < size(edge):
            S1.append(i)
        edge = X.get_edge_from(-i)
        if X.get_n_crossings(edge) < size(edge):
            S2.append(-i)
    
    polarity_labels = {}
    counter_pos = 0
    counter_neg = 0
    for i in S1:
        p = rho(i)
        if p == "+": 
            counter_pos += 1
            label = f"{p}{counter_pos}"
        if p == "-": 
            counter_neg += 1
            label = f"{p}{counter_neg}"

        polarity_labels[label] = [i]
    
    counter_pos = 0
    counter_neg = 0

    for i in S2:
        p = rho(i)
        if p == "+": 
            counter_pos += 1
            label = f"{p}{counter_pos}"
        if p == "-": 
            counter_neg += 1
            label = f"{p}{counter_neg}"

        polarity_labels[label].append(i)
    
    return polarity_labels

def tau(X : Tangle):
    polarity_labels = label_node_polarity(X)
    Z = []
    for i in range(1, X.N+1):
        edge = X.get_edge_from(i)
        if X.get_n_crossings(edge) >= size(edge):
            Z.append(tuple(edge))
    
    for x,y in polarity_labels.values():
        Z.append((x,y))
    
    return Tangle(Z)

def factorizeSN(X : Tangle):
  s = [-X.get_edge_from(i)[1] for i in range(1, X.N+1)]

  F = []
  for j in range(X.N-1, 0, -1):
    for i in range(0, j):
      if s[i] > s[i+1]:
        s[i], s[i+1] = s[i+1], s[i]
        F.append(f"T{i+1}")

  return F

def factorizeBN(X : Tangle):
    I = [int(f[1:]) for f in factorizeSN(tau(X))]
    if len(I) == 0: return ["I"]
    
    F = []
    for i in I:
        h = (i,i+1)
        if h in X:
            l = len(X)
            for edge in X.inv:
                if edge == h: continue
                X_new = X.copy()
                for d in X_new.inv:
                    if d == h or d == edge: continue
                    if X_new.are_intersecting_edges(d,edge):
                        X_new.n_crossings[d] -= 1
                
                new_edges = merge(X_new, h, edge)
                if new_edges is None: continue # merge is not defined

                e1,e2 = new_edges
                
                for d in X_new.inv:
                    if d == e1: continue
                    if X_new.are_intersecting_edges(d,e1):
                        X_new.n_crossings[e1] += 1
                        X_new.n_crossings[d] += 1
                
                for d in X_new.inv:
                    if d == e2: continue
                    if X_new.are_intersecting_edges(d,e2):
                        X_new.n_crossings[e2] += 1
                        X_new.n_crossings[d] += 1
                
                l_new = len(X_new)
                if l_new == l - 1:
                    F.append(f"U{i}")
                    X = X_new
                    break

        else:
            compose_with_T(i, X)
            F.append(f"T{i}")
            e1 = X.get_edge_from(i)
            e2 = X.get_edge_from(i+1)
            X.n_crossings[e1] -= 1
            X.n_crossings[e2] -= 1
    
    return F
