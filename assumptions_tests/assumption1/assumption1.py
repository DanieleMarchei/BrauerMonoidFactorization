import pickle as pkl
import random
from sys import argv
from tqdm import tqdm
import maude

def init():
    ''' Initilizes the maude package.
    '''
    
    maude.init()
    maude.load('brauer.maude')

def reduce(factorization, max_patience = "auto", dfs = True):
    ''' Reduces the input term as an equivalent and (possibly) smaller factorization. \n
    Arguments:
        factorization : the factor list to be reduced \n
        max_patience : integer, str (default: auto), how many steps inside a local minima are you willing to perform. If "auto" is equal to 1000 * length of term
    '''

    if len(factorization) == 0: return []
    if len(factorization) == 1: return factorization

    factorization = [factor for factor in factorization if factor != "I"]

    # covert factor list to Maude term
    term = []
    for factor in factorization:
        if factor == "I": continue

        term.append(f" {factor[0]} {factor[1:]} ")

    term = ",".join(term)


    # strategy_rule = "(delete ! ; move)*"
    strategy_rule = "(delete ! ; try(move) ) *"

    brauer = maude.getModule("REW-RULES")
    try:
        t = brauer.parseTerm(term)
    except:
        raise Exception("An error has occurred. Have you called init() before reduce()?")
        

    t_len = brauer.parseTerm(f"length({term})")
    t_len.reduce()
    t_len = t_len.toInt()

    if max_patience == "auto":
        max_patience = 1000 * t_len

    # calculate the max number of factors the term can have
    N = brauer.parseTerm(f"max_idx({term})")
    N.reduce()
    N = N.toInt() + 1
    factors_upper_bound = N * (N - 1) // 2

    strategy = brauer.parseStrategy(strategy_rule)
    
    # perform srew with depth first search
    results = t.srewrite(strategy, dfs)
    current_patience = max_patience
    for res in results:
        # calculate result length
        l = brauer.parseTerm(f"length({res[0]})")
        l.reduce()
        l = l.toInt()

        if l < t_len:
            # covert Maude term to factor list
            reduced_factorization = []
            for factor in str(res[0]).split(","):
                reduced_factorization.append(factor.replace(" ", ""))
            
            return reduce(reduced_factorization, max_patience)
        elif l > t_len:
            # stop looping because we have reached the bottom of the tree
            break
        elif l <= factors_upper_bound:
            # do not start to lose patience until we have reached the upper bound
            current_patience -= 1
        
        if current_patience <= 0:
            # stop looping if patience was lost
            break

    # no better solution was found

    return factorization

# This is a bare bones implementation because 
# the maude package has some conflicts with the sagemath package
# and cannot be loaded together apparently 
class Tangle:

    def __init__(self, edges) -> None:
        self.edges = edges
    
    def delete_edge_containing(self, i):
        edge = None
        for e in self.edges:
            if i in e:
                edge = e
                break

        if edge is not None: 
            self.edges.remove(edge)
        
        return edge

    def compose(self, prime_name):
        type, i = prime_name[0], int(prime_name[1:])
        
        e1 = self.delete_edge_containing(-i)
        e2 = self.delete_edge_containing(-(i+1))
        if e2 is None:
            # e2 is None because e1 was a lower hook of size 1
            self.edges.append(e1)
            return self
        
        if type == "T":
            if -i in e1:
                e1.remove(-i)
                e1.add(-(i+1))
                e2.remove(-(i+1))
                e2.add(-i)
            else:
                e1.remove(-(i+1))
                e1.add(-i)
                e2.remove(-i)
                e2.add(-(i+1))
            self.edges.append(e1)
            self.edges.append(e2)
        else:
            self.edges.append({-i, -(i+1)})
            other_edge = {d for d in e1.union(e2) if d not in (-i, -(i+1))}
            self.edges.append(other_edge)
        
        return self

    def base_diagram(self):
        # same code used by the SageMath's implementation
        return tuple(sorted(tuple(sorted(i)) for i in self.edges))
    
    def is_symmetric(self):
        is_transversal = lambda edge: sorted(edge)[0] < 0 and sorted(edge)[1] > 0
        return all(map(is_transversal, self.edges))
    
    def are_intersecting_edges(self, edge1, edge2):
        def phi(x):
            i = abs(x)
            if x < 0: return i

            return 2*len(self.edges) - i + 1

        _e1 = sorted((phi(edge1[0]), phi(edge1[1])))
        _e2 = sorted((phi(edge2[0]), phi(edge2[1])))

        b1 = _e1[0] < _e2[0] < _e1[1] < _e2[1]
        b2 = _e2[0] < _e1[0] < _e2[1] < _e1[1]

        return b1 or b2

    def n_intersecting_edges(self, edge):
        n = 0
        for e in tangle.base_diagram():
            if e == edge: continue
            if self.are_intersecting_edges(edge, e):
                n += 1
        return n

    def n_crossings(self):
        s = 0
        for edge in tangle.base_diagram():
            s += self.n_intersecting_edges(edge)
        
        return s // 2
    
    def is_TL(self):
        return self.n_crossings() == 0

    @staticmethod
    def t(i, n):
        edges = [{j, -j} for j in range(1, n+1) if j not in (i, i+1)]
        edges.append({i, -(i+1)})
        edges.append({i+1, -i})
        return Tangle(edges)
    
    @staticmethod
    def u(i, n):
        edges = [{j, -j} for j in range(1, n+1) if j not in (i, i+1)]
        edges.append({i, i+1})
        edges.append({-i, -(i+1)})
        return Tangle(edges)

    @staticmethod
    def id(n):
        edges = [{j, -j} for j in range(1, n+1)]
        return Tangle(edges)


    def __repr__(self) -> str:
        return str(self.base_diagram())

def generators(n):
    gens = []
    
    for i in range(1,n):
        name = f"U{i}"
        gens.append((name, Tangle.u(i, n)))
    
    for i in range(1,n):
        name = f"T{i}"
        gens.append((name, Tangle.t(i, n)))
    
    return gens

def random_factorization(n, gens, scale = 2):
    fact_len = random.randint(2, int(scale * n*(n-1)//2)+1)
    tangle = Tangle.id(n)
    fact = []
    for _ in range(fact_len):
        prime_name, prime = random.choice(gens)
        tangle.compose(prime_name)
        fact.append(prime_name)
    
    return tangle, fact




def double_fact(N):
    if N == 0 or N == 1:
        return 1
    return double_fact(N-2) * N

_, n = argv[0], int(argv[1]) 

init()

gens = generators(n)
generated = set()
bn_size = double_fact(2*n-1)
max_patience = 2000000
current_patience = max_patience

with open(f"../../BN/b{n}.pkl", "rb") as f:
    print(f"Loading BN/b{n}.pkl")
    bn = pkl.load(f)

counterexample_found = False

with tqdm() as pbar:
    while len(generated) < bn_size:
        tangle, factorization = random_factorization(n, gens)
        actual_factorization = bn[tangle.base_diagram()]
        if len(factorization) == len(actual_factorization): 
            # we want to check non minimal factorizations
            # therefore if we generate a minimal factorization, we discard it
            continue

        if tangle.base_diagram() not in generated:
            generated.add(tangle.base_diagram())

            predicted_factorization = reduce(factorization)
            predicted_factorization = [factor for factor in predicted_factorization if factor != "I"]

            predicted_len = len(predicted_factorization)
            actual_len = len(actual_factorization)

            if predicted_len != actual_len:
                counterexample_found = True
                print("Found counterexample")
                print(tangle)
                print(f"Predicted len={predicted_len} : ", predicted_factorization)
                print(f"Actual len={actual_len} : ", actual_factorization)
                break
        
            pbar.update(1)
            pbar.set_description(f"{len(generated)}/{bn_size} - {current_patience}/{max_patience}")
            current_patience = max_patience
        else:
            current_patience -= 1
            pbar.set_description(f"{len(generated)}/{bn_size} - {current_patience}/{max_patience}")

        
        if current_patience <= 0:
            break

if not counterexample_found:
    print(f"No counterexamples found in B{n} after {len(generated)} generated tangles")