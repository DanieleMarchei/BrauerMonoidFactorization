#all codes in here were adapted from sage.combinat.diagram_algebras
#because this allows the usage of PyPy without the need of installing sage on it

def perfect_matchings_iterator(n):
    if n == 0:
        yield []
        return

    i = 0
    x = 0
    y = 0
    g = 0
    j = 0
    J = 0
    e = [i for i in range(2*n)]
    
    f = [0 for _ in range(2*n)]
    for i in range(2*n):
        if i % 2 == 0:
            f[i] = i + 1
        else:
            f[i] = i - 1

    odd = False

    yield convert(f, n)
    while e[0] != n - 1:
        i = e[0]
        if odd:
            x = 2 * i
        else:
            x = i

        y = f[x]
        g = y - x - 1
        if g % 2 == odd:
            g += 1
            j = y + 1
        else:
            g -= 1
            j = y-1
        J = f[j]
        f[y] = J
        f[J] = y
        f[x] = j
        f[j] = x
        odd = not odd
        e[0] = 0
        if g == 0 or g == 2 * (n-i-1):
            e[i] = e[i+1]
            e[i+1] = i + 1

        yield convert(f, n)


def convert(f,n):
    """
    Convert a list ``f`` representing a fixed-point free involution
    to a set partition.
    """
    ret = []
    i = 0
    for i in range(2*n):
        if i < f[i]:
            ret.append((i, f[i]))
    return ret

def brauer_diagrams(k):
    r"""
    from sage.combinat.diagram_algebras
    """
    s = list(range(1,k+1)) + list(range(-k,0))
    for p in perfect_matchings_iterator(k):
        b = [(s[a],s[b]) for a,b in p]
        yield b
