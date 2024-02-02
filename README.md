# Brauer Monoid Factorization

Code for the paper Marchei D., Merelli E., Francis A. (2024). ``Factorizing the Brauer Monoid in polynomial time'' (to be submitted).

Example of usage:
```python
from brauermonoid import *

X = text_to_tangle("1:4,2:4',3:5,6:1',2':3',5':6'")
tau_X = tau(X)
print(factorizeSN(tau_X))
print(factorizeBN(X))
```
which outputs:

```bash
['T3', 'T4', 'T5', 'T2', 'T3', 'T1', 'T2']
['T3', 'U4', 'U5', 'T2', 'T3', 'U1', 'U2']
```

If you would like to run the code for testing the Assumptions, please install [SageMath](https://www.sagemath.org/) and Maude System for python:
```bash
pip install maude
```
the code is located in the folder ```assumption_tests```. Please unzip the file ```bn.zip```, which contains the databases we mentioned in the Appendix.

If you also would like to reproduce the tables in the Discussion section, please run the code in ```number_of_merges.py``` and ```number_of_tangles_with_k_factors.py```. We ***highly*** recommend you use the [PyPy](https://www.pypy.org/) interpreter, because enumerating all tangles is a demanding task from N = 8 onward.