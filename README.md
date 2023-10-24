# BrauerMonoidFactorization
Example of usage:
```python
from brauermonoid import *

X = text_to_tangle("1:4,2:4',3:5,6:1',2':3',5':6'")
tau_X = tau(X)
print(factorizeSN(tau_X))
print(factorizeBN(X))
```
