import numpy as np
import galois
import random

# Define a finite field
base_field = 131 #21888242871839275222246405745257275088548364400416034343698204186575808495617
Fp = galois.GF(base_field)
# 1. Prover --> prove the computation out <== x*y
# Define the matrices for constraints
O = Fp([[0, 1, 0, 0]])
L = Fp([[0, 0, 1, 0]])
R = Fp([[0, 0, 0, 1]])

# Provide the witness
x = Fp(random.randint(1, base_field))
y = Fp(random.randint(1, base_field))
print("x = ", x)
print("y = ", y)
out = x * y
witness = Fp([1, out, x, y])  # Cw = Aw*Bw
print("witness = ", witness)
result = np.dot(O, witness) == np.dot(L, witness) * np.dot(R, witness)

# Verifier
assert np.all(result), "Result contains an inequality"

# 2. Prover --> prove the computation out <== 3xˆ3 - 5xˆ2yˆ2 + 7xyˆ2 - 10y
# Constraints
"""
v1 = xx
v2 = yy
v3 = 3xv1
v4 = 5v1v2
out - v3 + v4 + 10y = 7xv2
"""
# From the above constraints, define the matrices for constraints
"""
witness = [1, out, x, y, v1, v2, v3, v4]
    | 1 out x  y  v1 v2 v3 v4 |  
    | 0  0  1  0  0  0  0  0  |
L = | 0  0  0  1  0  0  0  0  |
    | 0  0  3  0  0  0  0  0  |
    | 0  0  0  0  5  0  0  0  |
    | 0  0  7  0  0  0  0  0  | 
    
    | 1 out x  y  v1 v2 v3 v4 |  
    | 0  0  1  0  0  0  0  0  |
R = | 0  0  0  1  0  0  0  0  |
    | 0  0  0  0  1  0  0  0  |
    | 0  0  0  0  0  1  0  0  |
    | 0  0  0  0  0  1  0  0  | 
    
    | 1 out x  y  v1 v2 v3 v4 |  
    | 0  0  0  0  1  0  0  0  |
O = | 0  0  0  0  0  1  0  0  |
    | 0  0  0  0  0  0  1  0  |
    | 0  0  0  0  0  0  0  1  |
    | 0  1  0 10  0  0 -1  1  |            
"""

O = Fp([[0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1],
        [0, 1, 0, 10, 0, 0, Fp(base_field - 1), 1]])

L = Fp([[0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 3, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 5, 0, 0, 0],
        [0, 0, 7, 0, 0, 0, 0, 0]])
R = Fp([[0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0]])

# Provide the witness
v1 = x*x
v2 = y*y
v3 = 3*x*v1
v4 = 5*v1*v2
out = v3 - v4 - 10*y + 7*x*v2

witness = Fp([1, out, x, y, v1, v2, v3, v4])  # Cw = Aw*Bw
print("witness =", witness)

result = O.dot(witness) == L.dot(witness) * R.dot(witness)
# Verifier
assert result.all(), "Result contains an inequality 2"
