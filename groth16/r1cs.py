import numpy as np
import galois
import random

# Define a finite field Fp, where p is a prime number
p = 1327        # Let's take a small prime number
Fp = galois.GF(p)
# 1. Prover --> prove the computation out <== x*y
# Define the matrices for constraints
O = Fp([[0, 1, 0, 0]])
L = Fp([[0, 0, 1, 0]])
R = Fp([[0, 0, 0, 1]])

# Provide the witness
x = Fp(random.randint(1, p))
y = Fp(random.randint(1, p))
print("x = ", x)
print("y = ", y)
out = x * y
witness = Fp([1, out, x, y])  # Cw = Aw*Bw
print("witness = ", witness)
result = np.dot(O, witness) == np.dot(L, witness) * np.dot(R, witness)

# Verifier -> do nothing rather than very the computation result
assert np.all(result), "Result contains an inequality"

# 2. Prover --> prove the computation out <== 3xˆ3 - 5xˆ2yˆ2 + 7xyˆ2 - 21y + 11
# Constraints
"""
v1 = xx
v2 = yy
v3 = 3xv1
v4 = 5v1v2
out - 11 - v3 + v4 + 21y = 7xv2    -> as addition is 'free', we move them to the output side     
"""
# From the above constraints, define the matrices for constraints
"""
witness = [1, out, x, y, v1, v2, v3, v4]

L: Construct matrix L from left hand terms of constraints
    | 1 out x  y  v1 v2 v3 v4 |         
    | 0  0  1  0  0  0  0  0  |         -> left_hand_side x in v1 = xx
L = | 0  0  0  1  0  0  0  0  |         -> left_hand_side y in v2 = yy    
    | 0  0  3  0  0  0  0  0  |         -> left_hand_side 3x in v3 = 3xv1
    | 0  0  0  0  5  0  0  0  |         -> left_hand_side 5v1 in v4 = 5v1v2
    | 0  0  7  0  0  0  0  0  |         -> left_hand_side 7x in out - v3 + v4 + 10y = 7xv2
    
R: the matrix represents right hand side variables in constraints
    | 1 out x  y  v1 v2 v3 v4 |  
    | 0  0  1  0  0  0  0  0  |
R = | 0  0  0  1  0  0  0  0  |
    | 0  0  0  0  1  0  0  0  |
    | 0  0  0  0  0  1  0  0  |
    | 0  0  0  0  0  1  0  0  | 

O: the matrix represents output variables in constraints
    |  1 out x  y  v1 v2 v3 v4 |  
    |  0  0  0  0  1  0  0  0  |
O = |  0  0  0  0  0  1  0  0  |
    |  0  0  0  0  0  0  1  0  |
    |  0  0  0  0  0  0  0  1  |
    |-11  1  0 21  0  0 -1  1  |         -> constant -11 is under the term 1   
"""

O = Fp([[0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 1],
        [Fp(p - 11), 1, 0, 21, 0, 0, Fp(p - 1), 1]])    # in a finite field Fp, -1 = p - 1, and -11 = p - 11
                                                        # as p - 1 = -1 mod p and p - 11 = -11 mod p

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
out = v3 - v4 - 21*y + 7*x*v2 + Fp(11)

witness = Fp([1, out, x, y, v1, v2, v3, v4])  # Cw = Aw*Bw
print("witness =", witness)

result = O.dot(witness) == L.dot(witness) * R.dot(witness)

# Verifier  -> do nothing rather than very the computation result
assert result.all(), "Result contains an inequality 2"
