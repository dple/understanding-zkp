"""
This Python code demonstrate the content discussed in the following article by Rareskills:
R1CS to Quadratic Arithmetic Program over a Finite Field in Python
(https://www.rareskills.io/post/r1cs-to-qap)
"""
import galois
import numpy as np
from utils import *

# Define a finite field Fp, where p is a prime number
p = 7919        # Let's take a small prime number
Fp = galois.GF(p)

"""
1. Convert 3xˆ3 - 5xˆ2yˆ2 + 7xyˆ2 - 21y + 11 to r1cs constraints
2. Convert r1cs constraints to Quadratic Arithmetic Program
"""
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
    
Note: In a finite field Fp, -1 = p - 1, and -11 = p - 11 as p - 1 = -1 mod p and p - 11 = -11 mod p
"""

def transform_matrix2poly(mat, wit):
        # Interpolate the matrix mat to a polynomial
        xs = Fp(np.array([1, 2, 3, 4, 5]))
        sum = galois.Poly([0, 0, 0, 0, 0], field=Fp)
        for i in range(len(wit)):
                poly = galois.lagrange_poly(xs, mat[:,i])
                sum += poly*wit[i]
        return sum

if __name__ == '__main__':
        O = Fp([[0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0, 1],
                [Fp(p - 11), 1, 0, 21, 0, 0, Fp(p - 1), 1]])

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

        # Compute the witness
        x = Fp(sample_Zp(p))
        y = Fp(sample_Zp(p))
        v1 = x * x
        v2 = y * y
        v3 = 3 * x * v1
        v4 = 5 * v1 * v2
        out = v3 - v4 - 21 * y + 7 * x * v2 + Fp(11)

        witness = Fp([1, out, x, y, v1, v2, v3, v4])  # Cw = Aw*Bw

        assert all(np.equal(np.matmul(L, witness) * np.matmul(R, witness), np.matmul(O, witness))), "not equal"

        """
        Computing h(x), where U(x)*V(x) = W(x) + t(x)*h(x)
        t(x) = (x - 1)(x - 2)(x - 3)(x - 4)(x - 5) as there are 5 rows
        1. Computer inner (Hadamard) product of the polynomials and the witnesses
        2. Computer the h(x) = (U(x)*V(x) - W(x)) / t(x) 
        """
        U_poly = transform_matrix2poly(L, witness)
        V_poly = transform_matrix2poly(R, witness)
        W_poly = transform_matrix2poly(O, witness)
        t_poly = galois.Poly.Roots([1, 2, 3, 4, 5], field=Fp)
        LR_product = U_poly * V_poly
        h_poly = (LR_product - W_poly) // t_poly
        remainder = (LR_product - W_poly) % t_poly
        # Check if the remainder is equal to zero
        assert remainder == 0, "The remainder polynomial is not equal to 0!"

        # Verifier check if two sides are balanced
        assert LR_product == W_poly + t_poly * h_poly, "Not equal!"