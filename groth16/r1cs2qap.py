"""
This Python code demonstrate the content discussed in the following article by Rareskills:
        R1CS to Quadratic Arithmetic Program over a Finite Field in Python
                (https://www.rareskills.io/post/r1cs-to-qap)

Proving x, y are solutions of out = 3xˆ3 - 5xˆ2yˆ2 + 7xyˆ2 - 21y + 11

The task is to evaluate:
1. Convert out <== 3xˆ3 - 5xˆ2yˆ2 + 7xyˆ2 - 21y + 11 to r1cs constraints
        --> Results are three matrices L, R, and O
2. Transform the r1cs constraints to Quadratic Arithmetic Program
        --> move from vector space (including L, R, O and witness vector) to polynomial space
        --> Results are 5 polynomials: U(x), V(x), W(x), t(x) and h(x)
3. Verify the equality of polynomials
        --> U(x)*V(x) = W(x) + t(x)*h(x)


Step 1: Get constraints

v1 = xx
v2 = yy
v3 = 3xv1
v4 = 5v1v2
out - 11 - v3 + v4 + 21y = 7xv2    -> as addition is 'free', we move them to the output side

witness = [1, out, x, y, v1, v2, v3, v4]

Step 2: From the above constraints, define the matrices for constraints

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

import galois
import numpy as np
from utils import *

# Define a finite field Fp, where p is a prime number
p = 7919        # Let's take a small prime number
Fp = galois.GF(p)


"""
This function transforms a matrix (L, R, or O) and a witness to a polynomial
@:param mat: input matrix
@:param wit: witness
@:return a transformed polynomial
"""
def transform_matrix2poly(mat, wit):
        xs = Fp(np.array([1, 2, 3, 4, 5]))
        sum = galois.Poly([0, 0, 0, 0, 0], field=Fp)
        for i in range(len(wit)):
                # Interpolate each column in matrix mat to a polynomial
                poly = galois.lagrange_poly(xs, mat[:,i])
                # Perform scalar multiplication between the obtained polynomial with the i-th witness, then sum up
                sum += poly*wit[i]
        return sum

if __name__ == '__main__':
        # Define matrices from constraints
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

        witness = Fp([1, out, x, y, v1, v2, v3, v4])  # O*witness = (L*witness) * (R*witness)

        # Check if (L*witness)*(R*witness) = (O*witness)
        #print(np.matmul(L, witness) * np.matmul(R, witness))
        #print(np.matmul(O, witness))
        assert all(np.equal(np.matmul(L, witness) * np.matmul(R, witness), np.matmul(O, witness))), "not equal"

        """
        Computing h(x), where U(x)*V(x) = W(x) + t(x)*h(x) and
        t(x) = (x - 1)(x - 2)(x - 3)(x - 4)(x - 5)      (because there are 5 rows or constraints)
        1. Computer inner product of the polynomials and the witnesses
        2. Computer the h(x) = (U(x)*V(x) - W(x)) / t(x) 
        """
        U_poly = transform_matrix2poly(L, witness)
        V_poly = transform_matrix2poly(R, witness)
        W_poly = transform_matrix2poly(O, witness)
        #print("Degree of polynomial = ", U_poly.degree)
        t_poly = galois.Poly.Roots([1, 2, 3, 4, 5], field=Fp)
        LR_product = U_poly * V_poly
        h_poly = (LR_product - W_poly) // t_poly
        #print("h(x) = ", h_poly)
        remainder = (LR_product - W_poly) % t_poly
        # Check if the remainder is equal to zero
        assert remainder == 0, "The remainder polynomial is not equal to 0!"

        # Verifier checks if (U_poly)*(V_poly) = (O_poly) + t(x)*h(x)
        # We are checking if polynomials in both sides are equal. Note that, this
        # is not SUCCINCT. In a zk-snark system, this relation will be checked by
        # evaluating polys at a random value z thanks to the Schwartz-Zippel Lemma.
        assert LR_product == W_poly + t_poly * h_poly, "Not equal!"