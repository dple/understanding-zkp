"""
This python code consolidated and slightly modified the example shown by Rareskills
                Quadratic Arithmetic Programs
https://www.rareskills.io/post/quadratic-arithmetic-program

Instead of working with integers, I adapted the code in a finite field GF(p).

The equation under evaluation is: out <== x1^2 + 4x2^2*x1 - 2

Constraints:
x3 = x1*x1
x4 = x2*x2
out -x3 + 2 = 4x4*x1
"""

import numpy as np
import galois
from utils import *

p = 1327
Fp = galois.GF(p)

"""
Given a matrix (L, R, or O) and a witness, transform them to a polynomial 
"""
def transform(mat, wit):
        xs = Fp(np.array([1, 2, 3]))
        sum = galois.Poly([0, 0, 0], field=Fp)
        for i in range(len(wit)):
                poly = galois.lagrange_poly(xs, mat[:,i])
                sum += poly*wit[i]
        return sum

if __name__ == '__main__':
    L = Fp([[0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 4]])

    R = Fp([[0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 1, 0, 0, 0]])

    O = Fp([[0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
            [2, 1, 0, 0, Fp(p - 1), 0]])    # -1 = p - 1 in GF(p)

    # withness
    x1 = Fp(3) #Fp(sample_Zqstar(p))
    x2 = Fp(4) #Fp(sample_Zqstar(p))
    x3 = x1 * x1
    x4 = x2 * x2
    out = 4 * x4 * x1 + x3 - Fp(2)

    witness = Fp(np.array([1, out, x1, x2, x3, x4]))
    print(witness)

    L_poly = transform(L, witness)
    print(L_poly)
    R_poly = transform(R, witness)
    print(R_poly)
    O_poly = transform(O, witness)
    print(O_poly)
    t_poly = galois.Poly.Roots([1, 2, 3], field=Fp)
    print(t_poly)
    LR_product = L_poly*R_poly
    print(LR_product)
    h_poly = (LR_product - O_poly) // t_poly
    remainder = (LR_product - O_poly) % t_poly

    print(h_poly)
    print("Remainder = ", remainder)    # remainder must be Zero
    assert remainder == 0, "The remainder polynomial is not equal to zero!"

    assert LR_product == O_poly + t_poly*h_poly, "Not equal!"