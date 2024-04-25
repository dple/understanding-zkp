"""
    Python code demonstrating KZG polynomial commitment scheme
    Assume prover need to commit the polynomial f(X) = 5X^4 - 2X +3. Then he proves evaluate it at 2, i.e., f(2) = 79 
"""
from py_ecc.bn128 import G1, G2, multiply, add, neg, curve_order, Z1, pairing
import galois
from functools import reduce
import numpy as np
import random

def generate_powers_of_tau(tau, degree, point):
    return [multiply(point, int(tau ** i)) for i in range(degree + 1)]

def inner_product(ec_points, coeffs):    
    # Check if number of ec_points equal to the one of coeffs
    assert len(coeffs) == len(ec_points), "Check failed"

    return reduce(add, (multiply(point, int(coeff)) for point, coeff in zip(ec_points, coeffs)), Z1)

if __name__ == '__main__':
    print("Initializing a large field, this may take a while...")
    p = curve_order
    GF = galois.GF(p)
    print("Initialization completed. Start computation ...")

    # Trusted setup 
    print("Starting the trusted setup process ...")
    tau = random.randint(1, p)      # get a random number tau
    d = 4                           # here we are working with a poly of degree 4
    powers_of_tau_G1 = generate_powers_of_tau(tau, d, G1)   # Generate powers of tau [G1, [tau]G1, [tau^2]G1, [tau^3]G1, [tau^4]G1]
    tauG2 = multiply(G2, int(tau))
    
    print("Committing poly f(X) ...")
    fX = galois.Poly([5, 0, 0, - 2, 3], field = GF)
    com_f = inner_product(powers_of_tau_G1, fX.coeffs[::-1])

    # Open the commitment
    print("Opening the commitment ...")    
    u = 2               # Verifier chooses a random v, e.g., v = 2, then sends it the prover
    # Prover calculates v and the polynomial Q(X)
    v = fX(u)           
    qX = (fX - v) // galois.Poly([1, -u], field = GF)

    print("Generating proof ...")
    com_q = inner_product(powers_of_tau_G1[:d], qX.coeffs[::-1])

    # Verifier
    uG2 = multiply(G2, int(u))
    tau_minus_uG2 = add(tauG2, neg(uG2))
    vG1 = multiply(G1, int(v))
    com_f_minus_vG1 = add(com_f, neg(vG1))

    if pairing(tau_minus_uG2, com_q) == pairing(G2, com_f_minus_vG1):
        print("Proof for commitment is correct!")
    else:
        print("Failed to test the proof!")