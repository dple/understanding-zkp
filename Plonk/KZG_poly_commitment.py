"""
    Python code demonstrating KZG polynomial commitment scheme
    Assume prover need to commit the polynomial f(X) = 5X^4 - 2X +3
"""
from py_ecc.bn128 import G1, G2, multiply, add, neg, curve_order, Z1, pairing
import galois
from functools import reduce
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
    print("tau: ", tau)
    d = 4                           # here we are working with a poly of degree 4
    powers_of_tau_G1 = generate_powers_of_tau(tau, d, G1)   # Generate powers of tau [G1, [tau]G1, [tau^2]G1, [tau^3]G1, [tau^4]G1]
    tauG2 = multiply(G2, int(tau))
    
    print("Committing poly f(X) ...")
    fX = galois.Poly([5, 0, 0, - 2, 3], field = GF)
    com_f = inner_product(powers_of_tau_G1, fX.coeffs[::-1])
    print("Commitment com_f: ", com_f)

    # Open the commitment
    print("Opening the commitment ...")    
    u = random.randint(1, p)               # Verifier chooses a random u, then sends it the prover
    # Prover evaluates the polynomial fX at u, and calculates the polynomial Q(X) = (f(X) - f(u)) / (x - u)
    v = fX(u)           
    tX = galois.Poly([1, -u], field = GF)
    qX = (fX - v) // tX
    remainder = (fX - v) % tX
    assert remainder == 0, "The remainder polynomial is not equal to zero!"

    # Prover generates a proof of commitment com_f for the random challenge u, then sends it to the verifier
    print("Generating proof ...")
    proof = inner_product(powers_of_tau_G1[:d], qX.coeffs[::-1])
    print("Proof com_q: ", proof)

    # Verifier computes the below values and verify the proof using pairings
    uG2 = multiply(G2, int(u))
    tau_minus_uG2 = add(tauG2, neg(uG2))
    vG1 = multiply(G1, int(v))
    com_f_minus_vG1 = add(com_f, neg(vG1))

    if pairing(tau_minus_uG2, proof) == pairing(G2, com_f_minus_vG1):
        print("Proof for commitment is correct!")
    else:
        print("Failed to test the proof!")