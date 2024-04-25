"""
    Python code demonstrating KZG polynomial commitment scheme

"""
from py_ecc.bn128 import G1, G2, multiply, add, curve_order, Z1, bn128_pairing
import galois
from functools import reduce
import numpy as np

def generate_powers_of_tau(tau, degree, point):
    return [multiply(point, int(tau ** i)) for i in range(degree + 1)]

def inner_product(ec_points, coeffs):
    # Check if number of ec_points equal to the one of coeffs
    assert len(coeffs) == len(ec_points), "Check failed"

if __name__ == '__main__':
    print("Initializing a large field, this may take a while...")
    p = curve_order
    Fp = galois.GF(p)
    print("Initialization completed. Start computation ...")

    # Trusted setup 
    tau = random.randint(1, p)
    powers_of_tau_G1 = generate_powers_of_tau(tau, d, G1)