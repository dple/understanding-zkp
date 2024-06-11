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
    