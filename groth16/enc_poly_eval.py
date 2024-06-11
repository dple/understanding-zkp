"""
This Python code was copied from the following lecture by Rareskills
                Encrypted Polynomial Evaluation
(https://www.rareskills.io/post/encrypted-polynomial-evaluation)

Given a polynomial p(x) = a_0 + a_1x + a_2x^2 + ... + a_dx^d, evaluate p(x) at a value t is to compute:

p(t) = a_0 + a_1 * t + a_2 * t^2 + ... + a_d * t^d

However, in ZKP system, we are working with a group of points on defined elliptic curve E over a finite field Fp,
that is, we will convert

a_0 -> a_0 * G
a_1 * t -> [a_1 * t] * G
....
a_d * t^d -> [a_d * t^d] * G

where G is the generator of the group of points we are working on.
Those conversions require scalar multiplications of points on the elliptic curve E
And, evaluating a polynomial will be converted as:

p(t) * G = a_0 * G + [a_1 * t] * G + [a_2 * t^2] * G + ... + [a_d * t^d] * G

Calculating scalar multiplication is quite expensive, especially when the degree of
polynomial p(x), d, is getting big. Moreover, the proof is not succinct as its size
is linear to d. That's why we need cryptographic pairings e: G1 x G2 -> G12 that won't
be implemented in this file.
"""


from py_ecc.bn128 import G1, multiply, add, eq, curve_order, Z1
import galois
from functools import reduce
from utils import sample_Zp

print("Initializing a large field, this may take a while...")
Fp = galois.GF(curve_order)

def generate_powers_of_tau(tau, degree):
    return [multiply(G1, int(tau ** i)) for i in range(degree + 1)]

def inner_product(ec_points, coeffs):
    return reduce(add, (multiply(point, int(coeff)) for point, coeff in zip(ec_points, coeffs)), Z1)

if __name__ == '__main__':

    # The factored form of a polynomial to be evaluated, p = (x - 4) * (x + 2)
    p = galois.Poly.Roots([4, -2], field=Fp)

        # Trusted setup phase,
    # 1. Get a random value tau
    tau = Fp(sample_Zp(curve_order))
    # 2. Calculate ([tau^0*G], [tau^1*G], [tau^2*G], ..., [tau^d*G])
    powers_of_tau = generate_powers_of_tau(tau, p.degree)

    # Evaluate via encrypted evaluation coefficients need to be reversed to match the powers,
    # that is, compute a_0 * G + a_1 * [tau * G] + a_2 * [tau^2 * G] + ... + a_d * [tau^d * G]
    evaluate_on_ec = inner_product(powers_of_tau, p.coeffs[::-1])

    # Evaluate p(x) at tau then convert to the point of elliptic curve p(tau)*G
    evaluate_then_convert_to_ec = multiply(G1, int(p(tau)))

    """
    Assert the equality of two above calculations
        a_0 * G + a_1 * [tau * G] + a_2 * [tau^2 * G] + ... + a_d * [tau^d * G]
        = a_0 * G + (a_1 * tau) * G + (a_2 * tau^2) * G] + ... + (a_d * tau^d) * G
        = p(tau) * G
    """
    assert eq(evaluate_then_convert_to_ec, evaluate_on_ec), "Asserting encrypted polynomial evaluation failed"