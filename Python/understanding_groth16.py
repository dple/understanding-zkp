"""
This code demonstrates groth16 algorithm. Explanation of the algorithm is on this link
        https://www.rareskills.io/post/groth16

Let's reuse some code implemented in the file enc_eval_qap.py
"""

from py_ecc.bn128 import G1, G2, multiply, add, pairing, curve_order, Z1
from functools import reduce
import galois
import numpy as np
from utils import sample_Zp

print("Initializing a large field, this may take a while...")
p = curve_order
Fp = galois.GF(p)


def generate_powers_of_tau(tau, degree, point):
    return [multiply(point, int(tau ** i)) for i in range(degree + 1)]


def inner_product(ec_points, coeffs):
    # Check if number of ec_points equal to the one of coeffs
    assert len(coeffs) == len(ec_points), "Check failed"

    return reduce(add, (multiply(point, int(coeff)) for point, coeff in zip(ec_points, coeffs)), Z1)


def transform_matrix2poly(mat, wit):
    # Interpolate the matrix mat to a polynomial
    xs = Fp(np.array([1, 2, 3, 4, 5]))
    sum = galois.Poly([0, 0, 0, 0, 0], field=Fp)
    for i in range(len(wit)):
        poly = galois.lagrange_poly(xs, mat[:, i])
        sum += poly * wit[i]
    return sum


if __name__ == '__main__':
    # Create matrices from constraints
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

    witness = Fp([1, out, x, y, v1, v2, v3, v4])  # W*witness = (U*witness) * (V*witness)

    U_poly = transform_matrix2poly(L, witness)
    V_poly = transform_matrix2poly(R, witness)
    W_poly = transform_matrix2poly(O, witness)
    t_poly = galois.Poly.Roots([1, 2, 3, 4, 5], field=Fp)
    LR_product = U_poly * V_poly
    h_poly = (LR_product - W_poly) // t_poly


    """
        Trusted setup
    """
    # 1. Get a random values
    tau = Fp(sample_Zp(p))
    alpha = sample_Zp(p)
    beta = sample_Zp(p)

    alpha_G1 = multiply(G1, int(alpha))
    alpha_V_poly = V_poly * alpha
    beta_G1 = multiply(G1, int(beta))
    beta_G2 = multiply(G2, int(beta))
    beta_U_poly = U_poly * beta

    # 2. We need compute scalar multiplications for both G1 and G2
    d = len(L[:, 0]) - 1  # degree of polynomial = number of rows - 1 in matrices.
    # As each column will be interpolated to a polynomial of degree
    # Calculate tau*G1, tauˆ2*G1, ..., tauˆd*G1
    powers_of_tau_G1 = generate_powers_of_tau(tau, d, G1)
    # Calculate tau*G2, tauˆ2*G2, ..., tauˆd*G2
    powers_of_tau_G2 = generate_powers_of_tau(tau, d, G2)
    # Calculate tau*t(tau)*G1, tauˆ2*t(tau)*G1, ..., tauˆd*t(tau)*G1. Required to evaluate h(x)
    t_tau_G1 = multiply(G1, int(t_poly(tau)))
    powers_of_t_tau_G1 = generate_powers_of_tau(tau, d - 1, t_tau_G1)



    """
    U(x)*V(x) = W(x) + t(x)*h(x), thus when evaluated at tau, the following equation holds
                    U(tau)*V(tau) = W(tau) + t(tau)*h(tau),
    The straightforward to obtain U(tau), V(tau), ... is to evaluate those polynomials at tau,
    but remember that tau is generated in a trusted setup and must be kept secret. Only powers 
    of tau [tau]G1, [tau]G2, [tauˆ2]G1, [tauˆ2]G2, ..., are available to the prover. Given those 
    values, no one could recover tau unless he could solve the discrete logarithm problem.   

    Verifier will check that relation by checking if:
            pairing(U(tau)*G1, V(tau)*G2) = pairing(C*G1, G2), 
    where C = W(tau) + t(tau)*h(tau)
    """

    # Given powers of tau from the trusted setup phase, compute U(tau)*G1, V(tau)*G2, and C*G1
    evaluate_U_on_G1 = inner_product(powers_of_tau_G1, U_poly.coeffs[::-1])
    evaluate_V_on_G2 = inner_product(powers_of_tau_G2, V_poly.coeffs[::-1])
    """
        Secret shifting
        Introducing alpha and beta to prevent forgery
        U'(x) = alpha + U(x), V'(x) = beta + V(x)
    """
    alpha_U_G1 = add(alpha_G1, evaluate_U_on_G1)
    beta_V_G2 = add(beta_G2, evaluate_V_on_G2)

    # Compute W'= W + beta*U + alpha*V
    evaluate_W_on_G1 = inner_product(powers_of_tau_G1, W_poly.coeffs[::-1])
    evaluate_beta_U_on_G1 = inner_product(powers_of_tau_G1, beta_U_poly.coeffs[::-1])
    evaluate_alpha_V_on_G1 = inner_product(powers_of_tau_G1, alpha_V_poly.coeffs[::-1])
    W_beta_U_G1 = add(evaluate_W_on_G1, evaluate_beta_U_on_G1)
    W_beta_U_alpha_V_G1 = add(W_beta_U_G1, evaluate_alpha_V_on_G1)

    evaluate_ht_on_G1 = inner_product(powers_of_t_tau_G1, h_poly.coeffs[::-1])
    Wht_G1 = add(W_beta_U_alpha_V_G1, evaluate_ht_on_G1)

    # Verifying
    print("Proof Verification ...")
    assert pairing(beta_V_G2, alpha_U_G1) == pairing(beta_G2, alpha_G1) * pairing(G2, Wht_G1), "Failed to check pairings"