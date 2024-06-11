"""
This Python code was consolidated from the following lecture by Rareskills
        Encrypted Evaluation of a Quadratic Arithmetic Program
          (https://www.rareskills.io/post/elliptic-curve-qap)

The task is to evaluate out <== 3xˆ3 - 5xˆ2yˆ2 + 7xyˆ2 - 21y + 11

- Constraints

v1 = xx
v2 = yy
v3 = 3xv1
v4 = 5v1v2
out - 11 - v3 + v4 + 21y = 7xv2    -> as addition is 'free', we move them to the output side

witness = [1, out, x, y, v1, v2, v3, v4]

- From the above constraints, define the matrices for constraints

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

from py_ecc.bn128 import G1, G2, multiply, add, pairing, curve_order, Z1
from functools import reduce
import galois
import numpy as np
from utils import sample_Zp

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
    print("Initializing a large field, this may take a while...")
    p = curve_order
    Fp = galois.GF(p)

    print("Initialization completed. Starting computation ...")

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

    """
    Computing h(x), such that U(x)*V(x) = W(x) + t(x)*h(x), where
        t(x) = (x - 1)(x - 2)(x - 3)(x - 4)(x - 5) as there are 5 rows
        h(x) = (U(x)*V(x) - W(x)) / t(x) 
    """

    U_poly = transform_matrix2poly(L, witness)
    V_poly = transform_matrix2poly(R, witness)
    W_poly = transform_matrix2poly(O, witness)
    t_poly = galois.Poly.Roots([1, 2, 3, 4, 5], field=Fp)
    LR_product = U_poly * V_poly
    h_poly = (LR_product - W_poly) // t_poly
    remainder = (LR_product - W_poly) % t_poly

    assert remainder == 0, "The remainder polynomial is not equal to zero!"

    """
        Trusted setup phase, calculating ([tau^0*G], [tau^1*G], [tau^2*G], ..., [tau^d*G]), 
        where G could be G1 or G2. Recall that G1 is the group of points defined on the base field Fp,
        and G2 is the group of points defined on the extension field (i.e, Fpˆ12 in the BN curves), 
        thus calculations on G1 are much cheaper than ones on G2. 
        
    """
    # 1. Get a random value tau
    tau = Fp(sample_Zp(p))

    # 2. We need compute scalar multiplications for both G1 and G2
    d = len(L[:, 0]) - 1        # degree of polynomial = number of rows - 1 in matrices.
                                # As each column will be interpolated to a polynomial of degree
    # Calculate tau*G1, tauˆ2*G1, ..., tauˆd*G1
    powers_of_tau_G1 = generate_powers_of_tau(tau, d, G1)
    # Calculate tau*G2, tauˆ2*G2, ..., tauˆd*G2
    powers_of_tau_G2 = generate_powers_of_tau(tau, d, G2)

    # Calculate tau*t(tau)*G1, tauˆ2*t(tau)*G1, ..., tauˆd*t(tau)*G1. Required to evaluate h(x)
    # For this per circuit trusted setup, prover must send the polynomial t(x) to the trusted server
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
    where C*G1 = (W(tau) + t(tau)*h(tau))*G1, and
    - U(tau)*G1 = sum(u_i*[tauˆi*G1]), where u_i are coefficients of U(x)
    - V(tau)*G2 = sum(v_i*[tauˆi*G2]), where v_i are coefficients of V(x)
    - W(tau)*G1 = sum(w_i*[tauˆi*G2]), where w_i are coefficients of W(x)
    - t(tau)*h(tau)*G1 = sum(h_i*[tauˆî*t(tau)*G1]), where h_i are coefficients of h(x) 
        and tauî*t(tau)*G1 are precomputed in the 2nd phase of the trusted setup. 
    """

    # Given powers of tau from the trusted setup phase, compute U(tau)*G1, V(tau)*G2, and C*G1
    evaluate_U_on_G1 = inner_product(powers_of_tau_G1, U_poly.coeffs[::-1])
    evaluate_V_on_G2 = inner_product(powers_of_tau_G2, V_poly.coeffs[::-1])
    evaluate_W_on_G1 = inner_product(powers_of_tau_G1, W_poly.coeffs[::-1])
    evaluate_ht_on_G1 = inner_product(powers_of_t_tau_G1, h_poly.coeffs[::-1])
    Wht_G1 = add(evaluate_W_on_G1, evaluate_ht_on_G1)

    """
    Verifying. Prover sends (or publishes) three point elements to the veririfer, including:
    - evaluate_V_on_G2 = U(tau)*G1
    - evaluate_I_on_G1 = V(tau)*G2
    - Wht_G1 = (W(tau) + t(tau)*h(tau))*G1
    The verifier will be required to compute two cryptographic pairings, then compare their results. 
    - Accept the proof if the two pairings returns the same result
    - Otherwise, reject the proof
    """

    print("Proof Verification ...")
    assert pairing(evaluate_V_on_G2, evaluate_U_on_G1) == pairing(G2, Wht_G1),"Failed to check pairings"