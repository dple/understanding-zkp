"""
This code demonstrates groth16 algorithm. Explanation of the algorithm is on this link
        https://www.rareskills.io/post/groth16

Let's reuse some code implemented in the file enc_eval_qap.py.
You can get more explanation on why these calculations on the above source file.

The task is to evaluate out <== 3xˆ3 - 5xˆ2yˆ2 + 7xyˆ2 - 21y + 11

- Constraints

v1 = xx
v2 = yy
v3 = 3xv1
v4 = 5v1v2
out - 11 - v3 + v4 + 21y = 7xv2    -> as addition is 'free', we move them to the output side

witness = [1, out, x, y, v1, v2, v3, v4]

"""

from py_ecc.bn128 import G1, G2, multiply, add, neg, pairing, curve_order, Z1
from functools import reduce
import galois
import numpy as np
import secrets


def sample_Zp(p):
    return secrets.randbelow(p)

def inner_product(ec_points, coeffs):
    # Check if number of ec_points equal to the one of coeffs
    assert len(coeffs) == len(ec_points), "Check failed"

    return reduce(add, (multiply(point, int(coeff)) for point, coeff in zip(ec_points, coeffs)), Z1)

"""
Generate powers of tau used to evaluate polynomials U(x), V(x), W(x) and h(x)
"""
def generate_powers_of_tau(tau, degree, point):
    return [multiply(point, int(tau ** i)) for i in range(degree + 1)]

"""
Given three matrices L, R and O, generate powers of tau for public and private input, 
i,e., C = beta*U + alpha*V + W
    @:param d: degree of the polynomials formed from columns in the matrices, d = n - 1, where n = #rows
    @:param m: #columns of the matrices = len(witness)

    @:return two list of powers of tau (elliptic curve points):
     - one related to public inputs (i.e., first l values in the witness vector) will be for verifier
     - one related to private inputs (i.e, last m - l values in the witness vector) will be for prover
"""
def generate_powers_of_tau_for_inputs(powers_of_tau, d, m, L_mat, R_mat, O_mat, alpha, beta, ell, gamma_inv, delta_inv):
    taus_for_public_inputs = []
    taus_for_private_inputs = []
    # require degree + 1 points to interpolate a polynomial of degree d
    xs = Fp(np.array([i + 1 for i in range(d + 1)]))
    # Each i-th col of matrices L, R, and W will be converted to a polynomial U_i(x), V_i(x) and W_i(x)
    for i in range(m):
        # U_i(x) = interpolate the i-th column of matrix L
        poly_Ui = galois.lagrange_poly(xs, L_mat[:, i])
        # Perform a random shift by multiplying the poly with a random factor beta
        beta_Ui = poly_Ui * beta        # multiply U with beta
        # V_i(x) = interpolate the i-th column of matrix R
        poly_Vi = galois.lagrange_poly(xs, R_mat[:, i])
        # Perform a random shift by multiplying the poly with a random factor alpha
        alpha_Vi = poly_Vi * alpha
        # W_i(x) = interpolate the i-th column of matrix W
        poly_Wi = galois.lagrange_poly(xs, O_mat[:, i])

        sum_poly = beta_Ui + alpha_Vi + poly_Wi

        if i < ell:
            taus_for_public_inputs.append(inner_product(powers_of_tau, (sum_poly.coeffs[::-1]) * gamma_inv))
        else:
            taus_for_private_inputs.append(inner_product(powers_of_tau, (sum_poly.coeffs[::-1]) * delta_inv))

    return taus_for_public_inputs, taus_for_private_inputs

def transform_matrix2poly(mat, wit, d):
    # Interpolate the matrix mat to a polynomial
    xs = Fp(np.array([i + 1 for i in range(d + 1)]))    #xs = Fp(np.array([1, 2, 3, 4, 5]))
    sum = galois.Poly([0] * (d + 1), field=Fp)
    for i in range(len(wit)):
        poly = galois.lagrange_poly(xs, mat[:, i])
        sum += poly * wit[i]
    return sum


if __name__ == '__main__':
    print("Initializing a large field, this may take a while...")
    p = curve_order
    Fp = galois.GF(p)
    print("Initialization completed. Start computation ...")

    """
        Prover's computations
    """
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

    n = len(L[:, 0])        # number of rows
    d = n - 1               # degree polynomials interpolated from matrices' columns
    m = len(L[0, :])        # number of columns = len(witness)
    # Compute the witness
    x = Fp(sample_Zp(p))
    y = Fp(sample_Zp(p))
    v1 = x * x
    v2 = y * y
    v3 = 3 * x * v1
    v4 = 5 * v1 * v2
    out = v3 - v4 - 21 * y + 7 * x * v2 + Fp(11)

    witness = Fp([1, out, x, y, v1, v2, v3, v4])  # W*witness = (U*witness) * (V*witness)

    # U(x) = sum_{i = 0}^m (w_i u_i(x)), where w_i is i-th element of the witness,
    # u_i(x) is the poly interpolated from the i-th column of matrix L, m is no of elements of the witness
    U_poly = transform_matrix2poly(L, witness, d)
    V_poly = transform_matrix2poly(R, witness, d)   # V(x) = sum(w_i v_i(x))
    W_poly = transform_matrix2poly(O, witness, d)   # W(x) = sum(w_i w_i(x))
    # t(x) = (x - 1)(x - 2)...(x - (d + 1))
    t_poly = galois.Poly.Roots(np.array([i + 1 for i in range(d + 1)]), field=Fp)
    LR_product = U_poly * V_poly
    h_poly = (LR_product - W_poly) // t_poly


    """
        Trusted setup
        Phase 1: setup for all circuits
    """
    # 1. Get a random value tau
    tau = sample_Zp(p)
    # 2. Calculate tau*G1, tauˆ2*G1, ..., tauˆd*G1
    powers_of_tau_for_G1 = generate_powers_of_tau(tau, d, G1)
    # Calculate tau*G2, tauˆ2*G2, ..., tauˆd*G2
    powers_of_tau_for_G2 = generate_powers_of_tau(tau, d, G2)
    # Calculate powers of tau for evaluating h(x)t(x) at tau
    # Calculate tau*t(tau)*G1, tauˆ2*t(tau)*G1, ..., tauˆ{d - 1}*t(tau)*G1
    t_tau = t_poly(tau)
    powers_of_tau_for_ht = [multiply(powers_of_tau_for_G1[i], int(t_tau)) for i in range(d)]


    """
        Phase 2: trusted setup per circuit. This setup requires for each individual circuit
    """
    # Introduce random alpha, beta to prevent a potential cheating prover from making up
    # three curve points of proof A, B, and C
    alpha = sample_Zp(p)
    beta = sample_Zp(p)
    alpha_G1 = multiply(G1, int(alpha))     # alpha*G1
    beta_G1 = multiply(G1, int(beta))       # beta*G1
    beta_G2 = multiply(G2, int(beta))       # beta*G2


    """
    Prover: compute
    U(x)*V(x) = W(x) + t(x)*h(x), thus when evaluated at tau, the following equation holds
                    U(tau)*V(tau) = W(tau) + t(tau)*h(tau),
    The straightforward to obtain U(tau), V(tau), ... is to evaluate those polynomials at tau,
    but remember that tau is generated in a trusted setup and must be kept secret. Only powers 
    of tau [tau]G1, [tau]G2, [tauˆ2]G1, [tauˆ2]G2, ..., are available to the prover. Given those 
    values, no one could recover tau unless he could solve the discrete logarithm problem.   

    Verifier will check that relation by checking if:
            pairing(B, A) = pairing(G2, C), 
    where A = U(tau)*G1, B = V(tau)*G2, C = (W(tau) + t(tau)*h(tau))*G1
    """

    # Given powers of tau from the trusted setup phase, prover compute U(tau)*G1, V(tau)*G2, and C*G1
    A = inner_product(powers_of_tau_for_G1, U_poly.coeffs[::-1])        # U(tau)*G1
    B1 = inner_product(powers_of_tau_for_G1, V_poly.coeffs[::-1])       # V(tau)*G1
    B2 = inner_product(powers_of_tau_for_G2, V_poly.coeffs[::-1])       # V(tau)*G2
    """
        Secret shifting
        Introducing alpha and beta to prevent a malicious prover to make up values U(tau)*G1, V(tau)*G2 and C
        U'(x) = alpha + U(x), V'(x) = beta + V(x)
    """
    alpha_A = add(alpha_G1, A)          # random shift for A, [alpha + U(tau)]*G1
    beta_B1 = add(beta_G1, B1)          # random shift for B, [beta + V(tau)]*G1
    beta_B2 = add(beta_G2, B2)          # random shift for B, [beta + V(tau)]*G2

    # Check #1
    evaluate_ht_on_G1 = inner_product(powers_of_tau_for_ht, h_poly.coeffs[::-1])

    _, taus_for_C = generate_powers_of_tau_for_inputs(powers_of_tau_for_G1, d, m, L, R, O, alpha, beta, 0, 1, 1)
    C_taus = inner_product(taus_for_C, witness)
    C = add(C_taus, evaluate_ht_on_G1)

    print("Proof Verification ...")
    if pairing(beta_B2, alpha_A) == pairing(beta_G2, alpha_G1) * pairing(G2, C):
        print("Pass test #1, proof is correct after introducing alpha and beta!")
    else:
        print("Failed test #1")
        exit(1)

    """
    Separate public and private inputs with gamma and delta
    """
    ell = 2  # number of public inputs. In this example, public inputs are 1 and out
    taus_for_public_inputs, taus_for_private_inputs = \
        generate_powers_of_tau_for_inputs(powers_of_tau_for_G1, d, m, L, R, O, alpha, beta, ell, 1, 1)
    C_public = inner_product(taus_for_public_inputs, witness[:ell])
    C_private_taus = inner_product(taus_for_private_inputs, witness[ell:])
    C_private = add(C_private_taus, evaluate_ht_on_G1)

    # Check #2
    if pairing(beta_B2, alpha_A) == pairing(beta_G2, alpha_G1) * pairing(G2, C_private) * pairing(G2, C_public):
        print("Pass test #2, proof is correct after separating public and private inputs!")
    else:
        print("Failed test #2")
        exit(1)

    """
    Introducing gamma and delta to prevent forgeries with public inputs 
    """
    gamma = sample_Zp(p)
    gamma_G2 = multiply(G2, int(gamma))  # gamma*G2
    gamma_inv = Fp(1) / Fp(gamma)
    delta = sample_Zp(p)
    delta_G1 = multiply(G1, int(delta))  # delta*G1
    delta_G2 = multiply(G2, int(delta))  # delta*G2
    delta_inv = Fp(1) / Fp(delta)

    # Calculate powers of tau for evaluating (h(x)t(x))/delta at tau
    #      delta^{-1}*tau*t(tau)*G1, delta^{-1}*tauˆ2*t(tau)*G1, ..., delta^{-1}*tauˆ{d - 1}*t(tau)*G1
    t_tau_delta = Fp(t_poly(tau) * delta_inv)
    powers_of_tau_for_ht = [multiply(powers_of_tau_for_G1[i], int(t_tau_delta)) for i in range(d)]
    evaluate_ht_on_G1 = inner_product(powers_of_tau_for_ht, h_poly.coeffs[::-1])

    # Calculate alpha*V(tau) + beta*U(tau) + W(au)
    taus_for_public_inputs, taus_for_private_inputs = \
        generate_powers_of_tau_for_inputs(powers_of_tau_for_G1, d, m, L, R, O, alpha, beta, ell, gamma_inv, delta_inv)

    C_public = inner_product(taus_for_public_inputs, witness[:ell])
    C_private_taus = inner_product(taus_for_private_inputs, witness[ell:])
    C_private = add(C_private_taus, evaluate_ht_on_G1)

    # Check #3
    if pairing(beta_B2, alpha_A) == pairing(beta_G2, alpha_G1) * \
            pairing(delta_G2, C_private) * pairing(gamma_G2, C_public):
        print("Pass test #3, proof is correct after introducing gamma and delta to prevent forgeries with public inputs!")
    else:
        print("Failed test #3")
        exit(1)


    """
    Enforcing true Zero-Knowledge by introducing two random values r and s
    Prover samples random value r, and s to prevent verifier from guessing witness values
    """
    r = sample_Zp(p)
    s = sample_Zp(p)
    r_delta_G1 = multiply(delta_G1, int(r))
    s_delta_G1 = multiply(delta_G1, int(s))
    s_delta_G2 = multiply(delta_G2, int(s))
    rs_delta_G1 = multiply(r_delta_G1, int(s))

    # Update A, B
    A = add(alpha_A, r_delta_G1);
    B1 = add(beta_B1, s_delta_G1);
    B2 = add(beta_B2, s_delta_G2);
    sA = multiply(A, int(s))
    rB1 = multiply(B1, int(r))


    """
        Compute C = (W(x) + beta*U(x) + alpha*V(x))*w + h(x)t(x) + sA + rB1 + rs*delta*G1
    """
    C_ht_sa = add(C_private, sA);
    C_ht_sa_rb = add(C_ht_sa, rB1);
    C = add(C_ht_sa_rb, neg(rs_delta_G1))

    """
    Verifying: verifier obtains a proof from prover, consisting of three elliptic curve points:
        - [A, B2, C]
        and two points from trusted server:
        - alpha_G1
        - beta_G2
        - C_public 
        - gamma_G2, 
        - delta_G2
        
    Accept the proof if the below equation is true. Otherwise, reject

    Note that: pairing(beta_G2, alpha_G1) could be pre-computed to make the proof verifications faster
    """
    # Final check
    if pairing(B2, A) == pairing(beta_G2, alpha_G1) * \
            pairing(delta_G2, C) * pairing(gamma_G2, C_public):
        print("Pass final test after adding two random values (r, s) to ensure the zero-knowledge!")
    else:
        print("Failed final test")
        exit(1)
