"""
This code demonstrates groth16 algorithm. Explanation of the algorithm is on this link
        https://www.rareskills.io/post/groth16

Let's reuse some code implemented in the file enc_eval_qap.py.
You can get more explanation on why these calculations on the above source file.
"""

from py_ecc.bn128 import G1, G2, multiply, add, pairing, curve_order, Z1
from functools import reduce
import galois
import numpy as np
from utils import sample_Zp

print("Initializing a large field, this may take a while...")
p = curve_order
Fp = galois.GF(p)


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
Generate powers of tau used to compute C = beta*U + alpha*V + W
degree: degree of the polynomials formed from columns in the matrices, degree = n - 1, where n = #rows
m: #columns of the matrices = len(witness)

Return m elliptic curve points corresponding to m elements of the witness, or #columns of matrices L, R, or O 
"""
def generate_powers_of_tau_W(tau, d, m, point, L_mat, R_mat, O_mat, alpha, beta):
    list_taus = []
    # require degree + 1 points to interpolate a polynomial of degree d
    xs = Fp(np.array([i + 1 for i in range(d + 1)]))
    # Each col of matrices L, R, and W will be converted to a polynomial
    for i in range(m):
        # U_i(x) = interpolate the i-th column of matrix L
        poly_Ui = galois.lagrange_poly(xs, L_mat[:, i])
        beta_Ui = poly_Ui * beta        # multiply U with beta
        # V_i(x) = interpolate the i-th column of matrix R
        poly_Vi = galois.lagrange_poly(xs, R_mat[:, i])
        alpha_Vi = poly_Vi * alpha
        # W_i(x) = interpolate the i-th column of matrix W
        poly_Wi = galois.lagrange_poly(xs, O_mat[:, i])
        # Get beta*U_i(tau) + alpha*V_i(tau) + W_i(tau)
        sum_tau = beta_Ui(tau) + alpha_Vi(tau) + poly_Wi(tau)
        list_taus.append(multiply(point, int(sum_tau)))

    return list_taus

def transform_matrix2poly(mat, wit, d):
    # Interpolate the matrix mat to a polynomial
    xs = Fp(np.array([i + 1 for i in range(d + 1)]))    #xs = Fp(np.array([1, 2, 3, 4, 5]))
    sum = galois.Poly([0] * (d + 1), field=Fp)
    for i in range(len(wit)):
        poly = galois.lagrange_poly(xs, mat[:, i])
        sum += poly * wit[i]
    return sum


if __name__ == '__main__':
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

    U_poly = transform_matrix2poly(L, witness, d)
    V_poly = transform_matrix2poly(R, witness, d)
    W_poly = transform_matrix2poly(O, witness, d)
    t_poly = galois.Poly.Roots(np.array([i + 1 for i in range(d + 1)]), field=Fp)
    LR_product = U_poly * V_poly
    h_poly = (LR_product - W_poly) // t_poly


    """
        Trusted setup
        Phase 1: setup for all circuits
    """
    # 1. Get a random value tau
    tau = Fp(sample_Zp(p))
    # 2. Calculate tau*G1, tauˆ2*G1, ..., tauˆd*G1
    powers_of_tau_G1 = generate_powers_of_tau(tau, d, G1)
    # Calculate tau*G2, tauˆ2*G2, ..., tauˆd*G2
    powers_of_tau_G2 = generate_powers_of_tau(tau, d, G2)

    # Phase 2: trusted setup per circuit. This setup requires for each individual circuit
    alpha = sample_Zp(p)
    beta = sample_Zp(p)
    alpha_G1 = multiply(G1, int(alpha))     # alpha*G1
    beta_G1 = multiply(G1, int(beta))       # beta*G1
    beta_G2 = multiply(G2, int(beta))       # beta*G2
    # Calculate alpha*V(tau) + beta*U(tau) + W(au)
    powers_of_tau_for_C = generate_powers_of_tau_W(tau, d, m, G1, L, R, O, alpha, beta)


    # Calculate tau*t(tau)*G1, tauˆ2*t(tau)*G1, ..., tauˆd*t(tau)*G1. Required to evaluate h(x)
    t_tau_G1 = multiply(G1, int(t_poly(tau)))
    powers_of_t_tau_G1 = generate_powers_of_tau(tau, d - 1, t_tau_G1)

    """
    Prover: compute
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

    # Given powers of tau from the trusted setup phase, prover compute U(tau)*G1, V(tau)*G2, and C*G1
    evaluate_U_on_G1 = inner_product(powers_of_tau_G1, U_poly.coeffs[::-1])     # U(tau)*G1
    evaluate_V_on_G2 = inner_product(powers_of_tau_G2, V_poly.coeffs[::-1])     # V(tau)*G2
    """
        Secret shifting
        Introducing alpha and beta to prevent a malicious prover to make up values U(tau)*G1, V(tau)*G2 and C
        U'(x) = alpha + U(x), V'(x) = beta + V(x)
    """
    alpha_U_G1 = add(alpha_G1, evaluate_U_on_G1)        # [alpha + U(x)]*G1
    beta_V_G2 = add(beta_G2, evaluate_V_on_G2)          # [beta + V(x)]*G2
    evaluate_ht_on_G1 = inner_product(powers_of_t_tau_G1, h_poly.coeffs[::-1])

    """
        Compute W'= W + beta*U + alpha*V
            Option 1: Trusted server requires the prover to send polynomials U(x), V(x), W(x), t(x)
                Question: Can a trusted server deduce the witness from polynomials U(x), V(x) and W(x)?
        """
    """
    beta_U_tau_G1 = multiply(G1, int(beta*U_poly(tau)))
    alpha_V_tau_G1 = multiply(G1, int(alpha*V_poly(tau)))
    W_tau_G1 = multiply(G1, int(W_poly(tau)))
    pt_tmp = add(beta_U_tau_G1, alpha_V_tau_G1)
    W_beta_U_alpha_V_G1 = add(pt_tmp, W_tau_G1)
    """

    """
        Option 2: Prover only sends matrices L, R, and O to the trusted server, that are public to anyone
    """
    W_beta_U_alpha_V_G1 = inner_product(powers_of_tau_for_C, witness)
    Wht_G1 = add(W_beta_U_alpha_V_G1, evaluate_ht_on_G1)

    """
    Verifying: verifier obtains a proof from prover, consisting of three elliptic curve points:
        - beta_V_G2
        - alpha_U_G1     
        - Wht_G1
        and two points from trusted server:
        - alpha_G2
        - beta_G1
        
    Accept the proof if the below equation is true. Otherwise, reject
    """
    print("Proof Verification ...")
    assert pairing(beta_V_G2, alpha_U_G1) == pairing(beta_G2, alpha_G1) * pairing(G2, Wht_G1), "Failed to check pairings"