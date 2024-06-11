import numpy as np
from scipy.interpolate import lagrange

if __name__ == '__main__':
    x = np.array([1, 2, 3])
    y_v1 = np.array([4, 12, 8])     # vector v1
    y_v2 = np.array([2, 2, -2])     # vector v2

    poly_v1 = lagrange(x, y_v1)     # Using lagrange to compute a polynomial from v1
    print("Poly v1 = ", poly_v1)
    poly_v2 = lagrange(x, y_v2)     # Using lagrange to compute a polynomial from v2
    print("Poly v2 = ", poly_v2)

    # 1. Adding two vectors v1, v2
    y_vA = y_v1 + y_v2
    # Homomorphic to adding two polynomials
    poly_vA = lagrange(x, y_vA)
    print("Sum of polynomials v1 + v2:", poly_vA)
    result = poly_vA == poly_v1 + poly_v2

    assert result.all(), "Add not equal!"

    # 2. Multiplication (Hadamard product) of two vectors v1 and v2
    y_vM = np.multiply(y_v1, y_v2)      # this returns a vector of 3 elements
                                        # If you interpolate this vector to a polynomial,
                                        # it will be a polynomial of degree 2
    print("Element-wise (or Hadamard) product of v1 & v2", y_vM)
    # Homomorphic to multiplying two polynomials
    poly_vM = poly_v1 * poly_v2         # a polynomial of degree 4
    print("Product of two polynomials: ", poly_vM)
    print("Evaluate the poly product at 1, 2, 3 =", poly_vM(1), poly_vM(2),  poly_vM(3))
    result = y_vM == np.array([poly_vM(1), poly_vM(2),  poly_vM(3)]).astype(int)

    assert result.all(), "Multiply not equal!"

    # 3. Multiplying a vector by a scalar
    y_v3 = y_v1 * 3
    print("Scalar multiplication of vector v1: ", y_v3)
    # Homomorphic to multiplying the polynomial by the same scalar
    poly_v3 = poly_v1 * 3

    result = lagrange(x, y_v3) == poly_v3

    assert result.all()