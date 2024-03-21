"""
This Python code demonstrate the examples discussed in the following article by Rareskills:
    How arithmetic circuits are used to verify zero knowledge proofs
            https://www.rareskills.io/post/zk-circuits

"""

from py_ecc.bn128 import curve_order

# 1. Prove b is the inverse of a
# The order of the base field of the BN254 curve
p = curve_order
def compute_inverse(a):
    return pow(a, -1, p)

# Prover --> compute the inverse of a, then give it to the verifier
a = 15
b = compute_inverse(a)

# Verifier --> just verify, won't carry out any inversion
assert (a*b) % p == 1, "Given number b is not an inverse of a"

# 2. Prove x is the solution of a polynomial
def polyeval(x, a, b, c):
    return a*x**2 + b*x + c

# Prover --> proving that 11 is a solution of the poly 2xˆ2 - 23x + 11
x = 11
a = 2
b = -23
c = 11

# Verifier  --> just verify if x is the solution, don't  carry out any computations how to find x
assert 0 == polyeval(x, a, b, c), "Given x is not a solution of the polynomial"

# 3. Proving x is the maximum element in a list of elements
def maxList(lst):
    return max(lst)

# Prover --> compute the max of the list then given it to the verifier
lst = [5, 7, 1, 9, 3, 11]
m = maxList(lst)
x = 11

# Verifier --> just verify, don't execute any computations how to find the max value
assert x == m, "Given input is not the maximum of the list"

# 4. Prove if the given list was sorted
# Prover --> sort the list then give it to the verifier
lst = [1, 2, 3, 4, 5, 6]
sortedList = sorted(lst)

# Verifier
assert lst == sortedList, "List is not sorted"

# 5. Prove if the list has no duplicates
# Prover --> convert the list to a set, then compare their size
s = set(lst)

# Verifier
assert len(lst) == len(s), "List contains duplicates"

