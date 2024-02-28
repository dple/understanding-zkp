"""
Demonstrate some ZKP examples in Python
"""

# 1. Prove b is the inverse of a
# The order of the base field of the BN254 curve
base_field = 21888242871839275222246405745257275088548364400416034343698204186575808495617
def compute_inverse(a):
    return pow(a, -1, base_field)

# Prover --> compute the inverse of a, then give it to the verifier
a = 15
b = compute_inverse(a)

# Verifier --> just verify, won't carry out computation
assert (a*b) % base_field == 1, "Given number is not an inverse of a"

# 2. Prove x is the solution of a polynomial
def polyeval(x, a, b, c):
    return a*x**2 + b*x + c

# Prover --> proving that 11 is a solution of the poly 2xˆ2 - 23x + 11
x = 11
a = 2
b = -23
c = 11
s = polyeval(x, a, b, c)

# Verifier  --> just verify, won't carry out computation
assert s == 0, "Given x is not a solution of the polynomial"

# 3. Proving x is the maximum element in a list of elements
def maxList(lst):
    return max(lst)

# Prover --> compute the max of the list then given it to the verifier
lst = [5, 7, 1, 9, 3, 11]
m = maxList(lst)
x = 11

# Verifier --> just verify, won't carry out computation
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

