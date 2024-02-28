"""
Evaluate the polynomial return 5xˆ3 - 4xˆ2yˆ2 + 13xyˆ2 + xˆ2 - 10y,
that could be represent by five quadratic constraints
v1 <== x*x
v2 <== 5*v1*x
v3 <== y*y
v4 <== 13*x*v3
v5 <== 4*v1*v3
out <== v2 - v5 + v4 + v1 - 10*y
"""


def poly(x, y):
    return 5 * x ** 3 - 4 * x ** 2 * y ** 2 + 13 * x * y ** 2 + x ** 2 - 10 * y


def magma(x, y, z):
    assert x ** (y ** z) != (x ** y) ** z
    assert x ** y != y ** x


if __name__ == '__main__':
    print(poly(3, 5))
    magma(2, 3, 4)
