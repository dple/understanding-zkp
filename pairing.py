from py_ecc.bn128 import G1, G2, pairing, add, multiply, eq

if __name__ == '__main__':
    print(G1)
    print(G2)
    P = multiply(G1, 3)
    Q = multiply(G2, 2)
    #print(pairing(P, Q))
    A = multiply(G2, 5)
    B = multiply(G1, 6)
    print(pairing(A, B))