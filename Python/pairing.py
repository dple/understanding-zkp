from py_ecc.bn128 import G1, G2, pairing, add, multiply, eq, neg

if __name__ == '__main__':
    print(G1)
    print(G2)
    P = multiply(G1, 3)
    Q = add(G2, G2)
    """
    py-ecc library implemented Ate pairing, that switches the order of P & Q, 
    that is the first argument of pairing will be the element, Q, in the extension field F_pˆ12, 
    and the second argument, P, will be in the base field F_p   
    """
    print(pairing(Q, P))
    T = multiply(G2, 6)
    # As pairing(a*G2, b*G1) == pairing(a*b*G2, G1), the following return true
    print(eq(pairing(Q, P), pairing(T, G1)))
    # And this will print 1 as e(G2, G2)ˆ{-ab} * e(G2, G2)ˆ{ab} = 1
    print(pairing(neg(Q), P) * pairing(T, G1))