from py_ecc.bn128 import G1, G2, pairing, add, multiply, eq, neg

if __name__ == '__main__':
    #print(G1)
    #print(G2)
    P = multiply(G1, 3)
    Q = add(G2, G2)
    """
    py-ecc library implemented Optimal Ate pairing (https://eprint.iacr.org/2006/110.pdf), 
    that switches the order of P & Q, i.e., the first argument of pairing will be the element, Q, 
    in the extension field F_pˆ12, and the second argument, P, will be in the base field F_p   
    """
    print("pairing(Q, P) = ", pairing(Q, P))        # The result will be an element in the extension field, F_pˆ12,
                                # representing by 12 element in the base field.
                                # The computation takes a few seconds, pretty slow,
                                # perhaps the implementation of pairing in py-ecc is not optimized

    # As pairing(a*G2, b*G1) == pairing(a*b*G2, G1), the following return true
    T = multiply(G2, 6)
    print(eq(pairing(Q, P), pairing(T, G1)))
    # And this will print 1 as e(G2, G2)ˆ{-ab} * e(G2, G2)ˆ{ab} = e(G2, G1)ˆ0 = 1
    print(pairing(neg(Q), P) * pairing(T, G1))
