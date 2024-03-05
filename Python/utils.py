import secrets
import math

q = 1337    # A composite number
n = 32

def sample_Zq():
    return secrets.randbelow(q)

def sample_bits(k):
    return secrets.randbits(k)

def rejection_sampling_Zq():
    while True:
        e = sample_bits(n)
        if e >=0 and e < q:
            return e
        
def sample_Zqstar():
    while True:
        e = sample_Zq()
        if math.gcd(e, q) == 1:
            return e
