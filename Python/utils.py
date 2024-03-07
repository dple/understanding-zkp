import secrets
import math

def sample_Zp(p):
    return secrets.randbelow(p)

def sample_bits(k):
    return secrets.randbits(k)

def rejection_sampling_Zp(k):
    while True:
        e = sample_bits(k)
        if e >=0 and e < q:
            return e
        
def sample_Zqstar(p):
    while True:
        e = sample_Zp(p)
        if math.gcd(e, p) == 1:
            return e
