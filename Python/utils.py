import secrets
import math

def sample_Zp(p):
    return secrets.randbelow(p)

def sample_bits(k):
    return secrets.randbits(k)
        
def sample_Zqstar(p):
    while True:
        e = sample_Zp(p)
        if math.gcd(e, p) == 1:
            return e
