import logging
import timeit

from math import ceil, sqrt
from sage.all import Integer
from sage.all import inverse_mod
from sage.crypto.util import random_blum_prime
from sage.all import RR
from sage.all import ZZ
from sage.all import Integer

import small_roots

def modular_bivariate(f, N1, N2, m, s, t, X, Y, Z, W, beta1_bit_length, beta2_bit_length, roots_method="groebner"):
    """
    Computes small modular roots of a bivariate polynomial.
    More information: Herrmann M., May A., "Maximizing Small Root Bounds by Linearization and Applications to Small Secret Exponent RSA"
    :param f: the polynomial
    :param e: the modulus
    :param m: the amount of normal shifts to use
    :param t: the amount of additional shifts to use
    :param X: an approximate bound on the x roots
    :param Y: an approximate bound on the y roots
    :param roots_method: the method to use to find roots (default: "groebner")
    :return: a generator generating small roots (tuples of x and y roots) of the polynomial
    """
    f = f.change_ring(ZZ)

    pr = ZZ["x", "y", "z", "w"]
    x, y, z, w = pr.gens()
    qr = pr.quotient(N2 - z * w)

    logging.debug("Generating shifts...")

    shifts = []
    for i in range(m + 1):
        for j in range(m - i + 1):
            g = (y*z) ** j * f ** i * 2 ** (int(beta2_bit_length-beta1_bit_length)*(m-i)) * N1 ** max(t-i, 0)
            g = qr(g).lift()
            shifts.append(g)

    L, monomials = small_roots.create_lattice(pr, shifts, [X, Y, Z, W])
    L = small_roots.reduce_lattice(L)

    pr = f.parent()
    x, y = pr.gens()

    polynomials = small_roots.reconstruct_polynomials(L, f, None, monomials, [X, Y, U], preprocess_polynomial=lambda p: p(x, y, 1 + x * y))
    for roots in small_roots.find_roots(pr, polynomials, method=roots_method):
        yield roots[x], roots[y]

def attack(modulus_bit_length, N1, N2, alpha=0.25, gamma=0.75, m=1, s=None, t=None):
    """
    Recovers the prime factors if the private exponent is too small.
    This implementation exploits knowledge of least significant bits of prime factors, if available.
    More information: Boneh D., Durfee G., "Cryptanalysis of RSA with Private Key d Less than N^0.292"
    :param N1 and param N2: the modulus
    :param alpha: the prime number are bounded by N^alpha (default: 0.25)  
    :param gamma: a predicted bound on the shared bits (default: 0.75) 
    :param m: the m value to use for the small roots method (default: 1)
    :param t: the t value to use for the small roots method (default: automatically computed using m)
    :return: a tuple containing the prime factors
    """
    for beta1_bit_length in range(ceil(modulus_bit_length * (1-alpha-gamma))):
        for beta2_bit_length in range(beta1_bit_length, ceil(modulus_bit_length * (1-alpha-gamma))): # Suppose beta2>beta1, otherwise, we exchange the values of N1 and N2.
            x, y, z = ZZ["x", "y", "z"].gens()
            f = x*z+2**(ceil(modulus_bit_length * (gamma))+beta2_bit_length)*y*z+N2
            X = Integer(2 ** int(beta2_bit_length))
            Y = Integer(2 ** int(modulus_bit_length*(1-alpha-gamma)-beta1_bit_length))
            Z = Integer(2 ** int(modulus_bit_length*(alpha)))
            W = Integer(2 ** int(modulus_bit_length*(1-alpha)))
            s = int(sqrt(alpha) * m) if s is None else s
            t = int((1 - sqrt(alpha)) * m) if t is None else t
            logging.info(f"Trying m = {m}, s = {s}, t = {t}...")
            for x0, y0, z0 in modular_bivariate(f, N1, N2, m, s, t, X, Y, Z, W, beta1_bit_length, beta2_bit_length):
                result = int(f(x0, y0, z0))
                if result == 0:
                    q2=z0
                    p2=N2/z0
                    p1=(p2+y0+z0*2**(beta2_bit_length+ceil(modulus_bit_length * gamma)))/(2**(beta2_bit_length-beta1_bit_length))
                    q1=N1/p1
                    return p1, q1, p2, q2

    return None, None, None, None

def generate_RSA_instance(modulus_bit_length, alpha=0.25, gamma=0.75):
    """
    Generate an RSA instance with given bit-lengths of the modulus and private key d
    :param modulus_bit_length: the bit length of the modulus
    :param delta: a given size on the private exponent (d is roughly N^delta) (default: 0.25)
    :return: a tuple containing the public key (N, e)
    """
    # if modulus_bit_length < 1024:
    #     raise ValueError("RSA modulus length must be >= 1024")
    shared_bit = N1 = N2 = Integer(1)
    q_bit_length=ceil(modulus_bit_length * alpha)
    p_bit_length=modulus_bit_length - q_bit_length
    shared_bit_length = ceil(modulus_bit_length * gamma)
    logging.info(f"Generating RSA instance with {modulus_bit_length}-bit modulus and {shared_bit_length}-bit implicit shared bits for p...")
    while N1.nbits() != modulus_bit_length and N2.nbits() != modulus_bit_length and shared_bit.nbits() != shared_bit_length:
        p1 = random_blum_prime(2**(p_bit_length - 1), 2**p_bit_length - 1)
        q1 = random_blum_prime(2**(q_bit_length - 1), 2**q_bit_length - 1)
        p2 = random_blum_prime(2**(p_bit_length - 1), 2**p_bit_length - 1)
        q2 = random_blum_prime(2**(q_bit_length - 1), 2**q_bit_length - 1)
        N1 = p1 * q1
        N2 = p2 * q2

    return N1, N2

def attack_RSA_instance(modulus_bit_length, alpha, gamma, m=3):
    """
    Attack an RSA instance with given bit-lengths of the modulus and private key d
    :param modulus_bit_length: the bit length of the modulus
    :param delta: a given size on the private exponent (d is roughly N^delta) (default: 0.25)
    :param m: a given parameter for controlling the lattice dimension (default: 3)
    :return: 1 if attack succeeds else 0
    """
    N1, N2  = generate_RSA_instance(modulus_bit_length, alpha, gamma)
    p_bits = modulus_bit_length - ceil(modulus_bit_length * alpha)
    p1, q1, p2, q2 = attack(modulus_bit_length, N1, N2, alpha, gamma, m)
    if p1 is not None and q1 is not None and p1 * q1 == N1 and p2 is not None and q2 is not None and p2 * q2 == N2:
        logging.info(f"Succeeded!")
        logging.info(f"Found p1 = {p1}")
        logging.info(f"Found q1 = {q1}")
        logging.info(f"Found p2 = {p2}")
        logging.info(f"Found q2 = {q2}")
        return 1
    else:
        logging.info(f"Failed!")
        return 0

# Some logging so we can see what's happening.
# logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)

'''
p_bits = int(input("Input prime bit-length: "))
delta = float(input("Input delta: "))
test_times = int(input("Input attack times: "))
m = int(input("Input m: "))
modulus_bit_length = p_bits * 2
'''

modulus_bit_length = 1024
alpha = 0.25
gamma = 0.75
test_times = 5
m = 5

logging.info(f"Test for delta={delta} with {test_times} times:")
total_time = 0
results = []
for i in range(test_times):
    start_time = timeit.default_timer()
    result = attack_RSA_instance(modulus_bit_length, delta=delta, m=m)
    end_time = timeit.default_timer()
    test_time = end_time - start_time
    if result:
        total_time += test_time
        results.append(result)
    logging.info(f"Test {i+1} costs {test_time:.6f} seconds")

avg_time = total_time / len(results)
logging.info(f"The success rate for delta={delta} using m={m} is {sum(results)/test_times*100}% and average time is {avg_time:.6f} seconds")
