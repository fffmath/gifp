import logging
import timeit

from math import ceil
from sage.all import Integer
from sage.all import inverse_mod
from sage.crypto.util import random_blum_prime
from sage.all import RR
from sage.all import ZZ
from sage.all import Integer

import small_roots

def modular_bivariate(f, e, m, t, X, Y, roots_method="groebner"):
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

    pr = ZZ["x", "y", "u"]
    x, y, u = pr.gens()
    qr = pr.quotient(1 + x * y - u)
    U = X * Y

    logging.debug("Generating shifts...")

    shifts = []
    for k in range(m + 1):
        for i in range(m - k + 1):
            g = x ** i * f ** k * e ** (m - k)
            g = qr(g).lift()
            shifts.append(g)

    for j in range(1, t + 1):
        for k in range(m // t * j, m + 1):
            h = y ** j * f ** k * e ** (m - k)
            h = qr(h).lift()
            shifts.append(h)

    L, monomials = small_roots.create_lattice(pr, shifts, [X, Y, U])
    L = small_roots.reduce_lattice(L)

    pr = f.parent()
    x, y = pr.gens()

    polynomials = small_roots.reconstruct_polynomials(L, f, None, monomials, [X, Y, U], preprocess_polynomial=lambda p: p(x, y, 1 + x * y))
    for roots in small_roots.find_roots(pr, polynomials, method=roots_method):
        yield roots[x], roots[y]

def attack(N, e, factor_bit_length, delta=0.25, m=1, t=None):
    """
    Recovers the prime factors if the private exponent is too small.
    This implementation exploits knowledge of least significant bits of prime factors, if available.
    More information: Boneh D., Durfee G., "Cryptanalysis of RSA with Private Key d Less than N^0.292"
    :param N: the modulus
    :param e: the public exponent
    :param factor_bit_length: the bit length of the prime factors
    :param delta: a predicted bound on the private exponent (d < N^delta) (default: 0.25)
    :param m: the m value to use for the small roots method (default: 1)
    :param t: the t value to use for the small roots method (default: automatically computed using m)
    :return: a tuple containing the prime factors
    """
    A = N + 1
    x, y = ZZ["x", "y"].gens()
    f = x * (A + y) + 1
    X = Integer(RR(e) ** delta)
    Y = Integer(2 ** int(factor_bit_length + 1))
    t = int((1 - 2 * delta) * m) if t is None else t
    logging.info(f"Trying m = {m}, t = {t}...")
    for x0, y0 in modular_bivariate(f, e, m, t, X, Y):
        z = int(f(x0, y0))
        if z % e == 0:
            k = pow(x0, -1, e)
            s = (N + 1 + k) % e
            phi = N - s + 1
            factors = small_roots.factorize(N, phi)
            if factors:
                return factors

    return None, None

def generate_RSA_instance(modulus_bit_length, delta=0.25):
    """
    Generate an RSA instance with given bit-lengths of the modulus and private key d
    :param modulus_bit_length: the bit length of the modulus
    :param delta: a given size on the private exponent (d is roughly N^delta) (default: 0.25)
    :return: a tuple containing the public key (N, e)
    """
    # if modulus_bit_length < 1024:
    #     raise ValueError("RSA modulus length must be >= 1024")
    e = d = N = Integer(1)
    d_bit_length = ceil(modulus_bit_length * delta)
    logging.info(f"Generating RSA instance with {modulus_bit_length}-bit modulus and {d_bit_length}-bit private key d...")
    while N.nbits() != modulus_bit_length and d.nbits() != d_bit_length:
        prime_bit_length = modulus_bit_length // 2
        p = random_blum_prime(2**(prime_bit_length - 1), 2**prime_bit_length - 1)
        q = random_blum_prime(2**(prime_bit_length - 1), 2**prime_bit_length - 1)
        N = p * q
        phi = (p - 1) * (q - 1)
        d = random_blum_prime(2**(d_bit_length - 1), 2**d_bit_length - 1)
        e = inverse_mod(d, phi)

    return N, e

def attack_RSA_instance(modulus_bit_length, delta=0.25, m=3):
    """
    Attack an RSA instance with given bit-lengths of the modulus and private key d
    :param modulus_bit_length: the bit length of the modulus
    :param delta: a given size on the private exponent (d is roughly N^delta) (default: 0.25)
    :param m: a given parameter for controlling the lattice dimension (default: 3)
    :return: 1 if attack succeeds else 0
    """
    N, e  = generate_RSA_instance(modulus_bit_length, delta)
    p_bits = modulus_bit_length / 2
    p, q = attack(N, e, p_bits, delta, m)
    if p is not None and q is not None and p * q == N:
        logging.info(f"Succeeded!")
        logging.info(f"Found p = {p}")
        logging.info(f"Found q = {q}")
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

p_bits = 512
modulus_bit_length = p_bits * 2
delta = 0.25
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
