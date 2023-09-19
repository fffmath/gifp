import logging
import timeit

from math import ceil
from sage.all import Integer
from sage.all import inverse_mod
from sage.crypto.util import random_blum_prime

import boneh_durfee


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
    p, q = boneh_durfee.attack(N, e, p_bits, delta, m)
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

# p_bits = 512
# modulus_bit_length = p_bits * 2
# delta = 0.25
# test_times = 5
# m = 5

p_bits = int(input("Input prime bit-length: "))
delta = float(input("Input delta: "))
test_times = int(input("Input attack times: "))
m = int(input("Input m: "))
modulus_bit_length = p_bits * 2

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
