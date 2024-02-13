import random
import time
import logging
import sys
from sage.crypto.util import random_blum_prime
logging.basicConfig(filename='gifp.log', level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
def generate_gifp_instance(modulus_bit_length, alpha, gamma, beta1, beta2, seed=None, max_attempts=10):
    """
    Generate two RSA instances with the same modulus length and a shared bit. 
    :param modulus_bit_length: The bit length of the modulus.
    :param alpha: The ratio of the bit length of the larger prime to the modulus bit length.
    :param gamma: The ratio of the bit length of the shared bit to the modulus bit length.
    :param beta1: The ratio of the bit length of the smaller prime in the first RSA instance to the modulus bit length.
    :param beta2: The ratio of the bit length of the smaller prime in the second RSA instance to the modulus bit length.
    :param seed: The seed for the random number generator. If not provided, the current time's microsecond is used.
    :param delta: The probability of failure in generating RSA instances.
    :param max_attempts: The maximum number of attempts to generate RSA instances.
    :return: A list of the first RSA instance's primes, a list of the second RSA instance's primes, the shared bit, and the seed.
    """
    N1 = N2 = attempts = 0
    while attempts < max_attempts:
            # If the seed parameter is not provided, use the current time's microsecond as the seed.
        if seed is None:
            seed = int(time.time() * 1e6)
        set_random_seed(seed)
        if beta2<beta1:
            tmp=beta2
            beta2=beta1
            beta1=tmp
        else:
            beta1=beta1
            beta2=beta2
        share_bit_length = int(modulus_bit_length * gamma)
        beta1_bit_length = int(modulus_bit_length * beta1)
        beta2_bit_length = int(modulus_bit_length * beta2)
        share_bit = ZZ(randint(2**(share_bit_length - 1)+1, 2**share_bit_length - 1))
        p_bit_length = int(modulus_bit_length * (1-alpha))
        q_bit_length = int(modulus_bit_length * alpha)
        MSB4p1 = ZZ(randint(2**int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*gamma-modulus_bit_length*beta1-1)+1, 2**int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*gamma-modulus_bit_length*beta1) - 1))
        MSB4p2 = ZZ(randint(2**int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*gamma-modulus_bit_length*beta2-1)+1, 2**int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*gamma-modulus_bit_length*beta2) - 1))
        p1 = random_blum_prime(MSB4p1*2**int(modulus_bit_length*gamma+modulus_bit_length*beta1) + share_bit* 2**(beta1_bit_length)+2**(beta1_bit_length - 1), MSB4p1*2**int(modulus_bit_length*gamma+modulus_bit_length*beta1) + share_bit* 2**(beta1_bit_length)+2**beta1_bit_length - 1)
        p2 = random_blum_prime(MSB4p2*2**int(modulus_bit_length*gamma+modulus_bit_length*beta2) + share_bit* 2**(beta2_bit_length)+2**(beta2_bit_length - 1), MSB4p2*2**int(modulus_bit_length*gamma+modulus_bit_length*beta2) + share_bit* 2**(beta2_bit_length)+2**beta2_bit_length - 1)
        q1 = random_blum_prime(2**(q_bit_length - 1), 2**q_bit_length - 1)
        q2 = random_blum_prime(2**(q_bit_length - 1), 2**q_bit_length - 1)
        N1 = p1 * q1
        N2 = p2 * q2
        desired_solution = (int(bin(p1)[2+int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*beta1):],2)* 2**int(beta2*modulus_bit_length-beta1*modulus_bit_length)-int(bin(p2)[2+int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*beta2):],2),int(bin(MSB4p1)[2:],2)-int(bin(MSB4p2)[2:],2),q2,p2)
        if N1.nbits() != modulus_bit_length or N2.nbits() != modulus_bit_length:
            attempts += 1
            seed = int(time.time() * 1e6)
            logging.info(f"Regenerated seed: {seed}")
        else:
            logging.info(f"Geenerated gifp instance successfully with seed: {seed}")
            logging.debug(f'share_bit: {bin(share_bit)[2:]}')
            logging.debug(f'p1: {p1}')
            logging.debug(f'MSB4p1: {bin(MSB4p1)[2:]}')
            logging.debug(f'share_bit4p1: {bin(p1)[2+int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*gamma-modulus_bit_length*beta1):2+int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*beta1)]}')
            logging.debug(f'LSB4p1: {bin(p1)[2+int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*beta1):]}')
            logging.debug(f'q1: {q1}')
            logging.debug(f'N1: {N1}')
            logging.debug(f'p2: {p2}')
            logging.debug(f'MSB4p2: {bin(MSB4p2)[2:]}')
            logging.debug(f'share_bit4p2: {bin(p2)[2+int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*gamma-modulus_bit_length*beta2):2+int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*beta2)]}')
            logging.debug(f'LSB4p2: {bin(p2)[2+int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*beta2):]}')
            logging.debug(f'q2: {q2}')
            logging.debug(f'N2: {N2}')
            logging.debug(f"The desired solution is (x0, y0, z0, w0) = {(int(bin(p1)[2+int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*beta1):],2)* 2**int(beta2*modulus_bit_length-beta1*modulus_bit_length)-int(bin(p2)[2+int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*beta2):],2),int(bin(MSB4p1)[2:],2)-int(bin(MSB4p2)[2:],2),q2,p2)}")
            return [p1, q1, N1], [p2, q2, N2], share_bit, desired_solution
    logging.warning(f"Failed to generate RSA instances after {max_attempts} attempts.")
    return None

def create_lattice(pr, shifts, bounds, order="invlex", sort_shifts_reverse=False, sort_monomials_reverse=False):
    """
    Creates a lattice from a list of shift polynomials.
    :param pr: the polynomial ring
    :param shifts: the shifts
    :param bounds: the bounds
    :param order: the order to sort the shifts/monomials by
    :param sort_shifts_reverse: set to true to sort the shifts in reverse order
    :param sort_monomials_reverse: set to true to sort the monomials in reverse order
    :return: a tuple of lattice and list of monomials
    """
    logging.debug(f"Creating a lattice with {len(shifts)} shifts (order = {order}, sort_shifts_reverse = {sort_shifts_reverse}, sort_monomials_reverse = {sort_monomials_reverse})...")
    if pr.ngens() > 1:
        pr_ = pr.change_ring(ZZ, order=order)
        shifts = [pr_(shift) for shift in shifts]

    monomials = set()
    for shift in shifts:
        monomials.update(shift.monomials())

    shifts.sort(reverse=sort_shifts_reverse)
    monomials = sorted(monomials, reverse=sort_monomials_reverse)
    L = matrix(ZZ, len(shifts), len(monomials))
    for row, shift in enumerate(shifts):
        for col, monomial in enumerate(monomials):
            L[row, col] = shift.monomial_coefficient(monomial) * monomial(*bounds)

    monomials = [pr(monomial) for monomial in monomials]
    return L, monomials

def log_lattice(L):
    """
    Logs a lattice.
    :param L: the lattice
    """
    #for row in range(L.nrows()):
    #    logging.debug(L[row,:])
    for row in range(L.nrows()):
        r = ""
        for col in range(L.ncols()):
            if L[row, col] == 0:
                r += "_ "
            else:
                r += "X "
        logging.info(r)
def reduce_lattice(L, delta=0.8):
    """
    Reduces a lattice basis using a lattice reduction algorithm (currently LLL).
    :param L: the lattice basis
    :param delta: the delta parameter for LLL (default: 0.8)
    :return: the reduced basis
    """
    logging.debug(f"Reducing a {L.nrows()} x {L.ncols()} lattice...")
    return L.LLL(delta)

def reconstruct_polynomials(B, f, modulus, monomials, bounds):
    """
    Reconstructs polynomials from the lattice basis in the monomials.
    :param B: the lattice basis
    :param f: the original polynomial (if set to None, polynomials will not be divided by f if possible)
    :param modulus: the original modulus
    :param monomials: the monomials
    :param bounds: the bounds
    :param preprocess_polynomial: a function which preprocesses a polynomial before it is added to the list (default: identity function)
    :param divide_gcd: if set to True, polynomials will be pairwise divided by their gcd if possible (default: True)
    :return: a list of polynomials
    """
    logging.debug(f"Reconstructing polynomials (divide_original = {f is not None}, modulus_bound = {modulus is not None})...")
    logging.debug(f"Reconstructing polynomials with bounds: {bounds}")
    polynomials = []
    for row in range(B.nrows()):
        norm_squared = 0
        ww = 0
        polynomial = 0
        for col, monomial in enumerate(monomials):
            if B[row, col] == 0:
                continue
            norm_squared += B[row, col] ** 2
            ww += 1
            assert B[row, col] % monomial(*bounds) == 0
            polynomial += B[row, col] * monomial // monomial(*bounds)

        # Equivalent to norm >= modulus / sqrt(w)
        if modulus is not None and norm_squared * ww >= modulus ** 2:
            logging.debug(f"Row {row} is too large, ignoring...")
            continue

        if f is not None and polynomial % f == 0:
            logging.debug(f"Original polynomial divides reconstructed polynomial at row {row}, dividing...")
            polynomial //= f

        if polynomial.is_constant():
            logging.debug(f"Polynomial at row {row} is constant, ignoring...")
            continue

        polynomials.append(polynomial)

    logging.debug(f"Reconstructed {len(polynomials)} polynomials")
    return polynomials

def find_roots_univariate(x, polynomial):
    """
    Returns a generator generating all roots of a univariate polynomial in an unknown.
    :param x: the unknown
    :param polynomial: the polynomial
    :return: a generator generating dicts of (x: root) entries
    """
    if polynomial.is_constant():
        return

    for root in polynomial.roots(multiplicities=False):
        if root != 0:
            yield {x: int(root)}

def find_roots_groebner(N1, N2,pr,desired_solution, unknown_modular, polynomials,bounds):
    """
    Returns a generator generating all roots of a polynomial in some unknowns.
    Uses Groebner bases to find the roots.
    :param pr: the polynomial ring
    :param polynomials: the reconstructed polynomials
    :return: a generator generating dicts of (x0: x0root, x1: x1root, ...) entries
    """
    gens = pr.gens()
    x, y, z, w = pr.gens()
    #qr = pr.change_ring(QQ).quotient(z*w - N2)
    #tmp=[qr(g).lift() for g in polynomials]
    #polynomials=list(set(tmp))
    s = Sequence(polynomials, pr.change_ring(QQ,order='lex'))
    for ii in range(len(s)):
        logging.debug(f"Polynomials in Sequence: {s[ii]}")
        logging.debug(f'Test for Sequence: s(x0, y0, z0, w0) mod modular = {s[ii](desired_solution)% unknown_modular}')
        logging.debug(f'Test for Sequence: s(x0, y0, z0, w0) = {s[ii](desired_solution)} over Z')
    while len(s) > 0:
        start_time_gb=time.time()
        G = s.groebner_basis()
        logging.debug(f"Sequence length: {len(s)}, Groebner basis length: {len(G)}")
        if len(G) == len(gens):
            end_time_gb=time.time()
            logging.info(f'The time taken by gb is {end_time_gb-start_time_gb}')
            logging.debug(f"Found Groebner basis with length {len(gens)}, trying to find roots...")
            roots = {}
            for ii in range(len(G)):
                logging.debug(f"Groebner basis {G[ii]}")
            for polynomial in G:
                vars = polynomial.variables()
                if len(vars) == 1:
                    for root in find_roots_univariate(vars[0], polynomial.univariate_polynomial()):
                        roots |= root

            if len(roots) == pr.ngens():
                yield roots, G
                return

            return
        else:
            # Remove last element (the biggest vector) and try again.
            s.pop()

def eliminate_N2(f, modular):
    pr = ZZ["x", "y", "z", "w"]
    x, y, z, w = pr.gens()
    tmp_poly=0
    for mono in f.monomials():
        if f.monomial_coefficient(mono)%modular==0:
            tmp_poly+=mono*f.monomial_coefficient(mono)
        else:
            tmp_poly+=mono*(f.monomial_coefficient(mono)% modular)
    return tmp_poly

def modular_gifp(desired_solution, unknown_modular, f, M, N1, N2, m, s, t, X, Y, Z, W):
    """
    Computes small modular roots of a bivariate polynomial.
    More information: Herrmann M., May A., "Maximizing Small Root Bounds by Linearization and Applications to Small Secret Exponent RSA"
    :param f: the polynomial
    :param M: 2^(|\beta1-\beta2|*n)
    :param N1: the modulus of the first RSA instance
    :param N2: the modulus of the second RSA instance
    :param m: the number of shifts
    :param s: a parameter
    :param t: a parameter
    :param X: the bound for x
    :param Y: the bound for y
    :param Z: the bound for z
    :param W: the bound for w
    :return: a generator generating small roots (tuples of x and y roots) of the polynomial
    """
    f = f.change_ring(ZZ)

    pr = ZZ["x", "y", "z", "w"]
    x, y, z, w = pr.gens()
    qr = pr.quotient(z*w - N2)

    logging.debug("Generating shifts...")

    shifts = []

    modular = M^m*N1^t

    logging.debug(f'modular is {modular}')
    N2_inverse = inverse_mod(N2, modular)
    logging.debug(f'N2^{-1}: {N2_inverse}')
    for ii in range(m + 1):
        for jj in range(m - ii + 1):
            g = (y*z) ** jj * w ** s * f ** (ii) * M ** (m - ii) * N1 ** max(t - ii, 0) * N2_inverse ** min(ii + jj, s)
            logging.debug(f'Initial polynomial shift: {g}')
            g = qr(g).lift()
            g = eliminate_N2(g, modular)
            shifts.append(g)
            logging.debug(f'Add shift: {g}')
            logging.debug(f'Test for shift: g(x0, y0, z0, w0) mod modular = {g(desired_solution)% unknown_modular}')
    logging.info('Generating lattice for gifp')
    L, monomials = create_lattice(pr, shifts, [X, Y, Z, W])
    log_lattice(L)
    logging.info('Reduce the lattice')
    start_time_LLL=time.time()
    L = reduce_lattice(L)
    end_time_LLL=time.time()
    logging.info(f'The time taken by LLL is {end_time_LLL-start_time_LLL}')
    log_lattice(L)

    polynomials = reconstruct_polynomials(L, f, unknown_modular, monomials, [X, Y, Z, W])
    for poly in polynomials:
        logging.debug(f'Polynomials after reconstructing from short vectors: {poly}')
        logging.debug(f'Test for reconstruct_polynomial: g(x0, y0, z0, w0) mod modular = {poly(desired_solution)% unknown_modular}')
        logging.debug(f'Test for reconstruct_polynomial: g(x0, y0, z0, w0) = {poly(desired_solution)} over Z')
    #for roots in find_roots_variety(pr,desired_solution, unknown_modular, polynomials,[X, Y, Z, W]):
    #    yield roots[x], roots[y], roots[z], roots[w]
    for roots in find_roots_groebner(N1, N2,pr,desired_solution, unknown_modular, polynomials,[X, Y, Z, W]):
        yield roots[x], roots[y], roots[z], roots[w]

def attack_gifp_instance(modulus_bit_length, alpha, gamma, beta1, beta2, m=1,seed=None, s=None, t=None):
    """
    Attack an RSA instance with given bit-lengths of the modulus and private key d
    :param modulus_bit_length: the bit length of the modulus
    :param delta: a given size on the private exponent (d is roughly N^delta) (default: 0.25)
    :param m: a given parameter for controlling the lattice dimension (default: 3)
    :return: 1 if attack succeeds else 0
    """
    if beta2<beta1:
        tmp=beta2
        beta2=beta1
        beta1=tmp
    else:
        beta1=beta1
        beta2=beta2
    N1_list, N2_list, share_bit, desired_solution = generate_gifp_instance(modulus_bit_length, alpha, gamma, beta1, beta2, seed,max_attempts=10)
    N1 = N1_list[2]
    N2 = N2_list[2]
    x, y, z, w = ZZ["x", "y", "z", "w"].gens()
    logging.info('Generating polynimial f')
    f = x*z + 2**(int(beta2*modulus_bit_length)+int(gamma*modulus_bit_length))*y*z+N2
    logging.debug(f'Polynomial f: {f}')
    X = Integer(2 ** int(beta2*modulus_bit_length))
    Y = Integer(2 ** int(modulus_bit_length*1-modulus_bit_length*alpha-modulus_bit_length*gamma-modulus_bit_length*beta1))
    Z = Integer(2 ** int(alpha * modulus_bit_length))
    W = Integer(2 ** int(modulus_bit_length-modulus_bit_length*alpha))
    M = Integer(2 ** int(modulus_bit_length*beta2-modulus_bit_length*beta1))
    logging.debug(f'M = {M}')
    t = round((1 - sqrt(alpha)) * m) if t is None else t
    s = round(sqrt(alpha) * m) if s is None else s   
    unknown_modular = M ** m * N1_list[0] ** t
    logging.info(f"Trying m = {m}, t = {t}, s = {s}...")
    for x0, y0, z0, w0 in modular_gifp(desired_solution, unknown_modular, f, M, N1, N2, m, s, t, X, Y, Z, W):
        v = int(f(x0, y0, z0))
        if v == 0:
            p2 = w0
            q2 = z0
            p1 = w0+x0+y0*(2**int(gamma*modulus_bit_length+beta2*modulus_bit_length))
            q1 = N1/p1
            if p1 is not None and q1 is not None and p1 * q1 == N1 and p2 is not None and q2 is not None and p2 * q2 == N2:
                logging.info(f"Succeeded!")
                logging.debug(f"Found p1 = {p1}")
                logging.debug(f"Found q1 = {q1}")
                logging.debug(f"Found p2 = {p2}")
                logging.debug(f"Found q2 = {q2}")
                return 1
    else:
        logging.info(f"Failed!")
        return 0

if __name__ == "__main__":
    seed = None

    # Check if command-line arguments are provided
    if len(sys.argv) != 7:
        print("Usage: sage gifp.sage <modulus_bit_length> <alpha> <gamma> <beta1> <beta2> <m>")
        sys.exit("Sorry, try again.")

    # Parse command-line arguments
    modulus_bit_length = int(sys.argv[1])
    alpha, gamma, beta1, beta2 = map(RR, sys.argv[2:6])
    m = int(sys.argv[-1])
    result = attack_gifp_instance(modulus_bit_length, alpha, gamma, beta1, beta2, m, seed)
    print(result)
