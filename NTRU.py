import random
from sympy import symbols, Poly, invert, GF

x = symbols('x')

"""Note that the reference for the code is An introduction to mathematical cryptography by
Hoffstein et. al """

# Generate ternary polynomials with the specified number of 1's and -1's
def generate_ternary_poly(N, num_pos, num_neg):
    coeffs = [0] * N
    indices = random.sample(range(N), num_pos + num_neg)
    for i in indices[:num_pos]:
        coeffs[i] = 1
    for i in indices[num_pos:]:
        coeffs[i] = -1
    return Poly(coeffs, x)

# Generates a random polynomial with N/3 + 1 1's and N/3 -1's
def generate_f(N):
    return generate_ternary_poly(N, N // 3 + 1, N // 3)

# Generates a random polynomial with N/3 1's and N/3 -1's
def generate_g_or_r(N):
    return generate_ternary_poly(N, N // 3, N // 3)

# ------ Polynomial Operations -------

# The next two functions generates an invertible function
def check_invertibility(f, modulus, N):
    try:
        return invert(f, Poly(x**N - 1, x), domain=GF(modulus))
    except:
        return None

def generate_invertible_f(N, p, q):
    while True:
        f = generate_f(N)
        if check_invertibility(f, p, N) and check_invertibility(f, q, N):
            return f

# Given a polynomial, return polynomial (mod) q
def poly_mod_q(poly, q):
    coeffs_mod_q = [int(c) % q for c in poly.all_coeffs()]
    return Poly(coeffs_mod_q, poly.gen, modulus=q)

# Mu;tiply two polynomials mod (x^N-1) mod (modulus)
def poly_mul_mod(a, b, modulus, N):
    a_coeffs = list(reversed(a.all_coeffs()))
    b_coeffs = list(reversed(b.all_coeffs()))

    # Pad to length N
    a_coeffs =  a_coeffs + [0]*(N - len(a_coeffs))
    b_coeffs =  b_coeffs + [0]*(N - len(b_coeffs))

    result = [0] * (2 * N - 1)
    for i in range(N):
        for j in range(N):
            result[i + j] = (result[i + j] + a_coeffs[i] * b_coeffs[j]) % modulus

    # Reduce mod x^N - 1
    reduced = [0] * N
    for i in range(len(result)):
        reduced[i % N] = (reduced[i % N] + result[i]) % modulus

    return Poly(list(reversed(reduced)), x)

# --- NTRU Functions ---
def ntru_keygen(N, p, q):
    f = generate_invertible_f(N, p, q)
    # Optionally hardcode f and g here for testing:
    # f = Poly([1, 0, -1, 1, 1, 0, -1], x)
    # g = Poly([1, 0, 1, 0, -1, -1, 0], x)
    g = generate_g_or_r(N)
    f_inv_p = invert(f, Poly(x**N - 1, x), domain=GF(p))
    f_inv_q = invert(f, Poly(x**N - 1, x), domain=GF(q))
    h = poly_mul_mod(g, f_inv_q, q, N)

    return f, h, f_inv_p, f_inv_q

def ntru_encrypt(m_poly, h, N, p, q):
    r = generate_g_or_r(N)
    prh = poly_mul_mod(r, h, q, N) # prh = r*h mod q
    e = p * prh + m_poly
    return poly_mod_q(e, q)


N, p, q = 7, 3, 41

# Key Generation
f, h, f_inv_p, f_inv_q = ntru_keygen(N, p, q)

# Encrypt a message
m_poly = Poly([0, 2, 0, 1, 1, 2, 1], x)
ciphertext = ntru_encrypt(m_poly, h, N, p, q)

print("Encrypted message:", ciphertext)

def center_lift(poly, q):
    """Center-lift polynomial coefficients to range [-q//2, q//2]."""
    coeffs = poly.all_coeffs()
    centered = [((c + q//2) % q) - q//2 for c in coeffs]
    return Poly(centered, poly.gen)

def ntru_decrypt(ciphertext, f, f_inv_p, N, p, q):

    a = poly_mul_mod(f, ciphertext, q, N)
    
    # Center lift coefficients
    a_centered = center_lift(a, q)

    a_mod_p = Poly([int(c) % p for c in a_centered.all_coeffs()], x, modulus=p)
    m = poly_mul_mod(f_inv_p, a_mod_p, p, N)

    return m

msg = ntru_decrypt(ciphertext, f, f_inv_p, N, p, q)


print(f"Orginal message: {m_poly} and decrypted message: {msg}")

if m_poly == msg:
    print("Success!")

