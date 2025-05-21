import numpy as np


# Generate a full-rank integer basis for testing

def generate_full_rank_integer_basis(n, low=-10, high=10):
    while True:
        B = np.random.randint(low, high + 1, size=(n, n))
        if np.linalg.matrix_rank(B) == n:
            return B.tolist()  # convert to list of lists


# Basic vector operations

def dot(u, v):
    """Dot product of two vectors."""
    return sum(ui * vi for ui, vi in zip(u, v))

def scalar_mult(c, v):
    """Multiply a vector by a scalar."""
    return [c * vi for vi in v]

def vec_sub(u, v):
    """Subtract vector v from vector u."""
    return [ui - vi for ui, vi in zip(u, v)]



# Gram-Schmidt Orthogonalization

def gram_schmidt(B):
    """Given a basis B, returns: B_star: orthogonalized vectors, mu: projection coefficients"""
    n = len(B)
    B_star = []
    mu = [[0.0] * n for _ in range(n)]

    for i in range(n):
        b_star = B[i][:]  # make a copy of B[i]
        for j in range(i):
            denom = dot(B_star[j], B_star[j])
            mu[i][j] = dot(B[i], B_star[j]) / denom
            proj = scalar_mult(mu[i][j], B_star[j])
            b_star = vec_sub(b_star, proj)
        B_star.append(b_star)

    return B_star, mu


# Size Reduction Step

def size_reduce(B, mu, k):
    """Make coefficients mu[k][j] small (close to 0)."""
    """Subtract integer multiples of earlier vectors."""
    for j in reversed(range(k)):
        q = round(mu[k][j])  # nearest integer
        if q != 0:
            B[k] = vec_sub(B[k], scalar_mult(q, B[j]))  # adjust B[k]
            # Update mu[k] accordingly
            for i in range(j + 1):
                mu[k][i] -= q * mu[j][i]



# LLL Reduction Algorithm

def lll(B, delta=0.75):
    """
    Perform LLL basis reduction on basis B with parameter delta (0.5 < delta < 1).
    Returns the reduced basis.
    """
    n = len(B)
    B = [list(map(float, b)) for b in B]  # convert to float
    B_star, mu = gram_schmidt(B)
    k = 1

    while k < n:
        size_reduce(B, mu, k)  # reduce B[k] using earlier vectors

        B_star, mu = gram_schmidt(B)  # re-orthogonalize
        lhs = dot(B_star[k], B_star[k])
        rhs = (delta - mu[k][k - 1] ** 2) * dot(B_star[k - 1], B_star[k - 1])

        if lhs >= rhs:
            k += 1  # condition met, go to next vector
        else:
            # swap and reprocess
            B[k], B[k - 1] = B[k - 1], B[k]
            k = max(k - 1, 1)
            B_star, mu = gram_schmidt(B)

    return B



# Example usage

np.random.seed(43)

n = 10  # dimension
B = generate_full_rank_integer_basis(n)

print("Original 10D full-rank basis:")
for b in B:
    print(b)

reduced = lll(B)

print("\nLLL-reduced 10D basis:")
for b in reduced:
    print([round(x) for x in b])
