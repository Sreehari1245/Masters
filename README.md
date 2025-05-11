# Masters
# NTRU CPA Encryption & LLL Reduction

This repository contains a Python implementation of:

- A CPA-secure NTRU, which is a post-quantum lattice-based encryption scheme.
- The Lenstra–Lenstra–Lovász (LLL) lattice basis reduction algorithm.

The code and mathematical foundations follow the presentation in the textbook:

"An Introduction to Mathematical Cryptography"_  
by Jeffrey Hoffstein, Jill Pipher, and Joseph H. Silverman.

---

1. NTRU Encryption (CPA-secure)
- Key generation using truncated polynomial rings.
- Message encryption and decryption using public and private polynomials.
- Cleanly written for clarity and modularity.

2. LLL Reduction Algorithm
- Implementation of the LLL basis reduction algorithm.
- Uses Gram-Schmidt orthogonalization and size-reduction.

3. Examples

- Small toy example to check correctness.
- Random full-rank 10-dimensional basis generation using NumPy.
- Integration with lattice-based cryptography via NTRU.

---



Requirements

- Python 3.x
- NumPy


