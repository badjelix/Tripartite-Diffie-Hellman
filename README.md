# Tripartite Diffie Hellman

#### Cryptography and Security Protocols Project
This is our Cryptography project where we implemented the Tripartite Diffie Hellman protocol proposed by Antoine Joux in *A one round protocol for tripartite Diffie-Hellman. Journal of Cryptology, 17(4):263–276, September 2004*.

The project was implemented in python.

## Requirements
- Python3
- *SymPy*

### How to install *SymPy*
You can easily install *SymPy* by following the instructions on
their official website ([Link here](https://www.sympy.org/en/)).

## How to run
We implemented both Weil and Tate pairings to generate the shared key. You can run the Weil and Tate pairing approaches with the following commands, respectively:
```bash
python3 runWeil.py
python3 runTate.py
```

## Implementation
Miller's algorithm -> find functions with desired divisors: [2]

Weil Pairing -> [2]

Tate Pairing -> [3]

- fieldelement.py -> implements operations between elements of a Galois Field; adaptation from https://github.com/glassnotes/PyniteFields
- util.py -> implements operations in Weierstrass elliptic curves, as well as the Weil and Tate pairings. It also tests some pairing examples (taken from various literature, like [2] and [3]) if run directly (`python3 util.py`)
- participant.py -> file describing the behaviour of a participant in a Diffie Hellman exchange
- runWeil.py & runTate.py -> implement the tripartite exchange
- test.py -> contains the example shown in [1], for validation purposes (`python3 test.py`)
- runNIST.py -> here, we attempted to use known NIST curves to execute the protocol, instead of relying on the parameters of the example in [1]. However, we came across the following problem: NIST only provides a generator P (for a specific curve), but for the pairings to work, we need P and Q, linearly independent (from each other, such that e(P,Q) != 1) and of the same order. We attempted to solve this by performing a field extension (for example, from F_p to F_p^2) and finding Q, as explained in [3], and we also attempted to use distortion maps ([2] and [4]) to generate Q.
Unfortunately, the NIST curves did not have the desired properties for perfoming these operations, and we were unable to generate Q.

## References

[1] Antoine Joux, *A one round protocol for tripartite Diffie-Hellman. Journal of Cryptology*

[2] A. E. Aftuck, *The Weil pairing on elliptic curves and its cryptographic applications*

[3] C. Costello, *Pairings for beginners*

[4] Antoine Joux, *Separating Decision Diffie–Hellman from Computational Diffie–Hellman in Cryptographic Groups*

## Contributors
- Gabriel Figueira
- Miguel Grilo
