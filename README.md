# Tripartite Diffie Hellman

#### Cryptography and Security Protocols Project
This is our Cryptography project where we implemented the Tripartite Diffie Hellman protocol proposed by Antoine Joux in *A one round protocol for tripartite Diffie-Hellman. Journal of Cryptology, 17(4):263–276, September 2004*.

We implemented it in Python, and somes of our code is based on other libraries (for example the code to generate field elements of a given Galois Field).

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

## Contributors
- Gabriel Figueira
- Miguel Grilo
