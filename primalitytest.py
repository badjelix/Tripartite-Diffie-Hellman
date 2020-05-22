# Python3 program to implement Solovay-Strassen
# Primality Test
import random

def solovay_strassen(n, k=10):
	if n == 2:
		return True
	if not n & 1:
		return False

	def legendre(a, p):
		if p < 2:
			raise ValueError('p must not be < 2')
		if (a == 0) or (a == 1):
			return a
		if a % 2 == 0:
			r = legendre(a / 2, p)
			if p * p - 1 & 8 != 0:
				r *= -1
		else:
			r = legendre(p % a, a)
			if (a - 1) * (p - 1) & 4 != 0:
				r *= -1
		return r

	for i in xrange(k):
		a = randrange(2, n - 1)
		x = legendre(a, n)
		y = pow(a, (n - 1) / 2, n)
		if (x == 0) or (y != x % n):
			return False

	return True

# Driver Code
iterations = 50;
num2 = 2**224 * (2**32 - 1) + 2**192 + 2**96 - 1;
num3 = 2**256 - 2**32 - 2**9 - 2**8 - 2**7 - 2**6 - 2**4 - 1;
num = 5210644015679228794060694325390955853397149309953825381775591280356090833797121
num5 = 482677778157700435350444108563600470389539607291135742953085077414483299007817968457323051999107203153032937333023591271636050696817523671646492380723773419011;

if (solovay_strassen(num, iterations)):
    print(num, "is prime ");
else:
    print(num, "is composite");


# This code is contributed by mits
