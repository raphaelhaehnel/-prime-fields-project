from typing import List
import numpy as np


class FiniteField:
    """
    This class represents the field 'l' = k[x]/<f(x)>
    """

    def __init__(self, p: int, f_coeffs: List[int]):
        """
        We assume that p is indeed prime, and f is indeed irreducible.

        Parameters
        ----------
        p : a prime number defining the prime field
        f_coeffs : the coefficients for f(x) represented as [a_0, ..., a_{n-1}]

        """

        if f_coeffs[-1] == 0:
            error = f"Last coefficient of f cannot be 0"
            raise ValueError(error)

        # Define the type of the matrices (int or float)
        self.type = "int32" if f_coeffs[-1] == 1 else "float64"

        self.p = p
        self.f_coeffs = np.array(f_coeffs, dtype=self.type)

        # We want that coeffient a_{n-1} equals 1
        if self.f_coeffs[-1] != 1:
            self.f_coeffs /= self.f_coeffs[-1]

        self.n = len(f_coeffs) - 1  # degree of f(x)
        self.root = self.__find_root()
        self.residue = self.__find_residue()
        self.identity = np.identity(self.n, dtype=self.type)

        # If the degree of the polynom is 2 or 3, we want to check its roots and raise an error if the
        # polynom has roots
        if self.n == 2 or self.n == 3:
            if self.__has_root():
                error = f"f(x) is not irreducible in k"
                raise ValueError(error)

    def __has_root(self):
        """
        Check if the function f(x) is indeed irreducible in k
        """

        # Go through all the elements in k
        for i in range(self.p):

            res = 0
            res += (self.f_coeffs[0]) + (self.f_coeffs[1] * i) \
                + (self.f_coeffs[2] * i ** 2)
            if self.n == 3:
                res += self.f_coeffs[3] * i ** 3

            if res == 0:
                return True
        return False

    def __find_residue(self):
        """
        Find the residue of the polynomial equation
        If f(x) = x^2+bx+c = 0, returns (-bx-c) % p
        """

        return (- self.f_coeffs[:-1] + self.p) % self.p

    def __eq__(self, other):
        """
        Overloading the == operator
        """
        return self.f_coeffs == other.f_coeffs and self.p == other.p

    def __ne__(self, other):
        """
        Overloading the != operator
        """
        return not (self == other)

    def __find_root(self):
        """
        Find the root of f(x) (not in k)
        """
        return np.roots(self.f_coeffs[::-1])


pass
