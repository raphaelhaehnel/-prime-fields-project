from typing import List
import numpy as np
from finiteField import FiniteField
import math


class FiniteFieldElement:
    """
    This class represents an element 'alpha' from the field 'l' = k[x]/<f(x)>
    """

    def __init__(self, coeffs: List[int], field: FiniteField):
        """
        Generate an element from the extended finite field l

        Parameters
        ----------
        coeffs : the coefficients of the polynom alpha from l
        field : the extended finite field l
        """

        self.coeffs = np.array(coeffs)
        self.field = field

        # The order of the element
        self.n = len(self.coeffs) - 1

        # Check if the degree of the polynom corresponds to the degree of the field minus 1
        if self.field.n != self.n + 1:
            error = f"The degree of f(x) must be equal to the degree of the element - 1"
            raise ValueError(error)

        # Generate the matrix representing the polynom
        self.matrix = self.to_matrix()

    def __add__(self, other):
        """
        Add two elements according to the field logic
        """

        # Check if the other object is not a FiniteFieldElement object
        if not self.__isvalid(other):
            return

        # Add the coefficients together and make modulo p
        new_coeff = (self.coeffs + other.coeffs) % self.field.p

        # Returns a new object with the new coefficients and the same field
        return self.__class__(new_coeff, self.field)

    def to_matrix(self):
        """
        This method converts the polynomial form to a matrix form according to the method we learned in class
        It uses recursion to find the matrix for all degree n
        """

        # Initialize the helper matrix
        helper = np.zeros((self.n+1, 2*self.n+1), dtype=self.field.type)

        # Place the coefficients of the elements in the helper matrix
        for i in range(self.n+1):
            helper[i, i:i+self.n+1] += self.coeffs

        # Perform the mutiplications from x to x^{n-1}
        for i in range(self.n):
            helper = self.__to_matrix_helper(helper, i)

        return helper[:, :self.n+1] % self.field.p

    def __to_matrix_helper(self, helper, i):
        """
        We call this method to find the matrix corresponding to the polynom by calling this method recursively
        """

        if i < 0:
            return helper

        # Create a matrix where each line contains the coeffients of the residue
        a = np.repeat([self.field.residue], self.n+1, axis=0)

        # Extract the 'index' column from the helper matrix
        b = helper[:, self.n+1 + i]

        # Multiply the residue by the column extracted from the helper
        helper[:, i:self.n+1 + i] += a * b[:, np.newaxis]

        # Set the column to zero now we used it
        helper[:, self.n+1 + i] = 0

        # Move on the previouse column in the helper matrix
        return self.__to_matrix_helper(helper, i-1)

    def __isvalid(self, other):
        """
        Check if the other object is FiniteFieldElement, is p is the same and if f(x) is the same
        """

        if not isinstance(other, FiniteFieldElement):
            error = f"The second element is not a FiniteFieldElement object"
            raise TypeError(error)

        if self.field.p != other.field.p:
            error = f"Cannot perform the operation with two elements in different finiteField p={self.field.p} and p={other.field.p}"
            raise TypeError(error)

        if not np.all(self.field.f_coeffs == other.field.f_coeffs):
            error = f"Cannot perform the operation with two elements in different finiteField f={self.field.f_coeffs} and f={other.field.f_coeffs}"
            raise TypeError(error)

        return True

    def __sub__(self, other):
        """
        Substract two elements according to the field logic
        """

        # Check if the other object is not a FiniteFieldElement object
        if not self.__isvalid(other):
            return

        # Add the coefficients together and make modulo p
        new_coeff = (self.coeffs - other.coeffs) % self.field.p

        # Returns a new object with the new coefficients and the same field
        return self.__class__(new_coeff, self.field)

    def __mul__(self, other):
        """
        Overload the * operator
        """

        if not self.__isvalid(other):
            return

        result = (self.matrix @ other.matrix) % self.field.p
        output = self.__class__(result[0, :], self.field)

        return output

    def __truediv__(self, other):
        """
        Overload the / operator
        """

        if not self.__isvalid(other):
            return

        result = self.matrix @ other.inverse().matrix % self.field.p
        return self.__class__(result[0, :], self.field)

    def inverse(self):
        """
        Computes the inverse of the polynom by looking for the multiplication order of this polynom
        """

        inv = self ** (self.mult_order() - 1)
        return inv

    def __pow__(self, other):
        """
        Overload the ** operator
        """

        if not isinstance(other, int):
            error = f"The exponent has to be int (instead of {other.type})"
            raise TypeError(error)

        # If the exponent is negative, computes the inverse of the polynom and then computes with the positive exponent
        if other < 0:
            return self.inverse() ** (-other)

        # If the exponent is zero, returns the identity element
        if other == 0:
            return self.__class__(self.field.identity[0, :], self.field)

        # Performs the multiplications
        multiplicator = self.matrix
        for i in range(other - 1):
            multiplicator = (multiplicator @ self.matrix) % self.field.p

        return self.__class__(multiplicator[0, :], self.field)

    def mult_order(self):
        """
        Computes the multiplicative order of the element
        """

        multiplicator = self.matrix
        i = 1

        # Multiply the element with itself until it is equal to 1 (identity matrix)
        while not np.all(multiplicator == self.field.identity):

            multiplicator = (multiplicator @ self.matrix) % self.field.p
            i += 1

            if i > self.field.p ** self.field.n:
                error = f"The multiplication order is too big (> {self.field.p ** self.field.n} = self.field.p ** self.field.n)"
                raise TypeError(error)

        return i

    def __repr__(self):
        """
        The object is represented by the string of the form "ax^2+bx+c"
        """

        representation = ""
        for i in range(self.n+1):
            coeff = "" if self.coeffs[-1-i] == 1 else self.coeffs[-1-i]
            if self.coeffs[-1-i] == 0:
                continue
            if self.n-i == 1:
                representation += f"{coeff}x + "
            elif self.n-i == 0:
                representation += f"{self.coeffs[-1-i]} + "
            else:
                representation += f"{coeff}x^{self.n-i} + "

        if representation[-2:] == "+ ":
            representation = representation[:-2]

        return representation

    def __hash__(self):
        """
        Generate a hash code to represent our element (for the hashing table)
        """

        return hash(str(self.coeffs.tolist()))

    def __eq__(self, other):
        """
        Defines how to compare two elements
        """

        # The coefficients have to be the same, and the fields have to be the same
        return all(self.coeffs == other.coeffs) and self.field == other.field


def BSGS(g: FiniteFieldElement, h: FiniteFieldElement):
    """
    This method is using the BSGS algorithm to solve the discrete logarithm problem
    """

    m = math.ceil(math.sqrt(g.field.p ** g.field.n - 1))

    hash_table = {}
    iterator = FiniteFieldElement(g.field.identity[0, :], g.field)

    for i in range(m):  # Giant steps
        hash_table[h * g**(-m*i)] = i

    for j in range(m):  # Baby steps
        if (iterator in hash_table):
            i = hash_table[iterator]
            return (j + i*m) % (g.field.p ** g.field.n - 1)
        iterator = iterator * g
    return None


# We are running a test code
if __name__ == "__main__":

    # Exemple 1
    ff1 = FiniteField(47, [42, 3, 0, 1])
    a = FiniteFieldElement([5, 10, 15], ff1)
    matrix1 = a.to_matrix()

    # Exemple 2
    ff2 = FiniteField(3, [1, 2, 0, 1])
    b = FiniteFieldElement([2, 2, 1], ff2)
    matrix2 = b.to_matrix()

    for i in range(3):
        for j in range(3):
            for k in range(3):
                if i == 0 and j == 0 and k == 0:
                    continue
                alpha = FiniteFieldElement([k, j, i], ff2)
                print(f"{alpha} : order {alpha.mult_order()}")

    print(f"Element generator of ff2: {ff2.multiplicative_group()}")

    # Exemple 3
    c = FiniteFieldElement([3, 10, 1], ff1)
    d = FiniteFieldElement([0, 45, 40], ff1)

    result = c + d
    a.mult_order()

    g = c
    h = c ** 4

    for i in range(27):
        print(f"BSGS(b, b ** {i}) = {BSGS(b, b ** i)}")

    print("Hello world!")
