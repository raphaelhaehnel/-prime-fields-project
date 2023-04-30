from typing import List
import numpy as np
from finiteField import FiniteField


class FiniteFieldElement:
    """
    This class represents an element 'alpha' from the field 'l' = k[x]/<f(x)>
    """

    def __init__(self, coeffs: List[int], field: FiniteField):
        self.coeffs = np.array(coeffs)
        self.field = field
        self.n = len(self.coeffs) - 1

        if self.field.n != self.n + 1:
            error = f"The degree of f(x) must be equal to the degree of the element - 1"
            raise ValueError(error)

        self.matrix = self.to_matrix()

    def __add__(self, other):
        """
        TODO: must be replaced by the matrix operations
        """

        # Check if the other object is not a FiniteFieldElement object
        if not isinstance(other, FiniteFieldElement):
            error = f"Cannot add with non-FiniteFieldElement objects"
            raise ValueError(error)

        if self.field != other.field:
            error = f"Cannot add FiniteFieldElement objects from different fields"
            raise ValueError(error)

        new_coeff = (self.coeffs + other.coeffs) % self.field.p

        return self.__class__(new_coeff, self.field)

    def to_matrix(self):

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

    def isvalid(self, other):
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

    def __add__(self, other):

        if not self.isvalid(other):
            return

        result = (self.matrix + other.matrix) % self.field.p
        output = self.__class__(result[0, :], self.field)

        return output

    def __sub__(self, other):

        if not self.isvalid(other):
            return

        result = (self.matrix - other.matrix) % self.field.p
        output = self.__class__(result[0, :], self.field)

        return output

    def __mul__(self, other):

        if not self.isvalid(other):
            return

        result = (self.matrix @ other.matrix) % self.field.p
        output = self.__class__(result[0, :], self.field)

        return output

    def __truediv__(self, other):
        pass

    def __pow__(self, other):
        pass

    def mult_order(self):
        multiplicator = self.matrix

        i = 0
        while not np.all(multiplicator == self.field.identity):
            multiplicator = (multiplicator @ self.matrix) % self.field.p
            i += 1

            if i >= self.field.p ** self.field.n:
                error = f"The multiplication order is too big (> {self.field.p ** self.field.n} = self.field.p ** self.field.n)"
                raise TypeError(error)

        return i


if __name__ == "__main__":

    # Exemple 1
    ff1 = FiniteField(47, [42, 3, 0, 1])
    a = FiniteFieldElement([5, 10, 15], ff1)
    matrix1 = a.to_matrix()

    # Exemple 2
    ff2 = FiniteField(3, [2, 0, 0, 2, 1])
    b = FiniteFieldElement([1, 2, 0, 1], ff2)
    matrix2 = b.to_matrix()

    # Exemple 3
    c = FiniteFieldElement([3, 10, 1], ff1)
    d = FiniteFieldElement([0, 45, 40], ff1)

    result = c + d
    a.mult_order()
