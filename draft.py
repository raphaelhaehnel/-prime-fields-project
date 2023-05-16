# Python program to demonstrate working of extended
# Euclidean Algorithm

# function for extended Euclidean Algorithm
def gcdExtended(a, b):
    # Base Case
    if a == 0:
        return b, 0, 1

    gcd, x1, y1 = gcdExtended(b % a, a)

    # Update x and y using results of recursive
    # call
    x = y1 - (b//a) * x1
    y = x1

    return gcd, x, y


def inverse(a, p):
    gcd, x, y = gcdExtended(a, p)

    if gcd == p:
        error = f"{a} does not have an inverse in k"
        raise ValueError(error)

    if gcd == 1:
        return x % p

    error = f"The Euclidean algorithm didn't work with this element"
    raise ValueError(error)


# Driver code
a, b = 3, 7
g, x, y = gcdExtended(a, b)
print(f"{a}^-1 ={x % b}")
