import ast
from typing import List
from primeFieldElement import PrimeFieldElement
from finiteField import FiniteField
from finiteFieldElement import FiniteFieldElement

"""
This the main frame of the program
"""


# class syntax
class Choice:
    PRIMEFIELDELEMENT = 1
    FINITEFIELD = 2
    FINITEFIELDELEMENT = 3


WELCOME_MESSAGE = "Welcome to the user interface of the Advanced Algebra project !"
ASK_INPUT = "Please choose the action your want to do"
ACTION_LIST = "(1) - Create a new PrimeFieldElement\n\
(2) - Create a new FiniteField\n\
(3) - Create a new FiniteFieldElement"
NOT_DEFINED = "Not defined action"


def get_user_choice() -> int:
    valid = False
    while not valid:
        try:
            action = input()
            action = int(action)
            if action not in [Choice.FINITEFIELD, Choice.FINITEFIELDELEMENT, Choice.PRIMEFIELDELEMENT]:
                raise ValueError()

            valid = True
        except:
            print(NOT_DEFINED)
    return action


def get_integer() -> int:
    valid = False
    while not valid:
        try:
            val = int(input())
            valid = True
        except:
            print("Wrong value. Try again.")
            print("Enter number:", end=" ")
    return val


def get_prime(n) -> int:
    valid = False
    while not valid:
        if n > 1:
            for i in range(2, int(n/2)+1):
                if (n % i) == 0:
                    print(n, "is not a prime number")
                break
            return n
        # If the number is less than 1, its also not a prime number.
        else:
            print(n, "is not a prime number")
        n = get_integer()


def get_list() -> list:
    valid = False

    while not valid:
        # Get user input
        user_input = input(
            "Enter the array of integers in the form '[a1, a2, ... an]': ")

        # Validate input format
        try:
            # Convert input string to a list
            input_list = ast.literal_eval(user_input)

            # Check if the input is a list
            if isinstance(input_list, list):
                # Check if all elements are integers
                if all(isinstance(element, int) for element in input_list):
                    return input_list
                else:
                    print("Input should contain only integers.")
            else:
                print("Input should be a list.")
        except (SyntaxError, ValueError):
            print("Invalid input format.")


print(WELCOME_MESSAGE)


# Define arrays to store the elements the user will create
list_PF_element: List[PrimeFieldElement] = []
list_FF: List[FiniteField] = []
list_FF_element: List[FiniteFieldElement] = []

running = True
while running:
    print(ASK_INPUT)
    print(ACTION_LIST)

    # Get the action numer from the user
    choice = get_user_choice()

    # Creation of a PrimeFieldElement object
    if choice == Choice.PRIMEFIELDELEMENT:
        print("Enter number:", end=" ")
        a = get_integer()
        print("Enter prime number:", end=" ")
        p = get_prime(get_integer())
        alpha = PrimeFieldElement(a, p)
        list_PF_element.append(alpha)

    # Creation of a FiniteField object
    elif choice == Choice.FINITEFIELD:
        print("Enter prime number:", end=" ")
        p = get_prime(get_integer())
        print("Enter the function coefficients:", end=" ")
        f_coeffs = get_list()
        alpha = FiniteField(p, f_coeffs)
        list_FF_element.append(alpha)

    # Creation of a FiniteFieldElement object
    elif choice == Choice.FINITEFIELDELEMENT:
        if not len(list_FF):
            print("No Finite Field defined: you cannot define a FiniteFieldElement")
