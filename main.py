"""
This the main frame of the program
"""

WELCOME_MESSAGE = "Welcome to the user interface of the Advanced Algebra project !"
ASK_INPUT = "Please choose the action your want to do"
ACTION_LIST = "(1) - Create a new PrimeFieldElement\n\
(2) - Create a new FiniteField\n\
(3) - Create a new FiniteFieldElement"
NOT_DEFINED = "Not defined action"

print(WELCOME_MESSAGE)
print(ASK_INPUT)
print(ACTION_LIST)


def get_user_choice():
    correct_action = False
    while not correct_action:
        try:
            action = input()
            action = int(action)
            if action not in [1, 2, 3]:
                raise ValueError()

            correct_action = True
        except:
            print(NOT_DEFINED)
    return action


choice = get_user_choice()
print("USER CHOICE IS " + str(choice))
