from sympy.abc import _clash1 # So that sympify will be able to convert string to symbolic objects 
from statistics import mean # Average method 
from math import log10 # Log base 10 
from sympy import * # Module to work with symbolic variables 
from util import * # Auxiliar function in  util.py

# LINEAR ITERATIONS METHOD
# Receives all functinos and variables as lists, for it works for any number of them
def linearIterations(F_var: list, variables: list, var_0: list, epsilon = 10**-6, MAXITER = 100, displayConsoleLog = True, printTotalIterations = True):
    try: # Try and convert the functions from string to symbolic
        for func in F_var: func = sympify(func)
    except: # If it fails, print error message
        raise Exception("The entered strings (F_var) were formatted wrongly. Try something like: '3*x*exp(-7*y**2)'.")

    try: # Check if all initial approximations are floats
        for var in var_0: var = float(sympify(var).evalf())
    except: # If not, print error message
        raise Exception("The entered initial approximations (var_0) are invalid - must be floats.")

    try: # Check if all variables are valid
        for var in variables: var = sympify(var)
    except: raise Exception("The entered variables are invalid.")

    # Checks if there are as many variables as equations
    if len(F_var) != len(var_0) or len(F_var) != len(variables): raise Exception("F_x and x_0 must have the same length.")
    else: order = len(F_var) # System order (number of equations and variables)

    # Defined function to simplify the notation used to evaluate the received function F_var[i]
    def F(i: int, values: list):
        fun = F_var[i]
        for i in range(order): fun = funEvaluate(fun, variables[i], values[i])

        return float(fun)

    # No. of decimal places for the results
    precision = int(abs(log10(epsilon))) if epsilon > 0 and log10(epsilon).is_integer else int(6)

    # Only 2 iterations' values will be stored for each variable
    varIt = [[val,0] for val in var_0] # Variables iterations
    varValues = lambda e : [item[e] for item in varIt] # Get e (= j for current, = j-1 for previous) iteration value for all variables

    if displayConsoleLog: 
        printTotalIterations = True
        print("___________________________________________________________________________\n\n")
        print(f"LINEAR ITERATIONS METHOD ::: Iteration no.00:", end=" ")
        for k in range(order): print(f"F_{str(k).zfill(2)} = % .{precision-1}e" % F(k,varValues(0)), end=" ; ")
        for k in range(order):
            print(f"{variables[k]} = % .{precision}f" % varIt[k][0], end="")
            print(" ; ", end="") if k < order - 1 else print()

    # Auxiliar index variable correspondent to the current element in the variables lists
    j = 0

    # Loops through the algorithm until the maximum number of iterations is met
    for i in range(1,MAXITER):
        j += 1 if j < 1 else -j # j = 0 --> j = 1 --> j = 0 --> j = 1 ...

        # Calculates every variable's value in the current iteration
        for k in range(order): varIt[k][j] = F(k,varValues(j-1)) 

        if displayConsoleLog: 
            print(f"LINEAR ITERATIONS METHOD ::: Iteration no.{str(i).zfill(2)}:", end=" ")
            for k in range(order): print(f"F_{str(k).zfill(2)} = % .{precision-1}e" % F(k,varValues(j)), end=" ; ")
            for k in range(order):
                print(f"{variables[k]} = % .{precision}f" % varIt[k][j], end="")
                print(" ; ", end="") if k < order - 1 else print()


        # Conditions to check if it found either the actual solution or one within the specified error range, respectively
        isSolution = [( abs(var[j] - var[j-1]) < epsilon * max(1,abs(var[j])) ) for var in varIt]
        if isSolution.count(False) == 0:
            if printTotalIterations: 
                print("___________________________________________________________________________\n\n") if not displayConsoleLog else print("\n")
                print("LINEAR ITERATIONS METHOD ::: The solution was successfully found after %d iterations.\n" % (i+1) )
                print("___________________________________________________________________________\n\n")
            return varValues(j), precision

    # If it manages to get out of the loop, it means no root was found, so it prints an error message
    raise Exception("The maximum number of iterations was met and no solution was found.")

# NEWTON'S METHOD
def newtons(f_var: list, variables: list, var_0: list, epsilon = 10**-6, MAXITER = 100, displayConsoleLog = True, printTotalIterations = True):
    try: # Try and convert the functions from string to symbolic
        for func in f_var: func = sympify(func)
    except: # If it fails, print error message
        raise Exception("The entered strings (f_var) were formatted wrongly. Try something like: '3*x*exp(-7*y**2)'.")

    try: # Check if all initial approximations are floats
        for var in var_0: var = float(sympify(var).evalf())
    except: # If not, print error message
        raise Exception("The entered initial approximations (var_0) are invalid - must be floats.")

    try: # Check if all variables are valid
        for var in variables: var = sympify(var)
    except: raise Exception("The entered variables are invalid.")

    # Checks if there are as many variables as equations
    if len(f_var) != len(var_0) or len(f_var) != len(variables): 
        raise Exception("F_x and x_0 must have the same length.")
    else: order = len(f_var) # System order (number of equations and variables)

    # Jacobian Matrix
    J = Matrix([[diff(func,var) for var in variables] for func in f_var])
    # psi(x) analogue
    F_var = Matrix(variables)-J**-1*Matrix(f_var)

    if displayConsoleLog: 
        print("___________________________________________________________________________\n\n")
        print("NEWTON'S METHOD ::: The iterative formulas were calculated successfully: ")
        pprint(F_var)
        print()

    # Defined function to simplify the notation used to evaluate the received function F_var[i]
    def F(i: int, values: list):
        fun = F_var[i]
        for i in range(order): fun = funEvaluate(fun, variables[i], values[i])

        return float(fun)

    def f(i: int, values: list):
        fun = f_var[i]
        for i in range(order):
            fun = funEvaluate(fun, variables[i], values[i])

        return float(fun)

    # No. of decimal places for the results
    precision = int(abs(log10(epsilon))) if epsilon > 0 and log10(epsilon).is_integer else int(6)

    # Only 2 iterations' values will be stored for each variable
    varIt = [[val,0] for val in var_0] # Variables iterations
    varValues = lambda e : [item[e] for item in varIt] # Get e (= j for current, = j-1 for previous) iteration value for all variables

    if displayConsoleLog: 
        printTotalIterations = True
        print(f"NEWTON'S METHOD ::: Iteration no.00:", end=" ")
        for k in range(order): print(f"f_{str(k).zfill(2)} = % .{precision-1}e" % f(k,varValues(0)), end=" ; ")
        for k in range(order):
            print(f"{variables[k]} = % .{precision}f" % varIt[k][0], end="")
            print(" ; ", end="") if k < order - 1 else print()

    # Auxiliar index variable correspondent to the current element in the 'x' and 'y' lists
    j = 0

    # Loops through the algorithm until the maximum number of iterations is met
    for i in range(1,MAXITER):
        j += 1 if j < 1 else -j # j = 0 --> j = 1 --> j = 0 --> j = 1 ...

        # Calculates every variable's value in the current iteration
        for k in range(order): varIt[k][j] = F(k,varValues(j-1)) 

        if displayConsoleLog: 
            print(f"NEWTON'S METHOD ::: Iteration no.{str(i).zfill(2)}:", end=" ")
            for k in range(order): print(f"f_{str(k).zfill(2)} = % .{precision-1}e" % F(k,varValues(j)), end=" ; ")
            for k in range(order):
                print(f"{variables[k]} = % .{precision}f" % varIt[k][j], end="")
                print(" ; ", end="") if k < order - 1 else print()

        # Conditions to check if it found either the actual solution or one within the specified error range, respectively
        isSolution = [( abs(var[j] - var[j-1]) < epsilon * max(1,abs(var[j])) ) for var in varIt]
        if isSolution.count(False) == 0:
            if printTotalIterations is True: 
                print("___________________________________________________________________________\n\n") if not displayConsoleLog else print("\n")
                print("NEWTON'S METHOD ::: The solution was successfully found after %d iterations.\n" % (i+1) )
                print("___________________________________________________________________________\n\n")
            return varValues(j), precision

    # If it manages to get out of the loop, it means no root was found, so it prints an error message
    raise Exception("The maximum number of iterations was met and no solution was found.")


if __name__ == "__main__":
    f_xy = ["x**2 - y - 1", "y**2 - x - 1"]
    x_0 = 2
    y_0 = 2
    variables = ["x", "y"]
    var_0 = [x_0, y_0]
    error = 10**-5
    maxIt = 10**3
    discriminateIt = True
    totalIt = True

    sol, precision = newtons(f_xy,variables,var_0,error,maxIt,discriminateIt,False,totalIt) # NEWTON
    print("SOLUTION:")
    print(f"x = {sol[0]:.{precision}f}")
    print(f"y = {sol[1]:.{precision}f}\n")

    F_xy = ["sqrt(y+1)", "sqrt(x+1)"]

    sol, precision = linearIterations(F_xy,variables,var_0,error,maxIt,discriminateIt,totalIt) # LINEAR ITERATIONS
    print("SOLUTION:")
    print(f"x = {sol[0]:.{precision}f}")
    print(f"y = {sol[1]:.{precision}f}\n")
