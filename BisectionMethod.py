from sympy import * # Module to work with symbolic variables
from sympy.abc import _clash1 # So that sympify will be able to convert string to symbolic objects
from statistics import mean # Average method

# Receives a function as a string and evaluate it at a chosen value
def funEvaluate(function, variable, value: float):
    try: # Try and convert the function from string to symbolic
        function = sympify(function)
    except: # If it fails, print error message
        raise Exception("The entered function string was formatted wrongly. Try something like: '3*x*exp(-7*x**2)'.")

    try: # Same goes for for the variable
        variable = sympify(variable)
    except:
        raise Exception("The entered variable is invalid.")

    # Actually evaluate the function and return the result
    return (function.subs(variable,value)).evalf()

# Receives f(x) and searches for a root in the interval [a,b] with an error of epsilon. Goes as far as MAXITER iterations
def bisection(function_x, a_0: float, b_0: float, epsilon: float, MAXITER: int):
    try: # Try and convert the function from string to symbolic
        function_x = sympify(function_x)
    except: # If it fails, print error message
        raise Exception("The entered string was formatted wrongly. Try something like: '3*x*exp(-7*x**2)'.")

    def f(n): # Defined function to simplify the notation used to evaluate the received function (function_x)
        return funEvaluate(function_x,"x",n)

    # The intervals will be stored in lists, each iteration's being an element
    a = [a_0]
    b = [b_0]

    # Check if the interval's extremes are roots
    if f(a[0]) == 0:
        return a[0]

    elif f(b[0]) == 0:
        return b[0]

    # The method only works if it's guaranteed there's at least one root in the specified interval,
    # and if this condition is met, there is no guarantee
    elif f(a[0]) * f(b[0]) > 0:
        raise Exception("Not able to search root in the specified interval.")

    x = [a[0]]

    # Loops through the algorithm until the maximum number of iterations is met
    for i in range(1,MAXITER):
        x.append(mean([a[i-1],b[i-1]]))

        # Conditions to check if it found either the actual root or one within the specified error range, respectively
        if ( f(x[i]) == 0 ) or ( abs(x[i] - x[i-1]) < epsilon * max(1,abs(x[i])) ):
            return x[i]

        # Checks in which new subinterval is guaranteed to have a root
        elif f(a[i-1]) * f(x[i]) < 0:
            a.append(a[i-1])
            b.append(x[i])

        elif f(x[i]) * f(b[i-1]) < 0:
            a.append(x[i])
            b.append(b[i-1])

    # If it manages to get out of the loop, it means no root was found
    # So it prints an error message
    raise Exception("The maximum number of iterations met and no root was found.")

if __name__ == "__main__":
    fun = input("Enter function f(x) = ")
    print()

    print("Considering the interval [a,b], enter:")
    a = (float)(sympify(input("a = ")))
    b = (float)(sympify(input("b = ")))
    print()

    error = (float)(sympify(input("Enter the desired error: ")))
    maxIterations = (int) (sympify(input("Enter the maximum number os iterations desired: ")))
    print()

    root = bisection(fun,a,b,error,maxIterations)

    print("A ROOT WAS FOUND SUCCESSFULLY:")
    print("x = %f\n" % root)