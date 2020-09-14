from sympy.abc import _clash1 # So that sympify will be able to convert string to symbolic objects 
from statistics import mean # Average method 
from math import log10 # Log base 10 
from sympy import * # Module to work with symbolic variables 
from util import * # Auxiliar function in  util.py 

# BISECTION METHOD
# Receives f(x) and searches for a root in the interval [a,b], with an error margin of epsilon. Goes as far as MAXITER iterations
def bisection(function_x, a: float, b: float, epsilon = 10**-6, MAXITER = 100, displayConsoleLog = False, printTotalIterations = False):

    try: # Try and convert the function from string to symbolic
        function_x = sympify(function_x)
    except: # If it fails, print error message
        raise Exception("The entered string was formatted wrongly. Try something like: '3*x*exp(-7*x**2)'.")

    # Defined function to simplify the notation used to evaluate the received function (function_x)
    f = lambda n : funEvaluate(function_x,"x",n)

    if displayConsoleLog: printTotalIterations = True

    # No. of decimal places for the results
    precision = int(abs(log10(epsilon))) if epsilon > 0 and log10(epsilon).is_integer else 6

    # Check if the interval's extremes are roots
    if f(a) == 0:
        if printTotalIterations: 
            print("___________________________________________________________________________\n\n")
            print("\nBISECTION METHOD ::: No iterations were necessary, the root was successfuly found.\n")
            print("___________________________________________________________________________\n\n")
        return float(a), int(precision)

    elif f(b) == 0:
        if printTotalIterations: 
            print("___________________________________________________________________________\n\n")
            print("\nBISECTION METHOD ::: No iterations were necessary, the root was successfully found.\n")
            print("___________________________________________________________________________\n\n")
        return float(b), int(precision)

    # The method only works if it's guaranteed there's at least one root in the specified interval,
    # and if this condition is met, there is no guarantee
    elif f(a) * f(b) > 0:
        raise Exception("Not able to search root in the specified interval.")

    # Only 2 elements will be used, [0] for even iterations and [1] for odd ones
    x = [a,0]

    if displayConsoleLog:
        print("___________________________________________________________________________\n\n")

    # Loops through the algorithm until the maximum number of iterations is met
    for i in range(1,MAXITER+1):
        # Auxiliar index variables in order to access elements from the 'x' list
        if isOdd(i):
            j = 1 # Indexes the current x element
            k = 0 # Indexes the previous one
        else:
            j = 0
            k = 1

        x[j] = mean([a,b])

        if displayConsoleLog is True: print(f"BISECTION METHOD ::: Iteration no.{str(i).zfill(2)}: a_{str(i-1).zfill(2)} = {a:.{precision}f} ; b_{str(i-1).zfill(2)} = {b:.{precision}f} ; x_{str(i).zfill(2)} = {x[j]:.{precision}f} ; f(a_{str(i-1).zfill(2)}) = {f(a):.{precision}f} ; f(b_{str(i-1).zfill(2)}) = {f(b):.{precision}f}; f(x_{str(i).zfill(2)}) = {f(x[j]):.{precision}f}")

        # Conditions to check if it found either the actual root or one within the specified error range, respectively
        if ( f(x[j]) == 0 ) or ( abs(x[j] - x[k]) < epsilon * max(1,abs(x[j])) ):
            if printTotalIterations is True: 
                print("___________________________________________________________________________\n\n") if not displayConsoleLog else print("\n")
                print("BISECTION METHOD ::: The root was successfully found after %d iterations.\n" % i )
                print("___________________________________________________________________________\n\n")
            return float(x[j]), int(precision)

        # Checks in which new subinterval is guaranteed to have a root
        elif f(a) * f(x[j]) < 0:
            b = x[j]

        elif f(x[j]) * f(b) < 0:
            a = x[j]

    # If it manages to get out of the loop, it means no root was found, so it prints an error message
    raise Exception("The maximum number of iterations was met and no root was found.")


# LINEAR ITERATIONS METHOD
# Receives psi(x) (= x) and searches for f(x)'s root starting at a point x_0, with an error margin of epsilon. Goes as far as MAXITER iterations
def linearIterations(psi_x, x_0: float, epsilon = 10**-6, MAXITER = 100, confirmConvergence = False, infLim = 0.0, supLim = 0.0, displayConsoleLog = False, printTotalIterations = False, methodName = "LINEAR ITERATIONS"):

    try: # Try and convert the function from string to symbolic
        psi_x = sympify(psi_x)
    except: # If it fails, print error message
        raise Exception("The entered string was formatted wrongly. Try something like: '3*x*exp(-7*x**2)'.")

    if displayConsoleLog and methodName == "LINEAR ITERATIONS": 
        print("___________________________________________________________________________\n\n")

    if confirmConvergence: # Checks if the passed psi(x) meets |psi'(x)| < 1, x in [infLim, supLim] and if x_0 is in the interval
        if x_0 >= infLim and x_0 <= supLim: # Checks if x_0 is in the interval
            if displayConsoleLog: print(methodName + " METHOD ::: Confirmed x_0 is in [%.2f, %.2f] successfully." % (infLim, supLim))
        else: 
            if displayConsoleLog: print(methodName + " METHOD ::: x_0 is not in [%.2f, %.2f].\n\n" % (infLim, supLim))
            raise Exception("The chosen x_0 doesn't meet the necessary conditions ar the specified interval.")

        dpsi = diff(psi_x,'x') # Calculates psi'(x)
        if displayConsoleLog: print(methodName + " METHOD ::: Derivative function psi'(x) was calculated successfully.")

        maxValue, minValue = intervalMaxMin(dpsi,infLim,supLim,'x') # Calculates psi'(x) maximum and minimum values in the passed interval
        if displayConsoleLog: print(methodName + " METHOD ::: psi'(x)'s maximum (%.2f) and minimum (%.2f) values were calculated successfully." % (maxValue, minValue))

        if (maxValue > 1.0) or (minValue < -1.0): # Checks if the inequality is true
            if displayConsoleLog: print(methodName + " METHOD ::: |psi'(x)| !< 1, x in [%.2f,%.2f]. \n\n" % (infLim, supLim))
            raise Exception("The chosen psi(x) doesn't meet the necessary conditions at the specified interval.") # If it isn't, raise excpetion
        else: 
            if displayConsoleLog: print(methodName + " METHOD ::: This psi(x) meets the necessary conditions.\n")

    # Defined function to simplify the notation used to evaluate the received function (psi_x)
    psi = lambda n : funEvaluate(psi_x,"x",n)

    if displayConsoleLog: printTotalIterations = True

    # No. of decimal places for the results
    precision = int(abs(log10(epsilon))) if epsilon > 0 and log10(epsilon).is_integer else 6

    # Only 2 elements will be used, [0] for even iterations and [1] for odd ones
    x = [x_0,0]

    # Loops through the algorithm until the maximum number of iterations is met
    for i in range(1,MAXITER+1):
        # Auxiliar index variables in order to access elements from the 'x' list
        if isOdd(i):
            j = 1 # Indexes the current x element
            k = 0 # Indexes the previous one
        else:
            j = 0
            k = 1

        x[j] = psi(x[k])

        if displayConsoleLog is True: print(methodName + f" METHOD ::: Iteration no.{str(i).zfill(2)}: psi(x_{str(i).zfill(2)}) = {psi(x[j]):.{precision}f} ; x_{str(i).zfill(2)} = {x[j]:.{precision}f}")

        # Conditions to check if it found either the actual root or one within the specified error range, respectively
        if  abs(x[j] - x[k]) < epsilon * max(1,abs(x[j])):
            if printTotalIterations is True: 
                print("___________________________________________________________________________\n\n") if not displayConsoleLog else print("\n")
                print(methodName + " METHOD ::: The root was successfully found after %d iterations.\n" % i )
                print("___________________________________________________________________________\n\n")
            return float(x[j]), int(precision)

    # If it manages to get out of the loop, it means no root was found, so it prints an error message
    raise Exception("The maximum number of iterations was met and no root was found.")


# NEWTON'S METHOD
# Receives f(x), calculates the best possible psi(x) and then executes the Linear Iterations Method
def newtons(function_x, x_0: float, epsilon = 10**-6, MAXITER = 100, displayConsoleLog = False, printTotalIterations = False):

    try: # Try and convert the function from string to symbolic
        function_x = sympify(function_x)
    except: # If it fails, print error message
        raise Exception("The entered string was formatted wrongly. Try something like: '3*x*exp(-7*x**2)'.")

    dfun = diff(function_x,'x')

    psi_x = simplify(sympify('x') - (function_x/dfun))
    if displayConsoleLog: 
        print("___________________________________________________________________________\n\n")
        print("NEWTON'S METHOD ::: psi(x) was calculated successfully: %s.\n" % str(psi_x))

    return linearIterations(psi_x, x_0, epsilon, MAXITER, confirmConvergence=False, displayConsoleLog=displayConsoleLog, printTotalIterations=printTotalIterations, methodName="NEWTON'S")


# SECANT METHOD
# Receives f(x), estimates an expression for x_(k+1) and searches for f(x)'s root starting at a points x_0 and x_1, with an error margin of epsilon. Goes as far as MAXITER iterations
def secant(function_x, x_0: float, x_1: float, epsilon = 10**-6, MAXITER = 100, displayConsoleLog = False, printTotalIterations = False):
    try: # Try and convert the function from string to symbolic
        function_x = sympify(function_x)
    except: # If it fails, print error message
        raise Exception("The entered string was formatted wrongly. Try something like: '3*x*exp(-7*x**2)'.")

    if displayConsoleLog: 
        print("___________________________________________________________________________\n\n")

    # Defined function to simplify the notation used to evaluate the received function (psi_x)
    f = lambda n : funEvaluate(function_x,"x",n)

    if displayConsoleLog: printTotalIterations = True

    # No. of decimal places for the results
    precision = int(abs(log10(epsilon))) if epsilon > 0 and log10(epsilon).is_integer else 6

    # Only 3 elements will be used, [0] for even iterations and [1] for odd ones
    x = [x_0,x_1,0]

    # Auxiliar index variable in order to access elements from the 'x' list
    k = 1

    # Loops through the algorithm until the maximum number of iterations is met
    for i in range(1,MAXITER+1):
        k += 1 if k < 2 else -k # k = 1 --> k = 2 --> k = 0 --> k = 1 (...)

        x[k] = (f(x[k-1])*x[k-2] - f(x[k-2])*x[k-1])/(f(x[k-1]) - f(x[k-2]))

        if displayConsoleLog is True: print(f"SECANT METHOD ::: Iteration no.{str(i).zfill(2)}: f(x_{str(i+1).zfill(2)}) = {f(x[k]):.{precision}f} ; x_{str(i+1).zfill(2)} = {x[k]:.{precision}f}")

        # Conditions to check if it found either the actual root or one within the specified error range, respectively
        if  abs(x[k] - x[k-1]) < epsilon * max(1,abs(x[k])):
            if printTotalIterations is True: 
                print("___________________________________________________________________________\n\n") if not displayConsoleLog else print("\n")
                print("SECANT METHOD ::: The root was successfully found after %d iterations.\n" % i )
                print("___________________________________________________________________________\n\n")
            return float(x[k]), int(precision)

    # If it manages to get out of the loop, it means no root was found, so it prints an error message
    raise Exception("The maximum number of iterations was met and no root was found.")


# REGULA FALSI METHOD
# Mix between Bisection and Secant
def regulaFalsi(function_x, x_0: float, x_1: float, epsilon = 10**-6, MAXITER = 100, displayConsoleLog = False, printTotalIterations = False):
    try: # Try and convert the function from string to symbolic
        function_x = sympify(function_x)
    except: # If it fails, print error message
        raise Exception("The entered string was formatted wrongly. Try something like: '3*x*exp(-7*x**2)'.")

    if displayConsoleLog: 
        print("___________________________________________________________________________\n\n")

    # Defined function to simplify the notation used to evaluate the received function (psi_x)
    f = lambda n : funEvaluate(function_x,"x",n)

    if displayConsoleLog: printTotalIterations = True

    if f(x_0) * f(x_1) >= 0: raise Exception("x_0 and x_1 are not ideal.")

    # No. of decimal places for the results
    precision = int(abs(log10(epsilon))) if epsilon > 0 and log10(epsilon).is_integer else 6

    # Only 2 elements will be used, [0] for even iterations and [1] for odd ones
    s = [x_0,x_1]

    # Loops through the algorithm until the maximum number of iterations is met
    for i in range(1,MAXITER+1):
        x = (f(s[1]) * s[0] - f(s[0]) * s[1])/(f(s[1]) - f(s[0]))

        if displayConsoleLog is True: print(f"REGULA FALSI METHOD ::: Iteration no.{str(i).zfill(2)}: f(x_{str(i+2).zfill(2)}) = {f(x):.{precision}f} ; x_{str(i+2).zfill(2)} = {x:.{precision}f}")

        # Conditions to check if it found either the actual root or one within the specified error range, respectively
        if  abs(x - s[0]) < epsilon * max(1,abs(x)) or abs(x - s[1]) < epsilon * max(1,abs(x)):
            if printTotalIterations is True: 
                print("___________________________________________________________________________\n\n") if not displayConsoleLog else print("\n")
                print("REGULA FALSI METHOD ::: The root was successfully found after %d iterations.\n" % i )
                print("___________________________________________________________________________\n\n")
            return float(x), int(precision)
        else:
            if f(x) * f(s[1]) < 0:
                s[0] = s[1]
                s[1] = x
            else: s[1] = x

    # If it manages to get out of the loop, it means no root was found, so it prints an error message
    raise Exception("The maximum number of iterations was met and no root was found.")


if __name__ == "__main__":
    #fun = "E**x**2-4*x"
    fun = "x**3 + 4*x**2 - 10"
    a = 1
    b = 2
    error = 10**-9
    maxIt = 10**3
    discriminateIt = True
    totalIt = True

    root, precision = bisection(fun,a,b,error,maxIt,displayConsoleLog=discriminateIt,printTotalIterations=totalIt) # BISECTION
    print("ROOT:")
    print(f"x = {root:.{precision}f}\n")


    psi = "sqrt(10/(4+x))"
    checkCon = True
    x_0 = 1.5

    root, precision = linearIterations(psi,x_0,error,maxIt,checkCon,a,b,discriminateIt,totalIt) # LINEAR ITERATIONS
    print("ROOT:")
    print(f"x = {root:.{precision}f}\n")

    root, precision = newtons(fun,x_0,error,maxIt,discriminateIt,totalIt) # NEWTON
    print("ROOT:")
    print(f"x = {root:.{precision}f}\n")

    x_0 = 1.5
    x_1 = 1.7

    root, precision = secant(fun,x_0,x_1,error,maxIt,discriminateIt,totalIt) # SECANT
    print("ROOT:")
    print(f"x = {root:.{precision}f}\n")

    x_0 = a
    x_1 = b

    root, precision = regulaFalsi(fun,x_0,x_1,error,maxIt,discriminateIt,totalIt) # SECANT
    print("ROOT:")
    print(f"x = {root:.{precision}f}\n")