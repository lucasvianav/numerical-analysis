import numpy as np
from sympy import *
from sympy.abc import _clash1
from statistics import mean
from operator import itemgetter

def funEvaluate(function, variable, value: float):
    try:
        function = sympify(function)
    except:
        raise Exception("The entered function string was formatted wrongly. Try something like: '3*x*exp(-7*x**2)'.")

    try:
        variable = sympify(variable)
    except:
        raise Exception("The entered variable is invalid.")

    return (function.subs(variable,value)).evalf()


def bisection(function_x, a_0: float, b_0: float, epsilon: float, MAXITER: int):
    try:
        function_x = sympify(function_x)
    except:
        raise Exception("The entered string was formatted wrongly. Try something like: '3*x*exp(-7*x**2)'.")

    def f(n):
        return funEvaluate(function_x,"x",n)

    a = [a_0]
    b = [b_0]
    intervalExtremes = [a,b]

    if f(a[0]) == 0:
        return a[0]

    elif f(b[0]) == 0:
        return b[0]

    elif f(a[0]) * f(b[0]) > 0:
        raise Exception("Not able to search root in the specified interval.")

    x = [a[0]]

    for i in range(1,MAXITER):
        x.append(mean(map(itemgetter(i-1),intervalExtremes)))

        if ( f(x[i]) == 0 ) or ( abs(x[i] - x[i-1]) < epsilon * max(1,abs(x[i])) ):
            return x[i]

        elif f(a[i-1]) * f(x[i]) < 0:
            a.append(a[i-1])
            b.append(x[i])

        elif f(x[i]) * f(b[i-1]) < 0:
            a.append(x[i])
            b.append(b[i-1])

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