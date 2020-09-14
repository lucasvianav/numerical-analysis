from sympy import * # Module to work with symbolic variables
from sympy.abc import _clash1 # So that sympify will be able to convert string to symbolic objects
from findExtremums import findExtremums

# Checks a numbers parity
def isOdd(number: int):
    return False if number % 2 == 0 else True


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


# Return max and min values for a function in specified interval
def intervalMaxMin(func, infLim, supLim, var = 'x'):
    var = sympify(var)
    func = sympify(func)
    maxList, minList, inflection = findExtremums(func, var)

    # [infLim,supLim] --> a = func(infLim) , b = func(supLim)
    a = funEvaluate(func,var,infLim)
    b = funEvaluate(func,var,supLim)

    maxValue = max(a,b)
    minValue = min(a,b)

    for point in maxList:
        if point < supLim and point > infLim:
            aux = funEvaluate(func,var,point)
            if aux > maxValue: maxValue = aux

    for point in minList:
        if point < supLim and point > infLim:
            aux = funEvaluate(func,var,point)
            if aux < minValue: minValue = aux

    for point in inflection:
        if point < supLim and point > infLim:
            aux = funEvaluate(func,var,point)
            if aux > maxValue: maxValue = aux
            elif aux < minValue: minValue = aux

    return maxValue, minValue


if __name__ == "__main__":
    print(intervalMaxMin("E**x**2-4*x",0,1))
    print(intervalMaxMin("E**x**2-4*x",1,2))