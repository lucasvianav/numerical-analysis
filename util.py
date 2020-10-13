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


# Substitutes a matrix's row or column for another
def subMatrix(matrix: Matrix, sub: str, index: int, new: Matrix):
    if sub == "row":
        matrix.row_del(index)
        matrix = matrix.row_insert(index, new)

    elif sub == "col":
        matrix.col_del(index)
        matrix = matrix.col_insert(index, new)

    else:
        raise Exception("The argument 'sub' accepts only 'row' or 'col'.")

    return matrix


# Swap a matrix's rows or columns
def swapMatrix(matrix: Matrix, sub: str, i: int, j: int):
    if sub == "row":
        tmp = matrix.row(i)
        matrix = subMatrix(matrix, sub, i, matrix.row(j))
        matrix = subMatrix(matrix, sub, j, tmp)

    elif sub == "col":
        tmp = matrix.col(i)
        matrix = subMatrix(matrix, sub, i, matrix.col(j))
        matrix = subMatrix(matrix, sub, j, tmp)

    else:
        raise Exception("The argument 'sub' accepts only 'row' or 'col'.")

    return matrix


if __name__ == "__main__":
    print(intervalMaxMin("E**x**2-4*x",0,1))
    print(intervalMaxMin("E**x**2-4*x",1,2))