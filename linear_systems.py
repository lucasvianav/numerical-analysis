from sympy.abc import _clash1 # So that sympify will be able to convert string to symbolic objects 
from statistics import mean # Average method 
from math import log10 # Log base 10 
from sympy import * # Module to work with symbolic variables 
from util import * # Auxiliar function in  util.py

# EXACT METHODS _______________________________________________________________________________________

# GAUSSIAN ELIMINATION METHOD
def gaussianElimination(coeffs: Matrix, variables: Matrix, result: Matrix, displayConsoleLog = True):
    # Try and convert the variables to symbolic
    try: 
        for var in variables: var = sympify(var)
    except:
        raise Exception("The entered variables are invalid.")

    # Checks if all coefficients are numeric
    try: 
        for coef in coeffs: coef = float(coef)
    except: 
        raise Exception("The coefficients must be numeric.")

    # Checks if there are as many variables as equations
    if coeffs.shape[0] == coeffs.shape[1] == variables.shape[0] == result.shape[0] and variables.shape[1] == result.shape[1] == 1: 
        order = coeffs.shape[0] # System order (number of equations and variables)
    else:
        raise Exception("The system is invalid.")

    # Checks if it is an independent linear system
    if coeffs.det() != 0 and displayConsoleLog:
        print(f"GAUSSIAN ELIMINATION METHOD ::: The linear system is an independent one (has only 1 solution).")
    elif coeffs.det() == 0:
        if displayConsoleLog:
            print(f"GAUSSIAN ELIMINATION METHOD ::: The linear system has either no or infinite solutions.")
        
        return

    aug = coeffs.col_insert(order + 1, result) # Augmented Matrix

    if displayConsoleLog:
        print(f"GAUSSIAN ELIMINATION METHOD ::: Step no.{'0'.zfill(len(str(order)))}: augmented matrix [A|b] = ")
        pprint(aug)
        print()

    for j in range(order-1): # Columns
        while True:
            # Swaps elements if needed
            if aug.col(j)[j] == 0:
                i = j

                # Selects the first line below the current whose corresponding element isn't 0
                while aug.col(j)[i] == 0 and i < order: i += 1 

                # If no line was found, there is no solution
                if aug.col(j)[i] == 0: raise Exception("No solution was found (main diagonal element equals to 0).")

                aug = swapMatrix(aug, "row", i, j)
                if displayConsoleLog:
                    print(f"GAUSSIAN ELIMINATION METHOD ::: Step no.{str(j+1).zfill(len(str(order)))}: swap rows {j+1} and {i+1}.")

            for i in range(j+1, order): # Rows
                multiplier = aug.row(i)[j] / aug.row(j)[j]
                aug = subMatrix(aug, "row", i, aug.row(i) - multiplier*aug.row(j))

                if displayConsoleLog:
                    print(f"GAUSSIAN ELIMINATION METHOD ::: Step no.{str(j+1).zfill(len(str(order)))}: row{i+1} = row{i+1} - ({multiplier}) * row{j+1}.")

            # If there's no need for more swaping, breaks out of the while loop
            if aug.col(j)[j] != 0: break
        
        if displayConsoleLog:
            print(f"GAUSSIAN ELIMINATION METHOD ::: Step no.{str(j+1).zfill(len(str(order)))}: augmented matrix [A|b] = ")
            pprint(aug)
            print()

    result = aug.col(-1).copy()
    coeffs = aug.copy(); coeffs.col_del(-1)

    sol = solve(coeffs*variables - result, variables)
    sol = Matrix([sol.get(var) for var in variables])

    return coeffs, result, sol


# CHOLESKY METHOD
def cholesky(coeffs: Matrix, variables: Matrix, result = None, displayConsoleLog = True):
    # Try and convert the variables to symbolic
    try: 
        for var in variables: var = sympify(var)
    except:
        raise Exception("The entered variables are invalid.")

    # Checks if all coefficients are numeric
    try: 
        for coef in coeffs: coef = float(coef)
    except: 
        raise Exception("The coefficients must be numeric.")

    # Checks if there are as many variables as equations
    if coeffs.shape[0] == coeffs.shape[1] == variables.shape[0] and variables.shape[1] == 1: 
        order = coeffs.shape[0] # System order (number of equations and variables)
    else:
        raise Exception("The system is invalid.")

    # Checks if the result parameter is a matrix
    if result:
        if type(result) is not Matrix: raise Exception("The 'result' parameter must be a matrix.")
        elif result.shape[0] != order or result.shape[1] != 1: raise Exception("The system is invalid.")

    # Checks if it is an independent linear system
    if coeffs.det() != 0 and displayConsoleLog:
        print("CHOLESKY METHOD ::: The linear system is an independent one (has only 1 solution).")
    elif coeffs.det() == 0:
        if displayConsoleLog:
            print("CHOLESKY METHOD ::: The linear system has either no or infinite solutions.")
        
        return

    # Checks if the coefficient matrix is simetric
    if coeffs.is_symmetric() and displayConsoleLog:
        print("CHOLESKY METHOD ::: The coefficient matrix is symmetric.")
    elif not coeffs.is_symmetric():
        if displayConsoleLog:
            print("CHOLESKY METHOD ::: The coefficient matrix is not symmetric.")

        return

    G = []
    for i in range(order): # Rows
        row = []

        for j in range(order): # Columns
            if j < i: 
                el = (coeffs.row(i)[j] - sum([(row[k])*(G[j][k]) for k in range(j)]))/G[j][j]

            elif j == i: 
                el = sqrt(coeffs.row(j)[j] - sum([(row[k])**2 for k in range(j)]))

            else: # j > i
                el = 0

            if displayConsoleLog and j <= i:
                print(f"CHOLESKY METHOD ::: g_{i+1}{j+1} = % .4f." % el)

            row.append(el)

        G.append(row)

    G = Matrix(G)

    if displayConsoleLog:
        print(f"\nCHOLESKY METHOD ::: G = ")
        pprint(G)

        print(f"\nCHOLESKY METHOD ::: G^T = ")
        pprint(G.T)
        print()

    # If a result matrix was passed, it is possible to calculate the system's solution
    if result:
        # Solution:
        # A*x = b --> G*G_T*x = b --> G*y = b
        # G_T*x = y

        auxSol = solve(G*variables - result, variables)
        auxSol = Matrix([auxSol.get(var) for var in variables])

        if displayConsoleLog:
            print("CHOLESKY METHOD ::: A*x = b --> G*G^T*x = b --> G*y = b and G^T*x = y.")
            print("CHOLESKY METHOD ::: G*y = b --> y = ")
            pprint(auxSol)
            print()

        sol = solve(G.T*variables - auxSol, variables)
        sol = Matrix([sol.get(var) for var in variables])

        return G, sol

    else:
        return G


# ITERATIVE METHODS _______________________________________________________________________________________


if __name__ == "__main__":
    # A = Matrix([ [2, 1, 2], [1, 2, 1], [1, 1, 2] ])
    # # A = Matrix([ [1,2,0], [1,2,1], [0,1,1] ])
    # var = Matrix([ 'x', 'y', 'z' ])
    # res = Matrix([ 1, 2, 3 ])
    # # res = Matrix([ 3, 1, -4 ])

    # coeffs, result, sol = gaussianElimination(A, var, res)

    # print("\n")
    # pprint(coeffs)
    # print()
    # pprint(result)
    # print()
    # pprint(sol)
    # print()


    # A = Matrix([ [4,2,2], [2,2,1], [2,1,2] ])
    A = Matrix([ [4, 12, -16], [12, 37, -43], [-16, -43, 98] ])
    var = Matrix([ 'x', 'y', 'z' ])

    G = cholesky(A, var)