from sympy.abc import _clash1 # So that sympify will be able to convert string to symbolic objects 
from statistics import mean # Average method 
from math import log10 # Log base 10 
from sympy import * # Module to work with symbolic variables 
from util import * # Auxiliar function in  util.py

# EXACT METHODS _______________________________________________________________________________________

# GAUSSIAN ELIMINATION METHOD
def gaussianElimination(coeffs: Matrix, result: Matrix, displayConsoleLog = True):
    # Checks if all coefficients are numeric
    try: 
        for coef in coeffs: coef = float(coef)
    except: 
        raise Exception("The coefficients must be numeric.")

    # Checks the system's order
    if coeffs.shape[0] == coeffs.shape[1] == result.shape[0] and result.shape[1] == 1: 
        order = coeffs.shape[0] # System order (number of equations and variables)
        variables = Matrix([f"x{i}" for i in range(order)])
    else:
        raise Exception("The system is invalid.")

    # Checks if it is an independent linear system
    if coeffs.det() != 0 and displayConsoleLog:
        print("___________________________________________________________________________\n\n")
        print(f"GAUSSIAN ELIMINATION METHOD ::: The linear system is an independent one (has only 1 solution).")
    elif coeffs.det() == 0:
        if displayConsoleLog:
            print("___________________________________________________________________________\n\n")
            print(f"GAUSSIAN ELIMINATION METHOD ::: The linear system has either no or infinite solutions.\n")
            print("___________________________________________________________________________\n\n")
        
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

    if displayConsoleLog: print("___________________________________________________________________________\n\n")

    result = aug.col(-1).copy()
    coeffs = aug.copy(); coeffs.col_del(-1)

    sol = solve(coeffs*variables - result, variables)
    sol = Matrix([sol.get(var) for var in variables])

    return coeffs, result, sol


# CHOLESKY'S METHOD
def cholesky(coeffs: Matrix, result = None, displayConsoleLog = True):
    # Checks if all coefficients are numeric
    try: 
        for coef in coeffs: coef = float(coef)
    except: 
        raise Exception("The coefficients must be numeric.")

    # Checks if the coefficient matrix is a square one
    if coeffs.shape[0] == coeffs.shape[1]: 
        order = coeffs.shape[0] # System order
    else:
        raise Exception("The system is invalid.")

    # Checks if the result parameter is a matrix
    if result:
        if type(result) is not Matrix: raise Exception("The 'result' parameter must be a matrix.")
        elif result.shape[0] != order or result.shape[1] != 1: raise Exception("The system is invalid.")
        else:
            variables = Matrix([f"x{i}" for i in range(order)])

    # Checks if it is an independent linear system
    if coeffs.det() != 0 and displayConsoleLog:
        print("___________________________________________________________________________\n\n")
        print("CHOLESKY METHOD ::: The linear system is an independent one (has only 1 solution).")
    elif coeffs.det() == 0:
        if displayConsoleLog:
            print("___________________________________________________________________________\n\n")
            print("CHOLESKY METHOD ::: The linear system has either no or infinite solutions.\n")
            print("___________________________________________________________________________\n\n")
        
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

        if not result: print("___________________________________________________________________________\n\n")

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
            print("___________________________________________________________________________\n\n")

        sol = solve(G.T*variables - auxSol, variables)
        sol = Matrix([sol.get(var) for var in variables])

        return G, sol

    else:
        return G, False



# ITERATIVE METHODS _______________________________________________________________________________________

# JACOBI-RICHARDSON'S (JR) METHOD
def jacobiRichardson(coeffs: Matrix, variables: Matrix, result: Matrix, var0: Matrix, epsilon = 10**-6, MAXITER = 100, displayConsoleLog = True, printTotalIterations = True):
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
    if coeffs.shape[0] == coeffs.shape[1] == variables.shape[0] == result.shape[0]  and variables.shape[1] == result.shape[1] == 1: 
        order = coeffs.shape[0] # System order (number of equations and variables)
    else:
        raise Exception("The system is invalid.")

    # Checks if it is an independent linear system
    if coeffs.det() != 0 and displayConsoleLog:
        print("___________________________________________________________________________\n\n")
        print("JACOBI-RICHARDSON METHOD ::: The linear system is an independent one (has only 1 solution).")
    elif coeffs.det() == 0:
        if displayConsoleLog:
            print("___________________________________________________________________________\n\n")
            print("JACOBI-RICHARDSON METHOD ::: The linear system has either no or infinite solutions.\n")
            print("___________________________________________________________________________\n\n")
        
        return

    if displayConsoleLog: printTotalIterations = True

    L = []
    D = []
    R = []

    for row in range(order):
        rowL = []
        rowD = []
        rowR = []

        for col in range(order):

            if row > col: 
                rowL.append(coeffs.row(row)[col])
                rowD.append(0)
                rowR.append(0)

            elif row == col: 
                rowL.append(0)
                rowD.append(coeffs.row(row)[col])
                rowR.append(0)

            elif row < col: 
                rowL.append(0)
                rowD.append(0)
                rowR.append(coeffs.row(row)[col])

        L.append(rowL)
        D.append(rowD)
        R.append(rowR)

    L = Matrix(L)
    D = Matrix(D)
    R = Matrix(R)

    if displayConsoleLog:
        print("JACOBI-RICHARDSON METHOD ::: Ax = b --> A = L + D + R.")
        print("JACOBI-RICHARDSON METHOD ::: L = ")
        pprint(L)
        print()

        print("JACOBI-RICHARDSON METHOD ::: D = ")
        pprint(D)
        print()

        print("JACOBI-RICHARDSON METHOD ::: R = ")
        pprint(R)
        print()

    if D.det() == 0:
        if displayConsoleLog:
            print("JACOBI-RICHARDSON METHOD ::: As |D| = 0, the method cannot be used.")

        return

    # A x = b --> x_k = B x_(k-1) + g
    B = - D**-1 * (L + R)
    g = D**-1 * result

    if displayConsoleLog:
        print("JACOBI-RICHARDSON METHOD ::: A x = b --> x_k = B x_(k-1) + g.")
        print("JACOBI-RICHARDSON METHOD ::: B = ")
        pprint(B)
        print()

        print("JACOBI-RICHARDSON METHOD ::: g = ")
        pprint(g)
        print()

    # Checks convergence
    if B.norm(oo) < 1:
        if displayConsoleLog:
            print(f"JACOBI-RICHARDSON METHOD ::: ||B||_oo = {B.norm(oo)} < 1, so there is guarantee of convergence.")
    elif B.norm(1) < 1:
        if displayConsoleLog:
            print(f"JACOBI-RICHARDSON METHOD ::: ||B||_1 = {B.norm(1)} < 1, so there is guarantee of convergence.")
    elif B.norm(2) < 1:
        if displayConsoleLog:
            print(f"JACOBI-RICHARDSON METHOD ::: ||B||_2 = {B.norm(2)} < 1, so there is guarantee of convergence.")
    else:
        if displayConsoleLog:
            print(f"JACOBI-RICHARDSON METHOD ::: ||B||_oo = {B.norm(oo)}, ||B||_1 = {B.norm(1)}, ||B||_2 = {B.norm(2)} are all > 1, so there is no guarantee of convergence.")
    if displayConsoleLog: print()

    # No. of decimal places for the results
    precision = int(abs(log10(epsilon))) if epsilon > 0 and log10(epsilon).is_integer else int(6)

    # Only 2 elements will be used, [0] for even iterations and [1] for odd ones
    x = [var0, 0]

    if displayConsoleLog:
        print("JACOBI-RICHARDSON METHOD ::: Iteration no.00: ", end="")
        for i in range(len(variables)): print(f"{variables[i]} = % .{precision}f" % x[0][i], end = ", ")
        print("||vars_k - vars_(k-1)||_oo =  " + "-"*(precision + 2))

    # Auxiliar index variable correspondent to the current element in the 'x' list
    j = 0

    # Loops through the algorithm until the maximum number of iterations is met
    for i in range(1, MAXITER):
        j += 1 if j < 1 else -j # j = 0 --> j = 1 --> j = 0 --> j = 1 ...

        x[j] = B*x[j-1] + g
        
        if displayConsoleLog:
            print(f"JACOBI-RICHARDSON METHOD ::: Iteration no.{str(i).zfill(2)}: ", end="")
            for k in range(len(variables)): print(f"{variables[k]} = % .{precision}f" % x[j][k], end = ", ")
            print(f"||vars_k - vars_(k-1)||_oo = % .{precision}f" % (x[j] - x[j-1]).norm(oo))


        if( (x[j] - x[j-1]).norm(oo) < epsilon * max(1, x[j].norm(oo)) ):
            if printTotalIterations:
                print("___________________________________________________________________________\n\n") if not displayConsoleLog else print("\n")
                print(f"JACOBI-RICHARDSON METHOD ::: The solution was successfully found after {i+1} iterations.\n")
                print("___________________________________________________________________________\n\n")

            return x[j], int(precision)

    # If it manages to get out of the loop, it means no root was found, so it prints an error message
    raise Exception("The maximum number of iterations was met and no root was found.")


# GAUSS-SEIDEL'S (GS) METHOD
def gaussSeidel(coeffs: Matrix, variables: Matrix, result: Matrix, var0: Matrix, epsilon = 10**-6, MAXITER = 100, displayConsoleLog = True, printTotalIterations = True):
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
    if coeffs.shape[0] == coeffs.shape[1] == variables.shape[0] == result.shape[0]  and variables.shape[1] == result.shape[1] == 1: 
        order = coeffs.shape[0] # System order (number of equations and variables)
    else:
        raise Exception("The system is invalid.")

    # Checks if it is an independent linear system
    if coeffs.det() != 0 and displayConsoleLog:
        print("___________________________________________________________________________\n\n")
        print("GAUSS-SEIDEL METHOD ::: The linear system is an independent one (has only 1 solution).")
    elif coeffs.det() == 0:
        if displayConsoleLog:
            print("___________________________________________________________________________\n\n")
            print("GAUSS-SEIDEL METHOD ::: The linear system has either no or infinite solutions.\n")
            print("___________________________________________________________________________\n\n")
        
        return

    if displayConsoleLog: printTotalIterations = True

    L = []
    D = []
    R = []

    for row in range(order):
        rowL = []
        rowD = []
        rowR = []

        for col in range(order):

            if row > col: 
                rowL.append(coeffs.row(row)[col])
                rowD.append(0)
                rowR.append(0)

            elif row == col: 
                rowL.append(0)
                rowD.append(coeffs.row(row)[col])
                rowR.append(0)

            elif row < col: 
                rowL.append(0)
                rowD.append(0)
                rowR.append(coeffs.row(row)[col])

        L.append(rowL)
        D.append(rowD)
        R.append(rowR)

    L = Matrix(L)
    D = Matrix(D)
    R = Matrix(R)

    if displayConsoleLog:
        print("GAUSS-SEIDEL METHOD ::: Ax = b --> A = L + D + R.")
        print("GAUSS-SEIDEL METHOD ::: L = ")
        pprint(L)
        print()

        print("GAUSS-SEIDEL METHOD ::: D = ")
        pprint(D)
        print()

        print("GAUSS-SEIDEL METHOD ::: R = ")
        pprint(R)
        print()

    if D.det() == 0:
        if displayConsoleLog:
            print("GAUSS-SEIDEL METHOD ::: As |D| = 0, the method cannot be used.")

        return

    # A x = b --> x_k = B x_(k-1) + g
    B = - (D**-1 * L + eye(order))**-1 * (D**-1 * R)
    g = (D**-1 * L + eye(order))**-1 * (D**-1 * result)

    if displayConsoleLog:
        print("GAUSS-SEIDEL METHOD ::: A x = b --> x_k = B x_(k-1) + g.")
        print("GAUSS-SEIDEL METHOD ::: B = ")
        pprint(B)
        print()

        print("GAUSS-SEIDEL METHOD ::: g = ")
        pprint(g)
        print()

    # Checks convergence
    if B.norm(oo) < 1:
        if displayConsoleLog:
            print(f"GAUSS-SEIDEL METHOD ::: ||B||_oo = {B.norm(oo)} < 1, so there is guarantee of convergence.")
    elif B.norm(1) < 1:
        if displayConsoleLog:
            print(f"GAUSS-SEIDEL METHOD ::: ||B||_1 = {B.norm(1)} < 1, so there is guarantee of convergence.")
    elif B.norm(2) < 1:
        if displayConsoleLog:
            print(f"GAUSS-SEIDEL METHOD ::: ||B||_2 = {B.norm(2)} < 1, so there is guarantee of convergence.")
    else:
        if displayConsoleLog:
            print(f"GAUSS-SEIDEL METHOD ::: ||B||_oo = {B.norm(oo)}, ||B||_1 = {B.norm(1)}, ||B||_2 = {B.norm(2)} are all > 1, so there is no guarantee of convergence.")
    if displayConsoleLog: print()

    # No. of decimal places for the results
    precision = int(abs(log10(epsilon))) if epsilon > 0 and log10(epsilon).is_integer else int(6)

    # Only 2 elements will be used, [0] for even iterations and [1] for odd ones
    x = [var0, 0]

    if displayConsoleLog:
        print("GAUSS-SEIDEL METHOD ::: Iteration no.00: ", end="")
        for i in range(len(variables)): print(f"{variables[i]} = % .{precision}f" % x[0][i], end = ", ")
        print("||vars_k - vars_(k-1)||_oo =  " + "-"*(precision + 2))

    # Auxiliar index variable correspondent to the current element in the 'x' list
    j = 0

    # Loops through the algorithm until the maximum number of iterations is met
    for i in range(1, MAXITER):
        j += 1 if j < 1 else -j # j = 0 --> j = 1 --> j = 0 --> j = 1 ...

        x[j] = B*x[j-1] + g
        
        if displayConsoleLog:
            print(f"GAUSS-SEIDEL METHOD ::: Iteration no.{str(i).zfill(2)}: ", end="")
            for k in range(len(variables)): print(f"{variables[k]} = % .{precision}f" % x[j][k], end = ", ")
            print(f"||vars_k - vars_(k-1)||_oo = % .{precision}f" % (x[j] - x[j-1]).norm(oo))


        if (x[j] - x[j-1]).norm(oo) < epsilon * max(1, x[j].norm(oo)) :
            if printTotalIterations:
                print("___________________________________________________________________________\n\n") if not displayConsoleLog else print("\n")
                print(f"GAUSS-SEIDEL METHOD ::: The solution was successfully found after {i+1} iterations.\n")
                print("___________________________________________________________________________\n\n")

            return x[j], int(precision)

    # If it manages to get out of the loop, it means no root was found, so it prints an error message
    raise Exception("The maximum number of iterations was met and no root was found.")



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


    # # A = Matrix([ [4,2,2], [2,2,1], [2,1,2] ])
    # A = Matrix([ [4, 12, -16], [12, 37, -43], [-16, -43, 98] ])
    # var = Matrix([ 'x', 'y', 'z' ])

    # G = cholesky(A, var)


    A = Matrix([ [10,2,1], [1,5,1], [2,3,10] ])
    var = Matrix([ 'x', 'y', 'z' ])
    res = Matrix([ 7, -8, 6 ])
    var0 = Matrix([7/10, -8/5, 8/5])

    sol, precision = jacobiRichardson(A, var, res, var0, epsilon=10**-5)

    for k in range(len(var)): 
        print(f"{var[k]} = % .{precision}f" % sol[k], end = "")
        if k < len(var) - 1: print(", ", end="")
        else: print("\n")

    sol, precision = gaussSeidel(A, var, res, var0, epsilon=10**-5)

    for k in range(len(var)): 
        print(f"{var[k]} = % .{precision}f" % sol[k], end = "")
        if k < len(var) - 1: print(", ", end="")
        else: print("\n")