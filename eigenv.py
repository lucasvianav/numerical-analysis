from sympy.abc import _clash1 # So that sympify will be able to convert string to symbolic objects 
from statistics import mean # Average method 
from math import log10 # Log base 10 
from sympy import *
from sympy.solvers.diophantine.diophantine import norm # Module to work with symbolic variables 
from util import * # Auxiliar function in  util.py

# POWER METHOD
def powerMethod(A: Matrix, y_0: Matrix, epsilon = 10**-6, MAXITER = 100, assureConvergence = True, displayConsoleLog = True, printTotalIterations = True):
    # Checks if all A's elements are numeric
    try: 
        for e in A: e = float(e)
    except: 
        raise Exception("All of A's must be numeric.")

    # Checks if A is a square matrix
    if A.shape[0] != A.shape[1]:
        raise Exception("A must be a square matrix.")

    # Checks if y_0 is null
    if y_0 == 0:
        raise Exception("The initial approximation y_0 must be != 0.")

    # Checks if y_0 is a column vector with as much rows as A
    if y_0.shape[0] != A.shape[0] or y_0.shape[1] != 1:
        raise Exception("The initial approximation y_0 must be a column vector with as many rows as the matrix A.")

    # Checks if ||y_0||oo == 1
    if assureConvergence:
        if y_0.norm(oo) != 1:
            raise Exception("The initial approximation y_0 must have an infinity norm of 1.")

    # No. of decimal places for the results
    precision = int(abs(log10(epsilon))) if epsilon > 0 and log10(epsilon).is_integer else int(6)

    if displayConsoleLog: 
        printTotalIterations = True
        print("___________________________________________________________________________\n\n")
        print(f"POWER METHOD ::: Iteration no. 00: z^T = %s ; alpha = %s ; y^T = [ %s  ]"
            % ( '-'*((precision + 5)*len(y_0) + 3), '-'*(precision + 2), "  ".join([str(f"% .{precision}f" % e) for e in y_0]) )
        )

    z = [0, 0]
    y = [y_0, 0] # --> u (eigenvector)
    a = [0, 0] # alpha --> lambda (eigenvalue)

    # Auxiliar index variable correspondent to the current element in the above lists
    j = 0

    # When the solution was found, the algorithm'll go over one more
    # iteration in order to check the eigenvalue's signal
    # And this variable checks is that is the extra iterations or not
    hasFinished = False

    for i in range(1, MAXITER+1):
        j += 1 if j < 1 else -j # j = 0 --> j = 1 --> j = 0 --> j = 1 ...

        z[j] = (A*y[j-1]).evalf()
        a[j] = (z[j].norm(oo)).evalf()
        y[j] = (z[j]/a[j]).evalf()

        if displayConsoleLog and not hasFinished:
            print(f"POWER METHOD ::: Iteration no. %s: z^T = [ %s  ] ; alpha = %.{precision}f ; y^T = [ %s  ]"
                % ( str(i).zfill(2), "  ".join([str(f"% .{precision}f" % e) for e in z[j]]), a[j], "  ".join([str(f"% .{precision}f" % e) for e in y[j]]) )
            )

        if abs(a[j] - a[j-1]) < epsilon * max(1, abs(a[j])):
            # If the solution was just found, go over one more iteration
            if not hasFinished: hasFinished = True

            else:
                # Trucante the eigenvector's values to the requested precision
                y[j] = Matrix([float(f"% .{precision}f" % e) for e in y[j]])
                y[j-1] = Matrix([float(f"% .{precision}f" % e) for e in y[j-1]])

                # Checks if the lambda = |lambda| or lambda = -|lambda|
                if y[j] != y[j-1] and abs(y[j]) == abs(y[j-1]): a[j] = -a[j]

                if printTotalIterations:
                    print("___________________________________________________________________________\n\n") if not displayConsoleLog else print("\n")
                    print(f"POWER METHOD ::: The solution was successfully found after {i+1} iterations.\n")
                    print("___________________________________________________________________________\n\n")

                return a[j], y[j], precision # eigenvalue, eigenvector

    # If it manages to get out of the loop, it means no root was found, so it prints an error message
    raise Exception("The maximum number of iterations was met and no root was found.")

def inversePower(A: Matrix, y_0: Matrix, epsilon = 10**-6, MAXITER = 100, assureConvergence = True, displayConsoleLog = True, printTotalIterations = True):
    # Checks if all A's elements are numeric
    try: 
        for e in A: e = float(e)
    except: 
        raise Exception("All of A's elements must be numeric.")

    # Checks if A is a square matrix
    if A.shape[0] != A.shape[1]:
        raise Exception("A must be a square matrix.")

    # Checks if y_0 is null
    if y_0 == 0:
        raise Exception("The initial approximation y_0 must be != 0.")

    # Checks if y_0 is a column vector with as much rows as A
    if y_0.shape[0] != A.shape[0] or y_0.shape[1] != 1:
        raise Exception("The initial approximation y_0 must be a column vector with as many rows as the matrix A.")

    # Checks if ||y_0||oo == 1
    if assureConvergence:
        if y_0.norm(oo) != 1:
            raise Exception("The initial approximation y_0 must have an infinity norm of 1.")

    # Inverts A
    try:
        A = A**-1
    except:
        raise Exception("As A is not invertible, the method cannot be applied.")


    # No. of decimal places for the results
    precision = int(abs(log10(epsilon))) if epsilon > 0 and log10(epsilon).is_integer else int(6)

    if displayConsoleLog: 
        printTotalIterations = True
        print("___________________________________________________________________________\n\n")
        print(f"INVERSE POWER METHOD ::: Iteration no. 00: z^T = %s ; alpha = %s ; y^T = [ %s  ]"
            % ( '-'*((precision + 5)*len(y_0) + 3), '-'*(precision + 2), "  ".join([str(f"% .{precision}f" % e) for e in y_0]) )
        )


    z = [0, 0]
    y = [y_0, 0] # --> u (eigenvector)
    a = [0, 0] # alpha --> lambda (eigenvalue)

    # Auxiliar index variable correspondent to the current element in the above lists
    j = 0

    # When the solution was found, the algorithm'll go over one more
    # iteration in order to check the eigenvalue's signal
    # And this variable checks is that is the extra iterations or not
    hasFinished = False

    for i in range(1, MAXITER+1):
        j += 1 if j < 1 else -j # j = 0 --> j = 1 --> j = 0 --> j = 1 ...

        z[j] = (A*y[j-1]).evalf()
        a[j] = (z[j].norm(oo)).evalf()
        y[j] = (z[j]/a[j]).evalf()

        if displayConsoleLog and not hasFinished:
            print(f"INVERSE POWER METHOD ::: Iteration no. %s: z^T = [ %s  ] ; alpha = %.{precision}f ; y^T = [ %s  ]"
                % ( str(i).zfill(2), "  ".join([str(f"% .{precision}f" % e) for e in z[j]]), a[j], "  ".join([str(f"% .{precision}f" % e) for e in y[j]]) )
            )

        if abs(a[j] - a[j-1]) < epsilon * max(1, abs(a[j])):
            # If the solution was just found, go over one more iteration
            if not hasFinished: hasFinished = True

            else:
                # Trucante the eigenvector's values to the requested precision
                y[j] = Matrix([float(f"% .{precision}f" % e) for e in y[j]])
                y[j-1] = Matrix([float(f"% .{precision}f" % e) for e in y[j-1]])

                # Checks if the lambda = |lambda| or lambda = -|lambda|
                if y[j] != y[j-1] and abs(y[j]) == abs(y[j-1]): a[j] = -a[j]

                if printTotalIterations:
                    print("___________________________________________________________________________\n\n") if not displayConsoleLog else print("\n")
                    print(f"INVERSE POWER METHOD ::: The solution was successfully found after {i+1} iterations.\n")
                    print("___________________________________________________________________________\n\n")

                return a[j], y[j], precision # eigenvalue, eigenvector

    # If it manages to get out of the loop, it means no root was found, so it prints an error message
    raise Exception("The maximum number of iterations was met and no root was found.")

# is not working properly
def jacobi(A: Matrix, displayConsoleLog = True):
    # Checks if all A's elements are numeric
    try: 
        for e in A: e = float(e)
    except: 
        raise Exception("All of A's elements must be numeric.")

    # Checks if A is a square matrix
    if A.shape[0] != A.shape[1] or A.shape[0] <= 1: 
        raise Exception("A must be a square matrix.")
    else:
        order = A.shape[0]

    # Checks if the matrix is simetric
    if not A.is_symmetric(): raise Exception("As A is not symmetric, the method can't be applied.")

    D = [A, 0] # diagonal eigenvalues' matrix
    V = 1 # eigenvectors matrix

    # Auxiliar index variable correspondent to the current element in the 'D' list
    j = 0

    # Index variable correspondent to the current step
    i = 1

    # Total number of steps
    noSteps = len(str(int((len(A)-order)/2)))
    
    if displayConsoleLog:
        print("___________________________________________________________________________\n\n")
        print(f"JACOBI'S METHOD ::: Step no.{str(i).zfill(noSteps)}: A = ")
        pprint(A)
        print()

    while True:
        i += 1
        j += 1 if j < 1 else -j # j = 0 --> j = 1 --> j = 0 --> j = 1 ...

        a_pq = 0; p = q = 0
        for x in range(order):
            for y in range(x+1, order):
                if x == y: continue

                if abs(float(D[j-1].row(x)[y])) >= abs(a_pq):
                    a_pq = float(D[j-1].row(x)[y])
                    p = x
                    q = y

        if displayConsoleLog:
            print(f"JACOBI'S METHOD ::: Step no.{str(i).zfill(noSteps)}:  a_pq = % .3f ; p = {p+1} ; q = {q+1}"
                % a_pq
            )

        PHI = (D[j-1].row(q)[q] - D[j-1].row(p)[p])/(2*a_pq)
        t = 1/(PHI+(PHI/abs(PHI))*sqrt(1+PHI**2)) if PHI != 0 else 1
        cos_phi = 1/sqrt(1+t**2)
        sen_phi = t*cos_phi

        if displayConsoleLog:
            print(f"JACOBI'S METHOD ::: Step no.{str(i).zfill(noSteps)}: Φ = {PHI} = % .3f ; t = {t} = %.3f ; cos(φ) = {cos_phi} = %.3f ; sen(φ) = {sen_phi} = %.3f"
                % ( PHI, t, cos_phi, sen_phi )
            )

        U = []
        auxU = []
        for x in range(order):
            row = []
            auxRow = []

            for y in range(order):
                if x == y == p or x == y == q: 
                    row.append(cos_phi)
                    auxRow.append("cos(phi)")

                elif x == p and y == q:
                    row.append(sen_phi)
                    auxRow.append("sin(phi)")

                elif x == q and y == p:
                    row.append(-1*sen_phi)
                    auxRow.append("-sin(phi)")

                elif x == y:
                    row.append(1)
                    auxRow.append(1)

                else:
                    row.append(0)
                    auxRow.append(0)

            U.append(row)
            auxU.append(auxRow)

        U = Matrix(U)
        auxU = Matrix(auxU)

        if displayConsoleLog:
            print(f"JACOBI'S METHOD ::: Step no.{str(i).zfill(noSteps)}: U = ")
            pprint(auxU)
            print()
            pprint(U)
            print()

        D[j] = (U.T * D[j-1]) * U
        V = V*U

        if displayConsoleLog:
            print(f"JACOBI'S METHOD ::: Step no.{str(i).zfill(noSteps)}: A = ")
            pprint(D[j])
            print()

        # When D is diagonal, the proccess is finished
        if D[j].is_diagonal(): break
            
    return D[j], V


if __name__ == "__main__":
    A = Matrix([[-2, 6], [-6, 13]])
    y_0 = Matrix([1, 1])

    value, vector, precision = powerMethod(A, y_0, epsilon=10**-3)

    print(f"%.{precision}f" % value)
    pprint(vector)

    # A = Matrix([[2, 4, -2], [4, 2, 2], [-2, 2, 5]])

    # D, V = jacobi(A, displayConsoleLog=True)

    # pprint(D.evalf())
    # print()

    # pprint(V.evalf())
    # print()