import sys
from sympy import *
import nonlinear_equations as nonlinearEq
import nonlinear_systems as nonlinearSys
import linear_systems as linearSys

no = ["no", "não", "n", "false", "falso", "negative", "negativo", "denied", "nein",""]
yes = ["yes", "y", "sim", "ja", "ya", "si", "positivo", "true", "verdadeiro"]
exitNames = ["exit", "exit()", "close", "leave"]
goBackNames = ["back", "go back", "<-", "<", "back()", "up", "up()"]

nonlinearEqName = ["1", "equation", "single equation", "non-linear equation",
                   "nonlinear equations", "equação não-linear", "equação", "equação não linear"]
nonlinearSysName = ["2", "non-linear system", "nonlinear system", "sistema não-linear", "sistema de equações não-lineares",
                    "sistema de equações não lineares", "sistema não linear", "non-linear system of equations", "nonlinear system of equations"]
linearSysName = ["3", "linear system", "sistema linear", "sistema de equações lineares",
                 "sistema de equações lineares", "linear system of equations"]


while True:
    print("\nTHE AVALIABLE TYPES FOR NUMERICAL ANALYSIS ARE:")
    print("1. Non-Linear Equation")
    print("2. Non-Linear System of Equations")
    print("3. Linear System of Equations\n")

    slctdType = input("Select the desired type: ").lower()
    while True:
        if nonlinearEqName.count(slctdType) > 0:
            slctdType = "nonlinear equation"
            break
        elif nonlinearSysName.count(slctdType) > 0:
            slctdType = "nonlinear sys"
            break
        elif linearSysName.count(slctdType) > 0:
            slctdType = "linear sys"
            break
        elif exitNames.count(slctdType) > 0:
            print()
            sys.exit()
        else: slctdType = input("Invalid input, try again. Select the desired type: ").lower()

    print("\n___________________________________________________________________________\n")
    goBack = False

    if slctdType == "nonlinear equation":
        bisectionName = ['1', '1.', 'bisection method', 'bisection',
                         'bissecção', 'bisseção', 'método da bissecção', 'método da bisseção']
        linearIterationsName = ['2', '2.', 'linear iterations method', 'linear method', 'linear', 'iterations', 'iteration', 'iteração', 'iterações', 'iteração linear', 'iterações lineares', 'método das iterações lineares',
                                'método da iteração linear', 'ponto fixo', 'método do ponto fixo', 'método dos pontos fixos', 'pontos fixos', 'iteracao linear', 'iteracoes lineares', 'metodo da iteracao linear', 'metodo das iteracoes lineares']
        newtonsName = ['3', '3.', "newton's method", 'newton method',
                       'método de newton', 'newton', 'newtons method']
        secantName = ['4', '4.', 'secant method', "secant's method", 'secants method', 'secant', 'secants',
                      'secante', 'secantes', 'método da secante', 'método das secantes', 'método secante', 'método secantes']
        regulaFalsiName = ['5', '5.', 'regula falsi', 'regula', 'falsi', 'falsa', 'posição', 'falsa posição', 'posição falsa', 'método regula falsi', 'regula falsi method', 'method regula falsi',
                           'método da falsa posição', 'método da posição falsa', 'método posição falsa', 'método falsa posição', 'posição falsa método', 'regula falsi método', 'falsa posição método']

        while True:
            print("\nThe avaliable methods for non-linear equation's analysis are:")
            print("1. Bisection Method")
            print("2. Linear Iterations Method")
            print("3. Newton's Method")
            print("4. Secant Method")
            print("5. Regula Falsi Method\n")

            method = input("Select the desired method: ").lower()
            while True:
                if bisectionName.count(method) > 0:
                    method = "bisection"
                    break
                elif linearIterationsName.count(method) > 0:
                    method = "linear iterations"
                    break
                elif newtonsName.count(method) > 0:
                    method = "newton's"
                    break
                elif secantName.count(method) > 0:
                    method = "secant"
                    break
                elif regulaFalsiName.count(method) > 0:
                    method = "regula falsi"
                    break
                elif exitNames.count(method) > 0:
                    print()
                    sys.exit()
                elif goBackNames.count(method) > 0:
                    goBack = True
                    print("\n___________________________________________________________________________\n")
                    break
                else: method = input("Invalid input, try again. Select the desired method: ").lower()

            if goBack: break

            print("\n::: " + method.title() + " method selected.")

            if method != "linear iterations":
                fun = input("\nEnter function f(x) = ")
                while True:
                    try: fun = sympify(fun)
                    except: fun = input("Invalid input, try again. Enter function f(x) = ")
                    else: 
                        print("\n::: f(x) = ")
                        pprint(fun)
                        break
            else:
                fun = input("\nEnter function psi(x) = ")
                while True:
                    try: fun = sympify(fun)
                    except: fun = input("Invalid input, try again. Enter function psi(x) = ")
                    else: 
                        print("\n::: psi(x) = ")
                        pprint(fun)
                        break

                confirmConv = input("\nConfirm that the result converge for this psi(x)? ")
                while True:
                    if yes.count(confirmConv):
                        confirmConv = True
                        break
                    elif no.count(confirmConv):
                        confirmConv = False
                        break
                    else: confirmConv = input("Invalid input, try again. Confirm that the result converge for this psi(x)? ")

            if method == "bisection" or (method == "linear iterations" and confirmConv):
                print("\nConsidering the interval [a,b], enter:")
                a = input("a = ")
                while True:
                    try: a = float(sympify(a))
                    except: a = input("Invalid input, try again. a = ")
                    else: break
                b = input("b = ")
                while True:
                    try: b = float(sympify(b))
                    except: b = input("Invalid input, try again. b = ")
                    else: 
                        print("\n::: Interval [%.3f, %.3f]." % (a, b))
                        break

            if method != "bisection":
                x_0 = input("\nEnter x_0 = ")
                while True:
                    try: x_0 = float(sympify(x_0))
                    except: x_0 = input("Invalid input, try again. x_0 = ")
                    else: break

                if method == "secant" or method == "regula falsi":
                    x_1 = input("\nEnter x_1 = ")
                    while True:
                        try: x_1 = float(sympify(x_1))
                        except: x_1 = input("Invalid input, try again. x_1 = ")
                        else: 
                            print("\n::: x_0 = %f ; x_1 = %f." % (x_0, x_1))
                            break

                else: print("\n::: x_0 = %f." % x_0)

            error = input("\nEnter the desired error margin: ")
            while True:
                if no.count(error) > 0: 
                    error = False
                    break
                else:
                    try: error = float(sympify(error))
                    except: error = input("Invalid input, try again. Enter the desired error margin: ")
                    else:
                        if error >= 1 or error == 0: error = input("Invalid input, the error margin must be less than 1. Try again: ")
                        else: break

            maxIterations = input("\nEnter the maximum number of iterations desired: ")
            while True:
                if no.count(maxIterations) > 0: 
                    maxIterations = False
                    break
                try: maxIterations = int(sympify(maxIterations))
                except: maxIterations = input("Invalid input, try again. Enter the maximum number of iterations: ")
                else:
                    if maxIterations < 1: maxIterations = input("Invalid input, there must be at least 1 iterations. Try again: ")
                    else: break

            showLog = input("\nDisplay console log for all values? ")
            while True:
                if yes.count(showLog):
                    showLog = True
                    showTotalNo = True
                    break
                elif no.count(showLog):
                    showLog = False

                    showTotalNo = input("\nDisplay total number of iterations? ")
                    while True:
                        if yes.count(showTotalNo):
                            showTotalNo = True
                            break
                        elif no.count(showTotalNo):
                            showTotalNo = False
                            break
                        else: showTotalNo = input("Invalid input, try again. Display total number of iterations? ")

                    break
                else: showLog = input("Invalid input, try again. Display console log for all values? ")

            print("\n")

            if method == "bisection":
                if error and maxIterations: root, precision = nonlinearEq.bisection(fun,a,b,error,maxIterations,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                elif error: root, precision = nonlinearEq.bisection(fun,a,b,error,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                elif maxIterations: root, precision = nonlinearEq.bisection(fun,a,b,MAXITER=maxIterations,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                else: root, precision = nonlinearEq.bisection(fun,a,b,displayConsoleLog=showLog,printTotalIterations=showTotalNo)

            elif method == "linear iterations":
                if confirmConv:
                    if error and maxIterations: root, precision = nonlinearEq.linearIterations(fun,x_0,error,maxIterations,confirmConv,a,b,showLog,showTotalNo)
                    elif error: root, precision = nonlinearEq.linearIterations(fun,x_0,error,confirmConvergence=confirmConv,infLim=a,supLim=b,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                    elif maxIterations: root, precision = nonlinearEq.linearIterations(fun,x_0,MAXITER=maxIterations,confirmConvergence=confirmConv,infLim=a,supLim=b,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                    else: root, precision = nonlinearEq.linearIterations(fun,x_0,confirmConvergence=confirmConv,infLim=a,supLim=b,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                else:
                    if error and maxIterations: root, precision = nonlinearEq.linearIterations(fun,x_0,error,maxIterations,confirmConv,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                    elif error: root, precision = nonlinearEq.linearIterations(fun,x_0,error,confirmConvergence=confirmConv,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                    elif maxIterations: root, precision = nonlinearEq.linearIterations(fun,x_0,MAXITER=maxIterations,confirmConvergence=confirmConv,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                    else: root, precision = nonlinearEq.linearIterations(fun,x_0,confirmConvergence=confirmConv,displayConsoleLog=showLog,printTotalIterations=showTotalNo)

            elif method == "newton's":
                if error and maxIterations: root, precision = nonlinearEq.newtons(fun,x_0,error,maxIterations,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                elif error: root, precision = nonlinearEq.newtons(fun,x_0,error,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                elif maxIterations: root, precision = nonlinearEq.newtons(fun,x_0,MAXITER=maxIterations,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                else: root, precision = nonlinearEq.newtons(fun,x_0,displayConsoleLog=showLog,printTotalIterations=showTotalNo)

            elif method == "secant":
                if error and maxIterations: root, precision = nonlinearEq.secant(fun,x_0,x_1,error,MAXITER=maxIterations,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                elif error: root, precision = nonlinearEq.secant(fun,x_0,x_1,error,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                elif maxIterations: root, precision = nonlinearEq.secant(fun,x_0,x_1,MAXITER=maxIterations,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                else: root, precision = nonlinearEq.secant(fun,x_0,x_1,displayConsoleLog=showLog,printTotalIterations=showTotalNo)

            elif method == "regula falsi":
                if error and maxIterations: root, precision = nonlinearEq.regulaFalsi(fun,x_0,x_1,error,MAXITER=maxIterations,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                elif error: root, precision = nonlinearEq.regulaFalsi(fun,x_0,x_1,error,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                elif maxIterations: root, precision = nonlinearEq.regulaFalsi(fun,x_0,x_1,MAXITER=maxIterations,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                else: root, precision = nonlinearEq.regulaFalsi(fun,x_0,x_1,displayConsoleLog=showLog,printTotalIterations=showTotalNo)

            print("ROOT:")
            print(f"x = {root:.{precision}f}")
            print("\n___________________________________________________________________________\n")

    elif slctdType == "nonlinear sys":
        linearIterationsName = ['1', '1.', 'linear iterations method', 'linear method', 'linear', 'iterations', 'iteration', 'iteração', 'iterações', 'iteração linear', 'iterações lineares', 'método das iterações lineares',
                                'método da iteração linear', 'ponto fixo', 'método do ponto fixo', 'método dos pontos fixos', 'pontos fixos', 'iteracao linear', 'iteracoes lineares', 'metodo da iteracao linear', 'metodo das iteracoes lineares']
        newtonsName = ['2', '2.', "newton's method", 'newton method',
                       'método de newton', 'newton', 'newtons method']

        while True:
            print("\nThe avaliable methods for non-linear equations system's analysis are:")
            print("1. Linear Iterations Method")
            print("2. Newton's Method")

            method = input("Select the desired method: ").lower()
            while True:
                if linearIterationsName.count(method) > 0:
                    method = "linear iterations"
                    break
                elif newtonsName.count(method) > 0:
                    method = "newton's"
                    break
                elif exitNames.count(method) > 0:
                    print()
                    sys.exit()
                elif goBackNames.count(method) > 0:
                    goBack = True
                    print("\n___________________________________________________________________________\n")
                    break
                else:
                    method = input("Invalid input, try again. Select the desired method: ").lower()

            if goBack: break

            print("\n::: " + method.title() + " method selected.")

            order = input("\nWhat is the system's order (no. of equations, variables)? ")
            while True:
                try: order = int(order)
                except: order = input("Invalid input, try again. What is the system's order? ")
                else:
                    print(f"\n::: The system has order {order}. This must be the number of functions (equations), as well as variables.")
                    break

            while True:
                variables = input("\nWhat variables will be used? Separate them with a whitespace: ").split()
                hasInvalid = False
                for var in variables:
                    if var.isnumeric(): hasInvalid = True; break
                if hasInvalid: print("Invalid input, try again. A variable's 1st character can't be a number.")
                else: 
                    print("\n::: The variables used are: " + ", ".join(variables) + ".")
                    break

            fun = []
            if method == "newton's":
                print()
                for k in range(order):
                    fun.append(input(f"Enter function f_{k} = "))
                    while True:
                        try:
                            fun[k] = sympify(fun[k])
                        except:
                            fun[k] = input(f"Invalid input, try again. Enter function f_{k} = ")
                        else:
                            break
                for k in range(order):
                    print(f"\n::: f_{k}(" + ", ".join(variables) + ") = ")
                    pprint(fun[k])

            elif method == "linear iterations":
                for k in range(order):
                    print()
                    fun.append(input(f"Enter function F_{k} = "))
                    while True:
                        try:
                            fun[k] = sympify(fun[k])
                        except:
                            fun[k] = input(f"Invalid input, try again. Enter function f_{k} = ")
                        else:
                            break
                for k in range(order):
                    print(f"\n::: F_{k}(" + ", ".join(variables) + ") = ")
                    pprint(fun[k])

            var_0 = []
            print()
            for k in range(order):
                var_0.append(input(f"Enter {variables[k]}_0 = "))
                while True:
                    try:
                        var_0[k] = float(sympify(var_0[k]))
                    except:
                        var_0[k] = (input(f"Invalid input, try again. {variables[k]}_0 = "))
                    else: break
            print()
            for k in range(order): print(f"::: {variables[k]}_0 = {var_0[k]}")

            error = input("\nEnter the desired error margin: ")
            while True:
                if no.count(error) > 0:
                    error = False
                    break
                else:
                    try:
                        error = float(sympify(error))
                    except:
                        error = input("Invalid input, try again. Enter the desired error margin: ")
                    else:
                        if error >= 1 or error == 0:
                            error = input("Invalid input, the error margin must be less than 1. Try again: ")
                        else:
                            break

            maxIterations = input("\nEnter the maximum number of iterations desired: ")
            while True:
                if no.count(maxIterations) > 0:
                    maxIterations = False
                    break
                try:
                    maxIterations = int(sympify(maxIterations))
                except:
                    maxIterations = input("Invalid input, try again. Enter the maximum number of iterations: ")
                else:
                    if maxIterations < 1: maxIterations = input("Invalid input, there must be at least 1 iterations. Try again: ")
                    else: break

            showLog = input("\nDisplay console log for all values? ")
            while True:
                if yes.count(showLog):
                    showLog = True
                    showTotalNo = True
                    break

                elif no.count(showLog):
                    showLog = False
                    showTotalNo = input("\nDisplay total number of iterations? ")
                    while True:
                        if yes.count(showTotalNo):
                            showTotalNo = True
                            break
                        elif no.count(showTotalNo):
                            showTotalNo = False
                            break
                        else:
                            showTotalNo = input("Invalid input, try again. Display total number of iterations? ")
                    break

                else:
                    showLog = input("Invalid input, try again. Display console log for all values? ")

            print("\n")

            if method == "linear iterations":
                if error and maxIterations:
                    sol, precision = nonlinearSys.linearIterations(
                        fun, variables, var_0, error, maxIterations, confirmConv, displayConsoleLog=showLog, printTotalIterations=showTotalNo)
                elif error:
                    sol, precision = nonlinearSys.linearIterations(
                        fun, variables, var_0, error, confirmConvergence=confirmConv, displayConsoleLog=showLog, printTotalIterations=showTotalNo)
                elif maxIterations:
                    sol, precision = nonlinearSys.linearIterations(
                        fun, variables, var_0, MAXITER=maxIterations, confirmConvergence=confirmConv, displayConsoleLog=showLog, printTotalIterations=showTotalNo)
                else:
                    sol, precision = nonlinearSys.linearIterations(
                        fun, variables, var_0, confirmConvergence=confirmConv, displayConsoleLog=showLog, printTotalIterations=showTotalNo)

            elif method == "newton's":
                if error and maxIterations:
                    sol, precision = nonlinearSys.newtons(
                        fun, variables, var_0, error, maxIterations, displayConsoleLog=showLog, printTotalIterations=showTotalNo)
                elif error:
                    sol, precision = nonlinearSys.newtons(
                        fun, variables, var_0, error, displayConsoleLog=showLog, printTotalIterations=showTotalNo)
                elif maxIterations:
                    sol, precision = nonlinearSys.newtons(
                        fun, variables, var_0, MAXITER=maxIterations, displayConsoleLog=showLog, printTotalIterations=showTotalNo)
                else:
                    sol, precision = nonlinearSys.newtons(
                        fun, variables, var_0, displayConsoleLog=showLog, printTotalIterations=showTotalNo)

            print("SOLUTION:")
            for k in range(order): print(f"{variables[k]} = {sol[k]:.{precision}f}")
            print("\n___________________________________________________________________________\n")

    elif slctdType == "linear sys":
        gaussianEliminationName = ["1", "1.", "gaussian elimination", "gaussian", "elimination",
                                   "exact gaussian", "exact gauss", "gaussian method", "gaussian elimination method", "elimination method"]
        choleskyName = ["2", "2.", "cholesky",
                        "cholesky method", "exact cholesky"]
        jrName = ["3", "3.", "jr", "jacobi-richardson", "jacobirichardson", "jacobi", "richardson",
                  "jacobi method", "jacobi-richardson method", "jacobirichardson method", "richardson method", "jr method"]
        gsName = ["4", "4.", "gauss", "seidel", "gauss method", "gauss-seidel",
                  "gs", "gs method", "seidel method", "gauss-seidel method", "gaussseidel"]

        while True:
            print("\nThe avaliable methods for linear system equation's analysis are:")
            print("EXACT METHODS:")
            print(" 1. Gaussian Elimination Method")
            print(" 2. Cholesky Method")
            print("ITERATIVE METHODS:")
            print(" 3. Jacobi-Richardson Method")
            print(" 4. Gauss-Seidel Method\n")

            method = input("Select the desired method: ").lower()
            while True:
                if gaussianEliminationName.count(method) > 0:
                    method = "gaussian elimination"
                    break
                elif choleskyName.count(method) > 0:
                    method = "cholesky"
                    break
                elif jrName.count(method) > 0:
                    method = "jacobi-richardson"
                    break
                elif gsName.count(method) > 0:
                    method = "gauss-seidel"
                    break
                elif exitNames.count(method) > 0:
                    print()
                    sys.exit()
                elif goBackNames.count(method) > 0:
                    goBack = True
                    print("\n___________________________________________________________________________\n")
                    break
                else: method = input("Invalid input, try again. Select the desired method: ").lower()

            if goBack: break

            print("\n::: " + method.title() + " method selected.")

            print("\nConsider the system A x = b.\n")

            order = input("What is the system's order? (no. of equations, variables) ")
            while True:
                try: order = int(order)
                except: order = input("Invalid input, try again. What is the system's order? ")
                else: 
                    print(f"\n::: The system has order {order}. This must be the number rows and columns in the A matrix, as well as variables and results.")
                    break

            if method == "cholesky":
                print("\nThe Cholesky Method can calculate both the system's solution or only the G matrix (as in A = G G^T).")
                results = input("Do you wish to calculate the system's solution? ")
                while True:
                    if no.count(results) > 0:
                        results = None
                        variables = None
                        break
                    
                    elif yes.count(results) > 0:
                        results = True
                        break

            if method != "cholesky" or (method == "cholesky" and variables):
                variables = input("\nWhat variables will be used? Separate them with a whitespace: ").split()
                while True:
                    hasInvalid = False
                    for var in variables:
                        if var.isnumeric(): hasInvalid = True; break
                    if hasInvalid: variables = input("Invalid input, a variable's 1st character can't be a number. Try again: ").split()
                    else: 
                        print("\n::: The variables used are: " + ", ".join(variables) + ".")
                        variables = Matrix(variables)
                        break

            print("\nEnter A matrix (coefficients) - each row separated by a linebreak and each column separated by one whitespace.")
            while True:
                coeffs = []
                for i in range(order): coeffs.append(input().split())

                try: coeffs = Matrix(coeffs)

                except: 
                    print("Invalid input, try again. Each row must be separated by exactly one linebreak and each column must be separated by exactly one whitespace.")
                    print("e.g. order 3:")
                    print("1 2 3")
                    print("4 5 6")
                    print("7 8 9")
                    print()
                    
                else: 
                    hasInvalid = False
                    for coeff in coeffs:
                        try: float(sympify(coeff))

                        except: hasInvalid = True; break

                    if hasInvalid: 
                        print("Invalid input, the A matrix must be numeric. Try again: ")

                    elif coeffs.shape[0] != coeffs.shape[1] or coeffs.shape[0] != order: 
                        print(f"Invalid input, the A matrix must be a square matrix of order {order}. Try again: ")

                    else:
                        print("\n::: A = ")
                        pprint(coeffs)
                        break


            if method != "cholesky" or (method == "cholesky" and results):
                results = input("\nEnter the transpose of the b vector (results) - oneline with elements separated by one whitespace. ").split()
                while True:
                    hasInvalid = False

                    for res in results:
                        try: res = float(sympify(res))

                        except: 
                            results = input("Invalid input, the results must be numeric. Try again: ").split()
                            hasInvalid = True
                            break

                    if not hasInvalid:
                        results = Matrix(results)
                        print("\n::: b = ") 
                        pprint(results)
                        break

            if method == "jacobi-richardson" or method == "gauss-seidel":
                var_0 = input("\nEnter the transpose of the x_0 vector (variables approx. for iteration 0) - oneline with elements separated by one whitespace. ").split()
                while True:
                    hasInvalid = False

                    for var in var_0:
                        try: var = float(sympify(var))

                        except: 
                            var_0 = input("Invalid input, the approximations must be numeric. Try again: ").split()
                            hasInvalid = True
                            break

                    if not hasInvalid:
                        var_0 = Matrix(var_0)
                        print("\n::: x_0 = ") 
                        pprint(var_0)
                        break

                error = input("\nEnter the desired error margin: ")
                while True:
                    if no.count(error) > 0: 
                        error = False
                        break
                    else:
                        try: error = float(sympify(error))
                        except: error = input("Invalid input, try again. Enter the desired error margin: ")
                        else:
                            if error >= 1 or error == 0: error = input("Invalid input, the error margin must be less than 1. Try again: ")
                            else: break

                maxIterations = input("\nEnter the maximum number of iterations desired: ")
                while True:
                    if no.count(maxIterations) > 0: 
                        maxIterations = False
                        break
                    try: maxIterations = int(sympify(maxIterations))
                    except: maxIterations = input("Invalid input, try again. Enter the maximum number of iterations: ")
                    else:
                        if maxIterations < 1: maxIterations = input("Invalid input, there must be at least 1 iterations. Try again: ")
                        else: break

            showLog = input("\nDisplay console log for all values? ")
            while True:
                if yes.count(showLog):
                    showLog = True
                    showTotalNo = True
                    break

                elif no.count(showLog):
                    showLog = False

                    if method == "jacobi-richardson" or method == "gauss-seidel":
                        showTotalNo = input("\nDisplay total number of iterations? ")
                        while True:
                            if yes.count(showTotalNo):
                                showTotalNo = True
                                break
                            elif no.count(showTotalNo):
                                showTotalNo = False
                                break
                            else: showTotalNo = input("Invalid input, try again. Display total number of iterations? ")

                    break

                else: showLog = input("Invalid input, try again. Display console log for all values? ")

            print("\n")

            if method == "gaussian elimination":
                coeffs, result, sol = linearSys.gaussianElimination(coeffs, variables, results, displayConsoleLog=showLog)

                print("SOLUTION FOR 'A x = b': \n")

                print("A =")
                pprint(coeffs)
                print()

                print("b =")
                pprint(result)
                print()

                print("x =")
                pprint(sol)
                print()

            elif method == "cholesky":
                G, sol = linearSys.cholesky(coeffs, variables, result=results, displayConsoleLog=showLog)

                print("SOLUTION FOR 'A x = G G^T x = b': \n") if sol else print("DECOMPOSITION A = G G^T:")

                print("G =")
                pprint(G)
                print()

                print("G^T =")
                pprint(G.T)
                print()

                if sol:
                    print("x =")
                    pprint(sol)
                    print()

            elif method == "jacobi-richardson":
                if error and maxIterations: sol, precision = linearSys.jacobiRichardson(coeffs,variables,results,var_0,epsilon=error,MAXITER=maxIterations,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                elif error: sol, precision = linearSys.jacobiRichardson(coeffs,variables,results,var_0,error,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                elif maxIterations: sol, precision = linearSys.jacobiRichardson(coeffs,variables,results,var_0,MAXITER=maxIterations,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                else: sol, precision = linearSys.jacobiRichardson(coeffs,variables,results,var_0,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                

            elif method == "gauss-seidel":
                if error and maxIterations: sol, precision = linearSys.gaussSeidel(coeffs,variables,results,var_0,error,MAXITER=maxIterations,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                elif error: sol, precision = linearSys.gaussSeidel(coeffs,variables,results,var_0,error,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                elif maxIterations: sol, precision = linearSys.gaussSeidel(coeffs,variables,results,var_0,MAXITER=maxIterations,displayConsoleLog=showLog,printTotalIterations=showTotalNo)
                else: sol, precision = linearSys.gaussSeidel(coeffs,variables,results,var_0,displayConsoleLog=showLog,printTotalIterations=showTotalNo)

            if method == "jacobi-richardson" or method == "gauss-seidel":
                print("SOLUTION:")
                for k in range(order): print(f"{variables[k]} = % .{precision}f" % sol[k])

            print("\n___________________________________________________________________________\n")
