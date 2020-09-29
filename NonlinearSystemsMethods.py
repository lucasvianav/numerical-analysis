from sympy.abc import _clash1 # So that sympify will be able to convert string to symbolic objects 
from statistics import mean # Average method 
from math import log10 # Log base 10 
from sympy import * # Module to work with symbolic variables 
from util import * # Auxiliar function in  util.py 

# LINEAR ITERATIONS METHOD
# Receives F(x,y) (= x) and G(x,y) (= y) and searches for f(x) and g(x)'s root starting at a point (x_0,y_0), with an error margin of epsilon. Goes as far as MAXITER iterations
def linearIterations(F_xy, G_xy, x_0: float, y_0: float, epsilon = 10**-6, MAXITER = 100, displayConsoleLog = False, printTotalIterations = False):
    try: # Try and convert the functions from string to symbolic
        F_xy = sympify(F_xy)
        G_xy = sympify(G_xy)
    except: # If it fails, print error message
        raise Exception("The entered strings were formatted wrongly. Try something like: '3*x*exp(-7*y**2)'.")

    # Defined function to simplify the notation used to evaluate the received functions F(x,y) and G(x,y)
    F = lambda p,q : funEvaluate(funEvaluate(F_xy,"x",p),"y",q)
    G = lambda p,q : funEvaluate(funEvaluate(G_xy,"x",p),"y",q)

    # No. of decimal places for the results
    precision = int(abs(log10(epsilon))) if epsilon > 0 and log10(epsilon).is_integer else 6

    # Only 2 elements will be used for each list
    x = [x_0,0]
    y = [y_0,0]

    if displayConsoleLog is True: 
            printTotalIterations = True
            print("___________________________________________________________________________\n\n")
            print(f"LINEAR ITERATIONS METHOD ::: Iteration no.00: F(x_00,y_00) = % .{precision-1}e ; G(x_00,y_00) = % .{precision-1}e ; x_00 = % .{precision}f ; y_00 = % .{precision}f"
                % ( F(x[0],y[0]), G(x[0],y[0]), x[0], y[0] )
            )

    # Auxiliar index variable correspondent to the current element in the 'x' and 'y' lists
    j = 0

    # Loops through the algorithm until the maximum number of iterations is met
    for i in range(1,MAXITER):
        j += 1 if j < 1 else -j # j = 0 --> j = 1 --> j = 0 --> j = 1 ...

        x[j] = F(x[0],y[0])
        y[j] = G(x[0],y[0])

        if displayConsoleLog is True: 
            print(f"LINEAR ITERATIONS METHOD ::: Iteration no.%s: F(x_%s,y_%s) = % .{precision-1}e ; G(x_%s,y_%s) = % .{precision-1}e ; x_%s = % .{precision}f ; y_%s = % .{precision}f"
                % ( str(i).zfill(2), str(i).zfill(2), str(i).zfill(2), F(x[0],y[0]), str(i).zfill(2), str(i).zfill(2), G(x[0],y[0]), str(i).zfill(2), x[j], str(i).zfill(2), y[j] )
            )

        # Conditions to check if it found either the actual solution or one within the specified error range, respectively
        if  abs(x[j] - x[j-1]) < epsilon * max(1,abs(x[j])) and abs(y[j] - y[j-1]) < epsilon * max(1,abs(y[j])):
            if printTotalIterations is True: 
                print("___________________________________________________________________________\n\n") if not displayConsoleLog else print("\n")
                print("LINEAR ITERATIONS METHOD ::: The solution was successfully found after %d iterations.\n" % (i+1) )
                print("___________________________________________________________________________\n\n")
            return float(x[j]), float(y[j]), int(precision)

    # If it manages to get out of the loop, it means no root was found, so it prints an error message
    raise Exception("The maximum number of iterations was met and no solution was found.")

# NEWTON'S METHOD
# Receives f(x), calculates the best possible psi(x) and then executes the Linear Iterations Method
def newtons(f_xy, g_xy, x_0: float, y_0: float, epsilon = 10**-6, MAXITER = 100, displayConsoleLog = False, printPartialDerivatives = False, printTotalIterations = False):
    try: # Try and convert the functions from string to symbolic
        f_xy = sympify(f_xy)
        g_xy = sympify(g_xy)
    except: # If it fails, print error message
        raise Exception("The entered strings were formatted wrongly. Try something like: '3*x*exp(-7*y**2)'.")

    dfx = diff(f_xy,'x') # fx(x,y) = df/dx(x,y)
    dfy = diff(f_xy,'y') # fy(x,y) = df/dy(x,y)
    dgx = diff(g_xy,'x') # gx(x,y) = dg/dx(x,y)
    dgy = diff(g_xy,'y') # gy(x,y) = dg/dy(x,y)

    # F(x,y) = x - [(f*gy - g*fy)/(fx*gy - fy*gx)](x,y)
    F_xy = sympify('x') - ((f_xy*dgy - g_xy*dfy)/(dfx*dgy - dfy*dgx))
    # G(x,y) = y - [(g*fx - f*gx)/(fx*gy - fy*gx)](x,y)
    G_xy = sympify('y') - ((g_xy*dfx - f_xy*dfx)/(dfx*dgy - dfy*dgx))

    if displayConsoleLog: 
        print("___________________________________________________________________________\n\n")
        print("NEWTON'S METHOD ::: F(x,y) and G(x,y) were calculated successfully: ")
        print("NEWTON'S METHOD ::: F(x,y) = ")
        pprint(F_xy)
        print("NEWTON'S METHOD ::: G(x,y) = ")
        pprint(F_xy)
        print()

    # Defined function to simplify the notation used to evaluate the received functions F(x,y) and G(x,y)
    f = lambda p,q : funEvaluate(funEvaluate(f_xy,"x",p),"y",q)
    g = lambda p,q : funEvaluate(funEvaluate(g_xy,"x",p),"y",q)
    fx = lambda p,q : funEvaluate(funEvaluate(dfx,"x",p),"y",q)
    fy = lambda p,q : funEvaluate(funEvaluate(dfy,"x",p),"y",q)
    gx = lambda p,q : funEvaluate(funEvaluate(dgx,"x",p),"y",q)
    gy = lambda p,q : funEvaluate(funEvaluate(dgy,"x",p),"y",q)
    F = lambda p,q : funEvaluate(funEvaluate(F_xy,"x",p),"y",q)
    G = lambda p,q : funEvaluate(funEvaluate(G_xy,"x",p),"y",q)

    # No. of decimal places for the results
    precision = int(abs(log10(epsilon))) if epsilon > 0 and log10(epsilon).is_integer else 6

    # Only 2 elements will be used for each list
    x = [x_0,0]
    y = [y_0,0]

    if displayConsoleLog: 
            printTotalIterations = True
            print("___________________________________________________________________________\n\n")
            if printPartialDerivatives: 
                print(f"LINEAR ITERATIONS METHOD ::: Iteration no.00: f(x_00,y_00) = % .{precision-1}e ; g(x_00,y_00) = % .{precision-1}e ; fx(x_00,y_00) = % .{precision-1}e ; fy(x_00,y_00) = % .{precision-1}e ; gx(x_00,y_00) = % .{precision-1}e ; gy(x_00,y_00) = % .{precision-1}e"
                    % ( f(x[0],y[0]), g(x[0],y[0]), fx(x[0],y[0]), fy(x[0],y[0]), gx(x[0],y[0]), gy(x[0],y[0]) )
                )
            print(f"LINEAR ITERATIONS METHOD ::: Iteration no.00: F(x_00,y_00) = % .{precision-1}e ; G(x_00,y_00) = % .{precision-1}e ; x_00 = % .{precision}f ; y_00 = % .{precision}f"
                % ( F(x[0],y[0]), G(x[0],y[0]), x[0], y[0] )
            )

    # Auxiliar index variable correspondent to the current element in the 'x' and 'y' lists
    j = 0

    # Loops through the algorithm until the maximum number of iterations is met
    for i in range(1,MAXITER):
        j += 1 if j < 1 else -j # j = 0 --> j = 1 --> j = 0 --> j = 1 ...

        x[j] = F(x[j-1],y[j-1])
        y[j] = G(x[j-1],y[j-1])

        if displayConsoleLog: 
            if printPartialDerivatives:
                print(f"LINEAR ITERATIONS METHOD ::: Iteration no.%s: f(x_%s,y_%s) = % .{precision-1}e ; g(x_%s,y_%s) = % .{precision-1}e ; fx(x_%s,y_%s) = % .{precision-1}e ; fy(x_%s,y_%s) = % .{precision-1}e ; gx(x_%s,y_%s) = % .{precision-1}e ; gy(x_%s,y_%s) = % .{precision-1}e"
                    % ( str(i).zfill(2), str(i).zfill(2), str(i).zfill(2), f(x[j],y[j]), str(i).zfill(2), str(i).zfill(2), g(x[j],y[j]), str(i).zfill(2), str(i).zfill(2), 
                    fx(x[j],y[j]), str(i).zfill(2), str(i).zfill(2), fy(x[j],y[j]), str(i).zfill(2), str(i).zfill(2), gx(x[j],y[j]), str(i).zfill(2), str(i).zfill(2), gy(x[j],y[j]) )
                )
            print(f"LINEAR ITERATIONS METHOD ::: Iteration no.%s: F(x_%s,y_%s) = % .{precision-1}e ; G(x_%s,y_%s) = % .{precision-1}e ; x_%s = % .{precision}f ; y_%s = % .{precision}f"
                % ( str(i).zfill(2), str(i).zfill(2), str(i).zfill(2), F(x[j],y[j]), str(i).zfill(2), str(i).zfill(2), G(x[j],y[j]), str(i).zfill(2), x[j], str(i).zfill(2), y[j] )
            )

        # Conditions to check if it found either the actual solution or one within the specified error range, respectively
        if  abs(x[j] - x[j-1]) < epsilon * max(1,abs(x[j])) and abs(y[j] - y[j-1]) < epsilon * max(1,abs(y[j])):
            if printTotalIterations is True: 
                print("___________________________________________________________________________\n\n") if not displayConsoleLog else print("\n")
                print("LINEAR ITERATIONS METHOD ::: The solution was successfully found after %d iterations.\n" % (i+1) )
                print("___________________________________________________________________________\n\n")
            return float(x[j]), float(y[j]), int(precision)

    # If it manages to get out of the loop, it means no root was found, so it prints an error message
    raise Exception("The maximum number of iterations was met and no solution was found.")

if __name__ == "__main__":
    f = "x**2 - y - 1"
    g = "y**2 - x - 1"

    x_0 = 2
    y_0 = 2
    error = 10**-5
    maxIt = 10**3
    discriminateIt = True
    totalIt = True

    xSol, ySol, precision = newtons(f, g, x_0, y_0,error,maxIt,discriminateIt,False,totalIt) # NEWTON
    print("SOLUTION:")
    print(f"x = {xSol:.{precision}f}")
    print(f"y = {ySol:.{precision}f}\n")

    F_xy = "sqrt(y+1)"
    G_xy = "sqrt(x+1)"

    xSol, ySol, precision = linearIterations(F_xy,G_xy, x_0, y_0,error,maxIt,discriminateIt,totalIt) # LINEAR ITERATIONS
    print("SOLUTION:")
    print(f"x = {xSol:.{precision}f}")
    print(f"y = {ySol:.{precision}f}\n")