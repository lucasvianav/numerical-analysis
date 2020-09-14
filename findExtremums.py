from sympy.solvers import solve
from sympy.abc import x
from sympy import *

def getDeepDotQuality(func, arg, val, n = 3):
  dy = func.diff(arg)
  dyn = dy.subs(arg, val).evalf()
  if (dyn == 0):
    return getDeepDotQuality(dy, arg, val, n+1)
  elif (n % 2 == 1):
    return 'has an inflection point'
  elif (dyn > 0):
    return 'is min'
  else:
    return 'is max'
  return 'aaaaaa'

def getDotQuality(func, arg, val):
  dy = func.subs(arg, val).evalf()
  if (dy > 0):
    return 'is min'
  elif (dy < 0):
    return 'is max'
  else:
    return getDeepDotQuality(func, arg, val)

def findExtremums(func, arg):
  dy = func.diff(arg)
  ddy = dy.diff(arg)
  extremums = solve(dy, arg)

  maxValues = []
  minValues = []
  wat = []

  for val in extremums:
    if getDotQuality(ddy, arg, val) == "is max":
      val = float(simplify(sympify(val)).evalf())
      maxValues.append(val)
    elif getDotQuality(ddy, arg, val) == "is min":
      val = float(simplify(sympify(val)).evalf())
      minValues.append(val)
    elif getDotQuality(ddy, arg, val) == "has an inflection point":
      val = float(simplify(sympify(val)).evalf())
      wat.append(val)

  return maxValues, minValues, wat


if __name__ == "__main__":
  # findExtremums(x**2, x)
  # findExtremums(x**3 - 2*x**2 + x + 1, x)
  # findExtremums(2*x**4, x)
  # findExtremums(2*x**3, x)
  max1, min1 = findExtremums(E**x**2-4*x,x)
  print(max1)
  #print(type(max1[0]))
  print(min1)
  print(type(min1[0]))