""""
Created on Thu Feb  3 15:46:14 2022
@author: cecile
"""
""""
==============================================================================
Description : Implements a collection of basic numerical and algorithmic 
              exercises in Python, including:
              - number parity and perfect number check
              - variable swapping methods
              - maximum and sum of numbers
              - convolution of two lists
              - quadratic equation solver
              - generation of prime numbers (Sieve of Eratosthenes)
              - different power computation methods
==============================================================================
List of Functions and Exercises:
------------------------------------------------------------------------------
1. parite(n)
    Determines whether a number is even (returns 1) or odd (returns 0).

2. Variable Swapping (Two Methods)
    Method 1: Using a temporary variable.
    Method 2: Without using a temporary variable.

3. Perfect Number Check
    Determines whether a given integer is a perfect number.

4. Maximum Value in a List
    Finds and displays the largest number in a list of user inputs.

5. Discrete Convolution of Two Lists
    Performs the convolution (element-wise product sum) between lists A and B.

6. Sum of the First n Integers
    Computes the sum using both an iterative loop and the direct formula.

7. Quadratic Equation Solver
    Solves a second-degree equation ax² + bx + c = 0 and displays all roots.

8. eratosthene(n)
    Returns a list of all prime numbers less than or equal to n using 
    the Sieve of Eratosthenes algorithm.

9. puiss1(x, n)
    Computes xⁿ using an iterative method.

10. puiss2(x, r)
    Computes xʳ using a simple recursive approach.

11. puiss3(x, n)
    Computes xⁿ using fast recursive exponentiation.
==============================================================================
"""


# 1. Check if a number is even
def parite(n):
    """Return 1 if n is even, else 0."""
    if n % 2 == 0:
        return 1
    return 0

print(parite(22378))


# 2. Swap two numbers (two different methods)
# Method 1: Using a temporary variable
x = float(input('x: '))
y = float(input('y: '))
c = x
x = y
y = c
print('After swapping: x =', x, ', y =', y)

# Method 2: Without using a temporary variable
x = float(input('x: '))
y = float(input('y: '))
y = y + x
x = y - x
y = y - x
print('After swapping: x =', x, ', y =', y)



# 3. Check if a number is perfect
n = int(input("Enter an integer: "))
s = 0
for i in range(1, n):
    if n % i == 0:
        s += i

if s == n:
    print("%d is a perfect number" % n)
else:
    print("%d is not a perfect number" % n)




# 4. Find the maximum value in a list
T = []
N = int(input('Enter the number of elements N: '))
for i in range(N):
    T += [float(input('Enter next number: '))]
maxi = T[0]
for x in T:
    if x > maxi:
        maxi = x
print('The maximum value is:', maxi)



# 5. Discrete convolution of two lists
tau = int(input('Enter tau: '))
A = [int(input('Enter a value for A: ')) for i in range(tau)]
B = [int(input('Enter a value for B: ')) for i in range(tau)]
P = [0 for i in range(2 * tau - 1)]
for k in range(2 * tau - 1):
    for i in range(tau):
        j = k - i
        if 0 <= j < tau:
            P[k] += A[i] * B[j]
print("Convolution result:", P)



# 6. Sum of first n integers
n = int(input('Enter n: '))
S = 0
for i in range(n + 1):
    S += i
print("Sum (with loop):", S)

n = int(input('Enter n again: '))
S = n * (n + 1) / 2
print("Sum (with formula):", S)



# 7. Solve a quadratic equation
from math import sqrt

A = float(input("Enter A: "))
B = float(input("Enter B: "))
C = float(input("Enter C: "))

if A == 0:
    if B == 0:
        if C == 0:
            s = "Infinite number of solutions"
        else:
            s = "No solution"
    else:
        s = "x = " + str(-C / B)
else:
    delta = B**2 - 4*A*C
    if delta == 0:
        s = "x1 = x2 = " + str(-B / (2*A))
    elif delta < 0:
        s = "No real solutions"
    else:
        s = "x1 = " + str((-B - sqrt(delta)) / (2*A)) + "\\nx2 = " + str((-B + sqrt(delta)) / (2*A))

print("Solution:", s)



# 8. Sieve of Eratosthenes: list of prime numbers up to n
def eratosthene(n):
    """Return a list of all prime numbers ≤ n using the Sieve of Eratosthenes."""
    if n < 2:
        return []
    listt = list(range(n + 1))
    listt[0] = 0
    listt[1] = 0
    k = 2
    while k * k <= n:
        if listt[k] != 0:
            i = k * k
            while i <= n:
                listt[i] = 0
                i += k
        k += 1
    listpremier = [k for k in listt if k != 0]
    return listpremier

print(eratosthene(50))


# 9. Power functions
# Iterative method
def puiss1(x, n):
    """Return x^n using an iterative method."""
    p = 1
    for i in range(n):
        p = p * x
    return p

print(puiss1(5, 2))


# Recursive method (simple)
def puiss2(x, r):
    """Return x^r using a simple recursive method."""
    if r == 0:
        return 1
    else:
        return x * puiss2(x, r - 1)

print(puiss2(5, 2))


# Recursive method (fast exponentiation)
def puiss3(x, n):
    """Return x^n using fast recursive exponentiation."""
    if n == 0:
        return 1
    elif n % 2 == 0:
        return puiss3(x, n // 2) * puiss3(x, n // 2)
    else:
        return puiss3(x, (n - 1) // 2) * puiss3(x, (n - 1) // 2) * x

print(puiss3(5, 3))
