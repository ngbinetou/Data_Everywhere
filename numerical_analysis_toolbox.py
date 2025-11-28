#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 10:55:35 2022

@author: cecile

==============================================================================
Title       : Numerical Analysis — Root Finding, Interpolation, Quadrature,
              ODE Solvers, LU Factorization, and Iterative Linear Solvers
Description : Coordinated collection of basic numerical analysis routines,
              kept close to the original coding style, with minimal fixes.
==============================================================================
List of Functions:
1) Root finding
   - dichotomie(F, a, b, n)
   - newton(F, f, a0, n)
   - dichotomie2(F, a, b, e)

2) Lagrange interpolation
   - Polylagrange(X, x, i)
   - interLagrange(X, f, x)

3) Numerical quadrature
   - RectangleDroite(f, a, b, n)
   - RectangleGauche(f, a, b, n)
   - Trapeze(f, a, b, n)
   - Simpson(f, a, b, n)

4) ODE solvers (IVP)
   - EulerExplicite(f, a, b, y0, n)
   - Heun(f, a, b, y0, n)
   - RK4(f, a, b, y0, n)
   - Newton(F, f, a0, n)
   - EulerImplicite(f, ff, a, b, y0, n)
   - CranckNicholson(f, ff, a, b, y0, n)  # (name kept as in original)

5) Linear systems — LU and triangular solves
   - DecompositionLU(A)
   - systemeInferieur(L, b)
   - systemeSuperieur(U, b)
   - resolutionLU(A, b)

6) Linear systems — iterative methods
   - yjacobi(A, b, x)
   - jacobi(A, b, x0, k)
   - yGauss_Seidel(A, b, x)
   - Gauss_Seidel(A, b, x0, k)
   - Gradient(A, b, x0, alf, k)  # Richardson iteration
==============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si


# =============================================================================
# 1) ROOT FINDING
# =============================================================================

# Bisection method over [a, b] for n iterations.
# Note: For correctness, the initial bracket should satisfy F(a)*F(b) < 0.
def dichotomie(F,a,b,n):
    for i in range (n+1):
        d=(a+b)/2
        if F(d)> 0:
             b=d
        else:
            a=d
    return d     

# Newton's method starting from a0, for n iterations.
# Note: f(a0) should not be zero to avoid division-by-zero.
def newton(F,f,a0,n):
    for i in range (n):
        a0=a0-F(a0)/f(a0)
    return a0

# Bisection with target precision e (number of iterations computed by ceil).
def dichotomie2(F,a,b,e):
    n=int(np.ceil(np.log2((b-a)/e)))
    return dichotomie(F, a, b, n)


# --- Example for sqrt(2) ---
def F(x):
    return x**2-2
def f(x):
    return 2*x
x=2**0.5
a=1
b=3
a0=2
n=30
print("Newton sqrt(2):", newton(F, f, a0, n))
print("Bisection sqrt(2):", dichotomie(F, a, b, n))
print("Reference sqrt(2):", x)


# =============================================================================
# 2) LAGRANGE INTERPOLATION
# =============================================================================

# Lagrange basis polynomial ℓ_i(x) for nodes X.
def Polylagrange(X,x,i):
    q=1
    for j in range (len(X)):
      if j!=i:
          q=q* (x-X[j]) / (X[i]-X[j])
    return q

# Lagrange interpolation polynomial at x, using samples f(X[i]).
def interLagrange(X,f,x):
    p=0
    for i in range(len(X)):
        q=Polylagrange(X,x,i)
        p+=f(X[i])*q
    return p

# --- Small demo plot ---
a,b=0,2*np.pi
n=3
p=5
X=np.linspace(a,b,n+1)
T=np.linspace(a,b,p*(n+1))
images_sin=np.sin(T)
POLY=[interLagrange(X,np.sin,a) for a in T]
plt.figure()
plt.plot(T,images_sin,label="sin(x)")
plt.plot(T,POLY,label="Lagrange interp.", linestyle="--")
plt.legend()
plt.title("Lagrange interpolation of sin on [0, 2π]")


# =============================================================================
# 3) NUMERICAL QUADRATURE
# =============================================================================

# Right Riemann sum on n subintervals.
def RectangleDroite(f,a,b,n):
    h=(b-a)/n
    ak=a
    s=0
    for k in range (n):
        ak= ak+h
        s=s+f(ak)
    Inm=h*s
    return Inm

# Left Riemann sum on n subintervals.
def RectangleGauche(f,a,b,n):
    h=(b-a)/n
    ak=a
    s=0
    for k in range (n):
        s=s+f(ak)
        ak= ak+h
    Inp=h*s
    return Inp

# Trapezoidal rule on n subintervals.
def Trapeze(f,a,b,n):
    h=(b-a)/n
    ak=a
    s=0
    for k in range (n-1):
        ak= ak+h
        s=s+f(ak)
    Tn= h*((f(a)+f(b))/2+s)
    return Tn

# Simpson's rule on n subintervals (composite).
# Note: classical Simpson typically expects n even for best consistency.
def Simpson(f,a,b,n):
    h=(b-a)/n
    ak=a
    zk=a+h/2
    s=0
    s1=f(zk)
    for k in range (n-1):
        ak= ak+h
        zk= zk+h
        s= s+f(ak)
        s1= s1+f(zk)
    Sn=(h/6)*(f(a)+f(b)+2*s+4*s1)
    return Sn

# --- Quadrature demo vs scipy.integrate.quad ---
def f_gauss(t):
    return np.exp((-t**2)/2)
a=0
b=1
n=100
print('Rectangle Gauche:', RectangleGauche(f_gauss,a,b,n))
print('quad:', si.quad(f_gauss,a,b), '\n')

print('Rectangle Droite:', RectangleDroite(f_gauss,a,b,n))
print('quad:', si.quad(f_gauss,a,b), '\n')

print('Trapeze:', Trapeze(f_gauss,a,b,n))
print('quad:', si.quad(f_gauss,a,b), '\n')

print('Simpson:', Simpson(f_gauss,a,b,n))
print('quad:', si.quad(f_gauss,a,b), '\n')


# =============================================================================
# 4) ODE SOLVERS
# =============================================================================

# Explicit Euler method on [a,b] with n steps.
def EulerExplicite(f,a,b,y0,n):
    h=(b-a)/n
    yk=y0
    Y=[yk]
    T=np.linspace(a,b,n+1)
    for k in range(n):
        yk=yk+h*f(T[k],yk)
        Y+=[yk]
    return T,Y

# Heun's method (explicit trapezoidal) — fixed f calls (f(t,y)).
def Heun(f,a,b,y0,n):
    h=(b-a)/n
    yk=y0
    Y=[yk]
    T=np.linspace(a,b,n+1)
    for k in range(n):
        k1=yk+h*f(T[k],yk)
        yk=yk+(h/2)*(f(T[k],yk)+f(T[k+1],k1))
        Y+=[yk]
    return T,Y

# Classical Runge–Kutta 4 (RK4) — fixed f calls (f(t,y)).
def RK4 (f,a,b,y0,n):
    h=(b-a)/n
    yk=y0
    Y=[yk]
    T=np.linspace(a,b,n+1)
    for k in range(n):
        k1=f(T[k], yk)
        k2=f(T[k]+h/2, yk+(h/2)*k1)
        k3=f(T[k]+h/2, yk+(h/2)*k2)
        k4=f(T[k]+h, yk+h*k3)
        yk=yk+(h/6)*(k1+2*k2+2*k3+k4)
        Y+=[yk]
    return T,Y

# Newton helper for implicit schemes (kept from original style).
def Newton(F,f,a0,n):
    for i in range (n):
        a0=a0-F(a0)/f(a0)
    return a0

# Implicit Euler via Newton iteration; ff = ∂f/∂y (Jacobian or partial derivative).
def EulerImplicite(f,ff,a,b,y0,n):
    yk=y0
    Y=[yk]
    h=(b-a)/n
    T=np.linspace(a,b,n+1)
    for k in range(n):
        g=lambda x: x-yk-h*f(T[k+1],x)
        g1=lambda x: 1-h*ff(T[k+1],x)
        a0=yk
        yk=Newton(g,g1,a0,10)
        Y+=[yk]
    return T,Y

# Crank–Nicolson (name kept CranckNicholson).
def CranckNicholson(f,ff,a,b,y0,n):
    yk=y0
    Y=[yk]
    h=(b-a)/n
    T=np.linspace(a,b,n+1)
    for k in range(n):
        l=lambda x: x-yk-((h/2)*(f(T[k],yk)+f(T[k+1],x)))
        l1=lambda x: 1-(h/2)*ff(T[k+1],x)
        a0=yk
        yk=Newton(l,l1,a0,10)
        Y+=[yk]
    return T,Y


# =============================================================================
# 5) LINEAR SYSTEMS — LU & TRIANGULAR
# =============================================================================

# LU decomposition without pivoting (A = L U).
# Note: requires non-zero pivots U[i,i].
def DecompositionLU (A):
    n = A.shape[0]    
    U = np.copy(A)  
    L = np.eye(n)
    for i in range (n):
        c=U[i,i]
        for j in range (i+1,n):
            L[j,i] = U[j,i] / c
            U[j] = U[j]- U [i]*L [j,i]
    return L , U

# Forward substitution: solve L x = b.
def systemeInferieur (L,b):
    n= L.shape[0]
    x= [0 for i in range (n)]
    for i in range (n):
        s=0
        for j in range (i):
            s+=L [i,j]*x [j]
        x[i]= (b[i]-s)/L[i,i]
    return x

# Backward substitution: solve U x = b.
def systemeSuperieur (U,b):
    n= U.shape[0]
    x= [0 for i in range (n)]
    for i in range (1,n+1):
        s=0
        for j in range (n-i+1,n):
            s+=U[n-i,j]*x[j]
        x[n-i]=(b[n-i]-s)/U[n-i,n-i]
    return x

# Solve A x = b via LU factorization.
def resolutionLU(A,b):
    L,U = DecompositionLU(A)
    y = systemeInferieur(L, b)
    x = systemeSuperieur(U, y)
    return x

# --- small LU demo ---
A=np.array([[1,2,3,4],[2,3,4,1],[3,4,1,2],[4,1,2,3]])
b=  [1,2,3,4]
x= resolutionLU(A,b)
print("LU solve x:", x)
print("A*e1:", np.dot(A,np.array([[1], [0], [0], [0]])))
print("A*x:", np.dot(A, np.transpose(x)))


# =============================================================================
# 6) LINEAR SYSTEMS — ITERATIVE METHODS
# =============================================================================

# One Jacobi sweep: y = D^{-1}(b - (L+U) x)
def yjacobi(A,b,x):
    n=len(b)
    y=[0 for i in range (n)]
    for i in range (n):
        s=0
        for j in range (n):
            if (j!=i):  
                s+=A[i,j]*x[j]
        y[i]=1/A[i,i]*(b[i]-s)
    return y

# Jacobi iteration for k steps starting at x0.
def jacobi(A,b,x0,k):
    x=[a for a in x0]
    for i in range (k):
        x= yjacobi(A, b, x)
    return x

# One Gauss–Seidel sweep.
def yGauss_Seidel(A,b,x):
    n=len(b)
    y=[0 for i in range (n)]
    for i in range (n):
        s1=0
        s2=0
        for j in range(i):
            s1+=A[i,j]*y[j]
        for j in range (i+1,n):
            s2+=A[i,j]*x[j]
        y[i]=1/A[i,i]*(b[i]-s1-s2)
    return y

# Gauss–Seidel iteration for k steps starting at x0.
def Gauss_Seidel(A,b,x0,k):
    x=[a for a in x0]
    for i in range (k):
        x= yGauss_Seidel(A, b, x)
    return x

# Richardson / Gradient-like fixed-step iteration: x_{k+1} = (I - αA) x_k + α b
def Gradient(A,b,x0,alf,k):
    n=len(b)
    x=np.array(x0)
    B=np.eye(n)-alf*np.array(A)
    bvec=np.array(b)
    for i in range (k):
        x=np.dot(B,x)+alf*bvec
    return x
