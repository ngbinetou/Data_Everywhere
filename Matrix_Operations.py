'''
==============================================================================
Description : Implements common matrix operations such as:
              - addition, subtraction, multiplication
              - transposition, trace computation
              - triangular extraction
              - 3x3 determinant calculation
==============================================================================
List of Functions:
------------------------------------------------------------------------------
1. additiionmatrice(T1, T2)
    Adds two matrices of the same size.

2. soustractionmatrice(T1, T2)
    Subtracts the second matrix from the first, element-wise.

3. multiplicationmatrice(T1, T2)
    Multiplies two square matrices of the same size.

4. multiplicationmatrice_generale(T1, T2)
    Multiplies two matrices with compatible dimensions (general case).

5. transposeematrice(T1)
    Returns the transpose of the input matrix.

6. sommeElDiagonaux(T)
    Computes the sum of diagonal elements (trace) of a square matrix.

7. TriSup(T)
    Extracts the upper triangular part of a matrix.

8. TriInf(T)
    Extracts the lower triangular part of a matrix.

9. det3(m)
    Computes the determinant of a 3×3 matrix.

10. LU decomposition of a matrix without pivoting.
==============================================================================
'''


'''1. Addition of two matrices of the same size'''
def additiionmatrice(T1, T2):
    n = len(T1)
    m = len(T1[0])
    T3 = []
    for i in range(n):
        ligne = []
        for j in range(m):
            ligne.append(T1[i][j] + T2[i][j])
        T3.append(ligne)
    return T3


'''2. Subtraction of two matrices of the same size'''
def soustractionmatrice(T1, T2):
    n = len(T1)
    m = len(T1[0])
    T3 = []
    for i in range(n):
        ligne = []
        for j in range(m):
            ligne.append(T1[i][j] - T2[i][j])
        T3.append(ligne)
    return T3


'''3. Multiplication of two square matrices of the same size'''
def multiplicationmatrice(T1, T2):
    n = len(T1)
    T3 = []
    for i in range(n):
        ligne = []
        for j in range(n):
            s = 0
            for k in range(n):
                s += T1[i][k] * T2[k][j]
            ligne.append(s)
        T3.append(ligne)
    return T3


'''4. General matrix multiplication (compatible dimensions)'''
def multiplicationmatrice_generale(T1, T2):
    n = len(T1)
    m = len(T1[0])
    o = len(T2[0])
    T3 = []
    for i in range(n):
        ligne = []
        for j in range(o):
            s = 0
            for k in range(m):
                s += T1[i][k] * T2[k][j]
            ligne.append(s)
        T3.append(ligne)
    return T3


'''5. Transpose of a matrix'''
def transposeematrice(T1):
    n = len(T1)
    m = len(T1[0])
    T2 = []
    for j in range(m):
        ligne = []
        for i in range(n):
            ligne.append(T1[i][j])
        T2.append(ligne)
    return T2


'''6. Sum of diagonal elements (trace) of a square matrix'''
def sommeElDiagonaux(T):
    n = len(T)
    s = 0
    for i in range(n):
        s += T[i][i]
    return s


'''7. Extraction of the upper triangular matrix'''
def TriSup(T):
    n = len(T)
    m = len(T[0])
    T1 = []
    for i in range(n):
        ligne = []
        for j in range(m):
            if i > j:
                ligne.append(0)
            else:
                ligne.append(T[i][j])
        T1.append(ligne)
    return T1


'''8. Extraction of the lower triangular matrix'''
def TriInf(T):
    n = len(T)
    m = len(T[0])
    T1 = []
    for i in range(n):
        ligne = []
        for j in range(m):
            if i < j:
                ligne.append(0)
            else:
                ligne.append(T[i][j])
        T1.append(ligne)
    return T1


'''9. Determinant of a 3×3 matrix'''
def det3(m):
    return (m[0][0] * m[1][1] * m[2][2] +
            m[0][1] * m[1][2] * m[2][0] +
            m[0][2] * m[1][0] * m[2][1] -
            m[0][2] * m[1][1] * m[2][0] -
            m[0][1] * m[1][0] * m[2][2] -
            m[0][0] * m[1][2] * m[2][1])



import numpy as np
'''10. LU decomposition of a matrix without pivoting'''
def LU(A):
    n = A.shape[0]
    U = np.copy(A).astype(float)  # S'assurer que les divisions soient en flottant
    L = np.eye(n)

    for i in range(n):
        p = U[i, i]
        if np.isclose(p, 0):  # Éviter division par zéro
            raise ZeroDivisionError(f"Pivot nul à la position ({i}, {i})")
        for j in range(i + 1, n):
            L[j, i] = U[j, i] / p
            U[j] = U[j] - L[j, i] * U[i]

    return L, U
