#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 13:02:30 2022

@author: cecile
"""
"""
==============================================================================
Description : Implements basic array processing functions including:
              - counting elements below a threshold,
              - counting occurrences,
              - finding majority elements,
              - and performing linear and binary search.
==============================================================================
List of Functions:
------------------------------------------------------------------------------
- inferieur(x, T)
    Counts how many elements in T are strictly less than x.

- inferieurtrie(T, x)
    Counts how many elements are less than x in a sorted array using binary search.

- mult(T, x)
    Returns the number of times x appears in T.

- majoritaire(T)
    Checks if there is a majority element (appears ≥ n/2 + 1 times).

- majoritaireplus(T)
    Optimized version for sorted arrays to check for a majority element.

- recherchesequentielle(T, a, b)
    Performs a linear search for the first element between a and b.

- recherchedico(T, c, d)
    Searches for any element in range [c, d] in a sorted array using binary search.

- recherche_dicho_recursive(x, T)
    Recursively searches for x using binary search in a sorted array.

- recherche_dicho_iterative(x, T)
    Iteratively searches for x using binary search in a sorted array.
==============================================================================
"""

# Function that counts how many elements in T are strictly less than x
def inferieur(x, T):
    n = len(T)
    s = 0
    for i in range(n):
        if T[i] < x:
            s += 1
    return s


# Function that counts how many elements are less than x in a sorted array T (binary search method)
# T must be sorted in increasing order
def inferieurtrie(T, x):
    n = len(T)
    a = 0
    b = n - 1
    while a <= b:
        c = (a + b) // 2
        if T[c] < x:
            a = c + 1
        else:
            b = c - 1
    return a  # 'a' represents the number of elements < x


# Function that counts the multiplicity of an element x in array T
def mult(T, x):
    n = len(T)
    s = 0
    for i in range(n):
        if T[i] == x:
            s += 1
    return s


# Function that checks if there exists a majority element in T
# (appears at least n/2 + 1 times)
def majoritaire(T):
    n = len(T)
    for i in range(n):
        if mult(T, T[i]) >= (n // 2 + 1):
            return True
    return False


# Optimized function to determine if there is a majority element in a sorted array
# T must be sorted in increasing order
def majoritaireplus(T):
    n = len(T)
    if mult(T, T[n // 2]) >= (n // 2 + 1):
        return True
    return False


# Function that sequentially searches for the first element between a and b
def recherchesequentielle(T, a, b):
    n = len(T)
    for i in range(n):
        if T[i] >= a and T[i] <= b:
            return i
    return -1


'''Time complexity = O(n) comparisons'''


# Function that searches for an element between c and d in a sorted array (binary search)
# T must be sorted in increasing order
def recherchedico(T, c, d):
    n = len(T)
    a = 0
    b = n - 1
    if T[n - 1] < c or T[0] > d:
        return -1
    while a <= b:
        i = (a + b) // 2
        if T[i] >= c and T[i] <= d:
            return i
        elif T[i] < c:
            a = i + 1
        else:
            b = i - 1
    return -1


'''Time complexity = p because log2(2p) = p'''


# Function that performs a recursive binary search of x in T
# T must be sorted in increasing order
def recherche_dicho_recursive(x, T):
    n = len(T)
    if n == 0:
        return 0
    p = n // 2
    if T[p] == x:
        return 1
    if T[p] > x:
        return recherche_dicho_recursive(x, T[:p])
    else:
        return recherche_dicho_recursive(x, T[p + 1:])


# Function that performs an iterative binary search of x in T
# T must be sorted in increasing order
def recherche_dicho_iterative(x, T):
    n = len(T)
    if n == 0:
        return 0
    a = 0
    b = n - 1
    while a <= b:
        c = (a + b) // 2
        if T[c] == x:
            return 1
        if T[c] > x:
            b = c - 1
        else:
            a = c + 1
    return 0
