#!/usr/bin/env python3
# -*- coding: utf-8 -*-


""""
Created on Tue Nov 23 14:07:35 2021
@author: cecile
"""

"""
==============================================================================
Title       : Base Conversion Algorithms in Python
Description : Implements basic number system conversions between binary, 
              decimal, octal, hexadecimal, and other bases. 
              Includes general-purpose base conversion functions and
              examples of converting both integer and fractional values.
==============================================================================
List of Functions:
------------------------------------------------------------------------------
1. base_2_8(n)
    Converts a binary number (base 2) to octal (base 8).

2. base_2_10(n)
    Converts a binary number (base 2) to decimal (base 10).

3. base_5_10(n)
    Converts a base-5 number to decimal (base 10).

4. base_8_10(n)
    Converts an octal number (base 8) to decimal (base 10).

5. base_b_10(b, n)
    Converts a number from any base b to decimal (base 10).

6. base_16_10(n)
    Converts a hexadecimal number (base 16) to decimal (base 10).

7. base_10_2(n)
    Converts a decimal number (base 10) to binary (base 2).

8. base_10_b(n, b)
    Converts a decimal number (base 10) to another base b 
    (supports digits and letters A–F for bases > 10).

9. conversion(n)
    Demonstrates conversion of a decimal number with a fractional part 
    to binary (integer and fractional parts handled separately).
==============================================================================
"""

# 1. Convert from base 2 to base 8
# Returns a string representation of the octal number.
# Slightly approximate: padding logic avoids index errors but not fully robust for irregular binary lengths.
def base_2_8(n):
    ch=str(n) 
    l=len(ch)
    q,r=divmod(l,3)
    if r==1:
        ch="00"+ch
        q+=1
    elif r==2:
        ch="0"+ch
        q+=1
    w=""
    for i in range(q):
        ch1=int(ch[3*i])*4
        ch2=int(ch[1+3*i])*2
        ch3=int(ch[2+3*i])*1
        w=w+str(ch3+ch2+ch1)
    return (w)
        
print(base_2_8(10111111))


# 2. Convert from base 2 to base 10
# Converts a binary string to its decimal value.
def base_2_10(n):
    ch=str(n)
    l=len(ch)
    n2_10=0
    for i in range(l):
        n2_10+=int(ch[-1-i])*2**i
    return (n2_10)   

print(base_2_10(1011101001110101111))


# 3. Convert from base 5 to base 10
# Handles all valid base-5 inputs.
def base_5_10(n):
    ch=str(n)
    l=len(ch)
    n5_10=0
    for i in range (l):
        n5_10+=int(ch[-1-i])*5**i
    return (n5_10)
print(base_5_10(333))


# 4. Convert from base 8 to base 10
# Works for all valid octal numbers.
def base_8_10(n):
    ch=str(n)
    l=len(ch)
    n8_10=0
    for i in range (l):
        n8_10+=int(ch[-1-i])*8**i
    return (n8_10)
print(base_8_10(1631))


# 5. Convert from any base b to base 10
# Works correctly for bases ≤ 10.
# Limitation: does not support letters for bases greater than 10.
def base_b_10(b,n):
    ch=str(n)
    l=len(ch)
    nb_10=0
    for i in range (l):
        nb_10+=int(ch[-1-i])*b**i
    return (nb_10)
print(base_b_10(5,333))


# 6. Convert from base 16 (hexadecimal) to base 10
# Handles digits and letters A–F correctly.
def base_16_10(n):
    ch=str(n)
    l=len(ch)
    n16_10=0
    L=('0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F')
    for i in range(l):
        n16_10+=L.index(ch[i])*16**(l-i-1)
    return n16_10
print(base_16_10("FC2B"))


# 7. Convert from base 10 to base 2
# Returns the binary representation as an integer.
def base_10_2(n):
    r2=""
    while n!=0:
        q,r=divmod(n,2)
        n=q
        r2=str(r)+r2
    return int(r2)

print(base_10_2(9))


# 8. Convert from base 10 to base b
# Works for bases up to 16. For bases > 16, letters are not extended beyond 'F'.
def base_10_b(n,b):
    rb=""
    while n!=0:
        q,r=divmod(n,b)
        n=q
        if r>9:
            L=['A','B','C','D','E','F']
            r=L[r-10]
        rb=str(r)+rb
    return (rb)

print(base_10_b(30,16))


# 9. Conversion of a number with fractional part to binary (integer + decimal)
# Approximate method: converts digits of the fractional part, not the actual fraction value.
# Prints intermediate steps but does not perform true binary fractional conversion.
from math import floor
def conversion(n):
    e=floor(n)
    p=str(n)
    pp=p.split('.')
    print(pp)
    ppp=int(pp[1])
    che=str(e)
    chp=str(p)
    le=len(che)
    lp=len(chp)
    re=""
    while e!=0:
        q,r=divmod(e,2)
        e=q
        re=str(r)+re
    rp=""
    while ppp!=0:
        q,r=divmod(ppp,2)
        ppp=q
        rp=str(r)+rp
        print(rp)
    return re+','+rp

print(conversion(0.48))
