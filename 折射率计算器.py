# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 10:52:15 2022

@author: psypple
"""

lam=520E-9
h = 4.136E-15
c = 3E8
x=float(input('In content= '))
y=float(input('Al content= '))
EgAl=3.42*(1-y)+6.13*y-1.3*y*(1-y)
aAl=9.82661-8.21608*y+31.5902*y**2
bAl=2.73591+0.84249*y-6.29321*y**2
nalgan=(aAl*(h*c/lam/EgAl)**(-2)*(2-(1+h*c/lam/EgAl)**0.5-(1-h*c/lam/EgAl)**0.5)+bAl)**0.5

EgIn=3.42*(1-x)+0.77*x-1.43*x*(1-x)
aIn=9.82661
bIn=2.73591
ningan=(aIn*(h*c/(lam+(1.24/3.42*1E-6-1.24/EgIn*1E-6))/3.42)**(-2)*(2-(1+h*c/(lam+(1.24/3.42*1E-6-1.24/EgIn*1E-6))/3.42)**0.5-(1-h*c/(lam+(1.24/3.42*1E-6-1.24/EgIn*1E-6))/3.42)**0.5)+bIn)**0.5

print ('n_InGaN=%f'%ningan)
print ('n_AlGaN=%f'%nalgan)