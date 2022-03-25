# -*- coding: utf-8 -*-
"""
Created on Thu Feb 03 12:12:46 2017

@author: tony.withers@uwo.ca

Functions to calculate H2O molar volume (PSvolume) and fugacity (PSfug)
using the Pitzer and Sterner equation of state.

Pitzer, K.S. and Sterner, S.M., 1994. Equations of state valid
continuously from zero to extreme pressures for H2O and CO2.
Journal of Chemical Physics. 101: 3111-3116.
"""

import math
from scipy import optimize

coeff=[]
coeff.append([0,0,0.24657688e6,0.51359951e2,0,0])
coeff.append([0,0,0.58638965e0,-0.28646939e-2,0.31375577e-4,0])
coeff.append([0,0,-0.62783840e1,0.14791599e-1,0.35779579e-3,0.15432925e-7])
coeff.append([0,0,0,-0.42719875e0,-0.16325155e-4,0])
coeff.append([0,0,0.56654978e4,-0.16580167e2,0.76560762e-1,0])
coeff.append([0,0,0,0.10917883e0,0,0])
coeff.append([0.38878656e13,-0.13494878e9,0.30916564e6,0.75591105e1,0,0])
coeff.append([0,0,-0.65537898e5,0.18810675e3,0,0])
coeff.append([-0.14182435e14,0.18165390e9,-0.19769068e6,-0.23530318e2,0,0])
coeff.append([0,0,0.92093375e5,0.12246777e3,0,0])

def PSeos(volume, temperature, targetP):  # cc/mol, Kelvins, bars
    R=8314510  # Pa.cc/K/mol
    den=1/volume  # mol/cc
    c=[]
    for i in range(10):
        c.insert(i,coeff[i][0]*temperature**-4+coeff[i][1]*temperature**-2
                 +coeff[i][2]*temperature**-1+coeff[i][3]
                 +coeff[i][4]*temperature+coeff[i][5]*temperature**2)
    pressure=(den+c[0]*den**2-den**2*((c[2]+2*c[3]*den+3*c[4]*den**2
              +4*c[5]*den**3)/(c[1]+c[2]*den+c[3]*den**2+c[4]*den**3
              +c[5]*den**4)**2)+c[6]*den**2*math.exp(-c[7]*den)
              +c[8]*den**2*math.exp(-c[9]*den))*R*temperature/1e5
    return pressure-targetP  # bars

def PSvolume(pressure, temperature):  # bars, Kelvins
    volume=optimize.root(PSeos, 10, args = (temperature, pressure))
    return volume.x

def PSfugacity(pressure, temperature):  # bars, Kelvins
    R=8314510  # Pa.cc/K/mol
    c=[]
    for i in range(10):
        c.insert(i,coeff[i][0]*temperature**-4+coeff[i][1]*temperature**-2
                 +coeff[i][2]*temperature**-1+coeff[i][3]
                 +coeff[i][4]*temperature+coeff[i][5]*temperature**2)
    volume=PSvolume(pressure, temperature)
    den=1/volume  # mol/cc
    fug=math.exp(math.log(den)+c[0]*den+(1/(c[1]+c[2]*den+c[3]*den**2
                                         +c[4]*den**3+c[5]*den**4)-1/c[1])
                 -c[6]/c[7]*(math.exp(-c[7]*den)-1)
                 -c[8]/c[9]*(math.exp(-c[9]*den)-1)
                 +pressure*1e5/(den*R*temperature)
                 +math.log(R*temperature)-1)/1e5
    return fug  # bars
