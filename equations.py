#File containing equations for ICIC data analysis workshop 2023
#AUTHOR: ecarter | Aug 2023

#Import relevant modules
import numpy as np

#Equations defined as functions, intent is to feed one into another as needed
def calc_mu(z,h,D_l):
    #Equation in form c - a + b, where c is a constant.
    a = 5*(np.log10(h))
    b = 5*(np.log10(D_l))
    mu = 25 - a + b

    return mu

def calc_D_l(z,eta1, etaz):
    #Equation in form a*[b+c]
    a = 3000*(1+z)

    #print(a)

    #print(eta1-etaz)

    D_l = a*(eta1-etaz)

    return D_l

def calc_eta(s, a):
    #Equation in form d*[b]^-1/8
    d = 2*(np.sqrt(((s**3)+1)))

    #print(d)
    b = (1/(a**4)) - (0.1540*((s)/(a**3))) + (0.4304*((s**2)/(a**2))) + (0.19097*((s**3)/(a))) + (0.066941*(s**4))
    #print(b)
    eta = d*(b**(-(1/8)))

    return eta

def calc_s(omega):
    #Calc s cubed first as per (4) in script
    s3 = (1-omega)/omega #s**3

    s = s3**(1/3) #Cube root for s alone

    return s

def calc_mu_across_z(omega,h,z):

    distmods = []

    #Calculate once for each omega_m:
    s    = calc_s(omega)
    eta1 = calc_eta(s, 1)
    #For each value of z generated in range 0>z>2:
    for zj in z:
        #Calculate for every z:
        opz = 1+zj #To circumvent a bug in eta_opz
        eta_opz = calc_eta(s, 1/opz)
        D_l  = calc_D_l(zj, eta1,eta_opz)
        mu = calc_mu(zj,h,D_l)
        distmods.append(mu)

    return distmods
