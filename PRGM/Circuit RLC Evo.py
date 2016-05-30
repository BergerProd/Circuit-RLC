# Quentin Bergé et Benjamin Deprat
# Micro Projet Informatique
# Circuit RLC en Regime Libre et Forcé

import matplotlib.pyplot as plt
import scipy.integrate as spi
from math import sqrt,exp
import numpy as np
import time as tm

#Constantes Du Circuit
L = 0.01
C = 0.00001
R=[160,63,3]
W0 = 1/sqrt(L*C)
q0 = 0.0001


#########################
#Resolution Regime Libre#
#########################

#Méthode Euler
def Euler_Aperiodique(n):
    start = tm.time()
    h=(0.01-0)/n
    t,y,z=0,0.0001,0
    Teuler=[0]
    Yeuler_ape=[0.0001]
    Zeuler_ape=[0]
    R0 = R[0]
    Q = (1/R0)*sqrt(L/C) #Ici Q=0,2 environ, cas d'un regime Aperiodique
    for i in range (n):
        t,y,z=t+h,y+h*z,z+h*(-((W0/Q)*z)-(W0**2)*y)
        Teuler.append(t)
        Yeuler_ape.append(y)
        Zeuler_ape.append(z)
    end = tm.time()
    print ("Euler Aperiodique =",end - start)
    global Teuler,Yeuler_ape,Zeuler_ape    #globalisation des listes pour les tracer après

# W0 vaut environ 3162 rad/s, donc T = 0.002 s
# on met donc un intervalle de 5 périodes environ
# alors ti = 0, tf = 0,01
Euler_Aperiodique(10000)


def Euler_Critique(n):
    start = tm.time()
    h=(0.01-0)/n
    t,y,z=0,0.0001,0
    Teuler_crit=[0]
    Yeuler_crit=[0.0001]
    Zeuler_crit=[0]
    R1 = R[1]
    Q = (1/R1)*sqrt(L/C) #Ici Q=0,5, cas du regime Critique
    for i in range (n):
        t,y,z=t+h,y+h*z,z+h*(-((W0/Q)*z)-(W0**2)*y)
        Teuler_crit.append(t)
        Yeuler_crit.append(y)
        Zeuler_crit.append(z)
    end = tm.time()
    print("Euler Critique =",end-start)
    global Yeuler_crit,Zeuler_crit,Teuler_crit

Euler_Critique(10000)

def Euler_Pseudo(n):
    start = tm.time()
    h=(0.01-0)/n
    t,y,z=0,0.0001,0
    Teuler_pseudo=[0]
    Yeuler_pseudo=[0.0001]
    Zeuler_pseudo=[0]
    R2 = R[2]
    Q = (1/R2)*sqrt(L/C) #Q=10,5, Cas du regime Pseudo - Periodique
    for i in range (n):
        t,y,z=t+h,y+h*z,z+h*(-((W0/Q)*z)-(W0**2)*y)
        Teuler_pseudo.append(t)
        Yeuler_pseudo.append(y)
        Zeuler_pseudo.append(z)
    end = tm.time()
    print("Euler Pseudo - Periodique = ",end-start)
    global Yeuler_pseudo,Zeuler_pseudo,Teuler_pseudo

Euleur_Pseudo(10000)

#Méthode odeint

def FuncSpi(q,T):
    return (q[1],-(W0**2)*q[0]-(W0/Q)*q[1])
# définition de la fonction pour odeint selon discretisation du problème
# Eq differentielle du 2nd Ordre, alors retour d'un vecteur


def Ode_Aperiodique(n):
    start = tm.time()
    R0 = R[0]
    Q = (1/R0)*sqrt(L/C)
    global Q
    Tode_ape=np.linspace(0,0.01,n+1)
    ape=spi.odeint(FuncSpi,(0.0001,0),Tode_ape)
    Yode_ape=ape[:,0]
    end = tm.time()
    print("Odeint Critique = ",end-start)
    global Yode_ape,Tode_ape

#On defini n=10000 pour meilleure précision
Ode_Aperiodique(10000)

def Ode_Critique(n):
    start = tm.time()
    R1= R[1]
    Q = (1/R1)*sqrt(L/C)
    global Q
    Tode_crit=np.linspace(0,0.01,n+1)
    crit=spi.odeint(FuncSpi,(0.0001,0),Tode_crit)
    Yode_crit=crit[:,0]
    end = tm.time()
    print("Odeint Critique = ",end-start)
    global Yode_crit,Tode_crit

Ode_Critique(10000)

def Ode_Pseudo(n):
    start = tm.time()
    R2 = R[2]
    Q = (1/R2)*sqrt(L/C)
    global Q
    Tode_pseudo=np.linspace(0,0.01,n+1)
    pseudo=spi.odeint(FuncSpi,(0.0001,0),Tode_pseudo)
    Yode_pseudo=pseudo[:,0]
    end = tm.time()
    print("Odeint Pseudo Periodique =",end-start)
    global Yode_pseudo,Tode_pseudo

Ode_Pseudo(10000)

def Exacte_ape(n):
    start = tm.time()
    R0 = R[0]
    Q = (1/R0)*sqrt(L/C)
    r1=(-W0/(2*Q))+ (W0*sqrt((1/(2*Q)**2-1)))
    r2=(-W0/(2*Q)) - (W0*sqrt((1/(2*Q)**2-1)))
    A=-(q0*r2)/(r1-r2)
    B=(q0*r2)/(r1-r2)
    Yexa_ape=[]
    Texa_ape=np.linspace(0,0.01,n+1)
    for i in Texa_ape:
        q = A*(exp(r1*i))+B*(exp(r2*i))
        Yexa_ape.append(q)
    end = tm.time()
    print("Exacte Aperiodique = ",end-start)
    global Yexa_ape,Texa_ape

Exacte_ape(10000)

def Exacte_crit(n):
    start = tm.time()
    R1 = R[1]
    Q = (1/R1)*sqrt(L/C)
    A = q0
    B = W0*q0
    Yexa_crit=[]
    Texa_crit=np.linspace(0,0.01,n+1)
    for i in Texa_crit:
        q = (A+B*i)*(exp(-W0*i))
        Yexa_crit.append(q)
    end = tm.time()
    print("Exacte Critique = ",end-start)
    global Yexa_crit,Texa_crit

Exacte_crit(10000)

def Exacte_pseudo(n):
    start = tm.time()
    R2 = R[2]
    Q = (1/R2)*sqrt(L/C)
    W = W0 * sqrt(1-(1/(2*Q))**2)
    Yexa_pseudo=[]
    Texa_pseudo=np.linspace(0,0.01,n+1)
    for i in Texa_pseudo:
        q = q0*(exp(-W0/(2*Q)*i))*(np.cos(W*i)+(1/((2*Q)*sqrt(1-((1/(2*Q))**2))))*np.sin(W*i))
        Yexa_pseudo.append(q)
    end = tm.time()
    print ("Exacte Pseudo Periodique =",end-start)
    global Yexa_pseudo,Texa_pseudo

Exacte_pseudo(10000)


#Traçage Des Courbes
#Pour un n=10000 qui est assez élevé, les courbes se superposent quasiment

def Superpo_ape():
    p1=plt.plot(Teuler,Yeuler_ape,"r",label="Euler")
    p2=plt.plot(Tode_ape,Yode_ape,"b",label="Odeint")
    p3=plt.plot(Texa_ape,Yexa_ape,"g",label="Exacte")
    plt.xlabel("Temps")
    plt.ylabel("Charge")
    plt.legend()
    plt.grid(True)
    plt.title("Aperiodique Q<0.5")
    plt.savefig("Aperiodique Libre.png")

def Superpo_crit():
    p1=plt.plot(Teuler_crit,Yeuler_crit,"r",label="Euler")
    p2=plt.plot(Tode_crit,Yode_crit,"b",label="Odeint")
    p3=plt.plot(Texa_crit,Yexa_crit,"g",label="Exacte")
    plt.xlabel("Temps")
    plt.legend()
    plt.ylabel("Charge")
    plt.grid(True)
    plt.title("Critique Q=0.5")
    plt.savefig("Critique Libre.png")

def Superpo_pseu():
    p1=plt.plot(Teuler_pseudo,Yeuler_pseudo,"r",label="Euler")
    p2=plt.plot(Tode_pseudo,Yode_pseudo,"b",label="Odeint")
    p3=plt.plot(Texa_pseudo,Yexa_pseudo,"g",label="Exacte")
    plt.grid(True)
    plt.legend()
    plt.xlabel("Temps")
    plt.ylabel("Charge")
    plt.title("Pseudo Periodique Q>0.5")
    plt.savefig("Pseudo Periodique Libre.png")


##########################
#Résolution Régime Forcé#
#########################
