# Quentin Bergé et Benjamin Deprat
# Micro Projet Informatique
# Circuit RLC en Regime Libre et Forcé

import matplotlib.pyplot as plt
import scipy.integrate as spi
from math import sqrt,exp
import numpy as np
import time as tm

#Resolution Regime Libre

#Méthode Euler
def Qape(n):
    start = tm.time()
    h=(0.01-0)/n
    t,y,z=0,0.0001,0
    T1=[0]
    Y1=[0.0001]
    Z1=[0]
    L = 0.01
    C = 0.00001
    R = 160
    W0 = 1/sqrt(L*C)
    Q = (1/R)*sqrt(L/C) #Ici Q=0,2 environ, cas d'un regime Aperiodique
    for i in range (n):
        t,y,z=t+h,y+h*z,z+h*(-((W0/Q)*z)-(W0**2)*y)
        T1.append(t)
        Y1.append(y)
        Z1.append(z)
    end = tm.time()
    print (end - start)
    global T1,Y1,Z1     #globalisation des listes pour les tracer après

# W0 vaut environ 3162 rad/s, donc T = 0.002 s
# on met donc un intervalle de 5 périodes environ
# alors ti = 0, tf = 0,01
Qape(10000)

def Qcrit(n):
    start = tm.time()
    h=(0.01-0)/n
    t,y,z=0,0.0001,0
    T2=[0]
    Y2=[0.0001]
    Z2=[0]
    L = 0.01
    C = 0.00001
    R = 63
    W0 = 1/sqrt(L*C)
    Q = (1/R)*sqrt(L/C) #Ici Q=0,5, cas du regime Critique
    for i in range (n):
        t,y,z=t+h,y+h*z,z+h*(-((W0/Q)*z)-(W0**2)*y)
        T2.append(t)
        Y2.append(y)
        Z2.append(z)
    end = tm.time()
    print(end-start)
    global Y2,Z2

Qcrit(10000)

def Qpseu(n):
    start = tm.time()
    h=(0.01-0)/n
    t,y,z=0,0.0001,0
    T3=[0]
    Y3=[0.0001]
    Z3=[0]
    L = 0.01
    C = 0.00001
    R = 3
    W0 = 1/sqrt(L*C)
    Q = (1/R)*sqrt(L/C) #Q=10,5, Cas du regime Pseudo - Periodique
    for i in range (n):
        t,y,z=t+h,y+h*z,z+h*(-((W0/Q)*z)-(W0**2)*y)
        T3.append(t)
        Y3.append(y)
        Z3.append(z)
    end = tm.time()
    print(end-start)
    global Y3,Z3

# Test pour Q=10000
Qpseu(10000)

#Méthode odeint

def FuncSpi(q,T):
    return (q[1],-(W0**2)*q[0]-(W0/Q)*q[1])
# définition de la fonction pour odeint selon discretisation du problème
# Eq differentielle du 2nd Ordre, alors retour d'un vecteur

def Qape_ode(n):
    start = tm.time()
    L = 0.01
    C = 0.00001
    R = 160
    W0 = 1/sqrt(L*C)
    Q = (1/R)*sqrt(L/C)
    global W0,Q
    Ta=np.linspace(0,0.01,n+1)
    ape=spi.odeint(FuncSpi,(0.0001,0),Ta)
    Yape_ode=ape[:,0]
    end = tm.time()
    print(end-start)
    global Yape_ode,Ta

#On defini n=10000 pour meilleure précision
Qape_ode(10000)

def Qcrit_ode(n):
    start = tm.time()
    L = 0.01
    C = 0.00001
    R = 63
    W0 = 1/sqrt(L*C)
    Q = (1/R)*sqrt(L/C)
    global W0,Q
    Tc=np.linspace(0,0.01,n+1)
    crit=spi.odeint(FuncSpi,(0.0001,0),Tc)
    Ycrit_ode=crit[:,0]
    end = tm.time()
    print(end-start)
    global Ycrit_ode,Tc

Qcrit_ode(10000)


def Qpseudo_ode(n):
    start = tm.time()
    L = 0.01
    C = 0.00001
    R = 3
    W0 = 1/sqrt(L*C)
    Q = (1/R)*sqrt(L/C)
    global W0,Q
    Tp=np.linspace(0,0.01,n+1)
    pseudo=spi.odeint(FuncSpi,(0.0001,0),Tp)
    Ypseudo_ode=pseudo[:,0]
    end = tm.time()
    print(end-start)
    global Ypseudo_ode,Tp

Qpseudo_ode(10000)

#Solution Exacte
def Exacte_ape(n):
    start = tm.time()
    L = 0.01
    C = 0.00001
    R = 160
    q0 = 0.0001
    W0 = 1/sqrt(L*C)
    Q = (1/R)*sqrt(L/C)
    r1=(-W0/(2*Q))+ (W0*sqrt((1/(2*Q)**2-1)))
    r2=(-W0/(2*Q)) - (W0*sqrt((1/(2*Q)**2-1)))
    A=-(q0*r2)/(r1-r2)
    B=(q0*r2)/(r1-r2)
    Yape_exa=[]
    T=np.linspace(0,0.01,n+1)
    for i in T:
        q = A*(exp(r1*i))+B*(exp(r2*i))
        Yape_exa.append(q)
    end = tm.time()
    print(end-start)
    global Yape_exa

Exacte_ape(10000)

def Exacte_crit(n):
    start = tm.time()
    L = 0.01
    C = 0.00001
    R = 63
    q0 = 0.0001
    W0 = 1/sqrt(L*C)
    Q = (1/R)*sqrt(L/C)
    A = q0
    B = W0*q0
    Ycrit_exa=[]
    T=np.linspace(0,0.01,n+1)
    for i in T:
        q = (A+B*i)*(exp(-W0*i))
        Ycrit_exa.append(q)
    end = tm.time()
    print(end-start)
    global Ycrit_exa

Exacte_crit(10000)

def Exacte_pseudo(n):
    start = tm.time()
    L = 0.01
    C = 0.00001
    R = 3
    q0 = 0.0001
    W0 = 1/sqrt(L*C)
    Q = (1/R)*sqrt(L/C)
    W = W0 * sqrt(1-(1/(2*Q))**2)
    Ypseudo_exa=[]
    T=np.linspace(0,0.01,n+1)
    for i in T:
        q = q0*(exp(-W0/(2*Q)*i))*(np.cos(W*i)+(1/((2*Q)*sqrt(1-((1/(2*Q))**2))))*np.sin(W*i))
        Ypseudo_exa.append(q)
    end = tm.time()
    print (end-start)
    global Ypseudo_exa

Exacte_pseudo(10000)


#Traçage Des Courbes
#Pour un n=10000 qui est assez élevé, les courbes se superposent quasiment

def Superpo_ape():
    p1=plt.plot(T1,Y1,"r",label="Euler")
    p2=plt.plot(Ta,Yape_ode,"b",label="Odeint")
    p3=plt.plot(Ta,Yape_exa,"g",label="Exacte")
    plt.xlabel("Temps")
    plt.ylabel("Charge")
    plt.legend()
    plt.grid(True)
    plt.title("Aperiodique Q<0.5")
    plt.savefig("Aperiodique Libre.png")

def Superpo_crit():
    p1=plt.plot(T1,Y2,"r",label="Euler")
    p2=plt.plot(Tc,Ycrit_ode,"b",label="Odeint")
    p3=plt.plot(Tc,Ycrit_exa,"g",label="Exacte")
    plt.xlabel("Temps")
    plt.legend()
    plt.ylabel("Charge")
    plt.grid(True)
    plt.title("Critique Q=0.5")
    plt.savefig("Critique Libre.png")

def Superpo_pseu():
    p1=plt.plot(T1,Y3,"r")
    p2=plt.plot(Tp,Ypseudo_ode,"b")
    p3=plt.plot(Tp,Ypseudo_exa,"g")
    plt.grid(True)
    plt.legend()
    plt.xlabel("Temps")
    plt.ylabel("Charge")
    plt.title("Pseudo Periodique Q>0.5")
    plt.savefig("Pseudo Periodique Libre.png")



#Resolution Regime Forcé
