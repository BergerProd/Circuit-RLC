#Micro Projet Informatique - Etude RLC Régime Libre et Forcé
####*Quentin Bergé - Benjamin Deprat*

--
## Introduction



##Régime Libre

###Constantes
On choisis les constantes de notre circuit en cohérence avec les composants réels.

* L = 10mH 
* C = 10μF 
*  W0 = 1/((L*C)^0,5) = 3162 rad/s
*  Alors T = 0.002s, le regime transitoire est de 5T environ alors tf = 0.01s
*  On fera varier R = [160,63,3] Ω afin d'avoir Q = [0.2,0.5,10.5]
*  q0 = 0.1 mC

### Architecture du Programme

Importation des différentes bibliothèques qui vont nous être utiles et déclaration des constantes inhérentes à notre circuit

	import matplotlib.pyplot as plt
	import scipy.integrate as spi
	from math import sqrt,exp
	import numpy as np
	import time as tm
	L = 0.01
	C = 0.00001
	R=[160,63,3]
	W0 = 1/sqrt(L*C)
	q0 = 0.0001	
	
##### Méthode Euler

	def Euler_Aperiodique(n):
   		start = tm.time() 
    	h=(0.01-0)/n
   	 	t,y,z=0,0.0001,0
   	 	Teuler=[0]
   		Yeuler_ape=[0.0001]
    	Zeuler_ape=[0]
    	R0 = R[0]
    	Q = (1/R0)*sqrt(L/C) #ici Q=0,2 environ, cas d'un regime Aperiodique
    	for i in range (n):
        	t,y,z=t+h,y+h*z,z+h*(-((W0/Q)*z)-(W0**2)*y) #selon discrétisation du problème
        	Teuler.append(t)
        	Yeuler_ape.append(y)
        	Zeuler_ape.append(z)
    	end = tm.time()
    	print ("Euler Aperiodique =",end - start) #mesure du temps
    	global Teuler,Yeuler_ape,Zeuler_ape    #globalisation des listes pour les tracer après
    	
   On repète les mêmes instructions pour les régimes Critique, et Pseudos-Periodique, seul le Q change.
   On divise cependant en plusieurs programme pour avoir plus de facilité à tracer le tout à la fin de chaque execution de ces fonctions pour un pas d'intégration (n) différent.
   
  
##### Méthode Odeint

	def FuncSpi(q,T):
    	return (q[1],-(W0**2)*q[0]-(W0/Q)*q[1])
	# définition de la fonction pour odeint selon discretisation du problème
	# Eq differentielle du 2nd Ordre, alors retour d'un vecteur
	
	def Ode_Aperiodique(n):
   	 	start = tm.time()
    	R0 = R[0]
    	Q = (1/R0)*sqrt(L/C)
    	global Q #globalisation du Q pour l'utiliser dans la fonction ci dessus
    	Tode_ape=np.linspace(0,0.01,n+1) #création liste de temps
    	ape=spi.odeint(FuncSpi,(0.0001,0),Tode_ape) #recourt à la méthode odeint de scipy
    	Yode_ape=ape[:,0] #stockage des solutions qui seront sous formes d'un array
    	end = tm.time()
    	print("Odeint Critique = ",end-start) #information sur le temps d'execution
    	global Yode_ape,Tode_ape #globalisation pour traçage
    
#####Méthode Exacte

	def Exacte_ape(n):
   	 	start = tm.time()
    	R0 = R[0]
    	Q = (1/R0)*sqrt(L/C)
    	r1=(-W0/(2*Q))+ (W0*sqrt((1/(2*Q)**2-1))) 
    	r2=(-W0/(2*Q)) - (W0*sqrt((1/(2*Q)**2-1))) 
    	A=-(q0*r2)/(r1-r2)
    	B=(q0*r2)/(r1-r2)
   	 	Yexa_ape=[]
    	Texa_ape=np.linspace(0,0.01,n+1) #création liste de temps
    	for i in Texa_ape: #parcours de la liste de temps
        	q = A*(exp(r1*i))+B*(exp(r2*i)) #forme de la solution 
        	Yexa_ape.append(q)
    	end = tm.time()
    	print("Exacte Aperiodique = ",end-start)
    	global Yexa_ape,Texa_ape
    	
    




### Courbes
Voici les differentes courbes obtenues en Régime Libre en Fonction de Pas d'intégration (n).

![Alt text](http://img11.hostingpics.net/pics/894507ComparaisonReponsesLibres.png)

### Vitesse d'Exécution

![Alt Test](http://img11.hostingpics.net/pics/843642VitesseCalculfRegimeetPasLibre.png)

###Analyse des Résultats

On peut alors voir que plus le pas d'intégration est grand, plus les résultats des méthode d'Euler et d'Odeint sont en corrélation avec la solution exacte.
La méthode la plus influencé par le pas reste la méthode d'Euler.
Pour un pas de 100, la solution est complètement en désaccord avec les 2 autres.
Alors que pour un pas de 10000 par exemple, les courbes se superposent. 

Néanmoins, la vitesse d'éxécution de la fonction augmente elle aussi sensiblement avec le pas d'intégration. On peut observer que le régime pseudo-periodique nécéssitera lui le plus important temps de calcul pour quelque méthode que ce soit. 

Enfin là méthode Odeint est très peu influencé par le pas d'intégration par rapport aux méthodes d'Euler et de la Solution Exacte qui elle est la la plus influencé



##Regime Sinusoidal Forcé

###Courbes

###Vitesse d'Exécution

###Analyse des Résultats

### Difficultées Rencontrées

##Conclusion