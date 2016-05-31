#Micro Projet Informatique - Etude RLC Régime Libre et Forcé
####*Quentin Bergé - Benjamin Deprat*

--
## Introduction

Nous étudierons un circuit RLC d'abord dans le cas d'un régime libre puis d'un régime forcé. Le but de cette étude sera de montrer l'influence du pas d'intégration sur la résolution numérique d'une équation différentielle du second ordre par la méthode d'Euler et la méthode Odeint de SciPy.

Pour ceci, nous solverons l'équation différentielle de la charge de ce circuit par ces 2 méthodes en comparant nos résultats à la solution exacte dans les cas d'un régime apériodique, critique et pseudo-périodique.

Nous comparerons alors nos résultats avec différents pas d'intégration à la solution exacte connue et comparerons les vitesses d'exécution de ces méthode.

Dans un second temps nous étudierons le cas du circuit soumis à une excitation sinusoidale afin de conforter nos premières constatations.



##Régime Libre

###Constantes
On choisis les constantes de notre circuit en cohérence avec les composants réels.

* L = 10mH 
* C = 10μF 
*  W0 = 1/((L*C)^0,5) = 3162 rad/s
*  Alors T = 0.002s, le regime transitoire est de 5T environ alors tf = 0.01s
*  On fera varier R = [160,63,3] Ω afin d'avoir Q = [0.2,0.5,10.5]
*  q0 = 0.1 mC
*  Wf=(2*np.pi)/0.01 On choisi f=100Hz pour l'excitation sinusoidale
*  Amp=10 avec une amplitude de 10V

### Architecture du Programme

Importation des différentes bibliothèques qui vont nous être utiles et déclaration des constantes inhérentes à notre circuit.
	
##### Méthode Euler
    	
   On repète les mêmes instructions pour les régimes Critique, et Pseudos-Periodique, seul le Q change.
   On divise cependant en plusieurs programme pour avoir plus de facilité à tracer le tout à la fin de chaque exécution de ces fonctions pour un pas d'intégration (n) différent.
   
  
##### Méthode Odeint

	def FuncSpi(q,T):
    	return (q[1],-(W0**2)*q[0]-(W0/Q)*q[1])
	# définition de la fonction pour odeint selon discretisation du problème
	# Eq differentielle du 2nd Ordre, alors retour d'un vecteur
	
#####Méthode Exacte
On utilise les solutions connues de l'équation, puis on stocke les valeurs dans une liste.


### Courbes
Voici les differentes courbes obtenues en Régime Libre en Fonction de Pas d'intégration (n).

![Alt text](http://img11.hostingpics.net/pics/894507ComparaisonReponsesLibres.png)

### Vitesse d'Exécution

![Alt Test](http://img11.hostingpics.net/pics/708676VitesseCalculfRegimeetPasLibre.png)

###Analyse des Résultats

On peut alors voir que plus le pas d'intégration est grand, plus les résultats des méthodes d'Euler et d'Odeint sont en corrélation avec la solution exacte.
La méthode la plus influencé par le pas reste la méthode d'Euler.
Pour un pas de 100, la solution est complètement en désaccord avec les 2 autres.
Alors que pour un pas de 10000 par exemple, les courbes se superposent. 

Néanmoins, la vitesse d'éxécution de cette méthode augmente aussi sensiblement avec le pas d'intégration. On peut observer que le régime pseudo-periodique nécéssitera lui le plus important temps de calcul pour quelque méthode que ce soit. 

Enfin la méthode Odeint est très peu influencée par le pas d'intégration par rapport à la méthode d'Euler.



##Regime Sinusoidal Forcé

###Architecture du Programme

#####Constantes
On reprend les constantes précédentes mais enrichies de Wf, pulsation du générateur sinusoidal que l'on défini avec f=100 Hz et une amplitude de 10V.
#####Méthode Euler

On réalise la méthode d'Euler avec la nouvelle discrétisation de l'équation differentielle.
    	
##### Méthode Odeint
 
 Définition de la fonction qui retournera les valeurs pour la fonction odeint.
 	
 	def FuncSpiForce(q,T):
    	return (q[1],-(W0**2)*q[0]-(W0/Q)*q[1]+ Amp*np.cos(Wf*T))
		#nouvelle discrétisation du probleme selon régime forcé
	
Puis on intègre cette fonction dans la méthode Odeint.
 

###Courbes
![Alt Test](http://img11.hostingpics.net/pics/542712ComparaisonReponsesForc.jpg)

###Vitesse d'Exécution

![Alt Test](http://img11.hostingpics.net/pics/302370VitesseDeCalculFctRgimeetPasForc.png)

###Analyse des Résultats
On peut à nouveau observer que l'augmentation du pas d'intégration (n) entraine une augmentation de la précision.
Dans le cas du faible pas, la méthode Odeint à l'air de mieux se comporter. Néanmoins les deux méthodes se superposent graphiquement dans le cas d'un pas grand (10000).

La méthode d'Euler est encore une fois plus influencée par le pas d'intégration. Pour un pas faible le temps de calcul sera faible et vice versa. Tandis que la méthode Odeint reste environ stable selon le pas d'intégration. Elle est donc peu inflencée et son temps d'execution est sensiblement inferieur à la méthode d'Euler lorsque ce pas est grand.


### Difficultées Rencontrées

Dans le cas du circuit RLC en régime sinusoidal forcé, il est difficile de conclure sur la véracité des méthodes en l'absence de la solution exacte en méthode temporelle. 

##Conclusion

Avec cette illustration de cas, nous avons pu analyser l'influence du pas d'intégration sur les méthode de determination d'une equation différentielle du 2nd ordre dans ce cas ci.
Il parait les conclusions suivantes:

* La méthode d'Euler se comporte bien lorsque le pas d'intégration est grand voire très grand. Elle est bien plus facile à mettre en place. Néanmoins, elle nécéssite un temps d'exécution bien plus important que la méthode Odeint.
* La méthode Odeint est efficace des lors que le pas d'intégration est suffisamment grand (ici pour n=100). Le temps de calcul reste stable lors de la variation du pas. Elle est quant à elle un peu plus compliquée à mettre en place.

Il s'agira alors de choisir le meilleur compromis entre facilité à mettre en place, précision des résultats et temps de calcul.

