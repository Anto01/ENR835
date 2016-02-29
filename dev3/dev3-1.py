# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 14:04:12 2016

@author: Antoine
"""

from solar_mod import *
import numpy as np
#from collections import *

# donnees du probleme
#
T_pc = 52.0         # Temperature de la plaque (Celsius)
T_pk = T_pc+273     # Temperature de la plaque (Kelvin)
T_ac = 18           # Temperature ambiante (celsius)
T_ak = T_ac+273     # Temperature ambiante (kelvin)
h = 6               # coefficient de convection extérieur
Lair  = 0.030       # epaisseur d'air
alpha_n = 0.95      # Coef d'absortion solaire
rhog = 9.8          
KL = 0.0125         # Propriété du vitrage
n2=1.526            # indice de réfraction de la vitre
n1=1                # indice de réfraction de l'air
Emitt = 0.17        # Emissivité de la plaque
emig = 0.88         # Emissivitée du verre
N = 1               # Nombre de vitrage
beta = 0            # Angle à plat


# Importation des données du problème devoir 2 question 1
I = np.load('I.npy')
the = np.load('thez.npy')
Rb = np.load('Rb.npy')
Id = np.load('Id.npy')
Ib = np.load('Ib.npy')

# Pour la 12e tranche
tranche = 12
Rb = Rb[tranche]
Id = Id[tranche]
Ib = Ib[tranche]
I = I[tranche]
the = the[tranche]


the_g = angle_reflechi(the)
the_d = angle_diffus(the)
alpn = alp_alpn(the) # coeficient d,absorbeur faible longueur d'onde

 
#tap_b = alp_alpn(the_b)
#tap_g = alp_alpn(the_g)
#tap_d = alp_alpn(the_d)

tau_al_b = Calcul_tau_al(the,alpn,KL,n2,n1,N)
tau_al_g = Calcul_tau_al(the,alpn,KL,n2,n1,N)
tau_al_d = Calcul_tau_al(the,alpn,KL,n2,n1,N)


Ibt = Ib*Rb*tau_al_b
Idt = Id*tau_al_d*(1+cosd(beta))/2.0
Irt = I*rhog*tau_al_g*(1-cosd(beta))/2.0
S =  Ibt+Idt+Irt # ABSORBED SOLAR RADIATION (W/m2)

# Calcul des pertes par unités de surface
# Modèle empirique
pertes_empirique = U_Klein(T_pc,T_ac,beta,h,Emitt,emig,N) #W/m2 ◦C

q = pertes_empirique*(T_pc-T_ac)


# Calcul du rendement du capteur
Rendement = q/I

