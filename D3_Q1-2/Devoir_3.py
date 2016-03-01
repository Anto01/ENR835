
# coding: utf-8

# Devoir #3 ENR382
# ========

# ## Modules à importer

# In[26]:

# Python
import numpy as np
import scipy as sp
import scipy.constants as const
get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
import pandas as pd
import pylab
import sys

import plotly
from plotly.offline import plot
#init_notebook_mode()

# Externe
import solar_mod as sm
import properties_mod as pm


# 
# ## Question 1
# 
# ### Importation et apperçu des données du devoir #2

# In[27]:

data_d2 = pd.read_csv('data_devoir2.csv')
print(data_d2)
#data_d2.head()


# ### Données du problème

# In[28]:

T_pc = 52.0         # Temperature de la plaque (Celsius)
T_pk = T_pc+273     # Temperature de la plaque (Kelvin)
T_ac = 18           # Temperature ambiante (celsius)
T_ak = T_ac+273     # Temperature ambiante (kelvin)
h_w = 6.0             # coefficient de convection extérieur
Lair  = 0.030       # epaisseur d'air
alpha_n = 0.95      # Coef d'absortion solaire
g = 9.8  
rhog = 0.4
KL = 0.0125         # Propriété du vitrage
n2=1.526            # indice de réfraction de la vitre
n1=1                # indice de réfraction de l'air
eps_p = 0.17        # Emissivité de la plaque
eps_c = 0.88         # Emissivitée du verre
N = 1               # Nombre de vitrage
# Notes: Nous posons l'hyposthèse que les angles beta et gama du capteur sont les mêmes que ceux du devoir #2
beta = 60         # Angle à plat
gam = 0             # plein sud
phi  = 45.0 +30.0/60.0  # latitude  Montreal (45 deg 30 min nord)


# ### Sélection des données de la 12e tranche

# In[29]:

tranche = 11
I = data_d2.I[tranche]
Ibn = data_d2.Ibn[tranche]
n =sm.jour_mois_jour_annee(12,'juin')
thez = data_d2.thez[tranche]
delt  = sm.decl_solaire(n)
omen = data_d2.omen[tranche]
Ib = Ibn * sm.cosd(thez)
Id = I-Ib
Rb = sm.calcul_Rb(phi,n,omen,beta,gam)

It,Ib,Id,Ir = sm.modele_isotropique(I,Ib,Id,beta,Rb,rhog)

print('omen =',omen)
print('Rb = ',Rb)
print('It = ',It)
print('Ib = ',Ib)
print('Id = ',Id)
print('Ir =',Ir)
print('thez = ',thez)


# ### Calculs du problème
# #### Calcul des angles réfléchi et diffus

# In[30]:

the_g = sm.angle_reflechi(beta)
the_d = sm.angle_diffus(beta)
#the = abs(thez-beta)

the = sm.normale_solaire(delt,phi,omen,beta,gam)

print('the = ',the)
print('the_g = ',the_g)
print('the_d = ',the_d)


# <img  src="thetag.PNG"/>

# <img  src="thetad.PNG"/>

# #### Calcul du coéfficient absorbeur pour faible longueur d'ondes

# In[31]:

tau_al_b = sm.Calcul_tau_al(the,alpha_n,KL,n2,n1,N)
tau_al_g = sm.Calcul_tau_al(the_g,alpha_n,KL,n2,n1,N)
tau_al_d = sm.Calcul_tau_al(the_d,alpha_n,KL,n2,n1,N)

print('tau_al_b =', tau_al_b)
print('tau_al_g =', tau_al_g)
print('tau_al_d =', tau_al_d)


# #### Calcul de la radiation totale transmise

# In[32]:

Ibt = Ib*tau_al_b
Idt = Id*tau_al_d
Irt = Ir*tau_al_g
S =  Ibt+Idt+Irt

print('Radiation directe transmise =',Ibt)
print('Radiation diffuse transmise=',Idt)
print('Radiation réfléchie transmise =',Irt)
print('Radiation totale transmise (S) =',S)


# <img  src="S.PNG"/>

# #### Calcul des pertes par model itératif
# 
# Pour une plaque simple, le coéfficient de perte vers le haut est exprimé grâce à cette formule.
# <img  src="ut.PNG"/>

# Premièrement, il faut trouver chacun des termes du dénominateur. Les formules suivantes seront utilisées.
# Pour le coéfficient de radiation de la plaque interne jusqu'à la surface externe,
# <img  src="hrpc.PNG"/>
# Pour ce qui est de la radiation de la surface externe vers l'air ambiante,
# <img  src="hrca.PNG"/>
# La convection entre les deux surfaces ($ h_{c,p-c} $) sera déterminé grâce à la fonction suivante, 
#   <img  src="hcpc.PNG"/>
# Pour enfin trouver le coéfficient de température de la face externe du capteur.
# <img  src="Tc.PNG"/>

# In[33]:


T_cc = 40 # Température initiale (hypothèse)
T_old = 1
converg = []
x_it = 0
x_conv = []

while (abs(T_old-T_cc) > 0.0001):
    
    converg.append(T_cc)
    T_old = T_cc

    T_ck = T_cc + 273.0
    T_avg = (T_ck+T_pk)/2

    hr_pc = (const.sigma * (T_pk**2+T_ck**2) * (T_pk+T_ck))/((1/eps_p)+(1/eps_c)-1)
    
    hr_ca = eps_c*const.sigma*(T_ck**2+T_ak**2)*(T_ck+T_ak)

    nu = pm.air_prop('nu',T_avg)
    al = pm.air_prop('al',T_avg)
    k = pm.air_prop('k',T_avg)
    Ra = const.g*(1/T_avg)*abs(T_ck-T_pk)*Lair**3/(nu*al)

    f1 = max(0,1.0-1708.0/(Ra*sm.cosd(beta)))
    f2 = 1.0-1708*(sm.sind(1.8*beta))**1.6/(Ra*sm.cosd(beta))
    f3 = max(0,(Ra*sm.cosd(beta)/5830)**(1.0/3.0)-1.0)
    Nu = 1.0 + 1.44*f1*f2+f3 
    hc_pc = Nu * k/Lair

    T_p = T_avg-273.0
    
    
    Ut = ((1/(hc_pc + hr_pc))+(1/(h_w+hr_ca)))**-1
    
    T_cc = T_p - ((Ut*(T_p-T_ac))/(hc_pc+hr_pc))
    x_conv.append(x_it)
    x_it = x_it + 1

q = Ut*(T_pc-T_ac)
total_iteratif = S-q


#plotly.offline.iplot({
#"data": [{
#    "x": x_conv,
#    "y": converg
#}],
#"layout": {
#    "title": "Convergence de la température externe"
#}
#})



print('Le coéfficient de perte vers le haut =',Ut, 'W/m² °C')
print('Les pertes vers le haut =', q, 'W/m²')
print('Ce qui est capté',total_iteratif,'W/m²' )
print('La température de la surface extérieure =',T_cc,'°C')


# #### Calcul des pertes par modèle empirique Klein
# 
# En plus de calculer avec la méthode itérative, la méthode Klein est utilisé afin de comparer les résultats.

# <img  src="uklein.PNG"/>

# In[34]:

pertes_empirique_coef = sm.U_Klein(T_pc,T_ac,beta,h_w,eps_p,eps_c,N)
pertes_empirique = pertes_empirique_coef*(T_pc-T_ac)
total_empirique = S-pertes_empirique


print('Coefficient de perte vers le haut Klein',pertes_empirique_coef,'W/m² °C')
print('Pertes vers le haut =',pertes_empirique, 'W/m²')
print('Ce qui est capté =',total_empirique, 'W/m²')


# #### Calcul du rendement du capteur

# In[35]:

Rendement_empirique = total_empirique/It
Rendement_iteratif = total_iteratif/It
print('Le rendement du capteur est de',Rendement_empirique,'pour la méthode empirique')
print('Le rendement du capteur est de',Rendement_iteratif,'pour la méthode itérative')


# ### Résultats
# 
# Notes: Nous posons l'hyposthèse que les angles beta et gama du capteur sont les mêmes que ceux du devoir #2
# 

# In[36]:

print('a)',S,'J/m²')
print('b)',q, 'W/m²')
print('c)',Rendement_iteratif)


# ## Question 2

# ### Données du problème

# In[37]:

H= 0.8
N=8
Y=2
Ac=H*Y # Surface des capteurs
W=0.1 # Distance entre les tubes (m)
D = 0.007 # m
Cb=100.0 # W/mK
P=2*(H+Y) # Perimetre du capteur
deltaa= 0.001 #epaisseur de la plaque (m)
ka=400.0  # W/m2K
mpt=0.016 # Debit total Kg/s
hf=1100.0 #W/m2K
Cp=4180.0   #J/KgK
Rpjoint = 1/Cb    #hypoth?se Conductivite joint est infinie
Ti=18.0
UL=Ut


# #### Calcul du rendement d'ailette F

# <img  src="m.PNG"/>

# In[38]:

m = np.sqrt(UL/(ka*deltaa))


# <img  src="f.PNG"/>

# In[39]:

F = np.tanh((m*(W-D)/2))/(m*(W-D)/2)


# #### Calcul du rendement d'absorbeur F'
# <img  src="Fp.PNG"/>

# In[40]:

Fp = (1.0/UL)/(W*(1.0/(UL*(D+(W-D)*F))+Rpjoint+(1.0/(hf*const.pi*D))))
print('Fp =',Fp)


# #### Calcul du rendement d'absorbeur F''
# <img  src="fpp.PNG"/>

# In[41]:

Fpp1 = (mpt*Cp)/(Ac*UL*Fp)
Fpp2 =  1-np.exp(-(Ac*UL*Fp)/(mpt*Cp))
Fpp = Fpp1 * Fpp2
print(Fpp)


# #### Calcul du facteur FR

# <img  src="fr.PNG"/>

# In[42]:

Fr = Fpp*Fp
print(Fr)


# #### Calcul de la température à l'entrée du capteur

# <img  src="tfi.PNG"/>
# En isolant Tfi de la formule, on obtient,

# In[43]:

Qu = Ac*Fr*(S-UL*(Ti-T_ac))
Tfi = T_pc-((Ac)/(Fr*UL))*(1-Fpp)


# ### Résultats

# In[44]:

print('a)',Fr)
print('b)',Tfi,'°C')

