# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 15:44:48 2018

@author: alexh
"""

#import re
import seaborn as sns
import cobra 
from fonction_et_classes import *
import matplotlib.pyplot as plt
import re
#on charge les 3 models de bactéries dans une liste
BaA=Bacterie(cobra.io.read_sbml_model('Bacteroides_uniformis_ATCC_8492.xml'))
BaB=Bacterie(cobra.io.read_sbml_model('Bifidobacterium_longum_NCC2705.xml'))
BaC=Bacterie(cobra.io.read_sbml_model('Lactobacillus_iners_DSM_13335.xml'))
models=[BaA,BaB,BaC]
commu=Communaute(models)



C3=Interaction_milieu_communaute(commu)

commu.model.medium=commu.milieu_de_base
C3.milieu_avant=commu.milieu_de_base
#nutriments_indisponibles=["EX_g1p_e","EX_man1p_e"]
#for r in nutriments_indisponibles:
#    C3.milieu_avant[r]=0


#C3.deter_milieu_minimal()
#print(C3.milieu_minimal)
reaction=input("entrer reaction a étudier :")
if not re.search("EX",reaction) :
    reaction=traduire_reac_nom(reaction,commu.model)

borne=deter_zone_borne_interessant(commu,reaction)
bounds,T,R1,R2,R3=etudier_impact_reaction2(reaction,borne,C3,nb_de_test=10)

"""
dico={}
    
for r in C3.types_echanges.keys():
    if (2 in C3.types_echanges[r]) :
         dico[r]=C3.types_echanges[r]
         #print (r)
for r in C3.types_echanges.keys():
    if (-1 in C3.types_echanges[r]) :
        dico[r]=C3.types_echanges[r]
        #print (r)
for r in C3.types_echanges.keys():
     if r not in dico.keys():
        for i in range(len(bounds)-1):
            if C3.types_echanges[r][i]!=C3.types_echanges[r][i+1] :
                #print(r)
                dico[r]=C3.types_echanges[r]
                break
listreac=list(dico.keys())
nomreacs=traduire_nom_reacs(listreac,C3.bacterie.model)
df=DataFrame(dico,index=bounds)
df.columns=nomreacs
df=df.T
#changer nom
"""

#%matplotlib inline #spé au notebook

plt.plot(bounds,T,'r--',label='global')
plt.plot(bounds,[r1*100 for r1 in R1],label=BaA.nom[:-10])
plt.plot(bounds,[r1*100 for r1 in R2],label=BaB.nom[:-8])
plt.plot(bounds,[r1*100 for r1 in R3],label=BaC.nom[:-10])
plt.xlabel("quantité de {} disponible".format(traduire_nom_reac(reaction,commu.model)))
plt.ylabel("croissance global et proportion des espèces",fontsize= 10)
plt.legend()
plt.title("impact de l'accès en {}".format(traduire_nom_reac(reaction,commu.model)))
plt.show()

"""
plt.figure(figsize=(10,12))
sns.heatmap(df,cmap="YlGnBu",linecolor='black',linewidths=.5)
plt.show()"""