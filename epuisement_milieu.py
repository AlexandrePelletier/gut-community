# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 20:27:17 2018

@author: alexh
"""

import seaborn as sns
import cobra 
from fonction_et_classes import *
import matplotlib.pyplot as plt
from pandas import DataFrame

#on charge les 3 models de bactéries dans une liste
BaA=Bacterie(cobra.io.read_sbml_model('Bacteroides_uniformis_ATCC_8492.xml'))
BaB=Bacterie(cobra.io.read_sbml_model('Bifidobacterium_longum_NCC2705.xml'))
BaC=Bacterie(cobra.io.read_sbml_model('Lactobacillus_iners_DSM_13335.xml'))
models=[BaA,BaB,BaC]
commu=Communaute(models)

#remettre medium a 1000:
commu.model.medium=commu.milieu_de_base
C3=Interaction_milieu_communaute(commu)


temps=[]
T=[ 0 for i in range(4)]
#T1=[ 0 for i in range(6)]
#T2=[ 0 for i in range(6)]
#T3=[ 0 for i in range(6)]
R1=[ 0 for i in range(4)]
R2=[ 0 for i in range(4)]
R3=[ 0 for i in range(4)]

for i in range(4):
    C3.deter_les_croissances()
    print("ok pour croissances", i)
    C3.deter_essentialite()
    print("ok pour essentialité", i)
    T[i]=C3.pourcentageglobaldumax
    #T1[i]=C3.pourcentagesdumax['Growth']
    #T2[i]=C3.pourcentagesdumax['Growth_B1']
    #T3[i]=C3.pourcentagesdumax['Growth_B2']
    R1[i]=C3.ratio_croissances['Growth']
    R2[i]=C3.ratio_croissances['Growth_B1']
    R3[i]=C3.ratio_croissances['Growth_B2']
    print("ok pour save croissance", i)
    C3.echange_essentiel_par_fba()
    print("ok pour fba", i)
    C3.modif_milieu()
    print("ok pour modif milieu", i)
    C3.milieu_avant=C3.milieu_apres
    print("ok pour milieu avant change", i)
    print(C3.pourcentageglobaldumax)
    print(C3.ratio_croissances)
    #print(C3.results_fva)
    print(C3.milieu_apres)
    temps.append(i)
    if C3.croissance_global_max<0.1:
        print("fini")
        break
nutriments_limitants={}
i=0
for r, f in C3.milieu_apres.items():
    if f<1:
        print(r, "( {} )".format(traduire_nom_reac(r,commu.model)), " : ", f )
        liste=[1000]
        for qte in C3.changements_milieu[r] :
            if qte < 0 :
                qte=0
            liste.append(qte)
        nutriments_limitants[r]=liste
        plt.plot(liste, label=traduire_nom_reac(r,commu.model))
        if i == 4 :
            plt.xlabel("temps")
            plt.ylabel("disponibilité dans le milieu",fontsize= 10)
            plt.legend(fontsize=8)
            plt.show()
        i=i+1
plt.xlabel("temps")
plt.ylabel("disponibilité dans le milieu",fontsize= 10)
plt.legend(fontsize=8)
plt.show()

nutriments_diminuant={}
i=0
for r, f in C3.milieu_apres.items():
    if f<800 and f>1:
        print(r, "( {} )".format(traduire_nom_reac(r,commu.model)), " : ", f )
        liste=[1000]
        for qte in C3.changements_milieu[r] :
            if qte < 0 :
                qte=0
            liste.append(qte)
        nutriments_diminuant[r]=liste
        plt.plot(liste, label=traduire_nom_reac(r,commu.model))
        if i == 4 :
            plt.xlabel("temps")
            plt.ylabel("disponibilité dans le milieu",fontsize= 10)
            plt.legend(fontsize=8)
            plt.show()
        i=i+1
plt.xlabel("temps")
plt.ylabel("disponibilité dans le milieu",fontsize= 10)
plt.legend(fontsize=8)
plt.show()

nutriments_augmentant_beaucoup={}
i=0
for r, f in C3.milieu_apres.items():
    if f>1500:
        print(r, "( {} )".format(traduire_nom_reac(r,commu.model)), " : ", f )
        liste=[1000]
        for qte in C3.changements_milieu[r] :
            if qte < 0 :
                qte=0
            liste.append(qte)
        nutriments_limitants[r]=liste
        plt.plot(liste, label=traduire_nom_reac(r,commu.model))
        if i == 4 :
            plt.xlabel("temps")
            plt.ylabel("disponibilité dans le milieu",fontsize= 10)
            plt.legend(fontsize=8)
            plt.show()
        i=i+1
plt.xlabel("temps")
plt.ylabel("disponibilité dans le milieu",fontsize= 10)
plt.legend(fontsize=8)
plt.show()


nutriments_augmentant={}
i=0
for r, f in C3.milieu_apres.items():
    if f<1500 and f >1010:
        print(r, "( {} )".format(traduire_nom_reac(r,commu.model)), " : ", f )
        liste=[1000]
        for qte in C3.changements_milieu[r] :
            if qte < 0 :
                qte=0
            liste.append(qte)
        nutriments_limitants[r]=liste
        plt.plot(liste, label=traduire_nom_reac(r,commu.model))
        if i == 4 :
            plt.xlabel("temps")
            plt.ylabel("disponibilité dans le milieu",fontsize= 10)
            plt.legend(fontsize=8)
            plt.show()
        i=i+1
plt.xlabel("temps")
plt.ylabel("disponibilité dans le milieu",fontsize= 10)
plt.legend(fontsize=8)
plt.show()

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
        for i in range(3):
            if C3.types_echanges[r][i]!=C3.types_echanges[r][i+1] :
                #print(r)
                dico[r]=C3.types_echanges[r]
                break
            
listreac=list(dico.keys())
nomreacs=traduire_nom_reacs(listreac,C3.bacterie.model)
df=DataFrame(dico,index=temps)
df.columns=nomreacs
df=df.T
#changer nom
dico_intervalle_flux1={}
dico_intervalle_flux2={}
dico_intervalle_flux3={}
dico_intervalle_flux4={}
for r, liste in C3.results_fva.items():
    if r in dico.keys() and (liste[1][0] >0.1 or liste[1][1] <-0.1 or liste[2][0]>0.1 or liste[2][1]<-0.1 or liste[3][0]>0.1 or liste[3][1] < -0.1) :
        nom=traduire_nom_reac(r,commu.model)
        if len(nom)<20 :
            dico_intervalle_flux1[nom]=liste[0]
            dico_intervalle_flux2[nom]=liste[1]
            dico_intervalle_flux3[nom]=liste[2]
            dico_intervalle_flux4[nom]=liste[3]
print(len(dico_intervalle_flux2))
plot_intervalle_flux(dico_intervalle_flux1)
plot_intervalle_flux(dico_intervalle_flux2)
plot_intervalle_flux(dico_intervalle_flux3)
plot_intervalle_flux(dico_intervalle_flux4)


#%matplotlib inline #spé au notebook

plt.plot(temps,T,'r--',label='global')
plt.plot(temps,[r1*100 for r1 in R1],label=BaA.nom[:-10])
plt.plot(temps,[r1*100 for r1 in R2],label=BaB.nom[:-8])
plt.plot(temps,[r1*100 for r1 in R3],label=BaC.nom[:-10])
plt.xlabel("temps")
plt.ylabel("croissance global et proportion des espèces  ",fontsize= 10)
plt.legend()
plt.title("impact epuisement du milieu au cours du temps")
plt.show()
"""
plt.plot(temps,T,label='global')
plt.plot(temps,T1,label=BaA.nom)
plt.plot(temps,T2,label=BaB.nom)
plt.plot(temps,T3,label=BaC.nom)
plt.xlabel("temps")
plt.ylabel("proportion de l'espèce dans la communauté",fontsize= 8)
plt.legend()
plt.title("impact epuisement du milieu au cours du temps")
plt.show()
"""


plt.figure(figsize=(10,15))
sns.heatmap(df,linecolor='black',linewidths=.1)
plt.xlabel("temps")
#cmap="YlGnBu"
plt.show()