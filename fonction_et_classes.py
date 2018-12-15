# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 19:12:53 2018

@author: alexh
"""
def liste_les_reacs_echange(model):
    """ retourne la liste des réactions d'échange"""
    import re
    liste_reac_ext=[]
    for r in model.reactions :
        if re.search("^EX",r.id):
            liste_reac_ext.append(r.id)
    return liste_reac_ext

def traduire_nom_reacs(listereacs,model):
    listenom=[]
    for r in listereacs:
        reac=model.reactions.get_by_id(r)
        for m in reac.metabolites.keys():
            nom=m.name[:20]
        listenom.append(nom)
    return(listenom)
    
def traduire_nom_reac(reac,model):
    if type(reac) == str:
        reac=model.reactions.get_by_id(reac)
    for m in reac.metabolites.keys():
        nom=m.name    
    return nom    

def traduire_reac_nom(nom,model):
    import re
    for m in model.metabolites:
        if nom==m.name:
            for r in m.reactions :
                if re.search("^EX",r.id):
                    reac=r.id
    return reac
    
        
        
def etudier_impact_reaction(reaction,borne,bacterie):
    from pandas import DataFrame
    rname=reaction.id
    T=[ 0 for i in range(10)]
    T1=[ 0 for i in range(10)]
    T2=[ 0 for i in range(10)]
    T3=[ 0 for i in range(10)]
    bounds=[(borne/10)*i*2 for i in range(5)]
    bounds+=[borne*2+borne*i*2 for i in range(5)]
    
    bounds[0]=borne*0.01
    
    
    for i in range(len(bounds)) :
        reaction.lower_bound = -bounds[i]
        bacterie.interactions[(rname,bounds[i])]=Interaction_milieu_bacterie(bacterie)
        
        bacterie.interactions[(rname,bounds[i])].deter_essentialite()
        bacterie.deter_les_croissances()
        T[i]=bacterie.croissance_global_max
        if len(bacterie.croissances_max.keys())>1:
            T1[i]=bacterie.croissances_max['Growth']
            T2[i]=bacterie.croissances_max['Growth_B1']
            T3[i]=bacterie.croissances_max['Growth_B2']
    #return(bounds,T,T1,T2,T3)
    dico={}
    for r in bacterie.interactions[(rname,borne*2)].echanges_essentiels_fva.keys():
        dico[r]=[]   
        
    for r in bacterie.reactions_echanges:
        for i in range(len(bounds)-1):
            if bacterie.interactions[(rname,bounds[i])].types_echanges[r] != bacterie.interactions[(rname,bounds[i+1])].types_echanges[r] :
                dico[r]=[]
                break
                
    for r in dico.keys():
        for i in range(len(bounds)):
            dico[r].append(bacterie.interactions[(rname,bounds[i])].types_echanges[r])
    df=DataFrame(dico,index=bounds)
    df=df.T

    return(bounds,T,T1,T2,T3, df)
    
def etudier_impact_reaction2(reaction,borne,interaction,nb_de_test=10):
    #from pandas import DataFrame
    T=[ 0 for i in range(nb_de_test)]
    #T1=[ 0 for i in range(10)]
    #T2=[ 0 for i in range(10)]
    #T3=[ 0 for i in range(10)]
    R1=[ 0 for i in range(nb_de_test)]
    R2=[ 0 for i in range(nb_de_test)]
    R3=[ 0 for i in range(nb_de_test)]
    nb=int(nb_de_test/2)
    bounds=[(borne-i*borne/10) for i in range(nb)]
    bounds.sort()
    bounds+=[borne*2+borne*i*2 for i in range(nb)]
    
    if bounds[0]==0:
        bounds[0]=borne*0.01
    
    
    for i in range(len(bounds)) :
        interaction.milieu_avant[reaction] = bounds[i]
        interaction.deter_les_croissances()
        #interaction.deter_essentialite()
        
        print(interaction.pourcentageglobaldumax)
        T[i]=interaction.pourcentageglobaldumax
        if len(interaction.croissances_max.keys())>1:
            #T1[i]=interaction.pourcentagesdumax['Growth']
            #T2[i]=interaction.pourcentagesdumax['Growth_B1']
            #T3[i]=interaction.pourcentagesdumax['Growth_B2']
            R1[i]=interaction.ratio_croissances['Growth']
            R2[i]=interaction.ratio_croissances['Growth_B1']
            R3[i]=interaction.ratio_croissances['Growth_B2']
        #print(interaction.results_fva)
    #dico={}
    
    #for r in interaction.types_echanges.keys():
        #if (2 in interaction.types_echanges[r]) or (-1 in interaction.types_echanges[r]) \
        #or (1 in interaction.types_echanges[r] and 1.5 in interaction.types_echanges[r]) or \
         #(1 in interaction.types_echanges[r] and -0.5 in interaction.types_echanges[r]) :
             #dico[r]=interaction.types_echanges[r]
    #listreac=list(dico.keys())
    #nomreacs=traduire_nom_reacs(listreac,interaction.bacterie.model)
    #df=DataFrame(dico,index=bounds)
    #df.columns=nomreacs
    #df=df.T

    return(bounds,T,R1,R2,R3)
    
def deter_zone_borne_interessant(bacterie,reaction):
    if type(reaction) == str:
        reaction=bacterie.model.reactions.get_by_id(reaction)
        
    reaction.lower_bound=-1000
    fba=bacterie.model.optimize()
    base=fba.f
    essai=base
    i=0
    
    while (base - essai) <5 and i<7:
        
        reaction.lower_bound=-(1000/(10**i))
        essai=bacterie.model.optimize().f
        print(i)
        borne=(1000/(10**i))
        i=i+1
        if i==7:
            print("pas de changement de croissance !")
    reaction.lower_bound=-1000
    print("quand " + reaction.id +" est borné à "+str(borne) + ", on a une diminution de croissance de "+str(base-essai))
    return borne

class Bacterie:
    """ objet représentant une bacterie et ses reactions, et échanges possibles"""
    def __init__(self,model):
        """ prend un paramètre un model cobrapy permettand de construire la bacterie :
            son .nom, son .model, son .nb_metabolites, son .nb_reactions, ses
            .reactions_echanges, son .objectif, son .milieu_de_base"""
        self.nom=str(model)
        self.model=model
        self.nb_metabolites=len(model.metabolites)
        self.nb_reactions=len(model.reactions)
        self.objective=model.objective
        self.reactions_echanges=liste_les_reacs_echange(model)
        self.milieu_de_base=model.medium
        self.interactions={}
        self.croissances_max={}
        self.croissance_global_max=0
        self.milieu_minimal={}
    def deter_les_croissances(self,dict_medium_change=None):
        import re
        fba=self.model.optimize()
        self.croissance_global_max=fba.f
        for r,f in fba.fluxes.iteritems():
            if re.search("Growth",r):
                self.croissances_max[r]=f
    def deter_milieu_minimal(self):
        from cobra.medium import minimal_medium
        self.milieu_minimal=minimal_medium(self.model, self.croissance_global_max)
        
class Interaction_milieu_bacterie:
    def __init__(self,bacterie):
        """ objet contenant, la .bacterie, le .milieu_avant que la bactérie interagisse avec, 
        les .echanges_essentiels de la bactérie dans ce milieu, le .nb_echanges_bloques,
        .nb_echanges_alternatifs, .nb_echanges_essentiels dans ce milieu, la .croissance
        de la bactérie dans ce milieu, et le .milieu_apres changé par l'interaction de la bacterie"""
        self.bacterie=bacterie
        self.milieu_avant=bacterie.model.medium
        self.milieu_apres={}
        self.echanges_essentiels_fba={}
        self.echanges_essentiels_fva={}
        self.nb_echanges_bloques=0
        self.nb_echanges_alternatifs=0
        self.nb_echanges_essentiels=0
        self.croissance=0
        self.types_echanges={}
    def deter_essentialite(self):
        """Apres fva, pour chaque réaction d'échanges avec l'extérieur, si
        les intervalles des flux optimaux à 99% sont strictement négatifs ou positifs,
        on sauvegarde la réaction d'échange et son intervalle dans le dico .echanges_essentiels
        , et on ajoute 1 à  .nb_echanges_essentiels_fva on compte également le nb
        d'échanges bloqués (dans nb_echanges_bloques) et alternatif dans nb_echanges_alternatifs"""
        #import cobra
        #self.model.medium=self.bacterie.milieu_avant
        fva=cobra.flux_analysis.flux_variability_analysis(self.bacterie.model, fraction_of_optimum=0.99)
        #self.diffgrowth=abs(abs(fva.minimum.get("Growth")-abs(fva.maximum.get('Growth'))))
        #self.diffgrowth1=
        #self.diffgrowth2=
        for r in self.bacterie.reactions_echanges:
            if abs(fva.minimum.get(r))<0.0000001 and abs(fva.maximum.get(r))<0.0000001:
                self.nb_echanges_bloques+=1
                self.types_echanges[r]=0
            elif fva.minimum.get(r) > 0.0000001 or fva.maximum.get(r)< -0.0000001 :
                bornes=(fva.minimum.get(r),fva.maximum.get(r))
                self.echanges_essentiels_fva[r]=bornes
                self.nb_echanges_essentiels+=1
                #delta=abs(abs(bornes[1]-abs(bornes[0]))-diffgrowth)
                #if delta<0.01 and fva.maximum.get(r)< -0.0000001:
                    #self.types_echanges[r]=3
                if fva.minimum.get(r) > 0.0000001:
                    self.types_echanges[r]=-1
                else:
                    self.types_echanges[r]=2
            else:
                self.nb_echanges_alternatifs+=1
                self.types_echanges[r]=1
                
    def echange_essentiel_par_fba(self):
        """on fait une FBA, et on crée un dico .echanges_essentiels_FBA, contenant 
        les echanges essentiels du dico .echanges_essentiels_FVA, mais cette fois ci avec les 
        valeurs de FBA correspondant"""
        fba=self.bacterie.model.optimize()
        self.croissance=fba.f
        for r in self.echanges_essentiels_fva.keys() :
            string="fba.fluxes."+r
            self.echanges_essentiels_fba[r]=round(float(eval(string)),3)
            
    def fba(self):
        fba=self.bacterie.model.optimize()
        self.croissance=fba.f
        
            
    def modif_milieu(self):
        """on remplit le dico .milieu_apres en faisant pour chaque echanges essentiels (valeur dans .milieu_avant) + (flux calculés par fba)  
        (si fba donne -200 pour l'échange EX_o2_e (consommation d'oxygène),
        dans le dico .milieu_apres, valeur de EX_o2_e sera par exemple de 1000-200 =800 """
        for r, flux in self.echanges_essentiels_fba.items():
            self.milieu_apres[r]=int(self.milieu_avant[r])+flux
            
    def __repr__(self):
        return "{} dans ce milieu :\n - {} échanges essentiels,\n - {} échanges bloqués,\n - une croissance de {} ".format(self.bacterie.nom,self.nb_echanges_essentiels,self.nb_echanges_bloques,self.croissance)
    def __str__(self):
        return repr(self)
    
def plot_intervalle_flux(dictionnaire_intervalle_flux):
    
    import pylab as pl
    from matplotlib import collections  as mc
    #from matplotlib.pyplot import figure
    lines = []
    legends = []
    y=1
    for r, (l,u) in dictionnaire_intervalle_flux.items():
        
        if u < 0 :
            lines.append([(l,y),(u,y)])
            legends.append((l-len(r)*30,y,r))
        else:
            lines.append([(l,y),(u,y)])
            legends.append((u+len(r)*30,y,r))
            
        
        y+=2
    
    lc = mc.LineCollection(lines, linewidths=2)
    fig, ax = pl.subplots()
    for legend in legends :
        pl.text(legend[0],legend[1],legend[2],horizontalalignment='center',fontsize=11,)
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(x=0.3)
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_position(('data',0))
    ax.yaxis.set_ticks_position('none')
    ax.set_yticklabels([])
    pl.show()
    
def fusion (bacteries):
    import cobra
    import re
    if type(bacteries) ==list:
        model = cobra.Model('Communautés')
        model= bacteries[0].model.copy()
        
        for i in range(1,len(bacteries)) :
            A=bacteries[i].model.copy()
            for m in A.metabolites:
                if m.compartment != 'C_e':
                    m.compartment=m.compartment+"_B"+str(i)
                    m.id=m.id+"_B"+str(i)
                    
            for r in A.reactions:
                if not re.search("EX_",r.id):
                    r.id = r.id +"_B"+str(i)
                    set_compartiments=r.compartments.copy()
                    r.compartments.clear()
                    for c in set_compartiments :
                        c=c+"_B"+str(i)
                        r.compartments.add(c)        
            model.add_reactions(A.reactions.copy())
            exec("model.reactions.Growth_B"+str(i)+".objective_coefficient=1")
    else :
        model=bacteries.model
    return model

            
    
class Communaute:
    """ objet représentant une communité bacteriene et ses reactions, et échanges possibles"""
    def __init__(self,bacteries):
        """ prend un paramètre un model cobrapy permettand de construire la bacterie :
            son .nom, son .model, son .nb_metabolites, son .nb_reactions, ses
            .reactions_echanges, son .objectif, son .milieu_de_base"""
        
        self.model=fusion(bacteries)
        self.nb_metabolites=len(self.model.metabolites)
        self.nb_reactions=len(self.model.reactions)
        self.objective=self.model.objective
        self.reactions_echanges=liste_les_reacs_echange(self.model)
        self.milieu_de_base=self.model.medium
        self.interactions={}
        self.croissances_max={}
        self.croissance_global_max=0
        self.milieu_minimal={}
        
    def deter_les_croissances(self):
        import re
        self.model.medium=self.milieu_de_base
        fba=self.model.optimize()
        self.croissance_global_max=fba.f
        for r,f in fba.fluxes.iteritems():
            if re.search("Growth",r):
                self.croissances_max[r]=f
    def deter_milieu_minimal(self):
        from cobra.medium import minimal_medium
        self.deter_les_croissances()
        self.milieu_minimal=dict(minimal_medium(self.model, self.croissance_global_max))
        
class Interaction_milieu_communaute:
    def __init__(self,bacterie):
        """ objet contenant, la .bacterie, le .milieu_avant que la bactérie interagisse avec, 
        les .echanges_essentiels de la bactérie dans ce milieu, le .nb_echanges_bloques,
        .nb_echanges_alternatifs, .nb_echanges_essentiels dans ce milieu, la .croissance
        de la bactérie dans ce milieu, et le .milieu_apres changé par l'interaction de la bacterie"""
        self.bacterie=bacterie
        self.milieu_avant=bacterie.model.medium
        self.milieu_apres={}
        self.changements_milieu={}
        self.echanges_essentiels_fba={}
        self.echanges_essentiels_fva={}
        self.echanges_alternatifs_fva={}
        self.nb_echanges_bloques=0
        self.nb_echanges_alternatifs=0
        self.nb_echanges_essentiels=0
        self.croissance=0
        self.types_echanges={}
        self.croissances_max={}
        self.croissance_global_max=0
        self.milieu_minimal={}
        self.pourcentageglobaldumax=0
        self.pourcentagesdumax={}
        self.results_fva={}
        self.ratio_croissances={}
        
    def deter_essentialite(self):
        """Apres fva, pour chaque réaction d'échanges avec l'extérieur, si
        les intervalles des flux optimaux à 99% sont strictement négatifs ou positifs,
        on sauvegarde la réaction d'échange et son intervalle dans le dico .echanges_essentiels
        , et on ajoute 1 à  .nb_echanges_essentiels_fva on compte également le nb
        d'échanges bloqués (dans nb_echanges_bloques) et alternatif dans nb_echanges_alternatifs"""
        import cobra
        
        self.bacterie.model.medium=self.milieu_avant
        fva=cobra.flux_analysis.flux_variability_analysis(self.bacterie.model, fraction_of_optimum=0.99)
        #self.diffgrowth=abs(abs(fva.minimum.get("Growth")-abs(fva.maximum.get('Growth'))))
        #self.diffgrowth1=
        #self.diffgrowth2=
        if len(list(self.results_fva.keys()))<1:
            for r in self.bacterie.reactions_echanges:
                self.results_fva[r]=[]
        for r in self.bacterie.reactions_echanges:
            bornes=(fva.minimum.get(r),fva.maximum.get(r))
            self.results_fva[r].append(bornes)
            if r not in self.types_echanges.keys() :
                self.types_echanges[r]=[]
            if abs(fva.minimum.get(r))<0.0000001 and abs(fva.maximum.get(r))<0.0000001:
                #self.nb_echanges_bloques+=1
                self.types_echanges[r].append(0)
            elif fva.minimum.get(r) > 0.0000001 or fva.maximum.get(r)< -0.0000001 :
                self.echanges_essentiels_fva[r]=bornes
                
                #self.nb_echanges_essentiels+=1
                #delta=abs(abs(bornes[1]-abs(bornes[0]))-diffgrowth)
                #if delta<0.01 and fva.maximum.get(r)< -0.0000001:
                    #self.types_echanges[r]=3
                if fva.minimum.get(r) > 0.0000001:
                    self.types_echanges[r].append(-1)
                else:
                    self.types_echanges[r].append(2)
            else:
                #self.nb_echanges_alternatifs+=1
                if fva.minimum.get(r) < -0.1 and fva.maximum.get(r)< 0.001:
                    self.types_echanges[r].append(1.5)
                elif fva.minimum.get(r) > -0.001 and fva.maximum.get(r)> 0.1:
                    self.types_echanges[r].append(-0.5)
                else :
                    self.types_echanges[r].append(1)
                
                
                
    def echange_essentiel_par_fba(self):
        """on fait une FBA, et on crée un dico .echanges_essentiels_FBA, contenant 
        les echanges essentiels du dico .echanges_essentiels_FVA, mais cette fois ci avec les 
        valeurs de FBA correspondant"""
        fba=self.bacterie.model.optimize()
        for r in self.echanges_essentiels_fva.keys() :
            string="fba.fluxes."+r
            self.echanges_essentiels_fba[r]=round(float(eval(string)),2)
            
    def fba(self):
        fba=self.bacterie.model.optimize()
        self.croissance=fba.f
        
            
    def modif_milieu(self):
        """on remplit le dico .milieu_apres en faisant pour chaque echanges essentiels (valeur dans .milieu_avant) + (flux calculés par fba)  
        (si fba donne -200 pour l'échange EX_o2_e (consommation d'oxygène),
        dans le dico .milieu_apres, valeur de EX_o2_e sera par exemple de 1000-200 =800 """
        #self.milieu_apres=dict(self.milieu_avant)
        if len(list(self.changements_milieu.keys()))<1:
            for r in self.milieu_avant.keys():
                self.changements_milieu[r]=[self.milieu_avant[r]]
                self.milieu_apres[r]=self.milieu_avant[r]
                
        for r in self.milieu_avant.keys() :
            if r in self.echanges_essentiels_fba.keys():
                #if self.milieu_avant[r]+self.echanges_essentiels_fba[r] > 0.1:
                self.milieu_apres[r]=self.milieu_avant[r]+self.echanges_essentiels_fba[r]
                self.changements_milieu[r].append(self.milieu_avant[r]+self.echanges_essentiels_fba[r])
                #else:
                    #self.milieu_apres[r]=0.1
                    #self.changements_milieu[r].append(0.1)
            else:
                self.changements_milieu[r].append(self.milieu_avant[r])
            
            
    def __repr__(self):
        return "{} dans ce milieu :\n - {} échanges essentiels,\n - {} échanges bloqués,\n - une croissance de {} ".format(self.bacterie.nom,self.nb_echanges_essentiels,self.nb_echanges_bloques,self.croissance)
    def __str__(self):
        return repr(self)    
    def deter_les_croissances(self):
        import re
        if self.bacterie.croissance_global_max==0:
            self.bacterie.deter_les_croissances()
        
        self.bacterie.model.medium=self.milieu_avant
        fba=self.bacterie.model.optimize()
        self.croissance_global_max=fba.f
        self.pourcentageglobaldumax=self.croissance_global_max/self.bacterie.croissance_global_max*100
        for r,f in fba.fluxes.iteritems():
            if re.search("Growth",r):
                self.croissances_max[r]=f
                self.ratio_croissances[r]=f/self.croissance_global_max
                self.pourcentagesdumax[r]=self.croissances_max[r]/self.bacterie.croissances_max[r]*100
    def deter_milieu_minimal(self):
        from cobra.medium import minimal_medium
        self.deter_les_croissances()
        self.milieu_minimal=dict(minimal_medium(self.bacterie.model, self.croissance_global_max))
    
    #c = np.array([(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)])
    
    