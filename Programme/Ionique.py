#! /usr/bin/env python3

"""
Ce module contient des fonctions qui permettent de sélectionner 
les Carbone alpha des résidu chargés, de calculer leur distances
et de déterminer s'il y a ou non intéraction
"""

def Select_RESD_AromtqP(ListAtom, RESD_IonPos):
    """
    Sélectionne la liste des atomes associés à des residus chargés
    positivement

    Argument(s): (1) Liste des objets "atome"
                 (2) Liste des résidu positifs
    """
    List_ResIonP = []
    for atome in ListAtom:
        if atome.Resid in RESD_IonPos and atome.Nom=="CA":
            List_ResIonP.append(atome)

    return List_ResIonP

def Select_RESD_AromtqN(ListAtom, RESD_IonNeg):
    """
    Sélectionne la liste des atomes associés à des residus chargés
    négativement.

    Argument(s): (1) Liste des objets "atome"
                 (2) Liste des résidus négatifs
    """
    List_ResIonN = []
    for atome in ListAtom:
        if atome.Resid in RESD_IonNeg and atome.Nom=="CA":
            List_ResIonN.append(atome)

    return List_ResIonN

def Calcul_IntIon(List_ResIonP, List_ResIonN):
    """
    Calcule la distance entre un atome donneur et un atome accepteur.
    Renvoie une liste des atomes qui interagissent et leur distance.

    Argument(s): (1) Liste des CA des résidus chargés +
                 (2) Liste des CA des résidus chargés -
    """
    Ionq=[]
    for at1 in List_ResIonP:
        for at2 in List_ResIonN:
            
            if round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2)==0.00:
                continue

            elif at1.Nom == at2.Nom or at1.Num == at2.Num:
                continue

            elif round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2) < 6.00:
                Dist = round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2)
                Ionq.append([at1, at2, Dist])
    return Ionq

