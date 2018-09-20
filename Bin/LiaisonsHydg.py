#! /usr/bin/env python3

"""
Ce module contient des fonctions qui sélectionnent
les donneurs et les accepteurs potentiels de LH,
puis calculent leur distance.
"""

def Select_HbondD(ListAtom, ATOM_Dnn):
    """
    Sélectionne la liste des atomes donneurs de LH

    Argument(s): (1) Liste des objets "atome"
                 (2) Liste des atomes donneurs
    """
    List_HbondD = []
    for atome in ListAtom:
        if atome.Nom in ATOM_Dnn:
            if atome.Resid=="PRO" and atome.Nom=="N":
                continue
            List_HbondD.append(atome)

    return List_HbondD

def Select_HbondA(ListAtom, ATOM_Acc):
    """
    Sélectionne la liste des atomes accepteurs de LH

    Argument(s): (1) Liste des objets "atome"
                 (2) Liste des atomes accepteurs
    """
    List_HbondA = []
    for atome in ListAtom:
        if atome.Nom in ATOM_Acc:
            if atome.Resid=="PRO" and atome.Nom=="N":
                continue
            List_HbondA.append(atome)

    return List_HbondA


def Calcul_LH(List_HbondD, List_HbondA):
    """
    Calcule la distance entre un atome donneur et un atome accepteur.
    Renvoie une liste des atomes qui interagissent et leur distance.

    Argument(s): (1) Liste des potentiels donneurs
                 (2) Liste des potentiels accepteurs
    """
    LH=[]
    for at1 in List_HbondD:
        for at2 in List_HbondA:
            if at1.Resid=="PRO" and at1.Nom=="N":
                continue
            
            elif round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2)==0.00:
                continue

            elif at1.Nom == at2.Nom or at1.Num == at2.Num:
                continue

            elif round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2) < 3.00 \
                and (at1.Nom!="SG" or at2.Nom!="SG"):
                Dist = round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2)
                LH.append([at1, at2, Dist])

            elif round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2) < 4.00 \
                and (at1.Nom=="SG" or at2.Nom=="SG"):
                Dist = round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2)
                LH.append([at1, at2, Dist])
    return LH


