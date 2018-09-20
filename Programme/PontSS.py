#! /usr/bin/env python3

"""
Le module "PontSS" contient des fonctions qui permettent de
créer une liste d'objets "atome" qui ne contient que les atomes
de Soufre des Cystéines, puis de de calculer la distance entre ces 
derniers
"""

def Select_Disulfure(ListAtom):
    """
    Sélectionne la liste des atomes succeptibles de faire des ponts S-S.

    Argument(s): Liste des objets atomes
    """
    List_SS = []
    for atome in ListAtom:
        if atome.Nom=="SG":
            List_SS.append(atome)

    return List_SS


def Calc_Pts_SS(List_SS):
    """
    Calcule les distances entre les atomes de soufre des Cystéines.
    Renvoie une Liste de Liste contenant le couple d'atomes qui interagit,
    et leur distance.

    Argument(s): Liste des atomes de soufre de Cystéines.
    """
    PontSS = []
    for at1 in List_SS:
        for at2 in List_SS:
            if round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2) < 2.2 \
                and round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 \
                 + (float(at1.z)-float(at2.z))**2)**0.5, 2) > 0:

                dist = round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2)
                PontSS.append([at1, at2, dist])
    return PontSS


