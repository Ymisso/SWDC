#! /usr/bin/env python3

"""
Le module "Hydrophobe" contient 2 fonctions qui permettent
de sélectionner les résidus hydrophobes et ensuite de calculer
la distance qui existe entre leur carbone Apha.
"""

def Select_Hydro(ListAtom, RESD_Hydro):
    """
    Sélectionne la liste des atomes succeptibles de créer des 
    Liasons Hydrophobes

    Argument(s): (1)Liste des objets "atome"
                 (2)Liste des résidus hydrophobes
    """
    List_Hydro = []
    for atome in ListAtom:
        if atome.Resid in RESD_Hydro and atome.Nom=="CA":
            List_Hydro.append(atome)

    return List_Hydro



def Calc_Int_Hydro(List_Hydro):
    """
    Calcule les distances entre les carbones Alpha des résidus hydrophobes.
    Renvoie une Liste de Liste contenant le couple d'atomes qui interagit.

    Argument(s): Liste des objets "atome" hydrophobes
    """
    Hydrop = []
    for at1 in List_Hydro:
        for at2 in List_Hydro:
            if round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2) < 5.3 \
                and round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 \
                 + (float(at1.z)-float(at2.z))**2)**0.5, 2) > 0:
                Hydrop.append([at1, at2])
    return Hydrop

