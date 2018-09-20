#! /usr/bin/env python3

#creer une classe des atomes avec les infos tirées du pdb

"""
Le module "Atome" contient une classe qui permet de créer des objets "atome"
et une fonction qui parse un pdb et qui rempli les objets "atome" avec les
informations des atomes contenus dans le pdb
"""

class Atom:
    """
    La classe Atome permet de creer des objets atomes
    """
    
    def __init__(self, Num=0, Nom="CC", Resid="AAA", Chain="A", NumRes=0, 
                    x=0.0, y=0.0, z=0.0):
        """
        definition des objets atome
        """
        self.Num=Num
        self.Nom=Nom
        self.Resid=Resid
        self.Chain=Chain
        self.NumRes=NumRes
        self.x=x
        self.y=y
        self.z=z


def Parse_and_Constructor(pdbfile):
    """
	Récupère les informations de tous les atomes dans le pdbfile
    et renvoie une liste d'objets "atome" qui contiennent ces informations.

    Argument(s):  Fichier pdb
	"""
    ListAtom = []
    for line in pdbfile:
        infoATOM = line.split()
        if infoATOM[0] == "ATOM": 
            atom=Atom(infoATOM[1], infoATOM[2], infoATOM[3], infoATOM[4], infoATOM[5], \
                        infoATOM[6], infoATOM[7], infoATOM[8])
            ListAtom.append(atom)

    return ListAtom

