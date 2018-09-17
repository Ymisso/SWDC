#! /usr/bin/env python3

#creer une classe des atomes avec les infos tirées du pdb

class Atom:
    """
    Ceci est la classe Atome
    """
    
    def __init__(self, Num=0, Nom="CC", Resid="AAA", Chain="A", NumRes=0, 
                    x=0.0, y=0.0, z=0.0):
        """
        definition des objts atome
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
	Récupère les informations de tous les atomes et renvoie une liste
	"""
    ListAtom = []
    for line in pdbfile:
        infoATOM = line.split()
        if infoATOM[0] == "ATOM": 
            atom=Atom(infoATOM[1], infoATOM[2], infoATOM[3], infoATOM[4], infoATOM[5], \
                        infoATOM[6], infoATOM[7], infoATOM[8])
            ListAtom.append(atom)

    return ListAtom

