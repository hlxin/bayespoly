# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 13:38:00 2020

@author: Alexander van Teijlingen
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def smi2cm(smi, dimensions = 2, Hs = True, return_xyz = False):
    mol = Chem.MolFromSmiles(smi)
    if Hs:
        mol = Chem.AddHs(mol)
    else:
        mol = Chem.RemoveHs(mol)
    if int(dimensions) == 2:
        AllChem.Compute2DCoords(mol)
    elif int(dimensions) == 3:
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
    else:
        print("Invalid input for parameter: dimensions\nPlease only use either 2 or 3")
    xyz = Chem.MolToXYZBlock(mol)
    xyzmatrix = xyz_parse(xyz)
    if return_xyz:
        return xyzmatrix
    else:
        return gen_coulombmatrix(xyzmatrix)

def xyz_parse(xyz):
    nAtoms = int(xyz.split("\n")[0])
    xyzmatrix = np.ndarray((nAtoms, 4), dtype="float")
    ANs = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109}
    i = 0
    for line in xyz.split("\n"):
        line = line.split()
        if len(line) == 4:
            xyzmatrix[i][0] = float(ANs[line[0]])
            xyzmatrix[i][1] = float(line[1])
            xyzmatrix[i][2] = float(line[2])
            xyzmatrix[i][3] = float(line[3])
            i+=1
    return xyzmatrix

def gen_coulombmatrix(xyzmatrix):
    """
    From: Rupp, M.; Tkatchenko, A.; Müller, K. R.; Von Lilienfeld, O. A. Fast and Accurate Modeling of Molecular Atomization Energies with Machine Learning. Phys. Rev. Lett. 2012, 108 (5), 1–5. https://doi.org/10.1103/PhysRevLett.108.058301.
    """
    nAtoms = int(xyzmatrix.shape[0])

    cij = np.zeros((nAtoms, nAtoms))
    
    for i in range(nAtoms):
        for j in range(nAtoms):
            if i == j:
                cij[i][j] = 0.5 * xyzmatrix[i][0] ** 2.4  # Diagonal term described by Potential energy of isolated atom
            else:
                dist = np.linalg.norm(np.array(xyzmatrix[i][1:]) - np.array(xyzmatrix[j][1:]))

                cij[i][j] = xyzmatrix[i][0] * xyzmatrix[j][0] / dist  # Pair-wise repulsion
            #cij[i][j] = float(cij[i][j])
    return cij


def visualise(smi, showHs = True):
    from mayavi import mlab
    xyzmatrix = smi2cm(smi, 3, showHs, True)
    cij = smi2cm(smi, 3, showHs)
    x = xyzmatrix[:,1]
    y = xyzmatrix[:,2]
    z = xyzmatrix[:,3]
    s = [sum(cij[d]) for d in range(cij.shape[0])]
    pts = mlab.points3d(x, y, z, s)
    mlab.show()
    
