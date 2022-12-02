import smi2cm
from rdkit import Chem
from rdkit.Chem import EState
from rdkit.Chem.EState import Fingerprinter
import numpy as np
import pandas as pd


#generate CM descriptors for 16 rings
mols = [np.linalg.eig(smi2cm.smi2cm(smi, 3, True))[0] for smi in  ['c1cc(c(cc1)C)O',
                                                                'c1c(c(c(cc1)C)O)C',
                                                                'c1c(c(c(cc1C(C)(C)C)C)O)C(C)(C)C',
                                                                'c1cc(c(cc1OC)C)O',
                                                                'c1c(c(c(cc1N(=O)=O)C)O)C',
                                                                'c1c(c(c(cc1Cl)C)O)Cl',
                                                                'c1cc(c(cc1Cl)C)O',
                                                                'c1c(c(c(cc1C)C)O)[C@@H]1[C@@H]2C[C@@H]3C[C@H]1C[C@H](C2)C3',
                                                                'c1c(c(c(cc1)C)O)C(C)C',
                                                                'c1c(c(c(cc1)C)O)c1ccccc1',
                                                                'c1c(c(c(cc1C)C)O)C(C)(C)C',
                                                                'c1c(c(c(cc1C(C)(C)c1ccccc1)C)O)C(C)(C)c1ccccc1',
                                                                'c1c(c(c(cc1Cl)C)O)C(C)(C)C',
                                                                'c1c(c(c(cc1C)C)O)C',
                                                                'c1c(c(c(cc1)C)O)[Si](C)(C)C(C)(C)C',
                                                                'c1c(c(c(cc1Br)C)O)Br']]


mols_new = []
for mol in mols:

    mn = mol.tolist()+[0]*(80-len(mol))
    mnn = sorted(mn, reverse=True)
    mols_new.append(mnn)

data = pd.DataFrame(mols_new)
data.to_csv("CM_rings.csv")


#generate CM descriptors for 36 linkers
mols2 = [np.linalg.eig(smi2cm.smi2cm(smi, 3, True))[0] for smi in ['C(C/N=C/C)/N=C/C',
                                                                'C(CC/N=C/C)/N=C/C',
                                                                'CC(C)(C/N=C/C)C/N=C/C',
                                                                'CCC(CC)(C/N=C/C)C/N=C/C',
                                                                '[C@H]1([C@H](CCCC1)/N=C\C)/N=C/C',
                                                                'c1(c(cccc1)/N=C/C)/N=C\C',
                                                                'c1c(c2c(cc1)cccc2/N=C/C)/N=C/C',
                                                                'c1(c(cccc1)/N=C\C)c1c(cccc1)/N=C\C',
                                                                'c1cc(c(cc1)/N=C/C)CCc1ccccc1/N=C/C',
                                                                'c1(cccc(c1c1c(cccc1C)/N=C/C)C)/N=C/C',
                                                                'c1c(cccc1)[C@H]([C@@H](c1ccccc1)/N=C/C)/N=C\C',
                                                                'c1cc(ccc1)CC(Cc1ccccc1)(C/N=C/C)C/N=C/C',
                                                                'C(CN(CC)C)N(CC)C',
                                                                'C(CCN(C)CC)N(C)CC',
                                                                'CC(C)(CN(C)CC)CN(C)CC',
                                                                'CCC(CC)(CN(CC)C)CN(C)CC',
                                                                '[C@H]1([C@H](CCCC1)N(C)CC)N(C)CC',
                                                                'c1(c(cccc1)N(C)CC)N(C)CC',
                                                                'c1c(c2c(cc1)cccc2N(CC)C)N(CC)C',
                                                                'c1(c(cccc1)N(CC)C)c1c(cccc1)N(C)CC',
                                                                'c1cc(c(cc1)N(CC)C)CCc1ccccc1N(CC)C',
                                                                'c1(cccc(c1c1c(cccc1C)N(CC)C)C)N(CC)C',
                                                                'c1c(cccc1)[C@H]([C@@H](c1ccccc1)N(C)CC)N(C)CC',
                                                                'c1cc(ccc1)CC(Cc1ccccc1)(CN(CC)C)CN(CC)C',
                                                                'C(CN(CC)Cc1ccccc1)N(CC)Cc1ccccc1',
                                                                'C(CCN(Cc1ccccc1)CC)N(Cc1ccccc1)CC',
                                                                'CC(C)(CN(Cc1ccccc1)CC)CN(Cc1ccccc1)CC',
                                                                'CCC(CC)(CN(CC)Cc1ccccc1)CN(Cc1ccccc1)CC',
                                                                '[C@H]1([C@H](CCCC1)N(Cc1ccccc1)CC)N(Cc1ccccc1)CC',
                                                                'c1(c(cccc1)N(Cc1ccccc1)CC)N(Cc1ccccc1)CC',
                                                                'c1c(c2c(cc1)cccc2N(CC)Cc1ccccc1)N(CC)Cc1ccccc1',
                                                                'c1(c(cccc1)N(CC)Cc1ccccc1)c1c(cccc1)N(Cc1ccccc1)CC',
                                                                'c1cc(c(cc1)N(CC)Cc1ccccc1)CCc1ccccc1N(CC)Cc1ccccc1',
                                                                'c1(cccc(c1c1c(cccc1C)N(CC)Cc1ccccc1)C)N(CC)Cc1ccccc1',
                                                                'c1c(cccc1)[C@H]([C@@H](c1ccccc1)N(Cc1ccccc1)CC)N(Cc1ccccc1)CC',
                                                                'c1cc(ccc1)CC(Cc1ccccc1)(CN(CC)Cc1ccccc1)CN(CC)Cc1ccccc1']]


mols2_new = []
for mol in mols2:

    mn = mol.tolist()+[0]*(80-len(mol))
    mnn = sorted(mn, reverse=True)
    mols2_new.append(mnn)

data = pd.DataFrame(mols2_new)
data.to_csv("CM_linkers.csv")
