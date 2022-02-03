from utils import*
import gemmi

crambin = gemmi.read_pdb('data/1ejg.pdb')
residues = crambin[0][0]
atoms =[]
for res in residues:
    atoms += [[a.element.atomic_number] + a.pos.tolist() for a in res]
