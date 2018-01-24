from openeye.oechem import *
import sys

spec = """
!PARAMETER -in
  !TYPE string
  !REQUIRED true
  !BRIEF Input filename
!END
!PARAMETER -p
  !TYPE string
  !REQUIRED true
  !BRIEF Output protein filename
!END
!PARAMETER -l
  !TYPE string
  !REQUIRED true
  !BRIEF Output ligand filename
!END
"""

def CountPeptideBonds(mol):
    unique = True
    smarts = 'NC(=O)'
    ss = OESubSearch(smarts)
    return len(list(ss.Match(mol, unique)))

itf = OEInterface(spec, sys.argv)

iname = itf.GetString('-in')
ifs = oemolistream(iname)
imol = OEGraphMol()
assert ifs.IsValid()
OEReadMolecule(ifs, imol)
ifs.close()

sopt = OESplitMolComplexOptions()
sopt.ResetFilters(0)
sopt.SetWaterFilter(OEMolComplexFilterFactory(OEMolComplexFilterCategory_Nothing))

count = 0
mols = []
for mol in OEGetMolComplexComponents(imol, sopt):
    count += 1
    npep = CountPeptideBonds(mol)
    mol.SetData('npep', npep)
    mols.append(OEGraphMol(mol))

mols.sort(key=lambda mol: mol.GetData('npep'), reverse=True)

protein, ligand = mols[0], mols[1]

protein.SetTitle('Protein')
ligand.SetTitle('Ligand')

OEWriteMolecule(oemolostream(itf.GetString('-p')), protein)
OEWriteMolecule(oemolostream(itf.GetString('-l')), ligand)
