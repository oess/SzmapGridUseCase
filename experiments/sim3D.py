from openeye.oechem import *
from openeye.oegrid import *
from openeye.oeshape import *
import sys

shape_grid_name = sys.argv[1]
test_mols_fname = sys.argv[2]
out_fname = sys.argv[3]

shape_grid = OEScalarGrid()

OEReadGrid(shape_grid_name, shape_grid)

mols = [OEMol(mol) for mol in oemolistream(test_mols_fname).GetOEMols()]

prep = OEOverlapPrep()
overlay = OEOverlay()
overlay.SetupRef(shape_grid)

omols = []
for mol in mols:
    # prep.Prep(mol)
    scoreiter = OEBestOverlayScoreIter()
    OESortOverlayScores(scoreiter, overlay.Overlay(mol), OEHighestRefTversky())
    for score in scoreiter:
        title = mol.GetTitle()
        print('{} {:.3f} {}'.format(score.fitconfidx, score.GetRefTversky(), title))
        omol = OEGraphMol(mol.GetConf(OEHasConfIdx(score.fitconfidx)))
        omol.SetData('reftversky', score.GetRefTversky())
        score.Transform(omol)
        omol.SetTitle('{} ({:.3f})'.format(title, score.GetRefTversky()))
        omols.append(omol)
        break

omols.sort(key=lambda m: m.GetData('reftversky'), reverse=True)

ofs = oemolostream(out_fname)
for mol in omols:
    OEWriteMolecule(ofs, mol)
ofs.close()
