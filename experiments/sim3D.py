from openeye.oechem import *
from openeye.oegrid import *
from openeye.oeshape import *

grid = OEScalarGrid()

OEReadGrid('out.grd', grid)

mols = [OEMol(mol) for mol in oemolistream('testmc.oeb').GetOEMols()]

prep = OEOverlapPrep()
overlay = OEOverlay()
overlay.SetupRef(grid)

omols = []
for mol in mols:
    # prep.Prep(mol)
    scoreiter = OEBestOverlayScoreIter()
    OESortOverlayScores(scoreiter, overlay.Overlay(mol), OEHighestRefTversky())
    for score in scoreiter:
        print(score.fitconfidx, score.GetRefTversky())
        omol = OEGraphMol(mol.GetConf(OEHasConfIdx(score.fitconfidx)))
        omol.SetData('reftversky', score.GetRefTversky())
        score.Transform(omol)
        omol.SetTitle('{} ({:.3f})'.format(omol.GetTitle(), score.GetRefTversky()))
        omols.append(omol)
        break

omols.sort(key=lambda m: m.GetData('reftversky'), reverse=True)

ofs = oemolostream('out.oeb')
for mol in omols:
    OEWriteMolecule(ofs, mol)
ofs.close()
