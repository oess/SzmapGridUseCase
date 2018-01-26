from openeye.oechem import *
from openeye.oegrid import *
from algebraicgrid import AlgebraicGrid
import sys


spec = """
!PARAMETER -in
  !TYPE string
  !REQUIRED true
  !BRIEF input file which has szmap grids
!END
!PARAMETER -out
  !TYPE string
  !REQUIRED false
  !DEFAULT shape.grd
  !BRIEF output file of shape_grid
!END
!PARAMETER -weight
  !TYPE float
  !REQUIRED false
  !DEFAULT 2.0
  !BRIEF weight of unfavored water region
!END
"""


def ReadDataFromSzmap(iname):
    DIFF_GRID_NAME = 'neut_diff_apo_free_energy_grid'
    MASK_GRID_NAME = 'apo_mask_grid'
    META_NAME = 'szmap-meta-data'

    protein = OEGraphMol()
    ligand = OEGraphMol()

    ifs = oemolistream(iname)
    assert ifs.IsValid()
    OEReadMolecule(ifs, protein)
    OEReadMolecule(ifs, ligand)
    ifs.close()

    diff_grid = AlgebraicGrid(protein.GetData(DIFF_GRID_NAME))
    mask_grid = AlgebraicGrid(protein.GetData(MASK_GRID_NAME))
    assert AlgebraicGrid.IsCompatible(diff_grid, mask_grid)

    meta_data = {}
    for line in protein.GetData(META_NAME).strip().split('\n'):
        it = line.split(':', 1)
        meta_data[it[0].strip()] = it[1].strip()

    return protein, ligand, diff_grid, mask_grid, meta_data


def main(itf):
    """
    neutraldiff_grid: neut_diff_apoo_free_energy_grid
    mask_grid: apo_mask_grid (maskgrid)
    unfavor_grid: positive part of dgrid and normalized to [0, 1]
    openspace_grid: non-zero part of maskgrid multiplied by (1 - scaled distance)
    shape_grid: openspace_grid + weight*unfavor_grid
    """
    iname = itf.GetString('-in')
    weight = itf.GetFloat('-weight')
    oname = itf.GetString('-out')

    print('weight={}'.format(weight))

    protein, ligand, neutraldiff_grid, mask_grid, meta_data = ReadDataFromSzmap(iname)

    mask_grid.Ternarize()
    unfavor_grid = neutraldiff_grid.FilteredGrid(filter=lambda v: v > 0.0).Rescale()
    center = unfavor_grid.GetCentroid()
    print(center)

    unmask_grid = ((1 - mask_grid.DistanceGrid(center).Rescale())*mask_grid).Rescale()
    unfavor_grid = (unfavor_grid**(1.0/10)).Rescale()

    shape_grid = (unmask_grid + weight*unfavor_grid).Rescale()
    print(shape_grid.GetStats())

    OEWriteGrid(oname, shape_grid)

if __name__ == '__main__':
    itf = OEInterface(spec, sys.argv)
    main(itf)
