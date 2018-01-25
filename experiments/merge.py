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
  !DEFAULT out.grd
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
    diff_grid: neut_diff_apoo_free_energy_grid
    mask_grid: apo_mask_grid
    unfav_grid: positive part of diff_grid and normalized to [0, 1]
    undist_grid: non-zero part of mask_grid multiplied by (1 - normalized distance)
    shape_grid: undist_grid + weight*unfav_grid
    """
    iname = itf.GetString('-in')
    weight = itf.GetFloat('-weight')
    oname = itf.GetString('-out')

    protein, ligand, diff_grid, mask_grid, meta_data = ReadDataFromSzmap(iname)
    unfav_grid = diff_grid.GetFilteredGrid(filter=lambda v: 0.0 < v).GetNormalizedGrid()
    center = unfav_grid.GetCentroid()
    undist_grid = (1 - mask_grid.GetDistanceGrid(center).GetNormalizedGrid())*mask_grid.GetNormalizedGrid()
    shape_grid = (undist_grid + weight*unfav_grid).GetNormalizedGrid()
    print(shape_grid.GetMinMax())
    OEWriteGrid(oname, shape_grid)


if __name__ == '__main__':
    itf = OEInterface(spec, sys.argv)
    main(itf)
