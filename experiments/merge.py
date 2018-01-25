from openeye.oechem import *
from openeye.oegrid import *
import sys, math
from mathgrid import MathGrid


def ListNormalizeInPlace(v, reverse=False):
    vmin, vmax = min(v), max(v)
    for i in range(len(v)):
        v[i] = (v[i] - vmin)/(vmax - vmin)
        if reverse:
            v[i] = 1.0 - v[i]
        if v[i] < 0.0:
            v[i] = 0.0
        elif 1.0 < v[i]:
            v[i] = 1.0

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

    diff_grid = MathGrid(protein.GetData(DIFF_GRID_NAME))
    mask_grid = MathGrid(protein.GetData(MASK_GRID_NAME))
    meta_data = {}
    for line in protein.GetData(META_NAME).strip().split('\n'):
        it = line.split(':', 1)
        meta_data[it[0].strip()] = it[1].strip()

    return protein, ligand, diff_grid, mask_grid, meta_data


def GetUnmaskDistanceGrid(maskgrid, center):
    grid = MathGrid(maskgrid, initialize=True)
    N = grid.size
    cx, cy, cz = center
    buf = [0]*N
    for i in range(N):
        x, y, z = grid.ElementToSpatialCoord(i)
        buf[i] = math.sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz))
    ListNormalizeInPlace(buf, reverse=True)
    for i in range(N):
        grid[i] = buf[i]*maskgrid[i]
    return grid


protein, ligand, diff_grid, mask_grid, meta_data = ReadDataFromSzmap(sys.argv[1])
unfav_grid = diff_grid.GetFilteredGrid(filter=lambda v: 0.0 < v).GetNormalizedGrid()
center = unfav_grid.GetCentroid()
undist_grid = GetUnmaskDistanceGrid(mask_grid, center)
merged_grid = undist_grid + unfav_grid

OEWriteGrid('z.grd', merged_grid)
