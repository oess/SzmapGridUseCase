from openeye.oechem import *
from openeye.oegrid import *
import sys


def ClearGridContent(grid):
    for i in range(grid.GetSize()):
        grid.SetValue(i, 0.0)

def GetGridDimension(grid):
    D = {}
    D['dim'] = dict(x=grid.GetXDim(), y=grid.GetYDim(), z=grid.GetZDim())
    D['spacing'] = grid.GetSpacing()
    D['min'] = dict(x=grid.GetXMin(), y=grid.GetYMin(), z=grid.GetZMin())
    D['max'] = dict(x=grid.GetXMax(), y=grid.GetYMax(), z=grid.GetZMax())
    D['mid'] = dict(x=grid.GetXMid(), y=grid.GetYMid(), z=grid.GetZMid())
    return D

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

    diff_grid = protein.GetData(DIFF_GRID_NAME)
    mask_grid = protein.GetData(MASK_GRID_NAME)
    meta_data = {}
    for line in protein.GetData(META_NAME).strip().split('\n'):
        it = line.split(':', 1)
        meta_data[it[0].strip()] = it[1].strip()

    return protein, ligand, diff_grid, mask_grid, meta_data

def GetMoleculeBox(mol):
    xmin, xmax = float('inf'), float('-inf')
    ymin, ymax = float('inf'), float('-inf')
    zmin, zmax = float('inf'), float('-inf')
    for atom in mol.GetAtoms():
        x, y, z = mol.GetCoords(atom)
        if x < xmin:
            xmin = x
        elif xmax < x:
            xmax = x
        if y < ymin:
            ymin = y
        elif ymax < y:
            ymax = y
        if z < zmin:
            zmin = z
        elif zmax < z:
            zmax = z
    xmid = (xmin+xmax)/2
    ymid = (ymin+ymax)/2
    zmid = (zmin+zmax)/2
    return dict(min=dict(x=xmin, y=ymin, z=zmin), max=dict(x=xmax, y=ymax, z=zmax), mid=dict(x=xmid, y=ymid, z=zmid))

def MaskGridToVirtualMoleculeGrid(maskgrid):
    mol = OEGraphMol()
    for i in range(maskgrid.GetSize()):
        value = maskgrid.GetValue(i)
        if value == 0:
            continue
        x, y, z = maskgrid.ElementToSpatialCoord(i)
        atom = mol.NewAtom(OEElemNo_C)
        mol.SetCoords(atom, (x,y,z))
    grid = OEScalarGrid()
    OEMakeMolecularGaussianGrid(grid, mol, maskgrid.GetSpacing())
    return grid

def FilteredGrid(grid, filter):
    filtered_grid = OEScalarGrid(grid)
    for i in range(grid.GetSize()):
        value = grid.GetValue(i)
        filtered_grid.SetValue(i, value if filter(value) else 0)
    return filtered_grid

def MergeGrid(opengrid, unfavgrid, weight):
    """
    Merge grid to size of opengrid
    both grid values are normalized to [0, 1]
    weight is multiplied to unfav value
    """
    o = GetGridDimension(opengrid)
    u = GetGridDimension(unfavgrid)
    assert \
        o['min']['x'] < u['min']['x'] and \
        o['min']['y'] < u['min']['y'] and \
        o['min']['z'] < u['min']['z'] and \
        u['max']['x'] < o['max']['x'] and \
        u['max']['y'] < o['max']['y'] and \
        u['max']['z'] < o['max']['z']

    ovalues = opengrid.GetValues()
    uvalues = unfavgrid.GetValues()
    omin, omax = min(ovalues), max(ovalues)
    umin, umax = min(uvalues), max(uvalues)

    merged_grid = OEScalarGrid(opengrid)
    ClearGridContent(merged_grid)
    for i in range(merged_grid.GetSize()):
        x, y, z = merged_grid.ElementToSpatialCoord(i)
        v1 = opengrid.GetValue(i)
        v1 = (v1 - omin)/(omax - omin)
        try:
            j = unfavgrid.SpatialCoordToElement(x, y, z)
            v2 = unfavgrid.GetValue(j)
            v2 = (v2 - umin)/(umax - umin)
            merged_grid.SetValue(i, v1 + weight*v2)
        except IndexError:
            pass

    return merged_grid


protein, ligand, diff_grid, mask_grid, meta_data = ReadDataFromSzmap(sys.argv[1])
open_grid = MaskGridToVirtualMoleculeGrid(mask_grid)
unfav_grid = FilteredGrid(diff_grid, lambda v: 0. < v)

print(GetGridDimension(open_grid))
print(GetGridDimension(unfav_grid))

merged_grid = MergeGrid(open_grid, unfav_grid, 20)

OEWriteGrid('vmol.grd', open_grid)
OEWriteGrid('unfav.grd', unfav_grid)
OEWriteGrid('merged.grd', merged_grid)
