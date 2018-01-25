# def GetMoleculeBox(mol):
#     xmin, xmax = float('inf'), float('-inf')
#     ymin, ymax = float('inf'), float('-inf')
#     zmin, zmax = float('inf'), float('-inf')
#     for atom in mol.GetAtoms():
#         x, y, z = mol.GetCoords(atom)
#         if x < xmin:
#             xmin = x
#         elif xmax < x:
#             xmax = x
#         if y < ymin:
#             ymin = y
#         elif ymax < y:
#             ymax = y
#         if z < zmin:
#             zmin = z
#         elif zmax < z:
#             zmax = z
#     xmid = (xmin+xmax)/2
#     ymid = (ymin+ymax)/2
#     zmid = (zmin+zmax)/2
#     return dict(min=dict(x=xmin, y=ymin, z=zmin), max=dict(x=xmax, y=ymax, z=zmax), mid=dict(x=xmid, y=ymid, z=zmid))
#
# def MaskGridToVirtualMoleculeGrid(maskgrid):
#     mol = OEGraphMol()
#     for i in range(maskgrid.GetSize()):
#         value = maskgrid.GetValue(i)
#         if value == 0:
#             continue
#         x, y, z = maskgrid.ElementToSpatialCoord(i)
#         atom = mol.NewAtom(OEElemNo_C)
#         mol.SetCoords(atom, (x,y,z))
#     grid = OEScalarGrid()
#     OEMakeMolecularGaussianGrid(grid, mol, maskgrid.GetSpacing())
#     return grid
#
# def FilteredGrid(grid, filter):
#     filtered_grid = OEScalarGrid(grid)
#     for i in range(grid.GetSize()):
#         value = grid.GetValue(i)
#         filtered_grid.SetValue(i, value if filter(value) else 0)
#     return filtered_grid
#
# def MergeGrid(opengrid, unfavgrid, weight):
#     """
#     Merge grid to size of opengrid
#     both grid values are normalized to [0, 1]
#     weight is multiplied to unfav value
#     """
#     o = GetGridDimension(opengrid)
#     u = GetGridDimension(unfavgrid)
#     assert \
#         o['min']['x'] < u['min']['x'] and \
#         o['min']['y'] < u['min']['y'] and \
#         o['min']['z'] < u['min']['z'] and \
#         u['max']['x'] < o['max']['x'] and \
#         u['max']['y'] < o['max']['y'] and \
#         u['max']['z'] < o['max']['z']
#
#     ovalues = opengrid.GetValues()
#     uvalues = unfavgrid.GetValues()
#     omin, omax = min(ovalues), max(ovalues)
#     umin, umax = min(uvalues), max(uvalues)
#
#     merged_grid = OEScalarGrid(opengrid)
#     ClearGridContent(merged_grid)
#     for i in range(merged_grid.GetSize()):
#         x, y, z = merged_grid.ElementToSpatialCoord(i)
#         v1 = opengrid.GetValue(i)
#         v1 = (v1 - omin)/(omax - omin)
#         try:
#             j = unfavgrid.SpatialCoordToElement(x, y, z)
#             v2 = unfavgrid.GetValue(j)
#             v2 = (v2 - umin)/(umax - umin)
#             merged_grid.SetValue(i, v1 + weight*v2)
#         except IndexError:
#             pass
#
#     return merged_grid
