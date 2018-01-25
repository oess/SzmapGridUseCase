def GetMoleculeBoundingBox(mol):
    xs, ys, zs = [], [], []
    for atom in mol.GetAtoms():
        x, y, z = mol.GetCoords(atom)
        xs.append(x)
        ys.append(y)
        zs.append(z)
    xmid = (min(xs)+max(xs))/2
    ymid = (min(ys)+max(ys))/2
    zmid = (min(zs)+max(zs))/2
    return dict(min=dict(x=min(xs), y=min(ys), z=min(zs)),
                max=dict(x=max(xs), y=max(ys), z=max(zs)),
                mid=dict(x=xmid, y=ymid, z=zmid))
