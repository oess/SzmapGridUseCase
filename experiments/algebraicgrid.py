from openeye import oechem, oegrid
import numpy as np

class AlgebraicGrid(oegrid.OEScalarGrid):
    """
    Inherits OEScalarGrid to equip algebraic methods
    """

    def __init__(self, grid, initialize=False, default=0.0, title=''):
        oegrid.OEScalarGrid.__init__(self, grid)
        self.size = self.GetSize()
        self.xdim = self.GetXDim()
        self.ydim = self.GetYDim()
        self.zdim = self.GetZDim()
        self.xmin = self.GetXMin()
        self.xmax = self.GetXMax()
        self.ymin = self.GetYMin()
        self.ymax = self.GetYMax()
        self.zmin = self.GetZMin()
        self.zmax = self.GetZMax()
        self.spacing = self.GetSpacing()
        if initialize:
            if default == 0.0:
                self.ClearValues()
            else:
                self.SeValues(np.ones(self.size)*default)
            self.SetTitle(title)
        elif title:
            self.SetTitle(title)

    def GetDimension(self):
        return dict(size=self.size,
                    spacing=self.spacing,
                    min=dict(x=self.xmin, y=self.ymin, z=self.zmin),
                    max=dict(x=self.xmax, y=self.ymax, z=self.zmax),
                    dim=dict(x=self.xdim, y=self.ydim, z=self.zdim))

    def GetStats(self):
        v = self.Values()
        D = {}
        D['size'] = v.size
        D['min'] = v.min()
        D['max'] = v.max()
        D['mean'] = v.mean()
        D['sd'] = v.std()
        return D

    def ClearValues(self):
        return self.SetValues(np.zeros(self.size))

    def Values(self):
        return np.array(self.GetValues())

    def SetValues(self, v):
        oegrid.OEScalarGrid.SetValues(self, v)
        return self

    def GetCentroid(self):
        v = self.Values()
        vsum = v.sum()
        xs = np.zeros(self.size)
        ys = np.zeros(self.size)
        zs = np.zeros(self.size)
        for i in range(self.size):
            xs[i], ys[i], zs[i] = self.ElementToSpatialCoord(i)
        return (v*xs).sum()/vsum, (v*ys).sum()/vsum, (v*zs).sum()/vsum

    def Rescale(self, a=0.0, b=1.0):
        v = self.Values()
        vmin, vmax = v.min(), v.max()
        v = (v - vmin)/(vmax - vmin)*(b - a) + a
        return self.SetValues(v)

    def RescaledGrid(self, a=0, b=1):
        return AlgebraicGrid(self).Rescale(a, b)

    def Normalize(self, m=0.0, s=1.0):
        v = self.Values()
        a = v.mean()
        sd = v.std()
        v = ((v - a)/sd - m)*s
        return self.SetValues(v)

    def NormalizedGrid(self, m=0.0, s=1.0):
        return AlgebraicGrid(self).Normalize(m, s)

    def Filter(self, filter=lambda x: True):
        v = self.Values()
        for i in range(v.size):
            if not filter(v[i]):
                v[i] = 0
        return self.SetValues(v)

    def FilteredGrid(self, filter=lambda x: True):
        return AlgebraicGrid(self).Filter(filter)

    def Ternarize(self):
        v = np.sign(self.Values())
        return self.SetValues(v)

    def TernaryGrid(self):
        return AlgebraicGrid(self).Ternarize()

    def DistanceGrid(self, p, __square__=True):
        v = np.zeros(self.size)
        for i in range(v.size):
            x, y, z = self.ElementToSpatialCoord(i)
            v[i] = (x-p[0])*(x-p[0]) + (y-p[1])*(y-p[1]) + (z-p[2])*(z-p[2])
        if __square__:
            v = np.sqrt(v)
        return AlgebraicGrid(self, initialize=True).SetValues(v)

    def Distance2Grid(self, p):
        return self.DistanceGrid(p, __square__=False)

    def __neg__(self):
        v = -self.Values()
        return AlgebraicGrid(self).SetValues(v)

    def __add__(self, o):
        try:
            v = self.Values() + o.Values()
        except AttributeError:
            v = self.Values() + o
        return AlgebraicGrid(self).SetValues(v)

    def __radd__(self, o):
        return self + o

    def __sub__(self, o):
        try:
            v = self.Values() - o.Values()
        except AttributeError:
            v = self.Values() - o
        return AlgebraicGrid(self).SetValues(v)

    def __rsub__(self, o):
        return -self + o

    def __mul__(self, o):
        try:
            v = self.Values()*o.Values()
        except AttributeError:
            v = self.Values()*o
        return AlgebraicGrid(self).SetValues(v)

    def __rmul__(self, o):
        return self*o

    def __truediv__(self, o):
        try:
            v = self.Values()/o.Values()
        except AttributeError:
            v = self.Values()/o
        return AlgebraicGrid(self).SetValues(v)

    def __rtruediv__(self, o):
        v = o/self.Values()
        return AlgebraicGrid(self).SetValues(v)

    def __pow__(self, o):
        v = self.Values()**o
        return AlgebraicGrid(self).SetValues(v)

    def __IsCompatible__(self, o):
        return self.size == o.size and \
            self.xdim == o.xdim and \
            self.ydim == o.ydim and \
            self.zdim == o.zdim and \
            self.xmin == o.xmin and \
            self.xmax == o.xmax and \
            self.ymin == o.ymin and \
            self.ymax == o.ymax and \
            self.zmin == o.zmin and \
            self.zmax == o.zmax and \
            self.spacing == o.spacing

    @staticmethod
    def IsCompatible(a, b):
        return a.__IsCompatible__(b)
