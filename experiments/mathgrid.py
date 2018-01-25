from openeye import oechem, oegrid
import numpy as np

class MathGrid(oegrid.OEScalarGrid):
    """
    Inherits OEScalarGrid to equip arithmetic methods
    """

    def __init__(self, grid, initialize=False, v=0.0):
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
            self.ClearValues(v)

    def ClearValues(self, v=0.0):
        self.SetValues(np.ones(self.GetSize())*v)

    def GetMinMax(self):
        v = np.array(self.GetValues())
        return v.min(), v.max()

    def GetCentroid(self):
        v = np.array(self.GetValues())
        vsum = v.sum()
        xs = np.zeros(self.size)
        ys = np.zeros(self.size)
        zs = np.zeros(self.size)
        for i in range(self.size):
            xs[i], ys[i], zs[i] = self.ElementToSpatialCoord(i)
        return (v*xs).sum()/vsum, (v*ys).sum()/vsum, (v*zs).sum()/vsum

    def GetNormalizedGrid(self, n=1):
        grid = MathGrid(self, initialize=True)
        v = np.array(self.GetValues())
        vmin, vmax = v.min(), v.max()
        v = (v - vmin)/(vmax - vmin)*n
        grid.SetValues(v)
        return grid

    def GetFilteredGrid(self, filter):
        grid = MathGrid(self, initialize=True)
        v = np.array(self.GetValues())
        for i in range(v.size):
            if not filter(v[i]):
                v[i] = 0
        grid.SetValues(v)
        return grid

    def GetDistanceGrid(self, p):
        grid = MathGrid(self, initialize=True)
        v = np.zeros(self.size)
        for i in range(v.size):
            x, y, z = self.ElementToSpatialCoord(i)
            v[i] = (x-p[0])*(x-p[0]) + (y-p[1])*(y-p[1]) + (z-p[2])*(z-p[2])
        v = np.sqrt(v)
        grid.SetValues(v)
        return grid

    def GetDistance2Grid(self, p):
        grid = MathGrid(self, initialize=True)
        v = np.zeros(self.size)
        for i in range(v.size):
            x, y, z = self.ElementToSpatialCoord(i)
            v[i] = (x-p[0])*(x-p[0]) + (y-p[1])*(y-p[1]) + (z-p[2])*(z-p[2])
        grid.SetValues(v)
        return grid

    def __sub__(self, other):
        grid = MathGrid(self, initialize=True)
        if isinstance(other, MathGrid):
            v = np.array(self.GetValues()) - np.array(other.GetValues())
        else:
            v = np.array(self.GetValues()) - other
        grid.SetValues(v)
        return grid

    def __rsub__(self, other):
        grid = MathGrid(self, initialize=True)
        v = other - np.array(self.GetValues())
        grid.SetValues(v)
        return grid

    def __mul__(self, other):
        grid = MathGrid(self, initialize=True)
        if isinstance(other, MathGrid):
            v = np.array(self.GetValues())*np.array(other.GetValues())
        else:
            v = np.array(self.GetValues())*other
        grid.SetValues(v)
        return grid

    def __rmul__(self, other):
        return self*other

    def __add__(self, other):
        grid = MathGrid(self, initialize=True)
        if isinstance(other, MathGrid):
            v = np.array(self.GetValues()) + np.array(other.GetValues())
        else:
            v = np.array(self.GetValues()) + other
        grid.SetValues(v)
        return grid

    def __radd__(self, other):
        grid = MathGrid(self, initialize=True)
        v = other + np.array(self.GetValues())
        grid.SetValues(v)
        return grid

    def __div__(self, other):
        grid = MathGrid(self, initialize=True)
        if isinstance(other, MathGrid):
            v = np.array(self.GetValues())/np.array(other.GetValues())
        else:
            v = np.array(self.GetValues())/other
        grid.SetValues(v)
        return grid


    def __IsMathCompatible__(self, other):
        return self.size == other.size and \
            self.xdim == other.xdim and \
            self.ydim == other.ydim and \
            self.zdim == other.zdim and \
            self.xmin == other.xmin and \
            self.xmax == other.xmax and \
            self.ymin == other.ymin and \
            self.ymax == other.ymax and \
            self.zmin == other.zmin and \
            self.zmax == other.zmax and \
            self.spacing == other.spacing

    @staticmethod
    def IsMathCompatible(a, b):
        return a.__IsMathCompatible__(b)
