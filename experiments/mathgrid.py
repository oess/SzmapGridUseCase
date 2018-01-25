from openeye import oegrid

class MathGrid(oegrid.OEScalarGrid):
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
        for i in range(self.size):
            self[i] = v

    def GetMinMax(self):
        vs = self.GetValues()
        return min(vs), max(vs)

    def GetCentroid(self):
        xs, ys, zs = 0.0, 0.0, 0.0
        vs = 0
        for i in range(self.size):
            v = self[i]
            if v == 0.0:
                continue
            x, y, z = self.ElementToSpatialCoord(i)
            xs += v*x
            ys += v*y
            zs += v*z
            vs += v
        return xs/vs, ys/vs, zs/vs

    def GetNormalizedGrid(self, n=1):
        grid = MathGrid(self, initialize=True)
        buf = self.GetValues()
        bmin, bmax = min(buf), max(buf)
        for i in range(grid.size):
            v = (buf[i] - bmin)/(bmax - bmin)*n
            grid.SetValue(i, v)
        return grid

    def GetFilteredGrid(self, filter):
        grid = MathGrid(self, initialize=True)
        buf = self.GetValues()
        for i in range(grid.size):
            grid[i] = buf[i] if filter(buf[i]) else 0
        return grid

    def __mul__(self, n):
        grid = MathGrid(self, initialize=True)
        buf = self.GetValues()
        for i in range(grid.size):
            grid[i] = n*buf[i]
        return grid

    def __rmul__(self, n):
        return n*self

    def __add__(self, other):
        grid = MathGrid(self, initialize=True)
        buf = self.GetValues()
        for i in range(grid.size):
            grid[i] = self[i] + other[i]
        return grid

    def IsCompat(self, other):
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
