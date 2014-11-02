"""Really simplistic Landmark EzFault reader
JDK 08/03/09"""
# TODO: This whole thing needs documentation!
# TODO: Methods to convert to grid / triangle strips
# TODO: Write support... Eventually

import numpy as np

# Local files
import utilities
from volume import volume

class ezfault(object):
    """Simple geoprobe ezfault reader."""
    def __init__(self, filename):
        """Takes a filename as input"""
        self._infile = file(filename, 'r')
        self._readHeader()


    def _readHeader(self):
        """ Blindly read the first 6 lines and set properties based on them
        This doesn't check order or sanity of any sort...."""
        # This really is a bit of a dumb and brittle parser... Doesn't parse
        #   the {'s and }'s, just the position of the data in each line.

        self._infile.seek(0)

        #1st Line: "#GeoProbe Surface V1.0 ascii"
        # Identification string
        self.IDstring = self._infile.readline()

        #2nd Line: "FACE %i"
        #  Face Direction (1=X, 2=Y, 3=Z)
        self.face = self._infile.readline().strip().split()
        self.face = int(self.face[1])

        #3rd Line: "BOX {%xmin %ymin %zmin %xmax %ymax %zmax}"
        #  Bounding Box (for some weird reason, this is in volume indicies
        #  instead of model units!!)
        self.box = self._infile.readline().strip()[5:-1].split()
        self.box = [int(item) for item in self.box]

        #4th Line: "ORIENT  {%xcomp %ycomp %zcomp}"
        #  Orientation of ribs (?)
        self.orientation = self._infile.readline().strip()[9:-1].split()
        self.orientation = [float(item) for item in self.orientation]

        #5th Line: "COLOR %hexcode"
        self.color = self._infile.readline().strip()[6:]

        #6th Line: "TRANSPARENCY %value"
        self.transparency = self._infile.readline().strip()[13:]
        self.transparency = float(self.transparency)

    #Can't use self.box to set xmin, xmax, etc...
    # (self.box is in volume indicies instead of model units!)
    xmin = property(lambda self: self.points.x.min(),
            doc='Mininum X-coordinate (in inline/crossline)')
    ymin = property(lambda self: self.points.y.min(),
            doc='Mininum Y-coordinate (in inline/crossline)')
    zmin = property(lambda self: self.points.z.min(),
            doc='Mininum Z-coordinate (in inline/crossline)')
    xmax = property(lambda self: self.points.x.max(),
            doc='Maximum X-coordinate (in inline/crossline)')
    ymax = property(lambda self: self.points.y.max(),
            doc='Maximum Y-coordinate (in inline/crossline)')
    zmax = property(lambda self: self.points.z.max(),
            doc='Maximum Z-coordinate (in inline/crossline)')

    @property
    def ribs(self):
        """Reads points from a rib"""
        """Format looks like this:
        CURVE
        {
        SUBDIVISION 2
        VERTICES
          {
            { 2657.630371 5162.239258 2845.000000 },
            { 2652.424072 5187.238281 2845.000488 },
            { 2648.202637 5206.083496 2845.000488 },
            { 2645.693604 5213.675781 2845.000488 },
          }
        }
        """
        self._readHeader()
        while True:
            # Make sure that the next line is "CURVE"
            line = self._infile.readline()
            if line == '':
                raise StopIteration #End of File
            elif line.strip() != 'CURVE':
                raise ValueError('Expected next line to be "CURVE",'\
                                ' got %s' % line.strip())

            # Skip the next 4 lines, and store the 5th for processing
            for i in range(5): line = self._infile.readline()

            # Read points
            verts = []
            while line.strip() is not '}':
                line = line.strip().split()
                verts.append( tuple([float(item) for item in line[1:4]]) )
                line = self._infile.readline()
            # Read the other '}'
            line = self._infile.readline()

            yield verts

    @property
    def points(self):
        """Returns a numpy array of all points in the file"""
        # Have we already done this:
        try:
            return self._allPoints
        except AttributeError:
            dat = []
            self._readHeader()
            for rib in self.ribs:
                dat.extend(rib)
            self._allPoints = np.rec.fromrecords(dat, names='x,y,z')
            return self._allPoints

    def strikeDip(self, vol=None, velocity=None):
        """
        Returns a strike and dip of the ezfault following the Right-hand-rule.
        Input:
            vol (optional): A geoprobe volume object
                If specified, the x, y, and z units will be converted
                to world units based on the volume's header.
            velocity (optional): Velocity in meters/second
                If specified, the z units will be converted from time
                into depth using the velocity given.  Assumes the z
                units are milliseconds!!
        Output:
            strike, dip
        """
        return utilities.points2strikeDip(self.points.x, self.points.y,
                        self.points.z, vol=vol, velocity=velocity)

    def interpolate(self, xi, yi):
        from scipy import interpolate
        # Need to write a norm function that calculates distance from a rib...

        """
        def interp(x1,x2,x3, x1i, x2i):
            spline = interpolate.Rbf(x1, x2, x3, function='thin-plate', smooth=0)
            return spline(x1i,x2i)

        try:
            zi = interp(self.points.x, self.points.y, self.points.z, xi, yi)
        except np.linalg.linalg.LinAlgError:
            zi = interp(self.points.y, self.points.x, self.points.z, yi, xi)
        """

        # Segfaults... Problems with the way scipy is compiled?
        tck = interpolate.bisplrep(self.points.x, self.points.y, self.points.z)
        zi = interpolate.bisplev(yi, xi, tck)


        """
        spline = interpolate.Rbf(self.points.x, self.points.y, self.points.z,
                                 function='thin-plate', smooth=0)
        zi = spline(xi,yi)
        """

        return zi



    def grid(self, dx=None, dy=None, extents=None, vol=None):

        if extents is None:
            xstart, xstop = self.xmin, self.xmax
            ystart, ystop = self.ymin, self.ymax
        else:
            xstart, xstop, ystart, ystop = extents

        if vol is not None:
            # Interpolate ezfault at volume indicies
            if type(vol) == type('String'):
                vol = volume(vol)
            dx, dy = abs(vol.dx), abs(vol.dy)

            # Make sure we start at a volume index
            xstart, xstop = [vol.index2model(vol.model2index(item))
                            for item in [xstart, xstop] ]
            ystart, ystop= [vol.index2model(
                                vol.model2index(item, axis='y'),
                                axis='y')
                            for item in [ystart, ystop] ]

        else:
            if dx is None: dx = 1
            if dy is None: dy = dx

        X = np.arange(xstart, xstop, dx)
        Y = np.arange(ystart, ystop, dy)
        # Not needed w/ new way of interpolating?
#        X,Y = np.meshgrid(X,Y)
        grid = self.interpolate(X,Y)

#        grid = Z.reshape(X.shape)

        return grid




