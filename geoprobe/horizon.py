import numpy as np
import struct

#-- Imports from local files --------------------------------------
import utilities
from volume import volume
from common import BinaryFile, StaticCache

class horizon(object):
    """Reads a geoprobe horizon from disk.

    horizon.x, horizon.y, and horizon.z are the x,y, and z coordinates
    stored in the horizon file (in model coordinates, i.e. inline/crossline).  

    horizon.grid is a 2d numpy array of the z-values in the horizon filled
    with horizon.nodata in regions where there aren't any z-values

    horizon.{x,y,z}min and horizon.{x,y,z}max are the min's and max's 
    of horizon.grid
    """
    def __init__(self, input):
        """Takes either a filename or a numpy array"""

        if type(input) == type('String'):
            self._readHorizon(input)    
        else:
            # For the moment, just assume it's a numpy array
            # Possibly pass a list of recarrays for lines?
            raise TypeError('Only reading is suppored at this time.  Input must be a valid file name')

        # For gridding:
        self.nodata = -9999

        # Read horizon and initalize various properties
        self.xmin = self.x.min()
        self.ymin = self.y.min()
        self.zmin = self.z.min()
        self.xmax = self.x.max()
        self.ymax = self.y.max()
        self.zmax = self.z.max()
        # Need to make dx, dy, and dz properties... 
        # How do we determine spacing without a volume?
        #    d = np.abs(np.diff(self.x)); np.mean(d[d!=0]) (ideally, mode)?

    def _readHorizon(self,filename):
        self._file = HorizonFile(filename, 'r')

        self._header = self._file.readHeader()
        if self._header == "#GeoProbe Horizon V2.0 ascii\n":
            raise TypeError('Ascii horizons not currently supported')
        elif self._header != "#GeoProbe Horizon V2.0 binary\n":
            raise TypeError('This does not appear to be a valid geoprobe horizon')
        self.data = self._file.points

    @property
    def numpoints(self):
        return self.data.size

    #-- x,y,z properties -------------------------------------------------------
    def _get_coord(self, name):
        return self.data[name]
    def _set_coord(self, name, value):
        self.data[name] = value
    x = property(lambda self: self._get_coord('x'),
            lambda self, value: self._set_coord('x', value),
            'X-coordinates of all points stored in the horizon')
    y = property(lambda self: self._get_coord('y'),
            lambda self, value: self._set_coord('y', value),
            'Y-coordinates of all points stored in the horizon')
    z = property(lambda self: self._get_coord('z'),
            lambda self, value: self._set_coord('z', value),
            'Z-coordinates of all points stored in the horizon')
    #---------------------------------------------------------------------------

    #-- Grid Property ---------------------------------------------------------
    @StaticCache
    def _get_grid(self):
        """An nx by ny numpy array (dtype=float32) of the z values contained
        in the horizon file"""
        x, y, z = self.x, self.y, self.z
        grid = np.ones((y.ptp() + 1, x.ptp() +1 ), dtype=np.float32)
        grid *= self.nodata
        I = np.array(x - x.min(), dtype=np.int)
        J = np.array(y - y.min(), dtype=np.int)
        for k in xrange(I.size):
            i, j, d = I[k], J[k], z[k]
            grid[j,i] = d
        return grid
    def _set_grid(self, value):
        self._grid = value
    grid = property(_get_grid, _set_grid)
    #--------------------------------------------------------------------------

    def strikeDip(self, vol=None, velocity=None, independent='z'):
        """
        Returns a strike and dip of the horizon following the Right-hand-rule. 
        Input:
            vol (optional): A geoprobe volume object
                If specified, the x, y, and z units will be converted
                to world units based on the volume's header.
            velocity (optional): Velocity in meters/second
                If specified, the z units will be converted from time
                into depth using the velocity given.  Assumes the z
                units are milliseconds!!
        independent (optional): Independent variable to use when 
                fitting the plane. Defaults to 'z', as most horizons
                are assumed to be closer to horizontal than vertical
                Set to None to choose the best option (slower) or 'x'
                or 'y' when fitting near-vertical horizons.
        Output:
            strike, dip
        """
        return utilities.points2strikeDip(self.x, self.y, self.z, 
                                          vol=vol, velocity=velocity, 
                                          independent=independent)

    def toGeotiff(self, filename, vol=None, nodata=None, zscale=None):
        """
        Create and write a geotiff file from the geoprobe horizon object.
        The Z values in the output tiff will be stored as 32bit floats.
        Input:
            filename:  Output filename
            vol (optional): A geoprobe volume object or path to a geoprobe volume file. 
                If vol is specified, the geotiff will be georeferenced based on the
                data in the volume header (and will therefore be in same projection
                as the volume's world coordinates).  Otherwise the geotiff is created
                using the model coordinates stored in the geoprobe horizon file.
            nodata (default=self.nodata (-9999)): Value to use for NoData.
            zscale (optional): Scaling factor to use for the Z-values.  If vol is
                specified, and vol.dz is negative, this defaults to -1.  Otherwise
                this defaults to 1.
        """
        if vol is not None:
            if type(vol) == type('string'): 
                vol = volume(vol)
            Xoffset, Yoffset = vol.model2world(self.xmin, self.ymin)
            transform = vol
        else:
            Xoffset, Yoffset = 0,0
            transform = None

        if nodata==None: 
            nodata = self.nodata

        # Zscale is not 1 by default, as I want the default to be set by vol.dz 
        # and any specified value to override the default
        if zscale is None:
            if vol is None: zscale = 1
            elif vol.dz > 0: zscale = 1
            elif vol.dz < 0: zscale = -1

        data = self.grid
        if nodata != self.nodata: data[data == self.nodata] = nodata
        data[data != nodata] *= zscale

        utilities.array2geotiff(data, filename, nodata=nodata, extents=(Xoffset, Yoffset), transform=transform)



#-- This is currently very sloppy code... Need to clean up and document
class HorizonFile(BinaryFile):
    """Basic geoprobe horizon binary file format reader

        Disk layout of Geoprobe horizons
          Reverse engineered by JDK, Feb. 2009
        Format descrip:
            1 ascii line w/ version (terminated with newline)
            There are two "sections" in every file. 
            The first section contains x,y,z points making a "filled" surface
                (This is basically a sparse matrix)
            The second section contains lines (manual picks)
            Both section types have a 4 byte header (seems to be >I?)
                The actual values in this header seem to be complelely meaningless (??)
            subsections
                if surface:
                    info: (>I) Number of points
                if line:
                    info: (>4fI) xdir,ydir,zdir,ID,numPoints
            Point Format in all sections: (>4f3B)
                x,y,z,confidence,type,heridity,tileSize
    """
    _pointFormat = ('>f', '>f', '>f', '>f', '>B', '>B', '>B')
    _pointNames = ('x', 'y', 'z', 'conf', 'type', 'herid', 'tileSize')
    _lineHdrFmt = '>4f'

    def __init__(self, *args, **kwargs):
        """Accepts the same argument set as a standard python file object"""
        # Initalize the file object as normal
        file.__init__(self, *args, **kwargs)

        # Build a dtype definition 
        self.point_dtype = []
        for name, fmt in zip(self._pointNames, self._pointFormat):
            self.point_dtype.append((name,fmt))

        # Size in Bytes of a point (x,y,z,conf,type,...etc)
        self._pointSize = sum(map(struct.calcsize, self._pointFormat))

    def readHeader(self):
        self.seek(0)
        return self.readline()

    def readPoints(self):
        numPoints = self.readBinary('>I')
        points = np.fromfile(self, count=numPoints, dtype=self.point_dtype)
        return points

    def sectionType(self):
        secFmt = '>I'
        secID = self.readBinary(secFmt)
        if secID == 19: 
            sectype = 'Points'
        else:
            # otherwise, secID is the number of lines
            sectype = 'Lines'
        return sectype

    def lineInfo(self):
        xdir,ydir,zdir,ID = self.readBinary(self._lineHdrFmt)
        return xdir, ydir, zdir, ID

    @property
    @StaticCache
    def points(self):
        """A numpy array with the fields ('x', 'y', 'z', 'conf', 'type', 'herid', 'tileSize') 
        for each point in the horizon (regardless of whehter it's a manual pick (line) or 
        surface)"""
        # Note: The total number of points in the file is not directly stored 
        #   on disk. Therefore, we must read through the entire file, store 
        #   each section's points in a list, and then create a contigious array
        #   from them.  Using numpy.append is much simpler, but quite slow.

        # Jump to start of file, past header
        self.readHeader() 

        # Read points section
        secType = self.sectionType()
        temp_points = [self.readPoints()]

        # Read lines section
        try:
            secType = self.sectionType() 
            while True:
                lineInfo = self.lineInfo()
                currentPoints = self.readPoints()
                temp_points.append(currentPoints)
        except EOFError:
            pass

        # Create a single numpy array from the list of arrays (temp_points)
        numpoints = sum(map(np.size, temp_points))
        points = np.zeros(numpoints, dtype=self.point_dtype)
        i = 0
        for item in temp_points:
            points[i : i + item.size] = item
            i += item.size

        return points

