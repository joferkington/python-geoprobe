import numpy as np
import struct

#-- Imports from local files --------------------------------------
from volume import volume
from common import BinaryFile, array2geotiff

#-- Build dtype for points ----------------------------------------
_point_format = ('>f', '>f', '>f', '>f',  '>B',   '>B',    '>B')
_point_names =  ('x',  'y',  'z', 'conf', 'type', 'herid', 'tileSize')
_point_dtype = zip(_point_names, _point_format)

class horizon(object):
    """Reads a geoprobe horizon from disk.

    horizon.x, horizon.y, and horizon.z are the x,y, and z coordinates
    stored in the horizon file (in model coordinates, i.e. inline/crossline).  

    horizon.grid is a 2d numpy array of the z-values in the horizon filled
    with horizon.nodata in regions where there aren't any z-values

    horizon.{x,y,z}min and horizon.{x,y,z}max are the min's and max's 
    of horizon.grid
    """
    def __init__(self, *args, **kwargs):
        """Takes either a filename or a numpy array"""

        # If __init__ is just passed a string, assume it's a filename
        # and make a horizon object by reading from disk
        if (len(args) == 1) and isinstance(args[0], str):
            self._readHorizon(args[0])    

        # Otherwise, pass the args on to _make_horizon_from_data for
        # parsing
        else:
            self._parse_new_horizon_input(*args, **kwargs)

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

    def _readHorizon(self, filename):
        self._file = HorizonFile(filename, 'r')
        self._header = self._file.readHeader()

        if self._header == "#GeoProbe Horizon V2.0 ascii\n":
            raise TypeError('Ascii horizons not currently supported')
        elif self._header != "#GeoProbe Horizon V2.0 binary\n":
            raise TypeError('This does not appear to be a valid geoprobe horizon')
        
        self.data = self._file.readAll()

        # Surface and line attributes
        self.surface = self._file.surface
        self.lines = self._file.lines

        if self.data.size == 0:
            raise ValueError('This file does not contain any points!')

    def _parse_new_horizon_input(self, *args, **kwargs):
        """Parse input when given something other than a filename"""
        
        #-- Parse Arguments ---------------------------------------------------
        if len(args) == 0:
            pass
        elif len(args) == 1:
            # Assume argument is data (numpy array with dtype of _point_dtype)
            self.data = np.asarray(args[0], dtype=_point_dtype)
            self.surface = self.data
            
        elif len(args) == 2:
            # Assume arguments are surface + lines
            init_from_surface_lines(self, surface=args[0], lines=args[1])

        elif len(args) == 3:
            # Assume arguments are x, y, and z arrays
            init_from_xyz(self, *args)

        else:
            print args
            raise ValueError('Invalid number of input arguments.')

        #-- Parse keyword arguments -------------------------------------------
        if len(kwargs) == 0:
            pass
        elif ('x' in kwargs) and ('y' in kwargs) and ('z' in kwargs):
            self._init_from_xyz(kwargs['x'], kwargs['y'], kwargs['z'])
            
        elif ('surface' in kwargs) and ('lines' in kwargs):
            self._init_from_surface_lines(kwargs['surface'], kwargs['lines'])
            
        elif 'surface' in kwargs:
            self._init_from_surface_lines(surface=surface)

        elif 'lines' in kwargs:
            self._init_from_surface_lines(lines=lines)

        else:
            raise ValueError('Invalid keyword arguments. You must specify x,y,&z, surface, or lines')

    def _init_from_xyz(self, x, y, z):
        """Make a new horizon object from x, y, and z arrays"""
        x,y,z = [np.asarray(item, dtype=np.float32) for item in [x,y,z]]
        if x.size == y.size == z.size:
            self.data = np.zeros(x.size, dtype=_point_dtype)
            self.x = x
            self.y = y
            self.z = z
            self.surface = data
        else:
            raise ValueError('x, y, and z arrays must be the same length')

    def init_from_surface_lines(self, surface=None, lines=None):
        """Make a new horizon object from either a surface array or a list of line arrays"""
        if surface is not None:
            surface = np.asarray(surface, dtype=_point_dtype)

        # Calculate total number of points
        numpoints = surface.size if surface else 0
        if lines is not None:
            for info, line in lines:
                numpoints += line.size

        # Make self.data and make self.lines & self.surface views into self.data
        self.data = np.zeros(numpoints, dtype=_point_dtype)

        array_list = []
        if surface is not None:
            array_list.append((None, surface))
        if lines is not None:
            array_list.extend(lines)

        i = 0
        for info, item in array_list:
            self.data[i:i+item.size] = item
            if (surface is not None) and (i == 0):
                self.surface = self.data[i:i+item.size]
            else:
                self.lines.append(info, data[i:i+item.size])

    def write(self, filename):
        """
        Write the horizon to a new file ("filename")
        """
        self._file = HorizonFile(filename, 'w')

        # If self.lines isn't set, default to []
        try:
            self._file.lines = self.lines
        except AttributeError:
            self._file.lines = []

        # If self.surface isn't set, default to an empty numpy array
        try:
            self._file.surface = self.surface
        except AttributeError:
            self._file.surface = np.zeros(0, dtype=_point_dtype)

        self._file.writeAll()

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
    def _get_grid(self):
        """An nx by ny numpy array (dtype=float32) of the z values contained
        in the horizon file"""
        try:
            return self._grid
        except AttributeError:
            x, y, z = self.x, self.y, self.z
            grid = np.ma.masked_all((y.ptp() + 1, x.ptp() +1 ), dtype=np.float32)
            grid.fill_value = self.nodata
            I = np.array(x - x.min(), dtype=np.int)
            J = np.array(y - y.min(), dtype=np.int)
            grid[J,I] = z
            self._grid = grid
            return self._grid
    def _set_grid(self, value):
        self._grid = np.ma.asarray(value)
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
        data.fill_value = nodata
        data *= zscale

        array2geotiff(data, filename, nodata=nodata, extents=(Xoffset, Yoffset), transform=transform)



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
                The first section (surface) always (?) has a section header 
                value of '\x00\x00\x00\x13' (unpacks to 19)
                The section section (manual picks) contains the number of 
                manually picked lines in the file (packed as >I).
            Subsections
                The first section only has one subsection, a "filled" surface
                    Surface header: (>I) Number of points in the surface
                The second section contains "numlines" subsections containing
                manual picks (lines):
                    Line header: (>4f) xdir,ydir,zdir,ID
            Point Format in all sections: (>4f3B)
                x,y,z,confidence,type,heridity,tileSize
    """
    _sectionHdrFmt = '>I'
    _surfaceHdrFmt = '>I'
    _lineHdrFmt = '>4f'

    def __init__(self, *args, **kwargs):
        """Accepts the same argument set as a standard python file object"""
        # Initalize the file object as normal
        file.__init__(self, *args, **kwargs)

    def readHeader(self):
        self.seek(0)
        return self.readline()

    def readPoints(self):
        numPoints = self.readBinary(self._surfaceHdrFmt)
        points = np.fromfile(self, count=numPoints, dtype=_point_dtype)
        return points

    def readSectionHeader(self):
        return self.readBinary(self._sectionHdrFmt)

    def readLineHeader(self):
        xdir,ydir,zdir,ID = self.readBinary(self._lineHdrFmt)
        return xdir, ydir, zdir, ID

    def readAll(self):
        """
        Reads in the entire horizon file and returns a numpy array with the fields 
        ('x', 'y', 'z', 'conf', 'type', 'herid', 'tileSize') for each point in the 
        horizon.
        """
        # Note: The total number of points in the file is not directly stored 
        #   on disk. Therefore, we must read through the entire file, store 
        #   each section's points in a list, and then create a contigious array
        #   from them.  Using numpy.append is much simpler, but quite slow.

        # Jump to start of file, past header
        self.readHeader() 

        # Read points section
        self.readSectionHeader() # Should always return 19
        surface = self.readPoints()
        temp_points = [surface]

        # Read lines section
        line_info = [None]
        self.numlines = self.readSectionHeader() 
        for i in xrange(self.numlines):
            line_info.append(self.readLineHeader())
            currentPoints = self.readPoints()
            temp_points.append(currentPoints)

        # Create a single numpy array from the list of arrays (temp_points)
        numpoints = sum(map(np.size, temp_points))
        points = np.zeros(numpoints, dtype=_point_dtype)
        self.lines = []
        i = 0
        for info, item in zip(line_info, temp_points):
            points[i : i + item.size] = item
            # self.surface is a view into the first part of the points array
            if i == 0:
                self.surface = points[i:i+item.size]
            # self.lines is a list of tuples, the first item is a tuple of
            # (xdir,ydir,zdir,ID) where <xdir,ydir,zdir> form a vector in
            # the direction of the line. The second item is a view into the 
            # points array containg the relevant x,y,z,etc points.
            else:
                self.lines.append((info, points[i:i+item.size]))
            i += item.size

        return points

    def writeHeader(self):
        header = "#GeoProbe Horizon V2.0 binary\n"
        self.seek(0)
        self.write(header)

    def writePoints(self, points):
        numPoints = points.size
        self.writeBinary(self._surfaceHdrFmt, numPoints)
        points.tofile(self)

    def writeLineHeader(self, line_hdr):
        self.writeBinary(self._lineHdrFmt, line_hdr)

    def writeSectionHeader(self, sec_hdr):
        self.writeBinary(self._sectionHdrFmt, sec_hdr)

    def writeAll(self):
        self.writeHeader()
        self.writeSectionHeader(19)
        self.writePoints(self.surface)
        self.writeSectionHeader(len(self.lines))
        for (info, line) in self.lines:
            self.writeLineHeader(info)
            self.writePoints(line)



