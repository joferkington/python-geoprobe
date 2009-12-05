import numpy as np
import struct

from geoprobe import utilities
from volume import volume

"""
Geoprobe horizons
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
       
_pointFormat = 'f, f, f, f, B, B, B' # Big-endian byte-order!
_pointNames = 'x, y, z, conf, type, herid, tileSize'

class horizon(object):
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
        self.x = self.data.x
        self.y = self.data.y
        self.z = self.data.z
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
        self._file = _horizonFile(filename, 'r')

        self._header = self._file.readHeader()
        if self._header == "#GeoProbe Horizon V2.0 ascii\n":
            raise TypeError('Ascii horizons not currently supported')
        elif self._header != "#GeoProbe Horizon V2.0 binary\n":
            raise TypeError('This does not appear to be a valid geoprobe horizon')
        self.numpoints = self._file.numpoints
        self.data = self._file.readAllPoints()

    @property
    def grid(self):
        try: return self._grid
        except AttributeError:
            grid = self.nodata*np.ones((self.data.y.ptp()+1,self.data.x.ptp()+1),dtype=np.float32)
            I = np.array(self.x-self.xmin,np.int)
            J = np.array(self.y-self.ymin,np.int)
            for k in xrange(I.size):
                i,j,d = I[k],J[k],self.z[k]
                grid[j,i] = d
            self._grid = grid
            return grid

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
        try: import osgeo.gdal as gdal 
        except ImportError: 
            raise ImportError('Gdal not found! To use horizon.toGeotiff, gdal and the python bindings to gdal must be installed')
            
        if filename[-4:] in ['.tif','tiff']:
            format = 'GTiff'
        elif filename[-3:] == 'img':
            format = 'HFA'
        else: 
            # Assume geotiff and append ".tif"
            format = 'GTiff'
            outputFilename += '.tif'

        if vol is not None:
            if type(vol) == type('string'): vol = volume(vol)

        if nodata==None: nodata = self.nodata

        # Zscale is not 1 by default, as I want the default to be set by vol.dz 
        # and any specified value to override the default
        if zscale is None:
            if vol is None: zscale = 1
            elif vol.dz > 0: zscale = 1
            elif vol.dz < 0: zscale = -1

        ysize,xsize = self.grid.shape
        Xoffset, Yoffset = self.xmin, self.ymin

        #-- Create and write output file
        driver = gdal.GetDriverByName(format)
        dataset = driver.Create(filename,xsize,ysize,1,gdal.GDT_Float32) #One band, stored as floats

        # Georeference volume if vol is given
        if vol is not None:
            Xoffset, Yoffset = vol.model2world(Xoffset, Yoffset)
            tran = vol.transform
            dataset.SetGeoTransform( [Xoffset[0], tran[0,0], tran[0,1], 
                                      Yoffset[0], tran[1,0], tran[1,1]] )
        # If there's a nodata value, set it
        if nodata: 
            dataset.GetRasterBand(1).SetNoDataValue(nodata)
        # Write dataset
        data = self.grid
        if nodata != self.nodata: data[data == self.nodata] = nodata
        data[data != nodata] *= zscale
        dataset.GetRasterBand(1).WriteArray(data)



#-- This is currently very sloppy code... Need to clean up and document
class _horizonFile(file):
    """Basic geoprobe horizon binary file format reader"""
    def __init__(self, *args, **kwargs):
        """Accepts the same argument set as a standard python file object"""
        # Initalize the file object as normal
        file.__init__(self, *args, **kwargs)

        # Initalize attributes unique to this class...
        self._pointFormat = 'f, f, f, f, B, B, B' # Big-endian byte-order!
        self._pointNames = 'x, y, z, conf, type, herid, tileSize'
        self._lineHdrFmt = '>4f'

    def readBinary(self,fmt):
        size = struct.calcsize(fmt)
        dat = self.read(size)
        if len(dat) < size: #EOF Reached
            raise EOFError('End of File Reached')
        data = struct.unpack(fmt, dat)

        # Don't return a tuple if it only has one value
        if len(data) == 1: data = data[0]
        return data

    def readHeader(self):
        self.seek(0)
        return self.readline()

    def readPoints(self):
        numPoints = self.readBinary('>I')
        if numPoints > 0:
            points = np.rec.fromfile(self,
                    shape = numPoints,
                    formats = self._pointFormat,
                    names = self._pointNames,
                    byteorder = '>')
            return points
        # apparently, len(points) is not 0 when numPoints is 0...
        else: 
            return []

    def skipPoints(self):
        # Read number of points
        numPoints = self.readBinary('>I')
        # Slightly ugly bit to calculate the size (in bytes) of 1 point
        pointSize = struct.calcsize(self._pointFormat.replace(',',''))
        # Jump to next section
        self.seek(numPoints*pointSize, 1) 
        return numPoints

    def sectionType(self):
        # No idea what the difference between #34 and #28, #2, etc is... (pos, neg, 0, pick??)
        secFmt = '>I'
        secID = self.readBinary(secFmt)
        if secID == 19: 
            sectype = 'Points'
        else:
            sectype = 'Lines'
        return sectype

    @property
    def numpoints(self):
        try: 
            return self._numpoints
        except AttributeError:
            # Transverse file structure and sum up total number of points
            self.readHeader()
            numPoints = 0; secType = None
            # Read points section
            secType = self.sectionType() 
            numPoints += self.skipPoints() 
            # Read Lines section
            try:  
                secType = self.sectionType() 
                while True:
                    self.lineInfo()
                    numPoints += self.skipPoints()
            except EOFError:
                pass
            self._numpoints = numPoints
            return numPoints


    def lineInfo(self):
        xdir,ydir,zdir,ID = self.readBinary(self._lineHdrFmt)
        return xdir, ydir, zdir, ID

    def readAllPoints(self):
        self.readHeader()
        # Initalize an empty recarray to store things in
        points = np.recarray(shape=self.numpoints, formats=self._pointFormat, names=self._pointNames, byteorder = '>')
        lines = [] # To store line objects in
        i = 0; secType = None
        self.readHeader() # Jump to start of file, past header

        # Read points section
        secType = self.sectionType()
        currentPoints = self.readPoints()
        points[i:i+len(currentPoints)] = currentPoints
        i += len(currentPoints)

        # Read lines section
        try:
            secType = self.sectionType() 
            while True:
                lineInfo = self.lineInfo()
                currentPoints = self.readPoints()
                points[i:i+len(currentPoints)] = currentPoints
                lines.append((lineInfo, points[i:i+len(currentPoints)]))
                i += len(currentPoints)
        except EOFError:
                pass
        self.lines = lines
        return points

