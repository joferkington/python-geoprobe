#! /usr/bin/python
import struct
import textwrap
import numpy as np

# Local imports
import volume

#-- Miscellaneous -------------------------------------------------------------
def format_headerDef_docs(headerDef, initial_indent=8, subsequent_indent=12):
    """
    Format the attributes contained in a headerDef for pretty printing 
    (i.e. for use in docstrings)
    """
    attribute_docs = ''
    initial_indent *= ' '
    subsequent_indent *= ' '

    for key in sorted(headerDef.keys()):
        value = headerDef[key]
        default = value['default']
        if isinstance(default, basestring):
            default = default.strip()

        doc = '%s: %s (default=%s)' % (key, value['doc'], repr(default))
        doc = textwrap.fill(doc, initial_indent=initial_indent, 
                subsequent_indent=subsequent_indent)

        if not key.startswith('_'):
            attribute_docs += doc + '\n'

    return attribute_docs

class cached_property(object):
    """
    A decorator class that ensures that properties are only evaluated once.
    From: <http://code.activestate.com/recipes/363602-lazy-property-evaluation/>
    """
    def __init__(self, calculate_function):
        self._calculate = calculate_function

    def __get__(self, obj, _=None):
        if obj is None:
            return self
        value = self._calculate(obj)
        setattr(obj, self._calculate.func_name, value)
        return value

#-- Raw reading and writing ---------------------------------------------------
class BinaryFile(file):
    """
    Automatically packs or unpacks binary data according to a format
    when reading or writing.
    """
    def __init__(self, *args, **kwargs):
        """
        Initialization is the same as a normal file object
        %s""" % file.__doc__
        file.__init__(self, *args, **kwargs)

    def readBinary(self,fmt):
        """
        Read and unpack a binary value from the file based
        on string fmt (see the struct module for details).
        """
        size = struct.calcsize(fmt)
        data = self.read(size)
        # Reading beyond the end of the file just returns ''
        if len(data) != size:
            raise EOFError('End of file reached')
        data = struct.unpack(fmt, data)

        for item in data:
            # Strip trailing zeros in strings 
            if isinstance(item, basestring):
                item = item.strip('\x00')

        # Unpack the tuple if it only has one value
        if len(data) == 1: data = data[0]

        return data

    def writeBinary(self, fmt, dat):
        """Pack and write data to the file according to string fmt."""
        # Try expanding input arguments (struct.pack won't take a tuple)
        try: 
            dat = struct.pack(fmt, *dat) 
        except (TypeError, struct.error): 
            # If it's not a sequence (TypeError), or if it's a 
            # string (struct.error), don't expand.
            dat = struct.pack(fmt, dat) 
        self.write(dat)


#-- Functions in the utilities namespace---------------------------------------
def array2geotiff(data, filename, nodata=-9999, transform=None, extents=None):
    """
    Write a geotiff file called "filename" from the numpy array "data".
    Input:
        data: A 2D numpy array with shape (nx,ny)
        filename: The filename of the output geotiff
        nodata:  The nodata value of the array (defaults to -9999)
        transform (optional): Either a 3x2 numpy array describing 
            an affine transformation to georeference the geotiff 
            with, or an object that has a transform attribute 
            containing such an array (e.g. a geoprobe.volume object)
        extents (optional): Either a 2-tuple of the coords of the 
            lower left corner of "data" (xmin, ymin) or an object 
            with xmin, xmax, ymin, ymax attributes (e.g. a geoprobe 
            volume object). This overrides the x and y offsets in 
            "transform". If transform is not given, it is ignored.
            Note that if you're using a horizon object as input to
            "extents", you'll need to call volume.model2world on
            horizon.xmin & ymin, as a horizon's coordinates are
            stored in model (inline, crossline) space.
    """
    try: import osgeo.gdal as gdal 
    except ImportError: 
        raise ImportError('Gdal not found! To use horizon.toGeotiff, gdal and the python bindings to gdal must be installed')

    if transform is not None:
        try:
            transform = transform.transform
        except AttribueError:
            # Else, assume it's already 2x3 array containg an affine transformation
            pass
    if extents is not None:
        try:
            xmin,ymin = extents
        except:
            try:
                xmin, ymin = extents.xmin, extents.ymin
            except AttributeError:
                raise ValueError('Extents must be either a 2-tuple of xmin,ymin or an object with xmin,ymin properities')
    else:
        xmin, ymin = 0,0
            
    # Determine output format from filename
    if filename[-4:] in ['.tif','tiff']:
        format = 'GTiff'
    elif filename[-3:] == 'img':
        format = 'HFA'
    else: 
        # Assume geotiff and append ".tif"
        format = 'GTiff'
        outputFilename += '.tif'

    # Need to change horizon.grid, and change this when I do... 
    # Everything should probably expect a x by y array to maintain consistency with volume.data
    ysize,xsize = data.shape

    #-- Create and write output file
    driver = gdal.GetDriverByName(format)
    dataset = driver.Create(filename,xsize,ysize,1,gdal.GDT_Float32) #One band, stored as floats

    # Georeference volume if vol is given
    if transform is not None:
        dataset.SetGeoTransform( [xmin, transform[0,0], transform[0,1], 
                                  ymin, transform[1,0], transform[1,1]] )
    # Set the nodata value
    dataset.GetRasterBand(1).SetNoDataValue(nodata)
    # Write dataset
    dataset.GetRasterBand(1).WriteArray(data)


def points2strikeDip(x, y, z, vol=None, velocity=None, independent=None):
    """
    Takes a point cloud defined by 3 vectors and returns a strike and dip 
    following the Right-hand-rule. 
    Input:
        x, y, z: numpy arrays or lists containg the x, y, and z 
            coordinates, respectively
        vol (optional): A geoprobe volume object
            If specified, the x, y, and z units will be converted
            to world units based on the volume's header.
        velocity (optional): Velocity in meters/second
            If specified, the z units will be converted from time
            into depth using the velocity given.  Assumes the z
            units are milliseconds!!
        independent (optional): Independent variable to use when fitting
            the plane.  Defaults to None, which automatically chooses
            the best option.  Choose 'x', 'y', or 'z' (or 0,1,2) to 
            speed the operation up when constraints can be placed on 
            the orientation of the plane. See docstring in fitPlane
            for more detailed information.
    Output:
        strike, dip (in degrees)
    """

    # Get x,y, and z, converting to world coords if necessary
    if vol is not None:
        # If given a string, assume it's the filename of a volume
        if type(vol) == type('String'):
            vol = volume.volume(vol)
        # Convert coords
        x,y = vol.model2world(x, y)

    #Negate z if depth is positive (which it usually is)
    if vol is None or vol.dz < 0: 
        z = -z

    # Convert to depth using velocity (assumes velocity is in m/s and Z is in milliseconds)
    if velocity is not None:
        z = 0.5 * velocity * z / 1000  

    a,b,c,d = fitPlane(x,y,z, independent) 
    strike, dip = normal2SD(a,b,c)

    return strike, dip

def fitPlane(x,y,z,independent=None):
    """Fits a plane to a point cloud using least squares. Should handle 
    vertical and horizontal planes properly.
    Input:
        x,y,z: 
           numpy arrays of x, y, and z, respectively
        independent (optional): 
          specify the formulation used to solve for the plane. This is 
          useful to speed up the calculation when you have constraints
          on the result. 
          Values may be:
              None (default): Choose the best one automatically
            0 or 'x': Use x = ay + bz + c. Best when the result is 
                        nearly parallel to the z-y plane. Will fail
                    if used for near-horizontal planes.
            1 or 'y': Use y = ax + bz + c. Best when the result is
                    nearly parallel to the z-x plane. Will fail
                    if used for near-horizontal planes.
            2 or 'z': Use z = ax + by + c. Best when the result is 
                    nearly parallel to the x-y plane. Will fail 
                    if used for near-vertical planes.
    Returns:
        a,b,c,d where 0 = ax + by + cz +d 
        (The normal vector is < a, b, c >)
    """

    #---------------------------------------------------------------------------------
    #-- This function looks needlessly complex, but if you just try to fit z = ax+by+c,
    #-- you can't properly handle near-vertical planes, and if you use the other two,
    #-- you can't handle near-horizontal planes. Instead, we use all 3 methods of 
    #-- fitting a plane to points and choose the one with the lowest residual.
    #-- The kwarg independent allows one to specifiy which formulation to use.  This does
    #-- speed the calculation up, however, its main purpose is actually to allow horizons
    #-- to always be fit by z=ax+by+c. If the "best" fit is actually chosen, long narrow
    #-- horizons will be best fit by a vertical plane.  As it is not possible to store
    #-- vertical horizons in geoprobe's format, long, narrow horizons are far more likely
    #-- than near-vertical horizons (and in fact, this option was added because
    #-- of a seafloor horizon that was "best" fit by a vertical plane if the lowest 
    #-- residual was chosen)
    #---------------------------------------------------------------------------------------
    # Alternatively, one could calculate the condition number to choose the best G-matrix.
    # However, this does not consistently pick the best formulation when the G-matrix is very
    # large (millions of rows).  (I'm not sure why) Either way, you'll have to calculate the
    # SVD for each G-matrix, which is 90% of the time taken up in inverting. 
    # One other option is to use the eigenvectors of the point cloud to fit the point cloud
    # Unfortunately, this seems equally sensitive to the issues above. (Why?? No idea...)
    # At any rate, this function is a bit slow, but it works consistently. No need to optimize yet
    
    #-- Basic code to fit a plane to a cloud of points via least-squares
    def solveForPlane(x1, x2, independent):
        """Fits a plane to a point cloud, where independent = ax1 + bx2 + C
        returns (a,b,c) and the residual"""
        from numpy.linalg import lstsq
        # See explanation in volume.world2model code if 
        # you're not familar with inverse theory.
        n = x1.size
        G = np.ones((n,3))
        G[:,0], G[:,1] = x1,x2
        m,resid,rank,s = lstsq(G,independent)
        return m, resid

    #-- The three possible formulations for fitting a plane
    def indX(x,y,z):
        # Fit x = ay + bz + c
        (a,b,c), resid = solveForPlane(y,z,x)
        m = (-1, a, b, c) # 0 = -x + ay + bz + c
        return m, resid 
    def indY(x,y,z):
        # Fit y = ax + bz + c
        (a,b,c), resid = solveForPlane(x,z,y)
        m = (a, -1, b, c) # 0 = ax - y + bz + c
        return m, resid
    def indZ(x,y,z):
        # Fit z = ax + by + c
        (a,b,c), resid = solveForPlane(x,y,z)
        m = (a, b, -1, c) # 0 = ax + by - z + c
        return m, resid

    #-- Process input arrays (size check, array conversion, and shift to 1D)
    x,y,z = [np.squeeze(np.asarray(item)) for item in [x,y,z]]
    if not (x.size == y.size == z.size):
        raise ValueError('Inputs must be the same size!')

    # By default, fit all 3 and choose the best
    if independent is None:
        m1, resid1 = indX(x,y,z)
        m2, resid2 = indY(x,y,z)
        m3, resid3 = indZ(x,y,z)
        #-- Return the model with the lowest residual
        models, resids = [m1, m2, m3], np.hstack((resid1, resid2, resid3))
        return  models[resids.argmin()]

    elif independent in [0, 'x']: # Fit x = ay + bz + c
        return indX(x,y,z)[0]
    elif independent in [1, 'y']: # Fit y = ax + bz + c
        return indY(x,y,z)[0]
    elif independent in [2, 'z']: # Fit z = ax + by + c
        return indZ(x,y,z)[0]
    
    else:
        raise ValueError('"independent" if specified, must be one of 0,1,2,"x","y","z".')

def normal2SD(x,y,z):
    """Converts a normal vector to a plane (given as x,y,z)
    to a strike and dip of the plane using the Right-Hand-Rule.
    Input:
        x: The x-component of the normal vector
        y: The y-component of the normal vector
        z: The z-component of the normal vector
    Output:
        strike: The strike of the plane, in degrees clockwise from north
        dip: The dip of the plane, in degrees downward from horizontal
    """
    from math import acos, asin, atan2,  sqrt, degrees

    # Due to geologic conventions, positive angles are downwards
    z = -z

    # First convert the normal vector to spherical coordinates
    #  (This is effectively a plunge/bearing of the normal vector)
    r = sqrt(x*x + y*y + z*z)
    plunge = degrees(asin(z/r))
    bearing = degrees(atan2(y, x))


    # Rotate bearing so that 0 is north instead of east
    bearing = 90-bearing
    if bearing<0: bearing += 360

    # If the plunge angle is upwards, get the opposite end of the line
    if plunge<0: 
        plunge = -plunge
        bearing -= 180
        if bearing<0: 
            bearing += 360

    # Now convert the plunge/bearing of the pole to the plane that it represents
    strike = bearing+90
    dip = 90-plunge
    if strike > 360: strike -= 360

    return strike, dip


