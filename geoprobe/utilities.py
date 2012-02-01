"""Various utilities using the geoprobe format api's
07/2009"""

import numpy as np

import volume
import horizon

def extractWindow(hor, vol, upper=0, lower=None, offset=0, 
                region=None, masked=False):
    """Extracts a window around a horizion out of a geoprobe volume
    Input:
        hor: a geoprobe horizion object or horizion filename
        vol: a geoprobe volume object or volume filename
        upper: (default, 0) upper window interval around horizion in voxels 
        lower: (default, upper) lower window interval around horizion in voxels 
        offset: (default, 0) amount (in voxels) to offset the horizion by along
                the Z-axis
        region: (default, overlap between horizion and volume) sub-region to 
                use instead of full extent. Must be a 4-tuple of (xmin, xmax,
                ymin, ymax)
        masked: (default, False) if True, return a masked array where nodata 
                values in the horizon are masked. Otherwise, return an array
                where the nodata values are filled with 0.
    Output:
        returns a numpy volume "flattened" along the horizion
    """
    def bbox_overlap(bbox1, bbox2):
        """
        Input: 
            bbox1: 4-tuple of (xmin, xmax, ymin, ymax) for the first region
            bbox2: 4-tuple of (xmin, xmax, ymin, ymax) for the second region
        Output:
            Returns a 4-tuple of (xmin, xmax, ymin, ymax) for overlap between 
            the two input bounding-box regions
        """
        def intersects(bbox1, bbox2):
            # Check for intersection
            xmin1, xmax1, ymin1, ymax1 = bbox1
            xmin2, xmax2, ymin2, ymax2 = bbox2
            xdist = abs( (xmin1 + xmax1) / 2.0 - (xmin2 + xmax2) / 2.0 )
            ydist = abs( (ymin1 + ymax1) / 2.0 - (ymin2 + ymax2) / 2.0 )
            xwidth = (xmax1 - xmin1 + xmax2 - xmin2) / 2.0
            ywidth = (ymax1 - ymin1 + ymax2 - ymin2) / 2.0
            return (xdist <= xwidth) and (ydist <= ywidth)

        if not intersects(bbox1, bbox2):
            return None

        output = []
        for i, comparison in enumerate(zip(bbox1, bbox2)):
            # mininum x or y coordinate
            if i % 2 == 0:
                output.append(max(comparison))
            # maximum x or y coordinate
            elif i % 2 == 1:
                output.append(min(comparison))
        return output

    # If the lower window isn't specified, assume it's equal to the 
    # upper window
    if lower == None: lower=upper

    # If filenames are input instead of volume/horizion objects, create 
    # the objects
    if type(hor) == type('String'):
        hor = horizon.horizon(hor)
    if type(vol) == type('String'):
        vol = volume.volume(vol)

    #Gah, changed the way hor.grid works... Should probably change it back
    depth = hor.grid.T 

    # Find the overlap between the horizon and volume
    vol_extents = [vol.xmin, vol.xmax, vol.ymin, vol.ymax]
    hor_extents = [hor.xmin, hor.xmax, hor.ymin, hor.ymax]
    extents = bbox_overlap(hor_extents, vol_extents)

    # Raise ValueError if horizon and volume do not intersect
    if extents is None:
        raise ValueError('Input horizon and volume do not intersect!')

    # Find the overlap between the (optional) subregion current extent
    if region is not None: 
        extents = bbox_overlap(extents, region)
        if extents is None:
            raise ValueError('Specified region does not overlap with'\
                             ' horizon and volume')
        elif len(extents) != 4:
            raise ValueError('"extents" must be a 4-tuple of'\
                             ' (xmin, xmax, ymin, ymax)')

    xmin, xmax, ymin, ymax = extents

    # Convert extents to volume indicies and select subset of the volume 
    xstart, ystart = xmin - vol.xmin, ymin - vol.ymin
    xstop, ystop = xmax - vol.xmin, ymax - vol.ymin
    data = vol.data[xstart:xstop, ystart:ystop, :]

    # Convert extents to horizion grid indicies and select 
    # subset of the horizion 
    xstart, ystart = xmin - hor.xmin, ymin - hor.ymin
    xstop, ystop = xmax - hor.xmin, ymax - hor.ymin
    depth = depth[xstart:xstop, ystart:ystop]

    nx,ny,nz = data.shape
    
    # convert z coords of horizion to volume indexes
    depth -= vol.zmin
    depth /= abs(vol.dz)
    depth = depth.astype(np.int)

    # Initalize the output array
    window_size = upper + lower + 1
    # Not creating a masked array here due to speed problems when 
    # iterating through ma's
    subVolume = np.zeros((nx,ny,window_size), dtype=np.uint8)

    # Using fancy indexing to do this uses tons of memory...
    # As it turns out, simple iteration is much, much more memory 
    # efficient, and almost as fast
    mask = depth.mask      # Need to preserve the mask for later
    depth = depth.filled() # Iterating through masked arrays is much slower...
    for i in xrange(nx):
        for j in xrange(ny):
            if depth[i,j] != hor.nodata:
                # Where are we in data indicies
                z = depth[i,j] + offset
                top = z - upper 
                bottom = z + lower + 1

                # Be careful to extract the right vertical region in cases
                # where the window goes outside the data bounds (find region of
                # overlap)
                data_top = max([top, 0]) 
                data_bottom = min([bottom, nz]) 
                window_top = max([0, window_size - bottom])
                window_bottom = min([window_size, nz - top])

                # Extract the window out of data and store it in subVolume
                subVolume[i, j, window_top : window_bottom] \
                                = data[i, j, data_top : data_bottom]

    # If masked is True (input option), return a masked array
    if masked:
        nx,ny,nz = subVolume.shape
        mask = mask.reshape((nx,ny,1))
        mask = np.tile(mask, (1,1,nz))
        subVolume = np.ma.array(subVolume, mask=mask)

    # If upper==lower==0, (default) subVolume will be (nx,ny,1), so 
    # return 2D array instead
    subVolume = subVolume.squeeze()

    return subVolume

def wiggle(x, origin=0, posFill='black', negFill=None, lineColor='black', 
        resampleRatio=10, rescale=False, ymin=0, ymax=None, ax=None):
    """Plots a "wiggle" trace
    Input:
        x: input data (1D numpy array)
        origin: (default, 0) value to fill above or below (float)
        posFill: (default, black) color to fill positive wiggles with (string 
            or None)
        negFill: (default, None) color to fill negative wiggles with (string 
            or None)
        lineColor: (default, black) color of wiggle trace (string or None)
        resampleRatio: (default, 10) factor to resample traces by before 
            plotting (1 = raw data) (float)
        rescale: (default, False) If True, rescale "x" to be between -1 and 1
        ymin: (default, 0) The minimum y to use for plotting
        ymax: (default, len(x)) The maximum y to use for plotting
        ax: (default, current axis) The matplotlib axis to plot onto
    Output:
        a matplotlib plot on the current axes
    """
    from matplotlib import pyplot as plt
    from scipy.signal import cspline1d, cspline1d_eval

    if ymax is None:
        ymax = x.size

    # Rescale so that x ranges from -1 to 1
    if rescale:
        x = x.astype(np.float)
        x -= x.min()
        x /= x.ptp()
        x *= 2
        x -= 1

    # Interpolate at resampleRatio x the previous density
    y = np.linspace(0, x.size, x.size)
    interp_y = np.linspace(0, x.size, x.size * resampleRatio)
    cj = cspline1d(x)
    interpX = cspline1d_eval(cj,interp_y) #,dx=1,x0=0
    newy = np.linspace(ymax, ymin, interp_y.size)
    if origin == None: 
        origin = interpX.mean()

    # Plot
    if ax is None:
        ax = plt.gca()
        plt.hold(True)
    if posFill is not None: 
        ax.fill_betweenx(newy, interpX, origin,
                where=interpX > origin,
                facecolor=posFill)
    if negFill is not None:
        ax.fill_betweenx(newy, interpX, origin,
                where=interpX < origin,
                facecolor=negFill)
    if lineColor is not None:
        ax.plot(interpX, newy, color=lineColor)

def wiggles(grid, wiggleInterval=10, overlap=0.7, posFill='black', 
        negFill=None, lineColor='black', rescale=True, extent=None, ax=None):
    """Plots a series of "wiggle" traces based on a grid
    Input:
        x: input data (2D numpy array)
        wiggleInterval: (default, 10) Plot 'wiggles' every wiggleInterval 
            traces
        overlap: (default, 0.7) amount to overlap 'wiggles' by (1.0 = scaled 
            to wiggleInterval)
        posFill: (default, black) color to fill positive wiggles with (string 
            or None)
        negFill: (default, None) color to fill negative wiggles with (string 
            or None)
        lineColor: (default, black) color of wiggle trace (string or None)
        resampleRatio: (default, 10) factor to resample traces by before 
            plotting (1 = raw data) (float)
        extent: (default, (0, nx, 0, ny)) The extent to use for the plot. 
            A 4-tuple of (xmin, xmax, ymin, ymax)
        ax: (default, current axis) The matplotlib axis to plot onto.
    Output:
        a matplotlib plot on the current axes
    """
    from matplotlib import pyplot as plt

    # Rescale so that the grid ranges from -1 to 1
    if rescale:
        grid = grid.astype(np.float)
        grid -= grid.min()
        grid /= grid.ptp()
        grid *= 2
        grid -= 1

    if extent is None:
        xmin, ymin = 0, grid.shape[0]
        ymax, xmax = 0, grid.shape[1]
    else:
        xmin, xmax, ymin, ymax = extent

    ny,nx = grid.shape
    x_loc = np.linspace(xmin, xmax, nx)
    for i in range(wiggleInterval//2, nx, wiggleInterval):
        x = overlap * (wiggleInterval / 2.0) \
                    * (x_loc[1] - x_loc[0]) \
                    * grid[:,i]
        wiggle(x + x_loc[i], origin=x_loc[i], 
                posFill=posFill, negFill=negFill, 
                lineColor=lineColor, ymin=ymin, 
                ymax=ymax, ax=ax)


def roseDiagram(data, nbins=30, bidirectional=True, title='North'):
    """Plots a circular histogram or "rose diagram"
    Input: 
        data: A list or 1D array of orientation data that the histogram 
            will be made from. The data should be an in degrees clockwise
            from north.
        nbins (default: 30): The number of bins in the histogram
        bidirectional (default: True): Whether or not to treat the input data
            as bi-directional. (i.e. if True, the rose diagram will be 
            symmetric)
        title (default: 'North'): The title of the plot
    """
    # TODO: This needs to pass kwargs on to the plotting routines
    # TODO: (or just remove this function entirely? It shouldn't 
    # TODO: really be here)
    from matplotlib import pyplot as plt
    data = np.asarray(data)
    n = data.size

    if bidirectional:
        # Rather than doing some sort of fancy binning, just
        #  "double" the data with the complimentary end (+180)
        data = np.append(data, data+180)
        data[data>360] -= 360

    # Rotate the data so that north will plot at the top 
    # (90deg, in polar space)
    data = 90-data
    data[data<0] += 360

    # Make a figure with polar axes
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.7], polar=True, axisbg='#d5de9c')

    # Change the labeling so that north is at the top
    plt.thetagrids(range(0,360,45), 
            ['90','45','0','315','270','225','180','135'])

    # Plot a histogram on the polar axes
    data = np.radians(data)
    plt.hist(data, bins=nbins, axes=ax)
    plt.title(title + '\nn=%i'%n)


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
        except AttributeError:
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
        filename += '.tif'

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
    from math import asin, atan2,  sqrt, degrees

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


