"""Various utilities using the geoprobe format api's
07/2009"""

import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

from volume import volume

def extractWindow(hor, vol, upper=0, lower=None, offset=0, region=None):
    """Extracts a window around a horizion out of a geoprobe volume
    Input:
        hor: a geoprobe horizion object or horizion filename
        vol: a geoprobe volume object or volume filename
        upper: (default, 0) upper window interval around horizion
        lower: (default, upper) lower window interval around horizion
        offset: (default, 0) amount (in voxels) to offset the horizion by along the Z-axis
        region: (default, full extent of horizion) sub-region to use instead of full extent
    Output:
        returns a numpy volume "flattened" along the horizion
    """
    
    # If the lower window isn't specified, assume it's equal to the upper window
    if lower == None: lower=upper

    # If filenames are input instead of volume/horizion objects, create the objects
    if type(hor) == type('String'):
        hor = horizion(hor)
    if type(vol) == type('String'):
        vol = volume(vol)

    # Set nodata points to be something useable
    depth = hor.grid
    depth[depth==hor.nodata] = lower
    depth = depth.T #Gah, changed the way hor.grid works... Should probably change it back

    # If extents are not set, use the full extent of the horizion (assumes horizion is smaller than volume)
    if region == None: 
        xmin, xmax, ymin, ymax = hor.xmin, hor.xmax, hor.ymin, hor.ymax
    else:
        xmin, xmax, ymin, ymax = region

    # Convert start and stop coords to volume coords
    xstart, ystart = xmin - vol.xmin, ymin - vol.ymin
    xstop, ystop = xmax - vol.xmin, ymax - vol.ymin
    # Select only the volume data within the region of interest
    data = vol.data[xstart:xstop, ystart:ystop, :]

    # Convert start and stop coords to horizion grid coords
    xstart, ystart = xmin - hor.xmin, ymin - hor.ymin
    xstop, ystop = xmax - hor.xmin, ymax - hor.ymin
    # Select only the volume data within the region of interest
    depth = depth[xstart:xstop, ystart:ystop]
    nx, ny = depth.shape
    
    # convert z coords of horizion to volume indexes
    depth = depth - vol.zmin
    depth = depth / np.abs(vol.dz)
    depth = np.array(depth, dtype=int)

    # Make indicies to extract a window around the horizion
    idxI = np.arange(nx)[:,np.newaxis,np.newaxis]
    idxJ = np.arange(ny)[np.newaxis,:,np.newaxis]
    idxK = depth[:,:,np.newaxis] + np.arange(-lower, upper) + offset 

    # Extract the subVolume
    subVolume = data[idxI,idxJ,idxK]
    subVolume = subVolume.astype(np.float)

    return subVolume

def coherence(data, window, rescale):
    """Calculates a coherence volume from a 3D numpy array. This method
    of calculating coherence implicitly assumes that the array consists
    of perodic data with a mean of 0."""
    if data.dtype in [np.uint8, np.int8]:
        data = data.astype(np.int16)
    ndimage.gaussian_filter(data**2, window, output=data, mode='constant', cval=0)
    return data

def wiggle(x, origin=0, posFill='black', negFill=None, lineColor='black', resampleRatio=10, rescale=False):
    """Plots a "wiggle" trace
    Input:
        x: input data (1D numpy array)
        origin: (default, 0) value to fill above or below (float)
        posFill: (default, black) color to fill positive wiggles with (string or None)
        negFill: (default, None) color to fill negative wiggles with (string or None)
        lineColor: (default, black) color of wiggle trace (string or None)
        resampleRatio: (default, 10) factor to resample traces by before plotting (1 = raw data) (float)
    Output:
        a matplotlib plot on the current axes
    """

    # Rescale so that x ranges from -1 to 1
    if rescale:
        x = x.astype(np.float)
        x -= x.min()
        x /= x.ptp()
        x *= 2
        x -= 1

    # Interpolate at 10x the previous density
    from scipy.signal import cspline1d, cspline1d_eval
    y = np.arange(0,x.size,1)
    newy = np.arange(0,x.size,1/float(resampleRatio))
    cj = cspline1d(x)
    interpX = cspline1d_eval(cj,newy) #,dx=1,x0=0
    if origin == None: origin = interpX.mean()

    # Plot
    if posFill is not None: 
        plt.fill_betweenx(newy,interpX,origin,where=interpX>origin,hold=True,facecolor=posFill)
    if negFill is not None:
        plt.fill_betweenx(newy,interpX,origin,where=interpX<origin,hold=True,facecolor=negFill)
    if lineColor is not None:
        plt.plot(interpX,newy,color=lineColor,hold=True)

def wiggles(grid, wiggleInterval=10, overlap=0.7, posFill='black', negFill=None, lineColor='black', 
        rescale=True, extent=None):
    """Plots a series of "wiggle" traces based on a grid
    Input:
        x: input data (2D numpy array)
        wiggleInterval: (default, 10) Plot 'wiggles' every wiggleInterval traces
        overlap: (default, 0.7) amount to overlap 'wiggles' by (1.0 = scaled to wiggleInterval)
        posFill: (default, black) color to fill positive wiggles with (string or None)
        negFill: (default, None) color to fill negative wiggles with (string or None)
        lineColor: (default, black) color of wiggle trace (string or None)
        resampleRatio: (default, 10) factor to resample traces by before plotting (1 = raw data) (float)
    Output:
        a matplotlib plot on the current axes
    """
    # Rescale so that the grid ranges from -1 to 1
    if rescale:
        grid = grid.astype(np.float)
        grid -= grid.min()
        grid /= grid.ptp()
        grid *= 2
        grid -= 1

    # TODO: Finish this!! add support for plotting with proper coords...
    if extent is None:
        xmin, ymin = 0, 0
        xmax, ymax = grid.shape

    ny,nx = grid.shape
    for i in range(wiggleInterval//2, nx, wiggleInterval):
        x = overlap * wiggleInterval/2 * grid[:,i]
        wiggle(x+i,i,posFill,negFill,lineColor)


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
            vol = volume(vol)
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
    to a strike and dip of the plane using the Right-Hand-Rule"""
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


def roseDiagram(data, title='North', nbins=30, bidirectional=True):
    data = np.array(data)
    n = data.size

    if bidirectional:
        # Rather than doing some sort of fancy binning, just
        #  "double" the data with the complimentary end (+180)
        data = np.append(data, data+180)
        data[data>360] -= 360

    # Rotate the data so that north will plot at the top (90deg, in polar space)
    data = 90-data
    data[data<0] += 360

    # Make a figure with polar axes
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.7], polar=True, axisbg='#d5de9c')

    # Change the labeling so that north is at the top
    plt.thetagrids(range(0,360,45), ['90','45','0','315','270','225','180','135'])

    # Plot a histogram on the polar axes
    data = np.radians(data)
    plt.hist(data, bins=nbins, axes=ax)
    plt.title(title + '\nn=%i'%n)


