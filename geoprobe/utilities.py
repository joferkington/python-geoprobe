"""Various utilities using the geoprobe format api's
07/2009"""

import numpy as np

from volume import volume
from horizon import horizon

#Functions that were formerly in this file...
from common import array2geotiff, points2strikeDip, fitPlane, normal2SD

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
        hor = horizon(hor)
    if type(vol) == type('String'):
        vol = volume(vol)

    #Gah, changed the way hor.grid works... Should probably change it back
    depth = hor.grid.filled().T 

    # If extents are not set, use the full extent of the horizion (assumes horizion is smaller than volume)
    if region == None: 
        xmin, xmax, ymin, ymax = hor.xmin, hor.xmax, hor.ymin, hor.ymax
    else:
        xmin, xmax, ymin, ymax = region

    #-- Select region of overlap between volume and horizon
    # Convert region to volume indicies
    xstart, ystart = xmin - vol.xmin, ymin - vol.ymin
    xstop, ystop = xmax - vol.xmin, ymax - vol.ymin
    # Select only the volume data within the region of interest
    data = vol.data[xstart:xstop, ystart:ystop, :]

    # Convert region to horizion grid indicies
    xstart, ystart = xmin - hor.xmin, ymin - hor.ymin
    xstop, ystop = xmax - hor.xmin, ymax - hor.ymin
    # Select only the volume data within the region of interest
    depth = depth[xstart:xstop, ystart:ystop]

    nx,ny,nz = data.shape
    
    # convert z coords of horizion to volume indexes
    nodata = depth != hor.nodata
    depth[nodata] -= vol.zmin
    depth[nodata] /= abs(vol.dz)
    depth = depth.astype(np.int)

    # Using fancy indexing to do this uses tons of memory...
    # As it turns out, simple iteration is much, much more memory efficient, and almost as fast
    window_size = upper + lower + 1
    subVolume = np.zeros((nx,ny,window_size), dtype=np.uint8)
    for i in xrange(nx):
        for j in xrange(ny):
            if depth[i,j] != hor.nodata:
                # Where are we in data indicies
                z = depth[i,j] + offset
                top = z - upper 
                bottom = z + lower + 1

                # Be careful to extract the right vertical region in cases where the 
                # window goes outside the data bounds (find region of overlap)
                data_top = max([top, 0]) 
                data_bottom = min([bottom, nz]) 
                window_top = max([0, window_size - bottom])
                window_bottom = min([window_size, nz - top])

                # Extract the window out of data and store it in subVolume
                subVolume[i,j,window_top:window_bottom] = data[i,j,data_top:data_bottom]

    return subVolume.squeeze()


def coherence(data, window=(0.3, 0.3, 2.0)):
    """Calculates a coherence volume from a 3D numpy array using a
    gaussian-shaped moving window. This method of calculating coherence 
    implicitly assumes that the array consistsof perodic data with a 
    mean of 0. If the data consists of 8-bit values, it will be 
    converted to 16-bit to avoid overflow (and centered on 0, if the 
    data are unsigned integers)
    Input:
        data: Input data (a 3d numpy array)
        window: A tuple of (xlength, ylength, zlength) describing the
                size of the gaussian window. Fractional values _are_
                allowed.
    """
    from scipy import ndimage
    # To avoid overflow, if we're dealing with an 8-bit array, convert it to 16-bit
    if data.dtype == np.uint8:
        data = data.astype(np.int16) - 127
    elif data.dtype == np.int8:
        data = data.astype(np.int16)
    ndimage.gaussian_filter(data, window, output=data, mode='constant', cval=0)
    data = np.sqrt(data)
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
    from matplotlib import pyplot as plt
    from scipy.signal import cspline1d, cspline1d_eval

    # Rescale so that x ranges from -1 to 1
    if rescale:
        x = x.astype(np.float)
        x -= x.min()
        x /= x.ptp()
        x *= 2
        x -= 1

    # Interpolate at 10x the previous density
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
    from matplotlib import pyplot as plt

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
    # TODO: (also, take the title parameter out) 
    # TODO: (or just remove this function entirely? It shouldn't really be here)
    from matplotlib import pyplot as plt
    data = np.asarray(data)
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


