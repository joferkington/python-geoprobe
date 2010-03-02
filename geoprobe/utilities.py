"""Various utilities using the geoprobe format api's
07/2009"""

import numpy as np

from volume import volume
from horizon import horizon

#Functions that were formerly in this file...
from common import array2geotiff, points2strikeDip, fitPlane, normal2SD

def extractWindow(hor, vol, upper=0, lower=None, offset=0, region=None, masked=False):
    """Extracts a window around a horizion out of a geoprobe volume
    Input:
        hor: a geoprobe horizion object or horizion filename
        vol: a geoprobe volume object or volume filename
        upper: (default, 0) upper window interval around horizion
        lower: (default, upper) lower window interval around horizion
        offset: (default, 0) amount (in voxels) to offset the horizion by along the Z-axis
        region: (default, overlap between horizion and volume) sub-region to use instead of 
                full extent. Must be a 4-tuple of (xmin, xmax, ymin, ymax)
        masked: (default, False) if True, return a masked array where nodata values in the
                horizon are masked. Otherwise, return an array where the nodata values are
                filled with 0.
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

    # If the lower window isn't specified, assume it's equal to the upper window
    if lower == None: lower=upper

    # If filenames are input instead of volume/horizion objects, create the objects
    if type(hor) == type('String'):
        hor = horizon(hor)
    if type(vol) == type('String'):
        vol = volume(vol)

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
            raise ValueError('Specified region does not overlap with horizon and volume')
        elif len(extents) != 4:
            raise ValueError('"extents" must be a 4-tuple of (xmin, xmax, ymin, ymax)')

    xmin, xmax, ymin, ymax = extents

    # Convert extents to volume indicies and select subset of the volume 
    xstart, ystart = xmin - vol.xmin, ymin - vol.ymin
    xstop, ystop = xmax - vol.xmin, ymax - vol.ymin
    data = vol.data[xstart:xstop, ystart:ystop, :]

    # Convert extents to horizion grid indicies and select subset of the horizion 
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
    # Not creating a masked array here due to speed problems when iterating through ma's
    subVolume = np.zeros((nx,ny,window_size), dtype=np.uint8)

    # Using fancy indexing to do this uses tons of memory...
    # As it turns out, simple iteration is much, much more memory efficient, and almost as fast
    mask = depth.mask      # Need to preserve the mask for later
    depth = depth.filled() # Iterating through masked arrays is much slower, apparently
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

    # If masked is True (input option), return a masked array
    if masked:
        nx,ny,nz = subVolume.shape
        mask = mask.reshape((nx,ny,1))
        mask = np.tile(mask, (1,1,nz))
        subVolume = np.ma.array(subVolume, mask=mask)

    # If upper==lower==0, (default) subVolume will be (nx,ny,1), so return 2D array instead
    subVolume = subVolume.squeeze()

    return subVolume


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


