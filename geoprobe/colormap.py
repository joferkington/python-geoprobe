import six
import numpy as np

class colormap(object):
    """Reads a Geoprobe formatted colormap"""
    num_colors = 256
    def __init__(self, filename):
        self.filename = filename
        self._parse_infile()

    def _parse_infile(self):
        infile = open(self.filename, 'r')
        header = next(infile)
        if header.startswith('#'):
            _ = next(infile)
        self.num_keys = int(next(infile).strip())
        keys = []
        for i in range(self.num_keys):
            keys.append(next(infile).strip().split())
        self.keys = np.array(keys, dtype=np.float)
        num_colors = int(next(infile).strip())
        colors = []
        for i in range(num_colors):
            colors.append(next(infile).strip().split())
        self.lut = np.array(colors, dtype=np.float)
        dtype = {'names':['red', 'green', 'blue', 'alpha', 'keys'],
                 'formats':5 * [np.float]}
        self.lut = self.lut.view(dtype)

    @property
    def as_matplotlib(self):
        from matplotlib.colors import LinearSegmentedColormap
        cdict = dict(red=[], green=[], blue=[])

        # Make sure that there is a key at 0.0 and 1.0
        keys = self.keys.tolist()
        if keys[0][0] != 0:
            keys = [[0.0] + keys[0][1:]] + keys
        if keys[-1][0] != 1.0:
            keys.append([1.0] + keys[-1][1:])

        for stop_value, red, green, blue, alpha, in keys:
            for name, val in zip(['red', 'green', 'blue'], [red, green, blue]):
                cdict[name].append([stop_value, val, val])
        return LinearSegmentedColormap(self.filename, cdict, self.num_colors)

    @property
    def as_pil(self):
        return list((255 * self.lut.view(np.float)[:,:3]).astype(np.int).flat)

    @property
    def as_pil_rgba(self):
        return list((255 * self.lut.view(np.float)[:,:4]).astype(np.int).flat)







