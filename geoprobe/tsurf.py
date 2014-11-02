class tsurf(object):
    default_name = 'Undefined'
    default_color = (0, 1, 1, 1.0)
    def __init__(self, *args, **kwargs):
        """Accepts either a single filename or 4 arguments: x, y, z, triangles.
        keyword argumets are: "color" and "name"

        If a filename is given, the tsurf is read from the file.

        Otherwise:
        x, y, z are sequences of the x, y, and z coordinates of the vertices.
        triangles is a sequence of the indicies of the coord arrays making up
        each triangle in the mesh. E.g. [[0, 1, 2], [2, 1, 3], ...]"""
        if len(args) == 1:
            self._read_tsurf(args[0])
        elif len(args) == 4:
            self._init_from_xyz(*args)
        else:
            raise ValueError('Invalid input arguments')
        color = kwargs.get('color', None)
        name = kwargs.get('name', None)
        if color is not None:
            self.color = color
        if name is not None:
            self.name = name
        self.header['name'] = self.name
        self.header['color'] = self.color

    def _read_tsurf(self, filename):
        with open(filename, 'r') as infile:
            firstline = next(infile).strip()
            if not firstline.startswith('GOCAD TSurf'):
                raise IOError('This is not a valid TSurf file!')

            # Parse Header
            self.header = {}
            line = next(infile).strip()
            if line.startswith('HEADER'):
                line = next(infile).strip()
                while '}' not in line:
                    key, value = line.split(':')
                    self.header[key.lstrip('*')] = value
                    line = next(infile).strip()
            self.name = self.header.get('name', filename)
            try:
                self.color = [float(item) for item in self.header['color'].split()]
                self.color = tuple(self.color)
            except KeyError:
                self.color = self.default_color

            # Read vertices and triangles
            if not next(infile).startswith('TFACE'):
                raise IOError('Only "TFACE" format TSurf files are supported')
            self.vertices, self.triangles = [], []
            for line in infile:
                line = line.strip().split()
                if line[0] == 'VRTX':
                    self.vertices.append([float(item) for item in line[2:]])
                elif line[0] == 'TRGL':
                    self.triangles.append([int(item)-1 for item in line[1:]])
        self.x, self.y, self.z = zip(*self.vertices)

    def _init_from_xyz(self, x, y, z, triangles):
        self.vertices = zip(x, y, z)
        self.x, self.y, self.z = x, y, z
        self.triangles = triangles
        self.color = self.default_color
        self.name = self.default_name
        self.header = {'moveAs':'2', 'drawAs':'2', 'line':'3',
                'clip':'0', 'intersect':'0', 'intercolor':' 1 0 0 1'}

    def write(self, outname):
        with open(outname, 'w') as outfile:
            # Write Header...
            outfile.write('GOCAD TSurf 1\n')
            outfile.write('HEADER {\n')
            """
            for key in ['name', 'color', 'moveAs', 'drawAs', 'line', 'clip',
                        'intersect', 'intercolor']:
                value = self.header[key]
            """
            for key, value in self.header.iteritems():
                if not isinstance(value, basestring):
                    try:
                        value = ' '.join(repr(item) for item in value)
                    except TypeError:
                        value = repr(item)
                outfile.write('*{}:{}\n'.format(key, value))
            outfile.write('}\n')

            # Write data...
            outfile.write('TFACE\n')
            for i, (x, y, z) in enumerate(self.vertices, start=1):
                template = '\t'.join(['VRTX {}'] + 3*['{: >9.3f}']) + '\n'
                outfile.write(template.format(i, x, y, z))
            for a, b, c in self.triangles:
                outfile.write('TRGL {} {} {}\n'.format(a+1, b+1, c+1))
            outfile.write('END\n')


