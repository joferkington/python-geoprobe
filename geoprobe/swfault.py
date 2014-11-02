import numpy as np
import xml.etree.ElementTree as et
import xml.dom.minidom as minidom
import utilities
# Note: matplotlib.delaunay is required for interpolation and triangulation.
# If it isn't available

class swfault(object):
    def __init__(self, arg):
        if isinstance(arg, basestring):
            # Assume it's a filename
            self._read(arg)
        else:
            # Assume it's a sequence of segments
            self._init_from_data(arg)
        self.dx = 1.0
        self.dy = 1.0

    def _read(self, filename):
        with open(filename, 'r') as infile:
            reader = SwfaultXMLReader(infile)
            self.segments = reader.segments
            self.color = reader.color
            self.name = reader.name
            self.linewidth = reader.linewidth
            self.seisworks_ids = reader.seisworks_ids

    def write(self, filename):
        with open(filename, 'w') as outfile:
            writer = SwfaultXMLWriter(self)
            writer.write(outfile)

    def _init_from_data(self, segments, name='Undefined'):
        self.segments = segments
        self.name = name
        self.color = (1, 0, 0, 1)
        self.linewidth = 3.0
        self.seisworks_ids = ['2147483647'] * len(segments)

    def _get_coord(self, axis):
        try:
            return self.xyz[:,axis]
        except IndexError:
            return np.array([])
    def _set_coord(self, value, axis):
        self.xyz[:,axis] = value
    x = property(lambda self: self._get_coord(0),
                 lambda self, value: self._set_coord(value, 0))
    y = property(lambda self: self._get_coord(1),
                 lambda self, value: self._set_coord(value, 1))
    z = property(lambda self: self._get_coord(2),
                 lambda self, value: self._set_coord(value, 2))

    #-- Bounding values... ----------------------------------------------------
    xmin = property(lambda self: self.x.min(),
            doc='Mininum X-coordinate (in inline/crossline)')
    ymin = property(lambda self: self.y.min(),
            doc='Mininum Y-coordinate (in inline/crossline)')
    zmin = property(lambda self: self.z.min(),
            doc='Mininum Z-coordinate (in inline/crossline)')
    xmax = property(lambda self: self.x.max(),
            doc='Maximum X-coordinate (in inline/crossline)')
    ymax = property(lambda self: self.y.max(),
            doc='Maximum Y-coordinate (in inline/crossline)')
    zmax = property(lambda self: self.z.max(),
            doc='Maximum Z-coordinate (in inline/crossline)')

    @property
    def grid_extents(self):
        try:
            return self._grid_extents
        except AttributeError:
            return (self.xmin, self.xmax, self.ymin, self.ymax)
    @grid_extents.setter
    def grid_extents(self, value):
        xmin, xmax, ymin, ymax = value
        self._grid_extents = (xmin, xmax, ymin, ymax)

    @property
    def segments(self):
        return [self.xyz[item] for item in self._indices]
    @segments.setter
    def segments(self, value):
        def sequence(item):
            length = len(item)
            sli = slice(sequence.i, sequence.i + length)
            sequence.i += length
            return sli
        sequence.i = 0
        self.xyz = np.array([item for segment in value for item in segment])
        self._indices = [sequence(item) for item in value]

    @property
    def tri(self):
        try:
            from matplotlib.delaunay import Triangulation
        except ImportError:
            raise ImportError('Maplotlib is required for geoprobe.swfault '
                              'triangulation')
        try:
            return self._tri
        except AttributeError:
            x, y, z = self._internal_xyz.T
            self._tri = Triangulation(x, y)
            self._tri.x = self.x
            self._tri.y = self.y
            return self._tri

    @property
    def interp(self):
        try:
            from matplotlib.delaunay import LinearInterpolator
        except ImportError:
            raise ImportError('Maplotlib is required for geoprobe.swfault '
                              'interpolation')
        try:
            return self._interp
        except AttributeError:
            self._interp = LinearInterpolator(self.tri, self.z)
            return self._interp

    @property
    def grid(self):
        try:
            return self._grid
        except AttributeError:
            xmin, xmax, ymin, ymax = self.grid_extents
            #nx, ny = int((xmax-xmin) / self.dx), int((ymax-ymin) / self.dy)
            nx, ny = int((xmax-xmin) / self.dx)+1, int((ymax-ymin) / self.dy)+1
            self._grid = self.interp[ymin:ymax:ny*1j, xmin:xmax:nx*1j]
            return self._grid

    @property
    def _internal_xyz(self):
        try:
            return self._internal_xyz_data
        except:
            vecs, vals = utilities.principal_axes(self.x, self.y, self.z, True)
            rotated = self.xyz.dot(vecs)
            rotated -= rotated.mean(axis=0)
            rotated /= np.sqrt(vals)
            self._internal_xyz_data = rotated
            return rotated

    @property
    def _outline_order(self):
        rotated_segments = [self._internal_xyz[item] for item in self._indices]
        xyz = np.array(self._segment_endpoints(rotated_segments))
        x, y, z = xyz.T

        theta = np.arctan2(y, x)
        order = theta.argsort()
        return order

    def _segment_endpoints(self, segments):
        return [seg[0] for seg in segments] + [seg[-1] for seg in segments]

    @property
    def outline(self):
        try:
            return self._outline
        except:
            xyz = np.array(self._segment_endpoints(self.segments))
            outline = xyz[self._outline_order]
            return np.squeeze(outline)

    @property
    def _rotated_outline(self):
        rotated_segments = [self._internal_xyz[item] for item in self._indices]
        xyz = np.array(self._segment_endpoints(rotated_segments))
        return np.squeeze(xyz[self._outline_order])

class SwfaultXMLReader(object):
    def __init__(self, f):
        data = f.read()
        try:
            self.tree = et.fromstring(data)
        except et.ParseError:
            # Stupid geoprobe produces invalid xml...
            self.tree = et.fromstring(data + '</boost_serialization>\n')
        self._parse_xml_tree()

    def _parse_xml_tree(self):
        def get_xyz(element):
            coords = [element.findall('item/'+key) for key in ['x','y','z']]
            xyz = zip(*[[float(coord.text) for coord in item] for item in coords])
            return xyz

        fault = self.tree.find('swFault')
        self.name = fault.find('name').text
        color_fields = ['lineColor_' + item for item in ('r', 'g', 'b', 'a')]
        self.color = [float(fault.find(item).text) for item in color_fields]
        self.linewidth = float(fault.find('lineWidth').text)

        self.segments = [get_xyz(elem) for elem in
                            fault.findall('segments/item/Project_Coordinates')]

        self.seisworks_ids = [item.find('SWSegID').text for item in
                                fault.findall('segments/item')]

class SwfaultXMLWriter(object):
    def __init__(self, parent):
        self.parent = parent

        self.setup_document()
        self.write_metadata()
        self.write_all_segments()

    def write(self, f):
        f.write(self.doc.toprettyxml(indent='    '))

    def setup_document(self):
        imp = minidom.getDOMImplementation()
        doctype = imp.createDocumentType('boost_serialization', '', '')
        self.doc = imp.createDocument(None, 'boost_serialization', doctype)
        self.doc.documentElement.setAttribute('signature', 'serialization::archive')
        self.doc.documentElement.setAttribute('version', '3')

        self.root = self.add_node(self.doc.documentElement, 'swFault', None,
                dict(class_id='0', tracking_level='0', version='3'))

    def write_metadata(self):
        self.add_node(self.root, 'name', self.parent.name)
        color_fields = ['lineColor_' + item for item in ('r', 'g', 'b', 'a')]
        for color, field in zip(self.parent.color, color_fields):
            self.add_node(self.root, field, repr(color))
        self.add_node(self.root, 'lineWidth', repr(self.parent.linewidth))

    def write_all_segments(self):
        segments = self.add_node(self.root, 'segments', None,
                dict(class_id='1', tracking_level='0', version='0'))
        self.add_node(segments, 'count', repr(len(self.parent.segments)))

        items = zip(self.parent.segments, self.parent.seisworks_ids)
        for i, (seg, sw_id) in enumerate(items):
            self.write_segment(i, seg, sw_id, segments)

    def write_segment(self, i, seg, sw_id, segment_node):
        # First items are different than the rest...
        if i == 0:
            item_attrs = dict(class_id='2', tracking_level='1', version='2', object_id='_0')
            coords_attrs = dict(class_id='3', tracking_level='0', version='0')
        else:
            item_attrs = dict(class_id_reference='2', object_id='_{}'.format(i))
            coords_attrs = {}

        item = self.add_node(segment_node, 'item', None, item_attrs)
        self.add_node(item, 'SWSegID', sw_id)
        coords = self.add_node(item, 'Project_Coordinates', None, coords_attrs)
        self.add_node(coords, 'count', repr(len(seg)))

        firstrun = True
        for x, y, z in seg:
            if i == 0 and firstrun:
                point_attrs = dict(class_id='4', tracking_level='0', version='0')
                firstrun = False
            else:
                point_attrs = {}
            point = self.add_node(coords, 'item', None, point_attrs)
            for name, value in zip(['x', 'y', 'z'], [x, y, z]):
                self.add_node(point, name, repr(value))

    def add_node(self, root, name, text=None, attributes=None):
        """Create and append node "name" to the node "root"."""
        node = self.doc.createElement(name)
        if text is not None:
            text = self.doc.createTextNode(text)
            node.appendChild(text)
        if attributes is not None:
            for key, value in attributes.iteritems():
                node.setAttribute(key, value)
        root.appendChild(node)
        return node
