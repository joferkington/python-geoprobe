import os
import numpy as np
import geoprobe

datadir = os.path.join(os.path.dirname(__file__), '..', 'examples', 'data/')
faultdir = os.path.join(datadir, 'swFaults')

class TestBasic:
    def setup(self):
        self.normal = geoprobe.swfault(faultdir + '/example_normal.swf')
        self.strike_slip = geoprobe.swfault(faultdir + '/example_ss.swf')
        self.empty = geoprobe.swfault(faultdir + '/empty.swf')
    
    def test_empty(self):
        x, y, z = self.empty.x, self.empty.y, self.empty.z
        assert np.allclose(x, [])
        assert np.allclose(y, [])
        assert np.allclose(z, [])

