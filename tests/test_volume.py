import os
import numpy as np
import pytest
import geoprobe


class TestBase:
    def setup_method(self, method):
        self.basedir = os.path.dirname(geoprobe.__file__)
        self.datadir = os.path.join(self.basedir, '..', 'examples', 'data')
        volname = os.path.join(self.datadir, 'Volumes', 'example.vol')
        self.vol = geoprobe.volume(volname)

class TestExtents(TestBase):
    def test_bounds_sanity(self):
        assert self.vol.xmin <= self.vol.xmax
        assert self.vol.ymin <= self.vol.ymax
        assert self.vol.zmin <= self.vol.zmax

    def test_bounds(self):
        assert self.vol.xmin == 2168.0
        assert self.vol.xmax == 2537.0
        assert self.vol.ymin == 4900.0
        assert self.vol.ymax == 5154.0
        assert self.vol.zmin == 2720.0
        assert self.vol.zmax == 3990.0

    def test_dx_dy_dz(self):
        """Ensure that bounds stay min <= max even when d's are negative."""
        for item in ['dx', 'dy', 'dz']:
            setattr(self.vol, item, -getattr(self.vol, item))
            self.test_bounds_sanity()
            setattr(self.vol, item, -getattr(self.vol, item))


class TestConversions(TestBase):
    def test_model2index_mins(self):
        assert self.vol.model2index(self.vol.xmin, axis='x') == 0
        assert self.vol.model2index(self.vol.ymin, axis='y') == 0
        assert self.vol.model2index(self.vol.zmin, axis='z') == 0

    def test_model2index_maxs(self):
        assert self.vol.model2index(self.vol.xmax, axis='x') == self.vol.nx - 1
        assert self.vol.model2index(self.vol.ymax, axis='y') == self.vol.ny - 1
        assert self.vol.model2index(self.vol.zmax, axis='z') == self.vol.nz - 1

    def test_dx_dy_dz(self):
        """Ensure that bounds work even when d's are negative."""
        for item in ['dx', 'dy', 'dz']:
            setattr(self.vol, item, -getattr(self.vol, item))
            self.test_model2index_mins()
            self.test_model2index_maxs()
            setattr(self.vol, item, -getattr(self.vol, item))

class TestCrop(TestBase):
    def test_crop_in_memory(self):
        self.vol.load()
        v = self.vol.crop(xmin=2200, xmax=2300, ymin=5000, ymax=5100,
                          zmin=3000, zmax=3100)
        vdat = self.vol[2200:2300, 5000:5100, 3000:3100]
        assert np.allclose(v.data, vdat)

class TestExtraction(TestBase):
    def test_arbitrary_section(self):
        x = [self.vol.xmin, self.vol.xmax]
        y = [self.vol.ymin, self.vol.ymax]
        dist = np.hypot(np.diff(x), np.diff(y))
        section, xi, yi = self.vol.extract_section(x, y)
        assert section.shape == (int(dist), self.vol.nz)

    def test_arbitrary_section_interpolation_shape(self):
        x, y = self.vol.index2model([0, 2, 5, 10], [0, 0, 0, 0])
        section, xi, yi = self.vol.extract_section(x, y)
        assert section.shape[0] == 10

    def test_arbitrary_section_interpolation_non_increasing(self):
        x, y = self.vol.index2model([0, 5, 2, 10], [0, 0, 0, 0])
        section, xi, yi = self.vol.extract_section(x, y)
        assert section.shape[0] == 16

    def test_arbitrary_section_interpolation_zlimits(self):
        x, y = self.vol.index2model([0, 2, 5, 10], [0, 0, 0, 0])
        zmin = self.vol.zmin + 50
        zmax = self.vol.zmax - 50
        section, xi, yi = self.vol.extract_section(x, y, zmin, zmax)
        assert section.shape == (10, 235)

    def test_xslice(self):
        slice1 = self.vol.XSlice(self.vol.xmin)
        slice2 = self.vol[self.vol.xmin, :, :].T
        slice3 = self.vol.data[0, :, :].T
        assert np.all(slice1 == slice2)
        assert np.all(slice2 == slice3)

