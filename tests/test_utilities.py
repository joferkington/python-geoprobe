from geoprobe import utilities
import geoprobe
import numpy as np

class TestBboxMethods:
    def setup(self):
        self.bbox1 = [  0, 1,   1, 2]
        self.bbox2 = [0.5, 2, 1.5, 3]
        self.bbox3 = [10, 20, 30, 40]

    def test_intersects(self):
        assert utilities.bbox_intersects(self.bbox1, self.bbox2)
        assert not utilities.bbox_intersects(self.bbox1, self.bbox3)
    
    def test_intersection(self):
        overlap = utilities.bbox_intersection(self.bbox1, self.bbox2)
        assert overlap == [0.5, 1, 1.5, 2]
        assert utilities.bbox_intersection(self.bbox1, self.bbox3) is None

    def test_union(self):
        assert utilities.bbox_union(self.bbox1, self.bbox2) == [0, 2, 1, 3]
        assert utilities.bbox_union(self.bbox1, self.bbox3) == [0, 20, 1, 40]

class TestPlaneFitting:
    def generate_random_basis(self, state=None):
        if state is None:
            state = np.random.RandomState()
        a, b = state.rand(2,3)
        c = np.cross(a, b)
        b = np.cross(b, c)
        a, b, c = [item / np.sqrt(np.sum(item**2)) for item in (a,b,c)]
        return np.vstack([a, b, c]).T

    def setup(self):
        x, y = np.mgrid[40:50, 30:50]
        z = np.ones_like(x) + 8
        x, y, z = (item.astype(np.float).ravel() for item in [x,y,z])
        self.coords = np.vstack([x, y, z]).T
        self.x, self.y, self.z = x, y, z

    def test_xy_plane(self):
        norm = geoprobe.utilities.fit_plane(self.x, self.y, self.z)
        assert np.allclose(norm, [0, 0, 1, -9])

    def test_xz_plane(self):
        xz_basis = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]]).T
        x, y, z = self.coords.dot(np.linalg.inv(xz_basis)).T
        norm = geoprobe.utilities.fit_plane(x, y, z)
        assert np.allclose(norm, [0, 1, 0, -9])

    def test_yz_plane(self):
        yz_basis = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]).T
        x, y, z = self.coords.dot(np.linalg.inv(yz_basis)).T
        norm = geoprobe.utilities.fit_plane(x, y, z)
        assert np.allclose(norm, [1, 0, 0, -9])

    def test_assortment(self):
        state = np.random.RandomState(1977)
        for _ in range(10):
            basis = self.generate_random_basis(state)
            x, y, z = self.coords.dot(np.linalg.inv(basis)).T
            norm = np.array(geoprobe.utilities.fit_plane(x, y, z))
            a, b, c = basis[:,-1]
            expected = [a, b, c, -9]
            assert np.allclose(norm, expected) or np.allclose(-norm, expected)

def test_isopach():
    x, y = np.mgrid[50:100, 50:100]
    z1 = np.hypot(x - x.mean(), y - y.mean())
    z2 = 2 * z1

    hor1 = geoprobe.horizon(x=x.ravel(), y=y.ravel(), z=z1.ravel())
    hor2 = geoprobe.horizon(x=x.ravel(), y=y.ravel(), z=z2.ravel())

    iso = geoprobe.utilities.create_isopach(hor2, hor1)
    assert np.allclose(iso.x, x.ravel())
    assert np.allclose(iso.y, y.ravel())
    assert np.allclose(z2 - z1, iso.grid)
