import unittest
import numpy as np

import fastpli.model.sandbox as sb
import fastpli.objects as obj


class MainTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_triangular_grid(self):
        seeds = sb.seeds.triangular_grid(0, 0, 1, endpoint=False)
        self.assertTrue(seeds.size == 0)
        seeds = sb.seeds.triangular_grid(2, 0, 1, endpoint=False)
        self.assertTrue(seeds.size == 0)

        seeds = sb.seeds.triangular_grid(2, 0, 1, endpoint=True)
        self.assertTrue(seeds.shape[0] == 3)
        self.assertTrue(np.all(seeds[:, 1] == 0))
        self.assertTrue(np.all(seeds[:, 0] == np.array([0, 1, 2])))

        seeds = sb.seeds.triangular_grid(2, 2, 1, endpoint=True)
        self.assertTrue(seeds.shape[0] == 8)

    def test_crop_rectangle(self):
        seeds = sb.seeds.triangular_grid(2, 2, 1)
        new_seeds = sb.seeds.crop_rectangle(2, 2, seeds)
        self.assertTrue(np.array_equal(seeds, new_seeds))

        seeds = sb.seeds.triangular_grid(0, 0, 1)
        new_seeds = sb.seeds.crop_rectangle(2, 2, seeds, 1)
        self.assertTrue(new_seeds.size == 0)

        new_seeds = sb.seeds.crop_rectangle([-1, -1], [1, 1], seeds, 0)
        self.assertTrue(np.array_equal(new_seeds, [[0, 0]]))

        seeds = sb.seeds.triangular_grid(2, 2, 1)
        new_seeds = sb.seeds.crop_rectangle(2,
                                            2,
                                            seeds,
                                            radii=[1] * seeds.shape[0])
        self.assertTrue(new_seeds.size < seeds.size)

    def test_crop_circle(self):
        seeds = sb.seeds.triangular_grid(2, 2, 1)
        new_seeds = sb.seeds.crop_circle(100, seeds)
        self.assertTrue(np.array_equal(seeds, new_seeds))

        new_seeds = sb.seeds.crop_circle(1, seeds)
        self.assertTrue(
            np.array_equal(new_seeds,
                           [[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]]))

        new_seeds = sb.seeds.crop_circle(np.sqrt(2), seeds, [1, 1])
        self.assertTrue(np.array_equal(seeds, new_seeds))

        new_seeds = sb.seeds.crop_circle(1, seeds, center=[0, 0], radii=1)
        self.assertTrue(np.array_equal(new_seeds, [[0, 0]]))

        new_seeds = sb.seeds.crop_circle(1,
                                         seeds,
                                         center=[0, 0],
                                         radii=[1] * seeds.shape[0])
        self.assertTrue(np.array_equal(new_seeds, [[0, 0]]))

    def test_build_cylinder(self):
        seeds = sb.seeds.triangular_circle(10, 1)
        for m in ['p', 'r', 'c']:
            sb.build.cylinder(p=(0, 0, 0),
                              q=(10, 10, 10),
                              r_in=5,
                              r_out=8,
                              seeds=seeds,
                              radii=1,
                              alpha=np.deg2rad(20),
                              beta=np.deg2rad(160),
                              mode=m)
        self.assertTrue(True)

    def test_build_cuboid(self):
        p = np.array([0, 80, 50])
        q = np.array([40, 180, 100])
        d = np.max(np.abs(p - q)) * np.sqrt(3)
        seeds = sb.seeds.triangular_grid(a=d, b=d, spacing=5, center=True)
        sb.build.cuboid(p=p,
                        q=q,
                        phi=np.deg2rad(45),
                        theta=np.deg2rad(90),
                        seeds=seeds,
                        radii=1)
        self.assertTrue(True)

        seeds = sb.seeds.triangular_grid(300, 300, 2, center=True)
        fb = sb.build.cuboid(p=[-5] * 3,
                             q=[5] * 3,
                             phi=0,
                             theta=np.deg2rad(0),
                             seeds=seeds,
                             radii=1)
        cut_fb = obj.fiber_bundle.Cut(fb, [[-5] * 3, [5] * 3])
        for f0, f1 in zip(fb, cut_fb):
            self.assertTrue(np.array_equal(f0, f1))

    def test_build_bundle(self):
        traj = np.array([[0, 0, 0], [0, 0, 100]])
        seeds = np.array([[0, 0], [1, 1], [-1, 1]])
        fiber_bundle = sb.build.bundle(traj, seeds, 1)

        for i in range(len(fiber_bundle)):
            self.assertTrue(
                np.array_equal(fiber_bundle[i][0, :3],
                               [seeds[i][0], seeds[i][1], 0]))
            self.assertTrue(
                np.array_equal(fiber_bundle[i][1, :3],
                               [seeds[i][0], seeds[i][1], 100]))


if __name__ == '__main__':
    unittest.main()
