from unittest import TestCase

import numpy as np

from dimet.visualization import pca_plot


class TestPcaPlot(TestCase):

    def test_eigsorted(self):
        xdata = np.arange(-10, 10, 1)
        ydata = np.arange(-1, 1, 0.1)
        # get values to build the ellipse
        cov = np.cov(xdata, ydata)
        vals, vecs = pca_plot.eigsorted(cov)
        self.assertAlmostEqual(vals[0], 3.535e+01, 2)
        self.assertAlmostEqual(vals[1], 1.665e-16, 2)
        self.assertAlmostEqual(vecs[0][0], -0.995, 2)
        self.assertAlmostEqual(vecs[0][1],  0.099, 2)
        self.assertAlmostEqual(vecs[1][0], -0.099, 2)
        self.assertAlmostEqual(vecs[1][1], -0.995, 2)




