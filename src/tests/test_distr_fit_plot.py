#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski
"""

from unittest import TestCase

import helpers

import numpy as np

import pandas as pd

import scipy.stats as stats

import processing.fit_statistical_distribution as fit_statistical_distribution
import visualization.distr_fit_plot as distr_fit_plot


class TestDistrFitPlot(TestCase):
    def test_make_pdf(self):
        dist = getattr(stats, "gennorm")
        params_dict =  {'beta': 1.09, 'loc': 0.01, 'scale': 1.77}
        params = list(params_dict.values())
        result = distr_fit_plot.make_pdf(dist, params, size=10000)
        self.assertAlmostEqual(result.index[0],-5.9147, 2 )
        self.assertAlmostEqual(result.index[-3], 5.932, 2 )
        self.assertAlmostEqual(result.to_list()[1],  0.007, 3 )
        self.assertAlmostEqual(result.to_list()[-1], 0.00699, 3)

