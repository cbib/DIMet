#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski
"""

from unittest import TestCase

import numpy as np
import pandas as pd
from scipy import stats

from dimet.processing import fit_statistical_distribution

np.random.seed(123)


class TestFitStatisticalDistribution(TestCase):
    def test_compute_z_score(self):
        data = {
            "gmean_1": [15, 6, 3.8, 18.6, 16, 12],
            "gmean_2": [1, 20, 16, 12, 2, 1.6],
            "ratio": [15, 0.3, 0.23, 1.55, 8, 7.5],
        }
        df = pd.DataFrame(data)
        result = fit_statistical_distribution.compute_z_score(df, "ratio")
        self.assertTrue(any(np.array(result.loc[0, :]) == np.array(
            [15.0, 1.0, 15.0, 1.793223]))
                        )
        self.assertTrue(any(np.array(result.loc[1, :]) == np.array(
            [6.0, 20.0, 0.30, - 0.961258]))
                        )
        self.assertTrue(any(np.array(result.loc[2, :]) == np.array(
            [3.8, 16.0, 0.23, -0.974374]))
                        )
        self.assertTrue(any(np.array(result.loc[3, :]) == np.array(
            [18.6, 12.0, 1.55, -0.727033]))
                        )

    def test_best_fit_distribution(self):
        np.random.seed(123)
        data = {'zscore': np.random.laplace(loc=0.0, scale=1.6, size=500)}
        df = pd.DataFrame(data)
        MYDISTRIBUTIONS = [stats.laplace, stats.johnsonsu, stats.pareto]
        name, params = fit_statistical_distribution.best_fit_distribution(
            df, MYDISTRIBUTIONS, bins=200
        )
        self.assertIsInstance(name, str)
        self.assertTrue(name in ['laplace', 'johnsonsu'])
        self.assertAlmostEqual(params[0], -0.0257, places=3)
        self.assertAlmostEqual(params[1], 1.03437, places=3)
        self.assertAlmostEqual(params[2], -0.0834, places=3)
        self.assertAlmostEqual(params[3], 1.5379, places=3)

    def test_get_best_fit(self):
        data = {'zscore': np.random.laplace(loc=0.0, scale=1.6, size=500)}
        df = pd.DataFrame(data)
        dist = np.around(np.array((df["zscore"]).astype(float)), 5)
        MYDISTRIBUTIONS = [stats.laplace, stats.norm,
                           stats.johnsonsu, stats.pareto]
        best_dist, best_dist_name, best_fit_params = \
            fit_statistical_distribution.get_best_fit(dist,
                                                      MYDISTRIBUTIONS)
        self.assertTrue(best_dist_name in ['laplace', 'johnsonsu'])
        self.assertIsInstance(best_fit_params, str)

    def test_find_best_distribution(self):
        data = {'zscore': np.random.laplace(loc=0.0, scale=1.6, size=500)}
        df = pd.DataFrame(data)
        MYDISTRIBUTIONS = [stats.laplace, stats.norm,
                           stats.johnsonsu, stats.pareto]
        result, params = fit_statistical_distribution.find_best_distribution(
            df, MYDISTRIBUTIONS
        )
        self.assertTrue(result.name in ['laplace', 'johnsonsu'])
        self.assertIsInstance(params, dict)
        self.assertListEqual(list(params.keys()), ['loc', 'scale'])
        self.assertTrue(params['loc'] <= 0)
        self.assertTrue(params['scale'] >= 1.6)
        self.assertAlmostEqual(params['loc'], -0.03, places=1)
        self.assertAlmostEqual(params['scale'], 1.85, places=1)

    def test_compute_p_value(self):
        data = {'zscore': np.random.laplace(loc=0.0, scale=1.6, size=500)}
        df = pd.DataFrame(data)

        best_dist = stats.laplace
        args_param = {'loc': -0.03, 'scale': 1.85}
        result = fit_statistical_distribution.compute_p_value(
            df, "right-tailed", best_dist, args_param
        )
        self.assertEqual(len(result['pvalue']), len(df["zscore"]) )
        self.assertAlmostEqual(result['pvalue'][0], 0.7572, places=1)
        self.assertAlmostEqual(result['pvalue'][1], 0.5879, places=1)
        self.assertAlmostEqual(result['pvalue'][2], 0.1721, places=1)
        self.assertAlmostEqual(result['pvalue'][4], 0.1721, places=1)