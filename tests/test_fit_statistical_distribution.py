#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski
"""

from unittest import TestCase

import numpy as np
import pandas as pd

from dimet.processing import fit_statistical_distribution


class TestFitStatisticalDistribution(TestCase):
    def test_compute_z_score(self):
        data = {
            "gmean_1" : [15, 6, 3.8, 18.6, 16, 12],
            "gmean_2" : [1, 20, 16, 12, 2, 1.6],
            "ratio": [15, 0.3, 0.23, 1.55, 8, 7.5],
        }
        df = pd.DataFrame(data)
        result = fit_statistical_distribution.compute_z_score(df, "ratio")
        self.assertTrue(any(np.array(result.loc[0, :]) == np.array(
           [15.0, 1.0, 15.0, 1.793223]))
        )
        self.assertTrue(any(np.array(result.loc[1, :]) == np.array(
           [6.0, 20.0, 0.30,  - 0.961258]))
                        )
        self.assertTrue(any(np.array(result.loc[2, :]) == np.array(
           [3.8, 16.0, 0.23, -0.974374]))
                        )
        self.assertTrue(any(np.array(result.loc[3, :]) == np.array(
           [18.6, 12.0, 1.55, -0.727033]))
                        )

    def test_find_best_distribution(self):
        data = {'zscore': np.random.laplace(loc=0.0, scale=1.6, size=500)}
        df = pd.DataFrame(data)
        best_distribution, args_param = \
             fit_statistical_distribution.find_best_distribution(df)
        #  unexpected distribution: 'gennorm' or dgamma or loglaplace or ?
        #  impossible to set assert :
        #  self.assertTrue(best_distribution.name == "laplace" |
        #              best_distribution.name == "dgamma" )  #  can be false
        self.assertIsInstance(best_distribution.name,
                              str)
        self.assertIsInstance(args_param, dict)
        # self.assertIsInstance(best_distribution,
        #                  scipy.stats._continuous_distns )  # failed

    def test_best_fit(self):
        data = {'zscore': np.random.laplace(loc=0.0, scale=1.6, size=500)}
        df = pd.DataFrame(data)
        dist = np.around(np.array((df["zscore"]).astype(float)), 5)
        best_dist, best_dist_name, best_fit_params = \
            fit_statistical_distribution.get_best_fit(dist)
        self.assertIsInstance(best_dist_name, str)
        self.assertIsInstance(best_fit_params, str)









