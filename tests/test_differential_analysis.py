from unittest import TestCase

import numpy as np
import pandas as pd
from scipy import stats

from dimet.helpers import countnan_samples
from dimet.processing import differential_analysis


class TestDifferentialAnalysis(TestCase):

    def test_compute_span_incomparison(self):
        data = {
            "c1": [15, 22, 310],
            "c2": [8, 30, 220],
            "c3": [11, 25, 170],
            "c4": [9, 33, 100],
        }
        df = pd.DataFrame(data)
        result = differential_analysis.compute_span_incomparison(
            df, [["c1", "c2"], ["c3", "c4"]])

        self.assertTrue(result.at[0, 'span_allsamples'] == float(7))
        self.assertTrue(result.at[1, 'span_allsamples'] == float(11))
        self.assertTrue(result.at[2, 'span_allsamples'] == float(210))

    def test_distance_or_overlap(self):
        data = {
            "c1": [15, 22, 310],
            "c2": [8, 30, 220],
            "c3": [11, 25, 170],
            "c4": [9, 33, 100],
        }
        df = pd.DataFrame(data)
        result = differential_analysis.distance_or_overlap(
            df, [["c1", "c2"], ["c3", "c4"]]
        )
        self.assertTrue(result.at[0, "distance"] == float(-2))
        self.assertTrue(result.at[1, "distance"] == float(-5))
        self.assertTrue(result.at[2, "distance"] == float(50))

    def test_select_rows_with_sufficient_non_nan_values(self):
        data = {
            "c1": [15, 22, 310, np.nan], "c2": [8, np.nan, 2, 1],
            "c3": [11, np.nan, 70, np.nan], "c4": [9, 33, 100, np.nan],
            "c5": [3, np.nan, 5, 2], "c6": [np.nan, 8, 4, np.nan]
        }
        df = pd.DataFrame(data)
        groups = [["c1", "c2", "c3"], ["c4", "c5", "c6"]]
        df = countnan_samples(df, groups)
        result_good, result_bad = differential_analysis. \
            select_rows_with_sufficient_non_nan_values(df, groups)
        self.assertEqual(result_good.shape, (2, 8))
        self.assertEqual(result_bad.shape, (2, 8))
        self.assertTrue(np.any(np.array(result_good.loc[0, :]) ==
                               np.array(
                                   [15.0, 8.0, 11.0, 9.0, 3.0, np.nan, 0, 1]))
                        )

    def test_compute_mann_whitney_all_h0(self):
        array1 = np.array([722, 760, 750, 700])
        array2 = np.array([150, 177, 165, 110])
        result = differential_analysis.compute_mann_whitney_allH0(array1,
                                                                  array2)
        self.assertEqual(result[0], 16.0)
        self.assertAlmostEqual(result[1], 0.014285, 4)

    def test_run_statistical_test(self):
        data = {
            "c1": [15, 310], "c2": [8, 2],
            "c3": [11, 70], "c4": [9, 100],
            "c5": [3, 5], "c6": [np.nan, 4]
        }
        df = pd.DataFrame(data)
        df.index = ['met1', 'met3']
        groups = [["c1", "c2", "c3"], ["c4", "c5", "c6"]]
        result = differential_analysis.run_statistical_test(df,
                                                            groups,
                                                            'BrMu')
        self.assertAlmostEqual(
            result.loc[result['metabolite'] == "met1", "pvalue"].item(),
            0.131399, 6)
        self.assertAlmostEqual(
            result.loc[result['metabolite'] == "met3", "pvalue"].item(),
            0.436360, 6)

    def test_auto_detect_tailway(self):
        data = {'zscore': np.random.laplace(loc=0.0, scale=1.6, size=500)}
        df = pd.DataFrame(data)
        best_distribution = getattr(stats, 'gennorm')
        args_param = {'beta': 1.07, 'loc': 0.06, 'scale': 1.72}
        autoset_tailway = differential_analysis.auto_detect_tailway(
                 df, best_distribution, args_param
            )
        self.assertIsInstance(autoset_tailway, str)
        self.assertTrue(autoset_tailway in ["right-tailed", "two-sided"])


    def test_reorder_columns_diff_end(self):
        data = {
            "distance": [2, 1.5], "span_allsamples": [4, 6],
            "other_exta_column": [400, 500],
            "distance/span": [0.5, 0.25], "count_nan_samples_group1": [0, 0],
            "count_nan_samples_group2": [0, 0], "pvalue": [1e-3, 1e-4],
            "padj": [1e-3, 1e-4], "log2FC": [4., 5235],
            "FC": [8, 23], "compartment": ["med", "med"],
        }
        df = pd.DataFrame(data)
        df.index = ['met1', 'met3']
        result = differential_analysis.reorder_columns_diff_end(df)
        self.assertTrue(
            any(np.array(result.loc["met1", :]) == np.array(
                [4.0, 1e-3, 1e-3, 0.5, 8, 0, 0, 2.0, 4, 'med', 400]
            )))
        self.assertEqual(result.shape, (2, 11))

    def test_round_result_float_columns(self):
        data = {
            "distance": [2, 1.5], "span_allsamples": [4, 6],
            "distance/span": [0.59595959595959, 0.25],
            "count_nan_samples_group1": [0, 0],
            "count_nan_samples_group2": [0, 0],
            "pvalue": [0, 0],
            "padj": [1.54354343434e-3, 1.3543434354335e-4],
            "log2FC": [4.063030512104, 5.0202235],
            "FC": [8, 23], "compartment": ["med", "med"],
        }
        df = pd.DataFrame(data)
        df.index = ['met1', 'met3']
        result = differential_analysis.round_result_float_columns(df)
        self.assertTrue(
            any(np.array(result["padj"]) == np.array([0.001544, 0.000135])
                )
        )
        self.assertTrue(
            any(np.array(result["log2FC"]) == np.array([4.063031, 5.020224])
                )
        )

    def test_time_course_auto_list_comparisons(self):
        metadata = pd.DataFrame({
            'condition': ['cond1', 'cond1', 'cond1', 'cond1',
                          'cond2', 'cond2', 'cond2', 'cond2'],
            'timenum': [1, 2.7, 3, 1, 2.7, 3, 4, 4],
            'timepoint': ['1h', '2.7h', '3h', '1h', '2.7h', '3h', '4h', '4h']
        })
        result = differential_analysis.time_course_auto_list_comparisons(
            metadata
        )
        self.assertListEqual(result[0], [['cond2', '4h'], ['cond2', '3h']])
        self.assertListEqual(result[1], [['cond1', '3h'], ['cond1', '2.7h']])
        self.assertListEqual(result[2], [['cond2', '3h'], ['cond2', '2.7h']])
        self.assertListEqual(result[3], [['cond1', '2.7h'], ['cond1', '1h']])
