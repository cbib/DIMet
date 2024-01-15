from unittest import TestCase

import numpy as np
import pandas as pd
from scipy import stats

from dimet.processing import bivariate_analysis


class TestBivariateAnalysis(TestCase):

    def test_conditions_to_comparisons(self):
        conditions_list = ["A", "B", "C", "D"]
        result = bivariate_analysis.conditions_to_comparisons(
            conditions_list
        )
        self.assertListEqual(result[0], ['A', 'B'])
        self.assertListEqual(result[1], ['A', 'C'])
        self.assertListEqual(result[2], ['A', 'D'])
        self.assertListEqual(result[3], ['B', 'C'])
        self.assertListEqual(result[4], ['B', 'D'])

    def test_compute_time_ordered_gmeans(self):
        metadata_df =  pd.DataFrame({
            'condition': ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'],
            'timenum': [0, 0, 1.5, 1.5, 4, 4, 7, 7],
            'name_to_plot': ['A-T0-1', 'A-T0-2', 'A-T1-1', 'A-T1-2',
                         'A-T4-1', 'A-T4-1', 'A-T7-1', 'A-T7-1']})
        metadata_df_sorted = metadata_df.sort_values(by="timenum")
        curr_condition = 'A'
        tmp = pd.DataFrame({ 'metabolite' : [1,2],
            'A-T0-1': [4.6, 4.6],  'A-T0-2': [6.5, 6.5], 'A-T1-1': [2, 2],
            'A-T1-2': [10.1, 8.4],  'A-T4-1':[5.6, 3.6], 'A-T4-1':[1.6, 1.8],
              'A-T7-1':[7.6, 7.6], 'A-T7-1':[8.1, 9.3]
        })
        row = tmp.iloc[1]
        sorted_time_numeric = [0, 1.5, 4, 7]
        result = bivariate_analysis.compute_time_ordered_gmeans(
            sorted_time_numeric, metadata_df_sorted,
            curr_condition, row
        )
        self.assertListEqual(result, [5.468089, 4.09878, 1.8, 9.3]
)




