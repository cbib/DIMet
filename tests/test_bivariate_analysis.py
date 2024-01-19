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

    def test_metabolite_time_profiles_gmean_df_dict(self):
        df = pd.DataFrame({'metabolite': [1, 2],
                        'A-T0-1': [4.6, 4.6], 'A-T0-2': [6.5, 6.5],
                        'A-T1-1': [2, 2],
                        'A-T1-2': [10.1, 8.4], 'A-T4-1': [5.6, 3.6],
                        'A-T4-1': [1.6, 1.8],
                        'A-T7-1': [7.6, 7.6], 'A-T7-1': [8.1, 9.3],
                        'B-T0-1': [8.1, 9.3], 'B-T0-2': [1.6, 1.8],
                        'B-T1-1': [6, 7],
                        'B-T1-2': [10.1, 8.4], 'B-T4-1': [3.2, 3.6],
                        'B-T4-1': [6.5, 6.5],
                        'B-T7-1': [4.6, 4.6], 'B-T7-1': [7.6, 7.6]})
        metadata_df = pd.DataFrame({
            'condition': ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
                          'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B'],
            'timenum': [0, 0, 1.5, 1.5, 4, 4, 7, 7,
                        0, 0, 1.5, 1.5, 4, 4, 7, 7],
            'timepoint' : ['0h','0h', '1.5h', '1.5h','4h', '4h', '7h', '7h',
                           '0h','0h', '1.5h', '1.5h','4h', '4h', '7h', '7h'],
            'name_to_plot': ['A-T0-1', 'A-T0-2', 'A-T1-1', 'A-T1-2',
                             'A-T4-1', 'A-T4-1', 'A-T7-1', 'A-T7-1',
                             'B-T0-1', 'B-T0-2', 'B-T1-1', 'B-T1-2',
                             'B-T4-1', 'B-T4-1', 'B-T7-1', 'B-T7-1']})
        comparison = ["A", "B"]
        result = bivariate_analysis.metabolite_time_profiles_gmean_df_dict(
            df, metadata_df, comparison
        )  # result :
        # {'metabo_time_profile':
        #  metabolite       gmean_arr_1                 gmean_arr_2
        #  0    [5.468089, 4.494441, 1.6, 8.1]   [3.6, 7.7846, 6.5, 7.6]
        #  1    [5.468089, 4.09878, 1.8, 9.3]  [4.091455, 7.668116, 6.5, 7.6]}
        self.assertTrue(list(result.keys())[0] == 'metabo_time_profile')
        self.assertListEqual(
            result['metabo_time_profile'].iloc[0, 1].tolist(),
            [5.468089, 4.494441, 1.6, 8.1])

    def test_modify_gmean_by_sanity(self):
        df = pd.DataFrame({'A-T0-1': [4.6, np.nan],
                           'A-T0-2': [6.5, np.nan],
                           'A-T0-3': [np.nan, 2]})
        df.index = ["metabolite1", "metabolite2"]
        df = df.assign(gmean=df.apply(lambda x: stats.gmean(x.dropna()),
                                      axis=1))
        # df.gmean.to_numpy() :  [5.46808925,  2. ]
        result = bivariate_analysis.modify_gmean_by_sanity(df)
        # result : [5.46808925, nan] np array
        self.assertFalse(np.all(df["gmean"].to_numpy() == result))
        self.assertAlmostEqual(result[0], 5.46808925, places=5)
        self.assertTrue(np.isnan(result[1]))
