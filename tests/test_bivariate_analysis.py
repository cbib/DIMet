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
                           'A-T4-2': [1.6, 1.8],
                           'A-T7-1': [7.6, 7.6], 'A-T7-2': [8.1, 9.3],
                           'B-T0-1': [8.1, 9.3], 'B-T0-2': [1.6, 1.8],
                           'B-T1-1': [6, 7],
                           'B-T1-2': [10.1, 8.4], 'B-T4-1': [3.2, 3.6],
                           'B-T4-2': [6.5, 6.5],
                           'B-T7-1': [4.6, 4.6], 'B-T7-2': [7.6, 7.6]})
        metadata_df = pd.DataFrame({
            'condition': ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',
                          'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B'],
            'timenum': [0, 0, 1.5, 1.5, 4, 4, 7, 7,
                        0, 0, 1.5, 1.5, 4, 4, 7, 7],
            'timepoint': ['0h', '0h', '1.5h', '1.5h', '4h', '4h', '7h', '7h',
                          '0h', '0h', '1.5h', '1.5h', '4h', '4h', '7h', '7h'],
            'name_to_plot': ['A-T0-1', 'A-T0-2', 'A-T1-1', 'A-T1-2',
                             'A-T4-1', 'A-T4-2', 'A-T7-1', 'A-T7-2',
                             'B-T0-1', 'B-T0-2', 'B-T1-1', 'B-T1-2',
                             'B-T4-1', 'B-T4-2', 'B-T7-1', 'B-T7-2']})
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
            [5.468089, 4.494441, 2.993326, 7.846018])

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

    def test_inner_gmean_dict_filler(self):
        k = 0
        inner_gmean_dict = {"metabolite": ['CoA', 'Ala'],
                            "gmean_arr_1": list()}

        df_a_group = pd.DataFrame({
            'A-T0-1': [0.2, 0.3, 0.5, 0.1, 0.2, 0.3, 0.4],
            'A-T0-2': [0.17, 0.33, 0.5, 0.1, 0.22, 0.27, 0.41],
            'A-T0-3': [0.2, 0.28, 0.52, 0.11, 0.18, 0.29, 0.42]})
        df_a_group.index = ['CoA_m+0', 'CoA_m+1', 'CoA_m+2',
                            'Ala_m+0', 'Ala_m+1', 'Ala_m+2', 'Ala_m+3']

        clue_isotopologue_df = pd.DataFrame({
            'isotopologue_name': ['CoA_m+0', 'CoA_m+1', 'CoA_m+2',
                                  'Ala_m+0', 'Ala_m+1', 'Ala_m+2', 'Ala_m+3'],
            'metabolite': ['CoA', 'CoA', 'CoA', 'Ala', 'Ala', 'Ala', 'Ala'],
            'm+x': [0, 1, 2, 0, 1, 2, 3]})
        clue_isotopologue_df['m+x'] = clue_isotopologue_df['m+x'].astype(int)

        result = bivariate_analysis.inner_gmean_dict_filler(
            k, inner_gmean_dict, df_a_group, clue_isotopologue_df
        )  # result:
        # {'metabolite': ['CoA', 'Ala'],
        # ' gmean_arr_1': [array([0.189454, 0.302643, 0.50658 ]),
        #                 array([0.103228, 0.199331, 0.286392, 0.409919])]}
        self.assertTrue(np.allclose(np.array(
            [0.189454, 0.302643, 0.50658]),
            result['gmean_arr_1'][0], rtol=1e-6))
        self.assertTrue(np.allclose(np.array(
            [0.103228, 0.199331, 0.286392, 0.409919]),
            result['gmean_arr_1'][1], rtol=1e-6))

    def test_compute_test_for_df_dict(self):
        df = pd.DataFrame({
            'metabolite': ['CoA', 'Ala'],
            'gmean_arr_1': [np.array([0.1894, 0.3026, 0.506]),
                            np.array([0.1032, 0.1993, 0.2863, 0.4099])],
            'gmean_arr_2': [np.array([0.506, 0.3026, 0.1894]),
                            np.array([0.4099, 0.2863, 0.1993, 0.1032])],
        })
        df_dict = {'T0': df}  # df comparing two conditions A vs B, at T0
        test = "pearson"
        result = bivariate_analysis.compute_test_for_df_dict(df_dict, test)
        self.assertEqual('T0', list(result.keys())[0])
        self.assertTrue(
            set(list(['metabolite', 'pvalue', 'correlation_coefficient']
                     )).issubset(set(list(result['T0'].columns)))
        )

    def test_compute_statistical_correlation(self):
        df = pd.DataFrame({
            'metabolite': ['CoA', 'Ala'],
            'gmean_arr_1': [np.array([0.1894, 0.3026, 0.506]),
                            np.array([0.1032, 0.1993, 0.2863, 0.4099])],
            'gmean_arr_2': [np.array([0.506, 0.3026, 0.1894]),
                            np.array([0.4099, 0.2863, 0.1993, 0.1032])],
        })
        df.index = ['CoA', 'Ala']
        test = "pearson"
        result = bivariate_analysis.compute_statistical_correlation(df, test)

        self.assertAlmostEqual(result.loc['CoA', 'correlation_coefficient'],
                               -0.947312, places=5)
        self.assertAlmostEqual(result.loc['CoA', 'pvalue'],
                               0.207574, places=5)
        self.assertAlmostEqual(result.loc['Ala', 'correlation_coefficient'],
                               -0.992587, places=5)
        self.assertAlmostEqual(result.loc['Ala', 'pvalue'],
                               0.007413, places=5)
