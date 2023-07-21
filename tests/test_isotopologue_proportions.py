from unittest import TestCase

import numpy as np
import pandas as pd

from dimet.visualization import isotopologue_proportions


class TestIsotopologueProportionsPlot(TestCase):

    def test_isotopologue_proportions_2piled_df(self):
        data = {
            'b-1': [0.2, 0.5, 0.3],
            'b-2': [0.1, 0.7, 0.2],
            'ctl-1': [0.6, 0.1, 0.3],
            'ctl-2': [0.3, 0.2, 0.5]
        }
        df = pd.DataFrame(data)
        df.index = ['cit_m+0', 'cit_m+1', 'cit_m+2']
        metadata_df = pd.DataFrame({
            'name_to_plot': ['b-1', 'b-2', 'ctl-1', 'ctl-2'],
            'condition': ['b', 'b', 'ctl', 'ctl'],
            'timepoint': ['t0', 't0', 't0', 't0'],
            'timenum': [0, 0, 0, 0],
            'short_comp': ['ex', 'ex', 'ex', 'ex'],
            'original_name': ['', '', '', '']
        })
        result = isotopologue_proportions.isotopologue_proportions_2piled_df(
            df, metadata_df
        )
        self.assertListEqual(
            list(result.columns),
            ["timenum", "condition", "isotopologue_name",
             "Isotopologue Contribution (%)"]
        )
        self.assertTrue(any(np.array(result.loc[0, :]) == np.array(
            ['0', 'b', 'cit_m+0', 20.0]
        )))
        self.assertTrue(any(np.array(result.loc[3, :]) == np.array(
            ['0', 'ctl', 'cit_m+0', 30.0]
        )))

    def test_massage_isotopologues(self):
        df = pd.DataFrame({
            'timenum': ['0', '0', '0', '0', '0', '0',
                        '0', '0', '0', '0', '0', '0'],
            'condition': ['b', 'b', 'ctl', 'ctl', 'b', 'b',
                          'ctl', 'ctl', 'b', 'b', 'ctl', 'ctl'],
            'isotopologue_name': ['cit_m+0', 'cit_m+0', 'cit_m+0',
                                  'cit_m+0', 'cit_m+1', 'cit_m+1',
                                  'cit_m+1', 'cit_m+1', 'cit_m+2',
                                  'cit_m+2', 'cit_m+2', 'cit_m+2'],
            'Isotopologue Contribution (%)': [20, 10, 60, 30, 50, 70,
                                              10, 20, 30, 20, 30, 50]
        })
        result = isotopologue_proportions.massage_isotopologues(
            df
        )
        self.assertEqual(result.columns[-1], "m+x")
        self.assertEqual(result.columns[-2], "metabolite")
        self.assertEqual(result['metabolite'].unique().item(), "cit")

    def test_prepare_means_replicates(self):
        df = pd.DataFrame({
            'timenum': ['0', '0', '0', '0', '0', '0',
                        '0', '0', '0', '0', '0', '0'],
            'condition': ['b', 'b', 'ctl', 'ctl', 'b', 'b',
                          'ctl', 'ctl', 'b', 'b', 'ctl', 'ctl'],
            'isotopologue_name': ['cit_m+0', 'cit_m+0', 'cit_m+0',
                                  'cit_m+0', 'cit_m+1', 'cit_m+1',
                                  'cit_m+1', 'cit_m+1', 'cit_m+2',
                                  'cit_m+2', 'cit_m+2', 'cit_m+2'],
            'Isotopologue Contribution (%)': [20, 10, 60, 30, 50, 70,
                                              10, 20, 30, 20, 30, 50],
            'metabolite': ['cit' for i in range(12)],
            'm+x': ['m+0', 'm+0', 'm+0', 'm+0', 'm+1', 'm+1', 'm+1', 'm+1',
                    'm+2', 'm+2', 'm+2', 'm+2']
        })
        result_dict = isotopologue_proportions.prepare_means_replicates(
            df, metaboli_selected=['cit']
        )
        self.assertListEqual(list(result_dict.keys()), ['cit'])
        self.assertEqual(result_dict['cit'].shape, (6, 5))
        self.assertTrue(any(np.array(
            result_dict['cit']['Isotopologue Contribution (%)']) == np.array(
            [15.0, 60.0, 25.0, 45.0, 15.0, 40.0]
        )))

    def test_add_combined_conditime(self):
        df = pd.DataFrame({
            'timenum': ['0', '0', '0', '0', '0', '0'],
            'condition': ['b', 'ctl', 'b',
                          'ctl',  'b',  'ctl'],
            'isotopologue_name': ['cit_m+0', 'cit_m+0',
                                  'cit_m+1', 'cit_m+1',
                                   'cit_m+2', 'cit_m+2'],
            'Isotopologue Contribution (%)':   [15.0, 60.0, 25.0,
                                                45.0, 15.0, 40.0],
            'metabolite': ['cit' for i in range(6)],
            'm+x': ['m+0', 'm+0',  'm+1', 'm+1',
                    'm+2', 'm+2']})
        dfs_dict = {'cit': df}
        combined_tc_levels = ['0 : ctl', '0 : b', '0xemptyspace']
        result = isotopologue_proportions.add_combined_conditime(
            dfs_dict, combined_tc_levels)
        self.assertListEqual(
            result['cit']['time_and_condition'].tolist(),
            ['0 : b', '0 : ctl', '0 : b', '0 : ctl', '0 : b', '0 : ctl']
        )
        self.assertEqual(
            result['cit']['time_and_condition'].cat.categories[0],
            '0 : ctl'
        )
        self.assertEqual(
            result['cit']['time_and_condition'].cat.categories[1],
            '0 : b'
        )
        self.assertEqual(
            result['cit']['time_and_condition'].cat.categories[2],
            '0xemptyspace'
        )

    def test_add_categorical_time(self):
        df = pd.DataFrame({
            'timenum': ['2','3','1','3','2','1'],
            'condition': ['b', 'ctl', 'b',
                          'ctl',  'b',  'ctl'],
            'isotopologue_name': ['cit_m+0', 'cit_m+0',
                                  'cit_m+1', 'cit_m+1',
                                   'cit_m+2', 'cit_m+2'],
            'Isotopologue Contribution (%)':   [15.0, 60.0, 25.0,
                                                45.0, 15.0, 40.0],
            'metabolite': ['cit' for i in range(6)],
            'm+x': ['m+0', 'm+0',  'm+1', 'm+1',
                    'm+2', 'm+2']})
        dfs_dict = {'cit': df}
        levels_time_str = ['1', '2', '3']
        result = isotopologue_proportions.add_categorical_time(
            dfs_dict, levels_time_str
        )
        self.assertEqual(
            result['cit']['timenum'].cat.categories[0], '1')
        self.assertEqual(result['cit']['timenum'].cat.categories[1], '2')
        self.assertEqual(
            result['cit']['timenum'].cat.categories[2], '3'
        )

    def test_add_xemptyspace_tolabs(self):
        time_levels_list = [str(i) for i in sorted([1, 3, 2, 4])]
        conditions = ['b', 'ctl']
        result = isotopologue_proportions.add_xemptyspace_tolabs(
            conditions, time_levels_list
        )
        self.assertListEqual(result[0], ['b', 'ctl', 'xemptyspace'])
        self.assertListEqual(
            result[1],
            ['1 : b', '1 : ctl', '1xemptyspace',
             '2 : b', '2 : ctl', '2xemptyspace',
             '3 : b', '3 : ctl', '3xemptyspace',
             '4 : b', '4 : ctl', '4xemptyspace'])



