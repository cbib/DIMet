from unittest import TestCase

import numpy as np

import pandas as pd

import dimet.processing.pca_analysis as pca_analysis


class TestPcaAnalysis(TestCase):

    def test_handle_nan_values_before_pca(self):
        data = {
            "c1": [0.1, 0, np.nan, 2],
            "c2": [25, 20, np.nan, 35],
            "c3": [np.nan, 5, np.nan, 4],
            "c4": [10, 12, np.nan, 34],
        }
        df = pd.DataFrame(data)
        result = pca_analysis.handle_nan_values_before_pca(df)
        self.assertEqual(result.shape[0], 3)
        self.assertEqual(result.shape[1], 4)
        self.assertTrue(
            np.isnan(df.loc[0, 'c3']) and
            result.loc[0, 'c3'] == 0.1
        )

    def test_reduce_data_df(self):
        data = {
            "c1": [0.1, 0.1, 2.0],
            "c2": [25.0, 20.0, 35.0],
            "c3": [0.1, 5.0, 4.0],
            "c4": [10.0, 12.0, 34.0]
        }
        df = pd.DataFrame(data)
        result = pca_analysis.reduce_data_df(df)
        self.assertTrue(
            np.allclose(np.array(result.loc[0, :]),
                        np.array([0.0098, 2.4536, 0.0098, 0.9814]), 4)
        )
        self.assertTrue(
            np.allclose(np.array(result.loc[1, :]),
                        np.array([0.0133, 2.6672, 0.6668, 1.6003]), 4)
        )
        self.assertTrue(
            np.allclose(np.array(result.loc[2, :]),
                        np.array([0.1268, 2.2194, 0.2536, 2.1560]), 4)
        )

    def test_compute_pca(self):
        data = {
            'beta-1': [0.0098, 0.0133, 0.1268],
            'beta-2': [2.4536, 2.6672, 2.2194],
            'ctrl-1': [0.0098, 0.6668, 0.2536],
            'ctrl-2': [0.9814, 1.6003, 2.1560]
        }
        df = pd.DataFrame(data)
        metadata_df = pd.DataFrame({
            'name_to_plot': ['beta-1', 'beta-2', 'ctrl-1', 'ctrl-2'],
            'condition': ['beta-glu', 'beta-glu', 'control', 'control']
        })
        pc_df, var_explained_df = pca_analysis.compute_pca(df, metadata_df)

        self.assertTrue(
            np.allclose(np.array(pc_df['PC1']),
                        np.array([-1.81, 2.34, -1.35, 0.82]), 2)
        )
        self.assertTrue(
            np.allclose(np.array(pc_df['PC2']),
                        np.array([-0.10, -0.38, -0.13, 0.62]), 2)
        )
        self.assertTrue(
            np.allclose(np.array(pc_df['PC3']),
                        np.array([-0.2279, -0.0277, 0.2562, -0.0004]), 4)
        )
        self.assertListEqual(list(pc_df['name_to_plot']),
                             list(metadata_df['name_to_plot']))
        self.assertListEqual(list(pc_df['condition']),
                             list(metadata_df['condition']))
        # var_explained_df :
        self.assertTrue(
            np.allclose(
                np.array(var_explained_df['Explained Variance %']),
                np.array([94.26, 4.47, 0.98]), 2
            )
        )
        self.assertListEqual(list(var_explained_df['PC']),
                             ["PC1", "PC2", "PC3"])

    def test_pca_on_split_dataset(self):
        data = {
            'beta-1': [0.0098, 0.0133, 0.1268],
            'beta-2': [2.4536, 2.6672, 2.2194],
            'ctrl-1': [0.0098, 0.6668, 0.2536],
            'ctrl-2': [0.9814, 1.6003, 2.1560]
        }
        df = pd.DataFrame(data)
        metadata_df = pd.DataFrame({
            'name_to_plot': ['beta-1', 'beta-2', 'ctrl-1', 'ctrl-2'],
            'condition': ['beta-glu', 'beta-glu', 'control', 'control'],
            'timepoint': ['t0', 't0', 't0', 't0'],
            'short_comp': ['ex', 'ex', 'ex', 'ex']
        })
        description = ["my_file_name", "ex"]
        result_dict = pca_analysis.pca_on_split_dataset(
            df, metadata_df, chosen_column="condition",
            description=description
        )
        self.assertEqual(list(result_dict.keys())[0],
                         ('my_file_name', 'beta-glu', 'ex'))
        self.assertEqual(list(result_dict.keys())[1],
                         ('my_file_name', 'control', 'ex'))
        result_pc_df = result_dict[('my_file_name', 'beta-glu', 'ex')]['pc']
        self.assertAlmostEqual(np.array(result_pc_df['PC1'])[0],
                               -1.732, 3)
        self.assertAlmostEqual(np.array(result_pc_df['PC1'])[1],
                               1.732, 3)
        result_var_df = result_dict[('my_file_name', 'beta-glu', 'ex')]['var']
        self.assertAlmostEqual(np.array(result_var_df[
                                            'Explained Variance %'])[0],
                               100.0, 1)
        self.assertAlmostEqual(np.array(result_var_df[
                                            'Explained Variance %'])[1],
                               2.8e-31, 1)
        result_pc_df2 = result_dict[('my_file_name', 'control', 'ex')]['pc']
        self.assertAlmostEqual(np.array(result_pc_df2['PC2'])[1],
                               1.24e-16, 2)
        result_var_df2 = result_dict[('my_file_name', 'control', 'ex')]['var']
        self.assertAlmostEqual(np.array(result_var_df2[
                                            'Explained Variance %'])[1],
                               5.18e-31, 2)

    def test_pca_global_compartment_dataset(self):
        data = {
            'beta-1': [0.0098, 0.0133, 0.1268],
            'beta-2': [2.4536, 2.6672, 2.2194],
            'ctrl-1': [0.0098, 0.6668, 0.2536],
            'ctrl-2': [0.9814, 1.6003, 2.1560]
        }
        df = pd.DataFrame(data)
        metadata_df = pd.DataFrame({
            'name_to_plot': ['beta-1', 'beta-2', 'ctrl-1', 'ctrl-2'],
            'condition': ['beta-glu', 'beta-glu', 'control', 'control'],
            'timepoint': ['t0', 't0', 't0', 't0'],
            'short_comp': ['ex', 'ex', 'ex', 'ex']
        })
        description = ["my_file_name", "ex"]
        result_dict = pca_analysis.pca_global_compartment_dataset(
            df, metadata_df, description
        )
        self.assertEqual(list(result_dict.keys())[0],
                         ('my_file_name', 'ex'))
        result_pc_df = result_dict[('my_file_name', 'ex')]['pc']
        self.assertAlmostEqual(np.array(result_pc_df['PC1'])[0],
                               -1.814, 3)
        self.assertAlmostEqual(np.array(result_pc_df['PC1'])[1],
                               2.34, 2)
        self.assertAlmostEqual(np.array(result_pc_df['PC1'])[2],
                               -1.358, 3)
        self.assertAlmostEqual(np.array(result_pc_df['PC1'])[3],
                               0.8, 1)
        result_var_df = result_dict[('my_file_name', 'ex')]['var']
        self.assertAlmostEqual(np.array(result_var_df[
                                            'Explained Variance %'])[0],
                               94.266, 3)
        self.assertAlmostEqual(np.array(result_var_df[
                                            'Explained Variance %'])[1],
                               4.747, 3)
