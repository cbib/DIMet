import os
from unittest import TestCase

import numpy as np

import pandas as pd

import matplotlib.pyplot as plt

import warnings

from omegaconf import DictConfig

from dimet.visualization import abundance_bars


class TestAbundanceBars(TestCase):
    def test_pile_up_abundance(self):
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
            'compartment': ['ex', 'ex', 'ex', 'ex']
        })
        result = abundance_bars.pile_up_abundance(df, metadata_df)
        self.assertTrue(result.shape == (12, 4))
        self.assertListEqual(list(result.columns),
                             ['timepoint', 'condition',
                              'metabolite', 'abundance'])
        self.assertTrue(any(
            np.array(result['abundance']) == np.array(
                [0.0098, 2.4536, 0.0098, 0.9814, 0.0133, 2.6672, 0.6668,
                 1.6003, 0.1268, 2.2194, 0.2536, 2.156]
            )))

    def test_plot_one_metabolite(self):
        warnings.filterwarnings("ignore")
        df = pd.DataFrame({
            'timepoint': ['t0', 't0'],
            'condition': ['beta', 'alpha'],
            'metabolite': ['m1', 'm1'],
            'abundance': [200, 700]
        })
        fig_this_metabolite, axs_k = plt.subplots(
            nrows=1, ncols=1,
            figsize=(5, 5))
        result = abundance_bars.plot_one_metabolite(
            df, "m1", axisx_var="timepoint",
            hue_var="condition", axisx_labeltilt=30,
            palette_choice="dark", curr_ax=axs_k, do_stripplot=False)
        bar = result.get_children()[0].get_facecolor()
        # one of the bars color, rgba
        self.assertAlmostEqual(bar[0], 0.062, 2)  # r
        self.assertAlmostEqual(bar[1], 0.144, 2)  # g
        self.assertAlmostEqual(bar[2],  0.435, 2)  # b

    def test_plot_abundance_bars_no_grid(self):
        warnings.filterwarnings("ignore")
        df = pd.DataFrame({
            'timepoint': ['t0', 't0'],
            'condition': ['beta', 'alpha'],
            'metabolite': ['m1', 'm1'],
            'abundance': [200, 700]
        })
        try:
            os.makedirs("../__pycache__/")
        except FileExistsError:
            pass
        cfg_m = DictConfig({'analysis': {
            'method': {'palette': 'dark',
                       'do_stripplot': False,
                       'x_text_modify_as': None,
                       'figure_format': 'svg'}
            }
        })
        result = abundance_bars.plot_abundance_bars_no_grid(
            df, ["m1"], "med", "total_abundance",
            axisx_var="timepoint",
            hue_var="condition",
            output_directory="../__pycache__/",
            axisx_labeltilt=30,
            width_each_subfig=1.2,
            height_each_subfig=2.4,
            cfg=cfg_m)
        self.assertTrue(result is None)
