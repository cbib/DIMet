#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski
"""

import logging
import os, sys

import hydra

from dimet.method import Method
from omegaconf import DictConfig, OmegaConf

from dimet.data import Dataset

import click

logger = logging.getLogger(__name__)


@hydra.main(config_path="./config", config_name="config", version_base=None)
def main_run_analysis(cfg: DictConfig) -> None:
    logger.info(f"The current working directory is {os.getcwd()}")
    logger.info("Current configuration is %s", OmegaConf.to_yaml(cfg))

    dataset: Dataset = Dataset(
        config=hydra.utils.instantiate(cfg.analysis.dataset))
    dataset.preload()
    dataset.split_datafiles_by_compartment()
    method: Method = hydra.utils.instantiate(
        cfg.analysis.method).build()  # method factory

    method.run(cfg, dataset)


@click.command(context_settings=dict(
        help_option_names=['-h', '--help'])  # needed to recognize -h as well
    )
@click.option('-cd', type=str, required=True,
              help="Directory containing the general configuration "
                   "(see -cn)")
@click.option('-cn', type=str, required=True,
              help="File name of the general configuration "
                   "(only one .yaml file name)")
def main_run_analysis__or_display_help(cd, cn) -> None:
    main_run_analysis()  # run: hydra handles further arguments validation


if __name__ == "__main__":
    main_run_analysis__or_display_help()
