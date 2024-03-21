#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski
"""

import logging
import os, sys

import hydra

from omegaconf import DictConfig, OmegaConf

from dimet.data import Dataset
from dimet.method import Method


logger = logging.getLogger(__name__)


@hydra.main(config_path="./config", config_name="config", version_base=None)
def main_run_analysis(cfg: DictConfig) -> None:
    logger.info(f"The current working directory is {os.getcwd()}")
    logger.info("Current configuration and defaults, %s",
                OmegaConf.to_yaml(cfg))

    dataset: Dataset = Dataset(
        config=hydra.utils.instantiate(cfg.analysis.dataset))
    dataset.preload()
    dataset.split_datafiles_by_compartment()
    method: Method = hydra.utils.instantiate(
        cfg.analysis.method).build()  # method factory

    method.run(cfg, dataset)


def fully_empty_args_custom_message() -> None:
    custom_message: str = (
        "\nDIMet is a tool for the Differential analysis of "
        "targeted Isotope-labeled Metabolomics data. \n\n  Please type "
        "'python -m dimet --help' or 'python -m dimet -h' for usage.\n"
    )
    print(custom_message)


if __name__ == "__main__":
    if len(sys.argv) <= 1:
        fully_empty_args_custom_message()
    else:
        main_run_analysis()
