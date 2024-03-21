import logging
import os
from pathlib import Path
from typing import Dict, Literal, Optional, Set

import pandas as pd
from hydra.core.hydra_config import HydraConfig
from omegaconf import DictConfig, ListConfig
from pydantic import BaseModel as PydanticBaseModel

from dimet.constants import molecular_types_for_metabologram
from dimet.helpers import (df_to_dict_by_compartment,
                           drop_all_nan_metabolites_on_comp_frames,
                           extfind, set_samples_names,
                           verify_metadata_sample_not_duplicated,
                           message_bad_separator_input)


class BaseModel(PydanticBaseModel):
    class Config:
        arbitrary_types_allowed = True


logger = logging.getLogger(__name__)


class DatasetConfig(BaseModel):
    label: str
    name: str
    subfolder: str
    metadata: str
    # We allow for some default values for the following files
    # will be ignored if they do not exist

    # First condition is the reference condition (control)
    # Conditions should be a subset of the medata corresponding column
    conditions: ListConfig
    abundances: str = "AbundanceCorrected"
    mean_enrichment: str = "MeanEnrichment13C"
    isotopologue_proportions: str = "IsotopologuesProportions"
    isotopologues: str = "Isotopologues"  # isotopologue absolute values

    def build(self) -> "Dataset":
        return Dataset(config=self)


class Dataset(BaseModel):
    config: DatasetConfig
    sub_folder_absolute: str = None
    metadata_df: Optional[pd.DataFrame] = None
    abundances_df: Optional[pd.DataFrame] = None
    mean_enrichment_df: Optional[pd.DataFrame] = None
    isotopologue_proportions_df: Optional[pd.DataFrame] = None
    isotopologues_df: Optional[pd.DataFrame] = None
    available_datasets: Set[
        Literal[
            "metadata", "abundances", "mean_enrichment",
            "isotopologue_proportions", "isotopologues"]
    ] = set()
    compartmentalized_dfs: Dict[str, Dict[str, pd.DataFrame]] = {}

    def preload(self):
        # check if we have a relative or absolute path, compute the absolute
        # path then load the data using pandas
        # if the path is relative, we assume it is relative to the original
        # CWD (befre hydra changed it)
        # store the data in self.metadata
        original_cwd = HydraConfig.get().runtime.cwd
        logger.info("Current config directory is %s", original_cwd)
        if not self.config.subfolder.startswith("/"):
            self.sub_folder_absolute = os.path.join(Path(original_cwd),
                                                    "data",
                                                    self.config.subfolder)
            logger.info("looking for data in %s", self.sub_folder_absolute)
        else:
            self.sub_folder_absolute = self.config.subfolder

        # start loading the dataframes
        ext = self.get_files_extension_as_dict()  # extension str by file
        file_paths = [
            ("metadata", os.path.join(
                self.sub_folder_absolute,
                self.config.metadata + "." + ext['metadata'])),
            ("abundances", os.path.join(
                self.sub_folder_absolute,
                self.config.abundances + "." + ext['abundances'])),
            ("mean_enrichment", os.path.join(
                self.sub_folder_absolute,
                self.config.mean_enrichment + "." + ext['mean_enrichment'])),
            ("isotopologue_proportions", os.path.join(
                self.sub_folder_absolute,
                self.config.isotopologue_proportions + "." + ext['isotopologue_proportions'])),
            ("isotopologues", os.path.join(
                self.sub_folder_absolute,
                self.config.isotopologues + "." + ext['isotopologues'])),
        ]
        dfs = []
        for label, file_path in file_paths:
            try:
                if label != "metadata":
                    # the quantifications dfs take first column as index
                    # (metabolites), regardless the name of that column
                    dfs.append(pd.read_csv(file_path, sep="\t",
                                           header=0, index_col=0))
                else:
                    dfs.append(pd.read_csv(file_path, sep="\t", header=0))
                self.available_datasets.add(label)
            except FileNotFoundError:
                if file_path.endswith(self.config.isotopologues + "." + ext['isotopologues']):
                    message_detail = "isotopologue absolute values missing"
                    logger.critical(
                        "File %s not found (%s), continue"
                        % (file_path, message_detail))
                else:
                    logger.critical("File %s not found, continue",
                                    file_path)
                dfs.append(None)
            except Exception as e:
                logger.error(
                    "Failed to load file %s during preload, aborting",
                    file_path)
                raise e

        (
            self.metadata_df,
            self.abundances_df,
            self.mean_enrichment_df,
            self.isotopologue_proportions_df,
            self.isotopologues_df,
        ) = dfs

        # log the first 5 rows of the metadata
        logger.info("Loaded metadata: \n%s", self.metadata_df.head())
        logger.info(
            "Finished loading dataset %s, available dataframes are : %s",
            self.config.label, self.available_datasets
        )
        self.check_expectations()

    def check_expectations(self):
        # conditions should be a subset of the metadata corresponding column
        if not set(self.config.conditions).issubset(
                set(self.metadata_df["condition"].unique())):
            logger.error("Conditions %s are not a subset of "
                         "the metadata declared conditions",
                         self.config.conditions)
            raise ValueError(
                f"Conditions {self.config.conditions} are not a subset "
                f"of the metadata declared conditions"
            )

        verify_metadata_sample_not_duplicated(self.metadata_df)
        message_bad_separator_input(self.metadata_df, "metadata")
        message_bad_separator_input(self.abundances_df, "abundances")
        message_bad_separator_input(self.mean_enrichment_df, "enrichment")
        message_bad_separator_input(self.isotopologue_proportions_df,
                                    "isotopologue proportions")
        message_bad_separator_input(self.isotopologues_df, "isotopologues")

    def split_datafiles_by_compartment(self) -> None:
        frames_dict = {}
        for data_file_label in self.available_datasets:
            if 'metadata' in data_file_label:
                continue
            dataframe_label = data_file_label + "_df"  # TODO: this is fragile!
            tmp_co_dict = df_to_dict_by_compartment(
                getattr(self, dataframe_label),
                self.metadata_df)  # split by compartment
            frames_dict[data_file_label] = tmp_co_dict

        frames_dict = drop_all_nan_metabolites_on_comp_frames(
            frames_dict, self.metadata_df)
        frames_dict = set_samples_names(frames_dict, self.metadata_df)
        self.compartmentalized_dfs = frames_dict

    def get_file_for_label(self, label):
        if label == "abundances":
            return self.config.abundances
        elif label == "mean_enrichment":
            return self.config.mean_enrichment
        elif label == "isotopologue_proportions":
            return self.config.isotopologue_proportions
        elif label == "isotopologues":
            return self.config.isotopologues
        else:
            raise ValueError(f"Unknown label {label}")

    def get_files_extension_as_dict(self):
        """returns dictionary of file extensions, uses extfind (helpers)"""
        extension_dict: Dict[str, str] = dict()
        extension_dict['metadata'] = extfind(self.sub_folder_absolute,
                                             self.config.metadata)
        extension_dict['abundances'] = extfind(self.sub_folder_absolute,
                                               self.config.abundances)
        extension_dict['mean_enrichment'] = extfind(
            self.sub_folder_absolute, self.config.mean_enrichment)
        extension_dict['isotopologues'] = extfind(self.sub_folder_absolute,
                                                  self.config.isotopologues)
        extension_dict['isotopologue_proportions'] = extfind(
            self.sub_folder_absolute,
            self.config.isotopologue_proportions)
        return extension_dict


class DataIntegrationConfig(DatasetConfig):
    transcripts: ListConfig
    pathways: DictConfig


class DataIntegration(Dataset):
    config: DataIntegrationConfig
    deg_dfs: Dict[int, pd.DataFrame] = {}
    pathways_dfs: Dict[str, pd.DataFrame] = {}

    def set_dataset_integration_config(self):
        self.preload()
        self.split_datafiles_by_compartment()

        self.check_expectations()  # of the Dataset class
        self.check_expectations_integration_data()  # of this child class

    def check_expectations_integration_data(self):
        if not len(set(self.config.transcripts)) == \
               len(self.config.transcripts):
            logger.error(
                f"Duplicated names: {self.config.transcripts}"
                f" in transcripts")
            raise ValueError(
                f"Duplicated names: {self.config.transcripts}"
                f" in transcripts")

        if not len(set(molecular_types_for_metabologram).difference(
                set(self.config.pathways.keys())
        )) == 0:
            logger.error("Unrecognized pathways configuration"
                         " in dataset yaml file")
            raise ValueError(
                "Unrecognized pathways_file_names configuration"
                " in dataset yaml file"
            )

    def load_deg_dfs(self):
        # generate dictionary of transcripts dataframes (DEGs) :
        # the keys are integers, with the order of files in the dataset yml
        for i, file_name in enumerate(self.config.transcripts):
            try:
                file_extension = extfind(self.sub_folder_absolute, file_name)
                path_deg_file = os.path.join(
                    self.sub_folder_absolute,
                    f"{file_name}.{file_extension}")
                deg_df = pd.read_csv(path_deg_file, sep='\t', header=0)
                self.deg_dfs[i] = deg_df
            except FileNotFoundError:
                logger.info(f"{file_name}.{file_extension}: file not found")
            except Exception as e:
                logger.info(
                    f'Error while opening file {file_name}.{file_extension} '
                    f' \n {e}')

        logger.info("Finished loading transcripts dataframes: "
                    "%s", self.config.transcripts)

    def load_pathways_dfs(self):
        for k in self.config.pathways.keys():
            try:
                file_extension = extfind(
                    self.sub_folder_absolute, self.config.pathways[k])
                path_file = os.path.join(
                    self.sub_folder_absolute,
                    f"{self.config.pathways[k]}.{file_extension}")
                pathway_df = pd.read_csv(path_file, sep='\t', header=0)
                self.pathways_dfs[k] = pathway_df
            except FileNotFoundError:
                logger.info(
                    f"{self.config.pathways[k]}.{file_extension}: not found")
            except Exception as e:
                logger.info(
                    f'{e}. Error while opening file '
                    f'{self.config.pathways[k]}.{file_extension} \n {e}')

        logger.info("Finished loading pathways dataframes: "
                    "%s", self.config.pathways)

    def get_names_transcripts_files(self) -> Dict[int, str]:
        out_dict = {}
        for i, name in enumerate(self.config.transcripts):
            out_dict[i] = name
        return out_dict
