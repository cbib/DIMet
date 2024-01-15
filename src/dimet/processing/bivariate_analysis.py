#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Johanna Galvis, Florian Specque, Macha Nikolski
"""
import logging
import operator
import os
from functools import reduce
from typing import List, Dict

import numpy as np
import pandas as pd
import scipy.stats as stats
from dimet.constants import (assert_literal, availtest_methods_type,
                             data_files_keys_type)
from dimet.data import Dataset
from dimet.processing.differential_analysis import time_course_auto_list_comparisons
from dimet.helpers import (arg_repl_zero2value,
                           compute_padj,
                           row_wise_nanstd_reduction,
                           #split_rows_by_threshold
                           )
from omegaconf import DictConfig

logger = logging.getLogger(__name__)


def save_output(result_df, compartment, dataset, file_name, more_string,
                test, out_table_dir):

    result = result_df.sort_values(["padj"],
                                ascending=True)
    result['compartment'] = compartment

    base_file_name = dataset.get_file_for_label(file_name)
    base_file_name += f"--{compartment}-{more_string}-{test}"
    output_file_name = os.path.join(out_table_dir,
                                    f"{base_file_name}.tsv")
    result.to_csv(
        output_file_name,
        index_label="metabolite",
        header=True,
        sep="\t",
    )
    logger.info(f"Saved the result in {output_file_name}")


def bivariate_compute_current_comparison(
        file_name, df: pd.DataFrame, metadata_df,
        compartment: str,
        dataset: Dataset,
        cfg: DictConfig,
        comparison: List[str], behavior: str,
        test: str, out_table_dir: str) -> None:
    """
    Runs a bivariate analysis for blocks of values, handling two behaviors:
    a. comparing conditions, time block:
      - time course blocks, using total metabolite abundances
      - time course blocks, using mean enrichment
    b. comparing conditions, MDV block:
      -  MDV blocks, using isotopologue proportions
    c. comparing consecutive time points # TODO: not yet implemented
      - MDV blocks, using isotopologue proportions

    """
    df_dict = compute_bivariate_by_case(
         df, metadata_df, compartment, comparison, behavior, test
    )

    for akey in df_dict.keys():
        df = df_dict[akey]

        df.set_index("metabolite", inplace=True)
        result_df = compute_padj(df, 0.05,
                           cfg.analysis.method.correction_method)
        comparison_str = "-".join(comparison)
        more_string = f"{comparison_str}-{akey}"
        if akey == "ok":
            more_string = f"{comparison_str}"
        save_output(result_df, compartment, dataset, file_name,
                    more_string, test, out_table_dir)


def compute_time_ordered_gmeans(sorted_time_numeric, metadata_df_sorted,
                  curr_condition, row: pd.Series) -> List[float]:
    time_ordered_gmeans = []
    for k in sorted_time_numeric:
        samples_names = metadata_df_sorted.loc[
            (metadata_df_sorted['condition'] == curr_condition) &
            (metadata_df_sorted["timenum"] == k), 'name_to_plot'
        ].values
        x = list(row[samples_names].dropna())
        try:
            time_ordered_gmeans.append(np.around(
                stats.gmean(np.array(x)), decimals=6))
        except ValueError:
            time_ordered_gmeans.append(np.nan)
        except Exception as e:
            logger.info(e)
            time_ordered_gmeans.append(np.nan)

    return time_ordered_gmeans


def way_1(df, metadata_df , comparison) -> Dict[str, pd.DataFrame]:
    """
    Using mean enrichment or abundances,
    computes the arrays of geometric means, ordered by time.
    Outputs dict of df of arrays, for comparing 2 condition time blocks.
    """
    metadata_df = metadata_df.loc[metadata_df['condition'].isin(comparison), :]
    metadata_df['timenum'] = metadata_df['timenum'].astype(float)
    metadata_df_sorted = metadata_df.sort_values(by="timenum")

    df = df[metadata_df_sorted['name_to_plot']]
    sorted_time_numeric = sorted(metadata_df['timenum'].unique())
    init_val = [np.array([0 for j in range(len(sorted_time_numeric))])
                   for i in range(df.shape[0])]
    metabolites = list(df.index)
    tmp_dict: dict = {"metabolite": metabolites,
                      "arr_1": init_val.copy(), "arr_2": init_val.copy()}
    for i, metabolite in enumerate(metabolites):
        # time ordered array of values for this current condition
        row = df.loc[metabolite, :]
        for j, curr_condition in enumerate(comparison):
            time_ordered_gmeans = compute_time_ordered_gmeans(
                sorted_time_numeric, metadata_df_sorted,
                curr_condition, row
            )
            location = f'arr_{j + 1}'
            tmp_dict[location][i] = time_ordered_gmeans

    df_dict : dict = {"ok": pd.DataFrame(tmp_dict)}

    return df_dict


def compute_isotopologue_meaning(isotopologues_list):
    isotopologues_uniq = sorted(list(np.unique(np.array(isotopologues_list))))
    clue_isotopologue = pd.DataFrame({
        "isotopologue_name": isotopologues_uniq
    })
    clue_isotopologue[["metabolite", "m+x"]] = (
        clue_isotopologue["isotopologue_name"].str.split(
         "_m+", expand=True, regex=False))
    clue_isotopologue["m+x"] = clue_isotopologue["m+x"].astype(int)
    clue_isotopologue = clue_isotopologue.sort_values(
        by=["metabolite", "m+x"])
    return clue_isotopologue


def way_2(df, metadata_df, comparison):
    """
        Using isotopologue proportions,
        computes the arrays of geometric means, ordered by m+x.
        Outputs dict of df of arrays, for comparing 2 condition time blocks.
        (a separated computation is done by time-point)
    """
    df_dict = dict()
    clue_isotopologue = compute_isotopologue_meaning(list(df.index))
    metabolites_uniq = clue_isotopologue["metabolite"].unique()
    df_dict = dict()
    for timepoint in list(metadata_df["timepoint"].unique()):
        tmp_dict: dict = {"metabolite": metabolites_uniq,
        "arr_1" : list(),  "arr_2": list()}
        metadata_df_tp = metadata_df.loc[
                         metadata_df["timepoint"] == timepoint, :]
        df_timepoint = df[metadata_df_tp['name_to_plot']].copy()
        for j, condition in enumerate(comparison):
            metadata_df_tp_cond = metadata_df_tp.loc[
                                  metadata_df["condition"] == condition, :]
            df_compar = df_timepoint[
                metadata_df_tp_cond['name_to_plot']].copy()
            df_compar["gmean"] = np.around(df_compar.apply(
                lambda x: stats.gmean(x.dropna()), axis = 1), decimals=6)
            df_compar["isotopologue_name"] = list(df_compar.index)
            merged_df = pd.merge(df_compar, clue_isotopologue, how='left',
                                 on="isotopologue_name").sort_values(
                                 by=["metabolite", "m+x"])
            MDV_ordered_list_gmeans = list()
            for i, metabolite in enumerate(metabolites_uniq):
                mdv_arr = merged_df.loc[merged_df["metabolite"] == metabolite, 
                "gmean"].values
                MDV_ordered_list_gmeans.append(np.array(mdv_arr))

            location = f'arr_{j + 1}'
            tmp_dict[location] = MDV_ordered_list_gmeans

        df_dict[timepoint] =  pd.DataFrame(tmp_dict)

    return df_dict


def compute_test_for_df_dict(df_dict, test):
    """parses a dictionary of dataframes
       computes bivariate test for each dataframe"""
    for akey in df_dict.keys():
        df = df_dict[akey].copy()

        df.index = df['metabolite']
        stat_list = []
        pvalue_list = []
        for i, metabolite in enumerate(list(df['metabolite'])):
            # array of n-(timepoints or m+x) geometrical means values
            array_1 = df.loc[metabolite, "arr_1"]
            array_2 = df.loc[metabolite, "arr_2"]
            if test == "pearson":
                stat_res, pvalue = stats.pearsonr(array_1, array_2)

            elif test == "spearman":
                stat_res, pvalue = stats.spearmanr(array_1, array_2)

            stat_list.append(stat_res)
            pvalue_list.append(pvalue)

        df["correlation_coefficient"] = stat_list
        df["pvalue"] = pvalue_list
        df_dict[akey] = df.copy()

    return df_dict


def compute_bivariate_by_case(df, metadata_df,
                              compartment,  # TODO define if necesary here ?
                              comparison, behavior, test
                              ) -> pd.DataFrame:
    # format the data
    if behavior == "conditions_comparison_time_blocks":
        # as stated in config, only abundances or mean enrichment processed
        df_dict = way_1(df, metadata_df, comparison)

    elif behavior == "conditions_MDV_comparison":
        df_dict = way_2(df, metadata_df, comparison) # separately for each time point)

    elif behavior == "timepoints_MDV_comparison":
        import sys
        sys.exit()
        df_dict = way_3(df, metadata_df, comparison) # separately for each condition

    # call the computation
    df_dict = compute_test_for_df_dict(df_dict, test)
    # output the df

    return df_dict





def round_result_float_columns(df: pd.DataFrame) -> pd.DataFrame:
    result_float_columns = [
        "correlation_coefficient",
        "pvalue",
        "padj"]

    columns_gmean = [column for column in list(df.columns) if
                     column.startswith("gmean_")]  # also gmean columns
    result_float_columns = list(
        set(result_float_columns).union(set(columns_gmean)))

    for column in result_float_columns:
        if column in list(df.columns):
            df[column] = np.around(
                df[column].astype(float).to_numpy(),
                decimals=6
            )

    return df



def conditions_to_comparisons(conditions: List[str]) -> List[List[str]]:
    """
    for the bi-variate analysis
    builds a list of lists, where each array is a 1 to 1 comparison,
    of the elements of comparisons or timepoints list
    example:
    comparisons_list = [A, B, C]
    output : [[A, B], [A, C], [B, C]]
    Note: order is not relevant for bi-variate analysis
    """
    comparisons = list()
    for i in conditions:
        for j in conditions:
            if i != j:
                pair = np.sort(np.array([i, j]), axis=None)
                if list(pair) not in comparisons:
                    comparisons.append(list(pair))
    return comparisons


def set_comparisons_by_behavior(behavior: str, cfg: DictConfig,
                                metadata_df: pd.DataFrame) -> List[List[str]]:
    """
    for the bi-variate analysis
    builds a list of lists, where each array is a 1 to 1 comparison,
    whether for conditions or consecutive time points
    """
    if behavior in ["conditions_comparison_time_blocks",
                    "conditions_MDV_comparison"]:
        comparisons = conditions_to_comparisons(cfg.analysis.conditions)
    elif behavior == "timepoints_MDV_comparison":
        comparisons = time_course_auto_list_comparisons(
            metadata_df
        )
    return comparisons


def bi_variate_analysis(
        file_name: data_files_keys_type, dataset: Dataset, cfg: DictConfig,
        behavior: str, out_table_dir: str
) -> None:
    """
    Bi-variate analysis is performed on compartmentalized versions
    of data files
    Attention: we replace zero values using the provided method
    Writes the table with computed statistics in the relevant output directory
    """
    assert_literal(file_name, data_files_keys_type, "file name")
    assert behavior in ["conditions_comparison_time_blocks",
                        "conditions_MDV_comparison",
                        "timepoints_MDV_comparison"], "wrong behavior chosen"
    impute_value = cfg.analysis.method.impute_values[file_name]
    test = cfg.analysis.method[behavior][file_name]  # e.g. pearson

    for compartment, compartmentalized_df in \
            dataset.compartmentalized_dfs[file_name].items():
        df = compartmentalized_df
        metadata_df_subset = dataset.metadata_df.loc[
            dataset.metadata_df['compartment'] == compartment, :]
        metadata_df_subset = metadata_df_subset.loc[
            metadata_df_subset['condition'].isin(cfg.analysis.conditions), :]

        df = df[metadata_df_subset['name_to_plot']]

        val_instead_zero = arg_repl_zero2value(impute_value, df)
        df = df.replace(to_replace=0, value=val_instead_zero)
        # note: do not drop rows all zero or all nan, blocks can break !
        if file_name == "abundances":  # only reduce if abundances
            df = row_wise_nanstd_reduction(df)
        df = df.round(decimals=6)

        automatic_comparisons = set_comparisons_by_behavior(
            behavior, cfg, metadata_df_subset
        )

        for comparison in automatic_comparisons:
            bivariate_compute_current_comparison(
                file_name,
                df, metadata_df_subset,
                compartment, dataset,
                cfg, comparison,
                behavior,  # specifies the mode, coheret with the df
                test, out_table_dir)

