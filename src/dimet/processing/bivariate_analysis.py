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


def save_output(result, compartment, dataset, file_name, more_string,test,
                out_table_dir):

    result = result.sort_values(["padj"],
                                ascending=True)
    #comp = "-".join(map(lambda x: "-".join(x), comparison))
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


def bivariate_compute_this(
        file_name, df: pd.DataFrame, metadata_df,
        compartment: str,
        dataset: Dataset,
        cfg: DictConfig,
        comparison: List[str], behavior: str,
        test: str, out_table_dir:str) -> None:
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
    df_dict = compute_bivariate_fashion(
         df, metadata_df, compartment, comparison, behavior, test
    )

    for akey in df_dict.keys():
        df = df_dict[akey]
        df.set_index("metabolite", inplace=True)
        result = compute_padj(df, 0.05,
                           cfg.analysis.method.correction_method)
        more_string = f"{comparison}-{akey}"
        if akey == "ok":
            more_string = f"{comparison}"
        save_output(result, compartment, dataset, file_name,
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


def compute_test_for_df_dict(df_dict, test):
    """xxxx"""
    for akey in df_dict.keys():
        df = df_dict[akey].copy()
        df.index = df['metabolite']
        print(df.head())
        stat_list = []
        pvalue_list = []
        for i, metabolite in enumerate(list(df['metabolite'])):
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


def compute_bivariate_fashion(df, metadata_df,
                              compartment, # TODO define if necesary here ?
                              comparison, behavior, test
      ) -> pd.DataFrame:
    # format the data
    if behavior == "conditions_comparison_time_blocks":
        # as stated in config, only abundances or mean enrichment processed
        df_dict = way_1(df, metadata_df, comparison)

    elif behavior == "conditions_MDV_comparison":
        df_dict = way_2() # separately for each time point)

    elif behavior == "timepoints_MDV_comparison":
        df_dict = way_3() # separately for each condition

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
        df = row_wise_nanstd_reduction(df)
        df = df.round(decimals=6)

        automatic_comparisons = set_comparisons_by_behavior(
            behavior, cfg, metadata_df_subset
        )

        for comparison in automatic_comparisons:
            bivariate_compute_this(
                file_name,
                df, metadata_df_subset,
                compartment, dataset,
                cfg, comparison,
                behavior,  # specifies the mode, coheret with the df
                test, out_table_dir)






#
#
# def run_bivariate_test(df: pd.DataFrame, comparison: List,  # TODO fix argunemtns
#                          test: str) -> pd.DataFrame:
#     """
#     This is a switch function for computing statistics for a pairwise
#     differential analysis
#     The comparison is a list with 2 sublists that contain column names
#     """
#     metabolites = df.index.values
#     correlation_coeff_array = []
#     pval = []
#
#     for i in df.index:  # i is one  metabolite name
#         a1 = np.array(df.loc[i, comparison[0]], dtype=float)
#         a2 = np.array(df.loc[i, comparison[1]], dtype=float)
#         vInterest = a1[~np.isnan(a1)]
#         vBaseline = a2[~np.isnan(a2)]
#
#         if (len(vInterest) < 2) | (len(vBaseline) < 2):
#             return pd.DataFrame(
#                 data={
#                     "metabolite": metabolites,
#                     "stat": [float("nan")] * len(metabolites),
#                     "pvalue": [float("nan")] * len(metabolites),
#                 }
#             )
#
#         if test == "pearson":
#             stat_result, pval_result = 0
#
#         elif test == "spearman":
#             stat_result, pval_result = 0
#
#         correlation_coeff_array.append(stat_result)
#         pval.append(pval_result)
#
#     assert len(metabolites) == len(pval)
#     return pd.DataFrame(
#         data={
#             "metabolite": metabolites,
#             "pvalue": pval,
#         }
#     )