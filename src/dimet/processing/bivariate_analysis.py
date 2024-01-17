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
from dimet.helpers import (arg_repl_zero2value,
                           compute_padj,
                           row_wise_nanstd_reduction
                           )
from omegaconf import DictConfig

logger = logging.getLogger(__name__)


def compute_test_for_df_dict(df_dict, test: str) -> Dict[str, pd.DataFrame]:
    """parses a dictionary of dataframes
       computes bivariate test : e.g. pearson
       row-wise for each dataframe"""
    for akey in df_dict.keys():
        df = df_dict[akey].copy()

        df.index = df['metabolite']
        stat_list = []
        pvalue_list = []
        for i, metabolite in enumerate(list(df['metabolite'])):
            # array of n-(timepoints or m+x) geometrical means values
            array_1 = df.loc[metabolite, "gmean_arr_1"]
            array_2 = df.loc[metabolite, "gmean_arr_2"]
            if test == "pearson":
                try:
                    stat_res, pvalue = stats.pearsonr(array_1, array_2)
                except ValueError:
                    stat_res, pvalue = (np.nan, np.nan)

            elif test == "spearman":
                try:
                    stat_res, pvalue = stats.spearmanr(array_1, array_2)
                except ValueError:
                    stat_res, pvalue = (np.nan, np.nan)

            stat_list.append(stat_res)
            pvalue_list.append(pvalue)

        df["correlation_coefficient"] = stat_list
        df["pvalue"] = pvalue_list
        df_dict[akey] = df.copy()

    return df_dict


def compute_bivariate_by_behavior(df, metadata_df, comparison,
                                  behavior: str, test: str) -> pd.DataFrame:
    """
    performs two steps:
    1. calls functions to compute geometric means, obtaining df's in dict
    2. computes the bivariate statistical test (pearson by default)
    """
    if behavior == "conditions_MDV_comparison":
        df_dict = conditions_MDV_gmean_df_dict(df, metadata_df, comparison)

    elif behavior == "timepoints_MDV_comparison":
        df_dict = timepoints_MDV_gmean_df_dict(df, metadata_df, comparison)

    elif behavior == "conditions_metabolite_time_profiles":
        # only abundances or mean enrichment processed
        df_dict = metabolite_time_profiles_gmean_df_dict(
            df, metadata_df, comparison)

    # 2. compute the statistical test
    df_dict = compute_test_for_df_dict(df_dict, test)

    return df_dict


def conditions_MDV_gmean_df_dict(df, metadata_df, comparison: List[str]):
    """
    Note: e.g. comparison [Ctl, Treated1]
    Separately by time-point:
    Using isotopologue proportions, by metabolite and condition,
     computes the arrays of geometric means, ordered by m+x.
     Outputs dict of df, each df with columns arr_1 and arr_2 (np arrays);
             the keys of the dict are the time-points
    """
    clue_isotopologue_df = compute_isotopologue_meaning(list(df.index))
    metabolites_uniq = clue_isotopologue_df["metabolite"].unique()
    df_dict = dict()
    for timepoint in list(metadata_df["timepoint"].unique()):
        inner_gmean_dict = {"metabolite": metabolites_uniq,
                            "gmean_arr_1": list(),  "gmean_arr_2": list()}
        metadata_df_tp = metadata_df.loc[
                         metadata_df["timepoint"] == timepoint, :]
        df_timepoint = df[metadata_df_tp['name_to_plot']].copy()
        for k, condition in enumerate(comparison):
            metadata_df_tp_cond = metadata_df_tp.loc[
                                  metadata_df_tp["condition"] == condition, :]
            df_compar = df_timepoint[
                metadata_df_tp_cond['name_to_plot']].copy()

            inner_gmean_dict = inner_gmean_dict_filler(
                k, inner_gmean_dict, df_compar, clue_isotopologue_df)

        df_dict[timepoint] = pd.DataFrame(inner_gmean_dict)

    return df_dict


def timepoints_MDV_gmean_df_dict(df, metadata_df, comparison: List[str]):
    """
    Note: e.g. comparison [T1, T0]
    Separately by condition:
     Using isotopologue proportions, by metabolite and time-point,
     computes the arrays of geometric means, ordered by m+x.
     Outputs dict of df, the keys of the dict are the conditions
    """
    clue_isotopologue_df = compute_isotopologue_meaning(list(df.index))
    metabolites_uniq = clue_isotopologue_df["metabolite"].unique()
    df_dict = dict()
    for condition in list(metadata_df['condition'].unique()):
        metadata_df_sub1 = metadata_df.loc[
                         metadata_df["condition"] == condition, :]

        df_sub1 = df[metadata_df_sub1['name_to_plot']].copy()
        inner_gmean_dict = {"metabolite": metabolites_uniq,
                            "gmean_arr_1": list(), "gmean_arr_2": list()}

        for k, timepoint in enumerate(comparison):
            metadata_df_sub2 = metadata_df_sub1.loc[
                               metadata_df_sub1["timepoint"] == timepoint, :]
            df_sub2 = df_sub1[metadata_df_sub2['name_to_plot']].copy()
            inner_gmean_dict = inner_gmean_dict_filler(
                k, inner_gmean_dict, df_sub2, clue_isotopologue_df)

        df_dict[condition] = pd.DataFrame(inner_gmean_dict)

    return df_dict


def inner_gmean_dict_filler(k: int, inner_gmean_dict, df_sub2,
                            clue_isotopologue_df) -> Dict[str, list]:
    """
    For MDV, fills the dictionary (of one single condition or time point)
    inner_gmean_dict = {"metabolite": metabolites_uniq,
                          "gmean_arr_1": list(), "gmean_arr_2": list()}, where
    arr_1 represents the first str of comparison, and arr_2 the second.
    Each 'gmean_arr_-' list has as many np.arrays as metabolites; each 'gmean_arr_-' list
    respects the order of the metabolites in the 'metabolite' key.
    And each inner array respects the order of isotopologues m+x integer
    """
    assert k in [0, 1], "k can only take value 0 or 1"
    metabolites_uniq = clue_isotopologue_df["metabolite"].unique()
    # compute the arrays of geometric means
    df_sub2["gmean"] = np.around(df_sub2.apply(
        lambda x: stats.gmean(x.dropna()), axis=1), decimals=6)
    df_sub2["isotopologue_name"] = list(df_sub2.index)
    # order these geometric means by m+x
    merged_df = pd.merge(df_sub2, clue_isotopologue_df, how='left',
                         on="isotopologue_name").sort_values(
                         by=["metabolite", "m+x"])
    MDV_ordered_list_gmeans = list()
    for i, metabolite in enumerate(metabolites_uniq):
        mdv_arr = merged_df.loc[merged_df["metabolite"] == metabolite,
        "gmean"].values
        MDV_ordered_list_gmeans.append(np.array(mdv_arr))
    key_name = f'gmean_arr_{k + 1}'
    inner_gmean_dict[key_name] = MDV_ordered_list_gmeans

    return inner_gmean_dict


def metabolite_time_profiles_gmean_df_dict(
        df, metadata_df, comparison: List[str]) -> Dict[str, pd.DataFrame]:
    """
    Using mean enrichment or abundances,
    computes the arrays of geometric means, ordered by time.
    Outputs dict of df of arrays, for comparing 2 condition time profiles.
    """
    metadata_df = metadata_df.loc[
                  metadata_df['condition'].isin(comparison), :]
    metadata_df['timenum'] = metadata_df['timenum'].astype(float)
    metadata_df_sorted = metadata_df.sort_values(by="timenum")
    df = df[metadata_df_sorted['name_to_plot']]
    metabolites = list(df.index)
    inner_gmean_dict = {"metabolite": metabolites,
                        "gmean_arr_1": list(), "gmean_arr_2": list()}
    for k, condition in enumerate(comparison):
        metadata_sorted2 = metadata_df_sorted.loc[
                metadata_df_sorted["condition"] == condition, :].copy()
        tmp_gmean_time_dict = {}
        for curr_time in list(metadata_sorted2['timepoint'].unique()):
            metadata_sorted3 = metadata_sorted2.loc[
                metadata_sorted2['timepoint'] == curr_time, :].copy()
            data_time_df = df[metadata_sorted3['name_to_plot']].copy()
            data_time_df = data_time_df.assign(
                gmean=data_time_df.apply(lambda x: stats.gmean(x.dropna()),
                                         axis=1))
            # todo: add sanity check: if not at least 2 no NaN samples, replace gmean value by NaN
            tmp_gmean_time_dict[curr_time] = data_time_df["gmean"].tolist()
        tmp_df = pd.DataFrame(tmp_gmean_time_dict)
        #  tmp_df :   T0       T2h  ...
        # 0  0.471528  0.719920  ...
        # 1  3.692766  1.743812   ...
        key_name = f'gmean_arr_{k + 1}'
        inner_gmean_dict[key_name] = list(np.around(tmp_df.to_numpy(), 6))

    return {'metabo_time_profile': pd.DataFrame(inner_gmean_dict)}


def compute_isotopologue_meaning(isotopologues_list) -> pd.DataFrame:
    """output: df :
    isotopologue_name  metabolite  m+x
    Cit_m+0            Cit         0
    ...
    """
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


def conditions_to_comparisons(conditions: List[str]) -> List[List[str]]:
    """
    for the bi-variate analysis generates comparisons ONLY FOR CONDITIONS
    example:  input: conditions = [A, B, C]
    output : [[A, B], [A, C], [B, C]]
    The order of the comparison in each inner list, e.g. [A, B] or [B, A]
    does not affect the result of the bi-variate test (verified)
    """
    comparisons = list()
    for i in conditions:
        for j in conditions:
            if i != j:
                pair = np.sort(np.array([i, j]), axis=None)
                if list(pair) not in comparisons:
                    comparisons.append(list(pair))
    return comparisons


def timepoints_to_comparisons(metadata_df: pd.DataFrame) -> List[List[str]]:
    """generates comparisons specifically as consecutive time points
    output: e.g. [[T0, T1], [T1, T2]]"""
    # order the time by the numeric part
    time_df = metadata_df[["timepoint", "timenum"]].copy().drop_duplicates()
    time_df["timenum"] = time_df["timenum"].astype(float)
    time_df = time_df.sort_values(by='timenum')
    comparisons = list()
    for i in range(time_df.shape[0] - 1):
        j = i + 1
        comparisons.append([
            time_df["timepoint"].tolist()[j], time_df["timepoint"].tolist()[i]])

    return comparisons


def set_comparisons_by_behavior(behavior: str, cfg: DictConfig,
                                metadata_df: pd.DataFrame) -> List[List[str]]:
    """
    For the bi-variate analysis, builds a list of lists, each list is a
    comparison. The 'behavior' determines if conditions (1) or time-points (2)
    Output:
        e.g.(1): comparisons = [['Ctl', 'Treated1'], ['Ctl', 'Treated2']]
        or e.g.(2): comparisons = [[T1, T0], [T2, T1]]
    """
    if behavior in ["conditions_metabolite_time_profiles",
                    "conditions_MDV_comparison"]:
        comparisons = conditions_to_comparisons(cfg.analysis.conditions)
    elif behavior == "timepoints_MDV_comparison":
        comparisons = timepoints_to_comparisons(metadata_df)
    return comparisons


def save_output(result_df, compartment, dataset, file_name, out_file_name_str,
                test, out_table_dir) -> None:
    """saves result to tab delimited file, for one comparison"""
    result = result_df.sort_values(["padj"],
                                ascending=True)
    result['compartment'] = compartment

    base_file_name = dataset.get_file_for_label(file_name)
    base_file_name += f"--{compartment}-{out_file_name_str}-{test}"
    output_file_name = os.path.join(out_table_dir,
                                    f"{base_file_name}.tsv")
    result.to_csv(
        output_file_name,
        index_label="metabolite",
        header=True,
        sep="\t",
    )
    logger.info(f"Saved the result in {output_file_name}")


def bivariate_run_and_save_current_comparison(
        file_name, df: pd.DataFrame, metadata_df,
        compartment: str, dataset: Dataset,
        cfg: DictConfig, comparison: List[str],
        behavior: str, test: str, out_table_dir: str) -> None:
    """
    Runs a bivariate analysis for blocks of values, handling 3 behavior types:
    a - conditions_MDV_comparison:
        comparison of the MDV arrays between 2 user specified conditions (separately for each time point)
    b - timepoints_MDV_comparison:
        comparison of the MDV arrays between 2 consecutive time points (separately for each condition)
    c - conditions_metabolite_time_profiles:
        comparison of the time course profiles of the metabolites
        total abundances and mean enrichment, between two conditions
    """
    df_dict = compute_bivariate_by_behavior(
         df, metadata_df, comparison, behavior, test
    )
    for akey in df_dict.keys():
        df = df_dict[akey]

        df.set_index("metabolite", inplace=True)
        result_df = compute_padj(df, 0.05,
                                 cfg.analysis.method.correction_method)
        comparison_str = "-".join(comparison)
        out_file_name_str = f"-MDV-{comparison_str}--{akey}"
        if akey == "metabo_time_profile":
            out_file_name_str = f"{comparison_str}"
        save_output(result_df, compartment, dataset, file_name,
                    out_file_name_str, test, out_table_dir)


def bivariate_comparison(
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
    assert behavior in ["conditions_metabolite_time_profiles",
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
        if file_name == "abundances":  # only reduction values if abundances
            df = row_wise_nanstd_reduction(df)
        df = df.round(decimals=6)

        automatic_comparisons = set_comparisons_by_behavior(
            behavior, cfg, metadata_df_subset
        )
        for comparison in automatic_comparisons:
            # e.g. comparison = ['A', 'B'], list of exactly two elements
            bivariate_run_and_save_current_comparison(
                file_name, df, metadata_df_subset, compartment,
                dataset, cfg, comparison,
                behavior,  # specifies the type of comparison
                test, out_table_dir)

