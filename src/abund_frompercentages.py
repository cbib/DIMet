#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 12:58:52 2022

@author: johanna

Calculates and saves:
1. absolute abundances of isotopologues
using corrected isotopologues (%)

2. mean Enrichment as in Castano-cerezo, Bellvert et al, 2019

--
Note: Some platforms do not supply corrected absolute abundances,
but they do supply:
   a) corrected isotopologues (%)
   b) raw abundances
so this script is dedicated to data preparation in that case,
by it obtaining a proxy to corrected abundances using a) and b)

"""


from .fun_fm import *





def saveabundance(df, nameout):
    df.to_csv(nameout, sep="\t")



def calculate_meanEnri(
                datadi,
                tableAbund,
                tableIC,
                metadata,
                names_compartments,
                namesuffix,
                odir ):
    """
    output : two tables
    - calc_corrAbu_  :  "corrected" abundances, calculated
    - meanEnrich_ :  weighted mean as in Castano-cerezo, Bellvert et al, 2019.
    """
    if not os.path.exists(odir):
        os.makedirs(odir)

    for co in names_compartments.values():
        abun = pd.read_csv(
            datadi + tableAbund + "_" + namesuffix + "_" + co + ".tsv",
            sep="\t",
            index_col=0,
        )
        # note that pandas automatically transform any 99.9% in decimal 0.999

        propc = pd.read_csv(
            datadi + tableIC + "_" + namesuffix + "_" + co + ".tsv",
            sep="\t",
            index_col=0,
        )
        dicos = dict()
        dicos[co] = {}
        dicos[co]["metadata"] = metadata.loc[metadata.short_comp == co]
        selecols = dicos[co]["metadata"]["sample"]
        # order  columns (samples) , same order in both df:
        dicos[co][tableAbund] = abun[selecols]
        dicos[co][tableIC] = propc[selecols]
        prop_df = dicos[co][tableIC]
        prop_df = correction_prop_df(prop_df)
        isotosrowdata = yieldrowdataB(prop_df)  # this funciton is from fun_fm

        isotosrowdata = isotosrowdata.assign(
                            coefficient=isotosrowdata["m+x"].str.replace("m+", "", regex=False))
        isotosrowdata.coefficient = isotosrowdata.coefficient.astype('int')

        outdf = dicos[co][tableAbund].copy()
        proxy_correctabu = prop_df.copy()
        metabolites = outdf.index

        for met_i in metabolites:
            isos_rowdata_i = isotosrowdata.loc[isotosrowdata["metabolite"] == met_i, :]
            prop_df_i = prop_df.loc[isos_rowdata_i.isotopolgFull.tolist(), :]
            abu_i = dicos[co][tableAbund].loc[met_i, :]
            abus_isotops_i = prop_df_i.multiply(abu_i.T) / 100
            # take advantage of this corrected  abus_isotops_i, to save into proxy_...
            proxy_correctabu.loc[abus_isotops_i.index, :] = abus_isotops_i

            # prepare to apply formula (sum(coeffs*abus_isotops_i)) / n
            # from Castano-cerezo, Bellvert et al, 2019.
            nmax = max(isos_rowdata_i.coefficient)
            abu_coef_prod_df = abus_isotops_i.multiply(isos_rowdata_i.coefficient.tolist(),
                                                       axis="index")
            abu_coef_sum_df = abu_coef_prod_df.sum(axis="index")
            meanenrich_i = abu_coef_sum_df / nmax
            meanenrich_i = meanenrich_i[selecols]
            outdf.loc[met_i, ] = meanenrich_i

        proxyfile = odir + "calc_corrAbu_" + namesuffix + "_" + co + ".tsv"
        proxy_correctabu.to_csv(proxyfile, sep='\t', header =True, index=True)
        print("\nSaved calc_corrAbu_ file, to tmp/. Compartment:", co, "\n")

        ofile = odir + "meanEnrich_" + namesuffix + "_" + co + ".tsv"
        outdf.to_csv(ofile, sep='\t', header =True, index=True)
        print("\nSaved mean enrichment to tmp/. Compartment:", co, "\n")
    return 0


def split_mspecies_files(dirtmpdata, names_compartments, namesuffix,
               odir):
    if not os.path.exists(odir):
        os.makedirs(odir)

    for co in names_compartments.values():
        proxyfile = dirtmpdata + "calc_corrAbu_" + namesuffix + "_" + co + ".tsv"
        proxy_correctabu = pd.read_csv(proxyfile, sep='\t', index_col=0, header=0)
        isotosrowdata =  yieldrowdataB(proxy_correctabu)

        metabolites = isotosrowdata.metabolite.unique().tolist()
        samples = proxy_correctabu.columns.tolist()

        proxy_correctabu['isotopolgFull'] = proxy_correctabu.index
        workingdf = proxy_correctabu.merge(isotosrowdata, on="isotopolgFull")
        workingdf.index = proxy_correctabu.index


        # save total marked : sum only marked ( m+0 excluded in totmk)
        preoutmat = np.zeros(shape =(len(metabolites), len(samples)))
        outdf = pd.DataFrame(preoutmat)
        outdf.index = metabolites
        outdf.columns = samples
        # total marked : send to mktot files
        preptotmark = workingdf.loc[workingdf["m+x"] != "m+0", :]

        for met_i in metabolites:
           subi  = preptotmark.loc[preptotmark["metabolite"] == met_i ,: ]
           outdf.loc[met_i,] = subi[samples].sum(axis=0)
        nameout = f"{odir}abux_byProp_totmk_{co}.tsv"
        saveabundance(outdf, nameout)

        # splitting and saving dfs by m+x

        #print(max([int(i) for i in workingdf["m+x"].str.replace("m+","")]))
        maxmk = max([int(i) for i in
                     workingdf["m+x"].str.replace("m+", "", regex = False)])

        for k in range(maxmk+1):
            selmk = f"m+{k}"
            tmp = workingdf.loc[workingdf["m+x"] == selmk,:]
            if tmp.shape[0] > 0 :
                tmp = tmp.drop(columns=isotosrowdata.columns)
                rw = [i.split("_m+")[0] for i in tmp.index.tolist()]
                tmp.index = rw
                nameout = f"{odir}abux_byProp_{selmk}_{co}.tsv"
                saveabundance(tmp, nameout)
    return 0


# -----------------------------------------------------------------
# NOTE isotosrowdata is a df, guiding  metabolites&isotopologues grouping:
# 		  metabolite   m+x      isotopolgFull
# 	0          Hexose   m+0         Hexose_m+0
# 	1          Hexose   m+1         Hexose_m+1
#
# 	166  Stearic_acid  m+16  Stearic_acid_m+16
# 	167  Stearic_acid  m+17  Stearic_acid_m+17
# 	168  Stearic_acid  m+18  Stearic_acid_m+18
# -----------------------------------------------------------------

# def makematch_abund_rowdata(abund, isotosrowdata):
#     """
#     check lines that are in abundance but not in corrected %
#
#     """
#     set(abund.index) - set(isotosrowdata["metabolite"])
#     rtodrop = list(set(abund.index) - set(isotosrowdata["metabolite"]))
#     abund = abund.drop(rtodrop, axis=0)
#     return abund