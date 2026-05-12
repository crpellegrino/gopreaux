import sys
# from caat import CAAT, DataCube, SN, SNCollection, GP, GP3D, SNModel
# from caat.utils import WLE
import json
import matplotlib.pyplot as plt
import numpy as np
import logging
# from sklearn.gaussian_process.kernels import Matern, RBF, WhiteKernel
import glob
import pandas as pd
import re
from collections import defaultdict
import scipy as sp
from itertools import compress

def sample_with_min_spacing(df, n=20, min_gap=3.0, phase_col='time', seed=None):
    """
    Randomly select n rows such that no two selected points
    are within min_gap days of each other in phase.
    """
    rng = np.random.default_rng(seed)

    # Shuffle a copy to ensure random (not phase-ordered) selection
    shuffled_df = df.sample(frac=1, random_state=rng.integers(1e9)).reset_index(drop=True)

    selected = [] #list of single-row series
    selected_phases = []

    for _, row in shuffled_df.iterrows():
        p = row[phase_col]
        # Accept this point only if it's far enough from all already-selected phases
        if all(abs(p - sp) >= min_gap for sp in selected_phases):
            selected.append(row)
            selected_phases.append(p)
        if len(selected) == n:
            break

    if len(selected) < n:
        print(f"Warning: only found {len(selected)} points satisfying the constraint.")

    # randomly dropping ~10% of observations
    random_keep = np.array([sp.stats.binom.rvs(1,0.9) for i in range(len(selected))])
    final_selection = list(compress(selected, random_keep))

    return pd.DataFrame(final_selection).reset_index(drop=True)


def samples2dict(filepath: str, errs=False) -> dict:
    """
    translates txt file of format SN,sample into dictionary
    """
    result = defaultdict(list)
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            elif errs:
                sn,gerr,rerr = line.strip().split(",")
                result[sn].append(float(gerr))
                result[sn].append(float(rerr))
            else:
                sn, sample = line.strip().split(",")
                result[sn].append(int(sample))
    return dict(result)


def inspect_photometry(file, refit_file):
    """
    randomly samples photometry from GP fit and overplots onto GP
    for visual inspection of fake photometry / LC
    """
    
    # getting the sample number of the good GP fits for each SN 
    good_inds = samples2dict(file)
    sn_names = list(good_inds.keys())

    avg_errs = samples2dict('avg_lc_err.txt',errs=True) # format SN:[gerr, rerr]
    
    with open(refit_file, 'a') as file1:
        #randomly drawing samples from the good fits
        for sn in sn_names:
            #read in the GP model fit
            gp_data = pd.read_csv(f'nsample_GP_models/{sn}_GP_model.csv')
            g_inds = (gp_data['filt']=='g')
            r_inds = (gp_data['filt']=='r')

            #read in the original photometry datacube
            data = pd.read_csv(f'../data/SESNe/SNIIb/{sn}/{sn}_datacube_mangled.csv')
            datadf = data.loc[(data['Nondetection']==False)&((data['Filter']=='g')|(data['Filter']=='r'))&(data['Phase']>min(gp_data['time']))&(data['Phase']<max(gp_data['time'])),
                            ['Phase', 'Filter', 'ShiftedMag', 'Magerr']]
            
            #changing df column names to match expected names of fit_photometry function
            datadf.rename(columns={'Magerr': 'MagErr'}, inplace=True)
            # datadf.rename(columns={'ShiftedMag':'Mag'}, inplace=True)

            if sn in ['SN2022qzr', 'SN2021pb', 'SN2020sbw', 'SN2020rsc', 'SN2020ikq']: #some objs have inverse mags, so flipping them back
                datadf['Mag'] = -1*datadf['ShiftedMag']
            else:
                datadf.rename(columns={'ShiftedMag':'Mag'}, inplace=True)

            for ind in good_inds[sn]:
                # print(sn, ind)
                g_df = sample_with_min_spacing(gp_data.loc[g_inds&(gp_data['sample']==ind)], n=20, min_gap=2.0, phase_col='time', seed=9203)
                r_df = sample_with_min_spacing(gp_data.loc[r_inds&(gp_data['sample']==ind)], n=20, min_gap=2.0, phase_col='time', seed=15003)
                
                #random add mag-noise by sampling from uncertainty
                g_df['flux_rand'] = np.random.normal(g_df['flux'], avg_errs[sn][0]*2, len(g_df['flux']))
                r_df['flux_rand'] = np.random.normal(r_df['flux'], avg_errs[sn][1]*2, len(r_df['flux']))

                #append avg-err as col to the df's
                g_df['flux_err'] = np.ones(len(g_df))*avg_errs[sn][0]
                r_df['flux_err'] = np.ones(len(r_df))*avg_errs[sn][1]

                #create plot
                plt.figure()
                  #plotting the interpolated photometry
                plt.errorbar(g_df['time'], g_df['flux_rand'], g_df['flux_err'], ls='', marker='.', markerfacecolor='cyan', markeredgecolor='k', ecolor='cyan', markersize=10)
                plt.errorbar(r_df['time'], r_df['flux_rand'], r_df['flux_err'], ls='', marker='.', markerfacecolor='orange', markeredgecolor='k', ecolor='orange', markersize=10)
                #   #plotting the original photometry,
                # plt.errorbar(datadf.loc[datadf['Filter']=='g', 'Phase'], datadf.loc[datadf['Filter']=='g', 'Mag'], datadf.loc[datadf['Filter']=='g', 'MagErr'], 
                #              ls='', marker='.', color='teal', alpha=0.4, markersize=7,zorder=0)
                # plt.errorbar(datadf.loc[datadf['Filter']=='r', 'Phase'], datadf.loc[datadf['Filter']=='r', 'Mag'], datadf.loc[datadf['Filter']=='r', 'MagErr'], 
                #              ls='', marker='.', color='tomato', alpha=0.4,  markersize=7,zorder=0)
                  #plotting the GP model
                plt.plot(gp_data.loc[g_inds&(gp_data['sample']==ind),'time'], gp_data.loc[g_inds&(gp_data['sample']==ind),'flux'], color='cyan', zorder=1, alpha=0.1) #0.1
                plt.plot(gp_data.loc[r_inds&(gp_data['sample']==ind),'time'], gp_data.loc[r_inds&(gp_data['sample']==ind),'flux'], color='orange', zorder=1, alpha=0.1)
                
                plt.title(f"{sn}-sample {ind}")
                plt.xlabel('phase [days]')
                plt.ylabel('flux relative to peak')
                # plt.savefig(f'temp_plots/rand_added/{sn}_sample_{ind}_fake_LC_2sig.png', dpi=200)
                plt.show()

                #saving the sub-selected fake LCs that are good 
                #noting which need to be refit
                like = input("Did you like that fake lightcurve: ")
                if like.lower() in ["y", "yes"]:
                    #save the g_df and r_df as one dataframe and as a csv
                    print("good fake lc", f"{sn} sample {ind}")
                    combo_df = pd.concat([g_df,r_df],ignore_index=True)
                    combo_df.to_csv(f'fake_LCs/rand_added/{sn}_sample{ind}_fake_LC.csv', index=False)
                elif like.lower() in ['quit','exit']:
                    file1.write(f"{sn},{ind}" + '\n')
                    file1.flush()
                    return
                else:
                    #save the SN and sample number to txt file for refitting
                    file1.write(f"{sn},{ind}" + '\n')
                    file1.flush()
                    continue
    return 

if __name__=="__main__":
    assert len(sys.argv) == 3, "missing filename argument or you gave me too many arguments"
    inspect_photometry(sys.argv[1], sys.argv[2])
