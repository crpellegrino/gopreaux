import sys
from caat import CAAT, DataCube, SN, SNCollection, GP, GP3D, SNModel
from caat.utils import WLE
import json
import matplotlib.pyplot as plt
import numpy as np
import logging
# from sklearn.gaussian_process.kernels import Matern, RBF, WhiteKernel
import glob
import pandas as pd
import re
from collections import defaultdict



def samples2dict(filepath: str) -> dict:
    """
    translates txt file of format SN,sample into dictionary
    """
    result = defaultdict(list)
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                sn, sample = line.strip().split(",")
                result[sn].append(int(sample))
    return dict(result)


def inspect_gp(foldername, good_files):
    files = glob.glob(foldername+"/*")
    good_files = samples2dict(good_files)
    
    for i,file in enumerate(files): #iterating over all GP models which is fine b/c every model has at least one good realization/sample to plot
        sn = file[18:-13]
        
        # read in GP model 
        gp_df = pd.read_csv(file)

        # read in actual photometry points
        data = pd.read_csv(f'../data/SESNe/SNIIb/{sn}/{sn}_datacube_mangled.csv')
        datadf = data.loc[(data['Nondetection']==False)&((data['Filter']=='g')|(data['Filter']=='r'))&
                        (data['Phase']>-20)&(data['Phase']<50),
                        ['Phase', 'Filter', 'ShiftedMag', 'Magerr']] 
            #changing df column names to match expected names of fit_photometry function
        datadf.rename(columns={'Magerr': 'MagErr'}, inplace=True)

        if sn in ['SN2022qzr', 'SN2021pb', 'SN2020sbw', 'SN2020rsc', 'SN2020ikq']:
            datadf['Mag'] = -1*datadf['ShiftedMag']
        else:
            datadf.rename(columns={'ShiftedMag':'Mag'}, inplace=True)

        gband = (datadf['Filter']=='g')
        rband = ~gband
        
        for s in good_files[sn]: 
            smol_df = gp_df.loc[(gp_df['sample']==s)]

            gb = (smol_df['filt']=='g')
            rb = ~gb

            plt.figure()

            plt.plot(smol_df.loc[gb, 'time'], smol_df.loc[gb, 'flux'], color='cyan')
            plt.plot(smol_df.loc[rb, 'time'], smol_df.loc[rb, 'flux'], color='orange')

            plt.scatter(datadf.loc[gband,['Phase']].values, datadf.loc[gband,['Mag']].values, 
                        facecolor='cyan', edgecolor='k', label='g', s=15)
            plt.scatter(datadf.loc[rband,['Phase']].values, datadf.loc[rband,['Mag']].values, 
                        facecolor='orange', edgecolor='k', label='r', s=15)
        
            plt.xlabel('Phase [days]')
            plt.ylabel('Flux')
            plt.title(f"{sn} - sample {s}")
            plt.legend()
            plt.savefig(f'temp_plots/{sn}_sample_{s}.png', dpi=200)
            # plt.show()
    return

if __name__=="__main__":
    assert len(sys.argv) == 3, "missing folder name argument or you gave me too many arguments"
    inspect_gp(sys.argv[1], sys.argv[2])