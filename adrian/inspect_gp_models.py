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


# convert this to a python script
# load in gp model
# load in mangled csv for each object
# plot the gp model iteration and the shiftedflux values
# ask if gp model looks acceptable
# if yes, append filename to good_models = []
# if no, skip
# repeat for all objects

def inspect_gp(foldername):
    files = glob.glob(foldername+"/*")
    
    with open('searched_files.txt', 'w') as file1:
        with open('good_searched_files.txt', 'w') as file2:
            file1.write("# SN, sample\n")
            file2.write("# SN, sample\n")
            file1.flush()
            file2.flush()

            for i,file in enumerate(files): 
                sn = file[18:-13]

                good_fits = {}
                
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
                
                plt.figure()
                
                for s in gp_df['sample'].unique(): 
                    smol_df = gp_df.loc[(gp_df['sample']==s)]

                    gb = (smol_df['filt']=='g')
                    rb = ~gb

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
                    plt.show()

                    like = input("Did you like that lightcurve: ")
                    print(like)
                    good_fits = []
                    if like.lower() in ["y", "yes"]:
                        print("good lc", f"{sn} sample {s}") 
                        good_fits.append(s)
                        file1.write(f"{sn},{s}" + '\n')
                        file2.write(f"{sn},{s}" + '\n')
                        file1.flush()
                        file2.flush()
                    elif like.lower() in ['quit','exit']:
                        file1.write(f"{sn},{s}" + '\n')
                        file1.flush()
                        return good_fits
                    else:
                        file1.write(f"{sn},{s}" + '\n')
                        file1.flush()
                        continue
    return good_fits

if __name__=="__main__":
    assert len(sys.argv) == 2, "missing folder name argument or you gave me too many arguments"
    good_stuff = inspect_gp(sys.argv[1])
    print(good_stuff)