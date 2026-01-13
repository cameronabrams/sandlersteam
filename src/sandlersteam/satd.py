# Author: Cameron F Abrams <cfa22@drexel.edu>

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from importlib.resources import files

def merge_high_low_P_tables(lowP_df, highP_df):
    """ Merge two dataframes of steam table data at low and high pressures """
    lowP_df['P'] = lowP_df['P'] / 1000.0  # convert kPa to MPa
    merged_df = pd.concat([lowP_df, highP_df], axis=0)
    icolumn = merged_df.columns[0]
    merged_df.sort_values(by=icolumn, inplace=True)
    return merged_df

class svi:
    # wrap the interp1d function so that it returns a scalar
    def __init__(self, f):
        self.f = f
    def __call__(self, x):
        return self.f(x).item()

class SaturatedSteamTables:
    """ Saturated steam tables based on Sandler's Steam Tables
        Data from:
        Sandler, S. I. (2017). Chemical, biochemical, and engineering
        thermodynamics 5th ed. John Wiley & Sons.
    """
    data_path = files('sandlersteam') / 'resources' / 'data'
    tablesP = [data_path / 'SandlerSatdSteamTableP1.txt', 
               data_path / 'SandlerSatdSteamTableP2.txt']
    tablesT = [data_path / 'SandlerSatdSteamTableT1.txt', 
               data_path / 'SandlerSatdSteamTableT2.txt'] 

    # _sp = ['T', 'P', 'vL', 'vV', 'uL', 'uV', 'hL', 'hV', 'sL', 'sV']
    # _su = ['C','MPa','m3/kg','m3/kg','kJ/kg','kJ/kg','kJ/kg','kJ/kg','kJ/kg-K','kJ/kg-K']
    # _sfs = ['{: .1f}','{: .2f}','{: .6g}','{: .6g}','{: .6g}','{: .6g}','{: .6g}','{: .6g}','{: .6g}','{: .6g}']
    colorder = ['P', 'T', 'vL', 'vV', 'uL', 'uV', 'hL', 'hV', 'sL', 'sV']
    def __init__(self):
        self.DF = {'P': merge_high_low_P_tables(
                        pd.read_csv(self.tablesP[0], sep=r'\s+', header=0, index_col=None),
                        pd.read_csv(self.tablesP[1], sep=r'\s+', header=0, index_col=None))[self.colorder],
                   'T': merge_high_low_P_tables(
                        pd.read_csv(self.tablesT[0], sep=r'\s+', header=0, index_col=None),
                        pd.read_csv(self.tablesT[1], sep=r'\s+', header=0, index_col=None))[self.colorder]}
        self.lim = {'P':[self.DF['P']['P'].min(),self.DF['P']['P'].max()],
                    'T':[self.DF['T']['T'].min(),self.DF['T']['T'].max()]}
        for p in ['vL', 'vV', 'uL', 'uV', 'hL', 'hV', 'sL', 'sV']:
            self.lim[p] = [min(self.DF['P'][p].min(), self.DF['T'][p].min()),
                           max(self.DF['P'][p].max(), self.DF['T'][p].max())]
        self.interpolators = {}
        for bp, cp in zip(['P', 'T'], ['T', 'P']):
            self.interpolators[bp] = {}
            X = np.array(self.DF[bp][bp].to_list())
            for p in [cp, 'vL', 'vV', 'uL', 'uV', 'hL', 'hV', 'sL', 'sV']:
                Y = np.array(self.DF[bp][p].to_list())
                self.interpolators[bp][p] = svi(interp1d(X, Y, kind='linear'))
