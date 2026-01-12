# Author: Cameron F. Abrams <cfa22@drexel.edu>

import numpy as np
import pandas as pd
from importlib.resources import files
import io

"""
                                                                                                    1         1
          1         2         3         4         5         6         7         8         9         0         1
012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123    
 260  0.0012749  1127.9  1134.3   2.8830  0.0012645  1121.1   1133.7   2.8699  0.0012550  1114.6   1133.4   2.8576
(0,4),(6,15),    (17,24),(25,32), (34,40),(42,51),   (53,60), (62,69), (71,77),(79,88),(90,97),(99,106),(108,114)
"""
def my_split(data, hder, P, Tsat, fixw=False):
    ndfs = []
    with io.StringIO(data) as f:
        if fixw:
            df = pd.read_fwf(f, colspecs=((0,4),(6,15), (17,24),(25,32),(34,40),(42,51), (53,60),(62,69),(71,77),(79,88),(90,97),(99,106),(108,114)), header=None, index_col=None)
        else:
            df = pd.read_csv(f, sep=r'\s+', header=None, index_col=None)
        df.columns = hder
        i = 1
        for p, ts in zip(P, Tsat):
            ndf = pd.DataFrame({'T': df['T'].copy(), 'P': np.array([p for _ in range(df.shape[0])])})
            if ndf.iloc[0, 0] == 'Sat.':
                ndf.iloc[0, 0] = ts
            ndf['T'] = ndf['T'].astype(float)
            tdf = df.iloc[:, i:i+4].copy()
            i += 4
            ndf = pd.concat((ndf, tdf), axis=1)
            ndf.dropna(axis=0, inplace=True)
            ndf.sort_values(by='T', inplace=True)
            ndfs.append(ndf)
    mdf = pd.concat(ndfs, axis=0)
    return mdf

class UnsaturatedSteamTable:
    data_path = files('sandlersteam') / 'resources' / 'data'
    table_suph = data_path / 'SandlerSuphSteamTables.txt'
    table_subc = data_path / 'SandlerSubcSteamTables.txt'
    _p = ['T', 'P', 'u', 'v', 's', 'h', 'x']
    _u = ['C', 'MPa', 'kJ/kg', 'm3/kg', 'kJ/kg-K', 'kJ/kg', '']
    _fs = ['{: .1f}','{: .2f}','{: .6g}','{: .6g}','{: .6g}','{: .6g}','{: .2f}']

    def __init__(self, phase: str = 'V'):
        if phase == 'V':
            with open(self.table_suph,'r') as f:
                lines = f.read().split('\n')
        elif phase == 'L':
            with open(self.table_subc,'r') as f:
                lines = f.read().split('\n')
        # identify header
        hder = lines[0].split()
        # identify lines with pressures
        plines = []
        # save list of pressure values
        Pval = []
        for i, l in enumerate(lines):
            if 'P = ' in l:
                plines.append(i)
        # extract and format all blocks as concurrent dataframes
        DFS = []
        for i, (l, r) in enumerate(zip(plines[:-1], plines[1:])):
            # get pressures, and saturation temperatures, if present
            tokens = lines[l].split()
            kills = ['P', '=', 'MPa']
            for k in kills:
                while k in tokens:
                    tokens.remove(k)
            P = []
            Tsat = []
            for x in tokens:
                if x[0] == '(':
                    y = x[1:-1]
                    Tsat.append(float(y))
                else:
                    P.append(float(x))
            while len(Tsat) < len(P):
                Tsat.append(None)
            Pval.extend(P)
            data = '\n'.join(lines[l + 1 : r])
            ndf = my_split(data, hder, P, Tsat, fixw=(phase == 'L'))
            DFS.append(ndf)
        tokens = lines[plines[-1]].split()
        kills = ['P', '=', 'MPa']
        for k in kills:
            while k in tokens:
                tokens.remove(k)
        P = []
        Tsat = []
        for x in tokens:
            if x[0] == '(':
                y = x[1:-1]
                Tsat.append(float(y))
            else:
                P.append(float(x))
        while len(Tsat) < len(P):
            Tsat.append(None)
        Pval.extend(P)
        data = '\n'.join(lines[plines[-1] + 1 :])
        ndf = my_split(data, hder, P, Tsat, fixw=(phase == 'L'))
        DFS.append(ndf)
        self.data = pd.concat(DFS, ignore_index=True)
        dof = self.data.columns
        self.uniqs = {}
        for d in dof:
            self.uniqs[d] = np.sort(np.array(list(set(self.data[d].to_list()))))

    def TPBilinear(self, specdict: dict):
        """
        Bilinear interpolation given T and P
        """
        retdict = {}
        xn, yn = specdict.keys()
        assert [xn, yn] == ['T', 'P']
        xi, yi = specdict.values()
        df: pd.DataFrame = self.data
        dof = df.columns
        retdict = {}
        retdict['T'] = xi
        retdict['P'] = yi
        tdf = df[df['P'] == yi]
        if not tdf.empty:
            X = np.array(tdf['T'])
            for d in dof:
                if d not in ['T','P']:
                    Y = np.array(tdf[d])
                    retdict[d] = np.interp(xi, X, Y, left=np.nan, right=np.nan)
        else:
            for PL, PR in zip(self.uniqs['P'][:-1], self.uniqs['P'][1:]):
                if PL < yi < PR:
                    break
            else:
                raise Exception(f'P {yi} not between {self.uniqs["P"][0]} and {self.uniqs["P"][-1]} at {xi} C')
            X = np.array([PL, PR])
            LDF = df[df['P'] == PL]
            RDF = df[df['P'] == PR]
            LT = np.array(LDF['T'])
            RT = np.array(RDF['T'])
            CT = np.array([T for T in LT if T in RT])
            if xi in CT:
                for d in dof:
                    if d not in ['T','P']:
                        Y = np.array([LDF[LT == xi][d].values[0], RDF[RT == xi][d].values[0]])
                        retdict[d] = np.interp(yi, X, Y, left=np.nan, right=np.nan)
            else:
                for TL, TR in zip(CT[:-1], CT[1:]):
                    if TL < xi < TR:
                        break
                else:
                    raise Exception(f'T {xi} not between {CT[0]} and {CT[-1]} at {yi} MPa')
                LTDF = LDF[(LDF['T'] == TL) | (LDF['T'] == TR)].sort_values(by='T')
                RTDF = RDF[(RDF['T'] == TL) | (RDF['T'] == TR)].sort_values(by='T')
                iv = np.zeros(2)
                for p in dof:
                    if not p in ['T','P']:
                        for i in range(2):
                            Lp = LTDF[p].values[i]
                            Rp = RTDF[p].values[i]
                            Y = np.array([Lp, Rp])
                            iv[i] = np.interp(yi, X, Y, left=np.nan, right=np.nan)
                        retdict[p] = np.interp(xi, np.array([TL, TR]), iv)
        return retdict

    def TThBilinear(self,specdict):
        """
        Bilinear interpolation given T and another property 
        (v, u, s, h, but not P)
        """
        xn, yn = specdict.keys()
        assert xn == 'T'
        xi, yi = specdict.values()
        df = self.data
        dof = list(df.columns)
        dof.remove(yn)
        dof.remove(xn)
        retdict = {}
        retdict['T'] = xi
        retdict[yn] = yi
        LLdat = {}
        LLdat[yn] = []
        for d in dof:
            LLdat[d] = []
        for P in self.uniqs['P'][::-1]:  # VUSH properties decrease with increasing P
            tdf = df[df['P'] == P]
            X = np.array(tdf['T'])
            Y = np.array(tdf[yn])
            if Y.min() < yi < Y.max():
                LLdat['P'].append(P)
                for d in 'vush':
                    Y = np.array(tdf[d])
                    LLdat[d].append(np.interp(xi, X, Y, left=np.nan, right=np.nan))
        X = np.array(LLdat[yn])
        retdict['P'] = np.interp(yi, X, np.array(LLdat['P']))
        for d in 'vush':
            if d != yn:
                retdict[d] = np.interp(yi, X, LLdat[d])
        return retdict

    def PThBilinear(self,specdict):
        """
        Bilinear interpolation given P and another property 
        (v, u, s, h, but not T)
        """
        xn, yn = specdict.keys()
        xi, yi = specdict.values()
        assert xn == 'P'
        df = self.data
        dof = list(df.columns)
        dof.remove(yn)
        dof.remove(xn)
        retdict = {}
        retdict['T'] = xi
        retdict[yn] = yi
        if xi in self.uniqs['P']:
            tdf = df[df['P'] == xi]
            X = np.array(tdf[yn])
            for pp in dof:
                Y = np.array(tdf[pp])
                retdict[pp] = np.interp(yi, X, Y, left=np.nan, right=np.nan)
        else:
            for PL, PR in zip(self.uniqs['P'][:-1], self.uniqs['P'][1:]):
                if PL < xi < PR:
                    break
            else:
                raise Exception(f'Error: no two blocks bracket P={xi}')
            ldf, rdf = df[df['P'] == PL], df[df['P'] == PR]
            ldict, rdict = {}, {}
            ldict['P'], rdict['P'] = PL, PR
            ldict[yn], rdict[yn] = yi, yi
            for d, xdf in zip([ldict, rdict], [ldf, rdf]):
                X = xdf[yn]
                for pp in dof:
                    Y = np.array(xdf[pp])
                    d[pp] = np.interp(yi, X, Y, left=np.nan, right=np.nan)
            X = np.array([PL, PR])
            for pp in dof:
                Y = np.array([ldict[pp], rdict[pp]])
                retdict[pp] = np.interp(xi, X, Y, left=np.nan, right=np.nan)
        return retdict
    
    def ThThBilinear(self,specdict):
        """
        Bilinear interpolation given two properties from v, u, s, h 
        (not T or P)
        """
        xn, yn =specdict.keys()
        assert not xn in ['T','P'] and not yn in ['T','P']
        xi, yi = specdict.values()
        df = self.data
        dof = self.data.columns
        LLdat = {}
        for d in dof:
            if d not in ['T','P']:
                LLdat[d] = []
        LLdat['T'] = self.uniqs['T']
        for T in LLdat['T']:
            tdf = df[df['T'] == T]
            X = np.array(tdf[xn])
            for d in dof:
                if d!='T' and d!=xn:
                    Y = np.array(tdf[d])
                    if Y.min() < yi < Y.max():
                        LLdat[d] = np.interp(xi, X, Y, left=np.nan, right=np.nan)
        X = LLdat[yn]
        retdict = {}
        retdict[xn] = xi
        retdict[yn] = yi
        for d in dof:
            if d!=xn and d!=yn:
                retdict[d] = np.interp(yi, X, LLdat[d])
        return retdict

    def Bilinear(self, specdict):
        """
        General bilinear interpolation dispatcher
        """
        xn, yn = specdict.keys()
        if [xn, yn] == ['T', 'P']:
            return self.TPBilinear(specdict)
        elif xn == 'T':
            return self.TThBilinear(specdict)
        elif xn == 'P':
            return self.PThBilinear(specdict)
        else:
            return self.ThThBilinear(specdict)
