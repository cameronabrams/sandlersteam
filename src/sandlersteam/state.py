# Author: Cameron F. Abrams <cfa22@drexel.edu>

import numpy as np
from scipy.interpolate import interp1d
from .satd import SaturatedSteamTables
from .unsatd import UnsaturatedSteamTable
from sandlermisc.statereporter import StateReporter


SteamTables = dict(
    satd = SaturatedSteamTables(),
    suph = UnsaturatedSteamTable('V'),
    subc = UnsaturatedSteamTable('L')
)

def show_available_tables(args):
    print('Available steam tables:')
    print(f'  Saturated steam tables:')
    print(f'    T-sat: T from {SteamTables["satd"].lim["TC"][0]} to {SteamTables["satd"].lim["TC"][1]} C')
    print(f'             from {np.round(SteamTables["satd"].lim["TC"][0] + 273.15,2)} to {np.round(SteamTables["satd"].lim["TC"][1] + 273.15,2)} K')
    print(f'    P-sat: P from {SteamTables["satd"].lim["P"][0]} to {SteamTables["satd"].lim["P"][1]} MPa')
    print(f'             from {np.round(SteamTables["satd"].lim["P"][0]*10,2)} to {np.round(SteamTables["satd"].lim["P"][1]*10,2)} bar')
    print(f'  Superheated steam tables blocks:\nPressure (MPa) -> Temperatures (C):')
    for p in SteamTables["suph"].uniqs['P']:
        Tlist = SteamTables["suph"].data[SteamTables["suph"].data['P'] == p]['TC'].to_list()
        print(f'    {p:>5.2f} ->', ', '.join([f"{x:>7.2f}" for x in Tlist]))
    print(f'  Subcooled steam tables blocks:\nPressure (MPa) -> Temperatures (C):')
    for p in SteamTables["subc"].uniqs['P']:
        Tlist = SteamTables["subc"].data[SteamTables["subc"].data['P'] == p]['TC'].to_list()
        print(f'    {p:>5.2f} ->', ', '.join([f"{x:>6.2f}" for x in Tlist]))

def state_subcommand(args):
    state_kwargs = {}
    for p in State._p:
        val = getattr(args, p)
        if val is not None:
            state_kwargs[p] = val
    state = State(**state_kwargs)
    report = state.report()
    print(report)

class PHASE:
    def __init__(self):
        pass

LARGE=1.e99
NEGLARGE=-LARGE

G_PER_MOL = 18.01528  # g/mol for water
KG_PER_MOL = G_PER_MOL / 1000.000  # kg/mol for water
MOL_PER_KG = 1.0 / KG_PER_MOL  # mol/kg for water

class State:
    _p = ['TC','P','u','v','s','h','x']
    _u = ['C','MPa','kJ/kg','m3/kg','kJ/kg-K','kJ/kg','']
    _fs = ['{: .1f}','{: .2f}','{: .6g}','{: .6g}','{: .6g}','{: .6g}','{: .2f}']
    _sp = ['TC','P','VL','VV','UL','UV','HL','HV','SL','SV']
    _su = ['C','MPa','m3/kg','m3/kg','kJ/kg','kJ/kg','kJ/kg','kJ/kg','kJ/kg-K','kJ/kg-K']
    _sfs = ['{: .1f}','{: .2f}','{: .6g}','{: .6g}','{: .6g}','{: .6g}','{: .6g}','{: .6g}','{: .6g}','{: .6g}']

    def report(self):
        satd = hasattr(self, 'x') and self.x is not None
        msg = 'SATURATED ' if satd else 'UNSATURATED '
        reporter = StateReporter()
        for p, u, fs in zip(self._p, self._u, self._fs):
            if p == 'TC' and self.TC is not None:
                TK = self.TC + 273.15
                reporter.add_property('TC', self.TC, 'C', fstring=fs)
                reporter.add_value_to_property('TC', TK, 'K', fstring=fs)
            elif p != 'x':
                val = self.__dict__[p]
                if val is not None:
                    reporter.add_property(p, val, u, fstring=fs)
                    if p == 'P':
                        val_bar = val * 10.0  # convert MPa to bar
                        reporter.add_value_to_property('P', val_bar, 'bar', fstring=fs)
                    if p in 'hsuv':
                        val_spec = val * KG_PER_MOL
                        if 'kJ' in u:
                            val_spec = val_spec * 1000.0  # convert kJ to J
                            spec_u = u.replace('kJ/kg', 'J/mol')
                        else:
                            spec_u = u.replace('/kg', '/mol')
                        reporter.add_value_to_property(p, val_spec, spec_u, fstring=fs)
        if hasattr(self, 'x') and self.x is not None:
            reporter.add_property('x', self.x, 'kg vapor/kg', fstring=self._fs[-1])
            for sp, su, sfs in zip(self._p, self._u, self._sfs):
                if sp not in ['TC','P','x']:
                    valL = self.Liquid.__dict__[sp[0].lower()]
                    valV = self.Vapor.__dict__[sp[0].lower()]
                    reporter.add_property(f'{sp}L', valL, su, fstring=sfs)
                    reporter.add_property(f'{sp}V', valV, su, fstring=sfs)
                    if sp in 'hsuv':
                        valL_spec = valL * KG_PER_MOL
                        valV_spec = valV * KG_PER_MOL
                        if 'kJ' in su:
                            valL_spec = valL_spec * 1000.0  # convert kJ to J
                            valV_spec = valV_spec * 1000.0  # convert kJ to J
                            spec_u = su.replace('kJ/kg', 'J/mol')
                        else:
                            spec_u = su.replace('/kg', '/mol')
                        reporter.add_value_to_property(f'{sp}L', valL_spec, spec_u, fstring=sfs)
                        reporter.add_value_to_property(f'{sp}V', valV_spec, spec_u, fstring=sfs)
        return f'THERMODYNAMIC STATE OF {msg}STEAM/WATER:\n' + reporter.report()

    def _resolve(self):
        """ 
        Resolve the thermodynamic state of steam/water given specifications
        """
        self.spec = [p for p in self._p if self.__dict__[p]!=None]
        assert len(self.spec) == 2, f'Error: must specify two properties (of {self._p}) for steam; found {self.spec}'
        if self.spec[1] == 'x':
            ''' explicitly saturated '''
            self._resolve_satd()
        else:
            if self.spec == ['TC', 'P']:
                ''' T and P given explicitly '''
                self._resolve_subsup()
            elif 'TC' in self.spec or 'P' in self.spec:
                ''' T OR P given, along with some other property (v,u,s,h) '''
                self._resolve_TPTh()
            else:
                raise Exception('If not explicitly saturated, you must specify either TC or P')

    def _resolve_TPTh(self):
        ''' TC or P along with one other property (th) are specified '''
        p = self.spec[0]
        cp = 'P' if p == 'TC' else 'TC'
        th = self.spec[1]
        if self.satd.lim[p][0] < self.__dict__[p] < self.satd.lim[p][1]:
            ''' T or P is between saturation limits; may be a saturated state, so 
                check whether the second property value lies between its liquid
                and vapor phase values at this T or P '''
            thL = self.satd.interpolators[p][f'{th.upper()}L'](self.__dict__[p])
            thV = self.satd.interpolators[p][f'{th.upper()}V'](self.__dict__[p])
            self.__dict__[cp] = self.satd.interpolators[p][cp](self.__dict__[p])
            if thL < self.__dict__[th] < thV:
                ''' This is a saturated state! Use lever rule to get vapor fraction: '''
                self.x = (self.__dict__[th] - thL)/(thV - thL)
                self.Liquid = PHASE()
                self.Vapor = PHASE()
                self.Liquid.__dict__[th] = thL
                self.Vapor.__dict__[th] = thV
                for pp in self._sp:
                    if pp not in ['TC', 'P', f'{th.upper()}V', f'{th.upper()}L']:
                        ppp = self.satd.interpolators[p][pp](self.__dict__[p])
                        if pp[-1] =='V':
                            self.Vapor.__dict__[pp[0].lower()] = ppp
                        elif pp[-1] == 'L':
                            self.Liquid.__dict__[pp[0].lower()] = ppp
                for pp in self._p:
                    if pp not in [p, cp, 'x']:
                        self.__dict__[pp] = self.x * self.Vapor.__dict__[pp] + (1 - self.x) * self.Liquid.__dict__[pp]
            else:
                ''' even though T or P is between saturation limits, the other property is not '''
                specdict = {p.upper(): self.__dict__[p] for p in self.spec}
                if self.__dict__[th] < thL:
                    ''' Th is below its liquid-state value; assume this is a subcooled state '''
                    retdict = self.subc.Bilinear(specdict)
                    for p in self._p: 
                        if p not in self.spec and p != 'x':
                            self.__dict__[p] = retdict[p.upper()]
                else:
                    ''' Th is above its vapor-state value; assume this is a superheated state '''
                    retdict = self.suph.Bilinear(specdict)
                    for p in self._p: 
                        if p not in self.spec and p != 'x':
                            self.__dict__[p] = retdict[p.upper()]
        elif self.__dict__[p] > self.satd.lim[p][1]:
            ''' Th is above its vapor-state value; assume this is a superheated state '''
            specdict = {p.upper(): self.__dict__[p] for p in self.spec}
            retdict = self.suph.Bilinear(specdict)
            for p in self._p: 
                if p not in self.spec and p != 'x':
                    self.__dict__[p] = retdict[p.upper()]
        else:
            ''' Th is below its liquid-state value; assume this is a subcooled state '''
            specdict = {p.upper(): self.__dict__[p] for p in self.spec}
            retdict = self.subc.Bilinear(specdict)
            for p in self._p: 
                if p not in self.spec and p != 'x':
                    self.__dict__[p] = retdict[p.upper()]

    def _resolve_subsup(self):
        ''' T and P are both given explicitly.  Could be either superheated or subcooled state '''
        assert self.spec == ['TC', 'P']
        specdict = {'TC': self.TC, 'P': self.P}
        # print(f'at P {self.P}, checking T {self.T} between {self.satd.lim['TC'][0]} and {self.satd.lim['TC'][1]}')
        if self.satd.lim['TC'][0] < self.TC < self.satd.lim['TC'][1]:
            Psat = self.satd.interpolators['TC']['P'](self.TC)
            # print(f'Returns Psat of {Psat}')
        else:
            Psat = LARGE
        if self.P > Psat:
            ''' P is higher than saturation: this is a subcooled state '''
            retdict = self.subc.Bilinear(specdict)
        else:
            ''' P is lower than saturation: this is a superheated state '''
            retdict = self.suph.Bilinear(specdict)
        for p in self._p: 
            if p not in self.spec and p != 'x':
                self.__dict__[p] = retdict[p.upper()]

    def _resolve_satd(self):
        ''' This is an explicitly saturated state with vapor fraction (x) and one 
        other property (p) specified '''
        p = self.spec[0]
        self.Liquid = PHASE()
        self.Vapor = PHASE()
        if p == 'TC':
            ''' The other property is TC; make sure it lies between saturation limits '''
            if self.TC < self.satd.lim['TC'][0] or self.TC > self.satd.lim['TC'][1]:
                raise Exception(f'Cannot have a saturated state at TC = {self.TC} C')
            ''' Assign all other property values by interpolation '''
            for q in self._sp:
                if q != 'TC':
                    prop = self.satd.interpolators['TC'][q](self.TC)
                    if q == 'P': self.__dict__[q] = prop
                    if q[-1] == 'V':
                        self.Vapor.__dict__[q[0].lower()] = prop
                    elif q[-1] == 'L':
                        self.Liquid.__dict__[q[0].lower()] = prop
            for q in self._p:
                if not q in ['TC','P','x']:
                    self.__dict__[q] = self.x * self.Vapor.__dict__[q] + (1 - self.x) * self.Liquid.__dict__[q]
        elif p == 'P':
            ''' The other property is P; make sure it lies between saturation limits '''
            if self.P < self.satd.lim['P'][0] or self.P > self.satd.lim['P'][1]:
                raise Exception(f'Cannot have a saturated state at P = {self.P} MPa')
            ''' Assign all other property values by interpolation '''
            for q in self._sp:
                if q != 'P':
                    prop = self.satd.interpolators['P'][q](self.P)
                    if q == 'TC': self.__dict__[q] = prop
                    if q[-1] == 'V':
                        self.Vapor.__dict__[q[0].lower()] = prop
                    elif q[-1] == 'L':
                        self.Liquid.__dict__[q[0].lower()] = prop
            for q in self._p:
                if not q in ['TC','P','x']:
                    self.__dict__[q] = self.x * self.Vapor.__dict__[q] + (1 - self.x) * self.Liquid.__dict__[q]
        else:
            ''' The other property is neither T or P; must use a lever-rule-based interpolation '''
            self._resolve_satd_lever()

    def _resolve_satd_lever(self):
        p = self.spec[0]
        assert p != 'TC' and p != 'P'
        ''' Vapor fraction and one other property value (not T or P) is given '''
        th = self.__dict__[p]
        x = self.__dict__['x']
        ''' Build an array of V-L mixed properties based on given value of x '''
        Y = np.array(self.satd.DF['TC'][f'{p.upper()}V']) * x + np.array(self.satd.DF['TC'][f'{p.upper()}L']) * (1 - x)
        X = np.array(self.satd.DF['TC']['TC'])
        ''' define an interpolator '''
        f = interp1d(X, Y)
        try:
            ''' interpolate the Temperature '''
            self.T = f(x)
            ''' Assign all other property values '''
            for q in self._sp:
                if q != 'TC':
                    prop = self.satd.interpolators['TC'][q](self.TC)
                    if q == 'P': self.__dict__[q] = prop
                    if q[-1] == 'V':
                        self.Vapor.__dict__[q[0].lower()] = prop
                    elif q[-1] == 'L':
                        self.Liquid.__dict__[q[0].lower()] = prop
            for q in self._p:
                if not q in 'PTx':
                    self.__dict__[q] = self.x * self.Vapor.__dict__[q] + (1 - self.x) * self.Liquid.__dict__[q]
        except:
            raise Exception(f'Could not interpolate {p} = {th} at quality {x} from saturated steam table')

    def _scalarize(self):
        """ Convert all properties to scalars (not np.float64) """
        for p in self._p:
            val = self.__dict__[p]
            if isinstance(val, np.float64):
                self.__dict__[p] = val.item()
        if hasattr(self, 'x') and self.x is not None:
            if isinstance(self.x, np.float64):
                self.x = self.x.item()
        if hasattr(self, 'Liquid'):
            for k, p in self.Liquid.__dict__.items():
                if isinstance(p, np.float64):
                    self.Liquid.__dict__[k] = p.item()
        if hasattr(self, 'Vapor'):
            for k, p in self.Vapor.__dict__.items():
                if isinstance(p, np.float64):
                    self.Vapor.__dict__[k] = p.item()

    def __init__(self, **kwargs):
        self.satd = SteamTables['satd']
        self.suph = SteamTables['suph']
        self.subc = SteamTables['subc']
        for p in self._p:
            self.__dict__[p] = kwargs.get(p, None)
        self._resolve()
        self._scalarize()

class RandomSample(State):
    def __init__(self,phase='suph', satdDOF='TC', seed=None, Prange=None, Trange=None):
        if phase == 'satd':
            if satdDOF == 'TC':
                sample_this = SteamTables[phase].DF['TC']
            else:
                sample_this = SteamTables[phase].DF['P']
        elif phase == 'suph' or phase == 'subc':
            sample_this = SteamTables[phase].data
        abs_mins = {'TC': sample_this['TC'].min(), 'P': sample_this['P'].min()}
        abs_maxs = {'TC': sample_this['TC'].max(), 'P': sample_this['P'].max()}
        sample = sample_this.sample(n=1, random_state=seed)
        TC = sample['TC'].values[0]
        P = sample['P'].values[0]
        if Trange and not Prange:
            if Trange[0] < abs_mins['TC']:
                raise ValueError(f'Trange[0] ({Trange[0]}) is below the minimum TC in the data ({abs_mins["TC"]})')
            if Trange[1] > abs_maxs['TC']:
                raise ValueError(f'Trange[1] ({Trange[1]}) is above the maximum TC in the data ({abs_maxs["TC"]})')
            while not Trange[0] < TC < Trange[1]:
                sample = sample_this.sample(n=1)
                TC = sample['TC'].values[0]
        if Prange and not Trange:
            if Prange[0] < abs_mins['P']:
                raise ValueError(f'Prange[0] ({Prange[0]}) is below the minimum P in the data ({abs_mins["P"]})')
            if Prange[1] > abs_maxs['P']:
                raise ValueError(f'Prange[1] ({Prange[1]}) is above the maximum P in the data ({abs_maxs["P"]})')
            while not Prange[0] < P < Prange[1]:
                sample = sample_this.sample(n=1)
                P = sample['P'].values[0]
        if Prange and Trange:
            if Trange[0] < abs_mins['TC']:
                raise ValueError(f'Trange[0] ({Trange[0]}) is below the minimum TC in the data ({abs_mins["TC"]})')
            if Trange[1] > abs_maxs['TC']:
                raise ValueError(f'Trange[1] ({Trange[1]}) is above the maximum TC in the data ({abs_maxs["TC"]})')
            if Prange[0] < abs_mins['P']: 
                raise ValueError(f'Prange[0] ({Prange[0]}) is below the minimum P in the data ({abs_mins["P"]})')
            if Prange[1] > abs_maxs['P']: 
                raise ValueError(f'Prange[1] ({Prange[1]}) is above the maximum P in the data ({abs_maxs["P"]})')
            while not Prange[0] < P < Prange[1] or not Trange[0] < TC < Trange[1]:
                sample = sample_this.sample(n=1)
                TC = sample['TC'].values[0]
                P = sample['P'].values[0]
        if phase == 'satd':
            if satdDOF == 'TC':
                super().__init__(TC=TC, x=1.0)
            else:
                super().__init__(P=P, x=1.0)
        else:
            super().__init__(TC=TC, P=P)
    