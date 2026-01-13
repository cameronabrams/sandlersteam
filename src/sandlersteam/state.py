# Author: Cameron F. Abrams <cfa22@drexel.edu>
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field, fields
from scipy.interpolate import interp1d
from .satd import SaturatedSteamTables
from .unsatd import UnsaturatedSteamTables
from sandlermisc import statereporter

import logging

logger = logging.getLogger(__name__)

_instance = None

def get_tables():
    """Get or create the singleton tables instance."""
    global _instance
    if _instance is None:
        _instance = dict(
            satd = SaturatedSteamTables(),
            unsatd = UnsaturatedSteamTables()
        )
        _instance.update(
            suph = _instance['unsatd'].suph,
            subc = _instance['unsatd'].subc
        )
    return _instance

LARGE=1.e99
NEGLARGE=-LARGE

G_PER_MOL = 18.01528  # g/mol for water
KG_PER_MOL = G_PER_MOL / 1000.000  # kg/mol for water
MOL_PER_KG = 1.0 / KG_PER_MOL  # mol/kg for water
TCK = 647.3
TCC = TCK - 273.15
PCBAR = 221.20
PCMPA = PCBAR / 10.0

@dataclass
class State:
    
    _cache: dict = field(default_factory=dict, init=False, repr=False)
    """ Cache for computed properties """

    _input_state: dict | None = field(default=None, init=False, repr=False)
    """ Snapshot of input field values for cache validation """

    _STATE_VAR_ORDERED_FIELDS = ['T', 'P', 'v', 's', 'h', 'u']

    _STATE_VAR_FIELDS = frozenset(_STATE_VAR_ORDERED_FIELDS)
    """ Fields that define the input state for caching purposes """

    _UNIT_FIELDS = frozenset(['mass_unit', 'temperature_unit', 'pressure_unit', 'volume_unit'])
    """ Fields that define the units for caching purposes """

    def __setattr__(self, name, value):
        """ Clear cache on input field changes """
        if name in self._STATE_VAR_FIELDS.union(self._UNIT_FIELDS).union({'x'}) and hasattr(self, '_cache'):
            self._cache.clear()
        super().__setattr__(name, value)

    def _get_current_input_state(self):
        """ Snapshot of current input field values. """
        return {field: getattr(self, field) for field in self._STATE_VAR_FIELDS.union(self._UNIT_FIELDS).union({'x'})}

    def is_cache_stale(self, property_name: str = None):
        """
        Check if cache is stale.
        
        Parameters
        ----------
        property_name: str
            Specific property to check, or None for general staleness
        
        Returns
        -------
        bool or dict: True/False if property_name given, else dict of changed inputs
        """
        if property_name:
            return property_name not in self._cache
        
        # Show which inputs changed
        if self._input_state is None:
            return self._STATE_VAR_FIELDS.union(self._UNIT_FIELDS).union({'x'})  # All inputs "changed" (initial state)
        
        current = self._get_current_input_state()
        changed = {k for k in self._STATE_VAR_FIELDS.union(self._UNIT_FIELDS).union({'x'}) if current[k] != self._input_state.get(k)}
        return changed if changed else False

    # Units
    temperature_unit: str = 'C'
    pressure_unit: str = 'MPa'
    mass_unit: str = 'kg'
    volume_unit: str = 'm3'
    energy_unit: str = 'kJ'
    
    # State properties
    P: float = None
    T: float = None
    v: float = None
    u: float = None
    h: float = None
    s: float = None

    # Vapor fraction if two-phase saturated state    
    x: float = None
    L: State = None
    V: State = None

    def __post_init__(self):
        self.tables = get_tables()
        self.satd = self.tables['satd']
        self.suph = self.tables['suph']
        self.subc = self.tables['subc']

    def get_unit(self, field_name: str) -> str:
        """Get the unit for a given field"""
        unit_map = {
            'P': self.pressure_unit,
            'T': self.temperature_unit,
            'x': f'{self.mass_unit} vapor/{self.mass_unit} total',
            'v': f'{self.volume_unit}/{self.mass_unit}',
            'u': f'{self.energy_unit}/{self.mass_unit}',
            'h': f'{self.energy_unit}/{self.mass_unit}',
            's': f'{self.energy_unit}/{self.mass_unit}-K',
            'Pv': f'{self.energy_unit}/{self.mass_unit}',
        }
        return unit_map.get(field_name)
    
    def get_formatter(self, field_name: str) -> str:
        """Get the formatter for a given field"""
        formatter_map = {
            'P': '{: 5g}',
            'T': '{: 5g}',
            'x': '{: 5g}',
            'v': '{: 6g}',
            'u': '{: 6g}',
            'h': '{: 6g}',
            's': '{: 6g}',
            'Pv': '{: 6g}',
        }
        return formatter_map.get(field_name)

    def _to_table_units(self, x: float, units: str):
        """
        convert quantity x from specified units to table units (MPa, C, kJ/kg, m3/kg, kJ/kg-K)
        """
        match units:
            case 'C':
                return x
            case 'K':
                return x - 273.15
            case 'F':
                return (x - 32.0) * 5.0 / 9.0
            case 'MPa':
                return x
            case 'kPa':
                return x / 1000.0
            case 'bar':
                return x / 10.0
            case 'atm':
                return x / 10.1325
            case 'kJ/kg':
                return x
            case 'J/kg':
                return x / 1000.0
            case 'J/mol':
                return x / 1000.0 * MOL_PER_KG
            case 'm3/kg':
                return x
            case 'm3/mol':
                return x * MOL_PER_KG
            case 'kJ/kg-K':
                return x
            case 'J/kg-K':
                return x / 1000.0
            case 'J/mol-K':
                return x / 1000.0 * MOL_PER_KG
            case _:
                raise ValueError(f'Unsupported unit conversion from {units}')
    
    def _from_table_units(self, x: float, specified_unit: str):
        """ convert quantity x from table units (MPa, C, kJ/kg, m3/kg, kJ/kg-K) to specified units """
        match specified_unit:
            case 'C':
                return x
            case 'K':
                return x + 273.15
            case 'F':
                return x * 9.0 / 5.0 + 32.0
            case 'MPa':
                return x
            case 'kPa':
                return x * 1000.0
            case 'bar':
                return x * 10.0
            case 'atm':
                return x * 10.1325
            case 'kJ/kg':
                return x
            case 'J/kg':
                return x * 1000.0
            case 'J/mol':
                return x * 1000.0 / MOL_PER_KG
            case 'm3/kg':
                return x
            case 'm3/mol':
                return x / MOL_PER_KG
            case 'kJ/kg-K':
                return x
            case 'J/kg-K':
                return x * 1000.0
            case 'J/mol-K':
                return x * 1000.0 / MOL_PER_KG
            case _:
                raise ValueError(f'Unsupported unit conversion to {specified_unit}')

    @property
    def TC(self):
        """ Temperature in degrees Celsius """
        return self._to_table_units(self.T, self.temperature_unit)

    @property
    def PMPa(self):
        """ Pressure in MPa """
        return self._to_table_units(self.P, self.pressure_unit)

    @property
    def vm3KG(self):
        """ Specific volume in m3/kg """
        return self._to_table_units(self.v, f'{self.volume_unit}/{self.mass_unit}')

    @property
    def hkJKG(self):
        """ Specific enthalpy in kJ/kg """
        return self._to_table_units(self.h, f'{self.energy_unit}/{self.mass_unit}')
    
    @property
    def ukJKG(self):
        """ Specific internal energy in kJ/kg """
        return self._to_table_units(self.u, f'{self.energy_unit}/{self.mass_unit}')
    
    @property
    def skJKGK(self):
        """ Specific entropy in kJ/kg-K """
        return self._to_table_units(self.s, f'{self.energy_unit}/{self.mass_unit}-K')

    @property
    def Pv(self):
        """ Pressure * specific volume in kJ/kg """
        return self._from_table_units(self.PMPa * self.vm3KG * 1000.0, "kJ/kg")  # in kJ/kg

    def lookup(self):
        """ Lookup and compute all properties based on current inputs """
        if not 'lookup' in self._cache:
            self._cache['lookup'] = self._resolve()
        return self

    def _get_statespec(self):
        """ Get current state specification """
        vals = [getattr(self, p) for p in self._STATE_VAR_FIELDS]
        return_list = [p for p, x in zip(self._STATE_VAR_FIELDS, vals) if x is not None]
        logger.debug(f'Current state spec: {return_list} {vals}')

        if getattr(self, 'x') is not None:
            if (getattr(self, 'T') is None and getattr(self, 'P') is None):
                raise ValueError(f'Must specify T or P for an explicitly saturated state x={self.x:2g}')
            else:
                return_list.append('x')
                return return_list
        else: 
            # expect an unsaturated state       
            if sum(v is not None for v in vals) != 2:
                raise ValueError(f'Exactly two of {self._STATE_VAR_FIELDS} must be specified')
            return return_list 

    def _check_saturation(self, specs):
        """ Check if the specified state is saturated """
        if 'x' in specs:
            return True
        p = None
        op = None
        hasT = 'T' in specs
        hasP = 'P' in specs
        has_only_T_or_P = hasT ^ hasP
        if hasT and has_only_T_or_P:
            p = 'T'
            v = self.TC
            op = specs[0] if specs[1] == 'T' else specs[1]
            if v > TCC:
                logger.debug(f'Temperature {v}C exceeds critical temperature {TCC}C; cannot be saturated')
                return False
        elif hasP and has_only_T_or_P:
            p = 'P'
            v = self.PMPa
            op = specs[0] if specs[1] == 'P' else specs[1]
            if v > PCMPA:
                logger.debug(f'Pressure {v}MPa exceeds critical pressure {PCMPA} MPa; cannot be saturated')
                return False
        if p is not None and op is not None:
            logger.debug(f'Checking saturation at {p}={v} for property {op}={getattr(self, op)}')
            logger.debug(f'Saturation limits for {p}: {self.satd.lim[p]}')
            logger.debug(f'Between? {self.satd.lim[p][0] <= v <= self.satd.lim[p][1]}')
            if not (self.satd.lim[p][0] <= v <= self.satd.lim[p][1]):
                logger.debug(f'Out of saturation limits for {p}={v}')
                return False
            op_valV = self.satd.interpolators[p][f'{op}V'](v)
            op_valL = self.satd.interpolators[p][f'{op}L'](v)
            if op_valL < getattr(self, op) < op_valV or op_valV < getattr(self, op) < op_valL:
                return True
        return False

    def _resolve(self):
        """ 
        Resolve the thermodynamic state of steam/water given specifications
        """
        spec = self._get_statespec()

        if self._check_saturation(spec):
            self._resolve_satd(spec)
        else:
            self._resolve_unsatd(spec)
        self._scalarize()
        self._input_state = self._get_current_input_state()

    def report(self):
        reporter = statereporter.StateReporter()
        for p in self._STATE_VAR_ORDERED_FIELDS + ['Pv']:
            if getattr(self, p) is not None:
                reporter.add_property(p, self._from_table_units(getattr(self, p), self.get_unit(p)), self.get_unit(p), self.get_formatter(p))
        if self.x is not None:
            reporter.add_property('x', self.x, f'{self.mass_unit} vapor/{self.mass_unit} total')
            for phase, state in [('L', self.Liquid), ('V', self.Vapor)]:
                for p in self._STATE_VAR_ORDERED_FIELDS + ['Pv']:
                    if not p in 'TP':
                        if getattr(state, p) is not None:
                            reporter.add_property(f'{p}{phase}', self._from_table_units(getattr(state, p), self.get_unit(p)), self.get_unit(p), self.get_formatter(p))
        return reporter.report()

    def __repr__(self):
        self.lookup()
        return f'State(T={self.T}, P={self.P}, v={self.v}, u={self.u}, h={self.h}, s={self.s}, x={self.x})'

    def clone(self, **kwargs) -> State:
        """ Create a copy of this State instance """
        new_state = State()
        for f in self._STATE_VAR_FIELDS.union({'x'}).union(self._UNIT_FIELDS):
            setattr(new_state, f, getattr(self, f))
        for k, v in kwargs.items():
            if k in self._STATE_VAR_FIELDS.union({'x'}).union(self._UNIT_FIELDS):
                setattr(new_state, k, v)
        return new_state

    def _resolve_unsatd(self, specs: list[str]):
        """ 
        Resolve the thermodynamic state of steam/water given specifications
        """
        hasT = 'T' in specs
        hasP = 'P' in specs
        if hasT and hasP:
            """ T and P given explicitly """
            self._resolve_at_T_and_P()
        elif hasT or hasP:
            """ T OR P given, along with some other property (v,u,s,h) """
            self._resolve_at_TorP_and_Theta(specs)
        else:
            self._resolve_at_Theta1_and_Theta2(specs)

    def _resolve_at_T_and_P(self):
        """ T and P are both given explicitly.  Could be either superheated or subcooled state """
        specdict = {'T': self.TC, 'P': self.PMPa}

        if self.satd.lim['T'][0] < self.TC < self.satd.lim['T'][1]:
            Psat = self.satd.interpolators['T']['P'](self.TC)
            # print(f'Returns Psat of {Psat}')
        else:
            Psat = LARGE
        if self.P > Psat:
            ''' P is higher than saturation: this is a subcooled state '''
            retdict = self.subc.Bilinear(specdict)
        else:
            ''' P is lower than saturation: this is a superheated state '''
            retdict = self.suph.Bilinear(specdict)
        for p, v in retdict.items():
            if p not in specdict and p != 'x':
                setattr(self, p, v)
    
    def _resolve_at_TorP_and_Theta(self, specs: list[str]):
        """ T or P along with some other property (v,u,s,h) are specified """
        hasT = 'T' in specs
        hasP = 'P' in specs

        if not (hasT or hasP):
            raise ValueError('Either T or P must be specified along with another property')

        is_superheated = False
        is_subcooled = False
        supercritical = False
        if hasT:
            p = 'T'
            v = self.TC
            supercritical = v >= TCC
        else:
            p = 'P'
            v = self.PMPa
            supercritical = v >= PCMPA

        op = specs[0] if specs[1] == p else specs[1]
        th = getattr(self, op)

        if not supercritical:
            thL = self.satd.interpolators[p][f'{op}L'](v)
            thV = self.satd.interpolators[p][f'{op}V'](v)
            if th < thL:
                is_subcooled = True
            elif th > thV:
                is_superheated = True
            else:
                raise ValueError(f'Specified state is saturated based on {p}={v} and {op}={th}')
        else:
            if hasT:
                is_superheated = True
            else:
                is_subcooled = True

        if not is_superheated and not is_subcooled:
            raise ValueError(f'Specified state is saturated based on {p}={v} and {op}={th}')
        specdict = {p: v, op: th}
        if is_superheated:
            retdict = self.suph.Bilinear(specdict)
        else:
            retdict = self.subc.Bilinear(specdict)
        for p, v in retdict.items():
            if p not in specs and p != 'x':
                setattr(self, p, v)

    def _resolve_at_Theta1_and_Theta2(self, specs: list[str]):
        specdict = {specs[0]: getattr(self, specs[0]), specs[1]: getattr(self, specs[1])}
        try:
            sub_try = self.subc.Bilinear(specdict)
        except Exception as e:
            logger.debug(f'Subcooled Bilinear failed: {e}')
            sub_try = None
        try:
            sup_try = self.suph.Bilinear(specdict)
        except Exception as e:
            logger.debug(f'Superheated Bilinear failed: {e}')
            sup_try = None
        if sub_try and not sup_try:
            retdict = sub_try
        elif sup_try and not sub_try:
            retdict = sup_try
        elif sup_try and sub_try:
            raise ValueError(f'Specified state is ambiguous between subcooled and superheated states based on {specs}')
        else:
            raise ValueError(f'Specified state could not be resolved as either subcooled or superheated based on {specs}')
        logger.debug(f'Resolved state with {retdict}')
        for p, v in retdict.items():
            if p not in specs and p != 'x':
                setattr(self, p, v)

    def _resolve_satd(self, specs: list[str]):
        """
        Resolve an explicitly saturated state given specifications
        
        Parameters
        ----------
        specs: list[str]
            List of specified properties, must include 'x' and either 'T' or 'P'
        """

        if 'x' in specs:
            """ Vapor fraction is explicitly given """
            if 'T' in specs or 'P' in specs:
                """ Vapor fraction and one of T or P is given """
                p = 'T' if 'T' in specs else 'P'
                v = getattr(self, p)
                complement = 'P' if p == 'T' else 'T'
                prop = self.satd.interpolators[p][complement](v)
                setattr(self, complement, prop)
                self.Liquid = self.clone(x=0.0)
                self.Vapor = self.clone(x=1.0)
                for op in self._STATE_VAR_FIELDS - {'T','P'}:
                    for phase, state in [('L', self.Liquid), ('V', self.Vapor)]:
                        prop = self.satd.interpolators[p][f'{op}{phase}'](v)
                        setattr(state, op, prop)
                    setattr(self, op, self.x * getattr(self.Vapor, op) + (1 - self.x) * getattr(self.Liquid, op))
            else:
                """ Vapor fraction and one lever-rule-calculable property (u, v, s, h) is given """
                th = getattr(self, p)
                Y = np.array(self.satd.DF['T'][f'{p}V']) * self.x + np.array(self.satd.DF['T'][f'{p}L']) * (1 - self.x)
                X = np.array(self.satd.DF['T']['T'])
                f = svi(interp1d(X, Y))
                try:
                    self.T = f(th)
                    self.P = self.satd.interpolators['T']['P'](self.T)
                    self.Liquid = self.clone(x=0.0)
                    self.Vapor = self.clone(x=1.0)
                    for op in self._STATE_VAR_FIELDS - {'T','P'}:
                        for phase, state in [('L', self.Liquid), ('V', self.Vapor)]:
                            prop = self.satd.interpolators['T'][f'{op}{phase}'](self.T)
                            setattr(state, op, prop)
                        if op != th:
                            setattr(self, op, self.x * getattr(self.Vapor, op) + (1 - self.x) * getattr(self.Liquid, op))
                except:
                    raise Exception(f'Could not interpolate {p} = {th} at quality {self.x} from saturated steam table')
        else:
            """ T or P along with a lever-rule-calculable property (u, v, s, h) is given """
            p = 'T' if 'T' in specs else 'P'
            pcomp = 'P' if p == 'T' else 'T'
            v = getattr(self, p)
            op = specs[0] if specs[1] == p else specs[1]
            th = getattr(self, op)
            thL = self.satd.interpolators[p][f'{op}L'](v)
            thV = self.satd.interpolators[p][f'{op}V'](v)
            th_list = [thL, thV]
            th_list.sort()
            if th_list[0] < th < th_list[1]:
                """ This is a saturated state! Use interpolation to get saturation value of complement property and lever rule to get vapor fraction: """
                prop = self.satd.interpolators[p][pcomp](v)
                setattr(self, pcomp, prop)
                self.x = (th - thL) / (thV - thL)
                self.Liquid = self.clone(x=0.0)
                self.Vapor = self.clone(x=1.0)
                for phase, state in [('L', self.Liquid), ('V', self.Vapor)]:
                    prop = self.satd.interpolators[p][f'{op}{phase}'](v)
                    setattr(state, op, prop)
                for op2 in self._STATE_VAR_FIELDS - {'T','P', op}:
                    for phase, state in [('L', self.Liquid), ('V', self.Vapor)]:
                        prop = self.satd.interpolators[p][f'{op2}{phase}'](v)
                        setattr(state, op2, prop)
                    setattr(self, op2, self.x * getattr(self.Vapor, op2) + (1 - self.x) * getattr(self.Liquid, op2))
            else:
                raise ValueError(f'Specified property {op}={th} is not between saturated liquid ({thL}) and vapor ({thV}) values at {p}={v}')

    def _scalarize(self):
        """ Convert all properties to scalars (not np.float64) """
        for p in self._STATE_VAR_FIELDS.union({'x'}):
            val = getattr(self, p)
            if isinstance(val, np.float64):
                setattr(self, p, val.item())
        if hasattr(self, 'Liquid'):
            self.Liquid._scalarize()
        if hasattr(self, 'Vapor'):
            self.Vapor._scalarize()

    def delta(self, other: State) -> dict:
        """ Calculate property differences between this state and another state """
        self.lookup()
        other.lookup()
        delta_props = {}
        for p in self._STATE_VAR_FIELDS:
            val1 = getattr(self, p)
            val2 = getattr(other, p)
            if val1 is not None and val2 is not None:
                delta_props[p] = val2 - val1
        delta_props['Pv'] = other.Pv - self.Pv
        return delta_props

