# Author: Cameron F. Abrams <cfa22@drexel.edu>
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field, fields
from scipy.interpolate import interp1d
from .satd import SaturatedSteamTables
from .unsatd import UnsaturatedSteamTables
from sandlermisc import statereporter
import pint
import logging

logger = logging.getLogger(__name__)

ureg = pint.UnitRegistry(autoconvert_offset_to_baseunit = True)

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
TCK = 647.3 * ureg('K')  # critical temperature in K
TCC = TCK.to('degC')

PCBAR = 221.20 * ureg('bar')  # critical pressure in bar
PCMPA = PCBAR.to('MPa')

@dataclass
class State:
    """
    Thermodynamic state of steam/water.
    """
    P: pint.Quantity = None
    """ Pressure """
    T: pint.Quantity = None
    """ Temperature """
    v: pint.Quantity = None
    """ Specific volume """
    u: pint.Quantity = None
    """ Specific internal energy """
    h: pint.Quantity = None
    """ Specific enthalpy """
    s: pint.Quantity = None
    """ Specific entropy """

    x: float = None
    """ Vapor fraction (quality) """
    L: State = None
    """ Liquid phase state when saturated """
    V: State = None
    """ Vapor phase state when saturated """

    _cache: dict = field(default_factory=dict, init=False, repr=False)
    """ Cache for computed properties and input state snapshot """

    _STATE_VAR_ORDERED_FIELDS = ['T', 'P', 'v', 's', 'h', 'u']
    """ Ordered fields that define the input state """

    _STATE_VAR_FIELDS = frozenset(_STATE_VAR_ORDERED_FIELDS).union({'x'})
    """ Fields that define the input state for caching purposes """

    def __post_init__(self):
        """ Requests access to global steam tables singleton """
        self.tables = get_tables()
        self.satd = self.tables['satd']
        self.suph = self.tables['suph']
        self.subc = self.tables['subc']
   
    def _dimensionalize(self, value):
        # update value to carry default units if necessary
        if (isinstance(value, float) or isinstance(value, int)) and name in self._STATE_VAR_FIELDS:
            # apply default units to raw numbers
            default_unit = self.get_default_unit(name)
            if default_unit is not None:
                value = value * ureg(default_unit)
        elif isinstance(value, pint.Quantity) and name in self._STATE_VAR_FIELDS:
            # convert any incoming pint.Quantity to default units
            value = value.to(self.get_default_unit(name))
        return value

    def __setattr__(self, name, value):
        """ 
        Generalized setattr to handle units conversion for state variables.
        """
        if name in self._STATE_VAR_FIELDS:
            value = self._dimensionalize(value)
            # if this is a state variable, it could be being set as an input or as
            # a computed property. If value is none and it is in the input set, remove it and 
            # blank all computed properties.
            if not hasattr(self, '_cache'):
                self._initialize_cache(name, value)
            if self._INPUT_STATE_VARS is not None:
                if value is None and name in self._INPUT_STATE_VARS:
                    self._INPUT_STATE_VARS.remove(name)
                    """ Input state variable removed: blank all computed properties; 
                    also blanks this variable since it is not in the input set anymore """
                    for f in self._STATE_VAR_FIELDS:
                        if f not in self._INPUT_STATE_VARS:
                            super().__setattr__(f, None)
                    return
                if len(self._INPUT_STATE_VARS) < 2:
                    self._INPUT_STATE_VARS.add(name)
                if len(self._INPUT_STATE_VARS) == 2:
                    """ Input state is set: resolve the state now """
                    self._resolve()
            else:
                ### initialize the input state vars set
                if value is not None:
                    self._INPUT_STATE_VARS = {name}
        else:
            super().__setattr__(name, value)

        # if name in self._STATE_VAR_FIELDS and hasattr(self, '_cache'):
        #     self._cache.clear()
        # if value is None:
        #     super().__setattr__(name, value)

        # super().__setattr__(name, value)

    def _get_current_input_state(self):
        """ Snapshot of current input field values. """
        return {field: getattr(self, field) for field in self._STATE_VAR_FIELDS}

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
        if self._INPUT_STATE_VARS is None:
            return self._STATE_VAR_FIELDS.union  # All inputs "changed" (initial state)
        
        current = self._get_current_input_state()
        changed = {k for k in self._STATE_VAR_FIELDS if current[k] != self._input_state.get(k)}
        return changed if changed else False

    def get_default_unit(self, field_name: str) -> str:
        """
        Get the default unit for a given field
        
        Parameters
        ----------
        field_name: str
            Name of the field

        Returns
        -------
        str: Default unit as a string that is understandable by pint
        """
        default_unit_map = {
            'P': 'MPa',
            'T': 'degC',
            'v': 'm**3 / kg',
            'u': 'kJ / kg',
            'h': 'kJ / kg',
            's': 'kJ / (kg * K)',
            'Pv': 'kJ / kg',
        }
        return default_unit_map.get(field_name)
    
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

    @property
    def Pv(self):
        """ Pressure * specific volume in kJ/kg """
        return (self.P.to('kPa') * self.v.to('m**3 / kg')).to('kJ / kg')

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
            v = self.T
            op = specs[0] if specs[1] == 'T' else specs[1]
            if v > TCC:
                logger.debug(f'Temperature {v} exceeds critical temperature {TCC}; cannot be saturated')
                return False
        elif hasP and has_only_T_or_P:
            p = 'P'
            v = self.P
            op = specs[0] if specs[1] == 'P' else specs[1]
            if v > PCMPA:
                logger.debug(f'Pressure {v} exceeds critical pressure {PCMPA}; cannot be saturated')
                return False
        if p is not None and op is not None:
            logger.debug(f'Checking saturation at {p}={v} for property {op}={getattr(self, op)}')
            logger.debug(f'Saturation limits for {p}: {self.satd.lim[p]}')
            logger.debug(f'Between? {self.satd.lim[p][0] <= v.m <= self.satd.lim[p][1]}')
            if not (self.satd.lim[p][0] <= v.m <= self.satd.lim[p][1]):
                logger.debug(f'Out of saturation limits for {p}={v}')
                return False
            op_valV = self.satd.interpolators[p][f'{op}V'](v.m)
            op_valL = self.satd.interpolators[p][f'{op}L'](v.m)
            if op_valL < getattr(self, op).m < op_valV or op_valV < getattr(self, op).m < op_valL:
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
        for f in self._STATE_VAR_FIELDS.union({'x'}):
            setattr(new_state, f, getattr(self, f))
        for k, v in kwargs.items():
            if k in self._STATE_VAR_FIELDS.union({'x'}):
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
        specdict = {'T': self.T.m, 'P': self.P.m}

        if self.satd.lim['T'][0] < self.T.m < self.satd.lim['T'][1]:
            Psat = self.satd.interpolators['T']['P'](self.T.m)
            # print(f'Returns Psat of {Psat}')
        else:
            Psat = LARGE
        if self.P.m > Psat:
            ''' P is higher than saturation: this is a subcooled state '''
            retdict = self.subc.Bilinear(specdict)
        else:
            ''' P is lower than saturation: this is a superheated state '''
            retdict = self.suph.Bilinear(specdict)
        for p, v in retdict.items():
            if p not in specdict and p != 'x':
                # interpolators return scalars in default units, so we
                # put units on them here
                setattr(self, p, v * ureg(self.get_default_unit(p)))
    
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
            v = self.T
            supercritical = v >= TCC
        else:
            p = 'P'
            v = self.P
            supercritical = v >= PCMPA

        op = specs[0] if specs[1] == p else specs[1]
        th = getattr(self, op)

        if not supercritical:
            thL = self.satd.interpolators[p][f'{op}L'](v.m)
            thV = self.satd.interpolators[p][f'{op}V'](v.m)
            if th.m < thL:
                is_subcooled = True
            elif th.m > thV:
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
        specdict = {p: v.m, op: th.m}
        if is_superheated:
            retdict = self.suph.Bilinear(specdict)
        else:
            retdict = self.subc.Bilinear(specdict)
        for p, v in retdict.items():
            if p not in specs and p != 'x':
                setattr(self, p, v * ureg(self.get_default_unit(p)))

    def _resolve_at_Theta1_and_Theta2(self, specs: list[str]):
        specdict = {specs[0]: getattr(self, specs[0]).m, specs[1]: getattr(self, specs[1]).m}
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
                setattr(self, p, v * ureg(self.get_default_unit(p)))

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
                prop = self.satd.interpolators[p][complement](v.m)
                setattr(self, complement, prop * ureg(self.get_default_unit(complement)))
                self.Liquid = self.clone(x=0.0)
                self.Vapor = self.clone(x=1.0)
                for op in self._STATE_VAR_FIELDS - {'T','P'}:
                    for phase, state in [('L', self.Liquid), ('V', self.Vapor)]:
                        prop = self.satd.interpolators[p][f'{op}{phase}'](v.m)
                        setattr(state, op, prop * ureg(self.get_default_unit(op)))
                    setattr(self, op, self.x * getattr(self.Vapor, op) + (1 - self.x) * getattr(self.Liquid, op))
            else:
                """ Vapor fraction and one lever-rule-calculable property (u, v, s, h) is given """
                th = getattr(self, p)
                Y = np.array(self.satd.DF['T'][f'{p}V']) * self.x + np.array(self.satd.DF['T'][f'{p}L']) * (1 - self.x)
                X = np.array(self.satd.DF['T']['T'])
                f = svi(interp1d(X, Y))
                try:
                    self.T = f(th.m)
                    self.P = self.satd.interpolators['T']['P'](self.T.m)
                    self.Liquid = self.clone(x=0.0)
                    self.Vapor = self.clone(x=1.0)
                    for op in self._STATE_VAR_FIELDS - {'T','P'}:
                        for phase, state in [('L', self.Liquid), ('V', self.Vapor)]:
                            prop = self.satd.interpolators['T'][f'{op}{phase}'](self.T.m)
                            setattr(state, op, prop * ureg(self.get_default_unit(op)))
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
            thL = self.satd.interpolators[p][f'{op}L'](v.m)
            thV = self.satd.interpolators[p][f'{op}V'](v.m)
            th_list = [thL, thV]
            th_list.sort()
            if th_list[0] < th.m < th_list[1]:
                """ This is a saturated state! Use interpolation to get saturation value of complement property and lever rule to get vapor fraction: """
                prop = self.satd.interpolators[p][pcomp](v.m)
                setattr(self, pcomp, prop * ureg(self.get_default_unit(pcomp)))
                self.x = (th.m - thL) / (thV - thL)
                self.Liquid = self.clone(x=0.0)
                self.Vapor = self.clone(x=1.0)
                for phase, state in [('L', self.Liquid), ('V', self.Vapor)]:
                    prop = self.satd.interpolators[p][f'{op}{phase}'](v.m)
                    setattr(state, op, prop * ureg(self.get_default_unit(op)))
                for op2 in self._STATE_VAR_FIELDS - {'T','P', op}:
                    for phase, state in [('L', self.Liquid), ('V', self.Vapor)]:
                        prop = self.satd.interpolators[p][f'{op2}{phase}'](v.m)
                        setattr(state, op2, prop * ureg(self.get_default_unit(op2)))
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

