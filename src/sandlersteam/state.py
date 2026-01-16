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
    name: str = None
    """ state name assigned upon creation (optional) """

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
    Liquid: State = None
    """ Liquid phase state when saturated """
    Vapor: State = None
    """ Vapor phase state when saturated """

    _STATE_VAR_ORDERED_FIELDS = ['T', 'P', 'v', 's', 'h', 'u']
    """ Ordered fields that define the input state """

    _STATE_VAR_FIELDS = frozenset(_STATE_VAR_ORDERED_FIELDS).union({'x'})
    """ Fields that define the input state for caching purposes; includes vapor fraction 'x' """

    _cache: dict = field(default_factory=dict, init=False, repr=False)
    """ Internal cache for tracking input variables and state completeness """

    def __new__(cls, *args, **kwargs):
        """Set up _cache before any field assignments."""
        logger.debug(f'Creating new State instance with args: {args}, kwargs: {kwargs}')
        instance = object.__new__(cls)
        object.__setattr__(instance, '_cache', {
            '_input_vars': [],
            '_is_specified': False,
            '_is_complete': False
        })
        logger.debug(f'Initialized _cache for State instance: {instance._cache}')
        return instance

    def _dimensionalize(self, name, value):
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
        """Custom attribute setter with input tracking and auto-resolution."""

        logger.debug(f'State {self.name}: __setattr__ called for {name} with value {value} (current value: {getattr(self, name, None)})')
        if (value is None or value == {} or value == []) and getattr(self, name, None) is not None:
            logger.debug(f'State {self.name}: Attempt to set {name} to None ignored to preserve resolved value {getattr(self, name)}')
            return   # Don't overwrite resolved value with None

        # Non-state variables set normally
        if name not in self._STATE_VAR_FIELDS:
            logger.debug(f'State {self.name}: Setting non-state variable {name} to {value}')
            object.__setattr__(self, name, value)
            return

        current_inputs = self._cache.get('_input_vars', [])
        if name in current_inputs:
            if value is not None:
                value = self._dimensionalize(name, value)
                object.__setattr__(self, name, value)
                logger.debug(f'State {self.name}: Updated input variable {name} to {value}')
                self._cache['_is_complete'] = False
                self._cache['_is_specified'] = len(current_inputs) == 2
                if self._cache['_is_specified']:
                    logger.debug(f'__set_attr__: State {self.name}: Now fully specified with inputs: {current_inputs}, resolving state.')
                    self._resolve()
            else:
                current_inputs.remove(name)
                self._cache['_input_vars'] = current_inputs
                self._blank_computed_state_vars()
                logger.debug(f'State {self.name}: Removed input variable {name}')
                self._cache['_is_specified'] = False
                self._cache['_is_complete'] = False
                # remove from input vars
                # blank all computed state vars
        else:
            if value is not None:
                value = self._dimensionalize(name, value)
                object.__setattr__(self, name, value)
                if not self._cache.get('_is_specified', False):
                    logger.debug(f'__set_attr__: State {self.name}: Adding new input variable {name} with value {value}')
                    current_inputs.append(name)
                    self._cache['_input_vars'] = current_inputs
                    self._cache['_is_specified'] = len(current_inputs) == 2
                    if self._cache['_is_specified']:
                        logger.debug(f'__set_attr__: State {self.name}: Now fully specified with inputs: {current_inputs}, resolving state.')
                        self._resolve()
                elif self._cache.get('_is_complete', False):
                    logger.debug(f'__set_attr__: State {self.name}: Already fully specified; resetting state with first new input variable {name}')
                    # already fully specified - setting a new input var - blank all computed vars
                    current_inputs = [name]
                    self._cache['_input_vars'] = current_inputs
                    self._cache['_is_specified'] = False
                    self._cache['_is_complete'] = False
                    self._blank_computed_state_vars()
                # set the value
                # if not _is_specified, add to input vars
                # check for specification completeness, and if so, resolve
            else:
                logger.debug(f'__set_attr__: State {self.name}: Setting non-input variable {name} to None')
                # setting a non-input var to None - just set it
                object.__setattr__(self, name, value)

    def __post_init__(self):
        logger.debug(f'__post_init__: State {self.name}: __post_init__ called, checking specification completeness.')
        if self._cache['_is_specified'] and not self._cache['_is_complete']:
            current_inputs = self._cache.get('_input_vars', [])
            logger.debug(f'__post_init__: State {self.name}: Now fully specified with inputs: {current_inputs}, resolving state.')
            self._resolve()
            self._cache['_is_complete'] = True

    def _blank_computed_state_vars(self):
        """
        Clear all computed state variables
        """
        for var in self._STATE_VAR_FIELDS:
            if var not in self._cache.get('_input_vars', []):
                object.__setattr__(self, var, None)
        self._cache['_is_complete'] = False


    def _resolve(self):
        """
        Resolve all state variables from the two input variables.
        This is where you'd call your steam tables and calculate everything.
        """
        if not self._cache.get('_is_specified', False):
            logger.debug(f'State {self.name}: State not fully specified; cannot resolve.')
            return
        
        # try:
        states_speced = self._cache.get('_input_vars', [])
        if len(states_speced) != 2:
            return
        
        if self._check_saturation(states_speced):
            self._resolve_satd(states_speced)
        else:
            self._resolve_unsatd(states_speced)
        self._scalarize()
        logger.debug(f'_resolve: State {self.name}: Successfully resolved state with inputs: {states_speced}')
        self._cache['_is_complete'] = True
        logger.debug(f'_resolve: State {self.name}: State resolution complete {self._cache["_is_complete"]}')         
        # except:
        #     raise Exception(f'Could not resolve state {self.name} with inputs: {self._cache.get("_input_vars", [])}')

    def __repr__(self):
        """Show which variables are inputs vs computed."""
        inputs = self._cache.get('_input_vars', [])
        parts = []
        for var in self._STATE_VAR_ORDERED_FIELDS + ['x']:
            val = getattr(self, var)
            if val is not None:
                marker = '*' if var in inputs else ''
                parts.append(f"{var}{marker}={val}")
        return f"State({self.name}: {', '.join(parts)})"

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

    def _check_saturation(self, specs: dict[str, float]) -> bool:
        satd = get_tables()['satd']
        """ Check if the specified state is saturated """
        if 'x' in specs:
            return True
        p, op = None, None
        hasT, hasP = 'T' in specs, 'P' in specs
        has_only_T_or_P = hasT ^ hasP
        if hasT and has_only_T_or_P:
            p, v = 'T', self.T
            if v > TCC:
                logger.debug(f'Temperature {v} exceeds critical temperature {TCC}; cannot be saturated')
                return False
            op = specs[0] if specs[1] == p else specs[1]
        elif hasP and has_only_T_or_P:
            p, v = 'P', self.P
            if v > PCMPA:
                logger.debug(f'Pressure {v} exceeds critical pressure {PCMPA}; cannot be saturated')
                return False
            op = specs[0] if specs[1] == p else specs[1]
        if p is not None and op is not None:
            logger.debug(f'Checking saturation at {p}={v} for property {op}={getattr(self, op)}')
            logger.debug(f'Saturation limits for {p}: {satd.lim[p]}')
            logger.debug(f'Between? {satd.lim[p][0] <= v.m <= satd.lim[p][1]}')
            if not (satd.lim[p][0] <= v.m <=  satd.lim[p][1]):
                logger.debug(f'Out of saturation limits for {p}={v}')
                return False
            op_val_satd_vapor = satd.interpolators[p][f'{op}V'](v.m)
            op_val_satd_liquid = satd.interpolators[p][f'{op}L'](v.m)
            op_val = getattr(self, op).m
            if op_val_satd_liquid < op_val < op_val_satd_vapor or op_val_satd_vapor < op_val < op_val_satd_liquid:
                return True
        return False

    # def _resolve(self):
    #     """ 
    #     Resolve the thermodynamic state of steam/water given specifications
    #     """
    #     spec = self._get_statespec()

    #     if self._check_saturation(spec):
    #         self._resolve_satd(spec)
    #     else:
    #         self._resolve_unsatd(spec)
    #     self._scalarize()
    #     self._input_state = self._get_current_input_state()

    def report(self):
        reporter = statereporter.StateReporter()
        for p in self._STATE_VAR_ORDERED_FIELDS + ['Pv']:
            if getattr(self, p) is not None:
                reporter.add_property(p, getattr(self, p).m, self.get_default_unit(p), self.get_formatter(p))
        if self.x is not None:
            reporter.add_property('x', self.x, f'mass fraction vapor')
            if 0 < self.x < 1:
                for phase, state in [('L', self.Liquid), ('V', self.Vapor)]:
                    for p in self._STATE_VAR_ORDERED_FIELDS + ['Pv']:
                        if not p in 'TP':
                            if getattr(state, p) is not None:
                                reporter.add_property(f'{p}{phase}', getattr(state, p).m, self.get_default_unit(p), self.get_formatter(p))
        return reporter.report()

    def _resolve_unsatd(self, specs: list[str]):
        """ 
        Resolve the thermodynamic state of steam/water given specifications
        """
        logger.debug(f'Resolving unsaturated state with specs: {specs}')
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
        satd = get_tables()['satd']
        suph = get_tables()['suph']
        subc = get_tables()['subc']

        if satd.lim['T'][0] < self.T.m < satd.lim['T'][1]:
            Psat = satd.interpolators['T']['P'](self.T.m)
            # print(f'Returns Psat of {Psat}')
        else:
            Psat = LARGE
        if self.P.m > Psat:
            ''' P is higher than saturation: this is a subcooled state '''
            retdict = subc.Bilinear(specdict)
        else:
            ''' P is lower than saturation: this is a superheated state '''
            retdict = suph.Bilinear(specdict)
        for p, v in retdict.items():
            if p not in specdict and p != 'x':
                # interpolators return scalars in default units, so we
                # put units on them here
                setattr(self, p, v * ureg(self.get_default_unit(p)))
    
    def _resolve_at_TorP_and_Theta(self, specs: list[str]):
        satd = get_tables()['satd']
        suph = get_tables()['suph']
        subc = get_tables()['subc']

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
            thL = satd.interpolators[p][f'{op}L'](v.m)
            thV = satd.interpolators[p][f'{op}V'](v.m)
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
            retdict = suph.Bilinear(specdict)
        else:
            retdict = subc.Bilinear(specdict)
        for p, v in retdict.items():
            if p not in specs and p != 'x':
                setattr(self, p, v * ureg(self.get_default_unit(p)))

    def _resolve_at_Theta1_and_Theta2(self, specs: list[str]):
        suph = get_tables()['suph']
        subc = get_tables()['subc']
        specdict = {specs[0]: getattr(self, specs[0]).m, specs[1]: getattr(self, specs[1]).m}
        try:
            sub_try = subc.Bilinear(specdict)
        except Exception as e:
            logger.debug(f'Subcooled Bilinear failed: {e}')
            sub_try = None
        try:
            sup_try = suph.Bilinear(specdict)
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
            List of specified properties
        """
        satd = get_tables()['satd']

        if 'x' in specs:
            """ Vapor fraction is explicitly given """
            p = 'x'
            other_p = specs[0] if specs[1] == p else specs[1]
            if other_p in ['T', 'P']:
                """ Vapor fraction and one of T or P is given """
                other_v = getattr(self, other_p)
                complement = 'P' if other_p == 'T' else 'T'
                complement_value_satd = satd.interpolators[other_p][complement](other_v.m)
                setattr(self, complement, complement_value_satd * ureg(self.get_default_unit(complement)))
                exclude_from_lever_rule = {'T', 'P', 'x'}
                exclude_from_single_phase_saturated_resolve = {'T', 'P', 'x'}
                initialize_single_phase_saturated_with = {other_p: other_v}
            else:
                """ Vapor fraction and one lever-rule-calculable property (u, v, s, h) is given """
                other_v = getattr(self, other_p)
                Y = np.array(satd.DF['T'][f'{other_p}V']) * self.x + np.array(satd.DF['T'][f'{other_p}L']) * (1 - self.x)
                X = np.array(satd.DF['T']['T'])
                f = svi(interp1d(X, Y))
                try:
                    self.T = f(other_v.m)
                    self.P = satd.interpolators['T']['P'](self.T.m) * ureg(self.get_default_unit('P'))
                except:
                    raise Exception(f'Could not interpolate {other_p} = {other_v} at quality {self.x} from saturated steam table')
                exclude_from_lever_rule = {'T', 'P', 'x', other_p}
                exclude_from_single_phase_saturated_resolve = {'T', 'P', other_p, 'x'}
                initialize_single_phase_saturated_with = {'T': self.T, 'P': self.P}
            # we have for sure determined T; either it was given or we just calculated it

        else:
            """ x is not in specs -- expect that T or P along with a lever-rule-calculable property (u, v, s, h) is given """
            hasT, hasP = 'T' in specs, 'P' in specs
            has_only_T_or_P = hasT ^ hasP
            if hasT and has_only_T_or_P:
                p, complement = 'T', 'P'
            elif hasP and has_only_T_or_P:
                p, complement = 'P', 'T'
            else:
                raise ValueError('Either T or P must be specified along with another property for saturated state without explicit x')
            v = getattr(self, p)
            complement_value_satd = satd.interpolators[p][complement](v.m)
            setattr(self, complement, complement_value_satd * ureg(self.get_default_unit(complement)))
            other_p = specs[0] if specs[1] == p else specs[1]
            other_v = getattr(self, other_p)
            other_v_Lsat = satd.interpolators[p][f'{other_p}L'](v.m)
            other_v_Vsat = satd.interpolators[p][f'{other_p}V'](v.m)
            self.x = (other_v.m - other_v_Lsat) / (other_v_Vsat - other_v_Lsat)
            exclude_from_lever_rule = {'T', 'P', 'x', other_p}
            exclude_from_single_phase_saturated_resolve = {'T', 'P', 'x', other_p}
            initialize_single_phase_saturated_with = {p: v}
        if 0.0 < self.x < 1.0:
            # generate the two saturated single-phase substates and apply lever rule
            # to resolve remaining properties of the overall state
            self.Liquid = State(x=0.0, name=f'{self.name}_L' if self.name else 'Saturated Liquid', **initialize_single_phase_saturated_with)
            self.Vapor = State(x=1.0, name=f'{self.name}_V' if self.name else 'Saturated Vapor', **initialize_single_phase_saturated_with)
            for op in self._STATE_VAR_FIELDS - exclude_from_lever_rule:
                setattr(self, op, self.x * getattr(self.Vapor, op) + (1 - self.x) * getattr(self.Liquid, op))
        elif self.x == 0.0:
            # This is a saturated liquid state, need to resolve all properties not already set
            for op in self._STATE_VAR_FIELDS - exclude_from_single_phase_saturated_resolve:
                prop = satd.interpolators[other_p][f'{op}L'](other_v.m)
                setattr(self, op, prop * ureg(self.get_default_unit(op)))
        elif self.x == 1.0:
            # This is a saturated vapor state, need to resolve all properties not already set
            for op in self._STATE_VAR_FIELDS - exclude_from_single_phase_saturated_resolve:
                prop = satd.interpolators[other_p][f'{op}V'](other_v.m)
                setattr(self, op, prop * ureg(self.get_default_unit(op)))

    def _scalarize(self):
        """ Convert all properties to scalars (not np.float64) """
        for p in self._STATE_VAR_FIELDS.union({'x'}):
            val = getattr(self, p)
            if isinstance(val, np.float64):
                setattr(self, p, val.item())
        if hasattr(self, 'Liquid') and self.Liquid is not None:
            self.Liquid._scalarize()
        if hasattr(self, 'Vapor') and self.Vapor is not None:
            self.Vapor._scalarize()

    def delta(self, other: State) -> dict:
        """ Calculate property differences between this state and another state """
        delta_props = {}
        for p in self._STATE_VAR_FIELDS:
            val1 = getattr(self, p)
            val2 = getattr(other, p)
            if val1 is not None and val2 is not None:
                delta_props[p] = val2 - val1
        delta_props['Pv'] = other.Pv - self.Pv
        return delta_props

