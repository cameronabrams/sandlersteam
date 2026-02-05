# Author: Cameron F. Abrams <cfa22@drexel.edu>
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field, fields
from scipy.interpolate import interp1d
from .satd import SaturatedSteamTables
from .unsatd import UnsaturatedSteamTables
from sandlermisc import ThermodynamicState, ureg, R
import pint
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

LARGE = 1e20
G_PER_MOL = 18.01528  # g/mol for water
KG_PER_MOL = G_PER_MOL / 1000.000  # kg/mol for water
MOL_PER_KG = 1.0 / KG_PER_MOL  # mol/kg for water
TCK = 647.3 * ureg.K  # critical temperature in K
TCC = TCK.to('degC')

PCBAR = 221.20 * ureg.bar  # critical pressure in bar
PCMPA = PCBAR.to('MPa')

@dataclass
class State(ThermodynamicState):
    """
    Thermodynamic state of steam/water.
    """

    name: str = 'Sandler-Steam'
    
    description: str = 'Sandler Steam Table State'

    _PARAMETER_ORDERED_FIELDS = ['Tc', 'Pc', 'Molwt']
    _PARAMETER_FIELDS = frozenset(_PARAMETER_ORDERED_FIELDS)

    Tc: pint.Quantity = TCC
    """ Critical temperature """
    Pc: pint.Quantity = PCMPA
    """ Critical pressure """
    Molwt: float = G_PER_MOL
    """ Molar weight in g/mol """
    
    def get_default_unit(self, var: str) -> pint.Unit:
        """
        Get the default unit for a given state variable.
        
        Parameters
        ----------
        var : str
            State variable name (e.g., 'P', 'T', 'v', 'u', 'h', 's', 'x')
        
        Returns
        -------
        pint.Unit
            Default unit for the variable
        """
        _default_unit_map = {
            'P': ureg.MPa,
            'T': ureg.degC,
            'v': ureg.m**3 / ureg.kg,
            'u': ureg.kJ / ureg.kg,
            'h': ureg.kJ / ureg.kg,
            's': ureg.kJ / (ureg.kg * ureg.K),
            'Pv': ureg.kJ / ureg.kg,
        }
        return _default_unit_map.get(var, ureg.dimensionless)

    def resolve(self) -> bool:
        """
        Resolve all state variables from the two input variables.
        This is where you'd call your steam tables and calculate everything.
        """
        if not self._cache.get('_is_specified', False):
            logger.debug(f'State {self.name}: State not fully specified; cannot resolve.')
            return False
        
        # try:
        states_speced = self.get_input_varnames()
        if len(states_speced) != 2:
            return False
        
        if self._check_saturation(states_speced):
            self._resolve_satd(states_speced)
        else:
            self._resolve_unsatd(states_speced)
        self._scalarize()
        logger.debug(f'resolve: State {self.name}: Successfully resolved state with inputs: {states_speced}')
        self._cache['_is_complete'] = True
        logger.debug(f'resolve: State {self.name}: State resolution complete {self._cache["_is_complete"]}')
        return True
        
    def _check_saturation(self, specs: dict[str, float | pint.Quantity]) -> bool:
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
                setattr(self, p, v * self.get_default_unit(p))
    
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
        logger.debug(f'resolve_at_TorP_and_Theta: Checking saturation at {p}={v} for property {op}={th}')
        if not supercritical:
            thL = satd.interpolators[p][f'{op}L'](v.m)
            thV = satd.interpolators[p][f'{op}V'](v.m)
            logger.debug(f'Saturation limits for {p}: {thL} to {thV}')
            if th.m < thL:
                is_subcooled = True
            elif th.m > thV:
                is_superheated = True
            else:
                raise ValueError(f'Specified state is saturated based on {p}={v} and {op}={th}')
        logger.debug(f'is_superheated: {is_superheated}, is_subcooled: {is_subcooled}, supercritical: {supercritical}')
        if not is_superheated and not is_subcooled and not supercritical:
            raise ValueError(f'Specified state is saturated based on {p}={v} and {op}={th}')
        specdict = {p: v.m, op: th.m}
        if is_superheated:
            retdict = suph.Bilinear(specdict)
        elif is_subcooled:
            retdict = subc.Bilinear(specdict)
        else:
            # if the temperature is greater than the lowest temperature available in the superheated table
            # at the given pressure, we treat it as superheated
            if hasT:
                Th_min_suph = suph.minTh_at_T(op, self.T.m)
                if th.m >= Th_min_suph:
                    retdict = suph.Bilinear(specdict)
                else:
                    retdict = subc.Bilinear(specdict)
            elif hasP:
                Th_min_suph = suph.minTh_at_P(op, self.P.m)
                if th.m >= Th_min_suph:
                    retdict = suph.Bilinear(specdict)
                else:
                    retdict = subc.Bilinear(specdict)

        for p, v in retdict.items():
            if p not in specs and p != 'x':
                setattr(self, p, v * self.get_default_unit(p))

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
                setattr(self, p, v * self.get_default_unit(p))

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
                setattr(self, complement, complement_value_satd * self.get_default_unit(complement))
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
                    self.P = satd.interpolators['T']['P'](self.T.m) * self.get_default_unit('P')
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
            setattr(self, complement, complement_value_satd * self.get_default_unit(complement))
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
                setattr(self, op, prop * self.get_default_unit(op))
        elif self.x == 1.0:
            # This is a saturated vapor state, need to resolve all properties not already set
            for op in self._STATE_VAR_FIELDS - exclude_from_single_phase_saturated_resolve:
                prop = satd.interpolators[other_p][f'{op}V'](other_v.m)
                setattr(self, op, prop * self.get_default_unit(op))


