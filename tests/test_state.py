from unittest import TestCase
from sandlersteam.state import State
import pint
ureg = pint.UnitRegistry(autoconvert_offset_to_baseunit = True)
import logging
logger = logging.getLogger(__name__)

class TestSteamState(TestCase):

    def test_state_instance(self):
        state = State(name='test_state')
        self.assertIsInstance(state, State)
        self.assertIsNone(state.P)
        state.P = 1.0
        # default units are MPa
        self.assertEqual(state.P.dimensionality, ureg('MPa').dimensionality)

    def test_state_set_attributes(self):
        state = State(P=20.0*ureg('kPa'), T=150.0*ureg('degC'))
        # converts to default units
        self.assertEqual(state.P.m, 0.02)
        state.P = 3.0*ureg('bar')
        self.assertEqual(state.P.m, 0.3)
        another_state = State(P=5.0, T=200.0*ureg('degF'))  # default units
        self.assertAlmostEqual(another_state.T.m, 93.3333333)

    def test_state_check_spec(self):
        state = State(name='test_state')
        self.assertFalse(state._cache['_is_specified'])
        self.assertFalse(state._cache['_is_complete'])
        state.T = 300
        self.assertFalse(state._cache['_is_specified'])
        self.assertFalse(state._cache['_is_complete'])
        state.P = 5.0
        self.assertTrue(state._cache['_is_specified'])
        logger.debug(f'Test: State {state.name} is now fully specified, checking completeness.')
        logger.debug(f'cache: {state._cache}')
        self.assertTrue(state._cache['_is_complete'])
        self.assertIsNotNone(state.h)
        self.assertIsNotNone(state.u)
        self.assertIsNotNone(state.v)
        self.assertIsNotNone(state.s)
        self.assertIsNone(state.x)
        state.x = 0.5
        self.assertFalse(state._cache['_is_specified'])
        self.assertFalse(state._cache['_is_complete'])
        state.P = 3.6
        self.assertTrue(state._cache['_is_specified'])
        self.assertTrue(state._cache['_is_complete'])

    def test_state_resolve_superheated_exact_at_T_and_P(self):
        state = State(P=0.1, T=300, name='test_state')
        self.assertAlmostEqual(state.h.m, 3074.3, places=1)

        state = State(P=1.0, T=500, name='test_state')
        self.assertAlmostEqual(state.s.m, 7.7622, places=3)

        state = State(P=12.5, T=600, name='test_state')
        self.assertAlmostEqual(state.u.m, 3225.4, places=1)

    def test_state_resolve_superheated_interpolated_at_T_and_P(self):
        state = State(P=0.5, T=375, name='test_state')
        self.assertAlmostEqual(state.h.m, 3219.8, places=1)

        state = State(P=0.7, T=450, name='test_state')
        self.assertAlmostEqual(state.s.m, 7.787225, places=3)

    def test_state_resolve_superheated_at_T_and_h(self):
        state = State(T=400, h=3269, name='test_state')
        self.assertAlmostEqual(state.P.m, 0.68125, places=3)

    def test_state_resolve_superheated_at_P_and_s(self):
        state = State(P=2.0, s=6.5, name='test_state')
        self.assertAlmostEqual(state.T.m, 241.32848, places=2)

    def test_state_resolve_superheated_at_h_and_s(self):
        state = State(h=3500.0, s=7.5, name='test_state')
        self.assertAlmostEqual(state.T.m, 514.56, places=2)
        self.assertAlmostEqual(state.P.m, 2.0, places=2)

    def test_state_resolve_saturated_liquid(self):
        state = State(P=0.1, x=0, name='test_state')
        self.assertAlmostEqual(state.T.m, 99.63, places=2)

        state = State(T=150, x=0, name='test_state')
        self.assertAlmostEqual(state.P.m, 0.4758, places=2)

    def test_state_resolve_saturated_vapor(self):
        state = State(P=0.1, x=0.5, name='test_state')
        self.assertAlmostEqual(state.T.m, 99.63, places=2)
        self.assertAlmostEqual(state.Liquid.h.m, 417.46, places=2)
        self.assertAlmostEqual(state.Vapor.h.m, 2675.5, places=2)

        state = State(T=200, x=1, name='test_state')
        self.assertAlmostEqual(state.P.m, 1.5538, places=2)
        self.assertIsNone(state.Liquid)
        self.assertIsNone(state.Vapor)

    def test_state_resolve_saturated_at_T_and_h(self):
        state = State(T=120, h=1500.0, name='test_state')
        self.assertAlmostEqual(state.P.m, 0.19853, places=2)
        self.assertAlmostEqual(state.x, 0.452, places=2)