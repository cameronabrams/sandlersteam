from unittest import TestCase
from sandlersteam.state import State

class TestSteamState(TestCase):

    def test_state_instance(self):
        state = State()
        self.assertIsInstance(state, State)
        self.assertIsNone(state.P)
        state.P = 1.0
        self.assertEqual(state.get_unit('P'), 'MPa')
        state.pressure_unit = 'kPa'
        self.assertEqual(state.get_unit('P'), 'kPa')
        self.assertEqual(state.PMPa, 0.001)
        self.assertTrue(hasattr(state, 'tables'))

    def test_state_clone(self):
        state = State() 
        state.P = 2.0
        state.T = 150.0
        clone_state = state.clone()
        self.assertIsInstance(clone_state, State)
        self.assertEqual(clone_state.P, state.P)
        self.assertEqual(clone_state.T, state.T)

        clone_state_2 = state.clone(pressure_unit='kPa', temperature_unit='C')
        self.assertEqual(clone_state_2.get_unit('P'), 'kPa')

    def test_state_get_statespec(self):
        state = State()
        state.P = 1.0
        state.T = 100.0
        pspec = state._get_statespec()
        self.assertIn('P', pspec)
        self.assertIn('T', pspec)
        self.assertEqual(len(pspec), 2)
        self.assertEqual(type(pspec), list)

    def test_state_invalid_statespec(self):
        state = State()
        state.T = 100.0
        with self.assertRaises(ValueError):
            state._get_statespec()

    def test_state_check_saturation(self):
        state = State()

        state.T = 300
        state.P = 5.0
        spec = state._get_statespec()
        self.assertFalse(state._check_saturation(specs=spec))

        state.T = None
        state.P = 0.1
        state.x = 0.5
        spec = state._get_statespec()
        self.assertTrue(state._check_saturation(specs=spec))

        state.x = None
        state.P = 0.1
        state.h = 2000.0
        spec = state._get_statespec()
        self.assertTrue(state._check_saturation(specs=spec))

        state.x = None
        state.P = None
        state.h = None
        state.s = 5.0
        state.T = 130.0
        spec = state._get_statespec()
        self.assertTrue(state._check_saturation(specs=spec))

    def test_state_resolve_superheated_exact_at_T_and_P(self):
        state = State(P=0.1, T=300)
        state.lookup()
        self.assertAlmostEqual(state.h, 3074.3, places=1)

        state = State(P=1.0, T=500)
        state.lookup()
        self.assertAlmostEqual(state.s, 7.7622, places=3)

        state = State(P=12.5, T=600)
        state.lookup()
        self.assertAlmostEqual(state.u, 3225.4, places=1)

    def test_state_resolve_superheated_interpolated_at_T_and_P(self):
        state = State(P=0.5, T=375)
        state.lookup()
        self.assertAlmostEqual(state.h, 3219.8, places=1)

        state = State(P=0.7, T=450)
        state.lookup()
        self.assertAlmostEqual(state.s, 7.787225, places=3)

    def test_state_resolve_superheated_at_T_and_h(self):
        state = State(T=400, h=3269)
        state.lookup()
        self.assertAlmostEqual(state.P, 0.68125, places=3)

    def test_state_resolve_superheated_at_P_and_s(self):
        state = State(P=2.0, s=6.5)
        state.lookup()
        self.assertAlmostEqual(state.T, 241.32848, places=2)

    def test_state_resolve_superheated_at_h_and_s(self):
        state = State(h=3500.0, s=7.5)
        state.lookup()
        self.assertAlmostEqual(state.T, 514.56, places=2)
        self.assertAlmostEqual(state.P, 2.0, places=2)

    def test_state_resolve_saturated_liquid(self):
        state = State(P=0.1, x=0)
        state.lookup()
        self.assertAlmostEqual(state.T, 99.63, places=2)

        state = State(T=150, x=0)
        state.lookup()
        self.assertAlmostEqual(state.P, 0.4758, places=2)

    def test_state_resolve_saturated_vapor(self):
        state = State(P=0.1, x=1)
        state.lookup()
        self.assertAlmostEqual(state.T, 99.63, places=2)
        self.assertAlmostEqual(state.Liquid.h, 417.46, places=2)
        self.assertAlmostEqual(state.Vapor.h, 2675.5, places=2)

        state = State(T=200, x=1)
        state.lookup()
        self.assertAlmostEqual(state.P, 1.5538, places=2)

    def test_state_resolve_saturated_at_T_and_h(self):
        state = State(T=120, h=1500.0)
        state.lookup()
        self.assertAlmostEqual(state.P, 0.19853, places=2)
        self.assertAlmostEqual(state.x, 0.452, places=2)