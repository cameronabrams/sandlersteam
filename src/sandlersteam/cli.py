# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import os
import sys

import argparse as ap
import numpy as np

from importlib.metadata import version

from .state import State, get_tables
from sandlermisc.statereporter import StateReporter


banner = r"""
   ____             ____       
  / __/__ ____  ___/ / /__ ____
 _\ \/ _ `/ _ \/ _  / / -_) __/
/___/\_,_/_//_/\_,_/_/\__/_/   
      ______                   
     / __/ /____ ___ ___ _     
    _\ \/ __/ -_) _ `/  ' \    
   /___/\__/\__/\_,_/_/_/_/  v""" + version("sandlersteam") + """

"""

def show_available_tables_subcommand(args):
    SteamTables = get_tables()
    print(f'  Saturated steam:')
    print(f'    T-sat: T from {SteamTables["satd"].lim["T"][0]} to {SteamTables["satd"].lim["T"][1]} C')
    print(f'             from {np.round(SteamTables["satd"].lim["T"][0] + 273.15,2)} to {np.round(SteamTables["satd"].lim["T"][1] + 273.15,2)} K')
    print(f'    P-sat: P from {SteamTables["satd"].lim["P"][0]} to {SteamTables["satd"].lim["P"][1]} MPa')
    print(f'             from {np.round(SteamTables["satd"].lim["P"][0]*10,2)} to {np.round(SteamTables["satd"].lim["P"][1]*10,2)} bar')
    print(f'  Superheated steam blocks:\nPressure (MPa) -> Temperatures (C):')
    for p in SteamTables["suph"].uniqs['P']:
        Tlist = SteamTables["suph"].data[SteamTables["suph"].data['P'] == p]['T'].to_list()
        print(f'    {p:>5.2f} ->', ', '.join([f"{x:>7.2f}" for x in Tlist]))
    print(f'  Subcooled liquid blocks:\nPressure (MPa) -> Temperatures (C):')
    for p in SteamTables["subc"].uniqs['P']:
        Tlist = SteamTables["subc"].data[SteamTables["subc"].data['P'] == p]['T'].to_list()
        print(f'    {p:>5.2f} ->', ', '.join([f"{x:>6.2f}" for x in Tlist]))

def state_subcommand(args):
    state = State()
    for p in 'TPhsuvx':
        v = getattr(args, p, None)
        if v is not None:
            setattr(state, p, v)
    report = state.report()
    print(report)

def delta_subcommand(args):
    state1 = State()
    state2 = State()
    for p in 'TPhsuvx':
        v = getattr(args, f'{p}1', None)
        if v is not None:
            setattr(state1, p, v)
        v = getattr(args, f'{p}2', None)
        if v is not None:
            setattr(state2, p, v)
    state_1 = state1.report()
    state_2 = state2.report()
    print(f"State-change calculations for Steam:")
    if args.show_states:
        print()
        two_states = ["State 1:                                      State 2:"]
        for line1, line2 in zip(state_1.splitlines(), state_2.splitlines()):
            two_states.append(f"{line1:<41s}     {line2}")
        print("\n".join(two_states))
        print()
        print("Property changes:")
    delta = state1.delta(state2)
    for p in ['T', 'P', 'v', 'u', 'h', 's']:
        if p in delta:
            val = delta[p]
            print(f'Δ{p}  = {val: 6g}')

def cli():
    subcommands = {
        'avail': dict(
            func = show_available_tables_subcommand,
            help = 'show available steam tables locations (T, P ranges)'
        ),
        'state': dict(
            func = state_subcommand,
            help = 'display thermodynamic state for given inputs'
        ),
        'delta': dict(
            func = delta_subcommand,
            help = 'calculate property differences between two states'
        )
    }
    parser = ap.ArgumentParser(
        prog='sandlersteam',
        description='Interact with steam tables in Sandler\'s textbook',
        epilog="(c) 2025, Cameron F. Abrams <cfa22@drexel.edu>"
    )
    parser.add_argument(
        '-b',
        '--banner',
        default=False,
        action=ap.BooleanOptionalAction,
        help='toggle banner message'
    )
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=f'sandlercubics version {version("sandlercubics")}',
        help='show program version and exit'
    )
    subparsers = parser.add_subparsers(
        title="subcommands",
        dest="command",
        metavar="<command>",
        required=True,
    )
    command_parsers={}
    for k, specs in subcommands.items():
        command_parsers[k] = subparsers.add_parser(
            k,
            help=specs['help'],
            add_help=False,
            formatter_class=ap.RawDescriptionHelpFormatter
        )
        command_parsers[k].set_defaults(func=specs['func'])
        command_parsers[k].add_argument(
            '--help',
            action='help',
            help=specs['help']
        )

    state_args = [
        ('P', 'pressure', 'pressure in MPa', float, True),
        ('T', 'temperature', 'temperature in C', float, True),
        ('x', 'quality', 'vapor quality (0 to 1)', float, False),
        ('v', 'specific_volume', 'specific volume in m3/kg', float, False),
        ('u', 'internal_energy', 'internal energy in kJ/kg', float, False),
        ('h', 'enthalpy', 'enthalpy in kJ/kg', float, False),
        ('s', 'entropy', 'entropy in kJ/kg-K', float, False),]

    for prop, longname, explanation, tp, _ in state_args:
        command_parsers['state'].add_argument(
            f'-{prop}',
            f'--{longname}',
            dest=prop,
            type=tp,
            help=f'{explanation.replace("_"," ")}'
        )
        command_parsers['delta'].add_argument(
            f'-{prop}1',
            f'--{longname}1',
            dest=f'{prop}1',
            type=tp,
            help=f'{explanation.replace("_"," ")} for state 1'
        )
        command_parsers['delta'].add_argument(
            f'-{prop}2',
            f'--{longname}2',
            dest=f'{prop}2',
            type=tp,
            help=f'{explanation.replace("_"," ")} for state 2'
        )
    command_parsers['delta'].add_argument(
        '--show-states',
        default=False,
        action=ap.BooleanOptionalAction,
        help='also show full states for both state 1 and state 2'
    )
    args = parser.parse_args()
    if args.func == state_subcommand:
        nprops = 0
        for prop, _, _, _, _ in state_args:
            if hasattr(args, prop) and getattr(args, prop) is not None:
                nprops += 1
        if nprops > 2:
            parser.error('At most two of P, T, x, v, u, h, s, and x may be specified for "state" subcommand')

    if args.banner:
        print(banner)
    if hasattr(args, 'func'):
        args.func(args)
    else:
        my_list = ', '.join(list(subcommands.keys()))
        print(f'No subcommand found. Expected one of {my_list}')
    if args.banner:
        print('Thanks for using sandlersteam!')