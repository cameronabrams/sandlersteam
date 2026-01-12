# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import os
import sys

import argparse as ap
import numpy as np

from importlib.metadata import version

from .request import request_subcommand
from .state import State
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

def regularize_T(args):
    # args is a namespace in which TC and TK must agree
    if args.TK is not None:
        args.TC = args.TK - 273.15
    elif args.TC is not None:
        args.TK = args.TC + 273.15
    else:
        raise ValueError('Either temperature in K or C must be specified')

def show_available_tables_subcommand(args):
    print(f'  Saturated steam:')
    print(f'    T-sat: T from {SteamTables["satd"].lim["TC"][0]} to {SteamTables["satd"].lim["TC"][1]} C')
    print(f'             from {np.round(SteamTables["satd"].lim["TC"][0] + 273.15,2)} to {np.round(SteamTables["satd"].lim["TC"][1] + 273.15,2)} K')
    print(f'    P-sat: P from {SteamTables["satd"].lim["P"][0]} to {SteamTables["satd"].lim["P"][1]} MPa')
    print(f'             from {np.round(SteamTables["satd"].lim["P"][0]*10,2)} to {np.round(SteamTables["satd"].lim["P"][1]*10,2)} bar')
    print(f'  Superheated steam blocks:\nPressure (MPa) -> Temperatures (C):')
    for p in SteamTables["suph"].uniqs['P']:
        Tlist = SteamTables["suph"].data[SteamTables["suph"].data['P'] == p]['TC'].to_list()
        print(f'    {p:>5.2f} ->', ', '.join([f"{x:>7.2f}" for x in Tlist]))
    print(f'  Subcooled liquid blocks:\nPressure (MPa) -> Temperatures (C):')
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

def request_subcommand(args):
    R = Request()
    if args.satdP:
        R.register('satdP')
    if args.satdT:
        R.register('satdT')
    if args.suphP:
        for P in args.suphP:
            R.register(suphP=P)
    if args.subcP:
        for P in args.subcP:
            R.register(subcP=P)
    with args.output as f:
        f.write(R.to_latex())
    print(f'Request completed: wrote output to {args.output.name if args.output else "stdout"}')


def cli():
    subcommands = {
        'latex': dict(
            func = request_subcommand,
            help = 'make a latex steam table request'
        ),
        'avail': dict(
            func = show_available_tables_subcommand,
            help = 'show available steam tables locations (T, P ranges)'
        ),
        'state': dict(
            func = state_subcommand,
            help = 'display thermodynamic state for given inputs'
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

    command_parsers['latex'].add_argument(
        '-o',
        '--output',
        type=ap.FileType('w'),
        default=None,
        help='output file (default: stdout)'
    )
    command_parsers['latex'].add_argument(
        '--suphP',
        type=float,
        action='append',
        help='add superheated steam table at pressure P (MPa)'
    )
    command_parsers['latex'].add_argument(
        '--subcP',
        type=float,
        action='append',
        help='add subcooled liquid table at pressure P (MPa)'
    )
    command_parsers['latex'].add_argument(
        '--satdP',
        action='store_true',
        help='include saturated steam table by pressure'
    )
    command_parsers['latex'].add_argument(
        '--satdT',
        action='store_true',
        help='include saturated steam table by temperature'
    )

    state_args = [
        ('P', 'pressure', 'pressure in MPa', float, True),
        ('T', 'temperature', 'temperature in K', float, True),
        ('x', 'quality', 'vapor quality (0 to 1)', float, False),
        ('v', 'specific_volume', 'specific volume in m3/kg', float, False),
        ('u', 'internal_energy', 'internal energy in kJ/kg', float, False),
        ('h', 'enthalpy', 'enthalpy in kJ/kg', float, False),
        ('s', 'entropy', 'entropy in kJ/kg-K', float, False),]
    extra_args = [
        ('TC', 'temperatureC', 'temperature in C (if T not specified)', float, True),
    ]
    for prop, longname, explanation, tp, _ in state_args + extra_args:
        if prop=='T':
            prop = 'TK'
        command_parsers['state'].add_argument(
            f'-{prop}',
            f'--{longname}',
            dest=prop,
            type=tp,
            help=f'{explanation.replace("_"," ")}'
        )
    args = parser.parse_args()
    if args.func == state_subcommand:
        regularize_T(args)
        nprops = 0
        for prop, _, _, _, _ in state_args:
            if hasattr(args, prop) and getattr(args, prop) is not None:
                nprops += 1
        if nprops > 2:
            parser.error('At most two of P, T, x, v, u, h, and s may be specified for "state" subcommand')

    if args.banner:
        print(banner)
    if hasattr(args, 'func'):
        args.func(args)
    else:
        my_list = ', '.join(list(subcommands.keys()))
        print(f'No subcommand found. Expected one of {my_list}')
    if args.banner:
        print('Thanks for using sandlersteam!')