# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import argparse as ap

from .request import request_subcommand
from .state import state_subcommand, show_available_tables

def cli():
    subcommands = {
        'latex': dict(
            func = request_subcommand,
            help = 'make a latex steam table request'
        ),
        'avail': dict(
            func = show_available_tables,
            help = 'show available steam tables'
        ),
        'state': dict(
            func = state_subcommand,
            help = 'display thermodynamic state for given inputs'
        )
    }
    parser = ap.ArgumentParser(
        prog='sandlersteam',
        description='Interact with steam tables in Sandler\'s textbook',
    )
    parser.add_argument(
        '-b',
        '--banner',
        default=False,
        action=ap.BooleanOptionalAction,
        help='toggle banner message'
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
        command_parsers['state'].add_argument(
            f'-{prop}',
            f'--{longname}',
            dest=prop,
            type=tp,
            help=f'{explanation.replace("_"," ")}'
        )
    args = parser.parse_args()
    if args.TC == None and hasattr(args, 'T'):
        args.TC = args.T - 273.15
    elif args.T == None and hasattr(args, 'TC'):
        args.T = args.TC + 273.15
    nprops = 0
    for prop, _, _, _, _ in state_args:
        if hasattr(args, prop) and getattr(args, prop) is not None:
            nprops += 1
    if nprops > 2:
        parser.error('At most two of P, T, x, v, u, h, and s may be specified for "state" subcommand')

    if hasattr(args, 'func'):
        args.func(args)
    else:
        my_list = ', '.join(list(subcommands.keys()))
        print(f'No subcommand found. Expected one of {my_list}')