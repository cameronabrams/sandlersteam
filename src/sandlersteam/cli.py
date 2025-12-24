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
        description='Generate LaTeX-formatted steam tables'
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
            formatter_class=ap.RawDescriptionHelpFormatter
        )
        command_parsers[k].set_defaults(func=specs['func'])

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

    all_props = [('P','pressure in MPa'),
                 ('T','temperature in Celsius'),
                 ('x','quality (0 to 1)'),
                 ('v','specific volume in m3/kg'),
                 ('u','specific internal energy in kJ/kg'),
                 ('h','specific enthalpy in kJ/kg'),
                 ('s','specific entropy in kJ/kg-K')]
    for prop, explanation in all_props:
        command_parsers['state'].add_argument(
            f'--{prop}',
            type=float,
            help=f'{explanation}'
        )
    args = parser.parse_args()
    # only 2 of the 3 state variables P,T,x,V,U,H,S may be specified
    nprops = 0
    for prop, _ in all_props:
        if hasattr(args, prop) and getattr(args, prop) is not None:
            nprops += 1    
        if nprops > 2:
            parser.error('At most two of P,T,x,V,U,H,S may be specified for state command')

    if hasattr(args, 'func'):
        args.func(args)
    else:
        my_list = ', '.join(list(subcommands.keys()))
        print(f'No subcommand found. Expected one of {my_list}')