"""
Solubility & Miscibility calculations CLI
"""

import argparse
import cCOSMO
from cosmopharm import SLE, LLE, COSMOSAC
from cosmopharm.utils import create_components, read_params
try:
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError(
        "This script requires matplotlib. Please install it with "
        "'pip install matplotlib'."
    )

def parse_args():
    parser = argparse.ArgumentParser(description='Solubility & Miscibility Calculations CLI')
    parser.add_argument('--names', nargs=2, default=['SIM', 'PLGA50'], help='Names of components e.g. SIM PLGA50')
    parser.add_argument('--params_file', type=str, default="data/sle/table_params.xlsx", help='Path to parameters file')
    parser.add_argument('--profile_path', type=str, default="./profiles/_import_methods/UD/", help='Path to COSMO profiles')
    parser.add_argument('--fraction', type=str, default='w', choices=['w', 'x'], help='Fraction type for plotting (w = weight, x = mole)')
    parser.add_argument('--output', type=str, default='basic_usage.png', help='Output file name for the plot')
    return parser.parse_args()

def main():
    args = parse_args()

    # Create components and add parameters for SLE calculation
    parameters = read_params(args.params_file)
    components = create_components(args.names, parameters)
    api, polymer = components[:2]

    # Initialize COSMO-SAC model
    db = cCOSMO.DelawareProfileDatabase(
        f"{args.profile_path}/complist.txt",
        f"{args.profile_path}/sigma3/")

    for name in args.names:
        iden = db.normalize_identifier(name)
        db.add_profile(iden)
    COSMO = cCOSMO.COSMO3(args.names, db)

    # Add extended functionality to COSMO-SAC model (e.g. free volume, .. etc.)
    actmodel = COSMOSAC(COSMO, mixture=components)
    actmodel.combinatorial = 'FV'  # Free-Volume (FV), Staverman-Guggenheim (SG)
    actmodel.dispersion = True  # Turn on/off the dispersion (optional)

    # Calculate SLE (solubility)
    sle = SLE(actmodel=actmodel)
    ideal = sle.solubility(mix_type='ideal', show_progress=True)  # ideal mixture (gamma=1)
    real = sle.solubility(mix_type='real', show_progress=True)  # real mixture (gamma=COSMO)

    # Calculate LLE (miscibility)
    lle = LLE(actmodel=actmodel)
    options = dict(max_gap=0.1, dT=30, exponent=2.1, max_gap_type='weight')
    miscibility = lle.miscibility(T=310, **options)

    # =============================================================================
    # Plot results
    # =============================================================================
    plt.figure()
    plt.plot(*ideal[[args.fraction, 'T']].values.T, 'r.-', label='SLE (ideal)')
    plt.plot(*real[[args.fraction, 'T']].values.T, 'k.-', label='SLE (real)')
    plt.plot(miscibility[args.fraction+'L1'], miscibility['T'], 'k.--', mfc='w', label='LLE')
    plt.plot(miscibility[args.fraction+'L2'], miscibility['T'], 'k.--', mfc='w')

    plt.xlim(0, 1)
    plt.ylim(300, 500)
    plt.title(f"{api} + {polymer}")
    plt.ylabel("T / K")
    xlabel = {'w': 'Weight', 'x': 'Mole'}
    plt.xlabel(f"{xlabel[args.fraction]} fraction {api}")
    plt.legend()

    # Save the figure to a file
    plt.savefig(args.output)
    plt.show()

if __name__ == "__main__":
    main()
