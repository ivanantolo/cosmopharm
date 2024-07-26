"""
== COSMOPharm CALCULATION TOOL =============================================
This script is designed as a companion tool for the
COSMOPharm package, aimed at facilitating the replication of results presented
in our manuscript published in *Mol. Pharm.*. Through the
application of the COSMO-SAC model, it provides a streamlined approach for
predicting phase behavior in drug-polymer and drug-solvent systems, crucial
for the development of pharmaceutical amorphous solid dispersions and related
applications.

Utilizing COSMOPharm's capabilities, this script guides users through the
process of modeling solubility, miscibility, and phase behavior in various
chemical systems, offering a hands-on experience in reproducing the findings
from our publication and exploring new systems and configurations.

Features:
- Direct application of COSMOPharm for compatibility prediction, solubility
  calculation, and phase behavior analysis in drug formulation processes.
- Dynamic selection from a comprehensive list of systems with available
  experimental data, leveraging open-source COSMO-SAC models.
- Visualization of solubility curves and miscibility data, allowing for direct
  comparison with experimental or theoretical predictions.
- Modular design for ease of replication and extension, promoting further
  research and development within the COSMO-SAC community.

Prerequisites:
- Ensure COSMOPharm and its dependencies are installed. Follow installation
  instructions provided in the COSMOPharm package README for guidance.

Usage:
    The script is executed directly and is configured internally. Modifications
    to the script or the use of specific functions within COSMOPharm should
    follow the documentation and examples provided as part of the COSMOPharm
    package.

Examples:
    # Execute the script to replicate manuscript results or explore new systems.
    python example_usage.py

Developed by Ivan Antolovic (ivan.antolovic@tu-berlin.de) and Martin Klajmon
(martin.klajmon@vscht.cz) in 2024, this script and the accompanying
COSMOPharm package aim to advance the field of pharmaceutical science through
the innovative application of COSMO-SAC models.

All extensions and modifications are thoroughly documented within the code,
facilitating easy navigation and customization for specific research needs.

================================================================================
"""

import os
import pandas as pd
from pathlib import Path
from itertools import product
try:
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError(
        "This script requires matplotlib. Please install it with "
        "'pip install matplotlib' or install the 'examples' extra of this "
        "package using 'pip cosmopharm[examples]'."
    )

import cCOSMO
from cosmopharm import SLE, LLE, COSMOSAC
from cosmopharm.utils import read_params, create_components

DIR_SLE_EXP = Path("data/sle/experimental")
DIR_LLE_RES = Path("data/lle/calculated")
PARAMS_FILE = Path("data/sle/table_params.xlsx")
PROFILES_UD = Path("profiles/_import_methods/UD/sigma3")
PROFILES_API = Path("profiles/pharmaceuticals/sigma")
PROFILES_POLY = Path("profiles/polymers/sigma")

CONFIGS = {
    'ideal'       : dict(mix_type='ideal'),
    'CS'          : dict(mix_type='real', comb=False, dsp=False),
    'CS_sg'       : dict(mix_type='real', comb='sg',  dsp=False),
    'CS_fv'       : dict(mix_type='real', comb='fv',  dsp=False),
    'CS_dsp'      : dict(mix_type='real', comb=False, dsp=True),
    'CS_sg_dsp'   : dict(mix_type='real', comb='sg',  dsp=True),
    'CS_fv_dsp'   : dict(mix_type='real', comb='fv',  dsp=True),
}


# =============================================================================
# MAIN
# =============================================================================
def main():
    print(
        "NOTE: The data stored in the .xlsx files located in 'data/sle' and "
        "'data/lle' are fixed and will not be overwritten by this script."
    )

    LOAD_LLE = True # Toggle between (True) and (False)
    # LOAD_LLE = False # Comment/Uncomment this line for quick toggle
    USE_PROMPT = True # Toggle between (True) and (False)
    # USE_PROMPT = False # Comment/Uncomment this line for quick toggle

    # Select systems and model configurations
    if USE_PROMPT:
        selection = make_selection()
        if selection is None:
            print("Exiting the program.")
            return
        systems, configs = selection
    else:
        # Manually select system & configuration w/o using the Prompt
        systems = ['SIM_PLGA50']
        configs = ['CS_sg_dsp']

    # Loading experimental data set for quick lookup during processing
    systems_with_exp_sle = set(import_systems_with_exp())

    # Ensure the directory for saving plots exists
    output_dir = './figures'
    os.makedirs(output_dir, exist_ok=True)

    for system in systems:
        print(f"\n=== {system.upper()} ===")
        names = system.split('_')

        # Create components and add parameters for SLE calculation
        parameters = read_params(PARAMS_FILE)
        mixture = create_components(names, parameters)

        # Initialize model
        actmodel = initialize_model(names, mixture)

        # Initialize equilibria
        sle = SLE(actmodel=actmodel)
        lle = LLE(actmodel=actmodel)

        for config in configs:
            # Model configuration
            actmodel.config = config
            actmodel.configuration(**CONFIGS[config])

            # Calculate SLE (solubility)
            solubility = sle.solubility(**CONFIGS[config], show_progress=True)

            # Calculate LLE (miscibility)
            if LOAD_LLE:
                miscibility = import_lle_results(system, config)
            else: # LOAD_LLE == False
                # Note: The code within this else block is still experimental.
                # It's in a beta stage and may not function correctly for all systems.
                # Significant manual adjustments may be necessary as it is not
                # yet fully automated like the solubility calculations.
                # If in doubt, switch to LOAD_LLE=True.
                try:
                    res = import_lle_results(system, config)
                    init = res.loc[res['T'].idxmin()]
                    T0, *x0 = init[["T", "xL1", "xL2"]]
                except ValueError:
                    # Provide T0, x0 || x0=None will approxmiate x0=f(T0)
                    T0, x0 = 310, None
                options = dict(max_gap=0.1, max_gap_type='weight', dT=30, exponent=2.1)
                miscibility = lle.miscibility(T=T0, x0=x0, **options)

            # Plot
            c, l = COLORS[config], LABELS[config]
            plt.plot(*solubility[['w','T']].values.T,'.-', c=c, label=l)
            plt.plot(*miscibility[['wL1','T']].values.T,'.--', c=c, mfc='w')
            plt.plot(*miscibility[['wL2','T']].values.T,'.--', c=c, mfc='w')

        # Experimental Data
        if system in systems_with_exp_sle:
            exp_sle = load_sle_exp_data(system)
            plt.plot(exp_sle['w'], exp_sle['T'],'ko',mfc='gray',label='Exp. data')

        # Plot Settings
        plt.xlim(0,1)
        plt.ylim(300,500)
        plt.title(f"{names[0]} + {names[1]}")
        plt.ylabel("T / K")
        xlabel = {'w':'Weight', 'x':'Mole'}
        plt.xlabel(f"{xlabel['w']} fraction {names[0]}")
        plt.legend()

        # Save the plot as a .png file
        plot_filename = os.path.join(output_dir, f"{system}_{config}.png")
        plt.savefig(plot_filename)
        plt.close()  # Close the plot to free memory

    print(f"\nPlots have been saved to the folder {output_dir}.")


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
def import_components(profiles_api, profiles_poly):
    api = [f.stem for f in profiles_api.glob('*.sigma')]
    poly = [f.stem for f in profiles_poly.glob('*.sigma')]
    return api, poly

def import_systems_with_exp(exp_dir=None):
    """Import systems for which experimental data are available."""
    exp_dir = exp_dir or DIR_SLE_EXP
    return [f.stem for f in exp_dir.iterdir()
            if f.is_file() and not f.name.startswith('_')]

def import_systems_all(profiles_api, profiles_poly):
    "Import systems using all combinations of API and polymers."
    api, poly = import_components(profiles_api, profiles_poly)
    return [f"{a}_{p}" for a, p in product(api, poly)]

def import_systems_without_exp(systems_all=None, systems_exp=None):
    """Calculate systems for which no experimental data are available."""
    if systems_all is None:
        systems_all = import_systems_all(PROFILES_API, PROFILES_POLY)
    if systems_exp is None:
        systems_exp = import_systems_with_exp(DIR_SLE_EXP)
    return [system for system in systems_all if system not in systems_exp]

def initialize_model(names, mixture):
    path_to_sigma_profiles = PROFILES_UD  # Path to sigma profiles
    path_to_complist = PROFILES_UD.parent # Path to complist.txt
    # Initialize COSMO-SAC model
    db = cCOSMO.DelawareProfileDatabase(
        (path_to_complist / "complist.txt").as_posix(),
        path_to_sigma_profiles.as_posix())
    for name in names:
        iden = db.normalize_identifier(name)
        db.add_profile(iden)
    COSMO = cCOSMO.COSMO3(names, db)
    # Extend functionality of COSMO-SAC model (e.g. free volume, .. etc.)
    # by wrapping COSMO into the new 'cosmopharm.COSMOSAC' class
    return COSMOSAC(COSMO, mixture=mixture)

def load_sle_exp_data(system, exp_dir=None):
    exp_dir = exp_dir or DIR_SLE_EXP
    # Ensure the system name is formatted correctly to match file naming conventions
    file = '_'.join(system.split('_')[:2])
    file_path = exp_dir / f"{file}.dat"  # Use Path object for path manipulation
    # Load the data using pandas
    exp = pd.read_csv(file_path, sep='\t', header=None,
                      skipfooter=1, usecols=[0, 1], engine='python')
    exp.columns = ['T', 'w']  # Assign column names
    return exp

def load_lle_res_data():
    pass

# =============================================================================
# Prompts for User Input
# =============================================================================
def get_initial_user_choice():
    """Prompt the user to choose which systems they want to select from."""
    print("\nPlease choose which systems you want to select from:")
    print()
    print("01. Systems, for which experimental data are available (reflects manuscript)")
    print("02. Systems, for which no experimental data are available in 'data/sle/experimental'.")
    print("03. All Systems regardless of experimental availability.")
    print()
    print("Enter '3' or '(a)ll' to select all systems.")
    print("Enter '(e)xit' to quit the program.")

    while True:
        choice = input("\nYour choice: ").strip().lower()
        if choice in ['1', '2']:
            return int(choice)
        elif choice in ['3', 'a']:
            return 3  # Treat 'a' as equivalent to selecting '3' for all systems
        elif choice in ['exit', 'e']:
            return None  # Indicate the user chose to exit
        else:
            print("Invalid selection. Please enter 1, 2, 3, 'a' for all, or 'e' to exit.")


def prompt_for_options(options, title="options", ncols=3):
    """
    Prompt the user for input, handle exit, and print selected options.

    Parameters:
    - options: A list of items to choose from, or dictionary.keys().
    - title: A string title for what is being selected (e.g., "Available systems").
    - ncols: Number of columns to display the options in.
    """
    options = options if isinstance(options, list) else list(options)
    if title == 'systems':
        systems_with_res_lle = set(import_systems_with_lle())
        systems_with_LLE = [option for option in options if option in systems_with_res_lle]
        options_list = [option + " *(LLE)" if option in systems_with_res_lle
                        else option for option in options]
    else:
        options_list = options.copy()

    while True:
        print(f"\nAvailable {title}:")
        num_options = len(options_list)
        rows = (num_options + ncols - 1) // ncols

        for row in range(rows):
            for col in range(ncols):
                idx = row + col * rows
                if idx < num_options:
                    print(f"{idx+1:02d}. {options_list[idx]:<20}", end=' ')
            print()

        if title == 'systems' and systems_with_LLE:
            print("\n*LLE (AAPS) found in our study. Data available in 'data/lle/calculated'.")
            print("Note 1: LLE results only provided for 'CS_sg_dsp' and 'CS_fv_dsp' configuration.")
            print("Note 2: LLE might be at lower temperatures. You can adjust plt.xlim() to make it visible.")
        print(f"\nEnter the numbers of the {title} you want to select, separated by commas or specify a range with a hyphen (e.g., 1,3,5-7).")
        print("Enter '(a)ll' to select all options.")
        print("Enter '(e)xit' to quit the program.")

        user_input = input("\nYour choice: ").strip().lower()

        if user_input in ['exit', 'e']:
            selected_options = None
        elif user_input in ['all', 'a']:
            selected_options = options
        else:
            selected_options = []
            input_elements = user_input.split(',')
            for element in input_elements:
                element_cleaned = element.strip()
                if '-' in element_cleaned:  # Check if it's a range
                    start, end = element_cleaned.split('-', 1)
                    if start.isdigit() and end.isdigit():
                        start_index = int(start) - 1
                        end_index = int(end)
                        if 0 <= start_index < num_options and 0 < end_index <= num_options:
                            selected_options.extend(options[start_index:end_index])
                        else:
                            print(f"Range out of bounds: {element_cleaned}")
                    else:
                        print(f"Invalid range: {element_cleaned}")
                elif element_cleaned.isdigit():
                    index = int(element_cleaned) - 1
                    if 0 <= index < num_options:
                        selected_options.append(options[index])
                    else:
                        print(f"Index out of range: {element_cleaned}")
                else:
                    print(f"Option not found: {element_cleaned}")
        if selected_options:
            print("\nYou selected:")
            for option in selected_options:
                print(option)
                # print()
        return selected_options

def make_selection():
    user_choice = get_initial_user_choice()

    if user_choice is None:
        return None

    if user_choice == 1:
        systems = import_systems_with_exp(DIR_SLE_EXP)
    elif user_choice == 3:
        systems = import_systems_all(PROFILES_API, PROFILES_POLY)
    elif user_choice == 2:
        systems_exp = import_systems_with_exp(DIR_SLE_EXP)
        systems_all = import_systems_all(PROFILES_API, PROFILES_POLY)
        systems = import_systems_without_exp(systems_all, systems_exp)

    # Filter systems whose components are not listed in "data/sle/table_params.xslx"
    parameters = read_params(PARAMS_FILE)
    comps = set(parameters['Component'])
    systems = [sys for sys in systems if all(c in comps for c in sys.split('_'))]

    selected_systems = prompt_for_options(systems, title='systems')
    if selected_systems is None:
        return None

    selected_configs = prompt_for_options(CONFIGS.keys(), title='configurations', ncols=1)
    if selected_configs is None:
        return None

    return selected_systems, selected_configs


# =============================================================================
# SETTINGS
# =============================================================================
COLORS = {
    'ideal'     : 'k',
    'CS'        : 'C5',
    'CS_sg'     : 'C0',
    'CS_fv'     : 'C4',
    'CS_dsp'    : 'C2',
    'CS_sg_dsp' : '#333333',
    'CS_fv_dsp' : '#D45500'
    }

LABELS = {
    'ideal'       : 'Ideal',
    'CS'          : r'CS',
    'CS_sg'       : r'CS$^{SG}$',
    'CS_fv'       : r'CS$^{FV}$',
    'CS_dsp'      : r'CS$_{dsp}$',
    'CS_sg_dsp'   : r'CS$^{SG}_{dsp}$',
    'CS_fv_dsp'   : r'CS$^{FV}_{dsp}$',
    }


def import_systems_with_lle(dir_lle_res=None):
    # Ensure dir_lle_res is a Path object for consistent path operations
    dir_lle_res = Path(dir_lle_res or DIR_LLE_RES)
    file = dir_lle_res / 'lle_results.xlsx'

    # Safely open the Excel file and extract sheet names
    with pd.ExcelFile(file) as xlsx:
        sheet_names = xlsx.sheet_names

    # Return sheet names as a set
    return set(sheet_names)

def import_lle_results(system=None, config=None, dir_lle_res=None):
    dir_lle_res = Path(dir_lle_res or DIR_LLE_RES)
    file = dir_lle_res / 'lle_results.xlsx'
    empty_frame = pd.DataFrame(columns=['T', 'wL1', 'wL2'])
    try:
        res = pd.read_excel(file, sheet_name=system, header=[0,1], index_col=0)
    except ValueError:
        return empty_frame
    if system:
        configs = res.columns.get_level_values(0)
        if config in configs:
            return res[config].dropna()
        elif config is None:
            return res
        else:
            return empty_frame
    else:
        return res if config is None else {k:df[config].dropna() for k, df in res.items()}



# =============================================================================
# EXECUTE MAIN FUNCTION
# =============================================================================
if __name__ == "__main__":
    main()
