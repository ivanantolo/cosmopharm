import pandas as pd
import numpy as np
from ..components import Component

def read_params(params_filepath):
    params = pd.read_excel(params_filepath)
    cols = [col.split('/')[0] for col in params.columns]
    params.rename(columns=dict(zip(params.columns, cols)), inplace=True)
    return params

def get_params_for_component(component, parameters):
    """Retrieve the parameters for the given component from the dataframe."""
    params = parameters[parameters['Component'] == component.name]
    if params.empty:
        raise ValueError(f"No data found for component {component.name}")
    return params.iloc[0]

def add_parameters(c, params):
    """Add properties of a component based on provided parameters."""
    KILOJOULE_TO_JOULE = 1e3
    c.Mw = params['Mw']  # g/mol
    c.T_fus = params['T_fus'] if params['T_fus'] > 0 else np.nan  # K
    c.H_fus = params['H_fus'] * KILOJOULE_TO_JOULE  # J/mol
    c.Cp_fus_A = np.nan_to_num(params['Cp_fus_a_fit'])  # J/(mol K)
    c.Cp_fus_BT = np.nan_to_num(params['Cp_fus_bT_fit'])  # J/(mol K²)
    c.v_298 = params['v298']  # cm³/mol
    c.v_hc = params['v_hc']  # cm³/mol

def create_components(names, parameters):
    components = []
    for name in names:
        component = Component(name)
        params = get_params_for_component(component, parameters)
        add_parameters(component, params)
        components.append(component)
    return components
