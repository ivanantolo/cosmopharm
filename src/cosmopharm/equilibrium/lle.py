import numpy as np
import pandas as pd
from scipy.optimize import least_squares, root
from typing import Union, Optional, Type, List, Dict

from ..components import Component
from ..actmodels import ActModel
from ..utils.spacing import spacing
from ..utils.lle_scanner import estimate_lle_from_gmix


class LLE:
    def __init__(self,
                 actmodel: Union[ActModel, Type[ActModel]],
                 mixture: Optional[List[Component]] = None) -> None:
        self.actmodel = actmodel
        self.mixture = mixture
        self._validate_arguments()

    def miscibility(self,
                    T: float,
                    x0: np.ndarray = None,
                    x0_type: str = 'mole',
                    max_gap: float = 0.1,
                    max_gap_type: str = 'mole',
                    max_T: float = 1000,
                    dT: float = 10,
                    exponent: float = 1
                    ) -> pd.DataFrame:
        """ Calculate miscibility """
        self.config = getattr(self.actmodel, 'config', '')
        Mw = np.array([c.Mw for c in self.mixture])
        self.is_valid_Mw = self.is_valid_numpy_array(Mw)
        res = {'binodal':[], 'spinodal':[]}
        var = 'x' if max_gap_type == 'mole' else 'w'
        print()
        print("Calculating LLE...")

        # Define column names
        binodal_columns = ['T', 'xL1', 'xL2']
        if self.is_valid_Mw:
            binodal_columns += ['wL1', 'wL2']

        # Check for valid molar masses
        if x0_type == 'weight':
            if self.is_valid_Mw:
                x0 = self.convert_to_mole_fractions(x0, self.mixture.Mw)
            else:
                raise ValueError("Molar masses are not available for conversion from weight to mole fraction.")
# =============================================================================
# TODO: Implement all edge cases (no LLE, bad approximation, ....)
# TODO: Improve code structure
        # Approximate initial value for LLE
        if x0 is None:
            print("...searching for suitable initial value...")
            x0 = self.approx_init_x0(T)

            if any(x is None for x in x0):
                print("...no initial value at T0 was found. Try another T0.")
                return pd.DataFrame(columns=binodal_columns)

            # Check if initial guess is reasonable - otherwise increase T
            # TODO: Check whether it might be an LCST if isLower
            binodal = self.solve_lle(T, x0, show_output=False)
            isEqual = np.diff(binodal['x'])[0] < 1e-8 # check if both phases have equal composition
            if isEqual:
                print("...no initial value at T0 was found. Try another T0.")
                return pd.DataFrame(columns=binodal_columns)
            isLower = min(binodal['x']) < min(x0) # lower bound below min(x0)
            while isLower and T <= max_T:
                print('LLE: ', f"{T=:.2f}", "...no feasbible initial value found.")
                T += 10  # Increase T by 10
                x0 = self.approx_init_x0(T)
                binodal = self.solve_lle(T, x0, show_output=False)
            print("Suitable initial value found! Proceed with calculating LLE...")
# =============================================================================
        # First iteration step
        binodal = self.solve_lle(T, x0)
        gap = np.diff(binodal[var])[0]
        res['binodal'].append((T, *[val for vals in binodal.values() for val in vals]))

        # Subsequent iteration steps
        while gap > max_gap and T <= max_T:
            T += dT * gap**exponent
            x0 = binodal['x']
            binodal = self.solve_lle(T, x0)
            gap = np.diff(binodal[var])[0]
            res['binodal'].append((T, *[val for vals in binodal.values() for val in vals]))

        # Convert lists to DataFrames
        res = pd.DataFrame(res['binodal'], columns=binodal_columns)
        return res

# =============================================================================
# MATHEMATICS
# =============================================================================

    def solve_lle(self, T: float, x0: np.ndarray, show_output=True) -> Dict[str, np.ndarray]:
        """ Solve for liquid-liquid equilibrium (LLE) at a given temperature and initial composition. """
        binodal = {'x': self.binodal(T, x0)}
        output = [f"{k}={v:.4f}" for k,v in zip(['xL1', 'xL2'], binodal['x'])]

        if self.is_valid_Mw:
            binodal['w'] = self.actmodel._convert(binodal['x'])
            output += [f"{k}={v:.4f}" for k,v in zip(['wL1', 'wL2'], binodal['w'])]

        if show_output:
            prefix = f'LLE ({self.config})' if self.config else 'LLE'
            print(f'{prefix}: ', f"{T=:.2f}", *output)
        return binodal




# =============================================================================
# THERMODYNAMICS
# =============================================================================
    def fobj_binodal(self, x1, T):
        # Equilibrium: Isoactivity criterion (aL1 - aL2 = 0)
        x = np.array([x1, 1-x1])
        activity = self.actmodel.activity(T, x)
        equilibrium = np.diff(activity, axis=1)
        return equilibrium.ravel() # reshape from (2,1) --> (2,)

    def fobj_spinodal(self, x1):
        T = 0
        x = np.array([x1, 1-x1])
        return self.actmodel.thermofac(T, x)

    def binodal(self, T, x0=None):
        if x0 is None:
            x0 = [0.1, 0.999]
        kwargs = dict(bounds=(0,1), ftol=1e-15, xtol=1e-15)
        res = least_squares(self.fobj_binodal, x0, args=(T,), **kwargs)
        return res.x

    def spinodal(self, T, x0=None):
        if x0 is None:
            x0 = self.binodal(T, x0)
        kwargs = dict(bounds=(0,1), ftol=1e-15, xtol=1e-15)
        res = least_squares(self.fobj_spinodal, x0, args=(T,), **kwargs)
        return res.x


# =============================================================================
# AUXILLIARY FUNCTIONS
# =============================================================================
    def approx_init_x0(self, T):
        x1 = spacing(0,1,51,'poly',n=3)
        gmix = self.actmodel.gmix(T, x1)
        xL, xR, yL, yR = estimate_lle_from_gmix(x1, gmix, rough=True)
        return xL, xR

    def _validate_arguments(self):
        """Validate the arguments for the LLE class."""
        # TODO: Insert case where both actmodel and mixture are provided
        # (check if acmodel.mixture == mixture, if not raise warning)
        if isinstance(self.actmodel, ActModel):
            # If actmodel is an instance of ActModel
            self.mixture: List[Component] = self.mixture or self.actmodel.mixture
        elif isinstance(self.actmodel, type) and issubclass(self.actmodel, ActModel):
            # If actmodel is a class (subclass of ActModel)
            if self.mixture is None:
                raise ValueError("Please provide a valid mixture:Mixture.")
            self.actmodel: ActModel = self.actmodel(self.mixture)
        else:
            # If actmodel is neither an instance nor a subclass of ActModel
            err = "'actmodel' must be an instance or a subclass of 'ActModel'"
            raise ValueError(err)

    def is_valid_numpy_array(self, arr: np.ndarray) -> bool:
        """Check if a numpy array contains only numbers and no None values."""
        if not isinstance(arr, np.ndarray):
            return False
        if arr.dtype == object:  # Check if the array contains objects (which could include None)
            return not np.any(arr == None)
        else:
            return np.issubdtype(arr.dtype, np.number)
