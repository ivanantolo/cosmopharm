import numpy as np
import numbers
from numpy.typing import NDArray
from typing import List, Literal

from ..components import Component
from ..utils.convert import convert

class ActModel:

    def __init__(self, components: List[Component]):
        self.mixture = components

    def lngamma(self, T, x):
        raise NotImplementedError("lngamma() hasn't been implemented yet.")

    def activity(self, T, x):
        act = np.log(x) + self.lngamma(T, x)
        act[(act == np.inf) | (act == -np.inf)] = np.nan
        return act

    def gmix(self, T, x):
        is_scalar = np.isscalar(x)
        # Convert input as needed
        x = self._convert_input(x)
        # Create mask to identify columns that don't contain 0 or 1
        mask = np.any((x != 0) & (x != 1), axis=0)
        # Apply the mask to filter x
        _x = x[:, mask]
        # Calculate gmix for the  x values
        _gmix = _x * (np.log(_x) + self.lngamma(T, _x))
        _gmix = np.sum(_gmix, axis=0)
        # Initialize gmix array with zeros
        gmix = np.zeros(1 if x.ndim==1 else x.shape[1])
        # Fill gmix with calculated values where the mask is True
        gmix[mask] = _gmix
        return gmix[0] if is_scalar else gmix
    
    def thermofac(self, T, x):
        """ Approximate thermodynamic factor
        Simple derivative form, when no analytical equation is available. 
        """
        def f(x1):
            x = np.array([x1, 1-x1])
            return self.lngamma(T, x)[0]
        h, x = 0.0001, x[0]
        dy = (f(x+h)-f(x-h))/(2*h)
        return 1 + x * dy


# =============================================================================
# Wrapper functions (Decorators)
# =============================================================================
    @staticmethod
    def vectorize(func):
        ''' Intended vor ActModels where only single mole fractions can be
        handled, like e.g. COSMO-SAC. This function vectorizes the lngamma()
        to make it work with arrays of mole fractions.
        '''
        def wrapper(self, T, x):
            # Convert input to appropriate format
            x = self._convert_input(x)
            # Process based on the dimensionality of x
            if x.ndim == 1:
                return func(self, T, x)
            elif x.ndim == 2:
                results = [func(self, T, x[:, col]) for col in range(x.shape[1])]
                return np.array(results).T
            else:
                raise ValueError("Input must be either a scalar, 0D, 1D or 2D array")
        return wrapper


# =============================================================================
# Auxilliary functions
# =============================================================================
    def _convert_input(self, x):
        """Converts input to a 1-dim ndarray if it's a number or 0-dim ndarray."""
        if isinstance(x, numbers.Number) or (isinstance(x, np.ndarray) and x.ndim == 0):
            return np.array([float(x), 1 - float(x)])
        elif isinstance(x, np.ndarray) and x.ndim == 1 and len(x) != len(self.mixture):
            return np.array([x, 1 - x])
        return x

    def _convert(self, 
                 x : NDArray[np.float64], 
                 to : Literal['weight', 'mole'] ='weight'
                 ) -> NDArray[np.float64]:
        """
        Convert the fraction of a binary mixture between mole fraction and weight fraction.

        This method is designed for internal use with binary mixtures, where the mixture is defined by two components. 
        It uses the 'convert' function to perform the conversion by creating an array with the fractions of both 
        components and the molecular weights from the mixture's attributes.

        Parameters:
            x (NDArray[np.float64]): The mole or weight fraction of the first component of the mixture. 
                                    If converting 'to' weight, 'x' represents mole fractions; if converting 'to' mole, 
                                    'x' represents weight fractions. This should be a single value or a 1D array of values.
            to (Literal['weight', 'mole'], optional): The target type for the conversion. Defaults to 'weight'.
                                                    Use 'weight' to convert mole fractions to weight fractions,
                                                    and 'mole' to convert weight fractions to mole fractions.

        Returns:
            NDArray[np.float64]: The converted fraction(s) of the first component in the same shape as 'x'.
                                If 'x' is a single value, the return will be a single converted value;
                                if 'x' is a 1D array, the return will be a 1D array of converted values.

        Example:
            >>> mixture = Mixture(components=[component1, component2], Mw=np.array([18.01528, 46.06844]))
            >>> sle = SLE(mix=mixture)
            >>> x_mole_fraction = np.array([0.4])  # Mole fraction of the first component
            >>> x_weight_fraction = sle._convert(x_mole_fraction, to='weight')
            >>> print(x_weight_fraction)
            array([0.01373165])
        """
        Mw = np.array([c.Mw for c in self.mixture])
        return convert(x=np.array([x, 1-x], dtype=np.float64), Mw=Mw, to=to)[0]
