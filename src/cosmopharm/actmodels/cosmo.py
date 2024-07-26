import numpy as np
import cCOSMO

from typing import List, Union, Literal
from .actmodel import ActModel
from ..components import Component

class COSMOSAC(ActModel):

    # Handle invalid values for free volume calculation
    class InvalidFreeVolumeParametersException(Exception):
        pass

    def __init__(self,
                 COSMO: Union[cCOSMO.COSMO1, cCOSMO.COSMO3],
                 mixture: List[Component],
                 combinatorial: Union[Literal['sg', 'fv'], bool] = 'sg',
                 dispersion: bool = False,
                 ) -> None:
        self.COSMO = COSMO
        self.mixture = mixture
        # Flexible assignment of 'get_lngamma_comb' and 'get_lngamma_dsp'
        # that changes dynamically if the values for 'combinatorial' or
        # 'dispersion' are changed after initialization of an instance.
        self._combinatorial = combinatorial
        self._dispersion = dispersion

    @ActModel.vectorize
    def lngamma(self, T, x):
        resid = self.get_lngamma_resid(T, x)
        comb  = self.get_lngamma_comb(x)
        disp  = self.get_lngamma_disp(x)
        lngamma = resid + comb + disp
        return lngamma

    def get_lngamma_fv(self, x):
        """
        Calculates the free-volume term of the activity coefficient for a mixture.

        This implementation uses a formula to avoid numerical instability when
        `x_i` approaches zero, which is important in asymmetric API-polymer
        mixtures. The formula used is:

        ```
        phi_i^FV / x_i = v_i^F / sum_j(x_j * v_j^F)
        ```

        where
        - `phi_i^FV` is the free-volume fraction of component `i`,
        - `x_i` is the mole fraction of component `i`,
        - `v_i^F` is the free volume of component `i`,
        and the summation is over all components `j` in the mixture.

        Parameters
        ----------
        x : array_like
            Mole fractions of the components in the mixture.

        Returns
        -------
        np.ndarray
            Logarithm of the free-volume term of the activity coefficient.

        Note:
        Free-volume term of the activity coefficient according to Elbro et al.
        (can replace ln_gamma_comb of normal COSMO-SAC) - Kuo2013
        x, v_298, v_hc are 1D arrays (number of elements = number of components)
        """
        self.validate_free_volume_parameters()  # Ensure components are valid before proceeding
        v_298 = np.array([comp.v_298 for comp in self.mixture])
        v_hc = np.array([comp.v_hc for comp in self.mixture])
        vf = v_298-v_hc
        sum_vf = np.sum(x*vf)
        phix = vf/sum_vf
        return np.log(phix) + 1 - phix

    def get_lngamma_sg(self, x):
        return self.COSMO.get_lngamma_comb(0, x)

    def get_lngamma_resid(self, T, x):
        return self.COSMO.get_lngamma_resid(T, x)

    def get_lngamma_comb(self, x):
        if self._combinatorial is False:
            return np.zeros(len(x))
        elif self._combinatorial.lower() == 'sg':
            return self.get_lngamma_sg(x)
        elif self._combinatorial.lower() == 'fv':
            return self.get_lngamma_fv(x)

    def get_lngamma_disp(self, x):
        if self._dispersion:
            return self.COSMO.get_lngamma_disp(x)
        else:
            return np.zeros(len(x))

    @property
    def dispersion(self):
        return self._dispersion

    @dispersion.setter
    def dispersion(self, value):
        self._dispersion = value

    @property
    def combinatorial(self):
        return self._combinatorial

    @combinatorial.setter
    def combinatorial(self, value: Union[str, bool]):
        is_valid_string = isinstance(value, str) and value.lower() in ('sg', 'fv')
        is_False = value is False
        if is_valid_string or is_False:
            self._combinatorial = value
        else:
            msg = "Invalid value for combinatorial term. Please choose 'sg', 'fv', or set to False."
            raise ValueError(msg)

# =============================================================================
# Auxilliary functions
# =============================================================================
    def configuration(self,
                      comb: Union[Literal['sg', 'fv'], bool] = 'sg',
                      dsp: bool = False, **kwargs
                      ):
        """ Convenience function to quickly configure COSMO parameters """
        self._combinatorial = comb
        self._dispersion = dsp


    def validate_free_volume_parameters(self):
        # List of parameters to validate
        parameters_to_check = ["v_298", "v_hc"]

        for comp in self.mixture:
            invalid_params = []  # List to accumulate names of invalid parameters for this component
            for param in parameters_to_check:
                value = getattr(comp, param, None)
                # Check if value is None, not a number (np.nan), less than or equal to 0
                if value is None or np.isnan(value) or value <= 0:
                    invalid_params.append((param, value))  # Append parameter name and value tuple

            # Check if any errors were found for this component
            if invalid_params:
                # If errors were found, construct the warning message
                error_message = f"Invalid FV parameters for component {comp}: {invalid_params}"
                raise self.InvalidFreeVolumeParametersException(error_message)

            # Additionally check if v_298 and v_hc are equal
            if comp.v_298 == comp.v_hc:
                msg = f"v_298 and v_hc are equal for component {comp}: v_298={comp.v_298}, v_hc={comp.v_hc}"
                raise self.InvalidFreeVolumeParametersException(msg)
