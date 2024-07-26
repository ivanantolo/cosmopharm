import numpy as np
import pandas as pd
from scipy.optimize import fsolve, root
from numpy.typing import NDArray
from typing import Literal, Optional, Type, Union, List, Tuple, Generator, Dict

from ..components import Component
from ..actmodels import ActModel
from ..utils.spacing import spacing

NumericOrFrame = Union[float, List[float], Tuple[float, ...], NDArray[np.float64], pd.DataFrame]

class SLE:
    def __init__(self,
                 actmodel: Union[ActModel, Type[ActModel]],
                 mixture: Optional[List[Component]] = None) -> None:
        self.actmodel = actmodel
        self.mixture = mixture
        self._validate_arguments()
        # Assign 'solute' and 'solvent' based on order in 'mixture'
        # Default assignment can be changed in e.g. 'solubility()'
        self.solute, self.solvent = self.mixture

    def solubility(self,
                   solute: Optional[Component] = None,
                   solvent: Optional[Component] = None,
                   vary: Literal['T', 'w', 'auto'] = 'auto',
                   mix_type: Literal['ideal', 'real'] = 'real',
                   args: Optional[NumericOrFrame] = None,
                   init: Optional[NumericOrFrame] = None,
                   solver: Literal['root', 'fsolve'] = 'root',
                   show_progress=False, **kwargs):
        ''' Calculate solubility curve of solute in solvent.'''
        self.solute = solute or self.solute
        self.solvent = solvent or self.solvent
        self.vary, self.mix_type = vary, mix_type
        self.show_progress = show_progress
        self.config = getattr(self.actmodel, 'config', self.mix_type)
        if self.vary == 'auto':
            gen = self.auto_solve(solver)
        else:
            self._vary = self.vary
            args = self.set_args(args)
            init = self.set_x0(init)
            gen = self.solve_sle(args, init, solver)
        try:
            res = [k for k in gen]
            res = pd.DataFrame(res, columns=['T', 'x', 'vary', 'w'])
            res = res[['T', 'w', 'x', 'vary']]
            return res
        except self.actmodel.InvalidFreeVolumeParametersException as e:
            print(f"Warning: {e}")  # Inform the user
            return pd.DataFrame(columns=['T', 'w', 'x', 'vary'])


# =============================================================================
# MATHEMATICS
# =============================================================================
    def solve_sle(self, args: NDArray[np.float64], init: NDArray[np.float64],
                  solver: Literal['root', 'fsolve'] = 'root'
                  ) -> Generator[Dict[str, Union[float, str]], None, None]:
        # Check compatibility of the "init" values
        is_iterable = init.size > 1
        if is_iterable and not init.size == args.size:
            msg = 'The length of "init" must be the same as "args".'
            raise ValueError(msg)
        x0 = init
        # Setup solver and handle pure component case
        key, lock = ['T', 'x'] if self._vary == 'T' else ['x', 'T']
        solve = self.set_solver(solver=solver)
        args, pure_component = self._handle_pure_component(args)
        if pure_component: # no need to calculate pure component
            yield pure_component
        for i, arg in enumerate(args):
            x0 = init[i] if is_iterable else x0
            out = float(solve(x0, arg))
            x0 = out if not is_iterable else x0
            res = {key: arg, lock: out, 'vary': self._vary}
            res['w'] = self.actmodel._convert(res['x'])[0]
            text = (f"T={res['T']:.2f}", f"x={res['x']:.4f}", f"w={res['w']:.4f}")
            if self.show_progress:
                print(f'SLE ({self.config}): ', *text)
            yield res

    def auto_solve(self, solver: Literal['root', 'fsolve'] = 'root'):
        if self.show_progress:
            print()
            print(f"Calculating SLE ({self.config})...")
        # Start with varying 'w' until dTdw > THRESHOLD
        self._vary = 'w'
        args = self.set_args()
        x0 = self.set_x0()
        gen = self.solve_sle(args, x0, solver)
        previous = None
        for i, current in enumerate(gen):
            yield current
            if self._should_stop_generator(i, previous, current):
                break  # This will end the generator
            previous = current
        # Switch to varying 'T'
        self._vary = 'T'
        T0, x0 = current['T'], current['x']
        args = self.set_args(xmax=T0)[1:]  # exclude initial point
        gen = self.solve_sle(args, x0)
        yield from gen


# =============================================================================
# THERMODYNAMICS
# =============================================================================
    def ideal_mix(self, T):
        return np.exp(-self.gibbs_fusion(T))

    def real_mix(self, T, x):
        lngamma = self.actmodel.lngamma(T, x)[0]
        return np.log(x) + lngamma + self.gibbs_fusion(T)

    # Gibbs energy of fusion, i.e., the right-hand side of the solubility equation:
    def gibbs_fusion(self, T):
        T_fus = self.solute.T_fus
        H_fus = self.solute.H_fus
        Cp_fus_A = self.solute.Cp_fus_A
        Cp_fus_BT = self.solute.Cp_fus_BT

        R  = 8.314                 # J/(mol K)
        RT = R*T                   # J/mol
        A, B = Cp_fus_A, Cp_fus_BT
        G1 = H_fus*(1-T/T_fus)     # J/mol
        G2 = A * (T-T_fus) + 0.5*B*(T**2-T_fus**2)
        G3 = -T * (A * np.log(T/T_fus) + B*(T-T_fus))
        G_fus = G1 + G2 + G3       # J/mol
        return G_fus/RT


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================
    def set_args(self,
                 args: Optional[NumericOrFrame] = None,
                 xmin: Optional[float] = None,
                 xmax: Optional[float] = None,
                 dx: Optional[float] = None
                 ) -> NDArray[np.float64]:
        vary = self._vary
        # Determine argument values based on input data or generate
        # them based on range and type
        defaults = {
            'T': {'min': 310, 'max': self.solute.T_fus, 'step': 10},
            'w': {'min': 0.01, 'max': 1, 'step': 0.08}
        }
        mi = defaults[vary]['min'] if xmin is None else xmin
        ma = defaults[vary]['max'] if xmax is None else xmax
        dx = defaults[vary]['step'] if dx is None else dx

        if args is None:
            if self.vary != 'auto':  # auto_vary == False
                args = np.arange(ma, mi-dx, -dx)
                args[-1] = np.maximum(args[-1], mi)
            elif vary == 'T':  # auto_vary == True
                num, dT = 16, 175  # How many data points in this T-range
                num = int((ma-mi)/dT*num)  # fraction of points if dT smaller
                num = max(6, num)
                kwargs = dict(reverse=True, n=1.5)
                args = spacing(ma, mi, num, 'poly', **kwargs)
            else:  # vary == 'w'
                num = 16 if self.mix_type == 'ideal' else 21
                args = spacing(ma, mi, num, 'quadratic')
        args = np.asarray(args)
        args = args if vary != 'w' else self.actmodel._convert(args, to='mole')
        return args

    def set_x0(self, init: Optional[NumericOrFrame] = None) -> NDArray[np.float64]:
        vary = self._vary
        # Set up initial values based on the type of variable ('T' or 'w')
        if vary == 'T':
            x0 = 1. if init is None else self.actmodel._convert(init, to='mole')
        else:  # vary == 'w'
            x0 = self.solute.T_fus if init is None else init
        x0 = np.asarray(x0)
        return x0

    def set_solver(self, solver: Literal['root', 'fsolve'] = 'root'):
        vary, mix = self._vary, self.mix_type
        # Define the objective function (fobj) and the solver function (solve)
        # based on the mixture type (mix) and the variable type (vary)
        if mix == 'ideal' and vary == 'T':
            def fobj(x, T): return self.ideal_mix(T)
            def solve(x0, args): return fobj(x0, args)
        else:
            if mix == 'ideal':
                def fobj(T, x): return x - self.ideal_mix(T)
            elif vary == 'T':  # mix != 'ideal'
                def fobj(x, T): return self.real_mix(T, x)
            else:  # vary == 'w'
                def fobj(T, x): return self.real_mix(T, x)
            kwargs = dict(method='krylov', options={'maxiter': 5, 'xtol': 1e-3})
            if solver == 'fsolve':
                def solve(x0, args): return fsolve(fobj, x0, args)
            else:
                def solve(x0, args): return root(fobj, x0, args, **kwargs).x
        # Problem: "fsolve" and "root" return different types of np.arrays
        # (1) fsolve returns (1,) 1D array
        # (2) root returns () 0D array
        # Therefore, it is necessary to use float(solve(...)) to extract the
        # single value from the array, since solve()[0] does not work for root.
        return solve


# =============================================================================
# AUXILLIARY FUNCTIONS
# =============================================================================
    def _should_stop_generator(self, i, previous, current):
        THRESHOLD = 60
        if i > 1:  # ensuring there was a previous result
            dT = current['T'] - previous['T']
            dw = current['w'] - previous['w']
            return (dT / dw) > THRESHOLD
        return False  # If not enough elements, continue the generator

    def _handle_pure_component(self, args):
        res = {'T': self.solute.T_fus, 'x': 1, 'vary': self._vary, 'w': 1}
        if self._vary == 'T' and self.solute.T_fus in args:
            args = args[args != self.solute.T_fus]
            return args, res
        elif self._vary == 'w' and 1 in args:
            args = args[args != 1]
            return args, res
        return args, None

    def _validate_arguments(self):
        """Validate the arguments for the SLE class."""
        # TODO: Insert case where both actmodel and mixture are provided (check if acmodel.mixture == mixture, if not raise warning)
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
